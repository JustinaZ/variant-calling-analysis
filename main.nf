/*
========================================================================================
   ADRC DMS core Variant Calling Workflow with Nextflow, University of Washington
========================================================================================
   Github   : https://github.com/JustinaZ/variant-calling-analysis
   Contact  : jzurausk@uw.edu
----------------------------------------------------------------------------------------
*/

/*
========================================================================================
   P A R A M E T E R  block:
========================================================================================
*/



params.outdir = "${projectDir}/output"  // output directory
params.genome = "${projectDir}/reference_data/hg38.fa" // location where .fa is kept
params.dict = "${projectDir}/reference_data/hg38.dict" // location where .dict file will be created
params.bwt = "${projectDir}/reference_data/hg38.fa.bwt" // path to check if BWA index exists

params.db_path = "${projectDir}/output/genomics_db" // dir where genomic db to be created
params.tmp_dir = "${projectDir}/output/tmp" // dir where tmp files should be placed while creating genomic db
params.bed_file = "${projectDir}/data/bed_file/twist_plus_refseq_10-31-2018.grc38.bed" // location where the bed file is/should be kept

params.sample_gvcfs = "${projectDir}/data/gvcf_samples/1291919.merged.matefixed.sorted.markeddups.recal.g.vcf.gz"
 // this is single sample only!! To be used for creating initial Genomic DB, all other samples will be appended from sample mapping file

params.gvcf_files = "${projectDir}/data/gvcf_samples/*.g.vcf.gz" //  Path pattern to all GVCF files
params.sample_map = "${params.outdir}/sample_mapping_file/sample_map.txt" // Path for the sample map

// For joint variant calling
params.chromosomes = params.chromosomes = (1..22).collect { it.toString() } + ["X", "Y"]
params.joint_outdir = "${params.outdir}/joint_vcfs" // where to put chromosome based joind calls (for later merging)
params.memory = "2 GB"
params.cpus = 2



// Docker image to use for GATK
params.docker_image = 'broadinstitute/gatk:4.5.0.0' 

/*
========================================================================================
   C H A N N E L  definitions:
========================================================================================
*/

genome_ch = Channel.fromPath(params.genome, checkIfExists: true)
dict_ch = Channel.fromPath(params.dict, checkIfExists: true)
bed_file_ch = Channel.fromPath(params.bed_file, checkIfExists: true)


// input channel definition for gvcfs and their indices
gvcf_ch = Channel
    .from(params.sample_gvcfs)
    .map { gvcfPath ->
        def indexPath = "${gvcfPath}.tbi"  // Assumes the index has the same base name as the GVCF
        [file(gvcfPath), file(indexPath)]   // Creates a tuple of [GVCF file, index file]
    } // This channel should pass one GVCF and its corresponding index to the process, which will use that to create the genomic DB.


gvcf_files_ch = Channel.fromPath(params.gvcf_files, checkIfExists: true).collect() // Define the channel for all GVCF files (for sample mapping)


// gvcf files chANNEL to check if the .tbi file exists,
gvcf_files_ch_filtered = Channel
    .fromPath(params.gvcf_files, checkIfExists: true)
    .filter { gvcf_file -> !file("${gvcf_file}.tbi").exists() }  // to only include files without existing .tbi

//channel to split the processing by chromosome (aka parallel processing)
chr_ch = Channel.fromList(params.chromosomes) 
chr_ch.view { "Processing chromosome: $it" } // for test add .view() to check if all chrm are passed at once

// genomicDB channel
genomics_db_ch = Channel.fromPath(params.db_path, checkIfExists: true)

// Defining the VCF channel to emit file paths from the output directory
chr_vcfs_ch = Channel.fromPath("${projectDir}/output/joint_vcfs/*.vcf.gz", checkIfExists: true).collect()

println """\
         RUNNING ADRC DMS CORE V A R I A N T  C A L L I N G   P I P E L I N E
         ====================================================================
         genome        : ${params.genome}
         outdir        : ${params.outdir}
         gvcf_files    : ${params.gvcf_files}
         sample_map    : ${params.sample_map}
         """
         .stripIndent()



/*
========================================================================================
   Workflow block:
========================================================================================
*/

workflow {
    file(params.outdir).mkdirs()  // This line needs `params.outdir` to be initialized
    file(params.sample_map).parent.mkdirs() // Ensure the output directories for sample map exist
    file(params.joint_outdir).mkdirs()  // Ensure output directories exist for "per chrm" joint calls

    // Check if dictionary file exists
    if (!file(params.dict).exists()) {
        println "Dictionary file does not exist. Creating it now."
        createDictionary(genome_ch)
    } else {
        println "Dictionary file already exists: ${params.dict}. Skipping creation."
    }

    // Check if BWA index files exist
    if (!file(params.bwt).exists()) {
        println "BWA index files do not exist. Creating them now."
        BWA_INDEX(genome_ch)
    } else {
        println "BWA index files already exist. Skipping indexing."
    }
    
    // Check if .fai file exists
    if (!file("${params.genome}.fai").exists()) {
        println ".fai index file does not exist. Creating it now."
        createFAI(genome_ch)
    } else {
        println ".fai index file already exists: ${params.genome}.fai. Skipping creation."
    }

    //Create TBI index files for each GVCF file (.csi did not work for me for combining gcvf for joint calls)
    createTBI(gvcf_files_ch_filtered)

    // Check if DB files exist
    if (!file(params.db_path).exists()) {
        println "GenomicsDBImport: Creating genomic database now."
        GenomicsDBImport(genome_ch, bed_file_ch, gvcf_ch)
    } else {
        println "GenomicsDBImport: Database already exists at ${params.db_path}. Skipping creation."
    }


    // Always create a new sample mapping file (as samples will increase with time)
    println "Creating a new sample mapping file."
    createSampleMap(gvcf_files_ch)


    // Append new samples to the genomic DB using the created sample map
    appendSamplesToGenomicsDB(params.sample_map, params.tmp_dir, params.db_path)

    // Execute the joint variant calling process
    chr_ch
    .view { "Emitting chromosome: $it" }
    .map { chr -> tuple(chr, file(params.genome), file(params.db_path)) }
    | JointCallsPerChromosome                                         // Pass tuples to the process


   // Merge VCFs in a single process
    MergeVcfs(chr_vcfs_ch)
}



/*
========================================================================================
   Processes blocks:
========================================================================================
*/


/*
 * Process to create dictionary file
 */


process createDictionary {
    tag{"CREATE_DICT ${referenceFile}"}
    label 'process_low'
    
    container 'broadinstitute/gatk:4.5.0.0'

    input:
    path referenceFile 

    output:
    path "${referenceFile.baseName}.dict"

    script:
    """
    gatk CreateSequenceDictionary -R ${referenceFile} -O ${referenceFile.baseName}.dict
    cp ${referenceFile.baseName}.dict ${projectDir}/reference_data/
    """
}


/*
 * Process to index the reference genome using bwa
 */


process BWA_INDEX {
    tag{"BWA_INDEX ${referenceFile}"}
    label 'process_high'

    container 'biocontainers/bwa:v0.7.17_cv1'

    input:
    path referenceFile

    output:
    tuple path( referenceFile ), path( "*" ), emit: bwa_index

    publishDir "${projectDir}/reference_data", mode: 'copy'

    script:
    """
    echo "Running BWA index"
    bwa index ${referenceFile}
    """
}

/*
* Process to create FASTA index file (.fai) for the reference genome hg38.fa. 
*/
process createFAI {
    tag "CREATE_FAI ${referenceFile}"
    label 'process_low'

     container 'biocontainers/samtools:v1.9-4-deb_cv1'  // Use the samtools container

    input:
    path referenceFile

    output:
    path "${referenceFile}.fai"  // Outputs the .fai index file

    publishDir "${projectDir}/reference_data", mode: 'copy'  // Store output in the reference_data directory

    script:
    """
    samtools faidx ${referenceFile}
    """
}

/*
 * Create TBI index files for each GVCF file (p.s.: I had .csi, which was not suitable for combineGVCFs/DBImport, we need .tbi).
 */


process createTBI {
    tag "CREATE_TBI ${gvcf_file}"  // Tags each task with the file name for easier tracking
    label 'process_high' // set to high as I find this to be very intensive 

    container 'broadinstitute/gatk:4.5.0.0'

    publishDir("${projectDir}/data/gvcf_samples", mode: 'copy')

    maxForks 3  // Limits the number of concurrent processes to 2-3, as the memory when over the head and this process failed on my machine

    input:
    path gvcf_file  // Input is one GVCF file at a time, passed via the channel

    output:
    path "${gvcf_file}.tbi"  // Output is the .tbi file created by GATK

    script:
    """
    echo "Indexing file: ${gvcf_file}"
    gatk --java-options "-Xmx8g" IndexFeatureFile -I ${gvcf_file}

    # message to indicate completion
    if [ -f "${gvcf_file}.tbi" ]; then
        echo "Successfully created TBI file for ${gvcf_file}."
    else
        echo "Failed to create TBI file for ${gvcf_file}."
    fi
    """
}


/*
 * Process for creating initial genomic DB with GenomicsDBImport
 */


process GenomicsDBImport {
    tag "GenomicsDBImport ${gvcfFile}"
    label 'process_high'

    container 'broadinstitute/gatk:4.5.0.0'
    
    publishDir "${params.db_path}", mode: 'copy'

    input:
    	path genomeFile
    	path bedFile
    	tuple path(gvcfFile), path(indexFile)

    script:
    """
    echo "Creating Genomics data store for ${gvcfFile}"
    gatk --java-options "-Xmx80g" GenomicsDBImport \
        --genomicsdb-workspace-path ${params.db_path} \
        --tmp-dir ${params.tmp_dir} \
        --merge-input-intervals \
        -L ${bedFile} \
        -V ${gvcfFile} \
        --reader-threads 4
    """
}


/*
 * Process for create sample mapping file
 */


process createSampleMap {
    tag{"CREATE_SAMPLE_MAP"}
    label 'process_low'

    publishDir("${params.outdir}/sample_mapping_file", mode: 'copy')

    input:
    path gvcf_files

    output:
    path "sample_map.txt"

    script:
    """
    echo "Creating sample map"
    for gvcf in \$(ls *.g.vcf.gz | sort); do
        basename=\$(basename \${gvcf} .merged.matefixed.sorted.markeddups.recal.g.vcf.gz)
        full_path=\$(realpath \${gvcf})
        echo "\${basename}\t \${full_path}" >> temp_sample_map.txt
    done
    
    # Remove the line associated with the sample used for initial DB creation
    sed -i '' '/1291919.merged.matefixed.sorted.markeddups.recal.g.vcf.gz/d' temp_sample_map.txt
    

    # Ensure only one file is generated (recall -  had some issues )
    mv temp_sample_map.txt sample_map.txt
    """
}




/*
 *  Process to append new samples to the existing GenomicsDB using the sample map
 */

process appendSamplesToGenomicsDB {
    tag "AppendSamples ${sample_map}"
    label 'process_high'

    container params.docker_image

    // Use the memory and CPU parameters
    //memory params.memory
    //cpus params.cpus

    publishDir "${params.db_path}", mode: 'copy'

    input:
        path sample_map 
        path tmp_dir 
        path db_path 

    script:
    """
    echo "Appending new samples to GenomicsDB from ${sample_map}"
    gatk --java-options "-Xmx80g" GenomicsDBImport \
        --genomicsdb-update-workspace-path ${db_path} \
        --sample-name-map ${sample_map} \
        --tmp-dir ${tmp_dir} \
        --batch-size 10 \
        --reader-threads 4
    """
}

process JointCallsPerChromosome {
    tag "Joint Calls for Chromosome ${chromosome}"
    label 'process_high'
 
    container params.docker_image

    publishDir "${params.joint_outdir}", mode: 'copy'

    input:
    tuple val(chromosome), path(genome_file), path(db_path)  // Accept all inputs as a tuple
    //val chromosome  // Chromosome number or label
    //path genome_file  // Reference genome
    //path db_path  // Path to GenomicsDB for joint calling

    output:
    path "${chromosome}_joint.vcf.gz"  // The output VCF for this chromosome
    path "${chromosome}_joint.vcf.gz.tbi"  // The VCF index file

    script:
    """
    gatk --java-options "-Xmx2g -XX:ParallelGCThreads=2" GenotypeGVCFs \
        -R /Users/justinazurauskiene/Desktop/vc_docker/reference_data/hg38.fa  \    // this should be changed  
        -V gendb://${db_path} \
        -L chr${chromosome} \
        -O ${chromosome}_joint.vcf.gz
    """
    
}

process MergeVcfs {
    tag {"Generate VCF list and merge"}
    label 'process_low'

    container 'broadinstitute/picard:latest'

    input:
    path vcfs

    output:
    path "merged.vcf.gz"

    script:
    """
    java -jar /usr/picard/picard.jar GatherVcfs \
        I=1_joint.vcf.gz \
        I=2_joint.vcf.gz \
        I=3_joint.vcf.gz \
        I=4_joint.vcf.gz \
        I=5_joint.vcf.gz \
        I=6_joint.vcf.gz \
        I=7_joint.vcf.gz \
        I=8_joint.vcf.gz \
        I=9_joint.vcf.gz \
        I=10_joint.vcf.gz \
        I=11_joint.vcf.gz \
        I=12_joint.vcf.gz \
        I=13_joint.vcf.gz \
        I=14_joint.vcf.gz \
        I=15_joint.vcf.gz \
        I=16_joint.vcf.gz \
        I=17_joint.vcf.gz \
        I=18_joint.vcf.gz \
        I=19_joint.vcf.gz \
        I=20_joint.vcf.gz \
        I=21_joint.vcf.gz \
        I=22_joint.vcf.gz \
        I=X_joint.vcf.gz \
        I=Y_joint.vcf.gz \
        O=merged.vcf.gz
    """
}
