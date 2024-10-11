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
   Parameter block:
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
 // this is single sample only!! To be used for creating initial Genomic DB, all other samples will be appended


params.gvcf_files = "${projectDir}/data/gvcf_samples/*.g.vcf.gz" //  Path pattern to all GVCF files
params.sample_map = "${params.outdir}/sample_mapping_file/sample_map.txt" // Path for the sample map



/*
========================================================================================
   Channel definitions:
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

//gvcf_files_ch_flat = Channel.fromPath(params.gvcf_files, checkIfExists: true) // (for creating .tbi)
// Modifying the gvcf_files_ch_flat to check if the .tbi file exists
gvcf_files_ch_filtered = Channel
    .fromPath(params.gvcf_files, checkIfExists: true)
    .filter { gvcf_file -> !file("${gvcf_file}.tbi").exists() }  // Only include files without existing .tbi






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
    
    //Create TBI index files for each GVCF file (.csi does not work for me for combining gcvf for joint calls)
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
 * Process for uploading new samples to pre-existing DB
 */


process appendSamplesToGenomicsDB {
    tag "APPEND_SAMPLES_TO_GENOMICS_DB"
    label 'process_high'

    container 'broadinstitute/gatk:4.5.0.0'
    
    publishDir "${params.db_path}", mode: 'copy'

    input:
    path genomeFile
    path bedFile
    tuple path(gvcfFile), path(indexFile)  // Use the newly created GVCF and its index

    script:
    """
    echo "Appending sample: ${gvcfFile} to existing GenomicsDB"
    gatk --java-options "-Xmx80g" GenomicsDBImport \
        --genomicsdb-update-workspace-path ${params.db_path} \
        --tmp-dir ${params.tmp_dir} \
        -L ${bedFile} \
        -V ${gvcfFile} \
        --reader-threads 4
    """
}
