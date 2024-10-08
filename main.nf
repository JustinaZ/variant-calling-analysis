/*
========================================================================================
   ADRC DMS core Variant Calling Workflow with Nextflow, University of Washington
========================================================================================
   Github   :
   Contact  : jzurausk@uw.edu
----------------------------------------------------------------------------------------
*/

params.outdir = "${projectDir}/output"  // output directory
params.genome = "${projectDir}/reference_data/hg38.fa" // location where .fa is kept
params.dict = "${projectDir}/reference_data/hg38.dict" // location where .dict file will be created
params.bwt = "${projectDir}/reference_data/hg38.fa.bwt" // path to check if BWA index exists

params.db_path = "${projectDir}/output/genomics_db" // dir where genomic db to be created
params.tmp_dir = "${projectDir}/output/tmp" // dir where tmp files should be placed while creating genomic db
params.bed_file = "${projectDir}/data/bed_file/twist_plus_refseq_10-31-2018.grc38.bed" // location where the bed file should be kept

params.sample_gvcfs = "${projectDir}/data/gvcf_samples/1291919.merged.matefixed.sorted.markeddups.recal.g.vcf.gz"
 // single sample to be used for creating initial Genomic DB


genome_ch = Channel.fromPath(params.genome, checkIfExists: true)
dict_ch = Channel.fromPath(params.dict, checkIfExists: true)
bed_file_ch = Channel.fromPath(params.bed_file, checkIfExists: true)


// input channel definition for gvcfs and their indices
gvcf_ch = Channel
    .from(params.sample_gvcfs)
    .map { gvcfPath ->
        def indexPath = "${gvcfPath}.tbi"  // This assumes the index has the same base name as the GVCF sample
        [file(gvcfPath), file(indexPath)]   // Creates a tuple of [GVCF file, index file] (pairing together sample and it's index)
    } // This channel should pass one GVCF and its corresponding index to the process, which will use that to create the genomic DB.




workflow {
    file(params.outdir).mkdirs()  // This line needs `params.outdir` to be initialized

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

    // Check if DB files exist
    if (!file(params.db_path).exists()) {
        println "GenomicsDBImport: Creating genomic database now."
        GenomicsDBImport(genome_ch, bed_file_ch, gvcf_ch)
    } else {
        println "GenomicsDBImport: Database already exists at ${params.db_path}. Skipping creation."
    }
}





// createDictionary process
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

// BWA_INDEX process
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

// GenomicsDBImport process
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

