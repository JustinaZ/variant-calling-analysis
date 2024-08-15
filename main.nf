/*
========================================================================================
   ADRC DMS core Variant Calling Workflow with Nextflow, University of Washington
========================================================================================
   Github   :
   Contact  : jzurausk@uw.edu
----------------------------------------------------------------------------------------
*/


nextflow.enable.dsl=2

// Pipeline input parameters
params.outdir = "${projectDir}/output"
params.genome = "${projectDir}/data/hg38.fa"
params.dict = "${params.outdir}/dict_files/hg38.dict"
params.sample_map = "${params.outdir}/sample_map.txt"
params.gvcf_files = "${projectDir}/data/*.g.vcf.gz" 
params.bwt = "${params.outdir}/bwa_index/hg38.fa.bwt" // for checking that at least one relevant file exist


println """\
         RUNNING ADRC DMS CORE V A R I A N T  C A L L I N G   P I P E L I N E
         ====================================================================
         genome        : ${params.genome}
         outdir        : ${params.outdir}
         gvcf_files    : ${params.gvcf_files}
         sample_map    : ${params.sample_map}
         """
         .stripIndent()

workflow {
    // Ensure the output directories for dictionary, index, the sample map, etc. exist
    file(params.dict).parent.mkdirs()
    file(params.bwt).parent.mkdirs()
    file(params.sample_map).parent.mkdirs()

    // Define the channel for the reference genome
    genome_ch = Channel.fromPath(params.genome, checkIfExists: true)

    // Define the channel for GVCF files
    gvcf_files_ch = Channel.fromPath(params.gvcf_files, checkIfExists: true).collect()

    // Check if dictionary file exists; if not, create it
    if (!file(params.dict).exists()) {
        println "Dictionary file does not exist. Creating it now."
        createDictionary(genome_ch)
    } else {
        println "Dictionary file already exists: ${params.dict}. Skipping creation."
    }

    // Check if BWA index file exists; if not, create it
    if (!file(params.bwt).exists()) {
        println "BWA index files do not exist. Creating them now."
        BWA_INDEX(genome_ch)
    } else {
        println "BWA index files already exist: ${params.bwt}. Skipping creation."
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
 * Create dictionary file.
 */
process createDictionary {
    tag{"CREATE_DICT ${referenceFile}"}
    label 'process_low'

    publishDir("${params.outdir}/dict_files", mode: 'copy')

    input:
    path referenceFile

    output:
    path "${referenceFile.baseName}.dict"

    script:
    """
    echo "Running GATK CreateSequenceDictionary"
    gatk CreateSequenceDictionary -R ${referenceFile} -O ${referenceFile.baseName}.dict
    ls -lh ${referenceFile.baseName}.dict
    """
}

/*
 * Index the reference genome for use by bwa and samtools.
 */
process BWA_INDEX {
  tag{"BWA_INDEX ${referenceFile}"}
  label 'process_high'

  publishDir("${params.outdir}/bwa_index", mode: 'copy')

  input:
  path referenceFile

  output:
  tuple path( referenceFile ), path( "*" ), emit: bwa_index

  script:
  """
  echo "Running BWA index"
  bwa index ${referenceFile}
  """
} 


/*
 * Create sample mapping file.
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
        echo "\${basename} \${full_path}" >> temp_sample_map.txt
    done
    
    # Ensure only one file is generated (recall -  had some issues )
    mv temp_sample_map.txt sample_map.txt
    """
}


