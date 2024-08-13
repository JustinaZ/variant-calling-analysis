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

println """\
         V A R I A N T  C A L L I N G   P I P E L I N E
         ===================================================
         genome       : ${params.genome}
         outdir       : ${params.outdir}
         """
         .stripIndent()

workflow {
    // Ensure the output directory exists
    file(params.dict).parent.mkdirs()

    // Define the channel for the reference genome
    genome_ch = Channel.fromPath(params.genome, checkIfExists: true)

    // Check if dictionary file exists; if not, create it
    if (!file(params.dict).exists()) {
        println "Dictionary file does not exist. Creating it now."
        createDictionary(genome_ch)
    } else {
        println "Dictionary file already exists: ${params.dict}. Skipping creation."
    }
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
