/*
========================================================================================
   ADRC DMS core Variant Calling Workflow with Nextflow, University of Washington
========================================================================================
   Github   : https://github.com/JustinaZ/variant-calling-analysis
   Contact  : jzurausk@uw.edu
----------------------------------------------------------------------------------------
*/


// Make sure DSL2 is on - just in case
nextflow.enable.dsl = 2

/*
========================================
========================================
   P A R A M E T E R  block:
========================================
========================================
*/

//======================================================================================================
//================= 1. CORE OUTPUT & REFERENCE SETTINGS/PARAMETERS: ====================================
//======================================================================================================


params.outdir       = "${projectDir}/output"                       // output directory
params.genome       = "${projectDir}/reference_data/hg38.fa"       // location where .fa is kept
params.dict         = "${projectDir}/reference_data/hg38.dict"       // location where .dict file will be created
params.bwt          = "${projectDir}/reference_data/hg38.fa.bwt"     // path to check if BWA index exists
params.gvcf_files   = "${projectDir}/data/gvcf_samples/*.g.vcf.gz"   // path pattern to all GVCF files


//======================================================================================================
//================= 2. GENOMICSDB IMPORT / BUILD SETTINGS: =============================================
//======================================================================================================


params.db_path      = "${projectDir}/output/genomics_db"           // dir where genomic db to be created
params.tmp_dir      = "${projectDir}/output/tmp"                   // temporary directory for creating the genomic db
params.bed_file     = "${projectDir}/data/bed_file/twist_plus_refseq_10-31-2018.grc38.bed"  // bed file

params.sample_gvcfs = "${projectDir}/data/gvcf_samples/1291919.merged.matefixed.sorted.markeddups.recal.g.vcf.gz" // this is !!!only!!! a single sample given; =>
// To be used for creating initial Genomic DB; should be changed/coded-up as a user input type in the future & in case a new DB creation is needed
params.sample_map   = "${params.outdir}/sample_mapping_file/sample_map.txt"  // sample mapping file path


//======================================================================================================
//================= 3. JOINT GENOTYPING / VARIANT-CALLING SETTINGS: ====================================
//======================================================================================================


params.chromosomes  = (1..22).collect { it.toString() } + ['X','Y']  // list of chromosome IDs to process
params.joint_outdir = "${params.outdir}/joint_vcfs"                  // directory where per-chromosome & merged VCFs will be saved


//======================================================================================================
//================= 4. VARIANT RECALIBRATION (VQSR) SETTINGS: ==========================================
//======================================================================================================

// SNP public resources (previously downloaded files)
params.hapmap_resource = "${projectDir}/gatk_resources/hapmap_3.3.hg38.vcf.gz"
params.omni_resource   = "${projectDir}/gatk_resources/1000G_omni2.5.hg38.vcf.gz"
params.onekg_resource  = "${projectDir}/gatk_resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz"

// Shared SNP+Indel public resources (previously downloaded files)
params.dbsnp_resource  = "${projectDir}/gatk_resources/Homo_sapiens_assembly38.dbsnp138.vcf.gz"

// Indel public resource (previously downloaded files)
params.mills_resource  = "${projectDir}/gatk_resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

// VQSR outputs:
params.vqsr             = "${params.outdir}/vqsr"

params.snp_recal_file   = "${params.vqsr}/snp/snp_recalibration.recal"
params.snp_tranches_file= "${params.vqsr}/snp/snp_recalibration.tranches"
params.snp_plots_file   = "${params.vqsr}/snp/snp_recalibration.plots.R"

params.indel_recal_file    = "${params.vqsr}/indel/indel_recalibration.recal"
params.indel_tranches_file = "${params.vqsr}/indel/indel_recalibration.tranches"
params.indel_plots_file    = "${params.vqsr}/indel/indel_recalibration.plots.R"


//======================================================================================================
//================= 4. GENE SUBSETTING SETTINGS: =======================================================
//======================================================================================================

// Path to .bed files for subsetting AD related genes (R code to create it is within project directory)
params.genes_bed = "${projectDir}/data/bed_file/genes_of_interest.bed"


//======================================================================================================
//================= 5. SnpEff SETTINGS: ================================================================
//======================================================================================================


// 1. The local folder containing the SnpEff genome data.
// 2. The name of the genome as recognised by the SnpEff default config
// "hg38" is/should be recognised by the built-in config;
// other can be used by changing this param.
// 3. Number of threads to let SnpEff use (for faster annotation).

params.snpeff_db_dir = "${projectDir}/snpEff_db_resource"
params.snpeff_genome = "hg38"
params.snpeff_threads = 4
params.out_prefix     = "annotated_vc_subset"

//params.snpeff_allow_download = false

/*
========================================
========================================
   C H A N N E L definition block:
========================================
========================================
*/

//------------------------------------------------------------------------------------------------------
//----------------------- 1. Reference & targets -------------------------------------------------------
//------------------------------------------------------------------------------------------------------
// Single files that downstream processes consume (FASTA, dict, BED).


genome_ch = Channel.fromPath(params.genome, checkIfExists: true)
dict_ch   = Channel.fromPath(params.dict, checkIfExists: true)
bed_file_ch = Channel.fromPath(params.bed_file, checkIfExists: true)


//------------------------------------------------------------------------------------------------------
//----------------------- 2. Seed GVCF for initial GenomicsDB ------------------------------------------
//------------------------------------------------------------------------------------------------------
// One tuple: [ seed.g.vcf.gz , seed.g.vcf.gz.tbi ]; GVCFs + .tbi index


gvcf_ch = Channel
    .from(params.sample_gvcfs)
    .map { gvcfPath ->
        def indexPath = "${gvcfPath}.tbi"
        [file(gvcfPath), file(indexPath)]
    }


//------------------------------------------------------------------------------------------------------
//----------------------- 3. All GVCFs (for sample-map creation, etc.) ---------------------------------
//------------------------------------------------------------------------------------------------------
// `collect()` groups all matches into a single emission so Nextflow stages the
// entire set into the process work dir (needed by CREATE_SAMPLE_MAP loop).
// i.e. Gather all GVCFs once - this is also for checking if .tbi exist


gvcf_files_ch = Channel.fromPath(params.gvcf_files, checkIfExists: true).collect()


//------------------------------------------------------------------------------------------------------
//----------------------- 4. Only GVCFs that *lack* a Tabix index --------------------------------------
//--------------------------------------------------------------------------------------------------------
// Feeds CREATE_TBI so we index just what is missing.
// (It looks specifically for `.tbi`; files indexed with CSI won't match!!!)
// If this yields zero items, CREATE_TBI would not run


gvcf_files_ch_filtered = Channel
    .fromPath(params.gvcf_files, checkIfExists: true)
    .filter { gvcf_file -> !file("${gvcf_file}.tbi").exists() }


//--------------------------------------------------------------------------------------------------------
//---------------------- 5. GenomicsDB workspace ---------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
// Carries the directory where the GenomicsDB will be created/updated.
// It is primary: pass as a value so it works for create-or-update flows.


genomics_db_ch = Channel.value(params.db_path)


//--------------------------------------------------------------------------------------------------------
//---------------------- 6. Chromosome driver list for joint genotyping ----------------------------------
//--------------------------------------------------------------------------------------------------------


fai_ch = Channel.fromPath("${params.genome}.fai", checkIfExists: true)
chr_ch = Channel.fromList(params.chromosomes)
chr_vcfs_ch = Channel.fromPath("${projectDir}/output/joint_vcfs/*.vcf.gz").collect()


//--------------------------------------------------------------------------------------------------------
//---------------------- 7. Variant recalibration ----------------------------------
//--------------------------------------------------------------------------------------------------------


// Merged VCF channel
merged_vcf_ch = Channel
    .fromPath("${projectDir}/output/joint_vcfs/merged.vcf.gz", checkIfExists: true)
    .map { vcfPath -> [ file(vcfPath), file("${vcfPath}.tbi") ] }


// SNP resource VCFs
hapmap_ch = Channel.fromPath(params.hapmap_resource)
    .map { vcfPath -> [file(vcfPath), file("${vcfPath}.tbi")] }

omni_ch = Channel
    .fromPath(params.omni_resource)
    .map { vcfPath -> [file(vcfPath), file("${vcfPath}.tbi")] }
onekg_ch = Channel
    .fromPath(params.onekg_resource)
    .map { vcfPath -> [file(vcfPath), file("${vcfPath}.tbi")] }



// Indel resource VCFs
dbsnp_ch = Channel
    .fromPath(params.dbsnp_resource)
    .map { vcfPath -> [file(vcfPath), file("${vcfPath}.tbi")] }
mills_ch = Channel
    .fromPath(params.mills_resource)
    .map { vcfPath -> [file(vcfPath), file("${vcfPath}.tbi")] }


// Reference genome + .fai + .dict
fafai_ch = Channel
    .from(params.genome)
    .map { faPath ->
        def faiPath = "${faPath}.fai"
        def dictPath = faPath.replace(".fa", ".dict")
        [file(faPath), file(faiPath), file(dictPath)]
    }


// For SNP recall
snp_recal_file_ch = Channel.fromPath(params.snp_recal_file)
    .map { recalPath -> [ file(recalPath), file("${recalPath}.idx") ] }

snp_tranches_file_ch = Channel.fromPath(params.snp_tranches_file)



//--------------------------------------------------------------------------------------------------------
//---------------------- 8. Subsetting genes  ------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------


// Channel for subsetting vqsr final result to focus on AD40 gene set
genes_of_interest_ch = Channel.fromPath(params.genes_bed)



//--------------------------------------------------------------------------------------------------------
//---------------------- 9. SnpEff database  -------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------

// SnpEff database dir must contain a subfolder named like params.snpeff_genome, e.g. hg38

snpeff_db_dir_ch = Channel.value( file(params.snpeff_db_dir) )

snpeff_genome_ch = Channel.value( params.snpeff_genome )



/*
========================================
========================================
   F U N C T I O N  block:
========================================
========================================
*/

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

// should works for local and most provider-backed paths (use Files.*, not new File)
def needsIndex = { Path vcf ->
    Path tbi = Paths.get(vcf.toString() + '.tbi')
    if (!Files.exists(tbi)) return true
    return Files.getLastModifiedTime(vcf).toMillis() > Files.getLastModifiedTime(tbi).toMillis()
}





def _now = new Date().format("yyyy-MM-dd HH:mm:ss z")


println """\
         =========================================================================================
         RUNNING ADRC DMS • R A R E  V A R I A N T  I D E N T I F I C A T I O N  P I P E L I N E •
         -----------------------------------------------------------------------------------------
         Launch time   : ${_now}
         Project dir   : ${projectDir}
         Genome FASTA  : ${params.genome}
         Output dir    : ${params.outdir}
         GVCF glob     : ${params.gvcf_files}
         Sample map    : ${params.sample_map}
         GenomicsDB    : ${params.db_path}
         BED targets   : ${params.bed_file}
         =========================================================================================
         """.stripIndent()

/*
========================================================================================
   Workflow block:

   - Creates required reference indices if missing
   - Ensures sample map and GVCF tabix indexes exist
   - Creates (or reuses) the GenomicsDB workspace
   - Performs joint-calls per chromosome 
   - Merges per-chromosome VCFs and creates a Tabix index
   - Exposes the merged VCF via `merged_vcf_ch` for downstream VQSR
========================================================================================
*/



workflow {

//-------------------------------------------------------------------------------
// 0) Setup output directories: create required folders for pipeline outputs
//-------------------------------------------------------------------------------


    file(params.outdir).mkdirs()
    file(params.sample_map).parent.mkdirs()

    // Create tmp dir only if we’re also creating the DB
    if (!file(params.db_path).exists()) {
    	file(params.tmp_dir).mkdirs()
    }

  
//-------------------------------------------------------------------------------
// 1) Reference prep: .dict (GATK), BWA index, and .fai (samtools)
//-------------------------------------------------------------------------------


  // Create dictionary if needed
    if (!file(params.dict).exists()) {
        println "\nSTEP DICT: Dictionary  file (.dict) is missing → Creating it now."
        createDictionary(genome_ch)
    } else {
        println "\nSTEP DICT: Dictionary file already exists: ${params.dict} → Skipping creation."
    }
    
   
  // BWA indexing with more comprehensive missing file checking added upfront - needs testing:
    def bwaIndexFiles = [
    	"${params.genome}.amb",
    	"${params.genome}.ann",
    	"${params.genome}.bwt",
    	"${params.genome}.pac",
    	"${params.genome}.sa"
    ]

    def missingIndexFiles = bwaIndexFiles.findAll { !file(it).exists() }

    if (missingIndexFiles) {
        println "\nSTEP BWA: Missing BWA index files: ${missingIndexFiles.join(', ')} → Creating them now."
        BWA_INDEX(genome_ch)
    } else {
        println "\nSTEP BWA: All BWA index files present → Skipping creation."
    }


  // .fai

    if (!file("${params.genome}.fai").exists()) {
        println "\nSTEP FAI: The .fai file is missing → Creating it now."
        createFAI(genome_ch)
    } else {
        println "\nSTEP FAI: The .fai  file already exist ${params.genome}.fai → Skipping creation."
    }

//-------------------------------------------------------------------------------
// 2) GVCF preparation: create Tabix indexes only for files that are missing .tbi
//-------------------------------------------------------------------------------


   missing_or_stale_ch = gvcf_files_ch.flatMap { List gvcfs ->
    if (!gvcfs || gvcfs.isEmpty()) {
        println "[tbi] WARNING: no GVCF files matched: ${params.gvcf_files}"
        return []
    }
    def todo = gvcfs.findAll { f -> needsIndex(f) }
    if (todo) {
        println "\nSTEP TBI: ${todo.size()} Tabix indexes missing or stale — creating now:"
        println todo.collect { "   • ${it}" }.join('\n')
        return todo
    } else {
        println "\nSTEP TBI: All Tabix indexes present and up-to-date → Skipping indexing."
        return []
    }
   }

   // Launch indexing for those (0 tasks if none)
   createTBI(missing_or_stale_ch)

   // Pre-existing, up-to-date pairs
   have_tbi_ch = gvcf_files_ch.flatMap { List gvcfs ->
    gvcfs.findAll { f -> !needsIndex(f) }
   }.map { f -> tuple(file(f), file("${f}.tbi")) }


   indexed_gvcf_ch = have_tbi_ch.mix(
   createTBI.out.map { tbi ->
    def vcf = tbi.toString().replaceFirst(/\.tbi$/, '')
    tuple(file(vcf), file(tbi))
   }
   )


//-------------------------------------------------------------------------------
// 3) GenomicsDB: create if absent; otherwise reuse existing workspace
//-------------------------------------------------------------------------------
// p.s. run time of this was about 10mins to created db


    if (!file(params.db_path).exists()) {
        println "\nSTEP DB CREATE: Creating new workspace -- ${params.db_path}."
        GenomicsDBImport(genome_ch, bed_file_ch, gvcf_ch)
    } else {
        println "\nSTEP DB CREATE:  Workspace exists at ${params.db_path}  →  Skipping creation."
    }


//----------------------------------------------------------------------------
// 4) Sample map: rebuild each run to reflect current GVCF set (excludes seed sample)
//----------------------------------------------------------------------------


   if (!file(params.sample_map).exists()) {
    println "STEP SAMPLE_MAP: sample_map.txt not found → Creating a new one at ${params.sample_map}"
    sample_map_ch = createSampleMap(gvcf_files_ch)
   } else {
    println "STEP SAMPLE_MAP: sample_map.txt already exists → Rebuilding to reflect current GVCF set"
    sample_map_ch = createSampleMap(gvcf_files_ch)
   }


//-----------------------------------------------------------------------------
// 5) GenomicsDB workspace
//-----------------------------------------------------------------------------


// For now if you are just re-running with current samples (Darwas) the below 
// code line should be commented out to avoid crashing pipeline.
// Uncomment in case new samples/maps is present and needs uploading.
// In the future, implement if/else check based on new samples present and
// which are in current DB, so it is automatic. R script for pulling out samples
// from current DB is included - just needs incorporating in here for 
// cross- referencing.
// This can be prone to "hanging",i.e. adjust config file or process params 
// The run time for all 32 samples is ~  5h 42m 4s


// ------------- append everything (all samples) in that map -------------
//   appendSamplesToGenomicsDB(
//       sample_map_ch,                 // the path emitted by createSampleMap
//       Channel.value(params.tmp_dir), // tmp directory as a channel
//       Channel.value(params.db_path)  // existing GenomicsDB workspace
//   )


//-----------------------------------------------------------------------------
// 6) Joint variant calling per chromosome
//-----------------------------------------------------------------------------


// Commented out -  if not running new joint calls
// The run time for all 32 samples is ~ 4h 4m 59s

// chr_ch
//   .view { "Emitting chromosome: $it" }
//   .map { chr -> tuple(
//       chr,
//       file(params.genome),
//       file("${params.genome}.fai"),
//       file(params.dict),                 
//       file(params.db_path)
//   )}
//   | JointCallsPerChromosome


// Merge VCFs from per-chromosome calls
   merged_vcf_result = MergeVcfs(chr_vcfs_ch)
 

// Index the merged files if needed
if (!file("${projectDir}/output/joint_vcfs/merged.vcf.gz.tbi").exists()) {
    println "\nSTEP MERGED INDEX: Creating index for merged VCF."
    CreateIndex(merged_vcf_result)
} else {
    println "\nSTEP MERGED INDEX: Index for merged VCF exists → Skipping index creation."
}



//-----------------------------------------------------------------------------
// 7) Variant recalibration
//-----------------------------------------------------------------------------


// Run VQSRSNP for SNP 
// The run time for all 32 samples is ~ 14m 11s

    if (!file(params.snp_recal_file).exists()) {
        println "SNP recal file not found  → Running VQSRSNP."
        VQSRSNP_OUT = VQSRSNP(
            fafai_ch,            // [ fa, fa.fai, fa.dict ]
            merged_vcf_ch,       // [ merged.vcf.gz, merged.vcf.gz.tbi ]
            hapmap_ch,
            omni_ch,
            onekg_ch
        )
    } else {
        println "SNP recal file exists  → skipping VQSRSNP."
        VQSRSNP_OUT = Channel.value([
            file(params.snp_recal_file),
            file(params.snp_tranches_file),
            file(params.snp_plots_file)
        ])
    }


// Run VQSRIndel for INDEL 
// The run time for all 32 samples is ~  12m 25s

    if (!file(params.indel_recal_file).exists()) {
        println "Indel recal file not found  → Running VQSRIndel."
        VQSRINDEL_OUT = VQSRIndel(
            fafai_ch,
            merged_vcf_ch,
            mills_ch,
            dbsnp_ch
        )
    } else {
        println "Indel recal file exists  → skipping VQSRIndel."
        VQSRINDEL_OUT = Channel.value([
            file(params.indel_recal_file),
            file(params.indel_tranches_file),
            file(params.indel_plots_file)
        ])
    }



//Index the newly created .recal files for SNPs & Indels

  def snp_recal_file = file(params.snp_recal_file)
  def indel_recal_file = file(params.indel_recal_file)

  if (!file("${snp_recal_file}.idx").exists()) {
    println "STEP VQSR SNP INDEX: Missing index for ${snp_recal_file} → creating now."
    IndexSNPRecalFile(Channel.of(snp_recal_file))
  } else {
    println "STEP VQSR SNP INDEX: Index already exists for ${snp_recal_file} → skipping."
  }

  if (!file("${indel_recal_file}.idx").exists()) {
    println "STEP VQSR INDEL INDEX: Missing index for ${indel_recal_file} → creating now."
    IndexIndelRecalFile(Channel.of(indel_recal_file))
  } else {
    println "STEP VQSR INDEL INDEX: Index already exists for ${indel_recal_file} → skipping."
  }


//-----------------------------------------------------------------------------
// 7) Apply Variant recalibration
//-----------------------------------------------------------------------------
// The run time for this section is ~ 14m 11s2m 12s

//Build channels that contain [recalFile, recalFile.idx], plus the .tranches path
    
    snp_recal_file_ch = VQSRSNP_OUT.map { tuple( file(it[0]), file("${it[0]}.idx") ) }
    snp_tranches_file_ch = VQSRSNP_OUT.map { file(it[1]) }

    indel_recal_file_ch = VQSRINDEL_OUT.map { tuple( file(it[0]), file("${it[0]}.idx") ) }
    indel_tranches_file_ch = VQSRINDEL_OUT.map { file(it[1]) }


//Apply SNP recalibration to the merged VCF file, i.e. to produced variant calls

    snp_recalibrated_vcf = ApplyVQSRSNP(
        merged_vcf_ch,
        snp_recal_file_ch,
        snp_tranches_file_ch
    )


//Apply Indel recalibration to the SNP-recalibrated VCF
 
    indel_recalibrated_vcf = ApplyVQSRIndel(
        snp_recalibrated_vcf,
        indel_recal_file_ch,
        indel_tranches_file_ch
    )

//-----------------------------------------------------------------------------
// 8) Subset variant recalibrated VCF file
//-----------------------------------------------------------------------------


// We select high-confidence (i.e. PASS) variants and variants under gene list of interest

    subset_vcf_result = SubsetPassVariants(
        indel_recalibrated_vcf,
        genes_of_interest_ch
    )


//-----------------------------------------------------------------------------
// 9) Annotate with SnpEff
//-----------------------------------------------------------------------------

// Output will be a new file with annotations, ${params.out_prefix}.snpeff.vcf.gz, plus its .tbi.

    annotated_vcf = SnpEffAnnotation(subset_vcf_result)
    

//annotated_plain = SnpEffAnnotation(
//       subset_vcf_result,		// tuple (VCF.gz, VCF.gz.tbi)
//       snpeff_db_dir_ch,                // dir containing `data/<genome>/...`
//       snpeff_genome_ch 		// e.g., "hg38" 
// )

// then compress+index in case SnpEffAnnotation cannot handle anymore

//annotated_vcf = BgzipTabix(annotated_plain)

}


/*
========================================
========================================
   P R O C E S S definition block:
========================================
========================================
*/


//--------------------------------------------------------
// ---------- PROCESS TO CREATE DICTIONARY FILE ----------
//--------------------------------------------------------

 
/* INFO:
 * Build a .dict dictionary for the reference FASTA
 * (required by BWA and GATK). 
 *
 * Inputs
 *  • reference : FASTA file to index (hg38.fa)
 * Outputs
 *  • <ref>.dict :GATK sequence dictionary
 *
 */


process createDictionary {
    tag { "DICT_${referenceFile.baseName}" }   // should be seen as DICT_hg38
    label 'process_low'
    
    container 'broadinstitute/gatk:4.5.0.0'

    input:
      path referenceFile 

    output:
      path "${referenceFile.baseName}.dict"

    script:
    """
    set -euo pipefail
    echo "[\$(date)] Creating sequence dictionary for ${referenceFile} …"

    gatk CreateSequenceDictionary -R ${referenceFile} -O ${referenceFile.baseName}.dict
    cp ${referenceFile.baseName}.dict ${projectDir}/reference_data/
    """
}


//-----------------------------------------------------------------
// ---------- PROCESS TO INDEX REFERENCE GENOME WITH BWA ----------
//-----------------------------------------------------------------


/* INFO:
 * Create the BWA index files: (`.amb .ann .bwt .pac .sa`) for the
 * reference FASTA.  Needed by any downstream read-mapping or
 * GenotypeGVCFs steps that depend on BAM headers.
 *
 * Inputs
 *  • reference : FASTA file (e.g. hg38.fa)
 * Outputs
 *  • reference index tuple : [ reference FASTA , associated *.{amb,ann,bwt,pac,sa} ]
 *
 * NOTE: not consumed in the current workflow, but kept for
 *       future pipeline extensions that may need BWA-MEM.
 */


process BWA_INDEX {
    tag { "BWAIDX_${referenceFile.baseName}" }    // should be seen as BWAIDX_hg38
    label 'process_high'

    container 'biocontainers/bwa:v0.7.17_cv1'

    input:
      path referenceFile

    output:
      tuple path(referenceFile), path("*"), emit: bwa_index

    publishDir "${projectDir}/reference_data", mode: 'copy'

    script:
    """
    set -euo pipefail
    echo "[\$(date)] Building BWA index for ${referenceFile} …"
    bwa index ${referenceFile}
    """
}


//-----------------------------------------------------------------
// ---------- PROCESS TO CREATE FASTA INDEX FILE ------------------
//-----------------------------------------------------------------


/* INFO:
 * Generate the .fai FASTA index using samtools.
 *
 * Inputs
 *  • reference : FASTA file (e.g. hg38.fa)
 * Outputs
 *  • reference.fai : tab-delimited index
 *
 */


process createFAI {
    tag { "FAI_${referenceFile.baseName}" }                   // should be seen as FAI_hg38
    label 'process_low'

    container 'biocontainers/samtools:v1.9-4-deb_cv1'

    input:
      path referenceFile

    output:
      path "${referenceFile}.fai"

    publishDir "${projectDir}/reference_data", mode: 'copy'

    script:
    """
    set -euo pipefail
    echo "[\$(date)] Generating FASTA .fai index for ${referenceFile} …"
    samtools faidx ${referenceFile}
    """
}


//-----------------------------------------------------------------
// ---------- PROCESS TO CREATE TBI INDEX FILES  ------------------
//-----------------------------------------------------------------


/* INFO:
 * Create a .tbi tab-index for each gz-compressed GVCF sample file (if missing/needed)
 *
 * Inputs
 *  • gvcf : *.g.vcf.gz file
 * Outputs
 *  • gvcf.tbi : Tabix index
 *
 */


process createTBI {
  tag { "TBI_${gvcf_file.baseName}" }
  label 'process_low'
  container 'broadinstitute/gatk:4.5.0.0'
  // publishDir("${projectDir}/data/gvcf_samples", mode: 'copy', overwrite: true)

  maxForks 3

  input:
    path gvcf_file

  output:
    tuple path(gvcf_file), path("${gvcf_file}.tbi"), emit: indexed

  script:
  """
  set -euo pipefail
  echo "[\$(date)] Starting TBI indexing for ${gvcf_file}"
  gatk --java-options "-Xmx8g" IndexFeatureFile -I ${gvcf_file}
  test -f "${gvcf_file}.tbi"
  """
}


//-----------------------------------------------------------------
// ---------- PROCESS TO CREATE INITIAL GENOMIC DB INDEX FILES  ---
//-----------------------------------------------------------------


/*
 * Create the initial genomic DB with GenomicsDBImport;
 * i.e. Initialise the GenomicsDB workspace from one (seed) GVCF file/sample
 *
 * Inputs
 *  • genomeFile : reference FASTA
 *  • bedFile    : target intervals (BED)
 *  • (gvcf , tbi) tuple : seed sample GVCF + its index
 *
 * Outputs
 *  • <workspace-dir>      : the GenomicsDB directory (emitted as a `path`)
 *
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
    gatk --java-options "-Xmx80g" GenomicsDBImport \
        --genomicsdb-workspace-path ${params.db_path} \
        --tmp-dir ${params.tmp_dir} \
        --merge-input-intervals \
        -L ${bedFile} \
        -V ${gvcfFile} \
        --reader-threads 4
    """
}


//-----------------------------------------------------------------
// ---------- PROCESS TO CREATE SAMPLE MAPPING FILE  --------------
//-----------------------------------------------------------------


/* INFO:
 * Build a two-column -- "sample-name → absolute-path" -- file required by
 * GenomicsDBImport (called sample-map format).
 *
 * Inputs
 *   • gvcf_files : all *.g.vcf.gz files staged in work dir
 *     (this should not include previous single sample used to build initial genomicDB)
 * Outputs
 *   • sample_map.txt : tab-delimited mapping
 */


process createSampleMap {
    tag { "SAMPLE_MAP" }
    label 'process_low'

    publishDir("${params.outdir}/sample_mapping_file", mode: 'copy')

    input:
      path gvcf_files

    output:
      path "sample_map.txt"

    script:
    """
    set -euo pipefail
    echo "[\$(date)] Creating sample mapping file …"

    for gvcf in \$(ls *.g.vcf.gz | sort); do
        basename=\$(basename \${gvcf} .merged.matefixed.sorted.markeddups.recal.g.vcf.gz)
        full_path=\$(realpath \${gvcf})
        echo "\${basename}\t\${full_path}" >> temp_sample_map.txt
    done
    sed -i.bak '/1291919.merged.matefixed.sorted.markeddups.recal.g.vcf.gz/d' temp_sample_map.txt
    rm temp_sample_map.txt.bak
    mv temp_sample_map.txt sample_map.txt
    """
}


//-----------------------------------------------------------------
// ----- PROCESS TO APPEND NEW SAMPLES TO EXISTING GENOMIC DB  ----
//-----------------------------------------------------------------



/* INFO:
 * Append new samples to an existing GenomicsDB workspace using a
 * GATK sample-name map (two columns: <sample_id>\t<abs_path_to_gvcf>).
 *
 * Inputs
 *   • sample_map : tab-delimited sample-name map (new samples only)
 *   • tmp_dir    : temporary/scratch directory for GenomicsDBImport
 *   • db_path    : existing GenomicsDB workspace to update
 * Outputs
 *   • (none)     : workspace at db_path is updated in place
 *
 */


process appendSamplesToGenomicsDB {
    tag "AppendSamples ${sample_map}"
    label 'process_high'
    container 'broadinstitute/gatk:4.5.0.0'

    publishDir "${params.db_path}", mode: 'copy'

    input:
      path sample_map
      path tmp_dir
      path db_path

    script:
    """
    gatk --java-options "-Xmx80g" GenomicsDBImport \
        --genomicsdb-update-workspace-path ${db_path} \
        --sample-name-map ${sample_map} \
        --tmp-dir ${tmp_dir} \
        --batch-size 10 \
        --reader-threads 4
    """
}


//-----------------------------------------------------------------
// ---------- PROCESS FOR JOINT CALLING PER CHROMOSOME ------------
//-----------------------------------------------------------------


/* INFO:
 * Genotype GVCFs from an existing GenomicsDB workspace for a single
 * chromosome and write a bgzipped VCF + tabix index.
 *
 * Inputs
 *   • chromosome : e.g. "1", "2", …, "X", "Y"
 *   • genome_file: reference FASTA (must have .fai)
 *   • db_path    : GenomicsDB workspace
 * Outputs
 *   • <chr>_joint.vcf.gz and .tbi
 */


process JointCallsPerChromosome {
    tag { "JOINT_${chromosome}" }
    label 'process_high'
 
    container 'broadinstitute/gatk:4.5.0.0'
    publishDir "${params.joint_outdir}", mode: 'copy'

    input:
        tuple val(chromosome),
          path(genome_file),
          path(genome_fai),
          path(genome_dict),  
          path(db_path)
      //tuple val(chromosome), path(genome_file), path(genome_fai), path(genome_dict), path(db_path)


    output:
      path "${chromosome}_joint.vcf.gz"
      path "${chromosome}_joint.vcf.gz.tbi"

    script:
    """
    set -euo pipefail
    echo "[\$(date)] Joint genotyping: chromosome ${chromosome}"

    gatk --java-options "-Xmx2g -XX:ParallelGCThreads=2" GenotypeGVCFs \
        -R ${genome_file} \
        -V gendb://${db_path} \
        -L chr${chromosome} \
        -O ${chromosome}_joint.vcf.gz
    """
}



//-----------------------------------------------------------------
// ----------------- PROCESS TO MERGE VCFS ------------------------
//-----------------------------------------------------------------


/*
 * Process to merge VCFs
 */
process MergeVcfs {
    tag { "Generate VCF list and merge" }
    label 'process_low'

    container 'broadinstitute/picard:latest'
    publishDir "${params.joint_outdir}", mode: 'copy'

    input:
      path vcfs

    output:
      path "merged.vcf.gz"

    script:
    """
    set -euo pipefail
    echo "[\$(date)] Merging VCFs with Picard GatherVcfs …"

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

    echo "[\$(date)] Merge complete: merged.vcf.gz"
    """
}


//-----------------------------------------------------------------
// -------------- PROCESS TO INDEX MERGED VCF ---------------------
//-----------------------------------------------------------------


/* INFO:
 * Create a Tabix index (.tbi) for the merged bgzipped VCF.
 *
 * Inputs
 *   • merged_vcf : merged.vcf.gz (bgzip-compressed VCF)
 * Outputs
 *   • merged_vcf.tbi : Tabix index corresponding to the input VCF
 */


process CreateIndex {
    tag { 'INDEX_MERGED' }
    label 'process_low'

    container 'broadinstitute/gatk:4.5.0.0'
    publishDir "${params.joint_outdir}", mode: 'copy'

    input:
      path merged_vcf

    output:
      path "${merged_vcf}.tbi"

    script:
    """
    set -euo pipefail
    echo "[\$(date)] Creating index for merged file: ${merged_vcf}"

    gatk --java-options "-Xmx8g" IndexFeatureFile -I ${merged_vcf}

    echo "[\$(date)] Done: ${merged_vcf}.tbi"
    """
}


//-----------------------------------------------------------------
// ---------- PROCESS TO RUN VARIANT RECALIBRATION (SNPs) ---------
//-----------------------------------------------------------------


/* INFO:
 * Perform Variant Quality Score Recalibration (VQSR) on SNPs using
 * GATK VariantRecalibrator. Builds a statistical model from multiple
 * training resources (HapMap, Omni, 1000G) and annotates SNPs with
 * recalibrated quality scores and tranches.
 *
 * Inputs
 *   • genome_file      : reference FASTA (with .fai and .dict)
 *   • merged_vcfs      : cohort-level joint-called VCF (bgzip-compressed)
 *   • hapmap_vcf       : HapMap training set resource (+ index)
 *   • omni_vcf         : Omni 2.5M training resource (+ index)
 *   • onekg_vcf        : 1000G Phase 1 SNP training resource (+ index)
 *
 * Outputs
 *   • snp_recalibration.recal    : recalibration table
 *   • snp_recalibration.tranches : tranches file for filtering thresholds
 *   • snp_recalibration.plots.R  : R script for generating diagnostic plots
 *
 * Notes
 *   - This step is SNP-specific (uses `-mode SNP`).
 *   - Should be followed by ApplyVQSR to apply the recalibration model.
 */


process VQSRSNP {
    tag "Variant Recalibration for SNPs"
    label 'process_high'
    container 'broadinstitute/gatk:4.5.0.0'

    publishDir "${params.vqsr}/snp", mode: 'copy'

    input:
      tuple path(genome_file), path(fai_file), path(dict_file)
      tuple path(merged_vcfs), path(merged_vcfs_tbi)
      tuple path(hapmap_vcf), path(hapmap_tbi)
      tuple path(omni_vcf), path(omni_tbi)
      tuple path(onekg_vcf), path(onekg_tbi)

    output:
      tuple path("snp_recalibration.recal"),
            path("snp_recalibration.tranches"),
            path("snp_recalibration.plots.R")

    script:
    """
    gatk --java-options "-Xmx4G -XX:ParallelGCThreads=2" VariantRecalibrator \
      -tranche 100.0 -tranche 99.95 -tranche 99.9 \
      -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
      -tranche 95.0 -tranche 94.0 \
      -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
      -R ${genome_file} \
      -V ${merged_vcfs} \
      --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap_vcf} \
      --resource:omni,known=false,training=true,truth=false,prior=12.0 ${omni_vcf} \
      --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${onekg_vcf} \
      -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
      -mode SNP \
      -O snp_recalibration.recal \
      --tranches-file snp_recalibration.tranches \
      --rscript-file snp_recalibration.plots.R
    """
}


//-----------------------------------------------------------------
// -------- PROCESS TO RUN VARIANT RECALIBRATION (INDELs) ---------
//-----------------------------------------------------------------


/* INFO:
 * Perform Variant Quality Score Recalibration (VQSR) on INDELs using
 * GATK VariantRecalibrator. Builds a statistical model from trusted
 * training resources (Mills + dbSNP) and annotates INDELs with
 * recalibrated quality scores and tranches.
 *
 * Inputs
 *   • genome_file       : reference FASTA (with .fai and .dict)
 *   • merged_vcfs       : cohort-level joint-called VCF (bgzip-compressed)
 *   • mills_vcf         : Mills & 1000G gold standard indels (+ index)
 *   • dbsnp_vcf         : dbSNP known variants resource (+ index)
 *
 * Outputs
 *   • indel_recalibration.recal    : recalibration table
 *   • indel_recalibration.tranches : tranches file for filtering thresholds
 *   • indel_recalibration.plots.R  : R script for generating diagnostic plots
 *
 * Notes
 *   - This step is INDEL-specific (uses `-mode INDEL`).
 *   - dbSNP is used as a "known sites" resource, not for training.
 *   - Should be followed by ApplyVQSR to apply the recalibration model.
 */


process VQSRIndel {
    tag "Variant Recalibration for Indels"
    label 'process_high'
    container 'broadinstitute/gatk:4.5.0.0'

    publishDir "${params.vqsr}/indel", mode: 'copy'

    input:
      tuple path(genome_file), path(fai_file), path(dict_file)
      tuple path(merged_vcfs), path(merged_vcfs_tbi)
      tuple path(mills_vcf), path(mills_tbi)
      tuple path(dbsnp_vcf), path(dbsnp_tbi)

    output:
      tuple path("indel_recalibration.recal"),
            path("indel_recalibration.tranches"),
            path("indel_recalibration.plots.R")

    script:
    """
    gatk --java-options "-Xmx4G -XX:ParallelGCThreads=2" VariantRecalibrator \
      -tranche 100.0 -tranche 99.95 -tranche 99.9 \
      -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
      -tranche 95.0 -tranche 94.0 \
      -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
      -R ${genome_file} \
      -V ${merged_vcfs} \
      --resource:mills,known=false,training=true,truth=true,prior=12.0 ${mills_vcf} \
      --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp_vcf} \
      -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
      -mode INDEL \
      -O indel_recalibration.recal \
      --tranches-file indel_recalibration.tranches \
      --rscript-file indel_recalibration.plots.R
    """
}


//-----------------------------------------------------------------
// --------- PROCESS TO INDEX SNP RECALIBRATION TABLE -------------
//-----------------------------------------------------------------


/* INFO:
 * Create a GATK index (.idx) for the SNP recalibration table generated
 * by VariantRecalibrator. This index is required by downstream steps
 * that consume the recalibration file (e.g., ApplyVQSR).
 *
 * Inputs
 *   • snp_recal_file : recalibration table (snp_recalibration.recal)
 *
 * Outputs
 *   • snp_recal_file.idx : index for the recalibration table
 *
 * Notes
 *   - Must be run after VQSRSNP produces the recalibration table.
 *   - Ensures ApplyVQSR for SNPs can locate and use the .recal file efficiently.
 */


process IndexSNPRecalFile {
    tag "Index SNP recal file"
    label 'process_high'
    container 'broadinstitute/gatk:4.5.0.0'

    publishDir "${params.vqsr}/snp", mode: 'copy'

    input:
      path snp_recal_file

    output:
      path "${snp_recal_file}.idx"

    script:
    """
    gatk IndexFeatureFile -I ${snp_recal_file}
    """
}


//-----------------------------------------------------------------
// ------- PROCESS TO INDEX INDEL RECALIBRATION TABLE -------------
//-----------------------------------------------------------------


/* INFO:
 * Create a GATK index (.idx) for the INDEL recalibration table generated
 * by VariantRecalibrator. This index is required by downstream steps
 * that consume the recalibration file (e.g., ApplyVQSR).
 *
 * Inputs
 *   • indel_recal_file : recalibration table (indel_recalibration.recal)
 *
 * Outputs
 *   • indel_recal_file.idx : index for the recalibration table
 *
 * Notes
 *   - Must be run after VQSRIndel produces the recalibration table.
 *   - Ensures ApplyVQSR for INDELs can locate and use the .recal file efficiently.
 */


process IndexIndelRecalFile {
    tag "Index Indel recal file"
    label 'process_high'
    container 'broadinstitute/gatk:4.5.0.0'

    publishDir "${params.vqsr}/indel", mode: 'copy'

    input:
      path indel_recal_file

    output:
      path "${indel_recal_file}.idx"

    script:
    """
    gatk IndexFeatureFile -I ${indel_recal_file}
    """
}


//-----------------------------------------------------------------
// --------- PROCESS TO APPLY VQSR MODEL TO SNP VARIANTS ----------
//-----------------------------------------------------------------


/* INFO:
 * Apply the SNP recalibration model generated by VariantRecalibrator
 * to the cohort-level VCF. Produces a recalibrated SNP VCF at the
 * specified truth sensitivity threshold.
 *
 * Inputs
 *   • merged_vcfs         : cohort-level joint-called VCF (+ .tbi index)
 *   • snp_recal_file      : SNP recalibration table (.recal) + index (.idx)
 *   • snp_tranches_file   : SNP tranches file (defines filtering thresholds)
 *
 * Outputs
 *   • SNP.recalibrated_99.9.vcf.gz      : recalibrated SNP VCF (bgzipped)
 *   • SNP.recalibrated_99.9.vcf.gz.tbi  : Tabix index for the recalibrated SNP VCF
 *
 * Notes
 *   - Uses `--truth-sensitivity-filter-level 99.9` (default high sensitivity).
 *   - This step should follow VQSRSNP + IndexSNPRecalFile.
 *   - The resulting VCF is ready for downstream filtering or merging with INDELs.
 */


process ApplyVQSRSNP {
    tag "Apply VQSR for SNPs"
    label 'process_high'
    container 'broadinstitute/gatk:4.5.0.0'

    publishDir "${params.vqsr}/snp", mode: 'copy'

    input:
      tuple path(merged_vcfs), path(merged_vcfs_tbi)
      tuple path(snp_recal_file), path(snp_recal_file_idx)
      path snp_tranches_file

    output:
      tuple path("SNP.recalibrated_99.9.vcf.gz"), path("SNP.recalibrated_99.9.vcf.gz.tbi")

    script:
    """
    gatk --java-options "-Djava.io.tmpdir=/tmp -Xms2G -Xmx2G -XX:ParallelGCThreads=2" ApplyVQSR \
      -V ${merged_vcfs} \
      --recal-file ${snp_recal_file} \
      --tranches-file ${snp_tranches_file} \
      --truth-sensitivity-filter-level 99.9 \
      -mode SNP \
      -O SNP.recalibrated_99.9.vcf.gz
    """
}


//-----------------------------------------------------------------
// -------- PROCESS TO APPLY VQSR MODEL TO INDEL VARIANTS ---------
//-----------------------------------------------------------------


/* INFO:
 * Apply the INDEL recalibration model generated by VariantRecalibrator
 * to the SNP-recalibrated VCF. Produces a fully recalibrated VCF that
 * includes both SNP and INDEL quality adjustments.
 *
 * Inputs
 *   • snp_recalibrated_vcf : VCF already recalibrated for SNPs (+ .tbi index)
 *   • indel_recal_file     : INDEL recalibration table (.recal) + index (.idx)
 *   • indel_tranches_file  : INDEL tranches file (defines filtering thresholds)
 *
 * Outputs
 *   • indel.SNP.recalibrated_99.9.vcf.gz      : recalibrated SNP+INDEL VCF (bgzipped)
 *   • indel.SNP.recalibrated_99.9.vcf.gz.tbi  : Tabix index for the recalibrated VCF
 *
 * Notes
 *   - Uses `--truth-sensitivity-filter-level 99.9` (high sensitivity).
 *   - Must be run after ApplyVQSRSNP has produced the SNP-recalibrated VCF.
 *   - This is typically the final VQSR step, yielding the cohort’s
 *     fully recalibrated variant set.
 */


process ApplyVQSRIndel {
    tag "Apply VQSR for INDELs"
    label 'process_high'
    container 'broadinstitute/gatk:4.5.0.0'

    publishDir "${params.vqsr}/indel", mode: 'copy'

    input:
      tuple path(snp_recalibrated_vcf), path(snp_recalibrated_vcf_tbi)
      tuple path(indel_recal_file), path(indel_recal_file_idx)
      path indel_tranches_file

    output:
      tuple path("indel.SNP.recalibrated_99.9.vcf.gz"), path("indel.SNP.recalibrated_99.9.vcf.gz.tbi")

    script:
    """
    gatk --java-options "-Djava.io.tmpdir=/tmp -Xms2G -Xmx2G -XX:ParallelGCThreads=2" ApplyVQSR \
      -V ${snp_recalibrated_vcf} \
      --recal-file ${indel_recal_file} \
      --tranches-file ${indel_tranches_file} \
      --truth-sensitivity-filter-level 99.9 \
      -mode INDEL \
      -O indel.SNP.recalibrated_99.9.vcf.gz
    """
}


//------------------------------------------------------------------
// --------- PROCESS TO SUBSET PASS VARIANTS and TARGET GENES -----
//-----------------------------------------------------------------


/* INFO:
 * Use GATK SelectVariants to extract only high-confidence (PASS) variants
 * from the fully recalibrated VCF and restrict them to a user-provided
 * gene list (BED intervals). Produces a smaller, focused VCF for downstream
 * analysis or reporting.
 *
 * Inputs
 *   • finalVcf    : fully recalibrated cohort VCF (bgzipped) + Tabix index
 *   • genesBed    : BED file with target gene intervals (e.g., AD gene panel)
 *
 * Outputs
 *   • AD_subset_genes.vcf.gz      : subset VCF containing only PASS variants in target genes
 *   • AD_subset_genes.vcf.gz.tbi  : Tabix index for the subset VCF
 *
 * Notes
 *   - Filters out all variants that did not pass VQSR (`--exclude-filtered true`).
 *   - Intervals (`-L`) come from the BED file; allows flexible gene lists.
 *   - Useful for generating clinician-facing or project-specific variant reports.
 */


process SubsetPassVariants {
    tag "Subset + Filter PASS"
    label 'process_medium'
    
    container 'broadinstitute/gatk:4.5.0.0' 
    
    publishDir "${params.outdir}/subset_vcf", mode: 'copy'

    input:
    tuple path(finalVcf), path(finalVcfTbi)  // the final VQSRed, i.e. indel.SNP.recalibrated_99.9.vcf.gz;  VCF + .tbi files
    path genesBed                            // gene list of interest as .bed file

    output:
    tuple path("AD_subset_genes.vcf.gz"), path("AD_subset_genes.vcf.gz.tbi")

    script:
     """
    # GATK SelectVariants: keep only PASS & restrict to AD gene list
    gatk --java-options "-Xmx4g" SelectVariants \
        -V ${finalVcf} \
        --exclude-filtered true \
        -L ${genesBed} \
        -O AD_subset_genes.vcf.gz
    """
}


//------------------------------------------------------------------
// ------------- PROCESS TO RUN SNPEFF FUNCTIONAL ANNOTATION ------
//------------------------------------------------------------------


/* INFO:
 * Annotate variants with predicted functional effects using SnpEff.
 *  ---------- HAS A BUG, STOPPED WORKING & needs fixing -------------------
 * check container before running:
 * https://biocontainers.pro/tools/snpeff
 * Inputs
 *   • vcf     : bgzipped VCF to annotate (e.g., AD_subset_genes.vcf.gz as we have)
 *   • vcf_tbi : Tabix index for the input VCF
 *
 * Outputs
 *   • ${params.out_prefix}.snpeff.vcf.gz      : SnpEff-annotated VCF (bgzipped)
 *   • ${params.out_prefix}.snpeff.vcf.gz.tbi  : Tabix index for the annotated VCF
 *
 * Notes
 *   - Requires the target genome database to be present in `params.snpeff_db_dir`
 *     and referenced by `params.snpeff_genome` (e.g., we use hg38 locally).
 *   - 
 *   - for simplicity we are using -nodownload, so the DB !!must!! already exist locally, 
 *     and the folder name must match params.snpeff_genome.
 *   - keep an eye on the container as it failed
 */


process SnpEffAnnotation {
  tag 'SNPEFF_Annotation'
  label 'process_medium'

  container 'quay.io/biocontainers/snpeff:5.3.0a--hdfd78af_1'
  
  publishDir "${params.outdir}/snpeff_annotated_subset_vcf", mode: 'copy'

  input:
    tuple path(vcf), path(vcf_tbi)
    //path  snpeff_db_dir
    //val   snpeff_genome

  output:
    tuple path("${params.out_prefix}.snpeff.vcf.gz"), path("${params.out_prefix}.snpeff.vcf.gz.tbi")

  script:
  """

  snpEff -v \\
      -dataDir ${params.snpeff_db_dir} \\
      -threads ${params.snpeff_threads} \\
      ${params.snpeff_genome} \\
      ${vcf} \\
      | bgzip -c > ${params.out_prefix}.snpeff.vcf.gz

    tabix -p vcf ${params.out_prefix}.snpeff.vcf.gz
  """
}
