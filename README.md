# Rare-variant-identification
Variant Calling Pipeline + Rare variant identification with Nextflow.

It was built to run on Macs; specifically tested on:  MacBook Pro, M4Pro Chip, 48GB Memory 


Current target for the workflow is summarised in here:

<img width="11063" height="3848" alt="target_pipeline" src="https://github.com/user-attachments/assets/1c10fa3b-ba07-460d-a5a2-3a0fe043e459" />


NOTE: code after snpeff step should be cleaned

=============================================

Software needs (local instalation):
- Nexflow
- Docker
- Java (Nextflow itself requires Java to run)


Relevant Docker BioContainers:  
- GATK4
- BWA
- Picard Tools
- Samtools
- SnpEff
- bcftools
- bedtools

To run it using terminal/command line:

nextflow run main.nf

=============================================

Short Description:

Pipeline is developed for exome sequencing data and works with curated gene list (such as e.g. genes relevant to the Alzheimer's disease).
The goal is to find variants that are "rare", are in coding regions or within 50bp of coding region.
Currently it is set up as an “add-on” tool to operate on provided individual/intermediate genomic vcf files, assuming samples
were already pre-processed (e.g. merged, matefixed, sorted, marked duplicates, etc.)


=============================================

Input Files:

-	The hg38 reference genome, version: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
  
-	The .gvcf.gz sample files; and (optional) corresponding .tbi sample files
  (The .g.vcf.gz file where extension stands for “Genomic VCF” that is gzip-compressed. This file should contain genomic variant information in the VCF (Variant Call Format) specific to one sample. It includes all potential variants detected during the sequencing and analysis, with additional information about the confidence and quality of these variant calls.
These files are used in the pipeline for subsequent joint genotyping, where multiple .g.vcf.gz files from different samples are combined to make more accurate and robust variant calls.)

-	The .bed file that is matching samples (AD case study: NWGC provided the relevant one)
  
-	Other supporting files/resources (via gsutil tool):
  
  o	Download HapMap (for SNP):
  
    gsutil cp gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz 
    
    gsutil cp gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi

  o	Download Omni (for SNP):
  
    gsutil cp gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz 
    
    gsutil cp gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi

  o	Download 1000 Genomes High-Confidence SNPs (for SNP):
  
    gsutil cp gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
    
    gsutil cp gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi

  o	Download Mills and 1000G Gold Standard INDELs (for INDEL):
  
    gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz 
    
    gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi 

  o	Download dbSNP (for INDEL):
  
    gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf 
    
    gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

  o	Download genomAD (for allele frequency info: https://gnomad.broadinstitute.org/downloads): 
  
    gsutil cp gs://gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz

p.s. or manually download from here:

 https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?inv=1&invt=AbmUZw&prefix=&forceOnObjectsSortingFiltering=false

=============================================

Output Files:

All generated output (final and intermediate files) will appear in "output" directory tree. 

====================================================================





