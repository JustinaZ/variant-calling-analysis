# variant-calling-analysis
Variant Calling Pipeline with Nextflow

Software needs:
- GATK4
- BWA
- Picard Tools
- Samtools
- SnpEff
- R (dependency for some GATK steps)
- Nextflow

Also: 
- Java
- bcftools
- bedtools
  
=============================================

Short Description:

Pipeline is developed for exome sequencing data and works with curated gene list (such as e.g. relevant to Alzheimer's disease).
The goal is to find variants that are "rare", are in coding regions or within 50bp of coding region.
Currently it is set up as an “add-on” tool to operate on provided individual/intermediate genoic vcf files, assuming samples
were already pre-processed (e.g. merged, matefixed, sorted, marked duplicates, etc.)

=============================================

Accepted inputs are:

- The hg38 reference genome, version: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
- The .g.vcf.gz file where extension stands for “Genomic VCF” that is gzip-compressed. This file should contain genomic variant information in the VCF (Variant Call Format) specific to one sample. It includes all potential variants detected during the sequencing and analysis, with additional information about the confidence and quality of these variant calls.
These files are used in the pipeline for subsequent joint genotyping, where multiple .g.vcf.gz files from different samples are combined to make more accurate and robust variant calls.
- The .csi extension indicates 







