# Whole Exome or Genome Variant Calling
Example of Variant Calling using [GATK](https://gatk.broadinstitute.org/hc/en-us). These bash shell scripts were ran on HPC using [Slurm](https://slurm.schedmd.com/documentation.html) workload manager.  
  
The pipline performs [somatic](https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-) and [germline](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-) variant calling in tumor only and paired tumor-normal mode, depending on workflow.  
  
A generalized overview of the pipeline and tool used:  

 
  1. [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic): Trim Illumina pair-end reads in fastq format 
  2. [BWA-mem](https://github.com/lh3/bwa): Align reads to reference genome. __optionally build reference index with__ _bwa_index.sh_ 
  3. [GATK Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2): Call somatic mutations
  4. [HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller): Call germline mutations
  5. [GVCF](https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format): Perform joint germline mutation calling 
  6. [snpEff](http://pcingola.github.io/SnpEff/): Perform annotation 
  
_details of all steps can be found within the shell scripts_
