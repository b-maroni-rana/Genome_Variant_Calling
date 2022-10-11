# Whole Exome or Genome Variant Calling
Example of Variant Calling shell scripts using GATK.  
  
These bash shell scripts were ran on HPC using [Slurm](https://slurm.schedmd.com/documentation.html) workload manager.  
  
The pipline uses [GATK](https://gatk.broadinstitute.org/hc/en-us) to perform [somatic](https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-) and [germline](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-) variant calling in tumor only mode and paired tumor-normal mode, depending on workflow.  
  
The pipeline runs a number of tools but overall processes is below:  

 
  1. Trim reads with [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
  2. Align reads with [BWA-mem](https://github.com/lh3/bwa)
  3. Call somatic mutations with [GATK Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2)
  4. Call germline mutations with [HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller)
  5. Performs joint germline mutation calling with [GVCF](https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format)
  6. Performs annotation with [snpEff](http://pcingola.github.io/SnpEff/)
  
_details of all steps can be found within the shell scripts_
