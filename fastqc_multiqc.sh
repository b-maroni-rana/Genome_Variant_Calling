#!/bin/bash

# brett maroni-rana (adapted from Ting Want)
# 20211103

# runs fastqc then multiqc on the fastqc files
# run with sbatch 1_qc_bmr.sh

#SBATCH -J qc
#SBATCH --partition=campus-new
#SBATCH -c 4
#SBATCH -t 2-12
#SBATCH -o log-%j.out
#SBATCH -e error-%j.out

module load FastQC/0.11.9-Java-11

fastq_files=(`find /fh/fast/grady_w/users/bmaroni-rana/variant_calling/gatk/tempus/fastq/ -maxdepth 2 -name "*.fastq.gz" -type f | sort`)
for file1 in ${fastq_files[@]}; do fastqc $file1 -t ${SLURM_JOB_CPUS_PER_NODE} -o /fh/fast/grady_w/users/bmaroni-rana/variant_calling/gatk/tempus/qc/; done

module load MultiQC/1.9-foss-2019b-Python-3.7.4

multiqc -v -o /fh/fast/grady_w/users/bmaroni-rana/variant_calling/gatk/tempus/qc/multiqc/ /fh/fast/grady_w/users/bmaroni-rana/variant_calling/gatk/tempus/qc/

mv *.out out
