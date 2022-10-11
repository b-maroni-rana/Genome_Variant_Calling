#!/bin/bash

#SBATCH -J index
#SBATCH --partition=campus-new
#SBATCH -c 4
#SBATCH -t 2-12
#SBATCH -o index_log_%j.out
# #SBATCH -e index_error_%j.out

# run with sbatch test_index.sh

# reference sequence must be indexed for use in bwa mem
# index reference fasta with bwa index
# output is ref.fa{.amb, .ann, .bwt, .pac, .sa} and used in alignment

refgenome="/fh/fast/grady_w/users/twang23/Data/gatk_resource_b37/human_g1k_v37.fasta"

module load BWA/0.7.17-GCC-8.3.0
bwa index $refgenome