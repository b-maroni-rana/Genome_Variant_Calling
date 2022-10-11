#!/bin/bash

# brett maroni-rana (adopted from ting wang)
# date: 20211116

# run with merge_soma_anno.sh

# pipeline details
# 1. GenomicsDBimport - merges gvcfs from germline calling (HaplotypeCaller). must point to non-existant path
# 2. GVCF - performs joint genotyping
# 3. VariatnRecalibrator, snps - build a recalibration model to score variant quality for filtering purposes
# 4. ApplyVQSR, snps - apply a score cutoff to filter variants based on a recalibration table
# 5. VariatnRecalibrator, indels - build a recalibration model to score variant quality for filtering purposes
# 6. ApplyVQSR, indels - apply a score cutoff to filter variants based on a recalibration table
# 7. bcftools view - filter variants for 'PASS', extract exom specific regions with 'regions' option (usually specific to kit, especcially for exome)
# 8. snpEff - annotation. functional and predicted. uses GRCh37.p13, These are RefSeq transcripts from NCBI mapped to GRCh38/hg19 reference genome sequence

#SBATCH -J merge_soma
#SBATCH -p campus-new
#SBATCH -c 10
#SBATCH -t 7-12
#SBATCH --output merge_soma_log-%j.out

infolder="/fh/fast/grady_w/users/bmaroni-rana/variant_calling/gatk/tempus/data_out_hg19/fbxw7"
outfolder="/fh/fast/grady_w/users/bmaroni-rana/variant_calling/gatk/tempus/data_out_hg19/fbxw7/merged/somatic"
regions="/fh/fast/grady_w/users/twang23/Data/gatk_resource_b37/Broad.human.exome.b37.bed"

echo "=== merge vcf files of multiple samples into one ==="
module load VCFtools/0.1.16-foss-2019b-Perl-5.30.0
# move somatic vcfs into one file for use in vcf-merge
mkdir -p ${infolder}/vcf
vcf_files=(`find ${infolder} -maxdepth 2 -name "*soma.filtered.vcf.gz*" -type f | sort`)
for files in "${vcf_files[@]}"
do
    cp $files ${infolder}/vcf
done
vcf-merge ${infolder}/vcf/*.soma.filtered.vcf.gz ${infolder}/vcf/*.soma.filtered.vcf.gz > ${outfolder}/merge_soma_filtered.vcf
echo "=== done! ==="
echo ""

echo "=== vcf.gz and index ==="
module load tabix/0.2.6-GCCcore-8.3.0
bgzip ${outfolder}/merge_soma_filtered.vcf
tabix -p vcf ${outfolder}/merge_soma_filtered.vcf.gz
echo "=== done! ==="
echo ""

echo "=== extract specific chromosomes, regions and PASS variants ===" # "PASS" fails
module load BCFtools/1.9-GCC-8.3.0
bcftools view -f PASS -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y \
-Oz -o ${outfolder}/merge_soma_filtered_newChr_PASS.vcf.gz ${outfolder}/merge_soma_filtered.vcf.gz
tabix -p vcf ${outfolder}/merge_soma_filtered_newChr_PASS.vcf.gz

bcftools view -R ${regions} \
-Oz -o ${outfolder}/merge_soma_filtered_newChr_PASS_regionFilt.vcf.gz ${outfolder}/merge_soma_filtered_newChr_PASS.vcf.gz
tabix -p vcf ${outfolder}/merge_soma_filtered_newChr_PASS_regionFilt.vcf.gz
echo "=== done! ==="
echo ""

echo "=== annotation ==="
module load Java/11.0.2 
java -Xmx10G -jar /fh/fast/grady_w/users/bmaroni-rana/snpEff/snpEff.jar \
-v -stats ${outfolder}/merge_soma_filtered_newChr_PASS_ann.html \
GRCh37.75 ${outfolder}/merge_soma_filtered_newChr_PASS.vcf.gz > ${outfolder}/merge_soma_filtered_newChr_PASS_ann.vcf

module load tabix/0.2.6-GCCcore-8.3.0
bgzip ${outfolder}/merge_soma_filtered_newChr_PASS_ann.vcf
tabix -p vcf ${outfolder}/merge_soma_filtered_newChr_PASS_ann.vcf.gz

java -Xmx10G -jar /fh/fast/grady_w/users/bmaroni-rana/snpEff/snpEff.jar \
-v -stats ${outfolder}/merge_soma_filtered_newChr_PASS_regionFilt_ann.html \
GRCh37.75 ${outfolder}/merge_soma_filtered_newChr_PASS_regionFilt.vcf.gz > ${outfolder}/merge_soma_filtered_newChr_PASS_regionFilt_ann.vcf

bgzip ${outfolder}/merge_soma_filtered_newChr_PASS_regionFilt_ann.vcf
tabix -p vcf ${outfolder}/merge_soma_filtered_newChr_PASS_regionFilt_ann.vcf.gz
echo "=== done! ==="
echo ""