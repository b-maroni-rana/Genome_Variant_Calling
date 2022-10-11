#!/bin/bash

# brett maroni-rana (adopted from ting wang)
# date: 20211115

# run with merge_gvcfs_anno.sh

# pipeline details
# 1. GenomicsDBimport - merges gvcfs from germline calling (HaplotypeCaller). must point to non-existant path
# 2. GVCF - performs joint genotyping
# 3. VariatnRecalibrator, snps - build a recalibration model to score variant quality for filtering purposes
# 4. ApplyVQSR, snps - apply a score cutoff to filter variants based on a recalibration table
# 5. VariatnRecalibrator, indels - build a recalibration model to score variant quality for filtering purposes
# 6. ApplyVQSR, indels - apply a score cutoff to filter variants based on a recalibration table
# 7. bcftools view - filter variants for 'PASS', extract exom specific regions with 'regions' option (usually specific to kit, especcially for exome)
# 8. snpEff - annotation. functional and predicted. uses GRCh37.p13, These are RefSeq transcripts from NCBI mapped to GRCh38/hg19 reference genome sequence

#SBATCH -J merge_gvcfs
#SBATCH -p campus-new
#SBATCH -c 10
#SBATCH -t 7-12
#SBATCH --output merge_gvcfs_log-%j.out

infolder="/fh/fast/grady_w/users/bmaroni-rana/variant_calling/gatk/tempus/data_out_hg19"
outfolder="/fh/fast/grady_w/users/bmaroni-rana/variant_calling/gatk/tempus/data_out_hg19/fbxw7/merged/germline"
genomicsdb="/fh/fast/grady_w/users/bmaroni-rana/variant_calling/gatk/tempus/data_out_hg19/fbxw7/genomicsdb" 
tmpfolder="/fh/fast/grady_w/users/bmaroni-rana/variant_calling/gatk/tempus/data_out_hg19/fbxw7/tmp"
intervals="/fh/fast/grady_w/users/twang23/Data/gatk_resource_b37/Broad.human.exome.b37.interval_list"
samplemap="/fh/fast/grady_w/users/bmaroni-rana/variant_calling/gatk/tempus/data_out_hg19/fbxw7/merged/germline/file_germ.gvcf.map"
refgenome="/fh/fast/grady_w/users/twang23/Data/gatk_resource_b37/human_g1k_v37.fasta"
vcfHapmap="/fh/fast/grady_w/users/twang23/Data/gatk_resource_b37/hapmap_3.3.b37.vcf.gz"
vcfOmni="/fh/fast/grady_w/users/twang23/Data/gatk_resource_b37/1000G_omni2.5.b37.vcf.gz"
vcfGlk="/fh/fast/grady_w/users/twang23/Data/gatk_resource_b37/1000G_phase1.snps.high_confidence.b37.vcf.gz"
vcfDbsnp="/fh/fast/grady_w/users/twang23/Data/gatk_resource_b37/dbsnp_138.b37.vcf.gz"
vcfMills="/fh/fast/grady_w/users/twang23/Data/gatk_resource_b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
vcfAxiom="/fh/fast/grady_w/users/twang23/Data/gatk_resource_b37/Axiom_Exome_Plus.genotypes.all_populations.poly.vcf.gz"
regions="/fh/fast/grady_w/users/twang23/Data/gatk_resource_b37/Broad.human.exome.b37.bed"

module load GATK/4.1.8.1-GCCcore-8.3.0-Java-11

echo "=== GenomicsDBImport, super slow ==="
gatk --java-options "-Xms200G -Xmx200G -XX:+UseParallelGC -XX:ParallelGCThreads=${SLURM_JOB_CPUS_PER_NODE}" GenomicsDBImport \
	--reader-threads ${SLURM_JOB_CPUS_PER_NODE} \
	--tmp-dir ${tmpfolder} \
	--genomicsdb-workspace-path ${genomicsdb} \
	--batch-size 50 \
	--sample-name-map ${samplemap} \
	-L 1 \
	-L 2 \
	-L 3 \
	-L 4 \
	-L 5 \
	-L 6 \
	-L 7 \
	-L 8 \
	-L 9 \
	-L 10 \
	-L 11 \
	-L 12 \
	-L 13 \
	-L 14 \
	-L 15 \
	-L 16 \
	-L 17 \
	-L 18 \
	-L 19 \
	-L 20 \
	-L 21 \
	-L 22 \
	-L X \
	-L Y
echo "=== done! ==="
echo ""

echo "=== GenotypeGVCFs, slow ==="
gatk --java-options "-Xmx100G -XX:+UseParallelGC -XX:ParallelGCThreads=${SLURM_JOB_CPUS_PER_NODE}" GenotypeGVCFs \
	-R ${refgenome} \
	-V gendb://${genomicsdb} \
	-O ${outfolder}/merge_germ.vcf.gz \
	--tmp-dir ${tmpfolder}
echo "=== done! ==="
echo ""

echo "=== Variant Quality Score Recalibration for SNPs ==="
module load GATK/4.1.8.1-GCCcore-8.3.0-Java-11
module load R/4.0.2-foss-2019b
gatk --java-options "-Xmx100G -XX:+UseParallelGC -XX:ParallelGCThreads=${SLURM_JOB_CPUS_PER_NODE}" VariantRecalibrator \
	-R ${refgenome} \
	-V ${outfolder}/merge_germ.vcf.gz \
	--max-gaussians 6 \
	--tmp-dir ${tmpfolder} \
	--resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${vcfHapmap} \
	--resource:omni,known=false,training=true,truth=true,prior=12.0 ${vcfOmni} \
	--resource:1000G,known=false,training=true,truth=false,prior=10.0 ${vcfGlk} \
	--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${vcfDbsnp} \
	-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
	-mode SNP \
	-O ${outfolder}/merge_germ_SNP.recal \
	--tranches-file ${outfolder}/merge_germ_recal_SNP.tranches \
	--rscript-file ${outfolder}/merge_germ_recal_SNP.plots.R
echo "=== done! ==="
echo ""

echo "=== Apply variant recalibation to SNPs ==="
gatk --java-options "-Xmx100G -XX:+UseParallelGC -XX:ParallelGCThreads=${SLURM_JOB_CPUS_PER_NODE}" ApplyVQSR \
	-R ${refgenome} \
	-V ${outfolder}/merge_germ.vcf.gz \
	-O ${outfolder}/merge_germ_recal_SNP.vcf.gz \
	-mode SNP \
	--tmp-dir ${tmpfolder} \
	--recal-file ${outfolder}/merge_germ_SNP.recal \
	--tranches-file ${outfolder}/merge_germ_recal_SNP.tranches \
	--truth-sensitivity-filter-level 99.5
echo "=== done! ==="
echo ""

echo "=== Variant Quality Score Recalibration for Indels ==="
module load GATK/4.1.8.1-GCCcore-8.3.0-Java-11
module load R/4.0.2-foss-2019b
gatk --java-options "-Xmx100G -XX:+UseParallelGC -XX:ParallelGCThreads=${SLURM_JOB_CPUS_PER_NODE}" VariantRecalibrator \
	-R ${refgenome} \
	-V ${outfolder}/merge_germ_recal_SNP.vcf.gz \
	-O ${outfolder}/merge_germ_SNP_Indel.recal \
	-mode INDEL \
	--max-gaussians 4 \
	--tranches-file ${outfolder}/merge_germ_recal_SNP_Indel.tranches \
	-an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
	--rscript-file ${outfolder}/merge_germ_recal_SNP_Indel.plots.R \
	--tmp-dir ${tmpfolder} \
	--resource:mills,known=false,training=true,truth=true,prior=12.0 ${vcfMills} \
	--resource:axiomPoly,known=false,training=true,truth=false,prior=10.0 ${vcfAxiom} \
	--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${vcfDbsnp}
echo "=== done! ==="
echo ""

echo "=== Apply variant recalibation to Indels ==="
gatk --java-options "-Xmx100G -XX:+UseParallelGC -XX:ParallelGCThreads=${SLURM_JOB_CPUS_PER_NODE}" ApplyVQSR \
	-R ${refgenome} \
	-V ${outfolder}/merge_germ_recal_SNP.vcf.gz \
	-O ${outfolder}/merge_germ_recal_SNP_Indel.vcf.gz \
	-mode INDEL \
	--tmp-dir ${tmpfolder} \
	--recal-file ${outfolder}/merge_germ_SNP_Indel.recal \
	--tranches-file ${outfolder}/merge_germ_recal_SNP_Indel.tranches \
	--truth-sensitivity-filter-level 99.0
echo "=== done! ==="
echo ""

echo "=== extract specific chromosomes, regions and PASS variants ==="
module load BCFtools/1.9-GCC-8.3.0
bcftools view -f "PASS" -Oz -o ${outfolder}/merge_germ_recal_SNP_Indel_PASS.vcf.gz ${outfolder}/merge_germ_recal_SNP_Indel.vcf.gz
tabix -p vcf ${outfolder}/merge_germ_recal_SNP_Indel_PASS.vcf.gz

bcftools view -R ${regions} -Oz -o ${outfolder}/merge_germ_recal_SNP_Indel_PASS_regionFilt.vcf.gz ${outfolder}/merge_germ_recal_SNP_Indel_PASS.vcf.gz
tabix -p vcf ${outfolder}/merge_germ_recal_SNP_Indel_PASS_regionFilt.vcf.gz
echo "=== done! ==="
echo ""

echo "=== annotation ==="
module load Java/11.0.2
java -Xmx10G -jar /fh/fast/grady_w/users/bmaroni-rana/snpEff/snpEff.jar \
-v -stats ${outfolder}/merge_germ_recal_SNP_Indel_PASS_ann.html \
GRCh37.75 ${outfolder}/merge_germ_recal_SNP_Indel_PASS.vcf.gz > ${outfolder}/merge_germ_recal_SNP_Indel_PASS_ann.vcf

module load tabix/0.2.6-GCCcore-8.3.0
bgzip ${outfolder}/merge_germ_recal_SNP_Indel_PASS_ann.vcf
tabix -p vcf ${outfolder}/merge_germ_recal_SNP_Indel_PASS_ann.vcf.gz

java -Xmx10G -jar /fh/fast/grady_w/users/bmaroni-rana/snpEff/snpEff.jar \
-v -stats ${outfolder}/merge_germ_recal_SNP_Indel_PASS_regionFilt_ann.html \
GRCh37.75 ${outfolder}/merge_germ_recal_SNP_Indel_PASS_regionFilt.vcf.gz > ${outfolder}/merge_germ_recal_SNP_Indel_PASS_regionFilt_ann.vcf

bgzip ${outfolder}/merge_germ_recal_SNP_Indel_PASS_regionFilt_ann.vcf
tabix -p vcf ${outfolder}/merge_germ_recal_SNP_Indel_PASS_regionFilt_ann.vcf.gz
echo "=== done! ==="
echo ""