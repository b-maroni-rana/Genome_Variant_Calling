#!/bin/bash

# pipeline_vcf_tumor_hg19.sh
# author : brett maroni-rana (adopted from ting wang)
# title : Tumor-Only Mode, Single Patient (sample) Somatic Variant and Germline Caller - GRCh37
# 20211109

# run with sbatch pipeline_vcf_tumor_hg19.sh

# !script uses hg19 resources as input!
# full pipeline - change {fastqfolder} appropriately
# runs in tumor-only mode (not tumor-normal mode)
# does not process vcfs downstream
# run additional programs for group variant calling and annotation

# GATK pipeline: 
# 0. OPTIONAL and NOT INCLUDED. index reference gemone with bwa index (check dir for .amb, .ann, .bwt, .pac and .sa)
# 1. trim with trimmomatic
# 2. align reads with bwa mem and convert .sam to .bam
# 3. sort .bam and mark duplicates with picard, sorts alignments usually by coordinates and marks (usually) PCR duplicates
# 4. perform base call recalibration using prior snp data with GATK BaseRecalibrator
# 5. perform somatic mutation calling with GATK Mutect2
# 6. lean orientation bias artifacts with GATK LearnReadOrientationModel
# 7. calculate sample contamination with GATK GetPileUpSummaries and Calculate Contamination
# 8. filter variants with GATK FilterMutectCalls, runs statistics on different bases and generates F score threshold for calls
# 9. perform germline mutation calling with GATK HaplotypeCaller
# 10. perform joint mutation calling with GATK GenotypeGVCF, combines germline calls per-sample increase calling power

#SBATCH -J vcf_tumor19
#SBATCH --partition=campus-new
#SBATCH -c 8
#SBATCH -t 7-12
#SBATCH -o pipeline_vcf_tumor_hg19_log_%j.out
# #SBATCH -e pipeline_vcf_tumor_hg19_error_%j.out

fastqfolder="/fh/fast/grady_w/users/bmaroni-rana/variant_calling/gatk/tempus/fastq/tumor_only"
outfolder="/fh/fast/grady_w/users/bmaroni-rana/variant_calling/gatk/tempus/data_out_hg19/tumor_only"
adapters="/app/software/Trimmomatic/0.39-Java-11/adapters/TruSeq3-PE-2.fa"
refgenome="/fh/fast/grady_w/users/twang23/Data/gatk_resource_b37/human_g1k_v37.fasta"
dbsnp="/fh/fast/grady_w/users/twang23/Data/gatk_resource_b37/dbsnp_138.b37.vcf.gz"
vcfGlk="/fh/fast/grady_w/users/twang23/Data/gatk_resource_b37/1000G_phase1.snps.high_confidence.b37.vcf.gz"
germlineresource="/fh/fast/grady_w/users/twang23/Data/Mutect2_data/b37/af-only-gnomad.raw.sites.vcf.gz"
ponvcf="/fh/fast/grady_w/users/twang23/Data/Mutect2_data/b37/Mutect2-exome-panel.vcf.gz"
germlinecommon="/fh/fast/grady_w/users/twang23/Data/Mutect2_data/b37/small_exac_common_3.vcf.gz"

fastq_files=(`find ${fastqfolder} -maxdepth 1 -name "*.fastq.gz" -type f | sort`)
for ((i=0; i<${#fastq_files[@]}; i+=2)); do
		subfolder=$(echo ${fastq_files[$i]} | rev | cut -d "/" -f1 | rev | cut -d "." -f1 | cut -d "_" -f1,2 ) # works!
		#echo $subfolder # check files work
		mkdir -p ${outfolder}/${subfolder}
		echo "=== analyzing for sample ${subfolder} ==="
		echo "=== cut adapters and trim reads for sample ${subfolder} ==="
		module load Trimmomatic/0.39-Java-11
		# trim reads
		java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads ${SLURM_JOB_CPUS_PER_NODE} \
			-phred33 ${fastq_files[$i]} ${fastq_files[$i+1]} \
			${outfolder}/${subfolder}/${subfolder}.R1_P_T.fq.gz ${outfolder}/${subfolder}/${subfolder}.R1_UP_T.fq.gz \
			${outfolder}/${subfolder}/${subfolder}.R2_P_T.fq.gz ${outfolder}/${subfolder}/${subfolder}.R2_UP_T.fq.gz \
			ILLUMINACLIP:${adapters}:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
			2> ${outfolder}/${subfolder}/${subfolder}.trimmomatic.log
		echo "=== DONE! ==="
		echo ""

		echo "=== mapping to reference genome with BWA MEM for sample ${subfolder} ==="
		module load BWA/0.7.17-GCC-8.3.0
		module load SAMtools/1.10-GCCcore-8.3.0
		# alignment
		bwa mem -t ${SLURM_JOB_CPUS_PER_NODE} -M -T 0 -R @RG\\tID:${subfolder}\\tSM:${subfolder}\\tLB:${subfolder}\\tPL:ILLUMINA\\tPU:${subfolder} \
			${refgenome} \
			${outfolder}/${subfolder}/${subfolder}.R1_P_T.fq.gz \
			${outfolder}/${subfolder}/${subfolder}.R2_P_T.fq.gz | samtools view -Shb -o ${outfolder}/${subfolder}/${subfolder}.bam -
		echo "=== DONE! ==="
		echo ""
		
		echo "=== sort and mark duplicates in bam file for sample ${subfolder} ==="
		module load picard/2.21.6-Java-11
		mkdir -p ${outfolder}/${subfolder}/${subfolder}.tmp
		# create index
		java -Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=${SLURM_JOB_CPUS_PER_NODE} -jar $EBROOTPICARD/picard.jar \
			SortSam CREATE_INDEX=true \
			INPUT=${outfolder}/${subfolder}/${subfolder}.bam \
			OUTPUT=${outfolder}/${subfolder}/${subfolder}.sorted.bam \
			SORT_ORDER=coordinate VALIDATION_STRINGENCY=STRICT TMP_DIR=${outfolder}/${subfolder}/${subfolder}.tmp
		# mark dup
		java -Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=${SLURM_JOB_CPUS_PER_NODE} -jar $EBROOTPICARD/picard.jar \
			MarkDuplicates CREATE_INDEX=true \
			INPUT=${outfolder}/${subfolder}/${subfolder}.sorted.bam \
			OUTPUT=${outfolder}/${subfolder}/${subfolder}.sorted.markDup.bam \
			METRICS_FILE=${outfolder}/${subfolder}/${subfolder}.sorted.markDup.metrics \
			ASSUME_SORT_ORDER=coordinate \
			VALIDATION_STRINGENCY=STRICT \
			TMP_DIR=${outfolder}/${subfolder}/${subfolder}.tmp
		# align metrics
		java -Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=${SLURM_JOB_CPUS_PER_NODE} -jar $EBROOTPICARD/picard.jar \
			CollectAlignmentSummaryMetrics R=${refgenome} \
			I=${outfolder}/${subfolder}/${subfolder}.sorted.markDup.bam \
			O=${outfolder}/${subfolder}/${subfolder}.sorted.markDup.AlignSummary \
			TMP_DIR=${outfolder}/${subfolder}/${subfolder}.tmp
		echo "=== DONE! ==="
		echo ""
		
		echo "=== base quality score recalibration for sample ${subfolder} ==="
		module load GATK/4.1.8.1-GCCcore-8.3.0-Java-11
		# base recal
		gatk --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=${SLURM_JOB_CPUS_PER_NODE}" \
			BaseRecalibrator -R ${refgenome} \
			-I ${outfolder}/${subfolder}/${subfolder}.sorted.markDup.bam \
			--known-sites ${dbsnp} --known-sites ${vcfGlk} \
			-O ${outfolder}/${subfolder}/${subfolder}.sorted.markDup.recal.table \
			--tmp-dir ${outfolder}/${subfolder}/${subfolder}.tmp
		# apply base recal
		gatk --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=${SLURM_JOB_CPUS_PER_NODE}" \
			ApplyBQSR -R ${refgenome} \
			-I ${outfolder}/${subfolder}/${subfolder}.sorted.markDup.bam \
			--bqsr-recal-file ${outfolder}/${subfolder}/${subfolder}.sorted.markDup.recal.table \
			-O ${outfolder}/${subfolder}/${subfolder}.sorted.markDup.recal.bam \
			--tmp-dir ${outfolder}/${subfolder}/${subfolder}.tmp
		echo "=== DONE! ==="
		echo ""
		
		echo "=== call somatic variants for sample ${subfolder} ==="
		module load GATK/4.1.8.1-GCCcore-8.3.0-Java-11
		# mutect2 tumor-only mode
		gatk --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=${SLURM_JOB_CPUS_PER_NODE}" \
			Mutect2 -R ${refgenome} \
			-I ${outfolder}/${subfolder}/${subfolder}.sorted.markDup.recal.bam \
			-O ${outfolder}/${subfolder}/${subfolder}.soma.vcf.gz \
			--f1r2-tar-gz ${outfolder}/${subfolder}/${subfolder}.f1r2.tar.gz \
			--germline-resource ${germlineresource} \
			--panel-of-normals ${ponvcf} \
			--native-pair-hmm-threads ${SLURM_JOB_CPUS_PER_NODE} \
			--tmp-dir ${outfolder}/${subfolder}/${subfolder}.tmp
		echo "=== DONE! ==="
		echo ""
		
		echo "=== learn read orientation model for sample ${subfolder} ==="
		# read orientation
		gatk --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=${SLURM_JOB_CPUS_PER_NODE}" \
			LearnReadOrientationModel \
			-I ${outfolder}/${subfolder}/${subfolder}.f1r2.tar.gz \
			-O ${outfolder}/${subfolder}/${subfolder}.readOrientationModel.tar.gz \
			--tmp-dir ${outfolder}/${subfolder}/${subfolder}.tmp
		echo "=== DONE! ==="
		echo ""
		
		echo "=== get pileup summaries for sample ${subfolder} ==="
		# pileup
		gatk --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=${SLURM_JOB_CPUS_PER_NODE}" \
			GetPileupSummaries \
			-I ${outfolder}/${subfolder}/${subfolder}.sorted.markDup.recal.bam \
			-O ${outfolder}/${subfolder}/${subfolder}.sorted.markDup.recal.pileupSummary.table \
			-V ${germlinecommon} \
			-L ${germlinecommon} \
			--tmp-dir ${outfolder}/${subfolder}/${subfolder}.tmp
		echo "=== DONE! ==="
		echo ""
		
		echo "=== calculate contamination for sample ${subfolder} ==="
		# 
		gatk --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=${SLURM_JOB_CPUS_PER_NODE}" \
			CalculateContamination \
			-I ${outfolder}/${subfolder}/${subfolder}.sorted.markDup.recal.pileupSummary.table \
			-O ${outfolder}/${subfolder}/${subfolder}.sorted.markDup.recal.contamination.table \
			--tumor-segmentation ${outfolder}/${subfolder}/${subfolder}.sorted.markDup.recal.segments.table \
			--tmp-dir ${outfolder}/${subfolder}/${subfolder}.tmp
		echo "=== DONE! ==="
		echo ""
		
		echo "=== filter mutect2 variants for sample ${subfolder} ==="
		# filter calls
		gatk --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=${SLURM_JOB_CPUS_PER_NODE}" \
			FilterMutectCalls \
			-R ${refgenome} \
			-V ${outfolder}/${subfolder}/${subfolder}.soma.vcf.gz \
			-O ${outfolder}/${subfolder}/${subfolder}.soma.filtered.vcf.gz \
			--contamination-table ${outfolder}/${subfolder}/${subfolder}.sorted.markDup.recal.contamination.table \
			--tumor-segmentation ${outfolder}/${subfolder}/${subfolder}.sorted.markDup.recal.segments.table \
			--orientation-bias-artifact-priors ${outfolder}/${subfolder}/${subfolder}.readOrientationModel.tar.gz \
			--tmp-dir ${outfolder}/${subfolder}/${subfolder}.tmp
		echo "=== DONE! ==="
		echo ""
		
		echo "=== call germline variants for sample ${subfolder} ==="
		# germline caller
		gatk --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=${SLURM_JOB_CPUS_PER_NODE}" \
			HaplotypeCaller \
			-R ${refgenome} \
			-I ${outfolder}/${subfolder}/${subfolder}.sorted.markDup.recal.bam \
			-O ${outfolder}/${subfolder}/${subfolder}.germ.g.vcf.gz \
			-ERC GVCF \
			--native-pair-hmm-threads ${SLURM_JOB_CPUS_PER_NODE} \
			--tmp-dir ${outfolder}/${subfolder}/${subfolder}.tmp
		echo "=== DONE! ==="
		echo ""

		echo "=== GenotypeGVCFs for sample ${subfolder} ==="
		# convert to gvcf
		gatk --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=${SLURM_JOB_CPUS_PER_NODE}" \
			GenotypeGVCFs \
			-R ${refgenome} \
			-V ${outfolder}/${subfolder}/${subfolder}.germ.g.vcf.gz \
			-O ${outfolder}/${subfolder}/${subfolder}.germ.vcf.gz \
			--tmp-dir ${outfolder}/${subfolder}/${subfolder}.tmp # this runs per samples, loses group germline calling power. process downstream with GenomicsDBImport and GenotypeGVCFs next (different script)
		echo "=== Done! ==="
        rm ${outfolder}/${subfolder}/${subfolder}.bam
		rm ${outfolder}/${subfolder}/${subfolder}.sorted.ba*
		rm ${outfolder}/${subfolder}/*fq.gz
done