#!/bin/bash
#SBATCH --mem=40gb	 		# Job memory request

# $1 = $inter_dir/$prefix
# $2 = $input_dir/$ref.fa

# Checklist: Have you...
# 1) ...installed and added the paths to bwa, samtools, and GATK?
# 2) ...checked the version of Java available on the compute nodes? 
# 3) ...checked your FASTQ/BAM files for read-group (@RG) information?
# 4) ...altered the GATK variant-calling parameters, if necessary?
# 5) ...customized the amount of memory needed to run GATK HaplotypeCaller?

# Step 1.0) Standardize filenames. Make the pool part of the prefix of the input BAM and output gVCF files from the job array ID. 
prefix="${1}_p$SLURM_ARRAY_TASK_ID"
infastq="${prefix}_read"
interbam="${prefix}.rg.bam"
outgvcf="${prefix}.raw.g.vcf"

# Step 1.1.1) For each pool, align the simulated reads to a reference sequence, 
/path/to/bwa mem $2 ${infastq}1.fastq ${infastq}2.fastq | /path/to/samtools view -Shub - > ${prefix}.bam
/path/to/samtools sort ${prefix}.bam ${prefix}.srt

# Step 1.1.2) For each BAM file, add read-group information if it doesn't exist. GATK requires read group information. 
/path/to/java -jar /path/to/gatk AddOrReplaceReadGroups -I ${prefix}.srt.bam -O ${interbam} -ID ${SLURM_ARRAY_TASK_ID} -LB NPD -PL Illumina -PU NPD -SM ${SLURM_ARRAY_TASK_ID} 
/path/to/samtools index ${interbam}

# Step 1.2) For each pool, call segregating sites from the sequencing reads relative to the reference sequence.
/path/to/java -Xmx20g -XX:+UseConcMarkSweepGC -XX:ParallelGCThreads=4 -jar /gpfs/common/programs/gatk-4.0.5.1/gatk-package-4.0.5.1-local.jar HaplotypeCaller -R $2 -I ${interbam} -ERC GVCF -ploidy 150 --heterozygosity 0.01 --max-alternate-alleles 1 -O ${outgvcf}
