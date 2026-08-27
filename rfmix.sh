#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=64G
#$ -l h_rt=240:0:0
#$ -j y
#$ -N rfmix
#$ -o /data/scratch/hmy117
#$ -t 1:22

cd /data/scratch/hmy117/rfmix

# download reference haps
~/google-cloud-sdk/bin/gsutil cp gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2/hgdp1kgp_chr${SGE_TASK_ID}\.filtered.SNV_INDEL.phased.shapeit5.bcf* ./

# make ref vcf 
~/plink2 --bcf hgdp1kgp_chr${SGE_TASK_ID}\.filtered.SNV_INDEL.phased.shapeit5.bcf \
--recode vcf \
--geno 0.1 \
--maf 0.01 \
--hwe 1e-20 \
--mac 1 \
--out ref_chr${SGE_TASK_ID}

# make input vcf 
~/plink2 --vcf /data/scratch/hmy117/adams_imputed_severity_topmed/chr${SGE_TASK_ID}\.empiricalDose.vcf.gz \
--recode vcf \
--out chr${SGE_TASK_ID} \
--double-id \
--maf 0.05 \
--rm-dup exclude-all \
--snps-only just-acgt \
--geno 0.1 \
--hwe 1e-10 \
--mind 0.1

module load bcftools
~/rfmix/rfmix \
-r ref_chr${SGE_TASK_ID}\.vcf \
-f chr${SGE_TASK_ID}\.vcf \
-m sample_ref \
--chromosome=${SGE_TASK_ID} \
-g genetic_map \
-o rfmix_local_ancestry_chr${SGE_TASK_ID}

wait 
rm ref_chr${SGE_TASK_ID}\.vcf

