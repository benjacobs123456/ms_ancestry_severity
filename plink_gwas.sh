#!/bin/bash
#$ -pe smp 12
#$ -l h_vmem=12G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N plink_gwas
#$ -o /data/scratch/hmy117

cd /data/home/hmy117/ADAMS_severity/
for ancestry in CSA EUR AFR;
do
# plink
~/plink2 \
--bfile ./outputs/imputed_genotypes_$ancestry \
--pheno ./outputs/pheno_plink_$ancestry\.tsv \
--covar ./outputs/cov_plink_$ancestry\.tsv \
--glm hide-covar \
--covar-variance-standardize \
--out ./outputs/plink_gwas_$ancestry
done





