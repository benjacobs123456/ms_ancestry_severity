#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -l h_rt=240:0:0
#$ -j y
#$ -N tractor_gwas
#$ -o /data/scratch/hmy117
#$ -t 1:22

cd /data/scratch/hmy117/rfmix

module unload miniforge
module load R/4.4.1

#
for pheno in gARMSS_rint uGMSSS_rint edss_rint eq5d_vas_rint age_at_dx_rint msis_physical_normalised_rint;
  do
    ~/Tractor/scripts/run_tractor.R \
--hapdose /data/scratch/hmy117/rfmix/chr${SGE_TASK_ID} \
--phenofile /data/home/hmy117/ADAMS_severity/outputs/tractor_pheno.tsv \
--covarcollist ageatedss,Sex,batch,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
--phenocol $pheno \
--method linear \
--output /data/scratch/hmy117/tractor/tractor_gwas_$pheno\_${SGE_TASK_ID}
done

