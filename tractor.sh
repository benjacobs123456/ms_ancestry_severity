#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N tractor
#$ -o /data/scratch/hmy117
#$ -t 1:22

cd /data/scratch/hmy117/rfmix

module load miniforge
mamba activate tractor

# extract haps
cd /data/scratch/hmy117/rfmix



python3 ~/Tractor/scripts/extract_tracts.py \
--vcf chr${SGE_TASK_ID}\.vcf \
--msp /data/scratch/hmy117/rfmix/rfmix_local_ancestry_chr${SGE_TASK_ID}\.msp.tsv \
--num-ancs 7 \
--output-dir /data/scratch/hmy117/rfmix
