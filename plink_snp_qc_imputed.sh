#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -l h_rt=1:0:0
#$ -j y
#$ -N plink_qc
#$ -o /data/scratch/hmy117
#$ -t 1:22


cd /data/scratch/hmy117/adams_imputed_severity_topmed/

~/plink2 --vcf chr${SGE_TASK_ID}\.dose.vcf.gz dosage=HDS \
--extract-if-info "R2>0.7" \
--double-id \
--hard-call-threshold 0.3 \
--maf 0.01 \
--rm-dup exclude-all \
--make-bed \
--out ADAMS_imputed_tmp_qc_chr${SGE_TASK_ID}

~/plink2 \
--bfile ADAMS_imputed_tmp_qc_chr${SGE_TASK_ID} \
--new-id-max-allele-len 9999 \
--set-all-var-ids @:#:\$r:\$a \
--out ADAMS_imputed_qc_chr${SGE_TASK_ID} \
--make-bed
