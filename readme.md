# Preamble
This repository contains code used to generate the findings presented in 'Association between genetic ancestry and Multiple Sclerosis severity' published in _Annals Of Neurology_.

Per-ancestry GWAS summary statistics can be downloaded from the GWAS catalogue.

Date: 27-08-26.

Author: Ben Jacobs

Email: b.jacobs@qmul.ac.uk

Contents:
- [Imputation](#Imputation)
- [Genotype QC](#QC)
- [Phenotype QC & analysis](#Phenotype)
- [GWAS](#GWAS)

# Imputation
## Imputation prep
````unix
cd /data/home/hmy117/ADAMS_severity/

for i in {1..22};
do
~/plink2 --bfile /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ADAMS_geno_fid_iid \
--set-all-var-ids chr@:#:\$r\:\$a \
--make-bed \
--export vcf \
--snps-only just-acgt \
--chr $i \
--output-chr chrMT \
--rm-dup exclude-all \
--out ./scratch/adams_hg38_cpra_chr$i
done

# convert to vcf 

### Sort with bcf tools
module load bcftools
for i in {1..22};
do
  bcftools sort ./scratch/adams_hg38_cpra_chr$i\.vcf \
  -Oz -o ./imputation_raw_files/sorted_chr$i\.vcf.gz &
done
````

## Imputation 
Download files. 

Attempt imputation via TOPMED-r3 server.

TOPMED-r3 panel. 

R2 filter 0.001.

Eagle 2.4.

1st pass fails due to strand flips.

Download snps-excluded.txt & exclude SNPs as follows:

````unix
cut -f1 "/data/home/hmy117/ADAMS_severity/imputation_raw_files/snps_excluded.txt" > /data/home/hmy117/ADAMS_severity/imputation_raw_files/snps_to_exclude

for i in {1..22};
do
~/plink2 --bfile /data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ADAMS_geno_fid_iid \
--set-all-var-ids @:#:\$r\:\$a \
--make-bed \
--export vcf \
--snps-only just-acgt \
--chr $i \
--output-chr chrMT \
--rm-dup exclude-all \
--exclude /data/home/hmy117/ADAMS_severity/imputation_raw_files/snps_to_exclude \
--out ./scratch/adams_hg38_cpra_2ndpass_chr$i
done

# convert to vcf 

### Sort with bcf tools
module load bcftools
for i in {1..22};
do
  bcftools sort ./scratch/adams_hg38_cpra_2ndpass_chr$i\.vcf \
  -Oz -o ./imputation_raw_files/2nd_pass/sorted_chr$i\.vcf.gz &
done

````
Then resubmit the job and try again. 


## Download imputed data

````unix
mkdir /data/scratch/hmy117/adams_imputed_severity_topmed
cd /data/scratch/hmy117/adams_imputed_severity_topmed

curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1756372/a0ecd1d5174cad443fbbbbba32e2e5d67fd4e3bf96150204469c02e0656d4cc1 | bash

# copy wget commands 
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/b711773ac70528fd5b0ab90861fa32b22e8327e46478490cd17ecf2fb67be0a4/chr_1.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/bd758b6d65bed051edde01a08f51b4233efdf03ffdc3c5c66e435b4f1d7b7563/chr_10.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/1144fe1e81dcb260e266ccb6132ef2d9ea0b3669f7365bc735f76c36504ef1bf/chr_11.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/53a2f6f2a6dbffaebb1c38c3ae3002f2d0c69bf335144f92e97e1ec3711a1381/chr_12.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/e4a946942a9cafef38cb875b531f147b56a7f30f147362f4e4ee0a2940c36678/chr_13.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/9e7bfcf741cda71030e916f373d97e2233c1c0512b06cbdcc300253c7de23c6f/chr_14.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/d2e03e777057a8b4312c582cbd368e181a299a79ed619e3ea5a0bfafa259da06/chr_15.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/9f14641073f29d4e12f8d33d9d52e00398a3e404bafe48b7039519fd59df4ca9/chr_16.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/c548a3cd2786546f0b2ef52c724b275e6372489ef9e30b04c9b23128f3dc2b1a/chr_17.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/d38d659b696ab520787c192110a70e101fbff66c910c48a16ead3b470e2b5fb1/chr_18.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/cdfc31c27ced10f1a1019e1ee6c149d0ca00d504910d34c697cb77bdb134a8aa/chr_19.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/28b7ef367e24d88128d1fa9f8e0279d86a4ed105ecbdf0555ef2724c51fcc010/chr_2.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/650f18ee3b535db5ab0b129235177a6a39e75e5a78e1c5c9ed14d748afc5481b/chr_20.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/7ec5983654165c123bc4127e1df7b5eeb2897ad1c01a21c17cf08e8eb538c435/chr_21.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/1343bbe4a3f969d25931a40fa4f496b571ffcc612aba8c1b1a9d8dd1f2c76cf5/chr_22.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/898c496e852b94432ce330a3e9196774e047fc43f9a2a46857a96c466fb120a7/chr_3.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/865da255d2d73dda8b59a6860fe01246ba462cc14efac5f61d70b22845d56363/chr_4.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/47fa46685014c3455fc9f97736ae772558dc336b504fe7da3415b5c636e6cf0d/chr_5.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/fe1ba842f07167289316396c7979b7cf15502cad1d7d7e434697cf15443f36bd/chr_6.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/883db49f07538b7ce2bbc8521cee53c52efc6a9df39db4e506ca91e8446ef01c/chr_7.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/170f32b157d58dad336ebfcb4d5da86bc01521f5a5b125c8a1b0a9d20f5c1d8b/chr_8.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/581417180d0d69630ef374bf774a0cf018be8a2e992e1f2406bc4e1945e58e6f/chr_9.zip &
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/c86d148c881b194c61fee3662e114fc58e8229323dde85c421fcd90bfa69fe40/results.md5 &


## Unzip
for i in {1..22};
do
  unzip -P VtqN8svNP#2aEz -o chr_$i\.zip &
done

````

# QC
- SNP QC on individual VCFs
- Conversion back to plink, hard call with threshold 0.3

## Get imputation scores
````unix 

module load bcftools
echo "CHR POS ID MAF R2 ER2 TYPED IMPUTED" > /data/scratch/hmy117/snp_info_data
for i in {1..22};
do
bcftools query -f "%CHROM\t%POS\t%ID\t%INFO/MAF\t%INFO/R2\t%INFO/ER2\t%INFO/TYPED\t%INFO/IMPUTED" /data/scratch/hmy117/adams_imputed_severity_topmed/chr$i\.info.gz >> /data/scratch/hmy117/snp_info_data
done 
````

## Explore in R 
````R 
library(tidyverse)
dat = read_table("/data/scratch/hmy117/snp_info_data")

# filter 
dat = dat %>% filter(R2>0.7 & MAF > 0.01)

````

## Run QC

````unix

qsub /data/home/hmy117/ADAMS_severity/scripts/plink_snp_qc_imputed.sh
````

## Update fam files
- Modify IDs in R
````R
library(tidyverse)
setwd("/data/scratch/hmy117/adams_imputed_severity_topmed/")

for(i in c(1:22)){
df = read_table(paste0("ADAMS_imputed_qc_chr",i,".fam"),col_names=F) %>%
	mutate(X1 = str_remove(X1,"^0_")) %>%
	tidyr::separate(X1,sep = "_",into = c("part1","oragene")) %>%
	dplyr::select(X2,oragene) %>%
	mutate(oldfid = X2, oldiid = X2, newfid = oragene, newiid = oragene) %>%
	dplyr::select(oldfid,oldiid,newfid,newiid)
write_tsv(df,paste0("chr",i,"_newids.tsv"),col_names = F)
}

````

## Update IDs in PLINK
- Merge across chromosomes in PLINK
- Rename sample IDs
````unix
# rename IDs
cd /data/scratch/hmy117/adams_imputed_severity_topmed/
for i in {22..1}; do ~/plink --bfile ADAMS_imputed_qc_chr$i --update-ids chr$i\_newids.tsv --out chr$i\_combined_adams_imputed_newids --make-bed; done
````

## Merge chromosomes
````unix
cd /data/scratch/hmy117/adams_imputed_severity_topmed/
rm filelist_for_merge
for i in {2..22}; do echo chr$i\_combined_adams_imputed_newids >> filelist_for_merge; done
~/plink --bfile chr1_combined_adams_imputed_newids \
--merge-list filelist_for_merge \
--out combined_adams_imputed \
--make-bed

````

## Further QC
````unix
cd /data/scratch/hmy117/adams_imputed_severity_topmed/

~/plink --bfile combined_adams_imputed \
--snps-only just-acgt \
--make-bed \
--geno 0.1 \
--hwe 1e-10 \
--mind 0.1 \
--maf 0.01 \
--out combined_adams_imputed_qc \
--chr 1-22

wc combined_adams_imputed_qc.bim

````

## Individual QC 

### Heterozygosity
````unix
cd /data/scratch/hmy117/adams_imputed_severity_topmed/
~/plink --bfile combined_adams_imputed_qc \
--indep-pairwise 1000 500 0.1 \
--maf 0.05 \
--out pruned_for_het_and_kinship

~/plink --bfile combined_adams_imputed_qc \
--extract pruned_for_het_and_kinship.prune.in \
--make-bed \
--out pruned_for_het_and_kinship_genotypes

~/plink --bfile pruned_for_het_and_kinship_genotypes \
--het small-sample \
--out het_check
````


### Ancestry inference

#### Reference preparation

##### AJ ref 

Note this in in hg18

````unix
mkdir /data/scratch/hmy117/aj_ref 

cd /data/scratch/hmy117/aj_ref
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE23nnn/GSE23636/suppl/GSE23636%5FAJ%5F471%5F732k.map.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE23nnn/GSE23636/suppl/GSE23636%5FAJ%5F471%5F732k.ped.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE23nnn/GSE23636/suppl/GSE23636%5FAJ%5F471%5F732k%5Freadme.txt.gz
gunzip * 

# basic QC
~/plink --file /data/scratch/hmy117/aj_ref/GSE23636_AJ_471_732k \
--make-bed \
--chr 1-22 \
--geno 0.1 \
--hwe 1e-10 \
--mac 1 \
--mind 0.1 \
--out aj_qc_ref

# make a snp list 
cut -f2 aj_qc_ref.bim > aj_snp_list

cd /data/scratch/hmy117/aj_ref/

# liftover to hg38
awk '{print "chr"$1,$4-1,$4,$2}' aj_qc_ref.bim > hg18_bedfile

# run liftover
/data/Wolfson-UKBB-Dobson/liftover/liftOver \
hg18_bedfile \
/data/Wolfson-UKBB-Dobson/liftover/hg18ToHg38.over.chain.gz \
hg38_bedfile \
unmapped

awk '{print $4,$3}' hg38_bedfile > hg38_snp_positions

# Update SNP positions
~/plink --bfile aj_qc_ref \
--update-map hg38_snp_positions \
--make-bed \
--out aj_qc_ref_hg38

````

#### Download HGDP reference data
````unix
qsub /data/home/hmy117/ADAMS_severity/scripts/download_hgdp_1kg.sh

# metadata
~/google-cloud-sdk/bin/gsutil cp \
gs://gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg/data_intersection/hgdp_1kg_sample_info.unrelateds.pca_outliers_removed.with_project.tsv \
/data/scratch/hmy117/hgdp_1kg_genomes/

wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg/metadata_and_qc/gnomad_meta_v1.tsv
````

#### Filter HGDP to variants in AJ dataset
````unix
qsub /data/home/hmy117/ADAMS_severity/scripts/filter_hgdp_to_aj_vars.sh
````

#### Merge HGDP files
````unix 
grep 71452575 /data/scratch/hmy117/hgdp_1kg_genomes/filtered_1kg_hgdp_aj_1kg_chr2.bim

cd /data/scratch/hmy117/hgdp_1kg_genomes/
for i in {2..22};
  do
    echo filtered_1kg_hgdp_aj_1kg_chr$i >> merge_filelist
  done


# try merge
~/plink --bfile filtered_1kg_hgdp_aj_1kg_chr1 \
--merge-list merge_filelist \
--make-bed \
--biallelic-only \
--out combined_hgdp_1kg_aj

# remove duplicates & do QC
for i in {1..22};
  do
    ~/plink --bfile filtered_1kg_hgdp_aj_1kg_chr$i \
    --exclude combined_hgdp_1kg_aj-merge.missnp \
    --make-bed \
    --out filtered_1kg_hgdp_aj_1kg_nodups_chr$i \
    --maf 0.05 \
    --geno 0.01 \
    --hwe 1e-10
  done

# repeat merge 
rm merge_filelist
for i in {2..22};
  do
    echo filtered_1kg_hgdp_aj_1kg_nodups_chr$i >> merge_filelist
  done


# try merge
~/plink --bfile filtered_1kg_hgdp_aj_1kg_nodups_chr1 \
--merge-list merge_filelist \
--make-bed \
--out combined_hgdp_1kg_aj


````

#### Filter AJ variants to HGDP variants 
````unix 
~/plink --bfile /data/scratch/hmy117/aj_ref/aj_qc_ref_hg38 \
--extract combined_hgdp_1kg_aj.bim \
--out aj_for_merge \
--make-bed
````

#### Merge 
````unix 
# change to CPRA
~/plink2 --bfile combined_hgdp_1kg_aj \
--set-all-var-ids @:#:\$r\:\$a \
--out combined_hgdp_1kg_aj_cpra \
--make-bed 

~/plink2 --bfile aj_for_merge \
--set-all-var-ids @:#:\$r\:\$a \
--out aj_for_merge_cpra \
--make-bed 

# filter to perfect overlap 
~/plink --bfile combined_hgdp_1kg_aj_cpra \
--extract aj_for_merge_cpra.bim \
--make-bed \
--out hgdp_1kg_for_merge_qc

~/plink --bfile aj_for_merge_cpra \
--extract combined_hgdp_1kg_aj_cpra.bim \
--make-bed \
--out aj_for_merge_qc 

# merge
~/plink --bfile hgdp_1kg_for_merge_qc \
--bmerge aj_for_merge_qc \
--out hgdp_1kg_aj_combined \
--make-bed 

````

#### Combine with genotype data 
````unix 
cd /data/scratch/hmy117/hgdp_1kg_genomes/

#head /data/home/hmy117/ADAMS_severity/combined_adams_imputed.bim

# extract adams variants in ref panel 
~/plink2 --bfile hgdp_1kg_aj_combined \
--extract /data/scratch/hmy117/adams_imputed_severity_topmed/combined_adams_imputed_qc.bim \
--rm-dup exclude-all \
--out hgdp_1kg_aj_combined_cpra_in_adams \
--make-bed


# prune - changed r2 to 0.2 02-10
~/plink --bfile hgdp_1kg_aj_combined_cpra_in_adams \
--indep-pairwise 1000 100 0.01 \
--out snps_for_pca \
--make-bed

# pca 
~/plink2 --bfile hgdp_1kg_aj_combined_cpra_in_adams \
--out hgdp_kg_aj_pca \
--pca approx allele-wts 20 \
--freq \
--extract snps_for_pca.prune.in

# project ADAMS samples
~/plink2 --bfile /data/scratch/hmy117/adams_imputed_severity_topmed/combined_adams_imputed_qc \
--read-freq hgdp_kg_aj_pca.afreq \
--score hgdp_kg_aj_pca.eigenvec.allele 2 5 header-read no-mean-imputation variance-normalize list-variants \
--score-col-nums 6-25 \
--out adams_pcs

# project original dataset
~/plink2 --bfile hgdp_1kg_aj_combined_cpra_in_adams \
--read-freq hgdp_kg_aj_pca.afreq \
--extract adams_pcs.sscore.vars \
--score hgdp_kg_aj_pca.eigenvec.allele 2 5 header-read no-mean-imputation variance-normalize list-variants \
--score-col-nums 6-25 \
--out hgdp_kg_aj_pcs_rescored
````

### ADMIXTURE prep
````unix 
cd /data/scratch/hmy117/hgdp_1kg_genomes/
# filter
~/plink2 --bfile /data/scratch/hmy117/adams_imputed_severity_topmed/combined_adams_imputed_qc \
--extract snps_for_pca.prune.in \
--out /data/scratch/hmy117/adams_imputed_severity_topmed/combined_adams_imputed_qc_pruned \
--make-bed

# merge 
~/plink --bfile hgdp_1kg_aj_combined_cpra_in_adams \
--bmerge /data/scratch/hmy117/adams_imputed_severity_topmed/combined_adams_imputed_qc_pruned \
--make-bed \
--out /data/scratch/hmy117/adams_imputed_severity_topmed/combined_adams_imputed_qc_pruned_merge_for_admixture


````

#### Make admixture pop file 
````R 
library(tidyverse)
setwd("/data/scratch/hmy117/hgdp_1kg_genomes")

# read in data
fam = read_table("/data/scratch/hmy117/adams_imputed_severity_topmed/combined_adams_imputed_qc_pruned_merge_for_admixture.fam",col_names=F)
colnames(fam)[2]="IID"
colnames(fam)[1]="FID"
fam$fid_iid = paste0(fam$FID,"_",fam$IID)
aj_fam = read_table("aj_for_merge_cpra.fam",col_names=F)
colnames(aj_fam)[2]="IID"
colnames(aj_fam)[1]="FID"
aj_fam$fid_iid = paste0(aj_fam$FID,"_",aj_fam$IID)

meta = read_tsv("/data/scratch/hmy117/aj_ref/gnomad_meta_v1.tsv") %>%
  dplyr::select(1,2,hgdp_tgp_meta.Genetic.region)
colnames(meta) = c("FID","IID","superpop")
meta$fid_iid = paste0(meta$FID,"_",meta$IID)

# join with gnomad 
fam = fam %>% left_join(meta %>% dplyr::select(-fid_iid),by=c("FID","IID"))

# make unique fid_iid
fam = fam %>%
  mutate(pop = case_when(
    fid_iid %in% aj_fam$fid_iid ~ "ASJ",
    fid_iid %in% meta$fid_iid ~ superpop,
    .default = as.character("-")    
  ))
fam = fam %>% dplyr::select(pop)
write_tsv(fam,"/data/scratch/hmy117/adams_imputed_severity_topmed/combined_adams_imputed_qc_pruned_merge_for_admixture.pop",col_names=F)

````

#### Run admixture
````unix
cd /data/scratch/hmy117/hgdp_1kg_genomes/

# prune 
~/dist/admixture_linux-1.3.0/admixture \
/data/scratch/hmy117/adams_imputed_severity_topmed/combined_adams_imputed_qc_pruned_merge_for_admixture.bed 8 \
--supervised

````

#### Ancestry inference 
````R
library(tidyverse)
setwd("/data/scratch/hmy117/hgdp_1kg_genomes")

# read in data
adams = read_table("adams_pcs.sscore")

kg_hgdp = read_table("hgdp_kg_aj_pcs_rescored.sscore")
meta = read_tsv("/data/scratch/hmy117/aj_ref/gnomad_meta_v1.tsv") %>%
  dplyr::select(2,hgdp_tgp_meta.Population,hgdp_tgp_meta.Genetic.region)
colnames(meta) = c("IID","pop","superpop")
kg_hgdp = kg_hgdp %>%
left_join(meta,by="IID")
kg_hgdp = kg_hgdp %>%
mutate(superpop = ifelse(is.na(superpop),"ASJ",superpop))

# build RF on HGDP data
# filter out those with missing labels
kg_hgdp = kg_hgdp %>%
  filter(!is.na(superpop)) %>%
  dplyr::select(superpop,contains("PC"))
library(caret)

# Define the training control & tunegrid
trainControl = trainControl(method = "cv", number = 10)
tuneGrid = expand.grid(.mtry = c(2, 3, 4, 5))

rf_fit = train(superpop ~ .,
                      data=kg_hgdp,
                      method='rf',
                      metric='Accuracy',
                      trControl = trainControl, 
                      tuneGrid = tuneGrid)

# save rf
saveRDS(rf_fit,"rf_fit.rds")


# can reload here
# rf_fit = readRDS("rf_fit.rds")

# predict
adams$predicted_ancestry = predict(rf_fit,adams)

predicted_ancestry_confidence = predict(rf_fit,adams,type="prob")
predicted_ancestry_confidence$IID = adams$IID

predicted_ancestry_confidence = tibble(predicted_ancestry_confidence) %>%
	mutate(max_prob = pmax(AFR,AMR,CSA,EAS,EUR,OCE,MID,ASJ))

adams = adams %>% 
  mutate(max_prob = predicted_ancestry_confidence$max_prob) 

# save anc calls 
write_tsv(adams %>% dplyr::select(1,2,predicted_ancestry),"/data/home/hmy117/ADAMS_severity/outputs/ancestry_calls.tsv")

# compare with ADMIXTURE 
q = read_table("combined_adams_imputed_qc_pruned_merge_for_admixture.8.Q",col_names=F)
pop = read_table("/data/scratch/hmy117/adams_imputed_severity_topmed/combined_adams_imputed_qc_pruned_merge_for_admixture.pop",col_names=F) %>% 
filter(X1!="-") %>% 
distinct()
colnames(q) = pop$X1
fam = read_table("/data/scratch/hmy117/adams_imputed_severity_topmed/combined_adams_imputed_qc_pruned_merge_for_admixture.fam",col_names=F)
colnames(fam)[2]="IID"
q = q %>% 
mutate(IID = fam$IID) 
q = q %>% pivot_longer(-IID)
anc = read_tsv("/data/home/hmy117/ADAMS_severity/outputs/ancestry_calls.tsv")

q = q %>% inner_join(anc,by="IID") 

q = q %>% arrange(predicted_ancestry,value)
q$IID = factor(q$IID,levels=unique(q$IID),ordered=T)


````


## Kinship
````unix
cd /data/scratch/hmy117/adams_imputed_severity_topmed/
~/king -b combined_adams_imputed_qc.bed --duplicate
~/king -b combined_adams_imputed_qc.bed --related --degree 3
~/plink --bfile combined_adams_imputed_qc --missing --out missingness_report
````


# Phenotype
Cleaning of phenotypic data.  
````R 
library(tidyverse)

# read in data 
data = read_csv("/data/home/hmy117/ADAMS_severity/pheno/ADAMS_pheno.csv")
fam = read_table("/data/home/hmy117/ADAMS/genotypes/QMUL_Aug_23/outputs/ADAMS_geno_fid_iid.fam",col_names=F) %>% 
  dplyr::select(X2,X5) %>% 
  dplyr::rename("IID"=X2,"genetic_sex"=X5)

# light QC 

# remove people with all NAs
remove_na = function(dat, field){
  pre_n = nrow(dat)
  message("Showing removed observations")
  dat %>% filter(is.na(dat[[field]])) %>% print(n=1000)
  dat = dat %>% filter(!is.na(dat[[field]]))
  post_n = nrow(dat)
  message("N before filtering: ",pre_n)
  message("N after filtering: ",post_n)
  message("N remove: ",pre_n - post_n)
  dat
  
}

# remove missing oragenes
data = remove_na(data,"Oragene ID")

# filter to people in post-qc genetic dataset 
fam_file = read_table("/data/home/hmy117/ADAMS_severity/combined_adams_imputed.fam",col_names = F)

data = data %>% 
  filter(`Oragene ID` %in% fam_file$X2)


#############################################
# Clean EDSS 
#############################################

# visualise obvious data entry errors
data = data %>% mutate(obs_num = row_number())
og_data = data 

data %>% 
  filter(latest_edss == 0) %>% 
  dplyr::count(mobility_clean)

# edit EDSS where there is a clear mismatch with self-reported mobility to NA
data = data %>% 
  mutate(latest_edss = ifelse(latest_edss == 0 & !is.na(mobility_clean) & !(mobility_clean %in% 
  c("I can walk for 20 minutes without a stick or frame",
  "I can run or walk as far as I want without any walking aids",
  "I can walk for 5-10 minutes without a stick or frame")), NA,latest_edss)) %>% 
    arrange(latest_edss) %>%
  mutate(latest_edss = ifelse(latest_edss == 0 & !is.na(msis_physical_normalised) & msis_physical_normalised >=20, NA,latest_edss))  
  
og_data %>% filter(!is.na(latest_edss)) %>% nrow()
og_data %>% filter(!is.na(latest_edss)) %>% distinct(Token) %>% nrow()

a = data %>% filter(is.na(latest_edss)) %>% nrow()
b = og_data %>% filter(is.na(latest_edss)) %>% nrow()
a - b


#############################################
# Pheno correlations pre-cleaning 
#############################################

pheno_dat = data %>% 
dplyr::select(age_at_recruitment,latest_edss,eq5d_vas,msis_physical_normalised,msis_psych_normalised)
colnames(pheno_dat) = c("Age","EDSS","EQ5D","MSIS-Phys","MSIS-Psych")
cor.test(pheno_dat$EDSS,pheno_dat$`MSIS-Phys`,method="spearman")
corrmat = cor(pheno_dat,use="pairwise.complete.obs")

corrmat
# relevel mobility
data$mobility_clean = factor(data$mobility_clean,ordered=T,levels = 
c(
  "I am unable to mobilise out of bed",
  "I use a wheelchair",
  "I use two walking sticks or a frame",
  "I use a walking stick",
  "I can walk for 5-10 minutes without a stick or frame",
  "I can walk for 20 minutes without a stick or frame",
  "I can run or walk as far as I want without any walking aids"
))

library(ggcorrplot)  
png("/data/home/hmy117/ADAMS_severity/plots/outcome_corr.png",res=900,units="in",width=6,height=6)
ggcorrplot(corrmat, outline.color = "white",lab=T,legend.title="Correlation coefficient")
dev.off()

png("/data/home/hmy117/ADAMS_severity/plots/mobility_vs_edss_and_msis.png",res=900,units="in",width=10,height=6)
ggplot(data,aes(latest_edss,msis_physical_normalised,fill=mobility_clean))+
  geom_point(shape=21,size=3)+
  theme_bw()+
  scale_fill_brewer(palette="Set1")+
  labs(x="EDSS",y="MSIS",fill="Self-reported mobility")+
  scale_x_continuous(breaks = seq(0,10,by=1))
dev.off()

# pairs plot 
png("/data/home/hmy117/ADAMS_severity/plots/outcome_corr_full.png",res=900,units="in",width=10,height=10)
GGally::ggpairs(pheno_dat)+
theme_bw()
dev.off()

#############################################
# Impute missing EDSS 
#############################################

# add ancestry 
ancestry_calls = read_tsv("/data/home/hmy117/ADAMS_severity/outputs/ancestry_calls.tsv",col_types=c("ddc"))
data$IID = data$`Oragene ID`
data = data %>% left_join(ancestry_calls,by="IID")
data = data %>% mutate(simple_ancestry = ifelse(predicted_ancestry %in% c("EUR","CSA","AFR"),predicted_ancestry,"other"))

# add indicators for na
data = data %>% 
  mutate(has_edss = ifelse(!is.na(latest_edss),"Yes","No")) %>%
  mutate(has_msis = ifelse(!is.na(msis_physical_normalised),"Yes","No")) %>%
  mutate(has_eq5d = ifelse(!is.na(eq5d_vas),"Yes","No"))

# repeat with MSIS bin
data$edss_bin = Hmisc::cut2(data$latest_edss,cuts = seq(0,10,by=2))
data$msis_bin = Hmisc::cut2(data$msis_physical_normalised,cuts = seq(0,100,by=25))
data$eq5d_bin = Hmisc::cut2(data$eq5d_vas,cuts = seq(0,100,by=10))


# get median edss per mobility bin 
median_edss_per_mobility_bin = data %>% 
  group_by(mobility_clean) %>% 
  summarise(imputed_edss = median(latest_edss,na.rm=T)) %>% 
  filter(!is.na(mobility_clean))

median_edss_per_mobility_bin_full = data %>% 
  group_by(mobility_clean) %>% 
  summarise(imputed_edss = median(latest_edss,na.rm=T),
  iqr_edss = IQR(latest_edss,na.rm=T) ) %>% 
  filter(!is.na(mobility_clean))
write_csv(median_edss_per_mobility_bin_full,"/data/home/hmy117/ADAMS_severity/outputs/mobility_edss_imputation.csv")


data = data %>%
  left_join(median_edss_per_mobility_bin,by="mobility_clean")

data = data %>%
  mutate(edss_source = ifelse(is.na(latest_edss),"Imputed","Observed")) %>%
  mutate(latest_edss = ifelse(is.na(latest_edss),imputed_edss,latest_edss))

# count 
data %>% filter(edss_source == "Imputed" & !is.na(latest_edss)) %>% nrow()
data %>% filter(edss_source == "Imputed" & !is.na(latest_edss)) %>% distinct(IID) %>% nrow()

# MSIS bin
median_edss_per_msis_bin = data %>% 
  group_by(msis_bin) %>% 
  summarise(imputed_edss = median(latest_edss,na.rm=T)) %>% 
  filter(!is.na(msis_bin))
data %>% filter(!is.na(msis_physical_normalised) & is.na(latest_edss)) %>% nrow()

median_edss_per_msis_bin_full = data %>% 
  group_by(msis_bin) %>% 
  summarise(imputed_edss = median(latest_edss,na.rm=T),
  iqr_edss = IQR(latest_edss,na.rm=T) ) %>% 
  filter(!is.na(msis_bin))
write_csv(median_edss_per_msis_bin_full,"/data/home/hmy117/ADAMS_severity/outputs/msis_edss_imputation.csv")

data = data %>%
  dplyr::select(-imputed_edss) %>%
  left_join(median_edss_per_msis_bin,by="msis_bin")
data = data %>%
  mutate(latest_edss = ifelse(is.na(latest_edss),imputed_edss,latest_edss))
median_edss_per_msis_bin

data %>% filter(!is.na(eq5d_vas) & is.na(latest_edss)) %>% nrow()

# approximate age at edss if imputed  using age at recruitment 
data = data %>% 
  mutate(ageatedss = ifelse(edss_source=="Imputed",age_at_recruitment,age_at_latest_edss))

# N 
data %>% filter(!is.na(latest_edss)) %>% nrow()
data %>% filter(!is.na(latest_edss)) %>% distinct(IID) %>% nrow()

#############################################
# Exclude duplicates based on sample with no EDSS 
#############################################


dups = read_table("/data/scratch/hmy117/adams_imputed_severity_topmed/king.con",col_types = cols(.default = "d"))

just_edss =  data %>% dplyr::select(IID,latest_edss)
dups = dups %>% 
  dplyr::rename("IID" = ID1) %>% left_join(just_edss,by="IID") %>% 
  dplyr::rename("ID1" = IID) %>%
  dplyr::rename("IID" = ID2) %>% left_join(just_edss,by="IID") %>%  
  dplyr::rename("ID2" = IID) %>% 
  mutate(person_to_discard = ifelse(is.na(latest_edss.x),ID1,ID2)) 

# save to file 
dups = dups %>% 
  mutate(FID = person_to_discard,IID = person_to_discard) %>%
  dplyr::select(FID,IID) %>% 
  distinct()
write_tsv(dups,"/data/scratch/hmy117/adams_imputed_severity_topmed/dups_to_discard.tsv")


# exclude dups from phenotype data 
data = data %>% filter(!IID %in% dups$IID)

# arrange by missingness for full duplicate (oragene duplicates)
doubles = data %>% dplyr::count(IID) %>% filter(n>1)

data = data %>% 
  arrange(latest_edss) %>%
  distinct(IID,.keep_all=TRUE)

# N again 
data %>% filter(!is.na(latest_edss)) %>% nrow()
data %>% filter(!is.na(latest_edss)) %>% distinct(IID) %>% nrow()

#############################################
# Calculate ARMSS 
#############################################

# ARMSS 
# install.packages("https://cran.r-project.org/src/contrib/Archive/ms.sev/ms.sev_1.0.4.tar.gz",repos=NULL,type="source")

library(ms.sev)

# calculate ARMSS score
data = data %>%
  dplyr::rename("edss" = latest_edss)

# if age at edss is <18 truncate to 18
# if age at edss is >75 to 75
data %>% filter(!is.na(edss) & !is.na(ageatedss) & (ageatedss < 18 | ageatedss > 75) ) %>% nrow()
data = data %>%
  mutate(ageatedss = case_when(
    ageatedss < 18 ~ 18,
    ageatedss > 75 ~ 75,
    ageatedss >= 18 | ageatedss <= 75 ~ ageatedss
  ))
## calculate ARMSS
armss = ms.sev::global_armss(data)
data = armss$data

# calculate MSSS
# set age at dx to missing if <18 or if age at dx was after recruitment

data = data %>% mutate(age_at_dx = ifelse(age_at_dx<18,NA,age_at_dx))
data = data %>% mutate(age_at_dx = ifelse(age_at_dx > age_at_recruitment,NA,age_at_dx))
data = data %>% mutate(age_at_dx = ifelse(age_at_dx > ageatedss | is.na(ageatedss),NA,age_at_dx))


data = data %>% mutate(dd = ageatedss - age_at_dx)
msss = ms.sev::global_msss(data)
data = msss$data

data %>% filter(!is.na(edss) & !is.na(ageatedss)) %>% nrow()
data %>% filter(!is.na(gARMSS)) %>% nrow()
data %>% filter(!is.na(uGMSSS)) %>% nrow()
data %>% filter(!is.na(edss) & !is.na(ageatedss) & !is.na(age_at_dx)) %>% nrow()

#############################################
# Phenotype descriptives 
#############################################
# add genetic sex
data = data %>% left_join(fam %>% mutate(IID = as.numeric(IID)),by="IID")
table(data$genetic_sex,data$Sex)

data = data %>%
  mutate(Sex = ifelse(genetic_sex==2,"F","M"))

# overall 
data %>% nrow()
data %>% dplyr::count(Sex) %>% mutate(prop = n/sum(n))
data %>% summarise_at(.vars = "age_at_dx",.funs = c("median","IQR"),na.rm=T)
data %>% summarise_at(.vars = "age_at_recruitment",.funs = c("median","IQR"),na.rm=T)
data %>% filter(!is.na(subtype_clean)) %>% dplyr::count(subtype_clean) %>% mutate(prop = n/sum(n))
data %>% dplyr::count(Site) %>% mutate(prop = n/sum(n))
data %>% dplyr::count(simple_ancestry) %>% mutate(prop = n/sum(n))
data %>% summarise_at(.vars = "edss",.funs = c("median","IQR"),na.rm=T)
data %>% filter(!is.na(edss)) %>% dplyr::count(edss_source) 
data %>% group_by(edss_source) %>% summarise_at(.vars = "edss",.funs = c("median","IQR"),na.rm=T)
data %>% summarise_at(.vars = "eq5d_vas",.funs = c("median","IQR"),na.rm=T)
data %>% summarise_at(.vars = "msis_physical_normalised",.funs = c("median","IQR"),na.rm=T)
data %>% dplyr::count(predicted_ancestry)

# contrast ancestry with ethnicity
data %>%
  group_by(simple_ancestry) %>% 
  dplyr::count(ethnicity_clean) %>% 
  mutate(prop = n/sum(n))

# descriptives 
library(compareGroups)

# all cohort by ancestry 
tbl = compareGroups::compareGroups(data = data,
                                   simple_ancestry ~ age_at_dx + age_at_recruitment + Site +
                                     subtype_clean + Sex + edss + gARMSS + uGMSSS + msis_physical_normalised + 
                                     msis_psych_normalised + 
                                     mobility_clean +
                                     neurologist_clean +
                                     current_dmt_clean + ethnicity_clean + has_edss + has_msis + has_eq5d + edss_source,
                                   method = c(2,2,3,3,3,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3))
tbl2 = compareGroups::createTable(tbl)

export2csv(tbl2,"/data/home/hmy117/ADAMS_severity/outputs/demographics_table_by_ancestry.csv")

# by recruitment site
tbl = compareGroups::compareGroups(data = data,
                                   Site ~ age_at_dx + age_at_recruitment  + 
                                     subtype_clean + Sex + edss + gARMSS + uGMSSS + msis_physical_normalised + 
                                     msis_psych_normalised + 
                                     simple_ancestry + ethnicity_clean + has_edss + has_msis + has_eq5d+ edss_source,
                                   method = c(2,2,3,3,2,2,2,2,2,3,3,3,3,3,3))
tbl2 = compareGroups::createTable(tbl)

export2csv(tbl2,"/data/home/hmy117/ADAMS_severity/outputs/demographics_table_by_site.csv")


##############################################e 
# Add batch 
############################################### 

batch = read_csv("/data/home/hmy117/ADAMS_severity/modified_sample_sheet.csv",skip=10) %>% 
  dplyr::select(1,batch)
colnames(batch)[1] = "Oragene ID"  
batch$`Oragene ID` = as.numeric(batch$`Oragene ID`)
batch = batch %>% group_by(`Oragene ID`) %>% mutate(order = row_number()) %>% slice_max(order)
data = data %>% left_join(batch,by="Oragene ID")

###############################################
# Make pheno and covar files for genetics
###############################################

setwd("/data/home/hmy117/ADAMS_severity/")
## select pheno cols
pheno_file = data %>% 
  dplyr::select(`Oragene ID`,gARMSS,edss,msis_physical_normalised,eq5d_vas,age_at_dx,uGMSSS) %>%
  dplyr::rename("FID" = `Oragene ID`) %>%
  mutate("IID" = FID) %>%
  dplyr::select(FID,IID,everything()) %>%
  distinct(IID,.keep_all = T)
write_tsv(pheno_file,"./pheno/adams_pheno.tsv",col_names = T)

## make covar file
covar_file = data %>% 
  dplyr::select(`Oragene ID`,ageatedss,Sex,ethnicity_clean,Site,batch) %>%
  mutate(agesq = ageatedss^2) %>%
  dplyr::rename("FID" = `Oragene ID`) %>%
  mutate("IID" = FID) %>%
  dplyr::select(FID,IID,everything()) %>%
  distinct(IID,.keep_all = T)

write_tsv(covar_file,"./pheno/adams_covars.tsv",col_names = T)
saveRDS(data,"./outputs/cleaned_pheno_data.rds")

###############################################
# ARMSS validation
###############################################

# remake data 
library(tidyverse)
setwd("/data/home/hmy117/ADAMS_severity/")
data = readRDS("./outputs/cleaned_pheno_data.rds")

# derive correlations between phenotypes to get N effective tests 
pheno_dat = data %>% 
  dplyr::select(age_at_dx,edss,eq5d_vas,msis_physical_normalised,gARMSS,subtype_clean)
colnames(pheno_dat) = c("Age at dx","EDSS","EQ5D","MSIS-Phys","gARMSS","Subtype")
pheno_dat = pheno_dat %>% mutate(Subtype = case_when(
  Subtype == "PPMS" ~ 1,
  Subtype == "SPMS" ~ 2,
  Subtype == "RRMS" ~ 3  
))
prcomp(na.omit(pheno_dat))
data.frame(
  sd = prcomp(na.omit(pheno_dat))$sd
) %>%
mutate(cum_pve = cumsum(sd)/sum(sd))

  cor.test(pheno_dat$EDSS,pheno_dat$`MSIS-Phys`,method="spearman")
  corrmat = cor(pheno_dat,use="pairwise.complete.obs")



# add global ancestry proportions 
q = read_table("/data/scratch/hmy117/hgdp_1kg_genomes/combined_adams_imputed_qc_pruned_merge_for_admixture.8.Q",col_names=F)
pop = read_table("/data/scratch/hmy117/adams_imputed_severity_topmed/combined_adams_imputed_qc_pruned_merge_for_admixture.pop",col_names=F) %>% 
filter(X1!="-") %>% 
distinct()
colnames(q) = pop$X1
fam = read_table("/data/scratch/hmy117/adams_imputed_severity_topmed/combined_adams_imputed_qc_pruned_merge_for_admixture.fam",col_names=F)
colnames(fam)[2]="IID"
q = q %>% 
mutate(IID = fam$IID) 

data = data %>% inner_join(q %>% mutate(IID = as.numeric(IID)),by="IID")

# plot global ancestry
plot_dat = q %>% 
  arrange(desc(AFR),desc(AMR),desc(ASJ),desc(CSA),desc(EAS),desc(EUR),desc(MID),desc(OCE)) %>%
  pivot_longer(-IID) %>% 
  filter(IID %in% data$IID)  %>% 
  left_join(data %>% dplyr::select(IID,predicted_ancestry) %>% mutate(IID = as.character(IID)),by="IID") %>%
  arrange(predicted_ancestry) 

labs = plot_dat %>% distinct(IID,.keep_all=T) %>% 
mutate(row = row_number() ) %>%
group_by(predicted_ancestry) %>% 
summarise(min = min(row),max = max(row))

plot_dat$IID = factor(plot_dat$IID,levels = unique(plot_dat$IID),ordered=T)
png("/data/home/hmy117/ADAMS_severity/plots/ancestry_proportions.png",res=900,units="in",width=12,height=4)
ggplot(plot_dat,aes(IID,value,fill=name))+
geom_col(color="black",lwd=0.05)+
scale_fill_brewer(palette="Accent")+
theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.position="top")+
labs(fill="Global ancestry proportion",x="",y="Global ancestry proportion")+
guides(fill = guide_legend(nrow = 1))+
geom_segment(
    inherit.aes = FALSE, 
    data = labs,
    mapping = aes(x = min, xend = max, y = 1.05, yend = 1.05, color = predicted_ancestry),
    size = 2,show.legend=F
) +
scale_color_brewer(palette="Accent")

dev.off()

# compare categorical ancestry and global ancestry 
plot_dat %>% 
  group_by(predicted_ancestry,name) %>% 
  summarise(median(value))

plot_dat %>% 
  filter(predicted_ancestry == "AFR") %>%
  pivot_wider(names_from = name, values_from = value) %>%
  dplyr::count(EUR>0.1)  %>% 
  mutate(prop = n/sum(n), total = sum(n))

# ancestry vs ethnicity
png("/data/home/hmy117/ADAMS_severity/plots/ancestry_vs_ethnicity.png",res=900,units="in",width=6,height=4)
ggplot(data %>% filter(!is.na(ethnicity_clean)),aes(predicted_ancestry,fill=ethnicity_clean))+
  geom_bar(position="fill",color="black")+
  theme_bw()+
  labs(x="Genetic ancestry",y="Proportion",fill="Self-reported ethnicity")+
  scale_fill_brewer(palette="Accent")
dev.off()


# bootstrap regression models 
bootstrap_linreg = function(anc){
  all_coefs = list()
  for(i in c(1:1000)){
    dat = data %>% filter(!is.na(edss) & predicted_ancestry == anc)
    dat = sample_n(dat,size = nrow(dat),replace=T)
    coefs = data.frame(summary(
      lm(data = dat,edss ~ ageatedss)
    )$coefficients)
    all_coefs[[i]] = coefs$Estimate[2]
  }
  return(unlist(all_coefs))
}
ancestries = c("CSA","AFR","EUR","MID","ASJ")
res_df = data.frame(sapply(ancestries,bootstrap_linreg))
res_df = res_df %>% 
  pivot_longer(cols = everything()) %>%
  group_by(name) %>%
  summarise(mean= mean(value), 
            lower_ci = quantile(value,0.025),
            upper_ci = quantile(value,0.975))

# calculate estimated edss increment for 10 year increase in age 
res_df = res_df %>% 
  mutate(edss_increment = paste0(
    "+",
    round(mean*10,2)," (95% CI ",
    round(lower_ci*10,2)," - ",
    round(upper_ci*10,2),")" )
    ) 
res_df$simple_ancestry=res_df$name

png("./plots/garmss_age_linreg.png",res=900,units="in",width=12,height=5)
ggplot(data %>% filter(simple_ancestry %in% c("CSA","AFR","EUR")),
aes(ageatedss,edss,fill=gARMSS))+
geom_point(alpha=0.5,size=3,shape=21)+
geom_text(data = res_df%>% filter(simple_ancestry %in% c("CSA","AFR","EUR")),aes(x=50,y=10,label = paste0("Delta EDSS per 10y increase:\n",edss_increment),fill=NA))+
facet_wrap(~simple_ancestry)+
geom_smooth(method="lm",se=F,color="red",alpha=0.2,linetype="dashed",lwd=0.5)+
theme_bw()+
scale_y_continuous(breaks = seq(0,10,by=0.5))+
scale_x_continuous(breaks = seq(0,100,by=10))+
scale_fill_viridis_c(option="plasma",breaks = seq(0,10,by=1),limits = c(0,11))+
labs(x="Age at EDSS",y="EDSS",fill="ARMSS")
dev.off()


png("./plots/age_vs_edss_vs_armss.png",res=900,units="in",width=5,height=3)
ggplot(data,
aes(ageatedss,edss,fill=gARMSS))+
geom_point(alpha=0.5,size=3,shape=21)+
theme_bw()+
scale_y_continuous(breaks = seq(0,10,by=0.5))+
scale_x_continuous(breaks = seq(0,100,by=10))+
scale_fill_viridis_c(option="plasma",breaks = seq(0,10,by=1),limits = c(0,11))+
labs(x="Age at EDSS",y="EDSS",fill="ARMSS")
dev.off()

# repeat, adjusting for sex 
# bootstrap regression models 
bootstrap_linreg = function(anc){
  all_coefs = list()
  for(i in c(1:1000)){
    dat = data %>% filter(!is.na(edss) & predicted_ancestry == anc)
    dat = sample_n(dat,size = nrow(dat),replace=T)
    coefs = data.frame(summary(
      lm(data = dat,edss ~ Sex + ageatedss)
    )$coefficients)
    all_coefs[[i]] = coefs$Estimate[3]
  }
  return(unlist(all_coefs))
}
ancestries = c("CSA","AFR","EUR")
res_df = data.frame(sapply(ancestries,bootstrap_linreg))
res_df = res_df %>% 
  pivot_longer(cols = everything()) %>%
  group_by(name) %>%
  summarise(mean= mean(value), 
            lower_ci = quantile(value,0.025),
            upper_ci = quantile(value,0.975))


# calculate estimated edss increment for 10 year increase in age 
res_df = res_df %>% 
  mutate(edss_increment = paste0(
    "+",
    round(mean*10,2)," (95% CI ",
    round(lower_ci*10,2)," - ",
    round(upper_ci*10,2),")" )
    ) 
res_df$simple_ancestry=res_df$name

png("./plots/garmss_age_linreg_sensitivity_adjusted_for_sex.png",res=900,units="in",width=12,height=5)
ggplot(data %>% filter(simple_ancestry %in% c("CSA","AFR","EUR")),
aes(ageatedss,edss,fill=gARMSS))+
geom_point(alpha=0.5,size=3,shape=21)+
geom_text(data = res_df%>% filter(simple_ancestry %in% c("CSA","AFR","EUR")),aes(x=50,y=10,label = paste0("Delta EDSS per 10y increase:\n",edss_increment),fill=NA))+
facet_wrap(~simple_ancestry)+
geom_smooth(method="lm",se=F,color="red",alpha=0.2,linetype="dashed",lwd=0.5)+
theme_bw()+
scale_y_continuous(breaks = seq(0,10,by=0.5))+
scale_x_continuous(breaks = seq(0,100,by=10))+
scale_fill_viridis_c(option="plasma",breaks = seq(0,10,by=1),limits = c(0,11))+
labs(x="Age at EDSS",y="EDSS",fill="ARMSS")
dev.off()


###############################################
# Remove rare ancestral groups 
###############################################

prevalent_ancestries = data %>% dplyr::count(predicted_ancestry) %>% filter(n>20)
data = data %>% filter(predicted_ancestry %in% prevalent_ancestries$predicted_ancestry)

# impute missing armss
library(mice)
data = data %>%
  mutate(raw_edss = ifelse(edss_source=="Observed",edss,NA))

predictor_matrix = matrix(c(0,1,1,1,0,1,1,1,0),nrow=3)
rownames(predictor_matrix) = c("msis_physical_normalised","mobility_clean","raw_edss")
colnames(predictor_matrix) = c("msis_physical_normalised","mobility_clean","raw_edss")
just_edss = data %>% dplyr::select(raw_edss,msis_physical_normalised,mobility_clean)

# impute with n_iter iterations
n_iter = 100
imputed_data = mice(just_edss,m=n_iter,maxit=50,meth='pmm',seed=500,predictorMatrix = predictor_matrix)

# join with main data 
saveRDS(imputed_data,"./outputs/imputed_edss_scores.rds")
# imputed_data = readRDS("./outputs/imputed_edss_scores.rds")
imputed_edss = imputed_data$imp$raw_edss %>% mutate(row_id = rownames(imputed_data$imp$raw_edss))
colnames(imputed_edss)[c(1:n_iter)] = paste0("imputed_edss",seq(1:n_iter))
data$row_id = seq(1:nrow(data))
imputed_edss$row_id = as.numeric(imputed_edss$row_id)

# join 
data = data %>% left_join(imputed_edss,by="row_id")

# remove those with missing mobility & MSIS
model_dat = data %>%
 filter(!(
  is.na(mobility_clean) & is.na(raw_edss) & is.na(msis_physical_normalised)
 ))

# compare imputed vs observed values 
imputed_vs_observed = model_dat %>%
  dplyr::select(row_id,contains("imputed_edss")) %>%
  dplyr::select(-imputed_edss) %>%
  pivot_longer(cols = -c(1,2,3)) %>% 
  filter(!is.na(value)) %>%
  group_by(row_id) %>% 
  summarise(median_imp = median(value)) %>%
  full_join(
    model_dat %>% dplyr::select(mobility_clean,msis_physical_normalised,row_id,raw_edss),
    by="row_id"
  ) %>% 
  mutate(edss_source = ifelse(!is.na(raw_edss),"Observed","Imputed")) %>%
  mutate(edss_for_plot = ifelse(!is.na(raw_edss),raw_edss,median_imp))

png("/data/home/hmy117/ADAMS_severity/plots/multiple_imputation.png",res=900,units="in",width=10,height=4)

 ggplot(imputed_vs_observed %>% filter(!is.na(mobility_clean)),aes(msis_physical_normalised,edss_for_plot,fill=mobility_clean))+
  geom_point(shape=21,size=3)+
  facet_wrap(~edss_source)+
  theme_bw()+
  labs(x="MSIS Physical",y="EDSS",fill="Mobility")+
  scale_fill_brewer(palette="Set1")
dev.off()


# fill in observed values 
model_res = list()
for(i in c(1:n_iter)){
  model_dat = model_dat %>%
  mutate(model_edss = ifelse(is.na(raw_edss),.data[[paste0("imputed_edss",i)]],raw_edss))

  model_dat$edss = model_dat$model_edss

  # get armss 
  model_dat = ms.sev::global_armss(model_dat)$data
  model_dat = model_dat %>% filter(!is.na(gARMSS))
  
  # fit model 
  model_dat$armss_norm = RNOmni::RankNorm(model_dat$gARMSS)
  
  model_dat$predicted_ancestry = relevel(factor(model_dat$predicted_ancestry),ref="EUR")
  model_dat$Site = relevel(factor(model_dat$Site),ref="Website")

  armss_model = glm(data = model_dat, armss_norm ~ ageatedss + Site + Sex + predicted_ancestry) 
  model_res[[i]] = armss_model
}

pool(model_res) %>% summary() %>% 
  mutate(lower_ci = estimate - 1.96 * std.error) %>%
  mutate(upper_ci = estimate + 1.96 * std.error) %>%
  write_csv(file = "/data/home/hmy117/ADAMS_severity/outputs/imputed_model_mice.csv")


###############################################
# Ancestry vs age at dx 
###############################################

# plots
counts = data %>% dplyr::count(predicted_ancestry)
medians = data %>% group_by(predicted_ancestry) %>% summarise(median_age = median(age_at_dx,na.rm=T)) 

p1 = ggplot(data,aes(predicted_ancestry,age_at_dx,fill=predicted_ancestry))+
geom_violin(alpha=0.5)+
geom_boxplot(width=0.1,alpha=0.3)+
geom_jitter(width=0.1,alpha=0.1)+
theme_bw()+
labs(x="Ancestry",y="Age at diagnosis",fill="Ancestry")+
scale_fill_brewer(palette="Set1")+
geom_text(data = medians,aes(y=85,label = format(round(median_age,1),nsmall=1)))+
theme(legend.position="none")
png("/data/home/hmy117/ADAMS_severity/plots/age_at_dx.png",res=900,units="in",width=6,height=4)
p1
dev.off()

# repeat stratified by site 
counts = data %>% dplyr::count(predicted_ancestry,Site) %>% mutate(Site = ifelse(Site == "Other","Clinical site",Site))
medians = data %>% group_by(predicted_ancestry,Site) %>% summarise(median_age = median(age_at_dx,na.rm=T)) %>% mutate(Site = ifelse(Site == "Other","Clinical site",Site))

p1 = ggplot(data %>% filter(Site != "UKMSR")%>% mutate(Site = ifelse(Site == "Other","Clinical site",Site)),aes(predicted_ancestry,age_at_dx,fill=predicted_ancestry))+
geom_violin(alpha=0.5)+
facet_wrap(~Site)+
geom_boxplot(width=0.1,alpha=0.3)+
geom_jitter(width=0.1,alpha=0.1)+
theme_bw()+
labs(x="Ancestry",y="Age at diagnosis",fill="Ancestry")+
scale_fill_brewer(palette="Set1")+
geom_text(data = medians %>% filter(Site != "UKMSR"),aes(y=85,label = format(round(median_age,1),nsmall=1)))+
theme(legend.position="none")
png("/data/home/hmy117/ADAMS_severity/plots/age_at_dx_by_site.png",res=900,units="in",width=8,height=4)
p1
dev.off()

# subtype 
p1 = 
ggplot(data %>% filter(!is.na(subtype_clean)),aes(predicted_ancestry,fill=subtype_clean))+
geom_bar(position="fill",color="black")+
theme_bw()+
labs(x="Ancestry",y="Proportion",fill="Ancestry")+
scale_fill_brewer(palette="Accent")
png("/data/home/hmy117/ADAMS_severity/plots/subtypes.png",res=900,units="in",width=6,height=4)
p1
dev.off()



# full model for age at dx 
model_dat = data %>% filter(!is.na(age_at_dx))
model_dat %>% nrow()
model_dat %>% dplyr::count(predicted_ancestry)
model_dat$age_norm = RNOmni::RankNorm(model_dat$age_at_dx)
model_dat$predicted_ancestry = relevel(factor(model_dat$predicted_ancestry),ref="EUR")
glm(data = model_dat, age_norm ~ predicted_ancestry) %>% broom::tidy() 
glm(data = model_dat, age_norm ~ Sex + subtype_clean + predicted_ancestry) %>% broom::tidy() %>% arrange(p.value) 


# simple model bootstrap 
all_coefs = list()
for(i in c(1:10000)){
  dat = data %>% filter(!is.na(age_at_dx))
  dat = sample_n(dat,size = nrow(dat),replace=T)
  dat$predicted_ancestry = relevel(factor(dat$predicted_ancestry),ref="EUR")
  coefs = lm(data = dat,age_at_dx ~ Sex  + predicted_ancestry) %>% broom::tidy()
  all_coefs[[i]] = coefs
}
all_coefs = do.call("bind_rows",all_coefs)

# main model 
data$predicted_ancestry = relevel(factor(data$predicted_ancestry),ref="EUR")
main_model = lm(data = data,age_at_dx ~ Sex  + predicted_ancestry) %>% broom::tidy()

counts = data %>% filter(!is.na(age_at_dx)) %>% dplyr::count(predicted_ancestry)
model_res = all_coefs %>%
  group_by(term) %>% 
  summarise(sd_est = sd(estimate), lower = quantile(estimate,0.025), upper = quantile(estimate,0.975)) %>% 
  left_join(main_model %>%  dplyr::select(term,estimate) %>% dplyr::rename("beta" = estimate),by="term") %>%
  mutate(z = beta/sd_est) %>% 
  mutate(p = 2 * (1 - pnorm(abs(z)))) %>% 
  dplyr::select(term,beta,lower,upper,p)
write_csv(model_res,"./outputs/age_dx_model.csv")

# repeat with just abundant groups
all_coefs = list()
for(i in c(1:10000)){
  dat = data %>% filter(!is.na(age_at_dx) & predicted_ancestry %in% c("EUR","CSA","AFR"))
  dat = sample_n(dat,size = nrow(dat),replace=T)
  dat$predicted_ancestry = relevel(factor(dat$predicted_ancestry),ref="EUR")
  coefs = lm(data = dat,age_at_dx ~ Sex  + predicted_ancestry) %>% broom::tidy()
  all_coefs[[i]] = coefs
}
all_coefs = do.call("bind_rows",all_coefs)

# main model just major ancestries
model_dat = data %>% filter(!is.na(age_at_dx) & predicted_ancestry %in% c("EUR","CSA","AFR"))
model_dat$predicted_ancestry = relevel(factor(model_dat$predicted_ancestry),ref="EUR")
main_model = lm(data = model_dat,age_at_dx ~ Sex  + predicted_ancestry) %>% broom::tidy()

counts = data %>% filter(!is.na(age_at_dx)) %>% dplyr::count(predicted_ancestry)
model_res = all_coefs %>%
  group_by(term) %>% 
  summarise(sd_est = sd(estimate), lower = quantile(estimate,0.025), upper = quantile(estimate,0.975)) %>% 
  left_join(main_model %>%  dplyr::select(term,estimate) %>% dplyr::rename("beta" = estimate),by="term") %>%
  mutate(z = beta/sd_est) %>% 
  mutate(p = 2 * (1 - pnorm(abs(z)))) %>% 
  dplyr::select(term,beta,lower,upper,p)
write_csv(model_res,"./outputs/age_dx_model_just_abundant.csv")

# adjust for subtype 
data$predicted_ancestry = relevel(factor(data$predicted_ancestry),ref="EUR")
data$subtype_clean_binary = ifelse(data$subtype_clean=="PPMS","PMS","RMS")
all_coefs = list()
for(i in c(1:10000)){
  dat = data %>% filter(!is.na(age_at_dx))
  dat = sample_n(dat,size = nrow(dat),replace=T)
  dat$predicted_ancestry = relevel(factor(dat$predicted_ancestry),ref="EUR")
  coefs = lm(data = dat,age_at_dx ~ Sex + subtype_clean_binary + predicted_ancestry) %>% broom::tidy()
  all_coefs[[i]] = coefs
}
all_coefs = do.call("bind_rows",all_coefs)

# main model 
main_model = lm(data = data,age_at_dx ~ Sex  + subtype_clean_binary + predicted_ancestry) %>% broom::tidy()

model_res = all_coefs %>%
  group_by(term) %>% 
  summarise(sd_est = sd(estimate), lower = quantile(estimate,0.025), upper = quantile(estimate,0.975)) %>% 
  left_join(main_model %>%  dplyr::select(term,estimate) %>% dplyr::rename("beta" = estimate),by="term") %>%
  mutate(z = beta/sd_est) %>% 
  mutate(p = 2 * (1 - pnorm(abs(z)))) %>% 
  dplyr::select(term,beta,lower,upper,p)
write_csv(model_res,"./outputs/age_dx_model_subtype.csv")

# add site 
all_coefs = list()
for(i in c(1:10000)){
  dat = data %>% filter(!is.na(age_at_dx))
  dat = sample_n(dat,size = nrow(dat),replace=T)
  dat$predicted_ancestry = relevel(factor(dat$predicted_ancestry),ref="EUR")
  coefs = lm(data = dat,age_at_dx ~ Sex + subtype_clean_binary + Site + predicted_ancestry) %>% broom::tidy()
  all_coefs[[i]] = coefs
}
all_coefs = do.call("bind_rows",all_coefs)

# main model 
main_model = lm(data = data,age_at_dx ~ Sex  + subtype_clean_binary + Site + predicted_ancestry) %>% broom::tidy()

model_res = all_coefs %>%
  group_by(term) %>% 
  summarise(sd_est = sd(estimate), lower = quantile(estimate,0.025), upper = quantile(estimate,0.975)) %>% 
  left_join(main_model %>%  dplyr::select(term,estimate) %>% dplyr::rename("beta" = estimate),by="term") %>%
  mutate(z = beta/sd_est) %>% 
  mutate(p = 2 * (1 - pnorm(abs(z)))) %>% 
  dplyr::select(term,beta,lower,upper,p)
write_csv(model_res,"./outputs/age_dx_model_subtype_site.csv")

# stratified model - just website
# add site 
all_coefs = list()
for(i in c(1:10000)){
  dat = data %>% filter(!is.na(age_at_dx) & Site == "Website")
  dat = sample_n(dat,size = nrow(dat),replace=T)
  dat$predicted_ancestry = relevel(factor(dat$predicted_ancestry),ref="EUR")
  coefs = lm(data = dat,age_at_dx ~ Sex +  predicted_ancestry) %>% broom::tidy()
  all_coefs[[i]] = coefs
}
all_coefs = do.call("bind_rows",all_coefs)

# main model 
main_model = lm(data = data%>% filter(!is.na(age_at_dx) & Site == "Website"),age_at_dx ~ Sex  +  predicted_ancestry) %>% broom::tidy()
data%>% filter(!is.na(age_at_dx) & Site == "Website") %>% dplyr::count(predicted_ancestry)

model_res = all_coefs %>%
  group_by(term) %>% 
  summarise(sd_est = sd(estimate), lower = quantile(estimate,0.025), upper = quantile(estimate,0.975)) %>% 
  left_join(main_model %>%  dplyr::select(term,estimate) %>% dplyr::rename("beta" = estimate),by="term") %>%
  mutate(z = beta/sd_est) %>% 
  mutate(p = 2 * (1 - pnorm(abs(z)))) %>% 
  dplyr::select(term,beta,lower,upper,p)
write_csv(model_res,"./outputs/age_dx_model_subtype_stratified_by_site_website.csv")

# repeat for clinical sites
all_coefs = list()
for(i in c(1:10000)){
  dat = data %>% filter(!is.na(age_at_dx) & Site == "Other")
  dat = sample_n(dat,size = nrow(dat),replace=T)
  dat$predicted_ancestry = relevel(factor(dat$predicted_ancestry),ref="EUR")
  coefs = lm(data = dat,age_at_dx ~ Sex +  predicted_ancestry) %>% broom::tidy()
  all_coefs[[i]] = coefs
}
all_coefs = do.call("bind_rows",all_coefs)

# main model 
main_model = lm(data = data%>% filter(!is.na(age_at_dx) & Site == "Other"),age_at_dx ~ Sex   + predicted_ancestry) %>% broom::tidy()
data%>% filter(!is.na(age_at_dx) & Site == "Other") %>% dplyr::count(predicted_ancestry)
model_res = all_coefs %>%
  group_by(term) %>% 
  summarise(sd_est = sd(estimate), lower = quantile(estimate,0.025), upper = quantile(estimate,0.975)) %>% 
  left_join(main_model %>%  dplyr::select(term,estimate) %>% dplyr::rename("beta" = estimate),by="term") %>%
  mutate(z = beta/sd_est) %>% 
  mutate(p = 2 * (1 - pnorm(abs(z)))) %>% 
  dplyr::select(term,beta,lower,upper,p)
write_csv(model_res,"./outputs/age_dx_model_subtype_stratified_by_site_clinical.csv")



##################### 
# subtype 
##################### 
data$subtype_clean_binary = ifelse(data$subtype_clean=="PPMS",1,0)

main_model = glm(data = data, subtype_clean_binary ~ predicted_ancestry,family=binomial(link="logit")) %>% broom::tidy() %>% mutate(or = exp(estimate))
sex_model = glm(data = data, subtype_clean_binary ~ Sex + predicted_ancestry,family=binomial(link="logit")) %>% broom::tidy() %>% mutate(or = exp(estimate))
sex_site_model = glm(data = data, subtype_clean_binary ~ Site + Sex + predicted_ancestry,family=binomial(link="logit")) %>% broom::tidy() %>% mutate(or = exp(estimate))

major_anc_model = glm(data = data %>% filter(predicted_ancestry %in% c("EUR","CSA","AFR")), subtype_clean_binary ~ Sex + predicted_ancestry,family=binomial(link="logit")) %>% broom::tidy() %>% mutate(or = exp(estimate))


##################### 
# plots 
##################### 

# reorder ancestry for plots 
data$predicted_ancestry = factor(data$predicted_ancestry,levels = c("AFR","CSA","EUR","MID"))

# univariable comparisons
csa_armss = data[data$predicted_ancestry=="CSA",]$gARMSS
afr_armss = data[data$predicted_ancestry=="AFR",]$gARMSS
eur_armss = data[data$predicted_ancestry=="EUR",]$gARMSS
mid_armss = data[data$predicted_ancestry=="MID",]$gARMSS

wilcox.test(csa_armss,eur_armss)
wilcox.test(afr_armss,eur_armss)
wilcox.test(mid_armss,eur_armss)

# plot EDSS by age 
data$age_at_edss_bin = Hmisc::cut2(data$ageatedss,cuts = c(30,40,50))

ggplot(data = data %>% filter(!is.na(age_at_edss_bin)), aes(age_at_edss_bin,edss,fill=predicted_ancestry))+
  geom_boxplot()

# simple models - armss
## rank normalise 
data$Site = relevel(factor(data$Site),ref="Website")
data$subtype_clean_binary = ifelse(data$subtype_clean=="PPMS","PMS","RMS")
data$subtype_clean_binary = relevel(factor(data$subtype_clean_binary),ref="RMS")
data$predicted_ancestry = relevel(factor(data$predicted_ancestry),ref="EUR")

model_dat = data %>% filter(!is.na(gARMSS))
model_dat$armss_norm = RNOmni::RankNorm(model_dat$gARMSS)
unadjusted_model = glm(data = model_dat, armss_norm ~ predicted_ancestry) %>% broom::tidy() %>% mutate(model = "Unadjusted",outcome = "ARMSS") %>% mutate(lower_ci = estimate - 1.96 * `std.error`,upper_ci = estimate + 1.96 * `std.error`) %>%
mutate(n = nrow(model_dat), 
  n_eur = nrow(model_dat %>% filter(predicted_ancestry=="EUR")),
  n_afr = nrow(model_dat %>% filter(predicted_ancestry=="AFR")),
  n_csa = nrow(model_dat %>% filter(predicted_ancestry=="CSA")),
  n_mid = nrow(model_dat %>% filter(predicted_ancestry=="MID")))



# just adjust for sex 
sex_model = glm(data = model_dat, armss_norm ~ Sex + predicted_ancestry) %>% broom::tidy() %>% mutate(model = "Sex",outcome = "ARMSS") %>% mutate(lower_ci = estimate - 1.96 * `std.error`,upper_ci = estimate + 1.96 * `std.error`) %>%
mutate(n = nrow(model_dat), 
  n_eur = nrow(model_dat %>% filter(predicted_ancestry=="EUR")),
  n_afr = nrow(model_dat %>% filter(predicted_ancestry=="AFR")),
  n_csa = nrow(model_dat %>% filter(predicted_ancestry=="CSA")),
  n_mid = nrow(model_dat %>% filter(predicted_ancestry=="MID")))

# just adjust for age 
age_model = glm(data = model_dat, armss_norm ~ ageatedss + predicted_ancestry) %>% broom::tidy() %>% mutate(model = "Age",outcome = "ARMSS") %>% mutate(lower_ci = estimate - 1.96 * `std.error`,upper_ci = estimate + 1.96 * `std.error`) %>%
mutate(n = nrow(model_dat), 
  n_eur = nrow(model_dat %>% filter(predicted_ancestry=="EUR")),
  n_afr = nrow(model_dat %>% filter(predicted_ancestry=="AFR")),
  n_csa = nrow(model_dat %>% filter(predicted_ancestry=="CSA")),
  n_mid = nrow(model_dat %>% filter(predicted_ancestry=="MID")))

# just adjust for site
site_model = glm(data = model_dat, armss_norm ~ Site + predicted_ancestry) %>% broom::tidy() %>% mutate(model = "Site",outcome = "ARMSS") %>% mutate(lower_ci = estimate - 1.96 * `std.error`,upper_ci = estimate + 1.96 * `std.error`)%>%
mutate(n = nrow(model_dat), 
  n_eur = nrow(model_dat %>% filter(predicted_ancestry=="EUR")),
  n_afr = nrow(model_dat %>% filter(predicted_ancestry=="AFR")),
  n_csa = nrow(model_dat %>% filter(predicted_ancestry=="CSA")),
  n_mid = nrow(model_dat %>% filter(predicted_ancestry=="MID")))

# just adjust for age at diagnosis
age_at_dx_model = glm(data = model_dat, armss_norm ~ age_at_dx + predicted_ancestry) %>% broom::tidy() %>% mutate(model = "Age at diagnosis",outcome = "ARMSS") %>% mutate(lower_ci = estimate - 1.96 * `std.error`,upper_ci = estimate + 1.96 * `std.error`)%>%
mutate(n = nrow(model_dat), 
  n_eur = nrow(model_dat %>% filter(predicted_ancestry=="EUR")),
  n_afr = nrow(model_dat %>% filter(predicted_ancestry=="AFR")),
  n_csa = nrow(model_dat %>% filter(predicted_ancestry=="CSA")),
  n_mid = nrow(model_dat %>% filter(predicted_ancestry=="MID")))

# adjust for age, sex, and site
age_sex_site = glm(data = model_dat, armss_norm ~ ageatedss + Site + Sex + predicted_ancestry) %>% broom::tidy() %>% mutate(model = "Age + Sex + Site",outcome = "ARMSS") %>% mutate(lower_ci = estimate - 1.96 * `std.error`,upper_ci = estimate + 1.96 * `std.error`)%>%
mutate(n = nrow(model_dat), 
  n_eur = nrow(model_dat %>% filter(predicted_ancestry=="EUR")),
  n_afr = nrow(model_dat %>% filter(predicted_ancestry=="AFR")),
  n_csa = nrow(model_dat %>% filter(predicted_ancestry=="CSA")),
  n_mid = nrow(model_dat %>% filter(predicted_ancestry=="MID")))

# adjust for age & sex
age_sex = glm(data = model_dat, armss_norm ~ ageatedss + Sex + predicted_ancestry) %>% broom::tidy() %>% mutate(model = "Age + Sex",outcome = "ARMSS") %>% mutate(lower_ci = estimate - 1.96 * `std.error`,upper_ci = estimate + 1.96 * `std.error`)%>%
mutate(n = nrow(model_dat), 
  n_eur = nrow(model_dat %>% filter(predicted_ancestry=="EUR")),
  n_afr = nrow(model_dat %>% filter(predicted_ancestry=="AFR")),
  n_csa = nrow(model_dat %>% filter(predicted_ancestry=="CSA")),
  n_mid = nrow(model_dat %>% filter(predicted_ancestry=="MID")))

# stratified models 
just_clinical_sites = glm(data = model_dat %>% filter(Site == "Other"), armss_norm ~ predicted_ancestry) %>% broom::tidy() %>% mutate(model = "Unadjusted (clinical sites)",outcome = "ARMSS") %>% mutate(lower_ci = estimate - 1.96 * `std.error`,upper_ci = estimate + 1.96 * `std.error`)%>%
mutate(n = nrow(model_dat), 
  n_eur = nrow(model_dat %>% filter(predicted_ancestry=="EUR")),
  n_afr = nrow(model_dat %>% filter(predicted_ancestry=="AFR")),
  n_csa = nrow(model_dat %>% filter(predicted_ancestry=="CSA")),
  n_mid = nrow(model_dat %>% filter(predicted_ancestry=="MID")))


just_website = glm(data = model_dat %>% filter(Site == "Website"), armss_norm ~ predicted_ancestry) %>% broom::tidy() %>% mutate(model = "Unadjusted (website)",outcome = "ARMSS") %>% mutate(lower_ci = estimate - 1.96 * `std.error`,upper_ci = estimate + 1.96 * `std.error`)%>%
mutate(n = nrow(model_dat), 
  n_eur = nrow(model_dat %>% filter(predicted_ancestry=="EUR")),
  n_afr = nrow(model_dat %>% filter(predicted_ancestry=="AFR")),
  n_csa = nrow(model_dat %>% filter(predicted_ancestry=="CSA")),
  n_mid = nrow(model_dat %>% filter(predicted_ancestry=="MID")))

just_ukmsr = glm(data = model_dat %>% filter(Site == "UKMSR"), armss_norm ~ predicted_ancestry) %>% broom::tidy() %>% mutate(model = "Unadjusted (website)",outcome = "ARMSS") %>% mutate(lower_ci = estimate - 1.96 * `std.error`,upper_ci = estimate + 1.96 * `std.error`)%>%
mutate(n = nrow(model_dat), 
  n_eur = nrow(model_dat %>% filter(predicted_ancestry=="EUR")),
  n_afr = nrow(model_dat %>% filter(predicted_ancestry=="AFR")),
  n_csa = nrow(model_dat %>% filter(predicted_ancestry=="CSA")),
  n_mid = nrow(model_dat %>% filter(predicted_ancestry=="MID")))

# stratify by sex
just_men = glm(data = model_dat %>% filter(Sex == "M"), armss_norm ~ predicted_ancestry) %>% broom::tidy() %>% mutate(model = "Unadjusted (just males)",outcome = "ARMSS") %>% mutate(lower_ci = estimate - 1.96 * `std.error`,upper_ci = estimate + 1.96 * `std.error`)%>%
mutate(n = nrow(model_dat), 
  n_eur = nrow(model_dat %>% filter(predicted_ancestry=="EUR")),
  n_afr = nrow(model_dat %>% filter(predicted_ancestry=="AFR")),
  n_csa = nrow(model_dat %>% filter(predicted_ancestry=="CSA")),
  n_mid = nrow(model_dat %>% filter(predicted_ancestry=="MID")))

just_women = glm(data = model_dat %>% filter(Sex == "F"), armss_norm ~ predicted_ancestry) %>% broom::tidy() %>% mutate(model = "Unadjusted (just males)",outcome = "ARMSS") %>% mutate(lower_ci = estimate - 1.96 * `std.error`,upper_ci = estimate + 1.96 * `std.error`)%>%
mutate(n = nrow(model_dat), 
  n_eur = nrow(model_dat %>% filter(predicted_ancestry=="EUR")),
  n_afr = nrow(model_dat %>% filter(predicted_ancestry=="AFR")),
  n_csa = nrow(model_dat %>% filter(predicted_ancestry=="CSA")),
  n_mid = nrow(model_dat %>% filter(predicted_ancestry=="MID")))

# 5y of dx
model_dat = data %>% filter(!is.na(gARMSS))
model_dat = model_dat %>% 
  mutate(disease_duration_at_edss = ageatedss - age_at_dx) %>% 
  filter(disease_duration_at_edss < 5)

model_dat$armss_norm = RNOmni::RankNorm(model_dat$gARMSS)



age_sex_site_within_5y_of_dx = glm(data = model_dat, armss_norm ~ ageatedss + Site + Sex + predicted_ancestry) %>% broom::tidy() %>% mutate(model = "Age + Sex + Site (<5y from dx)",outcome = "ARMSS") %>% mutate(lower_ci = estimate - 1.96 * `std.error`,upper_ci = estimate + 1.96 * `std.error`)%>%
mutate(n = nrow(model_dat), 
  n_eur = nrow(model_dat %>% filter(predicted_ancestry=="EUR")),
  n_afr = nrow(model_dat %>% filter(predicted_ancestry=="AFR")),
  n_csa = nrow(model_dat %>% filter(predicted_ancestry=="CSA")),
  n_mid = nrow(model_dat %>% filter(predicted_ancestry=="MID")))

# just directly observed EDSS
model_dat = data %>% filter(!is.na(gARMSS))
model_dat = model_dat %>% 
  filter(edss_source == "Observed")

model_dat$armss_norm = RNOmni::RankNorm(model_dat$gARMSS)
age_sex_site_direct_edss = glm(data = model_dat, armss_norm ~ ageatedss + Site + Sex + predicted_ancestry) %>% broom::tidy() %>% mutate(model = "Age + Sex + Site (direct EDSS)",outcome = "ARMSS") %>% mutate(lower_ci = estimate - 1.96 * `std.error`,upper_ci = estimate + 1.96 * `std.error`)%>%
mutate(n = nrow(model_dat), 
  n_eur = nrow(model_dat %>% filter(predicted_ancestry=="EUR")),
  n_afr = nrow(model_dat %>% filter(predicted_ancestry=="AFR")),
  n_csa = nrow(model_dat %>% filter(predicted_ancestry=="CSA")),
  n_mid = nrow(model_dat %>% filter(predicted_ancestry=="MID")))

# primary analysis - EDSS itself 
model_dat = data %>% filter(!is.na(edss))
model_dat$edss_norm = RNOmni::RankNorm(model_dat$edss)
age_sex_site_edss = glm(data = model_dat, edss_norm ~ ageatedss + Site + Sex + predicted_ancestry) %>% broom::tidy() %>% mutate(model = "Age + Sex + Site",outcome = "EDSS") %>% mutate(lower_ci = estimate - 1.96 * `std.error`,upper_ci = estimate + 1.96 * `std.error`)%>%
mutate(n = nrow(model_dat), 
  n_eur = nrow(model_dat %>% filter(predicted_ancestry=="EUR")),
  n_afr = nrow(model_dat %>% filter(predicted_ancestry=="AFR")),
  n_csa = nrow(model_dat %>% filter(predicted_ancestry=="CSA")),
  n_mid = nrow(model_dat %>% filter(predicted_ancestry=="MID")))

# primary analysis - MSSS
model_dat = data %>% filter(!is.na(uGMSSS))
model_dat$msss_norm = RNOmni::RankNorm(model_dat$uGMSSS)
age_sex_site_msss = glm(data = model_dat, msss_norm ~ ageatedss + Site + Sex + predicted_ancestry) %>% broom::tidy() %>% mutate(model = "Age + Sex + Site",outcome = "MSSS") %>% mutate(lower_ci = estimate - 1.96 * `std.error`,upper_ci = estimate + 1.96 * `std.error`)%>%
mutate(n = nrow(model_dat), 
  n_eur = nrow(model_dat %>% filter(predicted_ancestry=="EUR")),
  n_afr = nrow(model_dat %>% filter(predicted_ancestry=="AFR")),
  n_csa = nrow(model_dat %>% filter(predicted_ancestry=="CSA")),
  n_mid = nrow(model_dat %>% filter(predicted_ancestry=="MID")))

# primary analysis - MSIS
model_dat = data %>% filter(!is.na(msis_physical_normalised))
model_dat$msis_phys_norm = RNOmni::RankNorm(model_dat$msis_physical_normalised)
age_sex_site_msis_phys_norm = glm(data = model_dat, msis_phys_norm ~ age_at_msis + Site + Sex + predicted_ancestry) %>% broom::tidy() %>% mutate(model = "Age + Sex + Site",outcome = "MSIS (Phys)") %>% mutate(lower_ci = estimate - 1.96 * `std.error`,upper_ci = estimate + 1.96 * `std.error`)%>%
mutate(n = nrow(model_dat), 
  n_eur = nrow(model_dat %>% filter(predicted_ancestry=="EUR")),
  n_afr = nrow(model_dat %>% filter(predicted_ancestry=="AFR")),
  n_csa = nrow(model_dat %>% filter(predicted_ancestry=="CSA")),
  n_mid = nrow(model_dat %>% filter(predicted_ancestry=="MID")))

# primary analysis - MSIS Psych
model_dat = data %>% filter(!is.na(msis_psych_normalised))
model_dat$msis_psych_norm = RNOmni::RankNorm(model_dat$msis_psych_normalised)
age_sex_site_msis_psych_norm = glm(data = model_dat, msis_psych_norm ~ age_at_msis + Site + Sex + predicted_ancestry) %>% broom::tidy() %>% mutate(model = "Age + Sex + Site",outcome = "MSIS (Psych)") %>% mutate(lower_ci = estimate - 1.96 * `std.error`,upper_ci = estimate + 1.96 * `std.error`)%>%
mutate(n = nrow(model_dat), 
  n_eur = nrow(model_dat %>% filter(predicted_ancestry=="EUR")),
  n_afr = nrow(model_dat %>% filter(predicted_ancestry=="AFR")),
  n_csa = nrow(model_dat %>% filter(predicted_ancestry=="CSA")),
  n_mid = nrow(model_dat %>% filter(predicted_ancestry=="MID")))

# primary analysis - EQ5D 
model_dat = data %>% filter(!is.na(eq5d_vas))

model_dat = model_dat %>%
  mutate(age_at_eq5d = age_at_recruitment + as.numeric(as.Date(date_of_eq5d_vas) - `Recruitment date`)/365.25) 

model_dat$eq5d_norm = RNOmni::RankNorm(model_dat$eq5d_vas)
age_sex_site_eq5d = glm(data = model_dat, eq5d_norm ~ age_at_eq5d + Site + Sex + predicted_ancestry) %>% broom::tidy() %>% mutate(model = "Age + Sex + Site",outcome = "EQ5D") %>% mutate(lower_ci = estimate - 1.96 * `std.error`,upper_ci = estimate + 1.96 * `std.error`)%>%
mutate(n = nrow(model_dat), 
  n_eur = nrow(model_dat %>% filter(predicted_ancestry=="EUR")),
  n_afr = nrow(model_dat %>% filter(predicted_ancestry=="AFR")),
  n_csa = nrow(model_dat %>% filter(predicted_ancestry=="CSA")),
  n_mid = nrow(model_dat %>% filter(predicted_ancestry=="MID")))


# combine all of the above models 
overall_model_dat = bind_rows(
  age_sex_site_eq5d,
  age_sex_site_msis_psych_norm,
  age_sex_site_msis_phys_norm,
  age_sex_site_msss,
  age_sex_site_edss,
  age_sex_site,
  just_men,just_women,just_ukmsr,just_website,just_clinical_sites,
age_sex,
age_sex_site_direct_edss,age_sex_site_within_5y_of_dx,
site_model,
age_model,
unadjusted_model
) %>%
  write_csv(file="./outputs/severity_models.csv")

# plot 
plot_dat = overall_model_dat %>% 
  filter(grepl("predicted",term)) %>%
  mutate(ancestry = str_remove_all(term,"predicted_ancestry")) %>% 
  dplyr::select(ancestry, model, outcome, p.value, estimate, lower_ci, upper_ci)


png("./plots/forest_plot.png",res=900,units="in",width=4,height=8)
ggplot(plot_dat %>% filter(model == "Age + Sex + Site"),aes(estimate,ancestry,fill=outcome))+
  geom_vline(xintercept = 0,alpha=0.9, linetype = "dashed", color="pink")+
  scale_fill_brewer(palette="Set3")+
  geom_errorbarh(mapping = aes(xmin = lower_ci, xmax = upper_ci, y = ancestry),height=0.1)+
  facet_wrap(~outcome,nrow=6)+
    geom_point(shape=21,color="black",size=3)+
  labs(x="Estimate (change in rank-normalised\noutcome vs European ancestry)",y="Ancestry ")+
  theme_bw()+
  theme(panel.grid = element_blank(),legend.position="none")
dev.off()


# plot armss
counts = data %>% dplyr::count(predicted_ancestry)
medians = data %>% group_by(predicted_ancestry) %>% summarise(median_armss = median(gARMSS,na.rm=T)) 

p1 = ggplot(data,aes(predicted_ancestry,gARMSS,fill=predicted_ancestry))+
geom_violin(alpha=0.5)+
geom_boxplot(width=0.1,alpha=0.3)+
geom_jitter(width=0.1,alpha=0.1)+
theme_bw()+
scale_y_continuous(limits=c(0,12),breaks=seq(0,10,by=1))+
labs(x="Ancestry",y="ARMSS",fill="Ancestry")+
scale_fill_brewer(palette="Set1")+
geom_text(data = medians,aes(y=11,label = format(round(median_armss,1),nsmall=1)))+
theme(legend.position="none",panel.grid = element_blank())

# plot msss
counts = data %>% dplyr::count(predicted_ancestry)
medians = data %>% group_by(predicted_ancestry) %>% summarise(median_msss = median(uGMSSS,na.rm=T)) 

p11 = ggplot(data,aes(predicted_ancestry,uGMSSS,fill=predicted_ancestry))+
geom_violin(alpha=0.5)+
geom_boxplot(width=0.1,alpha=0.3)+
geom_jitter(width=0.1,alpha=0.1)+
theme_bw()+
scale_y_continuous(limits=c(0,12),breaks=seq(0,10,by=1))+
labs(x="Ancestry",y="MSSS",fill="Ancestry")+
scale_fill_brewer(palette="Set1")+
geom_text(data = medians,aes(y=11,label = format(round(median_msss,1),nsmall=1)))+
theme(legend.position="none",panel.grid = element_blank())


# plot edss
counts = data %>% dplyr::count(predicted_ancestry)
medians = data %>% group_by(predicted_ancestry) %>% summarise(median_edss = median(edss,na.rm=T)) 

p12 = ggplot(data,aes(predicted_ancestry,edss,fill=predicted_ancestry))+
geom_violin(alpha=0.5)+
geom_boxplot(width=0.1,alpha=0.3)+
geom_jitter(width=0.1,alpha=0.1)+
theme_bw()+
scale_y_continuous(limits=c(0,12),breaks=seq(0,10,by=1))+
labs(x="Ancestry",y="EDSS",fill="Ancestry")+
scale_fill_brewer(palette="Set1")+
geom_text(data = medians,aes(y=11,label = format(round(median_edss,1),nsmall=1)))+
theme(legend.position="none",panel.grid = element_blank())


# plot msis
counts = data %>% dplyr::count(predicted_ancestry)
medians = data %>% group_by(predicted_ancestry) %>% summarise(median_msis = median(msis_physical_normalised,na.rm=T)) 

p2 = ggplot(data,aes(predicted_ancestry,msis_physical_normalised,fill=predicted_ancestry))+
geom_violin(alpha=0.5)+
geom_boxplot(width=0.1,alpha=0.3)+
geom_jitter(width=0.1,alpha=0.1)+
theme_bw()+
labs(x="Ancestry",y="MSIS (Physical)",fill="Ancestry")+
scale_fill_brewer(palette="Set1")+
scale_y_continuous(limits=c(0,120),breaks=seq(0,100,by=10))+
geom_text(data = medians,aes(y=110,label = format(round(median_msis,1),nsmall=1)))+
theme(legend.position="none",panel.grid = element_blank())

# plot msis-psych
counts = data %>% dplyr::count(predicted_ancestry)
medians = data %>% group_by(predicted_ancestry) %>% summarise(median_msis = median(msis_psych_normalised,na.rm=T)) 

p3 = ggplot(data,aes(predicted_ancestry,msis_physical_normalised,fill=predicted_ancestry))+
geom_violin(alpha=0.5)+
geom_boxplot(width=0.1,alpha=0.3)+
geom_jitter(width=0.1,alpha=0.1)+
theme_bw()+
labs(x="Ancestry",y="MSIS (Psychological)",fill="Ancestry")+
scale_fill_brewer(palette="Set1")+
scale_y_continuous(limits=c(0,120),breaks=seq(0,100,by=10))+
geom_text(data = medians,aes(y=110,label = format(round(median_msis,1),nsmall=1)))+
theme(legend.position="none",panel.grid = element_blank())

# plot eq5d
counts = data %>% dplyr::count(predicted_ancestry)
medians = data %>% group_by(predicted_ancestry) %>% summarise(median_eq5d = median(eq5d_vas,na.rm=T)) 

p4 = ggplot(data,aes(predicted_ancestry,eq5d_vas,fill=predicted_ancestry))+
geom_violin(alpha=0.5)+
geom_boxplot(width=0.1,alpha=0.3)+
geom_jitter(width=0.1,alpha=0.1)+
theme_bw()+
labs(x="Ancestry",y="EQ5D",fill="Ancestry")+
scale_fill_brewer(palette="Set1")+
scale_y_continuous(limits=c(0,120),breaks=seq(0,100,by=10))+
geom_text(data = medians,aes(y=110,label = format(round(median_eq5d,1),nsmall=1)))+
theme(legend.position="none",panel.grid = element_blank())

png("/data/home/hmy117/ADAMS_severity/plots/msis_and_garmss_vs_ancestry.png",res=900,units="in",width=6,height=12)
gridExtra::grid.arrange(p1,p12,p4,p2,p3,p11,nrow=6)

dev.off()


# armss vs global ancestry 
admixed = data %>% filter(EUR > 0.1 & AFR > 0.1)
ggplot(admixed,aes(EUR,gARMSS))+
  geom_point(shape=21,size=3)+
  geom_smooth(method="lm",se=F,linetype="dashed",color="pink")+
  theme_bw()+
  labs(y="ARMSS")
model_dat = admixed %>% filter(!is.na(gARMSS))
model_dat$armss_norm = RNOmni::RankNorm(model_dat$gARMSS)
glm(data = model_dat, armss_norm ~ ageatedss + Site + Sex + subtype_clean + AFR) %>% broom::tidy() 



````

# GWAS

GWAS of severity, and replication analysis of DYSF-ZNF638 SNPs.

## Remove duplicates & copy back to home
````unix 
cd /data/scratch/hmy117/adams_imputed_severity_topmed/
~/plink --bfile combined_adams_imputed_qc \
--remove dups_to_discard.tsv \
--keep /data/home/hmy117/ADAMS_severity/outputs/adams_pheno.tsv \
--mind 0.1 \
--out /data/home/hmy117/ADAMS_severity/combined_adams_imputed \
--make-bed

````

## PC plots 
````R 

library(tidyverse)

ancestry = read_tsv("/data/home/hmy117/ADAMS_severity/outputs/ancestry_calls.tsv",col_types="ccc")
pheno = read_tsv("/data/home/hmy117/ADAMS_severity/pheno/adams_pheno.tsv",col_types="ccddddd") %>% dplyr::select(-FID)
cov = read_tsv("/data/home/hmy117/ADAMS_severity/pheno/adams_covars.tsv",col_types="ccddddd") %>% dplyr::select(-FID)
fam = read_table("/data/home/hmy117/ADAMS_severity/combined_adams_imputed.fam",col_names=F)

# join 
dat = ancestry %>% 
  left_join(pheno,by="IID") %>% 
  left_join(cov,by="IID") %>% 
  filter(IID %in% fam$X2)

# add pcs 
# read in data
adams = read_table("/data/scratch/hmy117/hgdp_1kg_genomes/adams_pcs.sscore")
dat = dat %>% 
  left_join(adams,by="IID") 

# reference
kg_hgdp = read_table("/data/scratch/hmy117/hgdp_1kg_genomes/hgdp_kg_aj_pcs_rescored.sscore")
meta = read_tsv("/data/scratch/hmy117/aj_ref/gnomad_meta_v1.tsv") %>%
  dplyr::select(2,hgdp_tgp_meta.Population,hgdp_tgp_meta.Genetic.region)
colnames(meta) = c("IID","pop","superpop")
kg_hgdp = kg_hgdp %>%
left_join(meta,by="IID")
kg_hgdp = kg_hgdp %>%
mutate(superpop = ifelse(is.na(superpop),"ASJ",superpop))

# combine 
combo = bind_rows(kg_hgdp %>% mutate(study="Reference"),dat %>% mutate(study="ADAMS")) %>%
  mutate(ancestry = ifelse(!is.na(predicted_ancestry),predicted_ancestry,superpop))

# pc plots 
p = ggplot(data = combo,aes(PC1_AVG,PC2_AVG,fill=ancestry))+
  geom_point(shape=21,size=2,stroke=0.2)+
	theme_bw()+
  labs(x="PC1",y="PC2",fill="Ancestry")+
	scale_fill_brewer(palette="Accent")+
  theme(legend.position="top")+
  facet_wrap(~study)

png("/data/home/hmy117/ADAMS_severity/plots/all_pcs.png",res=900,units="in",width=6,height=4)
p
dev.off()


````

## Extract DYSF SNPs
````unix 

zgrep rs10191329 /data/scratch/hmy117/adams_imputed_severity_topmed/chr2.info.gz

echo rs7579497 > /data/scratch/hmy117/snps
echo rs10191329 >> /data/scratch/hmy117/snps
echo rs13384155 >> /data/scratch/hmy117/snps

~/plink2 --vcf /data/scratch/hmy117/adams_imputed_severity_topmed/chr2.dose.vcf.gz dosage=HDS \
--extract /data/scratch/hmy117/snps \
--recode A include-alt \
--out /data/scratch/hmy117/rs10191329_rawcalls

# ld 
for ancestry in CSA EUR AFR;
do
~/plink --bfile /data/home/hmy117/ADAMS_severity/outputs/imputed_genotypes_$ancestry \
--ld-snp '2:71449869:C:A' \
--ld-window-kb 9999999 \
--ld-window 9999999 \
--ld-window-r2 0 \
--r2 \
--out /data/home/hmy117/ADAMS_severity/outputs/rs10191329_ld_$ancestry
done

for ancestry in CSA EUR AFR;
do
~/plink --bfile /data/home/hmy117/ADAMS_severity/outputs/imputed_genotypes_$ancestry \
--snp '2:71449869:C:A' \
--freq \
--out /data/home/hmy117/ADAMS_severity/outputs/rs10191329_freq_$ancestry
done
````

qlogin -pe smp 1 -l h_vmem=256G -l h_rt=240:00:00
cd /data/home/hmy117/ADAMS_severity/

# Analysis in R
````R
library(tidyverse)

# plot ld 
res = list()
for(anc in c("EUR","CSA","AFR")){
  res[[length(res)+1]] = read_table(paste0("/data/home/hmy117/ADAMS_severity/outputs/rs10191329_ld_",anc,".ld")) %>%
    mutate(ancestry=anc)
}
res = do.call("bind_rows",res)

# plot 
p = ggplot(res %>% filter(BP_B > 70949869 & BP_B < 71949869),aes(BP_B,R2,fill=R2))+
geom_vline(xintercept=71449869,alpha=0.5)+
geom_point(shape=21,size=2)+
facet_wrap(~ancestry,nrow=3)+
theme_bw()+
scale_fill_viridis_c(option="plasma")+
labs(x="Position on chromosome 2",y=bquote(R^2),fill=bquote(R^2))

png("/data/home/hmy117/ADAMS_severity/plots/rs10191329_ld.png",res=900,units="in",width=6,height=6)
p
dev.off()


rs10191329 = read_table("/data/scratch/hmy117/rs10191329_rawcalls.raw")
ancestry = read_tsv("/data/home/hmy117/ADAMS_severity/outputs/ancestry_calls.tsv",col_types="ccc")
pheno = read_tsv("/data/home/hmy117/ADAMS_severity/pheno/adams_pheno.tsv",col_types="ccddddd") %>% dplyr::select(-FID)
cov = read_tsv("/data/home/hmy117/ADAMS_severity/pheno/adams_covars.tsv",col_types="ccdccccd") %>% dplyr::select(-FID)

# merge 
geno = rs10191329 %>% 
separate(IID,sep="_",into=c("IID","Other")) %>%
inner_join(pheno,by="IID") %>% 
inner_join(cov,by="IID") %>% 
inner_join(ancestry,by="IID") 

# additive alleles 
geno$rs10191329_A = 2 - geno$`rs10191329_C(/A)`
geno$rs7579497_A = 2 - geno$`rs7579497_G(/A)`
geno$rs13384155_C = 2 - geno$ `rs13384155_T(/C)`

# linear models for each phenotype, ancestry and SNP 
res = list()
for(phenotype in c("gARMSS","edss","msis_physical_normalised","eq5d_vas","uGMSSS")){

  model_dat = geno %>% filter(!is.na(geno[[phenotype]]))
  model_dat$pheno_norm = RNOmni::RankNorm(model_dat[[phenotype]])

  for(ancestry in c("CSA","AFR","EUR")){
    message(ancestry)
    message(phenotype)
      model_dat_anc = model_dat %>% filter(predicted_ancestry==ancestry)
      n = nrow(model_dat_anc)

      # add pcs 
      pcs = read_table(paste0("/data/home/hmy117/ADAMS_severity/outputs/pcs_",ancestry,".eigenvec"),col_types="ccdddd")

      # join 
      model_dat_anc = model_dat_anc %>% 
        inner_join(pcs,by="IID")


      rs10191329_counts = model_dat_anc %>% dplyr::count(rs10191329 = round(rs10191329_A,0)) %>% 
        pivot_wider(names_from = rs10191329, values_from = n) 
      colnames(rs10191329_counts) = paste0("rs10191329_A_",colnames(rs10191329_counts))
    
      rs7579497_counts = model_dat_anc %>% dplyr::count(rs7579497 = round(rs7579497_A,0)) %>% 
        pivot_wider(names_from = rs7579497, values_from = n) 
      colnames(rs7579497_counts) = paste0("rs7579497_A_",colnames(rs7579497_counts))
    
      rs13384155_counts = model_dat_anc %>% dplyr::count(rs13384155 = round(rs13384155_C,0)) %>% 
        pivot_wider(names_from = rs13384155, values_from = n) 
      colnames(rs13384155_counts) = paste0("rs13384155_C_",colnames(rs13384155_counts))
    
      # additive

    if(!is.null(rs10191329_counts[["rs10191329_A_1"]])){
    
    if(
      (rs10191329_counts[["rs10191329_A_1"]])
       > 2){
    message("additive")

      
      res[[length(res)+1]] = lm(data = model_dat_anc,pheno_norm ~ batch +ageatedss + Sex +  PC1 +PC2 + PC3 + PC4 +rs10191329_A) %>% broom::tidy() %>%
        mutate(lower_ci = estimate - 1.96*std.error, upper_ci = estimate + 1.96*std.error) %>% 
        filter(term == "rs10191329_A") %>% 
        dplyr::select(term,estimate,lower_ci,upper_ci,p.value) %>% 
        mutate(pheno = phenotype, anc = ancestry, N = n) %>% 
        bind_cols(rs10191329_counts) %>% mutate(model="Additive")
    }
    }

    # dominant 
 if(!is.null(rs10191329_counts[["rs10191329_A_1"]])){
    
    if(
      (rs10191329_counts[["rs10191329_A_1"]])
       > 2){
    message("dominant")

            model_dat_anc$rs10191329_carrier = ifelse(round(model_dat_anc$rs10191329_A,0)!=0,"carrier","non_carrier")
            model_dat_anc$rs10191329_carrier = relevel(factor(model_dat_anc$rs10191329_carrier),ref="non_carrier")

      res[[length(res)+1]] = lm(data = model_dat_anc,pheno_norm ~ batch +ageatedss + Sex + PC1 +PC2 + PC3 + PC4 +rs10191329_carrier) %>% broom::tidy() %>%
        mutate(lower_ci = estimate - 1.96*std.error, upper_ci = estimate + 1.96*std.error) %>% 
        filter(grepl("rs10191329_carrier",term)) %>% 
        dplyr::select(term,estimate,lower_ci,upper_ci,p.value) %>% 
        mutate(pheno = phenotype, anc = ancestry, N = n) %>% 
        bind_cols(rs10191329_counts) %>% mutate(model="Dominant")
    }
    }

    # recessive
 if(!is.null(rs10191329_counts[["rs10191329_A_2"]])){
    
    if(
      (rs10191329_counts[["rs10191329_A_2"]])
       > 2){
    message("recessive")

      model_dat_anc$rs10191329_hom = ifelse(round(model_dat_anc$rs10191329_A,0)==2,"hom","non_hom")
      model_dat_anc$rs10191329_hom = relevel(factor(model_dat_anc$rs10191329_hom),ref="non_hom")
      res[[length(res)+1]] = lm(data = model_dat_anc,pheno_norm ~ batch +ageatedss + Sex +  PC1 +PC2 + PC3 + PC4 +rs10191329_hom) %>% broom::tidy() %>%
        mutate(lower_ci = estimate - 1.96*std.error, upper_ci = estimate + 1.96*std.error) %>% 
        filter(grepl("rs10191329_hom",term)) %>% 
        dplyr::select(term,estimate,lower_ci,upper_ci,p.value) %>% 
        mutate(pheno = phenotype, anc = ancestry, N = n) %>% 
        bind_cols(rs10191329_counts) %>% mutate(model="Recessive")
    }
    }

  # repeat for rs7579497
      if(!is.null(rs7579497_counts[["rs7579497_A_1"]])){
    
    if(
      (rs7579497_counts[["rs7579497_A_1"]])
       > 2){
    message("additive")

      
      res[[length(res)+1]] = lm(data = model_dat_anc,pheno_norm ~ batch +ageatedss + Sex +  PC1 +PC2 + PC3 + PC4 +rs7579497_A) %>% broom::tidy() %>%
        mutate(lower_ci = estimate - 1.96*std.error, upper_ci = estimate + 1.96*std.error) %>% 
        filter(term == "rs7579497_A") %>% 
        dplyr::select(term,estimate,lower_ci,upper_ci,p.value) %>% 
        mutate(pheno = phenotype, anc = ancestry, N = n) %>% 
        bind_cols(rs7579497_counts) %>% mutate(model="Additive")
    }
    }

    # dominant 
 if(!is.null(rs7579497_counts[["rs7579497_A_1"]])){
    
    if(
      (rs7579497_counts[["rs7579497_A_1"]])
       > 2){
    message("dominant")

            model_dat_anc$rs7579497_carrier = ifelse(round(model_dat_anc$rs7579497_A,0)!=0,"carrier","non_carrier")
            model_dat_anc$rs7579497_carrier = relevel(factor(model_dat_anc$rs7579497_carrier),ref="non_carrier")

      res[[length(res)+1]] = lm(data = model_dat_anc,pheno_norm ~ batch +ageatedss + Sex +  PC1 +PC2 + PC3 + PC4 +rs7579497_carrier) %>% broom::tidy() %>%
        mutate(lower_ci = estimate - 1.96*std.error, upper_ci = estimate + 1.96*std.error) %>% 
        filter(grepl("rs7579497_carrier",term)) %>% 
        dplyr::select(term,estimate,lower_ci,upper_ci,p.value) %>% 
        mutate(pheno = phenotype, anc = ancestry, N = n) %>% 
        bind_cols(rs7579497_counts) %>% mutate(model="Dominant")
    }
    }

    # recessive
 if(!is.null(rs7579497_counts[["rs7579497_A_2"]])){
    
    if(
      (rs7579497_counts[["rs7579497_A_2"]])
       > 2){
    message("recessive")

      model_dat_anc$rs7579497_hom = ifelse(round(model_dat_anc$rs7579497_A,0)==2,"hom","non_hom")
      model_dat_anc$rs7579497_hom = relevel(factor(model_dat_anc$rs7579497_hom),ref="non_hom")
      res[[length(res)+1]] = lm(data = model_dat_anc,pheno_norm ~ batch +ageatedss + Sex +  PC1 +PC2 + PC3 + PC4 +rs7579497_hom) %>% broom::tidy() %>%
        mutate(lower_ci = estimate - 1.96*std.error, upper_ci = estimate + 1.96*std.error) %>% 
        filter(grepl("rs7579497_hom",term)) %>% 
        dplyr::select(term,estimate,lower_ci,upper_ci,p.value) %>% 
        mutate(pheno = phenotype, anc = ancestry, N = n) %>% 
        bind_cols(rs7579497_counts) %>% mutate(model="Recessive")
    }
    }

  
  # repeat
      if(!is.null(rs13384155_counts[["rs13384155_C_1"]])){
    
    if(
      (rs13384155_counts[["rs13384155_C_1"]])
       > 2){
    message("additive")

      
      res[[length(res)+1]] = lm(data = model_dat_anc,pheno_norm ~ batch +ageatedss + Sex +  PC1 +PC2 + PC3 + PC4 +rs13384155_C) %>% broom::tidy() %>%
        mutate(lower_ci = estimate - 1.96*std.error, upper_ci = estimate + 1.96*std.error) %>% 
        filter(term == "rs13384155_C") %>% 
        dplyr::select(term,estimate,lower_ci,upper_ci,p.value) %>% 
        mutate(pheno = phenotype, anc = ancestry, N = n) %>% 
        bind_cols(rs13384155_counts) %>% mutate(model="Additive")
    }
    }

    # dominant 
 if(!is.null(rs13384155_counts[["rs13384155_C_1"]])){
    
    if(
      (rs13384155_counts[["rs13384155_C_1"]])
       > 2){
    message("dominant")

            model_dat_anc$rs13384155_carrier = ifelse(round(model_dat_anc$rs13384155_C,0)!=0,"carrier","non_carrier")
            model_dat_anc$rs13384155_carrier = relevel(factor(model_dat_anc$rs13384155_carrier),ref="non_carrier")

      res[[length(res)+1]] = lm(data = model_dat_anc,pheno_norm ~ batch +ageatedss + Sex +  PC1 +PC2 + PC3 + PC4 +rs13384155_carrier) %>% broom::tidy() %>%
        mutate(lower_ci = estimate - 1.96*std.error, upper_ci = estimate + 1.96*std.error) %>% 
        filter(grepl("rs13384155_carrier",term)) %>% 
        dplyr::select(term,estimate,lower_ci,upper_ci,p.value) %>% 
        mutate(pheno = phenotype, anc = ancestry, N = n) %>% 
        bind_cols(rs13384155_counts) %>% mutate(model="Dominant")
    }
    }

    # recessive
 if(!is.null(rs13384155_counts[["rs13384155_C_2"]])){
    
    if(
      (rs13384155_counts[["rs13384155_C_2"]])
       > 2){
    message("recessive")

      model_dat_anc$rs13384155_hom = ifelse(round(model_dat_anc$rs13384155_C,0)==2,"hom","non_hom")
      model_dat_anc$rs13384155_hom = relevel(factor(model_dat_anc$rs13384155_hom),ref="non_hom")
      res[[length(res)+1]] = lm(data = model_dat_anc,pheno_norm ~ batch +ageatedss + Sex + PC1 +PC2 + PC3 + PC4 + rs13384155_hom) %>% broom::tidy() %>%
        mutate(lower_ci = estimate - 1.96*std.error, upper_ci = estimate + 1.96*std.error) %>% 
        filter(grepl("rs13384155_hom",term)) %>% 
        dplyr::select(term,estimate,lower_ci,upper_ci,p.value) %>% 
        mutate(pheno = phenotype, anc = ancestry, N = n) %>% 
        bind_cols(rs10191329_counts) %>% mutate(model="Recessive")
    }
    }

}
}

res = do.call("bind_rows",res)


# rename phenos 
res = res %>%
  mutate(pheno = case_when(
    pheno == "edss" ~ "EDSS",
    pheno == "gARMSS" ~ "ARMSS",
    pheno == "uGMSSS" ~ "MSSS",
    pheno == "msis_physical_normalised" ~ "MSIS-29 (Phys)",
    pheno == "eq5d_vas" ~ "EQ5D"
  )) 

# plots
p = ggplot(res %>% filter(model == "Additive"),aes(estimate,anc,fill=anc))+
  geom_vline(xintercept=0,alpha=0.5,linetype="dashed")+
  geom_errorbarh(mapping = aes(xmin = lower_ci,xmax=upper_ci,y=anc),height=0.1)+
  geom_point(shape=21,size=3)+
  theme_bw()+
  facet_grid(term~pheno)+
  theme(legend.position="none")+
  scale_fill_brewer(palette="Accent")+
  labs(x="Effect on normalised phenotype\nper severity allele",y="Ancestry",fill=NULL)


png("/data/home/hmy117/ADAMS_severity/plots/rs10191329_imsgc_snps_vs_phenos.png",res=900,units="in",width=8,height=4)
p
dev.off()

  
  res = res pivot_wider(id_cols = c(term,model,pheno),values_from = c(estimate,lower_ci,upper_ci,p.value,N,contains("rs")),
  names_from = anc) %>% 
  dplyr::select(1,2,3,contains("CSA"),contains("AFR"),contains("EUR"))
write_csv(res,"/data/home/hmy117/ADAMS_severity/outputs/imsgc_snps_vs_phenotypes.csv")


# power 
calc_power = function(n_wt,n_het,n_hom,beta=0.089){

  iters = list()
  for(i in c(1:1000)){
  # make genos 
  genos = c(rep(0,n_wt),rep(1,n_het),rep(2,n_hom))

  # simulate phenos 
  pheno_wt = rnorm(mean = 0,sd=1,n=n_wt)
  pheno_het = rnorm(mean = beta,sd=1,n=n_het)
  pheno_hom = rnorm(mean = beta*2,sd=1,n=n_hom)

  df = data.frame(genos,pheno = c(pheno_wt,pheno_het,pheno_hom))

  # regress
  iters[[i]] = lm(data = df, pheno ~ genos) %>% broom::tidy()
  }

  # combine 
  iters = do.call("bind_rows",iters)
  iters = iters %>% filter(term == "genos")

  power = iters %>% dplyr::count(hits = p.value < 0.05 & sign(estimate)==1)
  message("power = ",power[power$hits==T,]$n / 1000 * 100,"%")
}

calc_power(163,65,5)
calc_power(119,27,0)
calc_power(205,92,8)

# ld plots

# plots - cohort level
## age vs edss 
png("/data/home/hmy117/ADAMS_severity/plots/rs10191329_vs_age_vs_edss.png",res=900,units="in",width=6,height=6)
ggplot(geno %>% filter(predicted_ancestry %in% c("CSA","EUR","AFR")),aes(ageatedss,edss,color=rs10191329_genotype))+
geom_point(shape=21)+
geom_smooth(se=F,method="lm")+
theme_bw()+
facet_wrap(~predicted_ancestry)+
labs(color="rs10191329 genotype",x="Age",y="EDSS")+
scale_color_brewer(palette="Set1")
dev.off()


## armss
counts = geno %>% dplyr::count(rs10191329_genotype)
medians = geno %>% group_by(rs10191329_genotype) %>% summarise(median_armss = median(gARMSS,na.rm=T)) 

png("/data/home/hmy117/ADAMS_severity/plots/rs10191329_vs_garmss.png",res=900,units="in",width=6,height=6)
ggplot(geno,aes(rs10191329_genotype,gARMSS,fill=rs10191329_genotype))+
geom_violin(alpha=0.5)+
geom_boxplot(width=0.1,alpha=0.3)+
geom_jitter(width=0.1,alpha=0.1)+
theme_bw()+
labs(x="rs10191329 genotype",y="ARMSS at diagnosis")+
scale_fill_brewer(palette="Set1")+
geom_text(data = medians,aes(y=11,label = format(round(median_armss,1),nsmall=1)))+
theme(legend.position="none")
dev.off()

## msis
counts = geno %>% dplyr::count(rs10191329_genotype)
medians = geno %>% group_by(rs10191329_genotype) %>% summarise(median_armss = median(ugARMSS,na.rm=T)) 

png("/data/home/hmy117/ADAMS_severity/plots/rs10191329_vs_msss.png",res=900,units="in",width=6,height=6)
ggplot(geno,aes(rs10191329_genotype,gARMSS,fill=rs10191329_genotype))+
geom_violin(alpha=0.5)+
geom_boxplot(width=0.1,alpha=0.3)+
geom_jitter(width=0.1,alpha=0.1)+
theme_bw()+
labs(x="rs10191329 genotype",y="ARMSS at diagnosis")+
scale_fill_brewer(palette="Set1")+
geom_text(data = medians,aes(y=11,label = format(round(median_armss,1),nsmall=1)))+
theme(legend.position="none")
dev.off()

png("/data/home/hmy117/ADAMS_severity/plots/rs10191329_vs_eq5d.png",res=900,units="in",width=6,height=6)
ggplot(geno,aes(rs10191329_genotype,eq5d_vas))+
geom_boxplot()+
facet_wrap(~predicted_ancestry)
dev.off()


# models
## armss
model_dat = geno %>% filter(!is.na(gARMSS))
model_dat$gARMSS = RNOmni::RankNorm(model_dat$gARMSS)
res = list()
for(ancestry in c("CSA","EUR","AFR")){
  model_dat = geno %>% filter(predicted_ancestry == ancestry )
  res[[length(res)+1]] = lm(data = model_dat, gARMSS ~ ageatedss + Sex + rs10191329_A ) %>% 
    broom::tidy() %>%
    mutate(anc = ancestry)
}
res = do.call("bind_rows",res)
res

## edss
model_dat = geno %>% filter(!is.na(edss))
model_dat$edss = RNOmni::RankNorm(model_dat$edss)
res = list()
for(ancestry in c("CSA","EUR","AFR")){
  model_dat = geno %>% filter(predicted_ancestry == ancestry )
  res[[length(res)+1]] = lm(data = model_dat, edss ~ ageatedss + Sex + rs10191329_A ) %>% 
    broom::tidy() %>%
    mutate(anc = ancestry)
}
res = do.call("bind_rows",res)
res

## msis
model_dat = geno %>% filter(!is.na(msis_physical_normalised))
model_dat$msis_physical_normalised = RNOmni::RankNorm(model_dat$msis_physical_normalised)
res = list()
for(ancestry in c("CSA","EUR","AFR")){
  model_dat = geno %>% filter(predicted_ancestry == ancestry )
  res[[length(res)+1]] = lm(data = model_dat, msis_physical_normalised ~ ageatedss + Sex + rs10191329_A ) %>% 
    broom::tidy() %>%
    mutate(anc = ancestry)
}
res = do.call("bind_rows",res)
res

## eq5d
model_dat = geno %>% filter(!is.na(eq5d_vas))
model_dat$eq5d_vas = RNOmni::RankNorm(model_dat$eq5d_vas)
res = list()
for(ancestry in c("CSA","EUR","AFR")){
  model_dat = geno %>% filter(predicted_ancestry == ancestry )
  res[[length(res)+1]] = lm(data = model_dat, eq5d_vas ~ ageatedss + Sex + rs10191329_A ) %>% 
    broom::tidy() %>%
    mutate(anc = ancestry)
}
res = do.call("bind_rows",res)
res
````

## PCA on whole cohort 
````unix 
cd /data/home/hmy117/ADAMS_severity/
~/plink2 --bfile /data/home/hmy117/ADAMS_severity/combined_adams_imputed \
--indep-pairwise 1000 100 0.1 \
--geno 0.01 \
--maf 0.05 \
--hwe 1e-5 \
--out ./outputs/pruned_snps_for_pca_ALL

# run pca 
~/plink2 --bfile /data/home/hmy117/ADAMS_severity/combined_adams_imputed \
--extract ./outputs/pruned_snps_for_pca_ALL.prune.in \
--pca 10 \
--out ./outputs/pcs_ALL
````

## Prepare covariates 
````R 
library(tidyverse)
setwd("/data/home/hmy117/ADAMS_severity/")
pcs = read_table(paste0("./outputs/pcs_ALL.eigenvec"),col_types="ccdddddddddd")
# read in covars
cov = read_tsv("/data/home/hmy117/ADAMS_severity/pheno/adams_covars.tsv",col_types=cols(.default="c")) %>% dplyr::select(-FID)
pheno = read_tsv("/data/home/hmy117/ADAMS_severity/pheno/adams_pheno.tsv",col_types="ccddddd") %>% dplyr::select(-FID)

# join 
pcs = pcs %>% 
  inner_join(cov,by="IID")

# make cov file 
covars = pcs %>% dplyr::select(1,2,contains("PC"),ageatedss,Sex) 
colnames(covars)[1] = "FID"
# make pheno
pcs = pcs %>% 
  inner_join(pheno,by="IID")

pheno = pcs %>% dplyr::select(1,2,gARMSS,edss,eq5d_vas,age_at_dx,msis_physical_normalised)
colnames(pheno)[1] = "FID"

# save 
write_tsv(pheno,paste0("./outputs/pheno_ALL.tsv"))
write_tsv(covars,paste0("./outputs/cov_ALL.tsv"))

# make tractor pheno 
tractor_pheno = pcs %>% 
  dplyr::select(IID,gARMSS,edss,eq5d_vas,age_at_dx,msis_physical_normalised,ageatedss,Sex,contains("PC"),batch) %>%
  mutate(Sex = ifelse(Sex=="F",2,1))


for(i in c(2:6) ){
  name = colnames(tractor_pheno)[i]
  newname = paste0(name,"_rint")
  x = tractor_pheno %>% 
    dplyr::select(IID,all_of(i)) %>%
    na.omit() 
  x[[newname]] = RNOmni::RankNorm(x[[name]])
  x = x %>% dplyr::select(IID,all_of(newname))
  tractor_pheno <<- tractor_pheno %>% 
    left_join(x,by="IID")
}
tractor_pheno = tractor_pheno %>% 
  dplyr::select(IID,gARMSS_rint,edss_rint,eq5d_vas_rint,age_at_dx_rint,msis_physical_normalised_rint,ageatedss,Sex,batch,contains("PC"))
tractor_pheno$IID = paste0(tractor_pheno$IID,"_",tractor_pheno$IID,"_",tractor_pheno$IID,"_",tractor_pheno$IID)
batch_codex = data.frame(batch = unique(tractor_pheno$batch)) %>%
  mutate(batch_num = row_number())
tractor_pheno = tractor_pheno %>% left_join(batch_codex,by="batch") %>%
  dplyr::select(-batch) %>% dplyr::rename("batch"=batch_num)
# make batch numeric

write_tsv(tractor_pheno,paste0("./outputs/tractor_pheno.tsv"))


````

## Severity GWAS
````unix
cd /data/home/hmy117/ADAMS_severity/
head "/data/home/hmy117/ADAMS_severity/outputs/ancestry_calls.tsv"

for ancestry in CSA EUR AFR; 
do
# filter to ancestry
awk -v anc=$ancestry '{if($3==anc) print $1,$2}' "/data/home/hmy117/ADAMS_severity/outputs/ancestry_calls.tsv" > ./outputs/sample_to_keep_$ancestry

~/plink2 --bfile /data/home/hmy117/ADAMS_severity/combined_adams_imputed \
--keep ./outputs/sample_to_keep_$ancestry \
--out ./outputs/imputed_genotypes_$ancestry \
--make-bed \
--maf 0.01 \
--hwe 1e-10 \
--mind 0.1 \
--geno 0.1

# pca within ancestry 
~/plink2 --bfile ./outputs/imputed_genotypes_$ancestry \
--indep-pairwise 1000 100 0.1 \
--geno 0.01 \
--maf 0.05 \
--hwe 1e-5 \
--out ./outputs/pruned_snps_for_pca_$ancestry

# run pca 
~/plink2 --bfile ./outputs/imputed_genotypes_$ancestry \
--extract ./outputs/pruned_snps_for_pca_$ancestry\.prune.in \
--pca 4 \
--out ./outputs/pcs_$ancestry
done
````

## Make per-ancestry covariate files 
````R 
library(tidyverse)
setwd("/data/home/hmy117/ADAMS_severity/")
# transform pheno 
pheno = read_tsv("/data/home/hmy117/ADAMS_severity/pheno/adams_pheno.tsv",col_types="ccddddd") %>% dplyr::select(-FID)
for(i in c(2:7) ){
  name = colnames(pheno)[i]
  newname = paste0(name,"_rint")
  x = pheno %>% 
    dplyr::select(IID,all_of(i)) %>%
    na.omit() 
  x[[newname]] = RNOmni::RankNorm(x[[name]])
  x = x %>% dplyr::select(IID,all_of(newname))
  pheno <<- pheno %>% 
    left_join(x,by="IID")
}
pheno = pheno %>% dplyr::select(IID,contains("rint"))


for(ancestry in c("CSA","AFR","EUR")){
pcs = read_table(paste0("./outputs/pcs_",ancestry,".eigenvec"),col_types="ccdddd")
# read in covars
cov = read_tsv("/data/home/hmy117/ADAMS_severity/pheno/adams_covars.tsv",col_types=cols(.default="c")) %>% dplyr::select(-FID)

# join 
pcs = pcs %>% 
  inner_join(cov,by="IID")

# make cov file 
covars = pcs %>% dplyr::select(1,2,contains("PC"),ageatedss,Sex,batch) 
colnames(covars)[1] = "FID"
# make pheno
pcs = pcs %>% 
  inner_join(pheno,by="IID")


# save 
write_tsv(pheno,paste0("./outputs/pheno_",ancestry,".tsv"))
write_tsv(covars,paste0("./outputs/cov_",ancestry,".tsv"))


# plink file 
write_tsv(pheno %>% mutate(FID = IID) %>% dplyr::select(FID,IID,everything()),paste0("./outputs/pheno_plink_",ancestry,".tsv"))
covars = covars %>% na.omit()
write_tsv(covars,paste0("./outputs/cov_plink_",ancestry,".tsv"))

}

````

## Run GWAS per-ancestry
````unix
#  qsub plink_gwas.sh
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

````

## PLINK meta-analysis 
````unix 

# format input files 
cd /data/home/hmy117/ADAMS_severity/
awk 'BEGIN {OFS="\t"} NR==1 {print $0, "chr_pos", "A2"; next}{chr_pos = $1 "_" $2;if ($6 == $4) a2 = $5; else a2 = $4;print $0, chr_pos, a2;}' ./outputs/plink_gwas_CSA.gARMSS_rint.glm.linear > ./outputs_plink_gwas_csa_for_meta

awk 'BEGIN {OFS="\t"} NR==1 {print $0, "chr_pos", "A2"; next}{chr_pos = $1 "_" $2;if ($6 == $4) a2 = $5; else a2 = $4;print $0, chr_pos, a2;}' ./outputs/plink_gwas_EUR.gARMSS_rint.glm.linear > ./outputs_plink_gwas_eur_for_meta

awk 'BEGIN {OFS="\t"} NR==1 {print $0, "chr_pos", "A2"; next}{chr_pos = $1 "_" $2;if ($6 == $4) a2 = $5; else a2 = $4;print $0, chr_pos, a2;}' ./outputs/plink_gwas_AFR.gARMSS_rint.glm.linear > ./outputs_plink_gwas_afr_for_meta

# run meta
~/plink \
--meta-analysis \
./outputs_plink_gwas_afr_for_meta \
./outputs_plink_gwas_csa_for_meta \
./outputs_plink_gwas_eur_for_meta + qt weighted-z \
study \
--meta-analysis-snp-field ID \
--meta-analysis-chr-field "#CHROM" \
--meta-analysis-a1-field A1 \
--meta-analysis-a2-field A2 \
--meta-analysis-bp-field POS \
--meta-analysis-se-field SE \
--meta-analysis-ess-field OBS_CT \
--out ./outputs/random_effects_meta

# bring in IMSGC severity GWAS
awk 'BEGIN{OFS="\t"};NR==1{print "#CHROM","POS","chr_pos","A1","A2","BETA","SE","P"};NR>1{print $1,$10,$1"_"$10,$3,$4,$6,$7,$8}' "/data/home/hmy117/ADAMS/genotypes/IMSGC_GWAS/imsgc_mssev_discovery_hg38.tsv" > imsgc_sev_gwas_for_meta

~/plink \
--meta-analysis \
./outputs_plink_gwas_afr_for_meta \
./outputs_plink_gwas_csa_for_meta \
./outputs_plink_gwas_eur_for_meta \
imsgc_sev_gwas_for_meta + qt \
study report-all \
--meta-analysis-snp-field chr_pos \
--meta-analysis-chr-field "#CHROM" \
--meta-analysis-bp-field POS \
--out ./outputs/random_effects_meta_with_imsgc

````

## VEP 
````unix 

# get significant snps 
cd /data/home/hmy117/ADAMS_severity/
awk '{if($12<0.00005) print $3}' ./outputs/plink_gwas_CSA.gARMSS_rint.glm.linear > sigsnps
awk '{if($12<0.00005) print $3}' ./outputs/plink_gwas_EUR.gARMSS_rint.glm.linear >> sigsnps
awk '{if($12<0.00005) print $3}' ./outputs/plink_gwas_AFR.gARMSS_rint.glm.linear >> sigsnps
awk 'NR==1{print};NR>1{if($8<0.00005) print $3}' ./outputs/random_effects_meta.meta >> sigsnps


# repeat for imsgc
awk '{if($8 < 0.000005) print $1":"$2":"$4":"$5}' /data/home/hmy117/ADAMS_severity/imsgc_sev_gwas_for_meta > imsgc_sig_snps
awk '{if($8 < 0.000005) print $1":"$2":"$5":"$4}' /data/home/hmy117/ADAMS_severity/imsgc_sev_gwas_for_meta >> imsgc_sig_snps

# make one-sample vcf
head -n 1 /data/home/hmy117/ADAMS_severity/combined_adams_imputed.fam > sample_to_keep 
~/plink2 --bfile /data/home/hmy117/ADAMS_severity/combined_adams_imputed \
--extract sigsnps \
--keep sample_to_keep \
--export vcf \
--out site_only_vcf

head -n 1 /data/home/hmy117/ADAMS_severity/combined_adams_imputed.fam > sample_to_keep 
~/plink2 --bfile /data/home/hmy117/ADAMS_severity/combined_adams_imputed \
--extract imsgc_sig_snps \
--keep sample_to_keep \
--export vcf \
--out site_only_vcf_imsgc

module unload R
module load miniforge
mamba activate vep_env

# annotate with vep
vep -i site_only_vcf.vcf \
-o /data/home/hmy117/ADAMS_severity/outputs/snp_annotations \
--cache \
--dir_cache /data/scratch/hmy117/.vep \
--force_overwrite \
--nearest symbol \
--everything \
--pick \
--tab

# annotate with vep - imsgc
cd /data/home/hmy117/ADAMS_severity/
vep -i site_only_vcf_imsgc.vcf \
-o /data/home/hmy117/ADAMS_severity/outputs/snp_annotations_imsgc \
--cache \
--dir_cache /data/scratch/hmy117/.vep \
--force_overwrite \
--nearest symbol \
--pick \
--everything \
--tab
````

## Plots 
````R 
library(tidyverse)
setwd("/data/home/hmy117/ADAMS_severity/")

# read in 
eur = read_table("./outputs/plink_gwas_EUR.gARMSS_rint.glm.linear") %>% mutate(ancestry="EUR")%>%dplyr::rename("CHROM" = `#CHROM`,"GENPOS"=POS,"SNP"=ID) 
sas = read_table("./outputs/plink_gwas_CSA.gARMSS_rint.glm.linear")%>% mutate(ancestry="CSA")%>%dplyr::rename("CHROM" = `#CHROM`,"GENPOS"=POS,"SNP"=ID) 
afr = read_table("./outputs/plink_gwas_AFR.gARMSS_rint.glm.linear")%>% mutate(ancestry="AFR")%>%dplyr::rename("CHROM" = `#CHROM`,"GENPOS"=POS,"SNP"=ID) 
meta = read_table("./outputs/random_effects_meta.meta")%>% mutate(ancestry="META")%>%
dplyr::select(-BETA,-P) %>%
dplyr::rename("CHROM" = CHR,"GENPOS"=BP,"P"="P(R)","BETA"="BETA(R)") 

combo = bind_rows(eur,sas,afr,meta)

# compare with imsgc severity gwas
imsgc_sev = read_table("imsgc_sev_gwas_for_meta")
imsgc_sev_sig =  imsgc_sev %>% filter(P < 5e-6)

combo_imsgc = combo %>% 
  dplyr::select(CHROM,GENPOS,A1,REF,ALT,BETA,SE,P,ancestry) %>% 
  filter(ancestry != "Meta-analysis") %>%
  mutate(chr_pos = paste0(CHROM,"_",GENPOS)) %>% 
  filter(chr_pos %in% imsgc_sev_sig$chr_pos) %>% 
  pivot_wider(id_cols = c(chr_pos,A1,REF,ALT),names_from = ancestry, values_from = c(BETA,SE,P))
combo_imsgc %>% arrange(chr_pos)
imsgc_hits = imsgc_sev_sig %>%
  inner_join(combo_imsgc,by="chr_pos") %>% 
  arrange(`#CHROM`,POS)

# add annotations
anno = read_table("/data/home/hmy117/ADAMS_severity/outputs/snp_annotations_imsgc",skip=108,col_types=cols(.default="c")) %>% 
  filter(CANONICAL=="YES" & MANE == "MANE_Select") %>% 
  dplyr::rename("SNP" = `#Uploaded_variation`)

imsgc_hits = imsgc_hits %>%
  mutate(Location = paste0(`#CHROM`,":",POS)) %>% 
  left_join(anno,by="Location")

# concordance
imsgc_hits = imsgc_hits %>% 
  mutate(imsgc_beta_orientated_to_adams = ifelse(A1.x == A1.y,BETA,BETA*-1)) %>%
  mutate(concordant_AFR = ifelse(sign(imsgc_beta_orientated_to_adams) == sign(BETA_AFR),"Yes","No")) %>%
  mutate(concordant_CSA = ifelse(sign(imsgc_beta_orientated_to_adams) == sign(BETA_CSA),"Yes","No")) %>%
  mutate(concordant_EUR = ifelse(sign(imsgc_beta_orientated_to_adams) == sign(BETA_EUR),"Yes","No")) %>%
  mutate(replicated_AFR = ifelse(concordant_AFR == "Yes" & P_AFR < 0.05,"*","")) %>%
  mutate(replicated_CSA = ifelse(concordant_CSA == "Yes" & P_CSA < 0.05,"*","")) %>%
  mutate(replicated_EUR = ifelse(concordant_EUR == "Yes" & P_EUR < 0.05,"*","")) 


write.csv(imsgc_hits,"/data/home/hmy117/ADAMS_severity/outputs/imsgc_hits.csv")  

# add annotations
anno = read_table("/data/home/hmy117/ADAMS_severity/outputs/snp_annotations",skip=108) %>% 
  dplyr::rename("SNP" = `#Uploaded_variation`)

# plot 
dat = combo 
sig_hits = dat %>% filter(P < 5e-6)
sig_hits = sig_hits  %>% left_join(anno,by="SNP")

sig_hits = sig_hits %>% left_join(imsgc_sev %>% dplyr::rename("CHROM"=`#CHROM`,"GENPOS"=POS),by=c("CHROM","GENPOS"))

# save
write.csv(sig_hits,"/data/home/hmy117/ADAMS_severity/outputs/garmss_suggestive_hits.csv")  


# manhattan
dat = dat %>%
    mutate(colcode = case_when(
        P < 5e-6 ~ "sig",
        P >=5e-6 & CHROM %%2 == 0 ~ "even",
        P >=5e-6 & CHROM %%2 != 0 ~ "odd"
    ))

pal = c("blue","lavenderblush1","lavenderblush2")
names(pal) = c("sig","even","odd")

# define windows
chrcoords = dat %>% group_by(CHROM) %>% summarise(min_bp = min(GENPOS), max_bp = max(GENPOS)) %>%
    mutate(CHROM = CHROM + 1) %>%
    mutate(cumbp_total = cumsum(max_bp))
dat = dat %>%
    left_join(chrcoords,by="CHROM") %>%
    mutate(cumbp = ifelse(is.na(cumbp_total),GENPOS,GENPOS+cumbp_total))
midpoints = dat %>%
    group_by(CHROM) %>%
    summarise(midpoint = median(cumbp))
plot_dat = dat %>% filter(P<0.05)
plot_dat = plot_dat %>% left_join(anno %>% dplyr::select(SNP,NEAREST),by="SNP")

p = ggplot(plot_dat,aes(cumbp,-log10(P),color=colcode))+
    geom_point(data = plot_dat %>% filter(P < 0.1 & colcode != "sig"),alpha=0.3)+
    geom_point(data = plot_dat %>% filter(P < 0.1 & colcode == "sig" & P >1e-5),alpha=0.5)+
    geom_point(data = plot_dat %>% filter(P < 0.1 & colcode == "sig" & P <=1e-5),alpha=1)+
    scale_x_continuous(breaks = midpoints$midpoint,labels = midpoints$CHROM)+
    scale_color_manual(values = pal)+
    ggrepel::geom_text_repel(data = plot_dat %>% filter(P < 5e-6) %>% group_by(ancestry,CHROM) %>% slice_min(P,n=1),
    mapping = aes(label = NEAREST),min.segment.length=0,nudge_y=1,nudge_x=1000,color="black",segment.linetype = 5,segment.alpha=0.6)+
    theme_bw()+
    facet_wrap(~ancestry,nrow=4)+
    theme(legend.position="none")+
    labs(y=bquote(log[10]~P),x="Genomic co-ordinates")+
    geom_hline(yintercept=-log10(5e-6),linetype="dashed",alpha=0.5,color="blue")+
    geom_hline(yintercept=-log10(5e-8),linetype="dashed",alpha=0.5,color="red")

outplot = paste0(
    "/data/home/hmy117/ADAMS_severity/plots/manhattan.png"
)

png(outplot,units="in",res=900,width=10,height=10)
p
dev.off()

````


## RFMIX 
````unix
# mkdir /data/scratch/hmy117/rfmix
cd /data/scratch/hmy117/rfmix

# make sample ref 
cut -f1,177 /data/scratch/hmy117/aj_ref/gnomad_meta_v1.tsv > sample_ref

# get map 
wget https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz
gunzip genetic_map_hg38_withX.txt.gz
awk '{print $1,$2,$4}' genetic_map_hg38_withX.txt > genetic_map

cd /data/scratch/hmy117/rfmix

qsub /data/home/hmy117/ADAMS_severity/scripts/rfmix.sh 
````

## Tractor 
### Setup
````unix
# git clone https://github.com/Atkinson-Lab/Tractor.git
module load miniforge
mamba create -n tractor
mamba activate tractor

# extract haps
cd /data/scratch/hmy117/rfmix
qsub /data/home/hmy117/ADAMS_severity/scripts/tractor.sh
````

### Run tractor gwas
module unload miniforge
module load R/4.4.1

qsub  /data/home/hmy117/ADAMS_severity/scripts/tractor_gwas.sh 



### Explore in R
````R

library(tidyverse)

setwd("/data/scratch/hmy117/tractor/")

# compare ancestry calls vs admixture
rfmix = read_table("/data/scratch/hmy117/rfmix/rfmix_local_ancestry_chr1.rfmix.Q",skip=1)
# compare with ADMIXTURE 
q = read_table("/data/scratch/hmy117/hgdp_1kg_genomes/combined_adams_imputed_qc_pruned_merge_for_admixture.8.Q",col_names=F)
pop = read_table("/data/scratch/hmy117/adams_imputed_severity_topmed/combined_adams_imputed_qc_pruned_merge_for_admixture.pop",col_names=F) %>% 
filter(X1!="-") %>% 
distinct()
colnames(q) = pop$X1
fam = read_table("/data/scratch/hmy117/adams_imputed_severity_topmed/combined_adams_imputed_qc_pruned_merge_for_admixture.fam",col_names=F)
colnames(fam)[2]="IID"
q = q %>% 
mutate(IID = fam$IID) 
colnames(rfmix)[1]="IID"
rfmix = rfmix %>% pivot_longer(cols = -IID) %>% dplyr::rename("Ancestry"=name,"Proportion (RFMIX)"=value)
dat = q %>% pivot_longer(cols = -IID) %>% dplyr::rename("Ancestry"=name,"Proportion (ADMIXTURE)"=value) %>%
  mutate(IID = paste0(IID,"_",IID,"_",IID,"_",IID)) %>%
  left_join(rfmix,by=c("IID","Ancestry")) %>% 
  filter(Ancestry %in% c("CSA","AFR","EUR"))
png("/data/home/hmy117/ADAMS_severity/plots/lai_vs_admixture.png",units="in",res=900,width=6,height=3)
ggplot(dat,aes(`Proportion (ADMIXTURE)`,`Proportion (RFMIX)`))+
  geom_point()+
  facet_wrap(~Ancestry,nrow=1)+
  theme_bw()+
  geom_smooth(method="lm",linetype="dashed")+
  geom_abline(color="red",intercept=0,slope=1,alpha=0.3)
dev.off()


# armss
res = list()
for(i in c(1:22)){
  res[[i]] = read_tsv(paste0("/data/scratch/hmy117/tractor/tractor_gwas_edss_rint_",i))

}
res = do.call("bind_rows",res)

# fix colnames
colnames(res)  = str_replace_all(
  colnames(res),
  "anc0","AFR"
)

colnames(res) = str_replace_all(
  colnames(res),
  "anc1","AMR"
)

colnames(res) = str_replace_all(
  colnames(res),
  "anc2","CSA"
)

colnames(res) = str_replace_all(
  colnames(res),
  "anc3","EAS"
)

colnames(res) = str_replace_all(
  colnames(res),
  "anc4","EUR"
)

colnames(res) = str_replace_all(
  colnames(res),
  "anc5","MID"
)

colnames(res) = str_replace_all(
  colnames(res),
  "anc6","OCE"
)

# lengthen 
pval_res = res %>%
  pivot_longer(cols =c(7:ncol(res))) %>%
  filter(grepl("pval",name))

# 
pval_res = pval_res %>%
  separate(name,sep="_",into=c("variable","ancestry"))

# filter 
dat = pval_res %>% filter(value < 0.05) %>%
  dplyr::rename("CHROM"=CHR,"GENPOS"=POS,"P"=name)


# filter for plot
dat = dat %>%
    mutate(colcode = case_when(
        P < 5e-6 ~ "sig",
        P >=5e-6 & CHROM %%2 == 0 ~ "even",
        P >=5e-6 & CHROM %%2 != 0 ~ "odd"
    ))

pal = c("blue","lavenderblush1","lavenderblush2")
names(pal) = c("sig","even","odd")

# define windows
chrcoords = dat %>% group_by(CHROM) %>% summarise(min_bp = min(GENPOS), max_bp = max(GENPOS)) %>%
    mutate(CHROM = CHROM + 1) %>%
    mutate(cumbp_total = cumsum(max_bp))
dat = dat %>%
    left_join(chrcoords,by="CHROM") %>%
    mutate(cumbp = ifelse(is.na(cumbp_total),GENPOS,GENPOS+cumbp_total))
midpoints = dat %>%
    group_by(CHROM) %>%
    summarise(midpoint = median(cumbp))
plot_dat = dat %>% filter(P<0.05)
plot_dat = plot_dat %>% left_join(anno %>% dplyr::select(SNP,NEAREST),by="SNP")

p = ggplot(plot_dat,aes(cumbp,-log10(P),color=colcode))+
    geom_point(data = plot_dat %>% filter(P < 0.1 & colcode != "sig"),alpha=0.3)+
    geom_point(data = plot_dat %>% filter(P < 0.1 & colcode == "sig" & P >1e-5),alpha=0.5)+
    geom_point(data = plot_dat %>% filter(P < 0.1 & colcode == "sig" & P <=1e-5),alpha=1)+
    scale_x_continuous(breaks = midpoints$midpoint,labels = midpoints$CHROM)+
    scale_color_manual(values = pal)+
    ggrepel::geom_text_repel(data = plot_dat %>% filter(P < 5e-6) %>% group_by(ancestry,CHROM) %>% slice_min(P,n=1),
    mapping = aes(label = NEAREST),min.segment.length=0,nudge_y=1,nudge_x=1000,color="black",segment.linetype = 5,segment.alpha=0.6)+
    theme_bw()+
    facet_wrap(~ancestry,nrow=4)+
    theme(legend.position="none")+
    labs(y=bquote(log[10]~P),x="Genomic co-ordinates")+
    geom_hline(yintercept=-log10(5e-6),linetype="dashed",alpha=0.5,color="blue")+
    geom_hline(yintercept=-log10(5e-8),linetype="dashed",alpha=0.5,color="red")

outplot = paste0(
    "/data/home/hmy117/ADAMS_severity/plots/manhattan.png"
)

png(outplot,units="in",res=900,width=10,height=10)
p
dev.off()


# genome_wide plots
dat = read_tsv(paste0("/data/scratch/hmy117/rfmix/rfmix_local_ancestry_chr",i,".msp.tsv"),skip=1)
anc_labels = read_tsv(paste0("/data/scratch/hmy117/rfmix/rfmix_local_ancestry_chr",i,".msp.tsv"),n_max=1)


anc = read_tsv("/data/home/hmy117/ADAMS_severity/outputs/ancestry_calls.tsv")


# find afr-eur admixed 
admixed =  q %>% filter(AFR > 0.4 & EUR > 0.4)

# lengthen
dat = dat %>% 
  pivot_longer(c(7:ncol(dat))) %>%
  mutate(
    ancestry = case_when(
      value == 0 ~ "AFR",
      value == 1 ~ "AMR",
      value == 2 ~ "CSA",
      value == 3 ~ "EAS",
      value == 4 ~ "EUR",
      value == 5 ~ "MID",
      value == 6 ~ "OCE" 
      )
  ) 
  
dat = dat  %>%
  separate(name,sep="_",into=c("IID","IID2","IID3","fullhap")) %>%
  separate(fullhap,sep="\\.",into=c("trash","hap"))

# compare local and global ancestry 
q_long = q %>% pivot_longer(-IID)
local_ancestry_props = dat %>% group_by(IID,hap) %>% dplyr::count(ancestry) %>% mutate(prop = n/sum(n))
local_ancestry_props = local_ancestry_props %>%
  left_join(q_long,by="IID")

png("/data/home/hmy117/ADAMS_severity/plots/lai_calls_vs_global_anc.png",res=900,units="in",width=8,height=8)
ggplot(local_ancestry_props,aes(prop,value,color=hap))+facet_grid(ancestry~name)+geom_point()
dev.off()

i=1

png("/data/home/hmy117/ADAMS_severity/plots/lai_calls_admixed.png",res=900,units="in",width=8,height=8)
ggplot(dat %>% filter(IID %in% admixed$IID),aes(spos,ancestry,color=hap,group=as.character(hap)))+geom_line()+facet_wrap(~IID)
dev.off()



anc_labels = read_tsv(paste0("/data/scratch/hmy117/rfmix/rfmix_local_ancestry_chr",i,".msp.tsv"),n_max=1)


# read vcf directly 
vcf = read_table("rs10191329_vcf") 
vcf = vcf %>%
  pivot_longer(cols = c(10:ncol(vcf))) %>%
  dplyr::select(-INFO,-QUAL,-FILTER,-FORMAT) %>%
  separate(name,sep="_",into=c("IID","FID")) %>%
  separate(value,sep="\\|",into=c("0","1")) %>%
  pivot_longer(cols = c("0","1"))  %>%
  dplyr::rename("hap"=name) %>%
  mutate(allele=ifelse(value==0,"C","A"))

vcf = vcf %>% left_join(dat,by=c("IID","hap")) %>%
dplyr::select(IID,hap,allele,ancestry)

# freqs 
vcf %>% group_by(ancestry) %>% dplyr::count(allele) %>%
mutate(prop = n/sum(n))

````


## SUSIE
````unix
module load R/4.4.1

# mkdir 
# mkdir /data/scratch/hmy117/susiex/
cd /data/scratch/hmy117/susiex/


# get plink reference from https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg
# harmonise ref & gwas
~/plink2 --vcf /data/scratch/hmy117/hgdp_1kg_genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr2.vcf.bgz \
--set-all-var-ids @:#:\$r\:\$a \
--new-id-max-allele-len 9999 \
--extract-if-info 'AF > 0.05' \
--make-bed \
--geno 0.1 \
--maf 0.05 \
--hwe 1e-20 \
--mac 1 \
--from-bp 70449869 \
--to-bp 72449869 \
--chr 2 \
--out reference_kg_hg38_chr2

# get ancestry labels & filter to adams snps
awk 'BEGIN{FS="\t"};{print $1,$177}' /data/scratch/hmy117/aj_ref/gnomad_meta_v1.tsv | grep EUR | awk '{print 0,$1}' > eur_ids.txt
awk 'BEGIN{FS="\t"};{print $1,$177}' /data/scratch/hmy117/aj_ref/gnomad_meta_v1.tsv | grep CSA  | awk '{print 0,$1}' > csa_ids.txt
awk 'BEGIN{FS="\t"};{print $1,$177}' /data/scratch/hmy117/aj_ref/gnomad_meta_v1.tsv | grep AFR | awk '{print 0,$1}'  > afr_ids.txt

# filter gwas to these snps 
awk 'BEGIN{OFS="\t"};NR==1{print $0,"cpra"};NR>1{if($1==2 && $2 > 70449869 && $2 < 72449869) print $0,$1":"$2":"$6":"$15}' ~/ADAMS_severity/outputs_plink_gwas_csa_for_meta > csa_gwas_for_susie
awk 'BEGIN{OFS="\t"};NR==1{print $0,"cpra"};NR>1{if($1==2 && $2 > 70449869 && $2 < 72449869) print $0,$1":"$2":"$6":"$15 }' ~/ADAMS_severity/outputs_plink_gwas_afr_for_meta > afr_gwas_for_susie
awk 'BEGIN{OFS="\t"};NR==1{print $0,"cpra"};NR>1{if($1==2 && $2 > 70449869 && $2 < 72449869) print $0,$1":"$2":"$6":"$15 }' ~/ADAMS_severity/outputs_plink_gwas_eur_for_meta > eur_gwas_for_susie
awk 'BEGIN{OFS="\t"};NR==1{print $0,"ID"};NR>1{if($1==2 && $2 > 70449869 && $2 < 72449869) print $0,$1":"$2":"$4":"$5}' ~/ADAMS_severity/imsgc_sev_gwas_for_meta > imsgc_eur_gwas_for_susie



# filter    
~/plink2 --bfile reference_kg_hg38_chr2 \
--keep eur_ids.txt \
--set-all-var-ids @:#:\$a\:\$r \
--new-id-max-allele-len 9999 \
--make-bed \
--out g1000_eur_hg38_chrpos_chr2

~/plink2 --bfile reference_kg_hg38_chr2 \
--set-all-var-ids @:#:\$a\:\$r \
--new-id-max-allele-len 9999 \
--keep csa_ids.txt \
--make-bed \
--out g1000_csa_hg38_chrpos_chr2    

~/plink2 --bfile reference_kg_hg38_chr2 \
--keep afr_ids.txt \
--set-all-var-ids @:#:\$a\:\$r \
--new-id-max-allele-len 9999 \
--make-bed \
--out g1000_afr_hg38_chrpos_chr2



# run susiex
/data/home/hmy117/SuSiEx/bin_static/SuSiEx \
--sst_file=csa_gwas_for_susie,afr_gwas_for_susie,eur_gwas_for_susie,imsgc_eur_gwas_for_susie \
--n_gwas=230,141,304,12584 \
--ref_file=g1000_csa_hg38_chrpos_chr2,g1000_afr_hg38_chrpos_chr2,g1000_eur_hg38_chrpos_chr2,g1000_eur_hg38_chrpos_chr2 \
--ld_file=csa_ld,afr_ld,eur_ld,eur_ld \
--out_dir=/data/scratch/hmy117/susiex/ \
--out_name=dysferlin \
--chr=2 \
--bp=70449869,72449869 \
--chr_col=1,1,1,1 \
--snp_col=16,16,16,9 \
--bp_col=2,2,2,2 \
--a1_col=6,6,6,4 \
--a2_col=15,15,15,5 \
--eff_col=9,9,9,6 \
--se_col=10,10,10,7 \
--pval_col=12,12,12,8 \
--plink=~/plink \
--mult-step=T \
--threads=$NSLOTS \
--n_sig=10 \
--maf 0.05 \
--pval_thresh=1 \
--max_iter=100
````

## Plot fine-mapping results
````R 
library(tidyverse)
dat = read_table("/data/scratch/hmy117/susiex/dysferlin.snp") 

p = ggplot(dat %>% filter(BP > 70949869 & BP < 71949869),aes(BP,`PIP(CS1)`,fill=`PIP(CS1)`,label=SNP))+
ggrepel::geom_text_repel(data = dat %>% filter(`PIP(CS1)`>0.05 & BP > 70949869 & BP < 71949869),min.segment.length=0)+
theme_bw()+
geom_point(shape=21,size=2)+
scale_fill_viridis_c(option="plasma")+
labs(x="Position on chromosome 2",y="PIP")+
theme(legend.position="none")+

png("/data/home/hmy117/ADAMS_severity/plots/rs10191329_finemapping.png",res=900,units="in",width=6,height=4)
p
dev.off()


````

## GWAS catalogue prep
````unix
~/plink2 --bfile ~/ADAMS_severity/outputs/imputed_genotypes_AFR --freq --out ~/ADAMS_severity/outputs/freqs_AFR
~/plink2 --bfile ~/ADAMS_severity/outputs/imputed_genotypes_CSA --freq --out ~/ADAMS_severity/outputs/freqs_CSA
~/plink2 --bfile ~/ADAMS_severity/outputs/imputed_genotypes_AFR --freq --out ~/ADAMS_severity/outputs/freqs_EUR

Rscript gwas_catalogue_prep.R


````