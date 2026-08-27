library(tidyverse)

for(x in c("AFR","CSA","EUR")){
file = paste0("/data/home/hmy117/ADAMS_severity/outputs/plink_gwas_",x,".gARMSS_rint.glm.linear")
outfile = paste0("/data/home/hmy117/ADAMS_severity/outputs/sumstats_for_gwas_catalogue_ARMSS_",x,".tsv")
read_table(file) %>% 
  dplyr::rename("chromosome" = '#CHROM',
                      "base_pair_location" = POS,
                      "effect_allele" = A1) %>%
  mutate(other_allele = ifelse(REF == effect_allele,ALT,REF)) %>%
  dplyr::rename("beta" = BETA,"standard_error"=SE,"p_value" = P,"variant_id" = ID,"n" = OBS_CT,ref_allele = REF) %>%
  dplyr::select(chromosome,base_pair_location,effect_allele,other_allele,beta,standard_error,p_value,variant_id,n,ref_allele) %>%
  write_tsv(outfile)
}