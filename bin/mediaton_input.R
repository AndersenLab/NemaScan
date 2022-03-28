#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readr)
library(glue)
library(purrr)


# load arguments
args <- commandArgs(trailingOnly = TRUE)


gtrait = args[1]

# load phenotpe data and scale
phenotype_data <- readr::read_tsv(args[2]) %>%
  na.omit() %>%
  as.data.frame() %>% 
  dplyr::select(strain,ph=gtrait) %>%
  dplyr::mutate(newpheno = (ph - mean(ph, na.rm = T)) / sd(ph, na.rm = T)) %>%
  dplyr::select(strain, tr = newpheno)

# save mapping data set
readr::write_tsv(phenotype_data, 
                 path = glue::glue("{gtrait}_scaled_mapping_{args[8]}.tsv"),
                 col_names = T)


gwas_pchr = args[3] 

gwas_start = args[4] %>% as.numeric()

gwas_stop = args[5] %>% as.numeric()

gwas_p = args[6] %>% as.numeric()

 
eqtl_transcript <- read.delim(args[7], stringsAsFactors=FALSE)

# eQTL overlap trait GWA QTL
gene_qtl <- eqtl_transcript %>% 
  dplyr::filter(e_chr==gwas_pchr) %>% 
  dplyr::filter(!e_end<=(gwas_start-1e6)) %>% 
  dplyr::filter(!e_start>=(gwas_stop+1e6)) %>% 
  dplyr::mutate(gwtrait=gtrait, gwpeak=gwas_p)  

# save mapping data set
readr::write_tsv(gene_qtl, 
                 path = glue::glue("{gtrait}_{gwas_pchr}_{gwas_p}_eqtl_{args[8]}.tsv"),
                 col_names = T)

 