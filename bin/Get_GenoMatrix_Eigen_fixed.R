library(Rcpp)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readr)
library(glue)
library(purrr)
# library(correlateR)
library(coop)

library(parallel) #mjs

args <- commandArgs(trailingOnly = TRUE)

phenotyped_strain_snps <- readr::read_tsv(args[1]) %>%
  na.omit()

analysis_chromosome <- args[2]

strain_count <- (ncol(phenotyped_strain_snps) - 4)

chrom_geno <- list()
for(chrom in 1:length(unique(phenotyped_strain_snps$CHROM))){
  t_chrom <- unique(phenotyped_strain_snps$CHROM)[chrom]
  t_df <- dplyr::filter(phenotyped_strain_snps, CHROM == t_chrom)
  t_df <- t_df[,5:ncol(t_df)]
  
  keepMarkers <- data.frame(
    MAF = apply(t_df,
                MARGIN = 1,
                FUN = function(x){
                  x[x==-1] <- 0
                  return(sum(x, na.rm = T)/length(x))}
    )) 
  
  t_df <- dplyr::mutate(t_df, MAF = keepMarkers$MAF) %>%
    dplyr::filter(MAF >= 0.05, MAF <= 0.95) %>%
    dplyr::select(-MAF)
  
  chrom_geno[[chrom]] <- t_df
}

lapply(chrom_geno, nrow)

eigenvalues <- list()
for(chrom in 1:length(chrom_geno)){
  snpcor <- coop::pcor(t(chrom_geno[[chrom]]))
  snpeigen <- eigen(snpcor, nrow(snpcor), only.values=TRUE)
  snpeigen$values = abs(snpeigen$values)
  whole_snpeigen = ifelse(snpeigen$values > 1, 1, 0)
  partial_snpeigen = snpeigen$values - floor(snpeigen$values)
  eigenvalues[[chrom]] = sum(whole_snpeigen) + sum(partial_snpeigen)
} # end chrom loop 

evals <- unlist(eigenvalues)

independent_tests <- data.frame(independent_tests = sum(evals))

readr::write_csv(independent_tests, path = glue::glue("{analysis_chromosome}_independent_snvs.csv"))
