#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readr)
library(glue)
library(purrr)
library(data.table)
args <- commandArgs(trailingOnly = TRUE)

# Arguments
# [1]: LOCO mapping
# [2]: INBRED mapping
# args <- c("data/10_2_0.6_0.05_gamma_granitic_esok_lmm-exact.loco.mlma",
#           "data/10_2_0.6_0.05_gamma_granitic_esok_lmm-exact_inbred.fastGWA")
LOCO <- data.table::fread(args[1]) %>%
  dplyr::mutate(algorithm = "LOCO")
INBRED <- data.table::fread(args[2]) %>%
  dplyr::mutate(algorithm = "INBRED") %>%
  dplyr::select(-N)
colnames(LOCO) <- colnames(INBRED)

# chr.check <- FALSE %in% (LOCO$CHR == INBRED$CHR)
# snp.check <- FALSE %in% (LOCO$SNP== INBRED$SNP)
# pos.check <- FALSE %in% (LOCO$POS == INBRED$POS)
# a1check <- FALSE %in% (LOCO$A1 == INBRED$A1)
# a2check <- FALSE %in% (LOCO$A2 == INBRED$A2)
# af1check <- FALSE %in% (LOCO$AF1 == INBRED$AF1)
# nrow.check <- FALSE %in% nrow(LOCO) == nrow(INBRED)
# all.checks <- unlist(lapply(ls(pattern = "check"), get))


if(nrow(INBRED) > nrow(LOCO)){
  INBRED <- INBRED[which(INBRED$SNP %in% LOCO$SNP),]
  base <- INBRED %>%
    dplyr::select(-BETA,-SE,-P,-algorithm)
} else {
  LOCO <- LOCO[which(LOCO$SNP %in% INBRED$SNP),]
  base <- LOCO %>%
    dplyr::select(-BETA,-SE,-P,-algorithm)
}

inbred <- INBRED %>%
    dplyr::select(SNP, BETA, SE, P, algorithm)
  
loco <- LOCO %>%
    dplyr::select(SNP, BETA, SE, P, algorithm)
  
aggregate.mapping <- base %>%
    dplyr::mutate(BETA = if_else(condition = inbred$P < loco$P, 
                                 true = inbred$BETA, 
                                 false = loco$BETA),
                  SE = if_else(condition = inbred$P < loco$P, 
                               true = inbred$SE,
                               false = loco$SE),
                  P = if_else(condition = inbred$P < loco$P, 
                              true = inbred$P, 
                              false = loco$P),
                  algorithm = if_else(condition = inbred$P < loco$P,
                                      true = inbred$algorithm, 
                                      false = loco$algorithm))

readr::write_tsv(x = aggregate.mapping, "temp.aggregate.mapping.tsv")