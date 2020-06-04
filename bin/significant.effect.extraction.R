#!/usr/bin/env Rscript
require(tidyverse)
require(GenomicRanges)
args <- commandArgs(trailingOnly = TRUE)
setwd(paste(args[1],"Simulations",sep = "/"))
# setwd("~/Dropbox/AndersenLab/LabFolders/Sam/projects/NemaScan/")
effect.files <- list.files(pattern = ".par",recursive = T)
iterations <- purrr::map(effect.files, .f = function(x){
  paste(strsplit(strsplit(x,split = "/")[[1]][3],split = "_")[[1]][1:3], collapse = "_")
  # paste(strsplit(x,split = "_")[[1]][1:3], collapse = "_")
})
extract.effect.metrics <- function(x){
  print(x)
  effects <- read.table(paste(paste(strsplit(x, split = "_")[[1]][1],
                                    "Phenotypes", sep = "/"),
                              paste(x,"sims.par",sep = "_"), sep = "/"),
                        header = T) %>%
  # effects <- read.table(paste(x,"sims.par",sep = "_"), header = T) %>%
    tidyr::separate(col = QTL,
                    into = c("CHR","POS"),
                    remove = FALSE) %>%
    dplyr::mutate(SNP = QTL)

  
  fastGWA.inbred <- read.table(paste(paste(strsplit(x, split = "_")[[1]][1],
                                           "Mappings", sep = "/"),
                                     paste(x,"lmm-exact_inbred.fastGWA",sep = "_"), sep = "/"),
                               header = T) %>%
  # fastGWA.inbred <- read.table(paste(x,"lmm-exact_inbred.fastGWA",sep = "_"),
  #                                header = T) %>%
    dplyr::mutate(CHR = as.factor(CHR)) %>%
    dplyr::mutate(true.QTL = as.factor(SNP %in% effects$QTL)) %>%
    dplyr::filter(CHR != "Mt") %>%
    dplyr::mutate(log10p = -log(P)) %>%
    dplyr::mutate(sig = log10p > 7) %>%
    dplyr::select(-c("A1","A2","N","P")) %>%
    dplyr::mutate(algo = "fastGWA.inbred")
  
  merge(effects, fastGWA.inbred, "SNP") %>%
    dplyr::select(SNP, Frequency, Effect, CHR.x, POS.x, BETA, log10p, sig) %>%
    dplyr::mutate(trait = x) %>%
    tidyr::separate(trait, into = c("nQTL","rep","h2"), sep = "_")

}
effects <- purrr::map(iterations, extract.effect.metrics)
effects.extracted <- do.call(rbind, effects)
save(effects.extracted,
     file = "Extracted_Effects.RData")
