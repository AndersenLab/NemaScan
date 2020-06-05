#!/usr/bin/env Rscript
library(tidyverse)
library(reshape2)

simulations <- list.files(pattern="HSQ")
# Concatenation Functions
tidy.effects <- function(x){
      print(paste("Beginning",x,sep = " "))
      setwd(paste("/projects/b1059/projects/Sam/GWAS_sensitivity/",x, sep = ""))
      hsq <- stringr::str_split(x, pattern = "[:alpha:]")[[1]][4]
      maf <- stringr::str_split(x, pattern = "[:alpha:]")[[1]][7]
      LD <- stringr::str_split(x, pattern = "[:alpha:]")[[1]][9]
      
      architectures <- list.files(pattern="QTL")
      all.architectures <- list()
      for(i in 1:length(architectures)){
            qtl.effect.files <- list.files(path = paste(architectures[[i]],
                                                        "causal.QTLs",
                                                        sep = "/"))
            qtl.effect.data <- list()
            print(paste("Beginning", architectures[[i]]))
            for(j in 1:length(qtl.effect.files)){
                  A <- read.table(paste(getwd(), 
                                        architectures[[i]],
                                        "causal.QTLs",
                                        paste("causal.QTLs.rep.", j, sep = ""),
                                        sep = "/"), header = T)
                  
                  A$Rep <- paste("sim", j, sep = "_")
                  A$nQTL <- strsplit(stringr::str_split(architectures[[i]],pattern = "[:alpha:]")[[1]][1], 
                                     split = "[.]")[[1]]
                  A$hsq <- paste(strsplit(stringr::str_split(hsq,pattern = "[:alpha:]")[[1]][1], 
                                          split = "[.]")[[1]][2],
                                 strsplit(stringr::str_split(hsq,pattern = "[:alpha:]")[[1]][1], 
                                          split = "[.]")[[1]][3],
                                 sep = ".")
                  qtl.effect.data[[j]] <- A
            }
            simulated.effects <- do.call(rbind,qtl.effect.data)
            all.architectures[[i]] <- simulated.effects 
      }
      print(paste("Merging",x))
      all.causalQTL.effects <- do.call(rbind, all.architectures) %>%
            dplyr::mutate(LD = paste(strsplit(stringr::str_split(LD,pattern = "[:alpha:]")[[1]][1], 
                                              split = "[.]")[[1]][2],
                                     strsplit(stringr::str_split(LD,pattern = "[:alpha:]")[[1]][1], 
                                              split = "[.]")[[1]][3],
                                     sep = ".")) %>%
            dplyr::mutate(MAF = paste(strsplit(stringr::str_split(maf,pattern = "[:alpha:]")[[1]][1], 
                                               split = "[.]")[[1]][2],
                                      strsplit(stringr::str_split(maf,pattern = "[:alpha:]")[[1]][1], 
                                               split = "[.]")[[1]][3],
                                      sep = "."))
}

# Concatenating Causal QTL Effects
print("Concatenating Causal QTL Effects")
causal.QTL.effects <- purrr::map(simulations, tidy.effects)
setwd("/projects/b1059/projects/Sam/GWAS_sensitivity")
save(causal.QTL.effects,
     file = "simulation.effects.RData")