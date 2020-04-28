#!/usr/bin/env Rscript
library(tidyverse)
library(reshape2)

simulations <- list.files(pattern="HSQ")
# Concatenation Functions
tidy.phenotypes <- function(x){
  print(paste("Beginning",x,sep = " "))
  setwd(paste("/projects/b1059/projects/Sam/GWAS_sensitivity/",x, sep = ""))
  hsq <- stringr::str_split(x, pattern = "[:alpha:]")[[1]][4]
  maf <- stringr::str_split(x, pattern = "[:alpha:]")[[1]][7]
  LD <- stringr::str_split(x, pattern = "[:alpha:]")[[1]][9]
  
  architectures <- list.files(pattern="QTL")
  all.architectures <- list()
  for(i in 1:length(architectures)){
    phenotype.files <- list.files(path = paste(architectures[[i]],"simulated.phenotypes",sep = "/"))
    phenotype.data <- list()
    print(paste("Beginning", architectures[[i]]))
    for(j in 1:length(phenotype.files)){
      A <- read.table(paste(getwd(), 
                            architectures[[i]], 
                            "simulated.phenotypes", 
                            paste("sim.phenos.QTLs.rep.",j,sep = ""),
                            sep = "/")) %>%
        dplyr::select(-V2)
      colnames(A) <- c("strain", paste("sim",j,
                                       strsplit(stringr::str_split(architectures[[i]],pattern = "[:alpha:]")[[1]][1], 
                                                split = "[.]")[[1]],
                                       paste(strsplit(stringr::str_split(hsq,pattern = "[:alpha:]")[[1]][1], 
                                                split = "[.]")[[1]][2],
                                             strsplit(stringr::str_split(hsq,pattern = "[:alpha:]")[[1]][1], 
                                                split = "[.]")[[1]][3],
                                             sep = "."), 
                                       sep = "_"))
      phenotype.data[[j]] <- A
    }
    simulated.phenotypes <- Reduce(merge,phenotype.data)
    all.architectures[[i]] <- simulated.phenotypes 
  }
  print(paste("Merging",x))
  all.simulated.phenotypes <- Reduce(merge, all.architectures) %>%
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

# Concatenating Simulated Phenotypes
print("Concatenating Simulated Phenotypes")
phenotype.simulations <- purrr::map(simulations, tidy.phenotypes)
setwd("/projects/b1059/projects/Sam/GWAS_sensitivity")
save(phenotype.simulations,
     file = "simulated.phenotypes.RData")