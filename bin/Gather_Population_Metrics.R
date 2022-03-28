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

# [1] = project directory
# [2] = output directory
# setwd(paste(args[1],"Simulations",sep = "/"))
today <- format(Sys.time(), '%Y%m%d')
sweeps <- data.table::fread("/projects/b1059/projects/Sam/NemaScan/bin/sweep_summary_20210121.tsv") %>%
   dplyr::rename(strain = isotype)

# strainsets <- data.table::fread(paste0("/projects/b1059/projects/Sam/NemaScan/", 
#                         args[1])) %>%
strainsets <- read.table(paste0("/projects/b1059/projects/Sam/NemaScan/", 
                         args[1])) %>%
   dplyr::rename(strain_set = V1, strains = V2) %>%
   dplyr::group_by(strain_set) %>%
   tidyr::nest()
# 
# # Populations and Strains
# populations <- list.files(pattern = ".phen",recursive = T)
# population.list <- purrr::map(populations, .f = function(x){
#       paste(strsplit(strsplit(x,split = "/")[[1]][4],split = "_")[[1]][6], collapse = "_") # QUEST
# })

# x <- strainsets$data[[1]]
# y <- strainsets$strain_set[[1]]
population.metrics <- function(x,y){
   options(dplyr.summarise.inform = FALSE)
   strains <- data.frame(stringr::str_split(x$strains,pattern = ",")[[1]]) %>%
      `colnames<-`(c("strain"))
   
   
   strain.features <- strains %>%
         dplyr::left_join(., sweeps) %>%
         tidyr::pivot_longer(cols = -strain, names_to = "feature") %>%
         dplyr::group_by(feature) %>%
         dplyr::summarise(mean.swept = mean(value),
                          median.swept = median(value),
                          sd.swept = sd(value)) %>%
         tidyr::pivot_longer(cols = -feature)
      
   quant.features <- strain.features %>%
      tidyr::unite(col = "feature", c(name, feature)) %>%
      tidyr::pivot_wider(names_from = feature, values_from = value)%>%
      dplyr::mutate(strain_set = y)
   
   
   strains %>%
      dplyr::left_join(., sweeps) %>%
      tidyr::pivot_longer(cols = -strain, names_to = "CHROM", values_to = "swept") %>%
      dplyr::mutate(swept = if_else(swept > 0.5, true = "swept","unswept")) %>%
      dplyr::group_by(CHROM, swept) %>%
      dplyr::summarise(pct.population.swept = n()/nrow(strains)) %>%
      dplyr::filter(swept == "swept") %>%
      dplyr::select(-swept) %>%
      tidyr::pivot_wider(names_from = CHROM, values_from = pct.population.swept) %>%
      # dplyr::select(-value) %>%
      dplyr::mutate(strain_set = y) %>%
      dplyr::full_join(., quant.features) %>%
      dplyr::mutate(pop.size = nrow(strains)) %>%
      dplyr::select(strain_set, !contains("strain_set"))
      
      
}
population.metrics.df <- purrr::map2(strainsets$data, strainsets$strain_set, population.metrics) %>%
   Reduce(full_join,.)
save(population.metrics.df, file = paste0(args[2],"NemaScan_Population_Metrics_",today,".RData"))

# strain.summary <- function(x){
#    options(dplyr.summarise.inform = FALSE)
#    population.phens <- list.files(pattern=paste(x,"sims.phen",sep = "_"), recursive=T)
#    strains <- purrr::map(population.phens, function(x){
#       data.table::fread(x) %>%
#          `colnames<-`(c("strain","strain2","trait.value")) %>%
#          dplyr::select(-strain2,-trait.value)}) %>%
#       Reduce(rbind,.) %>%
#       dplyr::filter(!duplicated(strain))
#    
#    strains %>%
#       dplyr::mutate(population = x,
#                     sample.size = nrow(strains))
# }
# strain.summary.df <- purrr::map(population.list, strain.summary) %>%
#    Reduce(rbind,.)
# save(strain.summary.df, file = paste("NemaScan_Strain_Summary",args[1],today,"RData", sep = "."))
