#!/usr/bin/env Rscript
require(tidyverse)
require(data.table)

#############
# ARGUMENTS #
#############
# [1] = project directory
args <- commandArgs(trailingOnly = TRUE)
today <- format(Sys.time(), '%Y%m%d')

# Independent Tests
setwd(paste(args[1],"Genotype_Matrix",sep = "/"))
tests <- list.files(pattern = "total_independent_tests.txt",recursive = T)
gather.independent.tests <- function(x){
   data.table::fread(x) %>%
      data.frame() %>%
      dplyr::rename(independent_tests = V1) %>%
      dplyr::mutate(strain_set = strsplit(x, "_")[[1]][1])
}
independent.tests <- purrr::map(tests, gather.independent.tests) %>%
   Reduce(rbind,.)
write.csv(independent.tests, file = paste(args[1],today,"independent.tests.csv",sep = "_"), row.names = F)