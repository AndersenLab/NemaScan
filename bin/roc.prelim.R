#!/usr/bin/env Rscript
require(tidyverse)
require(tidymodels)
# setwd("~/Documents/AndersenLab/NemaScan/Simulations/")
args <- commandArgs(trailingOnly = TRUE)
setwd(paste(args[1],"Simulations",sep = "/"))
h2pal <- brewer_pal(palette = "Greens")(3)


# Bring in QTL effects
effects <- list.files(pattern = ".par",recursive = T)
iterations <- purrr::map(effects, .f = function(x){
   paste(strsplit(strsplit(x,split = "/")[[1]][3],split = "_")[[1]][1:3], collapse = "_")
})

# Define Assessments
fastGWA.assessment <- function(x){
   
   print(paste("Assessing fastGWA", x, sep = " "))
   nQTL <- strsplit(x,split = "_")[[1]][1]
   rep <- strsplit(x,split = "_")[[1]][2]
   h2 <- strsplit(x,split = "_")[[1]][3]
   
   effects <- read.table(paste(paste(strsplit(x, split = "_")[[1]][1], 
                                     "Phenotypes", sep = "/"),
                               paste(x,"sims.par",sep = "_"), sep = "/"),
                         header = T)
   
   lmm.exact.fastGWA <- read.table(paste(paste(strsplit(x, split = "_")[[1]][1], 
                                      "Mappings", sep = "/"),
                                paste(x,"lmm-exact.fastGWA",sep = "_"), sep = "/"),
                          header = T)
   tests <- nrow(lmm.exact.fastGWA)
   performance <- lmm.exact.fastGWA %>%
      dplyr::mutate(log10p = -log(P)) %>% 
      dplyr::mutate(truth = as.factor(SNP %in% effects$QTL)) %>%
      dplyr::mutate(BF = as.factor(P < (0.05/tests))) %>%
      dplyr::mutate(predicted = BF)
   levels(performance$truth) <- c("nonQTL","QTL")
   levels(performance$predicted) <- c("nonQTL","QTL")
   levels(performance$BF) <- c("0","1")
   
   
   class_metrics <- metric_set(accuracy,sensitivity,precision,specificity,ppv,mcc,detection_prevalence)
   a <- performance %>%
      class_metrics(truth, estimate = predicted)
      
   npv <- performance %>%
      yardstick::npv(truth, predicted)
   
   auc <- performance %>%
      yardstick::roc_auc(truth, log10p)
   
   gain.capture <- performance %>%
      yardstick::gain_capture(truth, log10p)
   
   a %>%
      dplyr::full_join(.,npv, by = c(".metric", ".estimator", ".estimate")) %>%
      dplyr::full_join(.,auc, by = c(".metric", ".estimator", ".estimate")) %>%
      dplyr::full_join(.,gain.capture, by = c(".metric", ".estimator", ".estimate")) %>%
      dplyr::mutate(nQTL = nQTL) %>%
      dplyr::mutate(rep = rep) %>%
      dplyr::mutate(h2 = h2) %>%
      dplyr::mutate(trait = x)

}
fastGWA.inbred.assessment <- function(x){
   
   print(paste("Assessing fastGWA (inbred)", x, sep = " "))
   nQTL <- strsplit(x,split = "_")[[1]][1]
   rep <- strsplit(x,split = "_")[[1]][2]
   h2 <- strsplit(x,split = "_")[[1]][3]
   
   effects <- read.table(paste(paste(strsplit(x, split = "_")[[1]][1], 
                                     "Phenotypes", sep = "/"),
                               paste(x,"sims.par",sep = "_"), sep = "/"),
                         header = T)
   
   lmm.exact.fastGWA <- read.table(paste(paste(strsplit(x, split = "_")[[1]][1], 
                                               "Mappings", sep = "/"),
                                         paste(x,"lmm-exact_inbred.fastGWA",sep = "_"), sep = "/"),
                                   header = T)
   tests <- nrow(lmm.exact.fastGWA)
   performance <- lmm.exact.fastGWA %>%
      dplyr::mutate(log10p = -log(P)) %>% 
      dplyr::mutate(truth = as.factor(SNP %in% effects$QTL)) %>%
      dplyr::mutate(BF = as.factor(P < (0.05/tests))) %>%
      dplyr::mutate(predicted = BF)
   levels(performance$truth) <- c("nonQTL","QTL")
   levels(performance$predicted) <- c("nonQTL","QTL")
   levels(performance$BF) <- c("0","1")
   
   
   class_metrics <- metric_set(accuracy,sensitivity,precision,specificity,ppv,mcc,detection_prevalence)
   a <- performance %>%
      class_metrics(truth, estimate = predicted)
   
   npv <- performance %>%
      yardstick::npv(truth, predicted)
   
   auc <- performance %>%
      yardstick::roc_auc(truth, log10p)
   
   gain.capture <- performance %>%
      yardstick::gain_capture(truth, log10p)
   
   a %>%
      dplyr::full_join(.,npv, by = c(".metric", ".estimator", ".estimate")) %>%
      dplyr::full_join(.,auc, by = c(".metric", ".estimator", ".estimate")) %>%
      dplyr::full_join(.,gain.capture, by = c(".metric", ".estimator", ".estimate")) %>%
      dplyr::mutate(nQTL = nQTL) %>%
      dplyr::mutate(rep = rep) %>%
      dplyr::mutate(h2 = h2) %>%
      dplyr::mutate(trait = x)
   
}
EMMA.assessment <- function(x){
   
   print(paste("Assessing EMMA", x, sep = " "))
   nQTL <- strsplit(x,split = "_")[[1]][1]
   rep <- strsplit(x,split = "_")[[1]][2]
   h2 <- strsplit(x,split = "_")[[1]][3]
   
   effects <- read.table(paste(paste(strsplit(x, split = "_")[[1]][1], 
                                     "Phenotypes", sep = "/"),
                               paste(x,"sims.par",sep = "_"), sep = "/"),
                         header = T)
   
   mapping <- readr::read_tsv(paste(paste(strsplit(x, split = "_")[[1]][1], 
                                               "Mappings", sep = "/"),
                                         paste(x,"processed_mapping.tsv",sep = "_"), sep = "/"), 
                                        col_names = TRUE, cols())

   performance <- mapping %>%
      dplyr::select(marker, CHROM, POS, log10p, trait, BF, aboveBF) %>%
      dplyr::mutate(marker = stringr::str_replace(marker, "_", ":")) %>%
      dplyr::mutate(truth = as.factor(marker %in% effects$QTL)) %>%
      dplyr::mutate(predicted = as.factor(aboveBF)) %>%
      dplyr::filter(!duplicated(marker))
   levels(performance$truth) <- c("nonQTL","QTL")
   levels(performance$predicted) <- c("nonQTL","QTL")
   
   
   class_metrics <- metric_set(accuracy,
                               sensitivity,
                               precision,
                               specificity,
                               ppv,
                               mcc,
                               detection_prevalence, 
                               bal_accuracy)
   a <- performance %>%
      class_metrics(truth, estimate = predicted)
   
   npv <- performance %>%
      yardstick::npv(truth, predicted)
   
   auc <- performance %>%
      yardstick::roc_auc(truth, log10p)
   
   gain.capture <- performance %>%
      yardstick::gain_capture(truth, log10p)
   
   
   a %>%
      dplyr::full_join(.,npv, by = c(".metric", ".estimator", ".estimate")) %>%
      dplyr::full_join(.,auc, by = c(".metric", ".estimator", ".estimate")) %>%
      dplyr::full_join(.,gain.capture, by = c(".metric", ".estimator", ".estimate")) %>%
      dplyr::mutate(nQTL = nQTL) %>%
      dplyr::mutate(rep = rep) %>%
      dplyr::mutate(h2 = h2) %>%
      dplyr::mutate(trait = x)
   
}
# Conduct Assessments
simulation.metrics <- purrr::map(iterations, fastGWA.assessment)
fastGWA.assessment.dat <- do.call(rbind, simulation.metrics) %>%
   dplyr::mutate(method = "fastGWA")
simulation.metrics <- purrr::map(iterations, fastGWA.inbred.assessment)
fastGWA.inbred.assessment.dat <- do.call(rbind, simulation.metrics) %>%
   dplyr::mutate(method = "fastGWA.inbred")
simulation.metrics <- purrr::map(iterations, EMMA.assessment) 
EMMA.assessment.dat <- do.call(rbind, simulation.metrics) %>%
   dplyr::mutate(method = "EMMA")

save(fastGWA.assessment.dat,
     fastGWA.inbred.assessment.dat,
     EMMA.assessment.dat,
     file = "Simulation.Performance.RData")
