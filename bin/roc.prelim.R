#!/usr/bin/env Rscript
require(tidyverse)
require(tidymodels)
require(RColorBrewer)
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

# print("Creating Final Data Frames and Plotting")
# 
# # Plotting fastGWA Performance
# ggplot(roc, mapping = aes(x = 1-specificity, y = sensitivity, group = rep)) +
#    theme_bw() +
#    geom_line(alpha = 0.5) +
#    geom_abline(slope = 1, intercept = 0) +
#    facet_grid(h2 ~ nQTL) +    
#    xlab("1 - Specificity") + 
#    ylab("Sensitivity") + 
#    theme(legend.position = "none") +
#    ggtitle("fastGWA ROC Curve")
# ggsave("fastGWA.ROC.pdf", width = 10, height = 8)
# 
# auc <- fast.GWA.performance %>%
#    dplyr::group_by(h2, nQTL, rep) %>%
#    yardstick::roc_auc(truth, log10p)
# auc$nQTL <- factor(auc$nQTL, levels = c("5","10","20","50","100","250"))
# ggplot(auc, mapping = aes(x = nQTL, y = .estimate, fill = h2)) +
#    theme_bw() +
#    geom_boxplot() +
#    ylim(0,1) +
#    ylab("Area Under ROC Curve") +
#    geom_hline(yintercept = 0.5, linetype = 3) +
#    scale_fill_brewer(type = "seq",palette = "Greens") +
#    theme(legend.position = "none") +
#    ggtitle("fastGWA AUC Comparisons")
# ggsave("fastGWA.AUC.pdf", width = 10, height = 8)
# 
# 
# 
# # Plotting fastGWA Performance
# print("Plotting fastGWA Performance (Inbred)")
# fast.GWA.inbred.performance <- do.call(rbind, temp.2)
# roc <- fast.GWA.inbred.performance %>%
#    dplyr::group_by(h2, nQTL, rep) %>%
#    yardstick::roc_curve(truth, log10p)
# roc$nQTL <- factor(roc$nQTL, levels = c("5","10","20","50","100","250"))
# ggplot(roc, mapping = aes(x = 1-specificity, y = sensitivity, group = rep)) +
#    theme_bw() +
#    geom_line(alpha = 0.5) +
#    geom_abline(slope = 1, intercept = 0) +
#    facet_grid(h2 ~ nQTL) + 
#    theme(legend.position = "none") +
#    xlab("1 - Specificity") + 
#    ylab("Sensitivity") + 
#    ggtitle("fastGWA Inbred ROC Curve")
# ggsave("fastGWA.inbred.ROC.pdf", width = 10, height = 8)
# 
# auc <- fast.GWA.inbred.performance %>%
#    dplyr::group_by(h2, nQTL, rep) %>%
#    yardstick::roc_auc(truth, log10p)
# auc$nQTL <- factor(auc$nQTL, levels = c("5","10","20","50","100","250"))
# ggplot(auc, mapping = aes(x = nQTL, y = .estimate, fill = h2)) +
#    theme_bw() +
#    geom_boxplot() +
#    ylim(0,1) +
#    ylab("Area Under ROC Curve") +
#    geom_hline(yintercept = 0.5, linetype = 3) +
#    scale_fill_brewer(type = "seq",palette = "Greens") +
#    theme(legend.position = "none") +
#    ggtitle("fastGWA Inbred AUC Comparisons")
# ggsave("fastGWA.inbred.AUC.pdf", width = 10, height = 8)
# print("Saving Results")




