#!/usr/bin/env Rscript
require(tidyverse)
require(tidymodels)
require(GenomicRanges)
# setwd("~/Documents/AndersenLab/NemaScan/Simulations/")
args <- commandArgs(trailingOnly = TRUE)
setwd(paste(args[1],"Simulations",sep = "/"))
markers <- as.numeric(args[2])

# Bring in QTL effects
effects <- list.files(pattern = ".par",recursive = T)
iterations <- purrr::map(effects, .f = function(x){
   paste(strsplit(strsplit(x,split = "/")[[1]][3],split = "_")[[1]][1:3], collapse = "_")
   # paste(strsplit(x,split = "_")[[1]][1:3], collapse = "_") # LOCAL
})
arms.centers <- readr::read_tsv("/projects/b1059/projects/Sam/NemaScan/bin/ARMS_CENTERS.tsv", 
                                col_names = c("CHR","START","STOP","TYPE"), 
                                cols())
arms.centers$CHR <- as.factor(arms.centers$CHR)
levels(arms.centers$CHR) <- c(1,2,3,4,5,6)

# Assessing Bin Performance
summarize.by.marker.bins <- function(x,y){
   n.markers.per.bin <- markers
   
   bins <- length(levels(as.factor(
      as.numeric(
         as.factor(rep(1:nrow(x)/n.markers.per.bin, each = n.markers.per.bin)[1:nrow(x)])
         )
      )))
   
   x %>%
      dplyr::mutate(SNP.bin = as.numeric(as.factor(rep(1:nrow(x)/n.markers.per.bin, each = n.markers.per.bin)[1:nrow(x)]))) %>%
      dplyr::group_by(SNP.bin) %>%
      dplyr::mutate(CHR = y) %>%
      dplyr::summarise(median(POS),
                       max(log10p),
                       TRUE %in% truth,
                       max(log10p) > -log(0.05/(bins*6)),
                       median(abs(BETA)),
                       median(AF1),
                       n(),
                       median(CHR),
                       min(POS),
                       max(POS),
                       max(POS)-min(POS)) %>%
      `colnames<-`(c("SNP.bin","median(POS)","max(log10p)","TRUE.HIT","DETECTED.HIT","median.SNP.effect","bin.median.MAF","marker.density","CHR","start","stop","bin.width.bp"))
}
fastGWA.inbred.SNP.bin.assessment <- function(x){
   
   print(paste("Assessing fastGWA (inbred) - SNP Bins", x, sep = " "))
   nQTL <- strsplit(x,split = "_")[[1]][1]
   rep <- strsplit(x,split = "_")[[1]][2]
   h2 <- strsplit(x,split = "_")[[1]][3]
   
   effects <- read.table(paste(paste(strsplit(x, split = "_")[[1]][1],
                                     "Phenotypes", sep = "/"),
                               paste(x,"sims.par",sep = "_"), sep = "/"),
                         header = T)
   # effects <- read.table("25_1_0.8_sims.par",header = T)
   
   
   # mapping <- read.table("25_1_0.8_lmm-exact_inbred.fastGWA",header = T) %>%
      mapping <- read.table(paste(paste(strsplit(x, split = "_")[[1]][1],
                                   "Mappings", sep = "/"),
                                  paste(x,"lmm-exact_inbred.fastGWA",sep = "_"), sep = "/"),
                            header = T) %>%
      dplyr::mutate(log10p = -log(P)) %>%
      dplyr::mutate(truth = as.factor(SNP %in% effects$QTL)) %>%
      dplyr::filter(CHR != 7) %>%
      droplevels() %>%
      dplyr::group_by(CHR) %>%
      tidyr::nest()
   binned.mapping <- purrr::map2(mapping$data, mapping$CHR, summarize.by.marker.bins)
   binned.mapping.df <- do.call(rbind,binned.mapping)
   
   gwa.results.gr <- GenomicRanges::GRanges(seqnames = binned.mapping.df$CHR,
                                            ranges = IRanges(start = binned.mapping.df$start,
                                                             end = binned.mapping.df$stop),
                                            POS = binned.mapping.df$`median(POS)`,
                                            BIN = binned.mapping.df$SNP.bin,
                                            log10p = binned.mapping.df$`max(log10p)`,
                                            truth = binned.mapping.df$TRUE.HIT,
                                            predicted = binned.mapping.df$DETECTED.HIT,
                                            med.effect = binned.mapping.df$median.SNP.effect,
                                            MAF = binned.mapping.df$bin.median.MAF,
                                            marker.density = binned.mapping.df$marker.density)
   
   arms.centers.gr <-  GenomicRanges::GRanges(seqnames = arms.centers$CHR,
                                              ranges = IRanges(start = arms.centers$START,
                                                               end = arms.centers$STOP),
                                              type = arms.centers$TYPE)
   
   arms.centers.labeled <- IRanges::findOverlapPairs(gwa.results.gr, arms.centers.gr) %>%
      as.data.frame() %>%
      dplyr::select(first.X.seqnames, first.X.POS, first.X.log10p,
                    first.X.truth, first.X.predicted,
                    first.X.med.effect, first.X.MAF, first.X.marker.density,
                    second.X.type, first.X.start, first.X.end) %>%
      `colnames<-`(c("CHR","POS","log10p",
                     "truth","predicted",
                     "median.effect","median.MAF","marker.density",
                     "type","bin.start","bin.end"))
   
   # For Later: Output Mappings and Bin Metrics
   # annotated.mapping <- arms.centers.labeled %>%
   #    tidyr::pivot_longer(cols = c("log10p","median.effect","median.MAF","marker.density"),
   #                        names_to = "bin.attributes",
   #                        values_to = "values")
   
   arms.centers.labeled$truth <- as.factor(arms.centers.labeled$truth)
   arms.centers.labeled$predicted <- as.factor(arms.centers.labeled$predicted)
   levels(arms.centers.labeled$predicted) <- c("FALSE","TRUE")
   
   # REQUIRED FOR ACCURATE MEASUREMENT
   arms.centers.labeled$truth <- factor(arms.centers.labeled$truth, levels = c("TRUE","FALSE"))
   arms.centers.labeled$predicted <- factor(arms.centers.labeled$predicted, levels = c("TRUE","FALSE"))
   
   class_metrics <- metric_set(accuracy,sensitivity,precision,specificity,mcc,detection_prevalence)
   a <- arms.centers.labeled %>%
      class_metrics(truth, estimate = predicted)
   npv <- arms.centers.labeled %>%
      yardstick::npv(truth, predicted)
   auc <- arms.centers.labeled %>%
      yardstick::roc_auc(truth, log10p)
   gain.capture <- arms.centers.labeled %>%
      yardstick::gain_capture(truth, log10p)
   wg <- a %>%
      dplyr::full_join(.,npv, by = c(".metric", ".estimator", ".estimate")) %>%
      dplyr::full_join(.,auc, by = c(".metric", ".estimator", ".estimate")) %>%
      dplyr::full_join(.,gain.capture, by = c(".metric", ".estimator", ".estimate")) %>%
      dplyr::mutate(nQTL = nQTL) %>%
      dplyr::mutate(rep = rep) %>%
      dplyr::mutate(h2 = h2) %>%
      dplyr::mutate(type = "genome-wide")
   
   b <- arms.centers.labeled %>%
      group_by(type) %>%
      class_metrics(truth, estimate = predicted)
   
   npv <- arms.centers.labeled %>%
      group_by(type) %>%
      yardstick::npv(truth, predicted)
   
   auc <- arms.centers.labeled %>%
      group_by(type) %>%
      yardstick::roc_auc(truth, log10p)
   
   gain.capture <- arms.centers.labeled %>%
      group_by(type) %>%
      yardstick::gain_capture(truth, log10p)
   
   ac <- b %>%
      dplyr::full_join(.,npv, by = c(".metric", ".estimator", ".estimate","type")) %>%
      dplyr::full_join(.,auc, by = c(".metric", ".estimator", ".estimate","type")) %>%
      dplyr::full_join(.,gain.capture, by = c(".metric", ".estimator", ".estimate","type")) %>%
      dplyr::mutate(nQTL = nQTL) %>%
      dplyr::mutate(rep = rep) %>%
      dplyr::mutate(h2 = h2)
   
   wg %>%
      full_join(., ac, by = c(".metric", ".estimator", ".estimate", "nQTL", "rep", "h2", "type"))
   
}
simulation.metrics.SNP.binned <- purrr::map(iterations, fastGWA.inbred.SNP.bin.assessment)
fastGWA.inbred.assessment.SNP.binned.dat <- do.call(rbind, simulation.metrics.SNP.binned)

save(fastGWA.inbred.assessment.SNP.binned.dat,
     file = paste("Sim.Perf.Bin",markers,"markers","RData", sep = "."))
