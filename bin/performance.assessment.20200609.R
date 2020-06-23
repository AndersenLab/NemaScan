#!/usr/bin/env Rscript
require(tidyverse)
require(tidymodels)
require(GenomicRanges)

# # LOCAL
# setwd("~/Dropbox/AndersenLab/LabFolders/Sam/projects/NemaScan/")
# markers <- 50

# QUEST
args <- commandArgs(trailingOnly = TRUE)
setwd(paste(args[1],"Simulations",sep = "/"))
markers <- as.numeric(args[2])


# Bring in QTL effects
effects <- list.files(pattern = ".par",recursive = T)
# effects <- list.files(pattern = "5_37_0.5_0.05_complete_sims.par")

iterations <- purrr::map(effects, .f = function(x){
   paste(strsplit(strsplit(x,split = "/")[[1]][3],split = "_")[[1]][1:5], collapse = "_") # QUEST
   # paste(strsplit(x,split = "_")[[1]][1:5], collapse = "_") # LOCAL
})
arms.centers <- readr::read_tsv("/projects/b1059/projects/Sam/NemaScan/bin/ARMS_CENTERS.tsv",
                                col_names = c("CHR","START","STOP","TYPE"),
                                cols())
# arms.centers <- readr::read_tsv("ARMS_CENTERS.tsv", 
#                                 col_names = c("CHR","START","STOP","TYPE"), 
#                                 cols())
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
   
   test <- x %>%
      dplyr::mutate(SNP.bin = as.numeric(as.factor(rep(1:nrow(x)/n.markers.per.bin, each = n.markers.per.bin)[1:nrow(x)]))) %>%
      dplyr::group_by(SNP.bin) %>%
      dplyr::mutate(CHR = y) %>%
      dplyr::summarise(median(POS),
                       max(log10p),
                       TRUE %in% truth,
                       max(log10p) > thresh, # thresh defined when all data read in, before chrs split
                       median(abs(BETA)),
                       median(AF1),
                       n(),
                       median(CHR),
                       min(POS),
                       max(POS),
                       max(POS)-min(POS)) %>%
      `colnames<-`(c("SNP.bin","median(POS)","max(log10p)","TRUE.HIT","DETECTED.HIT","median.SNP.effect","bin.median.MAF","marker.density","CHR","start","stop","bin.width.bp"))
}
summarize.by.marker.bins.RRBLUP <- function(x,y){
   
   n.markers.per.bin <- markers
   bins <- length(levels(as.factor(
      as.numeric(
         as.factor(rep(1:nrow(x)/n.markers.per.bin, each = n.markers.per.bin)[1:nrow(x)])
      )
   )))
   
   test <- x %>%
      dplyr::mutate(SNP.bin = as.numeric(as.factor(rep(1:nrow(x)/n.markers.per.bin, each = n.markers.per.bin)[1:nrow(x)]))) %>%
      dplyr::group_by(SNP.bin) %>%
      dplyr::mutate(CHROM = y) %>%
      dplyr::summarise(median(POS),
                       max(log10p),
                       TRUE %in% truth,
                       max(log10p) > thresh, # thresh defined when all data read in, before chrs split
                       n(),
                       median(CHROM),
                       min(POS),
                       max(POS),
                       max(POS)-min(POS)) %>% 
      `colnames<-`(c("SNP.bin","median(POS)","max(log10p)","TRUE.HIT","DETECTED.HIT","marker.density","CHR","start","stop","bin.width.bp"))
}
bin.assessment <- function(x){
   
   nQTL <- strsplit(x,split = "_")[[1]][1]
   rep <- strsplit(x,split = "_")[[1]][2]
   h2 <- strsplit(x,split = "_")[[1]][3]
   MAF <- strsplit(x,split = "_")[[1]][4]
   sample.population <- strsplit(x,split = "_")[[1]][5]
   
   print(x)
   
   effects <- read.table(paste(paste(strsplit(x, split = "_")[[1]][1],
                                     "Phenotypes", sep = "/"),
                               paste(x,"sims.par",sep = "_"), sep = "/"),
                         header = T)
   # effects <- read.table("5_37_0.5_0.05_complete_sims.par", header = T)
   
   #### GCTA: LMM-EXACT w/ INBRED GRM ####
   # mapping <- read.table("5_37_0.5_0.05_complete_lmm-exact_inbred.fastGWA",header = T)
   mapping <- read.table(paste(paste(strsplit(x, split = "_")[[1]][1],
                                   "Mappings", sep = "/"),
                                  paste(x,"lmm-exact.fastGWA",sep = "_"), sep = "/"),
                            header = T)
         
   bins <- length(levels(as.factor(as.numeric(as.factor(rep(1:nrow(mapping)/markers, 
                                                                  each = markers)[1:nrow(mapping)])))))
   # thresh <<- -log(0.05/(bins)) # BIN THRESH
   thresh <<- -log(0.05/(nrow(mapping))) # PURE BF
   mapping <- mapping %>%
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
                                            marker.density = binned.mapping.df$marker.density,
                                            bin.width = binned.mapping.df$bin.width.bp)
   arms.centers.gr <-  GenomicRanges::GRanges(seqnames = arms.centers$CHR,
                                              ranges = IRanges(start = arms.centers$START,
                                                               end = arms.centers$STOP),
                                              type = arms.centers$TYPE)
   
   arms.centers.labeled <- IRanges::findOverlapPairs(gwa.results.gr, arms.centers.gr) %>%
      as.data.frame() %>%
      dplyr::select(first.X.seqnames, first.X.POS, first.X.log10p,
                    first.X.truth, first.X.predicted,
                    first.X.med.effect, first.X.MAF, first.X.marker.density, first.X.bin.width,
                    second.X.type, first.X.start, first.X.end) %>%
      `colnames<-`(c("CHR","POS","log10p",
                     "truth","predicted",
                     "median.effect","median.MAF","marker.density","bin.width.bp",
                     "type","bin.start","bin.end"))
   
   arms.centers.labeled$truth <- as.factor(arms.centers.labeled$truth)
   arms.centers.labeled$predicted <- as.factor(arms.centers.labeled$predicted)
   levels(arms.centers.labeled$predicted) <- c("FALSE","TRUE")
   # REQUIRED FOR ACCURATE MEASUREMENT
   arms.centers.labeled$truth <- factor(arms.centers.labeled$truth, levels = c("TRUE","FALSE"))
   arms.centers.labeled$predicted <- factor(arms.centers.labeled$predicted, levels = c("TRUE","FALSE"))
   lmm.exact.inbred <- arms.centers.labeled %>% 
      dplyr::mutate(algo = as.factor("lmm.exact.inbred")) %>%
      dplyr::mutate(type = as.factor(type))
   #####
   
   #### GCTA: LMM-EXACT ####
   # mapping <- read.table("5_37_0.5_0.05_complete_lmm-exact.fastGWA",header = T)
   mapping <- read.table(paste(paste(strsplit(x, split = "_")[[1]][1],
                                   "Mappings", sep = "/"),
                                  paste(x,"lmm-exact_inbred.fastGWA",sep = "_"), sep = "/"),
                            header = T)
   
   bins <- length(levels(as.factor(as.numeric(as.factor(rep(1:nrow(mapping)/markers, 
                                                            each = markers)[1:nrow(mapping)])))))
   # thresh <<- -log(0.05/(bins)) # BIN THRESH
   thresh <<- -log(0.05/(nrow(mapping))) # PURE BF
   mapping <- mapping %>%
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
                                            marker.density = binned.mapping.df$marker.density,
                                            bin.width = binned.mapping.df$bin.width.bp)
   arms.centers.gr <-  GenomicRanges::GRanges(seqnames = arms.centers$CHR,
                                              ranges = IRanges(start = arms.centers$START,
                                                               end = arms.centers$STOP),
                                              type = arms.centers$TYPE)
   
   arms.centers.labeled <- IRanges::findOverlapPairs(gwa.results.gr, arms.centers.gr) %>%
      as.data.frame() %>%
      dplyr::select(first.X.seqnames, first.X.POS, first.X.log10p,
                    first.X.truth, first.X.predicted,
                    first.X.med.effect, first.X.MAF, first.X.marker.density, first.X.bin.width,
                    second.X.type, first.X.start, first.X.end) %>%
      `colnames<-`(c("CHR","POS","log10p",
                     "truth","predicted",
                     "median.effect","median.MAF","marker.density","bin.width.bp",
                     "type","bin.start","bin.end"))
   
   arms.centers.labeled$truth <- as.factor(arms.centers.labeled$truth)
   arms.centers.labeled$predicted <- as.factor(arms.centers.labeled$predicted)
   levels(arms.centers.labeled$predicted) <- c("FALSE","TRUE")
   # REQUIRED FOR ACCURATE MEASUREMENT
   arms.centers.labeled$truth <- factor(arms.centers.labeled$truth, levels = c("TRUE","FALSE"))
   arms.centers.labeled$predicted <- factor(arms.centers.labeled$predicted, levels = c("TRUE","FALSE"))
   lmm.exact <- arms.centers.labeled %>% 
      dplyr::mutate(algo = as.factor("lmm.exact")) %>%
      dplyr::mutate(type = as.factor(type))
   #####
   
   #### rrBLUP: EMMAx ####
   # mapping <- readr::read_tsv("5_37_0.5_0.05_complete_processed_mapping.tsv", col_names = TRUE, cols())
   mapping <- readr::read_tsv(paste(paste(strsplit(x, split = "_")[[1]][1],
                                          "Mappings", sep = "/"),
                                    paste(x,"processed_mapping.tsv",sep = "_"), sep = "/"),
                              col_names = TRUE, cols())
   
   bins <- length(levels(as.factor(as.numeric(as.factor(rep(1:nrow(mapping)/markers, 
                                                            each = markers)[1:nrow(mapping)])))))
   # thresh <<- -log(0.05/(bins)) # BIN THRESH
   thresh <<- -log(0.05/(nrow(mapping))) # PURE BF
   mapping <- mapping %>%
      dplyr::mutate(SNP = paste(CHROM,POS,sep = ":")) %>%
      dplyr::mutate(truth = as.factor(SNP %in% effects$QTL)) %>%
      dplyr::filter(CHROM != 7) %>%
      dplyr::filter(!duplicated(marker)) %>%
      dplyr::select(marker, CHROM, POS, log10p, trait, BF, aboveBF, SNP, truth) %>%
      droplevels() %>%
      dplyr::group_by(CHROM) %>%
      tidyr::nest()
   binned.mapping <- purrr::map2(mapping$data, mapping$CHROM, summarize.by.marker.bins.RRBLUP)
   binned.mapping.df <- do.call(rbind,binned.mapping)
   gwa.results.gr <- GenomicRanges::GRanges(seqnames = binned.mapping.df$CHR,
                                            ranges = IRanges(start = binned.mapping.df$start,
                                                             end = binned.mapping.df$stop),
                                            POS = binned.mapping.df$`median(POS)`,
                                            BIN = binned.mapping.df$SNP.bin,
                                            log10p = binned.mapping.df$`max(log10p)`,
                                            truth = binned.mapping.df$TRUE.HIT,
                                            predicted = binned.mapping.df$DETECTED.HIT,
                                            marker.density = binned.mapping.df$marker.density,
                                            bin.width = binned.mapping.df$bin.width.bp)
   arms.centers.gr <-  GenomicRanges::GRanges(seqnames = arms.centers$CHR,
                                              ranges = IRanges(start = arms.centers$START,
                                                               end = arms.centers$STOP),
                                              type = arms.centers$TYPE)
   
   arms.centers.labeled <- IRanges::findOverlapPairs(gwa.results.gr, arms.centers.gr) %>%
      as.data.frame() %>%
      dplyr::select(first.X.seqnames, first.X.POS, first.X.log10p,
                    first.X.truth, first.X.predicted,
                    first.X.marker.density, first.X.bin.width,
                    second.X.type, first.X.start, first.X.end) %>%
      `colnames<-`(c("CHR","POS","log10p",
                     "truth","predicted",
                     "marker.density","bin.width.bp",
                     "type","bin.start","bin.end"))
   
   arms.centers.labeled$truth <- as.factor(arms.centers.labeled$truth)
   arms.centers.labeled$predicted <- as.factor(arms.centers.labeled$predicted)
   levels(arms.centers.labeled$predicted) <- c("FALSE","TRUE")
   # REQUIRED FOR ACCURATE MEASUREMENT
   arms.centers.labeled$truth <- factor(arms.centers.labeled$truth, levels = c("TRUE","FALSE"))
   arms.centers.labeled$predicted <- factor(arms.centers.labeled$predicted, levels = c("TRUE","FALSE"))
   EMMA <- arms.centers.labeled %>% 
      dplyr::mutate(algo = as.factor("EMMA")) %>%
      dplyr::mutate(type = as.factor(type))
   #####
   

   all.assessments <- lmm.exact.inbred %>%
      dplyr::full_join(., lmm.exact, by = c("CHR", "POS", "log10p", "truth", "predicted", "median.effect", "median.MAF", "marker.density", "bin.width.bp", "type", "bin.start", "bin.end", "algo")) %>%
      dplyr::full_join(., EMMA, by = c("CHR", "POS", "log10p", "truth", "predicted", "marker.density", "bin.width.bp", "type", "bin.start", "bin.end", "algo"))
   
   
   class_metrics <- metric_set(accuracy,sensitivity,precision,specificity)
   a <- all.assessments %>%
      dplyr::group_by(algo) %>%
      class_metrics(truth, estimate = predicted)
   auc <- all.assessments %>%
      dplyr::group_by(algo) %>%
      yardstick::roc_auc(truth, log10p)
   wg <- a %>%
      dplyr::full_join(.,auc) %>%
      dplyr::mutate(nQTL = nQTL) %>%
      dplyr::mutate(rep = rep) %>%
      dplyr::mutate(h2 = h2) %>%
      dplyr::mutate(MAF = MAF) %>%
      dplyr::mutate(sample.population = sample.population) %>%
      dplyr::mutate(type = "genome-wide")
   
   
   
   b <- all.assessments %>%
      group_by(type, algo) %>%
      class_metrics(truth, estimate = predicted)
   auc <- all.assessments %>%
      group_by(type, algo) %>%
      yardstick::roc_auc(truth, log10p)
   
   ac <- b %>%
      dplyr::full_join(.,auc) %>%
      dplyr::mutate(nQTL = nQTL) %>%
      dplyr::mutate(rep = rep) %>%
      dplyr::mutate(h2 = h2) %>%
      dplyr::mutate(MAF = MAF) %>%
      dplyr::mutate(sample.population = sample.population)
   
   wg %>% 
      full_join(., ac, by = c("algo", ".metric", ".estimator", ".estimate", "type", "nQTL", "rep", "h2", "MAF","sample.population"))
   
}
simulation.metrics.SNP.binned <- purrr::map(iterations, bin.assessment)
assessment.SNP.binned.dat <- do.call(rbind, simulation.metrics.SNP.binned)

save(assessment.SNP.binned.dat,
     file = paste("performance",markers,"marker.bins","RData", sep = "."))
