#!/usr/bin/env Rscript
require(tidyverse)
require(tidymodels)
require(data.table)
require(GenomicRanges)

#############
# ARGUMENTS #
#############
# [1] = project directory
args <- commandArgs(trailingOnly = TRUE)
setwd(paste(args[1],"Simulations",sep = "/"))
today <- format(Sys.time(), '%Y%m%d')


poopy
# Simulated QTLs and Effects
effects <- list.files(pattern = ".par",recursive = T)
iterations <- purrr::map(effects, .f = function(x){
   paste(strsplit(strsplit(x,split = "/")[[1]][4],split = "_")[[1]][1:6], collapse = "_") # QUEST
})

# # Arms and Centers
# arms.centers <- readr::read_tsv("/projects/b1059/projects/Sam/NemaScan/bin/ARMS_CENTERS.tsv",
#                                 col_names = c("CHR","START","STOP","TYPE"),
#                                 cols())
# arms.centers$CHR <- as.factor(arms.centers$CHR)
# levels(arms.centers$CHR) <- c(1,2,3,4,5,6)


# Assessing Mapping Performance
simulation.metrics <- function(x){
   
   nQTL <- strsplit(x,split = "_")[[1]][1]
   rep <- strsplit(x,split = "_")[[1]][2]
   h2 <- strsplit(x,split = "_")[[1]][3]
   MAF <- strsplit(x,split = "_")[[1]][4]
   effect.range <- strsplit(x,split = "_")[[1]][5]
   sample.population <- strsplit(x,split = "_")[[1]][6]
   
   print(x)
   
   effects <- data.table::fread(paste(paste(effect.range,
                                     nQTL,
                                     "Phenotypes", sep = "/"),
                                     paste(x,"sims.par",sep = "_"), sep = "/"),
                                header = T)
   effects <- effects %>% 
      tidyr::separate(QTL, c("CHROM","POS"), sep = ":", remove = F)
   
   safe.lmm.exact.inbred <- purrr::safely(.f = function(x){data.table::fread(paste(paste(effect.range,nQTL,"Mappings", sep = "/"),
                                                                                   paste(x,"processed_LMM_EXACT_INBRED_mapping.tsv",sep = "_"), sep = "/"),header = T)}, 
                                          otherwise = "No successful mapping matching simulation parameters :( ")
   
   safe.lmm.exact <- purrr::safely(.f = function(x){data.table::fread(paste(paste(effect.range,nQTL,"Mappings", sep = "/"),
                                                                                   paste(x,"processed_LMM_EXACT_mapping.tsv",sep = "_"), sep = "/"),header = T)}, 
                                          otherwise = "No successful mapping matching simulation parameters :( ")
   
   safe.EMMA <- purrr::safely(.f = function(x){data.table::fread(paste(paste(effect.range,nQTL,"Mappings", sep = "/"),
                                                                       paste(x,"processed_mapping.tsv",sep = "_"), sep = "/"),header = T) %>%
         dplyr::mutate(marker = gsub(marker, pattern = "_", replacement = ":"))}, 
         otherwise = "No successful mapping matching simulation parameters :( ")
   
   #### GCTA: LMM-EXACT w/ INBRED GRM ####
   mapping.lmm.exact.inbred <- safe.lmm.exact.inbred(x)
   
   if(is.character(mapping.lmm.exact.inbred$result)){
      lmm.exact.inbred <- c("No mapping for simulation parameters")
   } else {
   peak.info <- mapping.lmm.exact.inbred$result %>%
      dplyr::mutate(causal.variant = as.factor(marker %in% effects$QTL)) %>%
      dplyr::filter(!is.na(peak_id)) %>%
      dplyr::select(CHROM, marker, POS, AF1, BETA, log10p, startPOS, peakPOS, endPOS, peak_id, interval_size, causal.variant, var.exp) %>%
      dplyr::filter(!duplicated(.)) %>%
      dplyr::mutate(detected.peak = marker)
   
   peaks <- GenomicRanges::GRanges(seqnames = peak.info$CHROM,
                                   ranges = IRanges::IRanges(start = peak.info$startPOS, 
                                                             end = peak.info$endPOS),
                                   peakPOS = peak.info$peakPOS,
                                   detected.peak = peak.info$detected.peak)
   real.effects <- GenomicRanges::GRanges(seqnames = effects$CHROM,
                                          ranges = IRanges::IRanges(start = as.numeric(effects$POS), 
                                                                    end = as.numeric(effects$POS)),
                                          QTL = effects$QTL)
   overlap <- IRanges::findOverlapPairs(real.effects, peaks) %>%
      as.data.frame() %>%
      dplyr::select(first.QTL, second.X.start, second.peakPOS, second.X.end, second.detected.peak) %>%
      `colnames<-`(c("QTL","startPOS","peakPOS","endPOS","detected.peak")) %>%
      dplyr::right_join(., peak.info) %>%
      dplyr::mutate(QTL = if_else(is.na(QTL), true = detected.peak, false = QTL))
   
   all.QTL <- data.frame(c(effects$QTL, overlap$QTL)) %>%
      `colnames<-`(c("QTL")) %>%
      dplyr::filter(!duplicated(QTL)) %>%
      dplyr::mutate(QTL = as.character(QTL),
                    Simulated = (QTL %in% effects$QTL),
                    Detected = (QTL %in% overlap$QTL)) %>%
      dplyr::full_join(.,effects, by = "QTL") %>%
      dplyr::full_join(.,overlap, by = "QTL") %>%
      dplyr::select(-CHROM.y, -marker, -POS.y, -AF1) %>%
      dplyr::mutate(sim = x,
                    algorithm = "LMM-EXACT-INBRED")

   all.QTL$Simulated <- factor(all.QTL$Simulated, levels = c("TRUE","FALSE"))
   all.QTL$Detected <- factor(all.QTL$Detected, levels = c("TRUE","FALSE"))
   lmm.exact.inbred <- all.QTL
   }
   #####
   
   #### GCTA: LMM-EXACT ####
   mapping.lmm.exact <- safe.lmm.exact(x)
   if(is.character(mapping.lmm.exact$result)){
      lmm.exact <- c("No mapping for simulation parameters")
   } else {
   peak.info <- mapping.lmm.exact$result %>%
      dplyr::mutate(causal.variant = as.factor(marker %in% effects$QTL)) %>%
      dplyr::filter(!is.na(peak_id)) %>%
      dplyr::select(CHROM, marker, POS, AF1, BETA, log10p, startPOS, peakPOS, endPOS, peak_id, interval_size, causal.variant, var.exp) %>%
      dplyr::filter(!duplicated(.)) %>%
      dplyr::mutate(detected.peak = marker)
   
   peaks <- GenomicRanges::GRanges(seqnames = peak.info$CHROM,
                                   ranges = IRanges::IRanges(start = peak.info$startPOS, 
                                                             end = peak.info$endPOS),
                                   peakPOS = peak.info$peakPOS,
                                   detected.peak = peak.info$detected.peak)
   real.effects <- GenomicRanges::GRanges(seqnames = effects$CHROM,
                                          ranges = IRanges::IRanges(start = as.numeric(effects$POS), 
                                                                    end = as.numeric(effects$POS)),
                                          QTL = effects$QTL)
   
   overlap <- IRanges::findOverlapPairs(real.effects, peaks) %>%
      as.data.frame() %>%
      dplyr::select(first.QTL, second.X.start, second.peakPOS, second.X.end, second.detected.peak) %>%
      `colnames<-`(c("QTL","startPOS","peakPOS","endPOS","detected.peak")) %>%
      dplyr::right_join(., peak.info) %>%
      dplyr::mutate(QTL = if_else(is.na(QTL), true = detected.peak, false = QTL))
   
   all.QTL <- data.frame(c(effects$QTL, overlap$QTL)) %>%
      `colnames<-`(c("QTL")) %>%
      dplyr::filter(!duplicated(QTL)) %>%
      dplyr::mutate(QTL = as.character(QTL),
                    Simulated = (QTL %in% effects$QTL),
                    Detected = (QTL %in% overlap$QTL)) %>%
      dplyr::full_join(.,effects, by = "QTL") %>%
      dplyr::full_join(.,overlap, by = "QTL") %>%
      dplyr::select(-CHROM.y, -marker, -POS.y, -AF1) %>%
      dplyr::mutate(sim = x,
                 algorithm = "LMM-EXACT")
   
   all.QTL$Simulated <- factor(all.QTL$Simulated, levels = c("TRUE","FALSE"))
   all.QTL$Detected <- factor(all.QTL$Detected, levels = c("TRUE","FALSE"))
   lmm.exact <- all.QTL
   }
   #####
   
   #### rrBLUP: EMMAx ####
   mapping.EMMA <- safe.EMMA(x)
   if(is.character(mapping.EMMA$result)){
      EMMA <- c("No mapping for simulation parameters")
   } else {
   peak.info <- mapping.EMMA$result %>%
      dplyr::mutate(causal.variant = as.factor(marker %in% effects$QTL)) %>%
      dplyr::filter(!is.na(peak_id)) %>%
      dplyr::select(CHROM, marker, POS, log10p, startPOS, peakPOS, endPOS, peak_id, interval_size, causal.variant, var.exp) %>%
      dplyr::filter(!duplicated(.)) %>%
      dplyr::mutate(detected.peak = marker)
   
   peaks <- GenomicRanges::GRanges(seqnames = peak.info$CHROM,
                                   ranges = IRanges::IRanges(start = peak.info$startPOS, 
                                                             end = peak.info$endPOS),
                                   peakPOS = peak.info$peakPOS,
                                   detected.peak = peak.info$detected.peak)
   real.effects <- GenomicRanges::GRanges(seqnames = effects$CHROM,
                                          ranges = IRanges::IRanges(start = as.numeric(effects$POS), 
                                                                    end = as.numeric(effects$POS)),
                                          QTL = effects$QTL)
   overlap <- IRanges::findOverlapPairs(real.effects, peaks) %>%
      as.data.frame() %>%
      dplyr::select(first.QTL, second.X.start, second.peakPOS, second.X.end, second.detected.peak) %>%
      `colnames<-`(c("QTL","startPOS","peakPOS","endPOS","detected.peak")) %>%
      dplyr::right_join(., peak.info) %>%
      dplyr::mutate(QTL = if_else(is.na(QTL), true = detected.peak, false = QTL))
   
   all.QTL <- data.frame(c(effects$QTL, overlap$QTL)) %>%
      `colnames<-`(c("QTL")) %>%
      dplyr::filter(!duplicated(QTL)) %>%
      dplyr::mutate(QTL = as.character(QTL),
                    Simulated = (QTL %in% effects$QTL),
                    Detected = (QTL %in% overlap$QTL)) %>%
      dplyr::full_join(.,effects, by = "QTL") %>%
      dplyr::full_join(.,overlap, by = "QTL") %>%
      dplyr::select(-CHROM.y, -marker, -POS.y) %>%
      dplyr::mutate(sim = x,
                 algorithm = "EMMA")
   
   all.QTL$Simulated <- factor(all.QTL$Simulated, levels = c("TRUE","FALSE"))
   all.QTL$Detected <- factor(all.QTL$Detected, levels = c("TRUE","FALSE"))
   EMMA <- all.QTL
   }
   #####
   
   print(paste("Joining All:",x))
   
   if(is.character(lmm.exact.inbred)){
      if(is.character(lmm.exact)){
         if(is.character(EMMA)){print("No Successful Mappings")
         } else {all.assessments <- EMMA}
         } 
      else {
         if(is.character(EMMA)){all.assessments <- lmm.exact
         } else {all.assessments <- lmm.exact %>% dplyr::full_join(., EMMA)}
         }
      } else {
         if(is.character(lmm.exact)){
            if(is.character(EMMA)){all.assessments <- lmm.exact.inbred
            } else {all.assessments <- lmm.exact.inbred %>% dplyr::full_join(., EMMA)}
         } else { 
               if(is.character(EMMA)){all.assessments <- lmm.exact.inbred %>% dplyr::full_join(.,lmm.exact)
               } else {all.assessments <- lmm.exact.inbred %>% dplyr::full_join(.,lmm.exact) %>% dplyr::full_join(., EMMA)}
         }
         }
}
simulation.metrics.list <- purrr::map(iterations, simulation.metrics)
simuation.metrics.df <- do.call(rbind, simulation.metrics.list)
save(simuation.metrics.df, 
     file = paste("NemaScan_Performance",args[1],today,"RData", sep = "."))
