#!/usr/bin/env Rscript
require(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
setwd(paste(args[1],"Simulations",sep = "/"))
effect.files <- list.files(pattern = ".par",recursive = T)
iterations <- purrr::map(effect.files, .f = function(x){
  paste(strsplit(strsplit(x,split = "/")[[1]][3],split = "_")[[1]][1:3], collapse = "_")
})
dir.create("Plots")
man.plot <- function(x){

effects <- read.table(paste(paste(strsplit(x, split = "_")[[1]][1], 
                                  "Phenotypes", sep = "/"),
                            paste(x,"sims.par",sep = "_"), sep = "/"),
                      header = T) %>%
  tidyr::separate(col = QTL, 
                  into = c("CHR","POS"), 
                  remove = FALSE) %>% 
  dplyr::mutate(SNP = QTL)

fastGWA.inbred <- read.table(paste(paste(strsplit(x, split = "_")[[1]][1], 
                                         "Mappings", sep = "/"),
                                   paste(x,"lmm-exact_inbred.fastGWA",sep = "_"), sep = "/"),
                             header = T) %>%
  dplyr::mutate(CHR = as.factor(CHR)) %>%
  dplyr::mutate(true.QTL = as.factor(SNP %in% effects$QTL)) %>%
  dplyr::filter(CHR != "Mt") %>%
  dplyr::mutate(log10p = -log(P)) %>%
  dplyr::select(-c("A1","A2","N","P")) %>%
  dplyr::mutate(algo = "fastGWA.inbred")

fastGWA <- read.table(paste(paste(strsplit(x, split = "_")[[1]][1], 
                                  "Mappings", sep = "/"),
                            paste(x,"lmm-exact.fastGWA",sep = "_"), sep = "/"),
                      header = T) %>%
  dplyr::mutate(CHR = as.factor(CHR)) %>%
  dplyr::mutate(true.QTL = as.factor(SNP %in% effects$QTL)) %>%
  dplyr::filter(CHR != "Mt") %>%
  dplyr::mutate(log10p = -log(P)) %>%
  dplyr::select(-c("A1","A2","N","P")) %>%
  dplyr::mutate(algo = "fastGWA")

EMMA <- readr::read_tsv(paste(paste(strsplit(x, split = "_")[[1]][1], 
                              "Mappings", sep = "/"),
                        paste(x,"processed_mapping.tsv",sep = "_"), sep = "/"), 
                  col_names = TRUE, cols()) %>% 
  dplyr::select(marker, CHROM, POS, log10p, trait, BF, aboveBF) %>%
  dplyr::mutate(marker = stringr::str_replace(marker, "_", ":")) %>%
  dplyr::mutate(true.QTL = as.factor(marker %in% effects$QTL)) %>%
  dplyr::select(-c(trait)) %>%
  `colnames<-`(c("SNP","CHR","POS","log10p","BF","aboveBF","true.QTL")) %>%
  dplyr::mutate(CHR = as.factor(CHR)) %>%
  dplyr::mutate(algo = "EMMA")


levels(fastGWA.inbred$CHR) <- c("I","II","III","IV","V","X","Mt")
levels(fastGWA$CHR) <- c("I","II","III","IV","V","X","Mt")
levels(EMMA$CHR) <- c("I","II","III","IV","V","X","Mt")

combined.results <- fastGWA %>%
  dplyr::full_join(., fastGWA.inbred) %>%
  dplyr::full_join(., EMMA)

crude.cutoff <- (median(as.numeric(EMMA$BF))*nrow(fastGWA))/nrow(EMMA)

combined.results %>%
  dplyr::mutate(sig = log10p > crude.cutoff) %>%
  dplyr::arrange(true.QTL) %>%
  ggplot(aes(x = POS/1000000, 
             y = log10p, 
             shape = true.QTL,
             fill = true.QTL,
             alpha = true.QTL,
             size = true.QTL)) + 
  theme_classic() + 
  geom_hline(yintercept = crude.cutoff, linetype = 3) + 
  geom_hline(yintercept = median(as.numeric(EMMA$BF), linetype = 1)) + 
  geom_point() +
  scale_fill_manual(values = c("black","red")) + 
  scale_shape_manual(values = c(1,21)) + 
  scale_alpha_manual(values = c(0.1,1)) + 
  scale_size_manual(values = c(0.5,3)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,max(combined.results$log10p) + 5)) +
  labs(x = "Genomic position (Mb)",
       y = expression(-log[10](italic(p)))) +
  theme(legend.position = "none") + 
  facet_grid(algo~CHR, scales = "free_x", space = "free") + 
  ggsave(paste("Plots/",x,".manhattan.plot.png",sep = ""), width = 10, height = 7)

}
purrr::map(iterations, man.plot)