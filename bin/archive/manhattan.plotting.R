#!/usr/bin/env Rscript
library(tidyverse)
# argument information - mapping
# 1 - nQTL
# 2 - Rep
# 3 - h2
# 4 - MAF
# 5 - effect range (chr)
# 6 - strain set
# 7 - path to mappings
# EXAMPLE
# Rscript bin/manhattan.plotting.R 5 1 0.1 0.05 gamma complete CeNDR2020_PowerPrecision_GAMMA
args <- commandArgs(trailingOnly = TRUE)
setwd(args[7])
file.prefix <- paste(args[1],args[2],args[3],args[4],args[5],args[6],sep = "_")
mapping.dir <- paste("Simulations", args[5], args[1], "Mappings", sep = "/")
real.effects <- data.table::fread(paste("Simulations", args[5], args[1], "Phenotypes", 
                                        paste(file.prefix,"sims.par",sep = "_"), sep = "/"))

gcta.inbred <- data.table::fread(paste(mapping.dir,"/",file.prefix,
                                       "_processed_LMM_EXACT_INBRED_mapping.tsv",sep = "")) %>%
  dplyr::mutate(algorithm = "LMM-EXACT-INBRED",
                strain = as.character(strain))
# gcta <- data.table::fread(paste(mapping.dir,"/",file.prefix,
#                                 "_processed_LMM_EXACT_mapping.tsv",sep = "")) %>%
#   dplyr::mutate(algorithm = "LMM-EXACT",
#                 strain = as.character(strain))
# gcta.loco.inbred <- data.table::fread(paste(mapping.dir,"/",file.prefix,
#                                 "_processed_LMM_EXACT_INBRED_LOCO_mapping.tsv",sep = "")) %>%
#   dplyr::mutate(algorithm = "LMM-EXACT-INBRED-LOCO",
#                 strain = as.character(strain))

gcta.loco <- data.table::fread(paste(mapping.dir,"/",file.prefix,
                                     "_processed_LMM_EXACT_LOCO_mapping.tsv",sep = "")) %>%
  dplyr::mutate(algorithm = "LMM-EXACT-LOCO",
                strain = as.character(strain))

# emmax <- data.table::fread(paste(mapping.dir,"/",file.prefix,
#                                 "_processed_mapping.tsv",sep = "")) %>%
#   dplyr::mutate(algorithm = "EMMAx",
#                 strain = as.character(strain),
#                 marker = gsub(marker, pattern = "_", replacement = ":"))

combined.fastGWA.results <- gcta.inbred %>%
  dplyr::full_join(., gcta.loco) %>%
  dplyr::mutate(CHROM = as.factor(CHROM),
                Simulated = marker %in% real.effects$QTL)

BF <- unique(combined.fastGWA.results$BF)[1]
if(max(combined.fastGWA.results$log10p) > BF){
  lim <- max(combined.fastGWA.results$log10p) + 2
} else {
  lim <- BF + 2
}
combined.fastGWA.results %>%
  dplyr::arrange(Simulated) %>%
  ggplot(.,aes(x = POS/1000000, 
             y = log10p,
             fill = as.factor(Simulated),
             colour = as.factor(aboveBF))) + 
  theme_classic() + 
  geom_rect(aes(xmin = startPOS/1000000, 
                xmax = endPOS/1000000, 
                ymin = 0, 
                ymax = Inf, 
                fill = "blue"),
            color = "blue",
            fill = "cyan",
            linetype = 2, 
            alpha=.3)+ 
  geom_point(shape = 21) +
  scale_colour_manual(values = c("black","red")) + 
  scale_fill_manual(values = c("black","yellow")) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,lim)) +
  geom_hline(yintercept = BF, linetype = 2) +
  labs(x = "Genomic position (Mb)",
       y = expression(-log[10](italic(p)))) +
  theme(legend.position = "none") + 
  facet_grid(algorithm~CHROM, scales = "free_x", space = "free") + 
  ggtitle(paste("Manhattan Plot:",file.prefix)) + 
  ggsave(paste(file.prefix,"manhattan.plot.png", sep = "."), width = 12, height = 6)
