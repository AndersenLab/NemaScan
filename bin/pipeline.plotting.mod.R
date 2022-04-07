#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readr)
library(glue)
library(purrr)
library(data.table)
library(ggrepel)
library(ggbeeswarm)
library(ggnewscale)
args <- commandArgs(trailingOnly = TRUE)

# Arguments
# [1]: Aggregated Mapping Result
# [2]: Isotype Sweep Metrics from CeNDR as of 20201203

pxg.plots <- function(trait, data){
  if(nrow(data) == 1){
    print("Mitochondria!")
  } else {
    
    data$allele <- as.factor(dplyr::recode(data$allele, "-1" = "REF", "1" = "ALT"))
    data$allele <- factor(data$allele, levels = c("REF","ALT"))
    
    CB.N2.code <- data %>%
      dplyr::select(strain, allele) %>%
      dplyr::distinct() %>%
      dplyr::filter(strain %in% c("CB4856","N2","PD1074")) %>%
      droplevels()
    
    if(length(levels(CB.N2.code$allele)) < 2){
      pal <- c("#726E75","#720E07")
    } else {
      pal <- c("#FFA500","#0000ff")
    }
    
    strains.of.interest <- c("PD1074", "N2", "CB4856", "RC301", "MY16", 
                             "ECA396", "ECA36", "XZ1516", "ECA248", "AB1", 
                             "CB4507", "CB4858", "CB4855", "CB4852", "MY1", 
                             "JU319", "JU345", "JU400", "PB306", "PX174", "PX179")
    
    pos <- position_beeswarm()
    
    peakPOS <- paste(unique(data$CHROM), unique(data$peakPOS), sep = ":")
    
    plot <- data %>%
      dplyr::filter(!is.na(allele)) %>%
      dplyr::mutate(SOI = strain %in% strains.of.interest,
                    SOI.3 = dplyr::case_when(strain %in% c("N2", "PD1074") ~ "N2",
                                                        strain == "CB4856" ~ "CB",
                                                        strain %in% strains.of.interest ~ "special",
                                                        TRUE ~ "other"),
                    SOI.2 = if_else(SOI == TRUE, true = strain, false = "")) %>%
      droplevels() %>%
      dplyr::arrange(SOI.2) %>%
      ggplot(mapping = aes(x = allele, y = value, text = SOI.2)) +
      theme_bw(base_size = 12) +
      geom_violin(aes(fill = allele), alpha = 0.5, scale = "count", draw_quantiles = c(0.25, 0.5, 0.75)) +
      ggplot2::scale_fill_manual(values = c("REF" = "#A79F92", "ALT" = "mediumpurple4"), guide = FALSE) +      
      ggnewscale::new_scale("fill") +
      ggplot2::geom_point(aes(fill = SOI.3, size = SOI), position = ggbeeswarm::position_beeswarm(), shape = 21) +
      ggplot2::scale_fill_manual(values = c("N2" = "orange", "CB" = "blue", "special" ="red", "other" ="grey50"), guide = FALSE) +
      ggplot2::scale_size_manual(values = c(1.5,2.5), guide = FALSE) +
      geom_text_repel(aes(label = SOI.2),
                      colour = "black", position = pos, max.overlaps = Inf) +
      theme(legend.position = "bottom") +
      ggtitle(peakPOS) +
      labs(y = "Trait Value",
           x = "Genotype")
    ggsave(paste(trait,"_", paste("CHR",unique(data$CHROM), sep = ""),"_",
                 paste(round(unique(data$peakPOS), digits = 2),"MB", sep = ""), "_effect_", args[3], ".plot.png",sep = ""), height = 5, width = 5)
  }
}

combined.mappings <- data.table::fread(args[1]) %>%
  dplyr::mutate(CHROM = factor(CHROM, levels = c("I","II","III","IV","V","X","MtDNA")))
# levels(combined.mappings$CHROM) <- c("I","II","III","IV","V","X","MtDNA")
combined.mappings <- combined.mappings %>%
  dplyr::select(-marker) %>%
  tidyr::unite("marker", CHROM, POS, sep = ":", remove = F)
tests <- args[2]

# do we have mito mapping?
mito_check <- combined.mappings %>%
    na.omit()

## MANHATTAN PLOTS ##
for.plot <- combined.mappings %>%
  dplyr::mutate(CHROM = as.factor(CHROM)) %>%
    {
        if(!("MtDNA" %in% mito_check$CHROM)) dplyr::filter(., CHROM != "MtDNA") else .
    }

BF <- combined.mappings %>% 
    dplyr::group_by(trait) %>% 
    dplyr::filter(log10p != 0) %>% 
    dplyr::distinct(marker, log10p) %>%
    dplyr::mutate(BF = -log10(0.05/sum(log10p > 0, na.rm = T))) %>%
    dplyr::ungroup() %>%
    dplyr::select(BF) %>%
    unique(.) %>%
  as.numeric()

ntests <- data.table::fread(tests) %>% 
  as.numeric()
EIGEN <- -log10(0.05/ntests)
BF.frame <- combined.mappings %>%
  dplyr::select(trait) %>%
  dplyr::filter(!duplicated(trait)) %>%
  dplyr::mutate(BF = BF, EIGEN  = EIGEN, user = unique(combined.mappings$BF))

# if user selected a different threshold, use that, otherwise plot BF and EIGEN
if(BF.frame$user %in% c(BF.frame$BF, BF.frame$EIGEN)) {
  for.plot.ann <- for.plot %>%
    dplyr::mutate(sig = case_when(log10p > BF.frame$BF ~ "BF",
                                  log10p > BF.frame$EIGEN ~ "EIGEN",
                                  TRUE ~ "NONSIG"))
  
  sig.colors <- c("red","#EE4266", "black")
  names(sig.colors) <- c("BF","EIGEN", "NONSIG")
} else {
  for.plot.ann <- for.plot %>%
    dplyr::mutate(sig = case_when(log10p > BF.frame$user ~ "user",
                                  TRUE ~ "NONSIG"))
  
  sig.colors <- c("red", "black")
  names(sig.colors) <- c("user", "NONSIG")
}

test <- BF.frame %>%
    tidyr::pivot_longer(BF:user) %>%
    dplyr::distinct() %>%
    dplyr::filter(name %in% names(sig.colors))

# are we plotting mito or no?
if("MtDNA" %in% unique(for.plot.ann$CHROM)) {
    facet_scales <- "fixed"
} else {
    facet_scales <- "free"
}

man.plot <- ggplot() + 
  theme_bw() + 
  geom_point(data = for.plot.ann, 
             mapping = aes(x = POS/1000000, 
                           y = log10p,
                           colour = sig,
                           alpha = sig)) +
  scale_alpha_manual(values = c("BF" = 1, "EIGEN" = 1, "user" = 1, "NONSIG" = 0.25)) +
  scale_colour_manual(values = sig.colors) + 
  scale_x_continuous(expand = c(0, 0), breaks = c(5, 10, 15, 20)) +
  geom_hline(data = test, aes(yintercept = value, linetype = name)) + 
  scale_linetype_manual(values = c("BF" = 1, "EIGEN" = 3, "user" = 2)) +
  labs(x = "Genomic position (Mb)",
       y = expression(-log[10](italic(p)))) +
  theme(legend.position = "none", 
        panel.grid = element_blank()) + 
  facet_grid(. ~ CHROM, scales = "free_x", space = facet_scales) + 
  ggtitle(BF.frame$trait)

ggsave(man.plot, filename = paste0(BF.frame$trait,"_manhattan_", args[3], ".plot.png"), width = 8, height = 4)


## SWEPTNESS & EFFECTS SUMMARY ##
QTLcheck <- combined.mappings %>%
  dplyr::filter(!is.na(peak_id)) %>%
  nrow()
if(QTLcheck > 0){
  nested.pxg.dat <- combined.mappings %>%
    dplyr::filter(!is.na(peak_id)) %>%
    dplyr::select(CHROM, marker, trait, startPOS, peakPOS, endPOS, AF1, value, strain, allele, peak_id) %>%
    dplyr::distinct() %>%
    dplyr::mutate(startPOS = startPOS/1000000,
                  peakPOS = peakPOS/1000000,
                  endPOS = endPOS/1000000) %>%
    dplyr::group_by(trait, peak_id) %>%
    tidyr::nest()
  
  
  purrr::pmap(.l = list(nested.pxg.dat$trait,
                        nested.pxg.dat$data), 
              .f = pxg.plots)
}

  
