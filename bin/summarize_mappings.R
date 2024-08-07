#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readr)
library(glue)
library(purrr)

args <- commandArgs(trailingOnly = TRUE)

# args[1] = qtl_peaks.tsv
# args[2] = chr_lens

# args <- c("~/Downloads/QTL_peaks.tsv", "~/Downloads/c_elegans_chr_lengths.tsv")

goodtraits <- readr::read_tsv(args[1], col_types = c('ccncnnnnn'))
chr_lens <- readr::read_tsv(args[2]) %>%
    dplyr::select(CHROM, startPOS = start, endPOS = stop) %>%
    dplyr::mutate(trait = goodtraits$trait[1])

# don't have log10p right now, so just color all QTL red.
goodtraits%>%
    dplyr::ungroup()%>%
    ggplot()+
    aes(x=peakPOS/1E6, y=trait)+
    theme_bw() +
    scale_fill_gradient(high = "#D7263D", low = "#0072B2",
                        name = expression(-log[10](italic(p))))+
    scale_color_gradient(high = "#D7263D", low = "#0072B2",
                         name = expression(-log[10](italic(p))))+
    geom_segment(aes(x = startPOS/1e6, y = trait, xend = endPOS/1e6, yend = trait, color = log10p), size = 2, alpha = 1) +
    geom_segment(data=chr_lens,aes(x = 0, y = trait, xend = endPOS/1e6, yend = trait), size = 2.5, alpha = 0) +
    geom_point(aes(fill = log10p), colour = "black",size = 2, alpha = 1, shape = 25)+
    xlab("Genomic Position (Mb)") + ylab("") +
    scale_x_continuous(expand = c(0, 0), breaks = c(5, 10, 15, 20)) +
    theme(strip.background = element_rect(colour = "black", fill = "white",
                                          size = 0.75, linetype = "solid")) +
    theme_bw(15) + 
    theme(strip.background = element_blank(),
          panel.background = element_blank(),
          strip.text.y = element_blank(),
          panel.spacing = unit(0.5, "lines"))+
    facet_grid(. ~ CHROM, scales = "free", space = "free")+
    ggplot2::labs(x = "Genomic Position (Mb)")

ggsave(glue::glue("Summarized_mappings_{args[3]}.pdf"), height = 12, width = 12)
