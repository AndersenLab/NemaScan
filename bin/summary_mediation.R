#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readr)
library(glue)
library(purrr)


args <- commandArgs(trailingOnly = TRUE)

tr = args[3]


## mediation with multiple test correction#####

med_analysis_multi <- read.delim(args[1], header=FALSE, stringsAsFactors=FALSE)



med_multi <- med_analysis_multi %>% dplyr::filter(!V1=="gene")

med_multi_head <- med_analysis_multi %>% dplyr::filter(V1=="gene")

colnames(med_multi) <- med_multi_head[1,]



## mediation with four different effect####


mediation_analysis_indi <- read.delim(args[2], header=FALSE, stringsAsFactors=FALSE)



med_indi <- mediation_analysis_indi %>% dplyr::filter(!V1=="var")

med_indi_head<-mediation_analysis_indi %>% dplyr::filter(V1=="var")

colnames(med_indi)<-med_indi_head[1,]



### join ####

summarized_med <- med_multi %>% left_join(med_indi) %>% 
  dplyr::select(gwtrait,
                gwchr=e_chr,
                gwpeak,
                mediator=gene,
                multi_abs_est=S,
                multi_padjust=p,
                var,
                estimate,
                ci_lower,
                ci_upper,
                prob,
                mediator_gene=ext_gene,
                e_peak,
                eQTL_classification
                )



# save  data set

write.table(summarized_med, paste(tr,"_mediation_", args[4], ".tsv",sep = ""), sep = "\t", row.names = F, quote = FALSE)

# plot

medprop_all <- summarized_med %>% 
  dplyr::mutate(e_peak=as.numeric(e_peak),
                estimate=as.numeric(estimate),
                prob=as.numeric(prob)) %>% 
  dplyr::filter(var == "prop_med") %>% 
  dplyr::filter(estimate >=0 & estimate<=1) %>% 
  dplyr::mutate(gene_gwas_qtl=paste(mediator,gwchr,gwpeak, sep = "_"))


med  <- summarized_med %>% 
  dplyr::select(-c("var","estimate","ci_lower","ci_upper","prob")) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(e_peak=as.numeric(e_peak),
                multi_abs_est=as.numeric(multi_abs_est),
                multi_padjust=as.numeric(multi_padjust)) %>% 
  dplyr::mutate(gene_gwas_qtl=paste(mediator,gwchr,gwpeak, sep = "_")) %>% 
  dplyr::mutate(mediator=paste(mediator,mediator_gene,sep = "_")) %>% 
  dplyr::mutate(q99 = quantile(multi_abs_est, probs = 0.99)[[1]]) %>% 
  dplyr::mutate(prop_med=ifelse(gene_gwas_qtl %in% medprop_all$gene_gwas_qtl, "0 ≤ prop_med ≤ 1","NA"),
                q99_mediator=ifelse(multi_abs_est>q99 & prop_med=="0 ≤ prop_med ≤ 1",mediator,NA)) %>% 
  dplyr::filter(multi_padjust<0.05)


if (length(unique(med$q99_mediator)) > 1) {
  
  
med_detailed_plot <- ggplot(med,aes(x=e_peak/1e6,y=multi_abs_est,color=q99_mediator )) +
  geom_point( size=2) +
  labs(x = "Genomic position (Mb)", y = "Mediation estimate") +
  scale_alpha_continuous(range = c(1, 0.1)) +
  geom_hline(aes(yintercept = q99), color = "grey") +
  theme_bw(10) +
  theme( panel.grid = element_blank(),
         panel.background = element_blank(), 
         #  legend.position = "none",
         axis.text = element_text(size=10,color="black"),
         axis.title = element_text(size=10, color="black"),
         strip.text = element_text(size=10, color = "black")) +
  facet_grid(eQTL_classification~gwchr,scales = "free_x")

ggsave(med_detailed_plot, filename = paste(tr,"_emed_detailed_plot_", args[4], ".png",sep = ""),  units = "in",height = 10, width =  10)


}








