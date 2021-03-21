#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(ggrepel)
library(ggbeeswarm)
library(ggnewscale)
args <- commandArgs(trailingOnly = TRUE)

# Arguments
# [1]: Aggregated Mapping Result
# [2]: Isotype Sweep Metrics from CeNDR as of 20201203
# local testing
compute.LD <- function(QTL, trait){
  if(length(unique(QTL$peak_id)) <= 1){
    data.frame(unique(QTL$marker),unique(QTL$marker),NA, trait) %>%
      `colnames<-`(c("QTL1","QTL2","r2","trait"))
  } else {
    QTLcombos <- data.frame(t(combn(x = unique(QTL$peak_id), m = 2))) %>%
      `colnames<-`(c("QTL1","QTL2"))
    
    LD <- list()
    for(i in 1:nrow(QTLcombos)){
      print(i)
      peaks.of.interest <- as.numeric(QTLcombos[i,])
      haps <- QTL %>%
        dplyr::select(strain, allele, peak_id) %>%
        dplyr::filter(peak_id %in% peaks.of.interest) %>%
        dplyr::mutate(allele = if_else(allele == -1, true = "REF", false = "ALT")) %>%
        dplyr::mutate(peak_id = as.factor(peak_id)) %>%
        tidyr::pivot_wider(names_from = peak_id, values_from = allele) %>%
        `colnames<-`(c("strain","A","B")) %>%
        tidyr::unite("hap", c(A,B), sep = "_", remove = F)
      
      P <- haps %>%
        dplyr::group_by(hap) %>%
        dplyr::summarise(n()/nrow(haps)) %>%
        `colnames<-`(c("P","freq")) %>%
        tidyr::pivot_wider(names_from = P, values_from = freq)
      
      n <- QTL %>%
        dplyr::select(peak_id) %>%
        dplyr::mutate(peak_id = as.factor(peak_id)) %>%
        dplyr::group_by(peak_id) %>%
        dplyr::summarise(n()) %>%
        `colnames<-`(c("peak_id","total"))
      
      p <- QTL %>%
        dplyr::select(strain, allele, peak_id) %>%
        dplyr::filter(peak_id %in% peaks.of.interest) %>%
        dplyr::mutate(allele = if_else(allele == -1, true = "REF", false = "ALT")) %>%
        dplyr::mutate(peak_id = as.factor(peak_id)) %>%
        dplyr::group_by(allele, peak_id) %>%
        dplyr::summarise(n()) %>%
        `colnames<-`(c("p","peak_id","n")) %>%
        dplyr::left_join(.,n) %>%
        dplyr::mutate(freq = n/total) %>%
        tidyr::unite("allele", c(p,peak_id), sep = "_") %>%
        dplyr::ungroup() %>%
        dplyr::select(allele, freq) %>%
        tidyr::pivot_wider(names_from = allele, values_from = freq)
      
      pApB <- p %>%
        dplyr::select(contains("REF")) %>%
        c(.) %>%
        unlist() %>%
        prod()
      
      pApBpapb <- p %>%
        unlist() %>%
        prod()
      
      D.AB <- P$REF_REF - pApB
      r2 <- (D.AB^2)/pApBpapb
      LD[[i]] <- QTL %>%
        dplyr::select(marker, peak_id) %>%
        dplyr::filter(peak_id %in% peaks.of.interest) %>%
        dplyr::distinct()  %>%
        tidyr::pivot_wider(names_from = peak_id, values_from = marker) %>%
        `colnames<-`(c("QTL1","QTL2")) %>%
        dplyr::mutate(r2 = r2,
                      trait = trait)
    }
    Reduce(rbind, LD)
  }
}
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
    
    plot <- data %>%
      dplyr::filter(!is.na(allele)) %>%
      dplyr::mutate(SOI = strain %in% strains.of.interest,
                    SOI.2 = if_else(SOI == TRUE, true = strain, false = "")) %>%
      droplevels() %>%
      dplyr::arrange(SOI.2) %>%
      ggplot(mapping = aes(x = allele, y = value, text = SOI.2)) +
      theme_bw(base_size = 12) +
      geom_violin(aes(fill = allele), alpha = 0.5, scale = "count", draw_quantiles = c(0.25, 0.5, 0.75)) +
      scale_fill_manual(values = pal, guide = FALSE) +
      
      ggnewscale::new_scale("fill") +
      
      geom_point(aes(fill = SOI), position = pos, size = 1.5, shape = 21) +
      # geom_point(aes(colour = sweep.share*100), size = 1.1, position = pos) +
      scale_fill_manual(values = c("#9297C4","#D33E43"), guide = FALSE) +
      # scale_colour_gradient(low = "black", high = "violetred", name = "Selective Sweep (% Chromosome)") +
      
      geom_text_repel(aes(label = SOI.2),
                      colour = "black", position = pos) +
      theme(legend.position = "bottom") +
      ggtitle(paste(trait, paste("CHR",unique(data$CHROM), sep = ""),
                    paste(round(unique(data$peakPOS), digits = 2),"MB", sep = ""), sep = ": ")) +
      labs(y = "Trait Value",
           x = "Genotype")
    print(plot)
    ggsave(paste(trait,"_", paste("CHR",unique(data$CHROM), sep = ""),"_",
                 paste(round(unique(data$peakPOS), digits = 2),"MB", sep = ""), "_effect.plot.png",sep = ""), height = 5, width = 5)
  }
}

# args <- c("~/Documents/projects/NemaScan_Performance/data/5_1_0.8_0.05_gamma_complete_processed_aggregate_mapping.tsv",
#           "~/Documents/projects/NemaScan_Performance/data/sweep_summary.tsv")

combined.mappings <- data.table::fread(args[1]) %>%
  dplyr::mutate(CHROM = factor(CHROM, levels = c("I","II","III","IV","V","X","MtDNA")))
# levels(combined.mappings$CHROM) <- c("I","II","III","IV","V","X","MtDNA")
combined.mappings <- combined.mappings %>%
  dplyr::select(-marker) %>%
  tidyr::unite("marker", CHROM, POS, sep = ":", remove = F)
sweeps <- data.table::fread(args[2])

## LD PLOTS ##
nested.QTL <- combined.mappings %>%
  as.data.frame() %>%
  dplyr::filter(!is.na(peak_id), !is.na(allele)) %>%
  dplyr::select(CHROM, marker, trait, algorithm, AF1, value, strain, allele, peak_id) %>%
  dplyr::distinct() %>%
  dplyr::group_by(trait) %>%
  tidyr::nest()

trait.LD <- purrr::pmap(.l = list(nested.QTL$data, 
                                  nested.QTL$trait), 
            .f = compute.LD) %>%
  Reduce(rbind,.)

options(scipen = 999)
if(!is.null(trait.LD)){
  check <- trait.LD %>%
    dplyr::filter(!is.na(r2)) %>%
    nrow()
  if(check != 0){
    LD.plot <- trait.LD %>%
      dplyr::filter(!is.na(r2)) %>%
      ggplot(., mapping = aes(x = QTL1, y = QTL2)) + 
      theme_classic() +
      geom_tile(aes(fill = r2),colour = "black", size = 3) + 
      geom_text(aes(label = round(r2, 4))) + 
      scale_fill_gradient(low="darkgreen", high="red", limits = c(0, 1), name = expression(italic(r^2))) + 
      theme(axis.title = element_blank(),
            axis.text = element_text(colour = "black")) + 
      labs(title = paste0("Linkage Disequilibrium: ",unique(trait.LD$trait)))
    ggsave(LD.plot, filename = paste0(unique(trait.LD$trait),"_LD.plot.png"), width = 7, height = 7)
  }
}


## MANHATTAN PLOTS ##
BF.frame <- combined.mappings %>%
  dplyr::select(trait, BF) %>%
  dplyr::filter(!duplicated(trait))
for.plot <- combined.mappings %>%
  dplyr::mutate(CHROM = as.factor(CHROM)) %>%
  dplyr::filter(CHROM != "MtDNA") %>%
  dplyr::mutate(algorithm = as.factor(algorithm))
  
man.plot <- ggplot() + 
  theme_bw() + 
  geom_point(data = for.plot[which(for.plot$aboveBF == 1),], 
             mapping = aes(x = POS/1000000, 
                           y = log10p,
                           fill = algorithm), shape = 21) +
  scale_fill_manual(values = c("blue","red"), name = "Algorithm") + 
  geom_point(data = for.plot[which(for.plot$aboveBF == 0),], 
             mapping = aes(x = POS/1000000, 
                           y = log10p), 
             alpha = 0.25) +
  scale_y_continuous(expand = c(0,0), limits = c(0,max(for.plot$log10p + 1))) +
  geom_hline(data = BF.frame, aes(yintercept = BF), linetype = 2) + 
  labs(x = "Genomic position (Mb)",
       y = expression(-log[10](italic(p)))) +
  theme(legend.position = "bottom", 
        panel.grid = element_blank()) + 
  facet_grid(. ~ CHROM, scales = "free_x", space = "free") + 
  ggtitle(BF.frame$trait)
ggsave(man.plot, filename = paste0(BF.frame$trait,"_manhattan.plot.png"), width = 8, height = 4)


## SWEPTNESS & EFFECTS SUMMARY ##
proc.sweeps <- sweeps %>%
  dplyr::select(c(isotype,contains("hapshare")))
colnames(proc.sweeps) <- gsub(colnames(proc.sweeps),pattern = "_hapshare", replacement = "")
sweep.chrom.pivot <- proc.sweeps %>%
  tidyr::pivot_longer(cols = -isotype, names_to = "CHROM", values_to = "sweep.share") %>%
  dplyr::rename(strain = isotype)

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
    dplyr::left_join(.,sweep.chrom.pivot) %>%
    dplyr::group_by(trait, peak_id) %>%
    tidyr::nest()
  
  
  purrr::pmap(.l = list(nested.pxg.dat$trait,
                        nested.pxg.dat$data), 
              .f = pxg.plots)
}

  
