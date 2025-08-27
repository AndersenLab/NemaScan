#!/usr/bin/env Rscript
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(stringr)
library(readr)
library(glue)
library(purrr)

# input arguments
# 1 = LD file
# 2 = phenotype file
# 3 = gene file
# 4 = annotation file
# 5 = algorithm
# 6 = species

args <- commandArgs(trailingOnly = TRUE)


pr_trait_ld <- data.table::fread(args[1]) %>%
    dplyr::mutate(CHR = case_when(CHR == 1 ~ "I",
                                  CHR == 2 ~ "II",
                                  CHR == 3 ~ "III",
                                  CHR == 4 ~ "IV",
                                  CHR == 5 ~ "V",
                                  CHR == 6 ~ "X",
                                  CHR == 7 ~ "MtDNA")) %>%
    dplyr::mutate(marker = paste(CHR,POS,sep = "_"),
                  log10p = -log(P))

output_sq <- sub("(.*)(\\.)(.*)(\\.)(.*)(\\.)(.*)(\\.)(prLD_df_)(.*.tsv)","\\5",args[1])
output_eq <- sub("(.*)(\\.)(.*)(\\.)(.*)(\\.)(.*)(\\.)(prLD_df_)(.*.tsv)","\\7",args[1])

phenotypes <- readr::read_tsv(args[2])

gene_ref_flat <- readr::read_tsv(args[3])

analysis_trait <- colnames(phenotypes)[2]

colnames(phenotypes) <- c("strain", "Phenotype_Value")

query_regions <- pr_trait_ld %>%
    dplyr::mutate(start_pos = min(POS),
                  end_pos = max(POS)) %>%
    dplyr::select(CHROM, start_pos, end_pos)%>%
    dplyr::mutate(CHROM = case_when(CHROM == 1 ~ "I",
                                    CHROM == 2 ~ "II",
                                    CHROM == 3 ~ "III",
                                    CHROM == 4 ~ "IV",
                                    CHROM == 5 ~ "V",
                                    CHROM == 6 ~ "X",
                                    CHROM == 7 ~ "MtDNA")) %>%
    dplyr::distinct()

query_regions

# query_regions == CHROM, start_pos, end_pos

# update 20210330 KSE: use impute vcf for genotypes and hard vcf annotation file for annotations
# this could be bcsq or snpeff
annotations <- data.table::fread(args[4])
names(annotations) = toupper(names(annotations))

# annotation type -snpeff or bcsq
if("CONSEQUENCE" %in% names(annotations)) {
    if("VARIANT_IMPACT" %in% names(annotations)) {
        ann_type <- "bcsq_old"
    } else {
        ann_type <- "bcsq"
    }
} else {
    ann_type <- "snpeff"
}

if( "STRAINS" %in% names(annotations)) {
    annotations <- annotations %>%
        dplyr::mutate(STRAIN = STRAINS) %>%
        dplyr::select(., -STRAINS)
}
if( "WORMBASE_ID" %in% names(annotations)) {
    annotations <- annotations %>%
        dplyr::mutate(WBGENE = WORMBASE_ID) %>%
        dplyr::select(., -WORMBASE_ID)
}
if( "GENE" %in% names(annotations)) {
    annotations <- annotations %>%
        dplyr::mutate(GENE_NAME = GENE) %>%
        dplyr::select(., -GENE)
}

species = args[6]

# do this for each QTL separately, then combine
annotation_out <- list()
for(r in 1:nrow(query_regions)){
    cq <- query_regions$CHROM[r]
    sq <- query_regions$start_pos[r]
    eq <- query_regions$end_pos[r]




# a function to check how many strains in the mapping set are in the column "Strains" of the annotation file

    alt_strain_check=function(df){

      df_split <- df %>% {
            if(ann_type == "bcsq") {
                tidyr::separate_rows( STRAIN, sep = " ")
             } else {
                tidyr::separate_rows( STRAIN, sep = ",")
             }} %>%
        tidyr::separate_rows( strains_used, sep = ",") %>%
        dplyr::filter(STRAIN==strains_used)

      return(nrow(df_split))

    }

# ALT strains in the mapping set at each position
    strain_list <- pr_trait_ld %>%
      dplyr::filter(CHR == cq ,
                    POS >= sq ,
                    POS <= eq ,
                    allele=="ALT") %>%
      dplyr::select(CHROM=CHR,POS,strains_used=strains) %>%
      dplyr::distinct()

# strain_list == CHROM, POS, strains_used

    annotations_region <- annotations %>%
      dplyr::filter(CHROM == cq,
                    POS >= sq,
                    POS <= eq)


    write_tsv(annotations_region,
            file = glue::glue("{query_regions$CHROM[r]}_{query_regions$start_pos[r]}_{query_regions$end_pos[r]}.test.tsv"))

# count and filter
    annotations_region <- annotations %>%
      dplyr::filter(CHROM == cq,
                    POS >= sq,
                    POS <= eq)  %>%
      dplyr::left_join(strain_list) %>%
      dplyr::filter(!STRAIN=="") %>%
      dplyr::filter(!is.na(strains_used)) %>%
      dplyr::select(CHROM,POS,STRAIN,strains_used) %>%
      dplyr::distinct() %>%
      dplyr::group_by(CHROM,POS,STRAIN,strains_used)
    annotations_region$x = 0   

    for (i in 1:nrow(annotations_region)){
        strainlist = strsplit(annotations_region$strains_used[i], ",")[[1]]
        if (ann_type == "bcsq"){
            strains = strsplit(annotations_region$STRAIN[i], " ")[[1]];
        } else { 
            strains = strsplit(annotations_region$STRAIN[i], ",")[[1]];
        }
        for (strain in strains) {
            if (strain %in% strainlist) annotations_region$x[i] = annotations_region$x[i] + 1
        }
    }

    annotations_region <- annotations_region %>%
      dplyr::filter(x>0) %>% #find overlap > 0
      dplyr::ungroup() %>%
      dplyr::group_by(CHROM,POS) %>%
      dplyr::filter(x==max(x))%>%   # find the max overlap for each position. Multiple max might exist
      dplyr::select(-strains_used,-x) %>%
      dplyr::left_join(annotations) %>% {
        if(ann_type != "snpeff") dplyr::select(., -STRAIN) else .
      }%>%
      tidyr::unite(marker, CHROM, POS, sep = "_")

    # pull variants from finemap impute -- don't need this?
    # impute_vcf <- data.table::fread(args[4]) %>%
    #     dplyr::select(marker:POS, REF:strains)

    # filter annotations to include variants within region
    # this excludes variants that are not annotated -- variants not within a gene
    annotation_out[[r]] <- annotations %>%
        dplyr::filter(CHROM == cq,
                      POS >= sq,
                      POS <= eq) %>% {
                          if(ann_type != "snpeff") dplyr::select(., -STRAIN) else .
                      } %>%
        # tidyr::separate_rows(Strains) %>% # is this going to crash R?
        # dplyr::filter(Strains %in% strains) %>%
        # dplyr::group_by_at(vars(-Strains)) %>%
        # dplyr::summarize(Strains = paste(Strains, collapse = ",")) %>%
        # dplyr::rename(strains = Strains) %>%
        tidyr::unite(marker, CHROM, POS, sep = "_") %>%
      dplyr::filter(!marker %in% annotations_region$marker) %>%
      dplyr::bind_rows(annotations_region)

    # add variants with no annotations???

}

# combine annotations for regions
annotation_out = dplyr::bind_rows(annotation_out)
if(ann_type != "snpeff") {
    annotation_out = tibble::add_column(annotation_out, dplyr::select(annotation_out, WBGENE) %>% dplyr::rename(WORMBASE_ID = WBGENE))
}

annotation_df <- annotation_out %>%
    dplyr::select(., -DIVERGENT) %>%
    dplyr::distinct() %>%
    dplyr::left_join(pr_trait_ld %>% filter(allele == "ALT"), ., by = c("marker", "REF", "ALT")) %>% {
        if(ann_type != "snpeff") dplyr::rename(., gene_id = WBGENE) else dplyr::select(., -GENE_NAME)
    }

genes_in_region <- gene_ref_flat %>%
    dplyr::filter(wbgene %in% annotation_df$gene_id) %>%
    dplyr::select(gene_id = wbgene, strand, txstart, txend, feature_id = gene, BIOTYPE = "biotype type") %>%
    dplyr::arrange(txstart, feature_id)%>%
    dplyr::distinct(gene_id, feature_id, .keep_all = TRUE)

ugly_genes_in_region <- genes_in_region %>%
    dplyr::left_join(annotation_df, ., by = "gene_id", relationship = "many-to-many") %>%
    dplyr::distinct(marker, CHR, POS, log10p, peak_marker, strains, .keep_all = T) %>% {
        if(ann_type == "snpeff") dplyr::rename(., VARIANT_IMPACT = impact) else .
    } %>%
    dplyr::mutate(start_pos = query_regions$start_pos,
                  end_pos = query_regions$end_pos)

tidy_genes_in_region <- if(ann_type != "snpeff") {
        # no feature_type
        ugly_genes_in_region %>%
            dplyr::select(MARKER = marker, CHROM = CHR, POS, REF, ALT, START_POS = start_pos, END_POS = end_pos, MAF_variant = maf_marker_b,
                      GENE_NAME, WBGeneID = gene_id,
                      WBFeature_ID = feature_id, TRANSCRIPT_BIOTYPE = BIOTYPE, CONSEQUENCE,
                      NUCLEOTIDE_CHANGE = DNACHANGE, AMINO_ACID_CHANGE=AA, BLOSUM = BLOSUM_SCORE, Grantham = GRANTHAM_SCORE, Percent_Protein=PERCENT_PROTEIN,
                      STRAND=strand, TRANSCRIPTION_START_POS = txstart, TRANSCRIPTION_END_POS = txend,
                      PEAK_MARKER = peak_marker, PEAK_MAF = peak_maf, STRAIN = strains,
                      VARIANT_LD_WITH_PEAK_MARKER = ld_r2, VARIANT_LOG10p = log10p, STRAIN_GENOTYPE = allele)
    } else {
        ugly_genes_in_region %>%
            dplyr::select(MARKER = marker, CHROM = CHR, POS, REF, ALT, START_POS = start_pos, END_POS = end_pos,MAF_variant = maf_marker_b,
                      GENE_NAME = GENE, WBGeneID = gene_id, WBFeature_TYPE = feature_type,
                      WBFeature_ID = feature_id.x, VARIANT_IMPACT,
                      NUCLEOTIDE_CHANGE = DNA_CHANGE, AMINO_ACID_CHANGE,
                      STRAND = strand, TRANSCRIPTION_START_POS = txstart, TRANSCRIPTION_END_POS = txend,
                      PEAK_MARKER = peak_marker, PEAK_MAF = peak_maf, STRAIN = strains, STRAIN_GENOTYPE = allele,
                      VARIANT_LD_WITH_PEAK_MARKER = ld_r2, VARIANT_LOG10p = log10p)
        }

#write_tsv(tidy_genes_in_region,
 #         path = glue::glue("{analysis_trait}_{cq}_{sq}-{eq}_{ann_type}_genes_{args[5]}.tsv"))

write_tsv(tidy_genes_in_region,
          path = glue::glue("{analysis_trait}_{cq}_{output_sq}-{output_eq}_{ann_type}_genes_{args[5]}.tsv"))


for(r in 1:length(unique(ugly_genes_in_region$start_pos))){

    gene_df <- ugly_genes_in_region %>%
        dplyr::filter(start_pos == unique(ugly_genes_in_region$start_pos)[r]) %>%
        dplyr::distinct(gene_id, peak_marker, CHR, strand, txstart, txend, start_pos, end_pos, log10p) %>%
        dplyr::group_by(gene_id) %>%
        dplyr::mutate(log10p = max(log10p, na.rm = T)) %>%
        dplyr::distinct()

    #peak_variant <- as.numeric(strsplit(unique(gene_df$peak_marker), split = ":")[[1]][2])
    peak_df <- gene_df %>%
    dplyr::filter(!is.na(peak_marker))

  peak_variant <- ifelse(identical(gsub("_", ":", unique(peak_df$peak_marker)), character(0)), NA,
                         as.numeric(strsplit(unique(peak_df$peak_marker), split = ":")[[1]][2]))

    variant_df <- ugly_genes_in_region %>%
        dplyr::filter(start_pos == unique(ugly_genes_in_region$start_pos)[r]) %>%
        dplyr::distinct(CHR, POS, log10p, CONSEQUENCE)

    if(ann_type == "bcsq") {
        variant_df$CONSEQUENCE[variant_df$CONSEQUENCE == "N/A"] <- "intergenic"
    } else {
        variant_df$VARIANT_IMPACT[is.na(variant_df$VARIANT_IMPACT)] <- "Intergenic"
    }

    xs <- unique(gene_df$start_pos)
    xe <- unique(gene_df$end_pos)
    cq <- unique(gene_df$CHR)

    max_logp <- unique(max(variant_df$log10p, na.rm = T))/150

    if(ann_type == "bcsq") {
      gene_plot <- ggplot(gene_df) +
        geom_vline(aes(xintercept = peak_variant/1e6),
                   linetype=3, color = "cyan")+
        geom_segment(aes(x = ifelse(strand == "+", txstart/1e6, txend/1e6),
                         xend = ifelse(strand == "+", txend/1e6, txstart/1e6),
                         y = log10p,
                         yend = log10p),
                     arrow = arrow(length = unit(5, "points")), size = 1) +
        geom_segment(aes(x = POS/1e6,
                         xend = POS/1e6,
                         y = log10p+max_logp,
                         yend = log10p-max_logp,
                         color = CONSEQUENCE), data = variant_df) +
        # remove moderate and modifier - we prob won't run snpeff anyways?
        scale_color_manual(values = c("synonymous" = "green",
                                      "missense" = "red",
                                      "splice_region" = "red",
                                      "intergenic" = "gray80",
                                      "intron" = "gray30",
                                      "5_prime_utr" = "orange",
                                      "3_prime_utr" = "orange",
                                      "Linker" = "gray80"),
                           breaks = c("intergenic", "synonymous", "intron", "3_prime_utr", "5_prime_utr", "missense", "splice_region"),
                           name = "EFFECT")+
        labs(x = "Genomic Position (Mb)",
             y = expression(-log[10](italic(p))))+
        theme_bw(18)+
        xlim(c(xs/1e6, xe/1e6)) +
        theme(legend.position = "top",
              panel.grid = element_blank())
    } else {
      gene_plot <- ggplot(gene_df) +
        geom_vline(aes(xintercept = peak_variant/1e6),
                   linetype=3, color = "cyan")+
        geom_segment(aes(x = ifelse(strand == "+", txstart/1e6, txend/1e6),
                         xend = ifelse(strand == "+", txend/1e6, txstart/1e6),
                         y = log10p,
                         yend = log10p),
                     arrow = arrow(length = unit(5, "points")), size = 1) +
        geom_segment(aes(x = POS/1e6,
                         xend = POS/1e6,
                         y = log10p+max_logp,
                         yend = log10p-max_logp,
                         color = VARIANT_IMPACT), data = variant_df) +
        # remove moderate and modifier - we prob won't run snpeff anyways?
        scale_color_manual(values = c("MODIFIER" = "gray50",
                                      "LOW" = "gray30",
                                      "MODERATE" = "orange",
                                      "HIGH" = "red",
                                      "Intergenic" = "gray80"),
                           breaks = c("HIGH", "MODERATE", "LOW", "MODIFIER", "Intergenic"),
                           name = "EFFECT")+
        labs(x = "Genomic Position (Mb)",
             y = expression(-log[10](italic(p))))+
        theme_bw(18)+
        xlim(c(xs/1e6, xe/1e6)) +
        theme(legend.position = "top",
              panel.grid = element_blank())
    }


    ggsave(gene_plot,
           filename = glue::glue("{analysis_trait}_{cq}_{xs}-{xe}_gene_plot_{ann_type}_{args[5]}.pdf"),
           height=10, width = 14)
}

