process summarize_mapping {

    label "xs"

    publishDir "${params.out}/INBRED/Plots", mode: 'copy', pattern: "*_inbred.pdf"
    publishDir "${params.out}/LOCO/Plots", mode: 'copy', pattern: "*_loco.pdf"

    input:
        tuple file(qtl_peaks_inbred), file(qtl_peaks_loco), file(chr_lens), file(summarize_mapping_file)

    output:
        file("Summarized_mappings*.pdf")

    """
    Rscript --vanilla ${summarize_mapping_file} ${qtl_peaks_inbred} ${chr_lens} inbred
    Rscript --vanilla ${summarize_mapping_file} ${qtl_peaks_loco} ${chr_lens} loco
    """
}

process generate_plots {

    label "xs"

    publishDir "${params.out}/INBRED/Plots/LDPlots", mode: 'copy', pattern: "*_LD_inbred.plot.png"
    publishDir "${params.out}/INBRED/Plots/EffectPlots", mode: 'copy', pattern: "*_effect_inbred.plot.png"
    publishDir "${params.out}/INBRED/Plots/ManhattanPlots", mode: 'copy', pattern: "*_manhattan_inbred.plot.png"

    publishDir "${params.out}/LOCO/Plots/LDPlots", mode: 'copy', pattern: "*_LD_loco.plot.png"
    publishDir "${params.out}/LOCO/Plots/EffectPlots", mode: 'copy', pattern: "*_effect_loco.plot.png"
    publishDir "${params.out}/LOCO/Plots/ManhattanPlots", mode: 'copy', pattern: "*_manhattan_loco.plot.png"

    input:
        tuple val(TRAIT), file(geno), file(pheno), file(tests), file(aggregate_mapping_inbred), \
              file(aggregate_mapping_loco), file(pipeline_plotting_mod) 

    output:
        tuple file(geno), file(pheno), file(aggregate_mapping_inbred), val(TRAIT), emit: maps_from_plot_inbred
        tuple file(geno), file(pheno), file(aggregate_mapping_loco), val(TRAIT), emit: maps_from_plot_loco
        file("*.png")

    """
    Rscript --vanilla ${pipeline_plotting_mod} ${aggregate_mapping_inbred} ${tests} inbred
    Rscript --vanilla ${pipeline_plotting_mod} ${aggregate_mapping_loco} ${tests} loco
    """
}


process LD_between_regions {

    label 'xs'

    publishDir "${params.out}/INBRED/Mapping/Processed", mode: 'copy', pattern: "*LD_between_QTL_regions_inbred.tsv"
    publishDir "${params.out}/LOCO/Mapping/Processed", mode: 'copy', pattern: "*LD_between_QTL_regions_loco.tsv"

    input:
        tuple val(TRAIT), file(geno), file(pheno), file(tests), file(aggregate_mapping_inbred), \
              file(aggregate_mapping_loco), file(ld_between_regions)

    output:
        tuple val(TRAIT), path("*LD_between_QTL_regions_inbred.tsv"), path("*LD_between_QTL_regions_loco.tsv") optional true

    """
    Rscript --vanilla ${ld_between_regions} ${geno} ${aggregate_mapping_inbred} ${TRAIT} inbred 
    Rscript --vanilla ${ld_between_regions} ${geno} ${aggregate_mapping_loco} ${TRAIT} loco
    """
}


/*
======================================
~ > *                            * < ~
~ ~ > *                        * < ~ ~
~ ~ ~ > *  RUN FINE MAPPING  * < ~ ~ ~
~ ~ > *                        * < ~ ~
~ > *                            * < ~
======================================
*/

/*
------------ Extract QTL interval genotype matrix
*/

process prep_ld_files {

    tag {TRAIT}
    label "prep_ld_files"
    errorStrategy 'retry'
    time { 20.minute * task.attempt }
    cpus = { 2 * task.attempt }
    memory = { 8.GB * task.attempt }

    publishDir "${params.out}/INBRED/Fine_Mappings/Data", mode: 'copy', pattern: "*ROI_Genotype_Matrix_inbred.tsv"
    publishDir "${params.out}/INBRED/Fine_Mappings/Data", mode: 'copy', pattern: "*LD_inbred.tsv"

    publishDir "${params.out}/LOCO/Fine_Mappings/Data", mode: 'copy', pattern: "*ROI_Genotype_Matrix_loco.tsv"
    publishDir "${params.out}/LOCO/Fine_Mappings/Data", mode: 'copy', pattern: "*LD_loco.tsv"

    input:
        tuple val(TRAIT), val(CHROM), val(marker), val(log10p), val(start_pos), val(peak_pos), \
              val(end_pos), val(peak_id), val(h2), val(algorithm), file(geno), file(pheno), \
              file(aggregate_mapping), file(imputed_vcf), file(imputed_index), file(phenotype), \
              file(num_chroms)

    output:
        tuple val(TRAIT), file(pheno), file("*ROI_Genotype_Matrix*.tsv"), file("*LD*.tsv"), file("*.bim"), file("*.bed"), file("*.fam"), val(algorithm), emit: finemap_preps
        tuple val(TRAIT), file("*ROI_Genotype_Matrix_${algorithm}.tsv"), emit: finemap_ROI
        tuple val(TRAIT), file("*LD_${algorithm}.tsv"), emit: finemap_LD

    """
    cat ${aggregate_mapping} |\\
    awk '\$0 !~ "\\tNA\\t" {print}' |\\
    awk '!seen[\$1,\$12,\$19,\$20,\$21]++' |\\
    awk 'NR>1{print \$1, \$11, \$18, \$19, \$20}' OFS="\\t" > ${TRAIT}_QTL_peaks.tsv
    filename='${TRAIT}_QTL_peaks.tsv'
    echo Start
    while read p; do 
        chromosome=`echo \$p | cut -f1 -d ' '`
        trait=`echo \$p | cut -f2 -d ' '`
        start_pos=`echo \$p | cut -f3 -d ' '`
        peak_pos=`echo \$p | cut -f4 -d ' '`
        end_pos=`echo \$p | cut -f5 -d ' '`
        if [ \$chromosome == "MtDNA" ]; then
            cat ${phenotype} | awk '\$0 !~ "strain" {print}' | cut -f1 > phenotyped_samples.txt
            bcftools view --regions \$chromosome:\$start_pos-\$end_pos ${imputed_vcf} \
            -S phenotyped_samples.txt |\\
            bcftools filter -i N_MISSING=0 |\\
            bcftools annotate --rename-chrs ${num_chroms} |\\
            awk '\$0 !~ "#" {print \$1":"\$2}' > \$trait.\$chromosome.\$start_pos.\$end_pos.txt
            bcftools view --regions \$chromosome:\$start_pos-\$end_pos ${imputed_vcf} \
            -S phenotyped_samples.txt |\\
            bcftools filter -i N_MISSING=0 -Oz --threads ${task.cpus} |\\
            bcftools annotate --rename-chrs ${num_chroms} -o finemap.vcf.gz
        else
            cat ${phenotype} | awk '\$0 !~ "strain" {print}' | cut -f1 > phenotyped_samples.txt
            bcftools view --regions \$chromosome:\$start_pos-\$end_pos ${imputed_vcf} \
            -S phenotyped_samples.txt |\\
            bcftools filter -i N_MISSING=0 |\\
            bcftools filter -e 'GT="het"' |\\
            bcftools annotate --rename-chrs ${num_chroms} |\\
            awk '\$0 !~ "#" {print \$1":"\$2}' > \$trait.\$chromosome.\$start_pos.\$end_pos.txt
            bcftools view --regions \$chromosome:\$start_pos-\$end_pos ${imputed_vcf} \
            -S phenotyped_samples.txt |\\
            bcftools filter -i N_MISSING=0 -Oz --threads ${task.cpus} |\\
            bcftools annotate --rename-chrs ${num_chroms} -o finemap.vcf.gz
        fi
        plink --vcf finemap.vcf.gz \\
            --threads ${task.cpus} \\
            --snps-only \\
            --maf ${params.maf} \\
            --biallelic-only \\
            --allow-extra-chr \\
            --set-missing-var-ids @:# \\
            --geno \\
            --make-bed \\
            --recode vcf-iid bgz \\
            --extract \$trait.\$chromosome.\$start_pos.\$end_pos.txt \\
            --out \$trait.\$chromosome.\$start_pos.\$end_pos
        nsnps=`wc -l \$trait.\$chromosome.\$start_pos.\$end_pos.txt | cut -f1 -d ' '`
        echo `cat ${num_chroms} | grep -w \$chromosome`
        chrom_num=`cat ${num_chroms} | grep -w \$chromosome | sed "s/\\t/ /" | cut -f2 -d ' '`
        plink --r2 with-freqs \\
            --threads ${task.cpus} \\
            --allow-extra-chr \\
            --snps-only \\
            --ld-window-r2 0 \\
            --ld-snp \$chrom_num:\$peak_pos \\
            --ld-window \$nsnps \\
            --ld-window-kb 6000 \\
            --chr \$chrom_num \\
            --out \$trait.\$chromosome:\$start_pos-\$end_pos.QTL \\
            --set-missing-var-ids @:# \\
            --vcf \$trait.\$chromosome.\$start_pos.\$end_pos.vcf.gz
        cut \$trait.\$chromosome:\$start_pos-\$end_pos.QTL.ld -f2-10 > \$trait.\$chromosome.\$start_pos.\$end_pos.LD_${algorithm}.tsv
        bcftools query --print-header -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' finemap.vcf.gz |\\
            sed 's/[[# 0-9]*]//g' |\\
            sed 's/:GT//g' |\\
            sed 's/0|0/-1/g' |\\
            sed 's/1|1/1/g' |\\
            sed 's/0|1/NA/g' |\\
            sed 's/1|0/NA/g' |\\
            sed 's/.|./NA/g'  |\\
            sed 's/0\\/0/-1/g' |\\
            sed 's/1\\/1/1/g'  |\\
            sed 's/0\\/1/NA/g' |\\
            sed 's/1\\/0/NA/g' |\\
            sed 's/.\\/./NA/g' |\\
            sed 's/^23/X/g' > \$trait.\$chromosome.\$start_pos-\$end_pos.ROI_Genotype_Matrix_${algorithm}.tsv
    done < \$filename
    """
}


process gcta_fine_maps {

    tag {TRAIT}
    label "gcta_fine_maps"
    errorStrategy 'retry'
    time { 40.minute * task.attempt }
    cpus { task.attempt }
    memory { 4.GB * task.attempt }

    publishDir "${params.out}/INBRED/Fine_Mappings/Data", mode: 'copy', pattern: "*inbred.fastGWA"
    publishDir "${params.out}/INBRED/Fine_Mappings/Data", mode: 'copy', pattern: "*_genes_inbred.tsv"
    publishDir "${params.out}/INBRED/Fine_Mappings/Plots", mode: 'copy', pattern: "*inbred.pdf"

    publishDir "${params.out}/LOCO/Fine_Mappings/Data", mode: 'copy', pattern: "*loco.fastGWA"
    publishDir "${params.out}/LOCO/Fine_Mappings/Data", mode: 'copy', pattern: "*_genes_loco.tsv"
    publishDir "${params.out}/LOCO/Fine_Mappings/Plots", mode: 'copy', pattern: "*loco.pdf"

    input:
        tuple val(TRAIT), file(pheno), file(ROI_geno), file(ROI_LD), file(bim), file(bed), file(fam), val(algorithm), \
              file(formatted_traits), file(trait_bed), file(trait_bim), file(trait_fam), \
              file(trait_map), file(trait_nosex), file(trait_ped), file(grm_bin), file(grm_id), \
              file(grm_nbin), file(grm_bin_inbred), file(grm_id_inbred), \
              file(annotation), file(genefile), file(finemap_qtl_intervals), file(plot_genes), val(species)

    output:
        tuple file("*.fastGWA"), val(TRAIT), file("*.prLD_df*.tsv"), file("*.pdf"), file("*_genes*.tsv"), val(algorithm)
        tuple file("*_genes*.tsv"), val(TRAIT), val(algorithm), emit: finemap_done
        tuple val(TRAIT), file("*inbred.fastGWA"), emit: finemap_GWA_inbred, optional: true
        tuple val(TRAIT), file("*.prLD_df_inbred.tsv"), emit: finemap_df_inbred, optional: true
        tuple val(TRAIT), file("*_genes_inbred.tsv"), emit: finemap_genes_inbred, optional: true
        tuple val(TRAIT), file("*loco.fastGWA"),  emit: finemap_GWA_loco, optional: true
        tuple val(TRAIT), file("*.prLD_df_loco.tsv"), emit: finemap_df_loco, optional: true
        tuple val(TRAIT), file("*_genes_loco.tsv"), emit: finemap_genes_loco, optional: true
    
    """
    tail -n +2 ${pheno} | awk 'BEGIN {OFS="\\t"}; {print \$1, \$1, \$2}' > plink_finemap_traits.tsv
    for i in *ROI_Genotype_Matrix_${algorithm}.tsv;
        do
        chr=`echo \$i | cut -f2 -d "."`
        start=`echo \$i | cut -f3 -d "." | cut -f1 -d "-"`
        stop=`echo \$i | cut -f3 -d "." | cut -f2 -d "-"`
        gcta64 --bfile ${TRAIT} \\
                --autosome \\
                --maf ${params.maf} \\
                --make-grm-inbred \\
                --out ${TRAIT}.FM_grm_inbred.${algorithm} \\
                --thread-num ${task.cpus}
        gcta64 --grm ${TRAIT}.FM_grm_inbred.${algorithm} \\
                --make-bK-sparse ${params.sparse_cut} \\
                --out ${TRAIT}.sparse_FM_grm_inbred.${algorithm}  \\
                --thread-num ${task.cpus}
        gcta64 --grm ${TRAIT}.FM_grm_inbred.${algorithm} \\
                --pca 1 \\
                --out ${TRAIT}.sparse_FM_grm_inbred.${algorithm}  \\
                --thread-num ${task.cpus}
        gcta64 --fastGWA-lmm-exact \\
                --grm-sparse ${TRAIT}.sparse_FM_grm_inbred.${algorithm} \\
                --bfile ${TRAIT}.\$chr.\$start.\$stop \\
                --qcovar ${TRAIT}.sparse_FM_grm_inbred.${algorithm}.eigenvec \\
                --out ${TRAIT}.\$chr.\$start.\$stop.finemap_inbred.${algorithm} \\
                --pheno plink_finemap_traits.tsv \\
                --maf ${params.maf} \\
                --thread-num ${task.cpus}
        
        Rscript --vanilla ${finemap_qtl_intervals} ${TRAIT}.\$chr.\$start.\$stop.finemap_inbred.${algorithm}.fastGWA \$i ${TRAIT}.\$chr.\$start.\$stop.LD_${algorithm}.tsv ${algorithm}
        
        Rscript --vanilla ${plot_genes} ${TRAIT}.\$chr.\$start.\$stop.prLD_df_${algorithm}.tsv ${pheno} ${genefile} ${annotation} ${algorithm} ${species}
    done
    """
}


/*
------ Slice out the QTL region for plotting divergent region and haplotype data.
*/


process divergent_and_haplotype {

    label "xs"

    publishDir "${params.out}/INBRED/Divergent_and_haplotype", mode: 'copy', pattern: '*inbred*'
    publishDir "${params.out}/LOCO/Divergent_and_haplotype", mode: 'copy', pattern: '*loco*'

    input:
        tuple file(QTL_peaks), val(algorithm), file(divergent_bins), file(divergent_df_isotype), \
              file(haplotype_df_isotype), file(div_isotype_list)

    output:
    tuple file("all_QTL_bins_inbred.bed"), file("all_QTL_div_inbred.bed"), file("haplotype_in_QTL_region_inbred.txt"), file("div_isotype_list2_inbred.txt"), emit: div_hap_table_inbred, optional: true 
    tuple file("all_QTL_bins_loco.bed"), file("all_QTL_div_loco.bed"), file("haplotype_in_QTL_region_loco.txt"), file("div_isotype_list2_loco.txt"), emit: div_hap_table_loco, optional: true 

    """
    awk NR\\>1 ${QTL_peaks} | awk -v OFS='\t' '{print \$1,\$5,\$7}' > QTL_region_${algorithm}.bed
    bedtools intersect -wa -a ${divergent_bins} -b QTL_region_${algorithm}.bed | sort -k1,1 -k2,2n | uniq > all_QTL_bins_${algorithm}.bed
    bedtools intersect -a ${divergent_df_isotype} -b QTL_region_${algorithm}.bed | sort -k1,1 -k2,2n | uniq > all_QTL_div_${algorithm}.bed
    bedtools intersect -a ${haplotype_df_isotype} -b QTL_region_${algorithm}.bed -wo | sort -k1,1 -k2,2n | uniq > haplotype_in_QTL_region_${algorithm}.txt
    cp ${div_isotype_list} ./div_isotype_list2_${algorithm}.txt
    """

}

// generate trait-specific html reports
process html_report_main {

    tag {"${TRAIT} - HTML REPORT" }
    label "html_report_main"
    errorStrategy 'retry'
    time { 40.minute * task.attempt }
    cpus { 2 * task.attempt }
    memory { 8.GB * task.attempt }

    publishDir "${params.out}/Reports/scripts/", pattern: "*.Rmd", overwrite: true, mode: 'copy'
    publishDir "${params.out}/Reports", pattern: "*.html", overwrite: true, mode: 'copy'

    input:
        // there is no way this will work for every trait... UGH
        tuple   val(TRAIT), \
                file(pheno), \
                file(qtl_peaks_inbred), \
                file(qtl_peaks_loco), \
                file(strain_issues), \
                file(tests), \
                file(geno), \
                file(ns_report_md), \
                file(ns_report_template_md), \
                file(render_markdown), \
                file(ns_report_alg), \
                val(mediate), \
                val(species), \
//divergent_and_haplotypes
                file(qtl_bins_inbred), \
                file(qtl_div_inbred), \
                file(haplotype_qtl_inbred), \
                file(div_isotype_inbred), \
//divergent_and_haplotypes
                file(qtl_bins_loco), \
                file(qtl_div_loco), \
                file(haplotype_qtl_loco), \
                file(div_isotype_loco), \
//gcta_interval_maps
                    file(pmap_inbred), \
                    file(pmap_loco), \
//gcta_find_maps - may be missing
                    file(fastGWA_inbred), \
                    file(prLD_inbred), \
                    file(bcsq_genes_inbred), \
//gcta_find_maps - may be missing
                    file(fastGWA_loco), \
                    file(prLD_loco), \
                    file(bcsq_genes_loco), \
//prep_ld_files - may be missing
                    file(roi_geno_inbred), \
                    file(roi_ld_inbred), \
//prep_ld_files - may be missing
                    file(roi_geno_loco), \
                    file(roi_ld_loco),\
 //summary_mediation - may be missing
                    file(mediation_summary_inbred), \
                    file(mediation_summary_loco)

    output:
        tuple file("NemaScan_Report_*.Rmd"), file("NemaScan_Report_*.html")


    """
    # edit the file paths for generating these reports
    cat ${ns_report_md} | \\
    sed "s+TRAIT_NAME_HOLDER+${TRAIT}+g" | \\
    sed "s+SPECIES+${species}+g" | \\
    sed "s+MEDIATION+${mediate}+g" | \\
    sed "s+../../Phenotypes/strain_issues.txt+${strain_issues}+g" | \\
    sed "s+../../Genotype_Matrix/total_independent_tests.txt+${tests}+g" | \\
    sed 's+paste0("../../INBRED/Mapping/Processed/processed_",trait_name,"_AGGREGATE_mapping_inbred.tsv")+"${pmap_inbred}"+g' | \\
    sed 's+paste0("../../LOCO/Mapping/Processed/processed_",trait_name,"_AGGREGATE_mapping_loco.tsv")+"${pmap_loco}"+g' | \\
    sed "s+inbred <- knitr::knit_child('NemaScan_Report_algorithm_template.Rmd')+inbred <- knitr::knit_child('NemaScan_Report_INBRED_template.Rmd')+g" | \\
    sed "s+loco <- knitr::knit_child('NemaScan_Report_algorithm_template.Rmd')+loco <- knitr::knit_child('NemaScan_Report_LOCO_template.Rmd')+g" > NemaScan_Report_${TRAIT}_main.Rmd

    # this is so complicating... ugh
    cat ${ns_report_alg} | \\
    sed 's+glue::glue("../../{alg}/Mapping/Processed/QTL_peaks_{stringr::str_to_lower(alg)}.tsv")+"${qtl_peaks_inbred}"+g' | \\
    sed "s+../../Genotype_Matrix/Genotype_Matrix.tsv+${geno}+g" | \\
    sed "s+NemaScan_Report_region_template.Rmd+NemaScan_Report_region_INBRED_${TRAIT}.Rmd+g" > NemaScan_Report_INBRED_template.Rmd

    cat ${ns_report_alg} | \\
    sed 's+glue::glue("../../{alg}/Mapping/Processed/QTL_peaks_{stringr::str_to_lower(alg)}.tsv")+"${qtl_peaks_loco}"+g' | \\
    sed "s+../../Genotype_Matrix/Genotype_Matrix.tsv+${geno}+g" | \\
    sed "s+NemaScan_Report_region_template.Rmd+NemaScan_Report_region_LOCO_${TRAIT}.Rmd+g" > NemaScan_Report_LOCO_template.Rmd

    cat "${ns_report_template_md}" | \\
    sed 's+../../{alg}/Fine_Mappings/Data/++g' | \\
    sed 's+glue::glue("../../{alg}/Divergent_and_haplotype/all_QTL_bins_{stringr::str_to_lower(alg)}.bed")+"${qtl_bins_inbred}"+g' | \\
    sed 's+glue::glue("../../{alg}/Divergent_and_haplotype/all_QTL_div_{stringr::str_to_lower(alg)}.bed")+"${qtl_div_inbred}"+g' | \\
    sed 's+../../{alg}/Mediation/++g' | \\
    sed 's+glue::glue("../../{alg}/Divergent_and_haplotype/haplotype_in_QTL_region_{stringr::str_to_lower(alg)}.txt")+"${haplotype_qtl_inbred}"+g' > NemaScan_Report_region_INBRED_${TRAIT}.Rmd

    cat "${ns_report_template_md}" | \\
    sed 's+../../{alg}/Fine_Mappings/Data/++g' | \\
    sed 's+glue::glue("../../{alg}/Divergent_and_haplotype/all_QTL_bins_{stringr::str_to_lower(alg)}.bed")+"${qtl_bins_loco}"+g' | \\
    sed 's+glue::glue("../../{alg}/Divergent_and_haplotype/all_QTL_div_{stringr::str_to_lower(alg)}.bed")+"${qtl_div_loco}"+g' | \\
    sed 's+../../{alg}/Mediation/++g' | \\
    sed 's+glue::glue("../../{alg}/Divergent_and_haplotype/haplotype_in_QTL_region_{stringr::str_to_lower(alg)}.txt")+"${haplotype_qtl_loco}"+g' > NemaScan_Report_region_LOCO_${TRAIT}.Rmd

    Rscript --vanilla ${render_markdown} NemaScan_Report_${TRAIT}_main.Rmd
    """
}
