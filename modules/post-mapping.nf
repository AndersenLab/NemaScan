process summarize_mapping {

    publishDir "${params.out}/Plots", mode: 'copy'

    input:
        tuple file(qtl_peaks), file(chr_lens), file(summarize_mapping_file)

    output:
        file("Summarized_mappings.pdf")

    """
    echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${summarize_mapping_file} > summarize_mapping_file
    Rscript --vanilla summarize_mapping_file ${qtl_peaks} ${chr_lens}

    """
}


process generate_plots {


    publishDir "${params.out}/Plots/LDPlots", mode: 'copy', pattern: "*_LD.plot.png"
    publishDir "${params.out}/Plots/EffectPlots", mode: 'copy', pattern: "*_effect.plot.png"
    publishDir "${params.out}/Plots/ManhattanPlots", mode: 'copy', pattern: "*_manhattan.plot.png"

    input:
        tuple file(geno), file(pheno), val(TRAIT), file(tests), file(aggregate_mapping), file(pipeline_plotting_mod)

    output:
        tuple file(geno), file(pheno), file(aggregate_mapping), val(TRAIT),  emit: maps_from_plot
        file("*.png")

    """
    echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${pipeline_plotting_mod} > pipeline.plotting.mod
    Rscript --vanilla pipeline.plotting.mod ${aggregate_mapping} ${tests}
    """
}


process LD_between_regions {

  publishDir "${params.out}/Mapping/Processed", mode: 'copy', pattern: "*LD_between_QTL_regions.tsv"

  input:
        tuple file(geno), file(pheno), val(TRAIT), file(tests), file(aggregate_mapping), file(ld_between_regions)

  output:
        tuple val(TRAIT), path("*LD_between_QTL_regions.tsv") optional true
        val TRAIT, emit: linkage_done

  """
    echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${ld_between_regions} > LD_between_regions 
    Rscript --vanilla LD_between_regions ${geno} ${aggregate_mapping} ${TRAIT}
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

    // machineType 'n1-standard-4'
    label 'med'

    tag {TRAIT}

    publishDir "${params.out}/Fine_Mappings/Data", mode: 'copy', pattern: "*ROI_Genotype_Matrix.tsv"
    publishDir "${params.out}/Fine_Mappings/Data", mode: 'copy', pattern: "*LD.tsv"

    input:
        tuple val(TRAIT), val(CHROM), val(marker), val(log10p), val(start_pos), val(peak_pos), val(end_pos), val(peak_id), val(h2), file(geno), \
        file(pheno), file(aggregate_mapping), file(imputed_vcf), file(imputed_index), file(phenotype), file(num_chroms)

    output:
        tuple val(TRAIT), file(pheno), file("*ROI_Genotype_Matrix.tsv"), file("*LD.tsv"), file("*.bim"), file("*.bed"), file("*.fam"), emit: finemap_preps
        tuple val(TRAIT), file("*ROI_Genotype_Matrix.tsv"), file("*LD.tsv"), emit: finemap_LD

    """
        echo "HELLO"
        cat ${aggregate_mapping} |\\
        awk '\$0 !~ "\\tNA\\t" {print}' |\\
        awk '!seen[\$1,\$12,\$19,\$20,\$21]++' |\\
        awk 'NR>1{print \$1, \$11, \$18, \$19, \$20}' OFS="\\t" > ${TRAIT}_QTL_peaks.tsv
        filename='${TRAIT}_QTL_peaks.tsv'
        echo Start
        while read p; do 
            chromosome=`echo \$p | cut -f1 -d' '`
            trait=`echo \$p | cut -f2 -d' '`
            start_pos=`echo \$p | cut -f3 -d' '`
            peak_pos=`echo \$p | cut -f4 -d' '`
            end_pos=`echo \$p | cut -f5 -d' '`
        
        cat ${phenotype} | awk '\$0 !~ "strain" {print}' | cut -f1 > phenotyped_samples.txt
        bcftools view --regions \$chromosome:\$start_pos-\$end_pos ${imputed_vcf} \
        -S phenotyped_samples.txt |\\
        bcftools filter -i N_MISSING=0 |\\
        bcftools annotate --rename-chrs ${num_chroms} |\\
        awk '\$0 !~ "#" {print \$1":"\$2}' > \$trait.\$chromosome.\$start_pos.\$end_pos.txt
        bcftools view --regions \$chromosome:\$start_pos-\$end_pos ${imputed_vcf} \
        -S phenotyped_samples.txt |\\
        bcftools filter -i N_MISSING=0 -Oz --threads 5 |\\
        bcftools annotate --rename-chrs ${num_chroms} -o finemap.vcf.gz
        plink --vcf finemap.vcf.gz \\
            --threads 5 \\
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
        nsnps=`wc -l \$trait.\$chromosome.\$start_pos.\$end_pos.txt | cut -f1 -d' '`
        chrom_num=`cat ${num_chroms} | grep -w \$chromosome | cut -f2 -d' '`
        plink --r2 with-freqs \\
            --threads 5 \\
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
        cut \$trait.\$chromosome:\$start_pos-\$end_pos.QTL.ld -f2-10 > \$trait.\$chromosome.\$start_pos.\$end_pos.LD.tsv
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
            sed 's/^23/X/g' > \$trait.\$chromosome:\$start_pos-\$end_pos.ROI_Genotype_Matrix.tsv
        done < \$filename
    """
}


process gcta_fine_maps {
    // machineType 'n2-highmem-8'
    label "highmem"

    publishDir "${params.out}/Fine_Mappings/Data", mode: 'copy', pattern: "*.fastGWA"
    publishDir "${params.out}/Fine_Mappings/Data", mode: 'copy', pattern: "*_genes.tsv"
    publishDir "${params.out}/Fine_Mappings/Plots", mode: 'copy', pattern: "*.pdf"

    input:
        tuple val(TRAIT), file(pheno), file(ROI_geno), file(ROI_LD), file(bim), file(bed), file(fam), file(annotation), file(genefile), file(finemap_qtl_intervals), file(plot_genes)

    output:
        tuple file("*.fastGWA"), val(TRAIT), file("*.prLD_df.tsv"), file("*.pdf"), file("*_genes.tsv")
        //val true, emit: finemap_done
        tuple file("*_genes.tsv"), val(TRAIT), emit: finemap_done
        tuple val(TRAIT), file("*.fastGWA"), file("*.prLD_df.tsv"), file("*_genes.tsv"), emit: finemap_html

    """
    tail -n +2 ${pheno} | awk 'BEGIN {OFS="\\t"}; {print \$1, \$1, \$2}' > plink_finemap_traits.tsv
    for i in *ROI_Genotype_Matrix.tsv;
      do
      chr=`echo \$i | cut -f2 -d "." | cut -f1 -d ":"`
      start=`echo \$i | cut -f2 -d "." | cut -f2 -d ":" | cut -f1 -d "-"`
      stop=`echo \$i | cut -f2 -d "." | cut -f2 -d ":" | cut -f2 -d "-"`
      gcta64 --bfile ${TRAIT}.\$chr.\$start.\$stop \\
              --autosome \\
              --maf ${params.maf} \\
              --make-grm-inbred \\
              --out ${TRAIT}.\$chr.\$start.\$stop.FM_grm_inbred \\
              --thread-num 9
      gcta64 --grm ${TRAIT}.\$chr.\$start.\$stop.FM_grm_inbred \\
              --make-bK-sparse ${params.sparse_cut} \\
              --out ${TRAIT}.\$chr.\$start.\$stop.sparse_FM_grm_inbred  \\
              --thread-num 9
      gcta64 --grm ${TRAIT}.\$chr.\$start.\$stop.FM_grm_inbred \\
              --pca 1 \\
              --out ${TRAIT}.\$chr.\$start.\$stop.sparse_FM_grm_inbred  \\
              --thread-num 9
      gcta64 --fastGWA-lmm-exact \\
              --grm-sparse ${TRAIT}.\$chr.\$start.\$stop.sparse_FM_grm_inbred \\
              --bfile ${TRAIT}.\$chr.\$start.\$stop \\
              --qcovar ${TRAIT}.\$chr.\$start.\$stop.sparse_FM_grm_inbred.eigenvec \\
              --out ${TRAIT}.\$chr.\$start.\$stop.finemap_inbred \\
              --pheno plink_finemap_traits.tsv \\
              --maf ${params.maf} \\
              --thread-num 9
      
      echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${finemap_qtl_intervals}  > Finemap_QTL_Intervals
      Rscript --vanilla Finemap_QTL_Intervals ${TRAIT}.\$chr.\$start.\$stop.finemap_inbred.fastGWA \$i ${TRAIT}.\$chr.\$start.\$stop.LD.tsv
      
      echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${plot_genes} > plot_genes 
      Rscript --vanilla plot_genes ${TRAIT}.\$chr.\$start.\$stop.prLD_df.tsv ${pheno} ${genefile} ${annotation}
    done
    """
}


/*
------ Slice out the QTL region for plotting divergent region and haplotype data.
*/


process divergent_and_haplotype {

  //executor 'local'

  publishDir "${params.out}/Divergent_and_haplotype", mode: 'copy'


  input:
    tuple file("QTL_peaks.tsv"), file("divergent_bins"), file(divergent_df_isotype), file(haplotype_df_isotype), file(div_isotype_list)

  output:
    tuple file("all_QTL_bins.bed"), file("all_QTL_div.bed"), file("haplotype_in_QTL_region.txt"), file("div_isotype_list2.txt"), emit: div_hap_table
    val true, emit: div_done


  """
  awk NR\\>1 QTL_peaks.tsv | awk -v OFS='\t' '{print \$1,\$5,\$7}' > QTL_region.bed
  bedtools intersect -wa -a ${divergent_bins} -b QTL_region.bed | sort -k1,1 -k2,2n | uniq > all_QTL_bins.bed
  bedtools intersect -a ${divergent_df_isotype} -b QTL_region.bed | sort -k1,1 -k2,2n | uniq > all_QTL_div.bed
  bedtools intersect -a ${haplotype_df_isotype} -b QTL_region.bed -wo | sort -k1,1 -k2,2n | uniq > haplotype_in_QTL_region.txt
  cp ${div_isotype_list} ./div_isotype_list2.txt
  """

}

// generate trait-specific html reports
process html_report_main {

  //executor 'local'
  tag {"${TRAIT} - HTML REPORT" }

  // machineType 'n1-highmem-2'
  label "highmem"

  publishDir "${params.out}/Reports", pattern: "*.Rmd", overwrite: true, mode: 'copy'
  publishDir "${params.out}/Reports", pattern: "*.html", overwrite: true, mode: 'copy'

  input:
    tuple val(TRAIT), file(qtl_peaks), file(pheno), file(strain_issues), file(tests), file(geno), file(ns_report_md), \
    file(ns_report_template_md), file(render_markdown), val(algorithm), val(mediate), val(species), file(qtl_bins), file(qtl_div), \
    file(haplotype_qtl), file(div_isotype), file(pmap), file(fastGWA), file(prLD), file(bcsq_genes), file(roi_geno), file(roi_ld), \
    file(mediation_summary)

  output:
    tuple file("NemaScan_Report_*.Rmd"), file("NemaScan_Report_*.html")


  """
    # edit the file paths for generating these reports
    cat ${ns_report_md} | \\
    sed "s+TRAIT_NAME_HOLDER+${TRAIT}+g" | \\
    sed "s+SPECIES+${species}+g" | \\
    sed "s+MEDIATION+${mediate}+g" | \\
    sed "s+Phenotypes/strain_issues.txt+${strain_issues}+g" | \\
    sed "s+Genotype_Matrix/total_independent_tests.txt+${tests}+g" | \\
    sed 's+paste0("Mapping/Processed/processed_",trait_name,"_AGGREGATE_mapping.tsv")+"${pmap}"+g' | \\
    sed "s+Mapping/Processed/QTL_peaks.tsv+${qtl_peaks}+g" | \\
    sed "s+Genotype_Matrix/Genotype_Matrix.tsv+${geno}+g" | \\
    sed "s+ALGORITHM+${algorithm}+" | \\
    sed "s+NemaScan_Report_region_template.Rmd+NemaScan_Report_region_${TRAIT}.Rmd+g" > NemaScan_Report_${TRAIT}_main.Rmd

    cat "${ns_report_template_md}" | \\
    sed 's+Fine_Mappings/Data/++g' | \\
    sed "s+Divergent_and_haplotype/div_isotype_list2.txt+${div_isotype}+g" | \\
    sed "s+Divergent_and_haplotype/all_QTL_bins.bed+${qtl_bins}+g" | \\
    sed "s+Divergent_and_haplotype/all_QTL_div.bed+${qtl_div}+g" | \\
    sed 's+glue::glue("Mediation/file_summary/{trait_name}_mediation.tsv")+"${mediation_summary}"+g' | \\
    sed "s+Divergent_and_haplotype/haplotype_in_QTL_region.txt+${haplotype_qtl}+g" > NemaScan_Report_region_${TRAIT}.Rmd
    
    echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" > .Rprofile
    
    echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${render_markdown}  > render_markdown 
    Rscript --vanilla render_markdown NemaScan_Report_${TRAIT}_main.Rmd
  """
}
