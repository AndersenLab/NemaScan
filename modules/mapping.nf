/*
======================================
~ > *                            * < ~
~ ~ > *                        * < ~ ~
~ ~ ~ > *  RUN GWAS MAPPING  * < ~ ~ ~
~ ~ > *                        * < ~ ~
~ > *                            * < ~
======================================
*/

/*
------------ GCTA
*/

process prepare_gcta_files {
    label "prepare_gcta_files"
    errorStrategy 'retry'
    time { 20.minute * task.attempt }
    cpus = { 4 * task.attempt }
    memory = { 16.GB * task.attempt }

    input:
        tuple file(strains), val(TRAIT), file(traits), file(vcf), file(index), file(num_chroms)

    output:
        tuple val(TRAIT), file("plink_formatted_traits.tsv"), file("${TRAIT}.bed"), file("${TRAIT}.bim"), file("${TRAIT}.fam"), file("${TRAIT}.map"), file("${TRAIT}.nosex"), file("${TRAIT}.ped")

    """
    bcftools annotate --rename-chrs ${num_chroms} ${vcf} |\\
    bcftools view -S ${strains} -Ou |\\
    bcftools filter -i N_MISSING=0 -Oz --threads 5 -o renamed_chroms.vcf.gz
    tabix -p vcf renamed_chroms.vcf.gz
    plink --vcf renamed_chroms.vcf.gz \\
          --threads ${task.cpus} \\
          --snps-only \\
          --biallelic-only \\
          --maf ${params.maf} \\
          --set-missing-var-ids @:# \\
          --indep-pairwise 50 10 0.8 \\
          --geno \\
          --allow-extra-chr
    tail -n +2 ${traits} | awk 'BEGIN {OFS="\\t"}; {print \$1, \$1, \$2}' > plink_formatted_traits.tsv
    plink --vcf renamed_chroms.vcf.gz \\
          --threads ${task.cpus} \\
          --make-bed \\
          --snps-only \\
          --biallelic-only \\
          --maf ${params.maf} \\
          --set-missing-var-ids @:# \\
          --extract plink.prune.in \\
          --geno \\
          --recode \\
          --out ${TRAIT} \\
          --allow-extra-chr \\
          --pheno plink_formatted_traits.tsv
    """
}

process gcta_grm {
 
    label "gcta_grm"
    errorStrategy 'retry'
    time { 20.minute * task.attempt }
    cpus = { 4 * task.attempt }
    memory = { 16.GB * task.attempt }

    input:
        tuple val(TRAIT), file(traits), file(bed), file(bim), file(fam), file(map), file(nosex), \
              file(ped)

    output:
        tuple val(TRAIT), file(traits), file(bed), file(bim), file(fam), file(map), file(nosex), \
              file(ped), file("${TRAIT}_gcta_grm.grm.bin"), file("${TRAIT}_gcta_grm.grm.id"), \
              file("${TRAIT}_gcta_grm.grm.N.bin"), file("${TRAIT}_gcta_grm_inbred.grm.bin"), \
              file("${TRAIT}_gcta_grm_inbred.grm.id")

    when:
        params.mapping

    """
    gcta64 --bfile ${TRAIT} \\
           --autosome \\
           --maf ${params.maf} \\
           --make-grm-inbred \\
           --out ${TRAIT}_gcta_grm_inbred \\
           --thread-num ${task.cpus}
    gcta64 --bfile ${TRAIT} \\
           --autosome \\
           --maf ${params.maf} \\
           --make-grm \\
           --out ${TRAIT}_gcta_grm \\
           --thread-num ${task.cpus}
    """
}


process gcta_lmm_exact_mapping {

    label "gcta_lmm_exact_mapping"
    errorStrategy 'retry'
    time { 20.minute * task.attempt }
    cpus = { 4 * task.attempt }
    memory = { 16.GB * task.attempt }

    publishDir "${params.out}/INBRED/Mapping/Raw", pattern: "*inbred_pca.fastGWA", overwrite: true
    publishDir "${params.out}/LOCO/Mapping/Raw", pattern: "*loco_pca.mlma", overwrite: true

    input:
    tuple val(TRAIT), file(traits), file(bed), file(bim), file(fam), file(map), \
          file(nosex), file(ped), file(grm_bin), file(grm_id), file(grm_nbin), \
           file(grm_bin_inbred), file(grm_id_inbred)
          

    output:
    tuple val(TRAIT), file("${TRAIT}_lmm-exact_inbred_pca.fastGWA"), file("${TRAIT}_lmm-exact_pca.loco.mlma")


    """
    gcta64 --grm ${TRAIT}_gcta_grm \\
           --make-bK-sparse ${params.sparse_cut} \\
           --out ${TRAIT}_sparse_grm \\
           --thread-num ${task.cpus}
    gcta64 --grm ${TRAIT}_gcta_grm \\
           --pca 1 \\
           --out ${TRAIT}_sparse_grm \\
           --thread-num ${task.cpus}\
    gcta64 --mlma-loco \\
           --grm ${TRAIT}_sparse_grm \\
           --bfile ${TRAIT} \\
           --qcovar ${TRAIT}_sparse_grm.eigenvec \\
           --out ${TRAIT}_lmm-exact_pca \\
           --pheno ${traits} \\
           --maf ${params.maf} \\
           --thread-num ${task.cpus}

    gcta64 --grm ${TRAIT}_gcta_grm_inbred \\
           --make-bK-sparse ${params.sparse_cut} \\
           --out ${TRAIT}_sparse_grm_inbred \\
           --thread-num ${task.cpus}
    gcta64 --grm ${TRAIT}_gcta_grm_inbred \\
           --pca 1 \\
           --out ${TRAIT}_sparse_grm_inbred \\
           --thread-num ${task.cpus}
    gcta64 --fastGWA-lmm-exact \\
           --grm-sparse ${TRAIT}_sparse_grm_inbred \\
           --bfile ${TRAIT} \\
           --qcovar ${TRAIT}_sparse_grm_inbred.eigenvec \\
           --out ${TRAIT}_lmm-exact_inbred_pca \\
           --pheno ${traits} \\
           --maf ${params.maf} \\
           --thread-num ${task.cpus}
    """
}

process gcta_lmm_exact_mapping_nopca {

    label "gcta_lmm_exact_mapping_nopca"
    errorStrategy 'retry'
    time { 20.minute * task.attempt }
    cpus = { 4 * task.attempt }
    memory = { 16.GB * task.attempt }

    publishDir "${params.out}/INBRED/Mapping/Raw", pattern: "*fastGWA", overwrite: true
    publishDir "${params.out}/LOCO/Mapping/Raw", pattern: "*.mlma", overwrite: true

    input:
    tuple val(TRAIT), file(traits), file(bed), file(bim), file(fam), file(map), \
    file(nosex), file(ped), file(grm_bin), file(grm_id), file(grm_nbin), \
    file(grm_bin_inbred), file(grm_id_inbred)

    output:
    //tuple val(TRAIT), file("${TRAIT}_lmm-exact_inbred.fastGWA"), file("${TRAIT}_lmm-exact.loco.mlma")
    tuple val(TRAIT), file("${TRAIT}_lmm-exact_inbred.fastGWA"), file("${TRAIT}_lmm-exact.loco.mlma")

    """
    gcta64 --grm ${TRAIT}_gcta_grm \\
           --make-bK-sparse ${params.sparse_cut} \\
           --out ${TRAIT}_sparse_grm \\
           --thread-num ${task.cpus}
    gcta64 --mlma-loco \\
           --grm ${TRAIT}_sparse_grm \\
           --bfile ${TRAIT} \\
           --out ${TRAIT}_lmm-exact \\
           --pheno ${traits} \\
           --maf ${params.maf} \\
           --thread-num ${task.cpus}

    gcta64 --grm ${TRAIT}_gcta_grm_inbred \\
           --make-bK-sparse ${params.sparse_cut} \\
           --out ${TRAIT}_sparse_grm_inbred \\
           --thread-num ${task.cpus}
    gcta64 --fastGWA-lmm-exact \\
           --grm-sparse ${TRAIT}_sparse_grm_inbred \\
           --bfile ${TRAIT} \\
           --out ${TRAIT}_lmm-exact_inbred \\
           --pheno ${traits} \\
           --maf ${params.maf} \\
           --thread-num ${task.cpus}
    """
}


process gcta_intervals_maps {

    label "gcta_intervals_maps"
    errorStrategy 'retry'
    time { 20.minute * task.attempt }
    cpus = { 4 * task.attempt }
    memory = { 16.GB * task.attempt }

    publishDir "${params.out}/INBRED/Mapping/Processed", mode: 'copy', pattern: "*_inbred.tsv"
    publishDir "${params.out}/LOCO/Mapping/Processed", mode: 'copy', pattern: "*_loco.tsv" 

    input:
        tuple val(TRAIT), file(pheno), file(tests), file(geno), val(sig_thresh), \
              val(qtl_grouping_size), val(qtl_ci_size), file(lmmexact_inbred), file(lmmexact_loco), \
              file(find_aggregate_intervals_maps)

    output:
        // tuple file(geno), file(pheno), val(TRAIT), file(tests), file("*AGGREGATE_mapping_inbred.tsv"), file("*AGGREGATE_mapping_loco.tsv"), emit: maps_to_plot
        path "*AGGREGATE_qtl_region_inbred.tsv", emit: qtl_peaks_inbred
        path  "*AGGREGATE_qtl_region_loco.tsv", emit: qtl_peaks_loco
        tuple val(TRAIT), file("*AGGREGATE_mapping_inbred.tsv"), file("*AGGREGATE_mapping_loco.tsv"), emit: for_html
        tuple val(TRAIT), file(geno), file(pheno), file(tests), file("*AGGREGATE_mapping_inbred.tsv"), file("*AGGREGATE_mapping_loco.tsv"), emit: maps_to_plot
        // path  "*AGGREGATE_qtl_region_inbred.tsv", emit: qtl_peaks_inbred, optional: true
        // tuple  val(TRAIT), file("*AGGREGATE_mapping_inbred.tsv"), emit: for_html_inbred, optional: true
        // tuple  val(TRAIT), file("*AGGREGATE_mapping_loco.tsv"), emit: for_html_loco, optional: true

    """
    Rscript --vanilla ${find_aggregate_intervals_maps} ${geno} ${pheno} ${lmmexact_inbred} \\
                      ${tests} ${qtl_grouping_size} ${qtl_ci_size} ${sig_thresh} \\
                      ${TRAIT}_AGGREGATE inbred
    Rscript --vanilla ${find_aggregate_intervals_maps} ${geno} ${pheno} ${lmmexact_loco} \\
                      ${tests} ${qtl_grouping_size} ${qtl_ci_size} ${sig_thresh} \\
                      ${TRAIT}_AGGREGATE loco
    """
}

