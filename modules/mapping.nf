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

    label "ml"

    input:
        tuple file(strains), val(TRAIT), file(traits), file(vcf), file(index), file(num_chroms)

    output:
        tuple val(TRAIT), file("plink_formated_trats.tsv"), file("${TRAIT}.bed"), file("${TRAIT}.bim"), file("${TRAIT}.fam"), file("${TRAIT}.map"), file("${TRAIT}.nosex"), file("${TRAIT}.ped"), file("${TRAIT}.log")

    """
    bcftools annotate --rename-chrs ${num_chroms} ${vcf} |\\
    bcftools view -S ${strains} -Ou |\\
    bcftools filter -i N_MISSING=0 -Oz --threads 5 -o renamed_chroms.vcf.gz
    tabix -p vcf renamed_chroms.vcf.gz
    plink --vcf renamed_chroms.vcf.gz \\
          --threads 5 \\
          --snps-only \\
          --biallelic-only \\
          --maf ${params.maf} \\
          --set-missing-var-ids @:# \\
          --indep-pairwise 50 10 0.8 \\
          --geno \\
          --allow-extra-chr
    tail -n +2 ${traits} | awk 'BEGIN {OFS="\\t"}; {print \$1, \$1, \$2}' > plink_formated_trats.tsv
    plink --vcf renamed_chroms.vcf.gz \\
          --threads 5 \\
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
          --pheno plink_formated_trats.tsv
    """
}

process gcta_grm {

    label "lg"

    input:
        tuple val(TRAIT), file(traits), file(bed), file(bim), file(fam), file(map), file(nosex), \
              file(ped), file(log), val(algorithm)

    output:
        tuple val(TRAIT), file(traits), file(bed), file(bim), file(fam), file(map), file(nosex), \
              file(ped), file(log), file("${TRAIT}_gcta_grm_${algorithm}.grm.bin"), \
              file("${TRAIT}_gcta_grm_${algorithm}.grm.id"), \
              file("${TRAIT}_gcta_grm_${algorithm}.grm.N.bin"), val(algorithm)

    when:
        params.mapping

    """
    if [[ ${algorithm} == "inbred" ]];
    then
        option="--make-grm-inbred"
    else
        options="--make-grm"
    fi

    gcta64 --bfile ${TRAIT} \\
           --autosome \\
           --maf ${params.maf} \\
           \${option} \\
           --out ${TRAIT}_gcta_grm_${algorithm} \\
           --thread-num 5
    """
}


process gcta_lmm_exact_mapping {

    label "lg"

    publishDir "${params.out}/INBRED/Mapping/Raw", pattern: "*inbred_pca.fastGWA", overwrite: true
    publishDir "${params.out}/LOCO/Mapping/Raw", pattern: "*loco_pca.mlma", overwrite: true

    input:
    tuple val(TRAIT), file(traits), file(bed), file(bim), file(fam), file(map), \
          file(nosex), file(ped), file(log), file(grm_bin), file(grm_id), file(grm_nbin),
          val(algorithm)

    output:
    tuple val(TRAIT), file("${TRAIT}_lmm-exact_${algorithm}_pca.*", arity=1)


    """
    if [[ ${algorithm} == "inbred" ]];
    then
        options="--fastqGWA-lmm-exact --grm-sparse"
    else
        options="--mlma-loco --grm"
    fi
    
    gcta64 --grm ${TRAIT}_gcta_grm_${algorithm} \\
           --make-bK-sparse ${params.sparse_cut} \\
           --out ${TRAIT}_sparse_grm_${algorithm} \\
           --thread-num ${task.cpus}
    gcta64 --grm ${TRAIT}_gcta_grm_${algorithm} \\
           --pca 1 \\
           --out ${TRAIT}_sparse_grm_${algorithm} \\
           --thread-num ${task.cpus}
    gcta64 \${options} ${TRAIT}_sparse_grm \\
           --bfile ${TRAIT} \\
           --out ${TRAIT}_lmm-exact_${algorithm} \\
           --pheno ${traits} \\
           --maf ${params.maf} \\
           --thread-num ${task.cpus}
    gcta64 \${options} ${TRAIT}_sparse_grm \\
           --bfile ${TRAIT} \\
           --qcovar ${TRAIT}_sparse_grm_${algorithm}.eigenvec \\
           --out ${TRAIT}_lmm-exact_${algorithm}_pca \\
           --pheno ${traits} \\
           --maf ${params.maf} \\
           --thread-num ${task.cpus}
    """
}

process gcta_lmm_exact_mapping_nopca {

    // machineType 'n1-highmem-4'
    label "lg"

    publishDir "${params.out}/INBRED/Mapping/Raw", pattern: "*fastGWA", overwrite: true
    publishDir "${params.out}/LOCO/Mapping/Raw", pattern: "*loco.mlma", overwrite: true

    input:
    tuple val(TRAIT), file(traits), file(bed), file(bim), file(fam), file(map), \
    file(nosex), file(ped), file(log), file(grm_bin), file(grm_id), file(grm_nbin), \
    val(algorithm)

    output:
    tuple val(TRAIT), file("${TRAIT}_lmm-exact_${algorithm}.*", arity=1)


    """
    if [[ ${algorithm} == "inbred" ]];
    then
        options="--fastqGWA-lmm-exact --grm-sparse"
    else
        options="--mlma-loco --grm"
    fi
    
    gcta64 --grm ${TRAIT}_gcta_grm_${algorithm} \\
           --make-bK-sparse ${params.sparse_cut} \\
           --out ${TRAIT}_sparse_grm_${algorithm} \\
           --thread-num ${task.cpus}
    gcta64 --grm ${TRAIT}_gcta_grm_${algorithm} \\
           --pca 1 \\
           --out ${TRAIT}_sparse_grm_${algorithm} \\
           --thread-num ${task.cpus}
    gcta64 \${options} ${TRAIT}_sparse_grm \\
           --bfile ${TRAIT} \\
           --out ${TRAIT}_lmm-exact_${algorithm} \\
           --pheno ${traits} \\
           --maf ${params.maf} \\
           --thread-num ${task.cpus}
    gcta64 \${options} ${TRAIT}_sparse_grm \\
           --bfile ${TRAIT} \\
           --qcovar ${TRAIT}_sparse_grm_${algorithm}.eigenvec \\
           --out ${TRAIT}_lmm-exact_${algorithm}_pca \\
           --pheno ${traits} \\
           --maf ${params.maf} \\
           --thread-num ${task.cpus}
    """
}


process gcta_intervals_maps {

    label "lg"

    publishDir "${params.out}/INBRED/Mapping/Processed", mode: 'copy', pattern: "*_inbred.tsv"
    publishDir "${params.out}/LOCO/Mapping/Processed", mode: 'copy', pattern: "*_loco.tsv" 

    input:
        tuple val(TRAIT), file(pheno), file(tests), file(geno), val(P3D), val(sig_thresh), \
              val(qtl_grouping_size), val(qtl_ci_size), file(lmmexact), \
              file(find_aggregate_intervals_maps), val(algorithm)

    output:
        tuple file(geno), file(pheno), val(TRAIT), file(tests), file("*AGGREGATE_mapping_${algorithm}.tsv"), emit: maps_to_plot
        path "*AGGREGATE_qtl_region${algorithm}.tsv", emit: qtl_peaks
        tuple file("*AGGREGATE_mapping_${algorithm}.tsv"), val(TRAIT), val(algorithm), emit: for_html

    """
    Rscript --vanilla ${find_aggregate_intervals_maps} ${geno} ${pheno} ${lmmexact} \\
                      ${tests} ${qtl_grouping_size} ${qtl_ci_size} ${sig_thresh} \\
                      ${TRAIT}_AGGREGATE ${algorithm}
    """
}

