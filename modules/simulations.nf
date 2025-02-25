/*
======================================
~ > *                            * < ~
~ ~ > *                        * < ~ ~
~ ~ ~ > *  GWAS SIMULATIONS  * < ~ ~ ~
~ ~ > *                        * < ~ ~
~ > *                            * < ~
======================================
*/


process prepare_simulation_files {
    
    label "ml"

    input:
        tuple val(strain_set), val(strains), file(vcf), file(index), file(num_chroms), val(MAF)

    output:
        tuple val(strain_set), val(strains), file("TO_SIMS.bed"), file("TO_SIMS.bim"), file("TO_SIMS.fam"), file("TO_SIMS.map"), file("TO_SIMS.nosex"), file("TO_SIMS.ped"), file("TO_SIMS.log"), file("${strain_set}_${MAF}_Genotype_Matrix.tsv"), val(MAF), emit: sim_geno
        tuple val(strain_set), val(strains), val(MAF), file("renamed_chroms.vcf.gz"), file("renamed_chroms.vcf.gz.tbi"), emit: renamed_chrom_vcf_to_ld


    """
    export > node_name.txt
    bcftools view --force-samples -s ${strains} ${vcf} | \\
    bcftools annotate --rename-chrs ${num_chroms} ${vcf} |\\
    bcftools filter -i N_MISSING=0 -Oz -o renamed_chroms.vcf.gz
    tabix -p vcf renamed_chroms.vcf.gz
    plink --vcf renamed_chroms.vcf.gz \\
    --snps-only \\
    --biallelic-only \\
    --maf ${MAF} \\
    --set-missing-var-ids @:# \\
    --indep-pairwise 50 10 0.8 \\
    --geno \\
    --not-chr MtDNA \\
    --allow-extra-chr
    plink --vcf renamed_chroms.vcf.gz \\
    --make-bed \\
    --snps-only \\
    --biallelic-only \\
    --maf ${MAF} \\
    --set-missing-var-ids @:# \\
    --extract plink.prune.in \\
    --geno \\
    --recode \\
    --out TO_SIMS \\
    --allow-extra-chr
    awk -F":" '\$1=\$1' OFS="\\t" plink.prune.in | \\
    sort -k1,1d -k2,2n > markers.txt
    bcftools query -l renamed_chroms.vcf.gz |\\
    sort > sorted_samples.txt
    bcftools view -v snps \\
    -S sorted_samples.txt \\
    -R markers.txt \\
    renamed_chroms.vcf.gz |\\
    bcftools query --print-header -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' |\\
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
    sed 's/.\\/./NA/g' > ${strain_set}_${MAF}_Genotype_Matrix.tsv
    """
}


/*
------------ Decomposition per chromosome
*/

process chrom_eigen_variants_sims {

    tag { CHROM }

    label "ml"

    input:
        tuple val(CHROM), \
        val(strain_set), \
        val(strains), \
        file(bed), \
        file(bim), \
        file(fam), \
        file(map), \
        file(sex), \
        file(ped), \
        file(log), \
        file(geno), \
        val(MAF), \
        file(get_genomatrix_eigen)

    output:
        tuple val(strain_set), val(strains), val(MAF), file(bed), file(bim), file(fam), file(map), file(sex), file(ped), file(log), file(geno), emit: sim_geno_meta
        tuple val(strain_set), val(strains), val(MAF), file("${CHROM}_${strain_set}_${MAF}_independent_snvs.csv"), emit: sim_geno_eigen_join


    """
    export > node_name.txt
    cat ${geno} |\\
    awk -v chrom="${CHROM}" '{if(\$1 == "CHROM" || \$1 == chrom) print}' > ${CHROM}_gm.tsv
    Rscript --vanilla ${get_genomatrix_eigen} ${CHROM}_gm.tsv ${CHROM}
    mv ${CHROM}_independent_snvs.csv ${CHROM}_${strain_set}_${MAF}_independent_snvs.csv
    """

}

/*
------------ Sum independent tests for all chromosomes
*/

process collect_eigen_variants_sims {

    executor 'local'
    container null

    publishDir "${params.out}/Genotype_Matrix", mode: 'copy'

    input:
        tuple val(strain_set), val(strains), val(MAF), file(tests), file(bed), file(bim), file(fam), file(map), file(sex), file(ped), file(log), file(geno)


    output:
        tuple val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(sex), file(ped), file(log), file(geno), val(MAF), file("${strain_set}_${MAF}_total_independent_tests.txt")

    """
    cat *independent_snvs.csv |\\
    grep -v inde |\\
    awk '{s+=\$1}END{print s}' > ${strain_set}_${MAF}_total_independent_tests.txt
    """

}


process simulate_effects_loc {

    tag {NQTL}

    label "md"

    input:
        tuple val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(MAF), file(n_indep_tests), val(NQTL), file(qtl_loc_bed), val(effect_range), val(SIMREP), file(create_causal_qtls)

    output:
        tuple val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(MAF), file(n_indep_tests), val(NQTL), val(SIMREP), val(effect_range), file("causal.variants.sim.${NQTL}.${SIMREP}.txt")

    """
    export > node_name.txt
    Rscript --vanilla ${create_causal_qtls} ${bim} ${NQTL} ${effect_range} ${qtl_loc_bed}
    mv causal.variants.sim.${NQTL}.txt causal.variants.sim.${NQTL}.${SIMREP}.txt
    """
}


process simulate_effects_genome {

    tag {NQTL}
    
    label "xs"

    input:
        tuple val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(MAF), file(n_indep_tests), val(NQTL), val(effect_range), val(SIMREP), file(create_causal_qtls)

    output:
        tuple val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(MAF), file(n_indep_tests), val(NQTL), val(SIMREP), val(effect_range), file("causal.variants.sim.${NQTL}.${SIMREP}.txt")


    """
    export > node_name.txt
    Rscript --vanilla ${create_causal_qtls} ${bim} ${NQTL} ${effect_range}
    mv causal.variants.sim.${NQTL}.txt causal.variants.sim.${NQTL}.${SIMREP}.txt
    """
}


process simulate_map_phenotypes {

    tag {"${NQTL} - ${SIMREP} - ${H2} - ${MAF}"}

    label "md"

    publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Mappings", pattern: "*fastGWA", overwrite: true
    publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Mappings", pattern: "*loco.mlma", overwrite: true
    publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Phenotypes", pattern: "*.phen", overwrite: true
    publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Phenotypes", pattern: "*.par", overwrite: true

    input:
        tuple val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(MAF), file(n_indep_tests), val(NQTL), val(SIMREP), val(effect_range), file(loci), val(H2), file(check_vp)

    output:
        tuple file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.bed"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.bim"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.fam"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.map"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.nosex"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.ped"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.log"), val(NQTL), val(SIMREP), file(loci), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.par"), emit: sim_phen_output
        tuple file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred.fastGWA"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact.loco.mlma"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact.log"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred.log"), emit: sim_GCTA_mapping_results
        path "${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred.fastGWA", emit: lmm_exact_inbred_analyze_sims
        path "${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred_pca.fastGWA", emit: lmm_exact_inbred_pca_analyze_sims
        path "${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact.loco.mlma", emit: lmm_exact_loco_analyze_sims
        path "${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_pca.loco.mlma", emit: lmm_exact_loco_pca_analyze_sims
        path "${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen", emit: simphen_analyze_sims
        path "${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.par", emit: simgen_analyze_sims
        tuple val(strain_set), val(strains), val(NQTL), val(SIMREP), val(H2), file(loci), file(gm), val(effect_range), file(n_indep_tests), val(MAF), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred_pca.fastGWA"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact.loco.mlma"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.par"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen"), emit: gcta_intervals

    """
    export > node_name.txt
    gcta64 --bfile TO_SIMS \\
         --simu-qt \\
         --simu-causal-loci ${loci} \\
         --simu-hsq ${H2} \\
         --simu-rep 1 \\
         --autosome-num 6 \\
         --thread-num ${task.cpus} \\
         --out ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims
    echo "a"
    plink --bfile TO_SIMS \\
        --make-bed \\
        --snps-only \\
        --biallelic-only \\
        --maf ${MAF} \\
        --set-missing-var-ids @:# \\
        --geno \\
        --recode \\
        --out TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set} \\
        --allow-extra-chr \\
        --allow-no-sex \\
        --pheno ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen
    echo "b"
    gcta64 --bfile TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set} \\
            --autosome --maf ${MAF} --make-grm \\
            --out TO_SIMS_${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_gcta_grm \\
            --thread-num ${task.cpus}
    echo "c"
    gcta64 --bfile TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set} \\
            --autosome --maf ${MAF} --make-grm-inbred \\
            --out TO_SIMS_${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_gcta_grm_inbred \\
            --thread-num ${task.cpus}
    echo "d"
    gcta64 --grm TO_SIMS_${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_gcta_grm_inbred \\
            --pheno ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen \\
            --reml --out check_vp \\
            --thread-num ${task.cpus}
    
    echo "e"
    python ${check_vp} --check_vp check_vp.hsq --simulated_phenos ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen 

    if [[ -f "new_phenos.temp" ]]
    then
        rm ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen
        mv new_phenos.temp ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen
    fi

    echo "f"
    gcta64 --grm TO_SIMS_${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_gcta_grm \\
           --make-bK-sparse ${params.sparse_cut} \\
           --out ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sparse_grm \\
           --thread-num ${task.cpus}
    echo "g"
    gcta64 --grm TO_SIMS_${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_gcta_grm \\
           --pca 1 \\
           --out ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sparse_grm \\
           --thread-num ${task.cpus}
    echo "h"
    gcta64 --mlma-loco \\
           --bfile TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set} \\
           --grm ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sparse_grm \\
           --out ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact \\
           --pheno ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen \\
           --maf ${MAF} \\
           --thread-num ${task.cpus}
    echo "i"
    gcta64 --mlma-loco \\
           --bfile TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set} \\
           --grm ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sparse_grm \\
           --qcovar ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sparse_grm.eigenvec \\
           --out ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_pca \\
           --pheno ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen \\
           --maf ${MAF} \\
           --thread-num ${task.cpus}
    echo "j"

    gcta64 --grm TO_SIMS_${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_gcta_grm_inbred \\
          --make-bK-sparse ${params.sparse_cut} \\
          --out ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sparse_grm_inbred \\
          --thread-num ${task.cpus}
    echo "k"
    gcta64 --grm TO_SIMS_${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_gcta_grm_inbred \\
          --pca 1 \\
          --out ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sparse_grm_inbred \\
          --thread-num ${task.cpus}
    echo "l"
    gcta64 --fastGWA-lmm-exact \\
          --grm-sparse ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sparse_grm_inbred \\
          --bfile TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set} \\
          --out ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred \\
          --pheno ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen \\
          --maf ${MAF} \\
          --thread-num ${task.cpus}
    echo "m"
    gcta64 --fastGWA-lmm-exact \\
          --grm-sparse ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sparse_grm_inbred \\
          --qcovar ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sparse_grm_inbred.eigenvec \\
          --bfile TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set} \\
          --out ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred_pca \\
          --pheno ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen \\
          --maf ${MAF} \\
          --thread-num ${task.cpus}
    """
}


process get_gcta_intervals {

    tag {"${NQTL} - ${SIMREP} - ${H2}"}

    publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Mappings", mode: 'copy', pattern: "*processed_LMM-EXACT-INBRED_PCA_mapping.tsv"
    publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Mappings", mode: 'copy', pattern: "*processed_LMM-EXACT-LOCO_PCA_mapping.tsv"
    publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Mappings", mode: 'copy', pattern: "*qtl_region.tsv"

    memory '50GB'
    time '20min'
    cpus 1
    errorStrategy 'ignore'

    input:
    tuple val(strain_set), val(strains), val(NQTL), val(SIMREP), val(H2), file(loci), file(gm), val(effect_range), file(n_indep_tests), val(MAF), file(lmmexact_inbred), file(lmmexact_loco), \
    file(var_effects), file(phenotypes) ,val(THRESHOLD), val(QTL_GROUP_SIZE), val(QTL_CI_SIZE), file(aggregate_mappings), file(find_aggregate_intervals), file(find_gcta_intervals), file(find_gcta_intervals_loco)

    output:
    tuple val(strain_set), val(strains), val(NQTL), val(SIMREP), val(H2), val(MAF), val(effect_range), path(var_effects), path(phenotypes), path(gm), file("*processed_LMM-EXACT-INBRED_PCA_mapping.tsv"), val("LMM-EXACT-INBRED_PCA"), emit: assess_data_inbred_pca 
    tuple val(strain_set), val(strains), val(NQTL), val(SIMREP), val(H2), val(MAF), val(effect_range), path(var_effects), path(phenotypes), path(gm), file("*processed_LMM-EXACT-LOCO_PCA_mapping.tsv"), val("LMM-EXACT-LOCO_PCA"), emit: assess_data_loco_pca 

    script:
    """
    export > node_name.txt
    Rscript --vanilla ${find_gcta_intervals} ${gm} ${phenotypes} ${lmmexact_inbred} ${n_indep_tests} ${NQTL} ${SIMREP} ${QTL_GROUP_SIZE} ${QTL_CI_SIZE} ${H2} ${params.maf} ${THRESHOLD} ${strain_set} ${MAF} ${effect_range} LMM-EXACT-INBRED_PCA
    Rscript --vanilla ${find_gcta_intervals_loco} ${gm} ${phenotypes} ${lmmexact_loco} ${n_indep_tests} ${NQTL} ${SIMREP} ${QTL_GROUP_SIZE} ${QTL_CI_SIZE} ${H2} ${params.maf} ${THRESHOLD} ${strain_set} ${MAF} ${effect_range} LMM-EXACT-LOCO_PCA
    """
}

process assess_sims_INBRED {

    container 'mckeowr1/asess_sims:1.1'
 
    executor 'local'   
    publishDir "${params.out}/scored_sims", mode: 'copy', pattern: "*_mapping.tsv"
    
    input:
        tuple val(strain_set), val(strains), val(NQTL), val(SIMREP), val(H2), val(MAF), val(effect_range), path(var_effects), path(phenotypes), path(gm), path(mapping_processed), val(algorithm_id), path(R_assess_sims)
    output: 
        file("*_mapping.tsv")        
    
    script:
    """
    export > node_name.txt
    Rscript --vanilla ${R_assess_sims} ${mapping_processed} ${gm} ${var_effects} ${phenotypes} ${NQTL} ${SIMREP} ${H2} ${MAF} ${effect_range} ${strain_set} ${algorithm_id}
    """
    
}
process assess_sims_LOCO {

    container 'mckeowr1/asess_sims:1.1'
    executor = 'local'

    publishDir "${params.out}/scored_sims", mode: 'copy', pattern: "*_mapping.tsv"
    
    input:
        tuple val(strain_set), val(strains), val(NQTL), val(SIMREP), val(H2), val(MAF), val(effect_range), path(var_effects), path(phenotypes), path(gm), path(mapping_processed), val(algorithm_id), path(R_assess_sims)
    output: 
        file("*_mapping.tsv")        
    
    script:
    """
    export > node_name.txt
    Rscript --vanilla ${R_assess_sims} ${mapping_processed} ${gm} ${var_effects} ${phenotypes} ${NQTL} ${SIMREP} ${H2} ${MAF} ${effect_range} ${strain_set} ${algorithm_id}
    """
    
}
/*
======================================
~ > *                            * < ~
~ ~ > *                        * < ~ ~
~ ~ ~ > *         LD         * < ~ ~ ~
~ ~ > *                        * < ~ ~
~ > *                            * < ~
======================================
    gcta_qtl_to_ld
        .combine(emmma_qtl_to_ld, by: [0,1,2,3,4,5,6])
        .combine(simulated_phenotypes, by: [0,1,2,3,4,5,6])
        .combine(renamed_chrom_vcf_to_ld, by: [0,1,2])
        .into{peaks_to_ld}
    process extract_qtl_ld {
    cpus 1
    input:
    set val(strain_set), val(strains), val(MAF), val(NQTL), val(SIMREP), val(H2), val(effect_range), file(lmmexact_peaks), file(lmmexact_inbred_peaks), file(emma_peaks), file(loci), file(phenotypes), file(vcf), file(vcfindex) from peaks_to_ld
    output:
    set val(strain_set), val(strains), val(MAF), val(NQTL), val(SIMREP), val(H2), val(effect_range), file(lmmexact_peaks), file(lmmexact_inbred_peaks), file(emma_peaks), file(loci), file(phenotypes), file(vcf), file(vcfindex), file("*ROI_Genotype_Matrix.tsv"), file("*LD.tsv") into LD_files_to_finemap
    """
    cat *_qtl_region.tsv |\\
    awk '\$0 !~ "CHROM" {print}' > QTL_peaks.tsv
    filename='QTL_peaks.tsv'
    echo Start
    while read p; do
        chromosome=`echo \$p | cut -f1 -d' '`
        trait=`echo \$p | cut -f3 -d' '`
        start_pos=`echo \$p | cut -f4 -d' '`
        peak_pos=`echo \$p | cut -f5 -d' '`
        end_pos=`echo \$p | cut -f6 -d' '`
        map_algo=`echo \$p | cut -f8 -d' '`
    bcftools view --regions \$chromosome:\$start_pos-\$end_pos -s ${strains} ${vcf} |\\
    bcftools filter -i N_MISSING=0 |\\
    awk '\$0 !~ "#" {print \$1":"\$2}' > \$trait.\$chromosome.\$start_pos.\$end_pos.txt
    bcftools view --regions \$chromosome:\$start_pos-\$end_pos -s ${strains} ${vcf} |\\
    bcftools filter -i N_MISSING=0 -Oz -o finemap.vcf.gz
    plink --vcf finemap.vcf.gz \\
        --snps-only \\
        --maf ${MAF} \\
        --biallelic-only \\
        --allow-extra-chr \\
        --set-missing-var-ids @:# \\
        --geno \\
        --make-bed \\
        --recode vcf-iid bgz \\
        --extract \$trait.\$chromosome.\$start_pos.\$end_pos.txt \\
        --out \$trait.\$chromosome.\$start_pos.\$end_pos
    nsnps=`wc -l \$trait.\$chromosome.\$start_pos.\$end_pos.txt | cut -f1 -d' '`
    plink --r2 with-freqs \\
        --allow-extra-chr \\
        --snps-only \\
        --ld-window-r2 0 \\
        --ld-snp \$chromosome:\$peak_pos \\
        --ld-window \$nsnps \\
        --ld-window-kb 6000 \\
        --chr \$chromosome \\
        --out \$trait.\$chromosome:\$start_pos-\$end_pos.QTL \\
        --set-missing-var-ids @:# \\
        --vcf \$trait.\$chromosome.\$start_pos.\$end_pos.vcf.gz
    sed 's/  \\t/g' \$trait.\$chromosome:\$start_pos-\$end_pos.QTL.ld |\\
    cut -f2-10 |\\
    sed 's/^23/X/g' | sed 's/\\t23\\t/\\tX\\t/g' > \$trait.\$chromosome.\$start_pos.\$end_pos.\$map_algo.LD.tsv
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
        sed 's/^23/X/g' > \$trait.\$chromosome:\$start_pos-\$end_pos.\$map_algo.ROI_Genotype_Matrix.tsv
    done < \$filename
    """
}
*/

/*
------------ Run fine mapping
process sim_fine_maps {
    echo 'true'
    input:
        set val(strain_set), val(strains), val(MAF), val(NQTL), val(SIMREP), val(H2), val(effect_range), file(lmmexact_peaks), file(lmmexact_inbred_peaks), file(emma_peaks), file(loci), file(phenotypes), file(vcf), file(vcfindex), file(roi_geno_matrix), file(roi_ld) from LD_files_to_finemap
    output:
    """
        for i in *ROI_Genotype_Matrix.tsv;
        do
            echo \$i
        done
    """
}
*/
