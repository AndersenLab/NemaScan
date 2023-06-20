

process prepare_repeated_simulation_files {
    container 'andersenlab/nemascan:20220407173056db3227'
    cpus 5
    memory 30.GB
    time '30m'

    input:
        tuple val(sp), val(strain_set), val(strains), file(vcf), file(index), file(plink_dir), file(num_chroms), val(MAF)

    output:
        tuple val(sp), val(strain_set), val(strains), file("TO_SIMS.bed"), file("TO_SIMS.bim"), file("TO_SIMS.fam"), file("TO_SIMS.map"), file("TO_SIMS.nosex"), file("TO_SIMS.ped"), file("TO_SIMS.log"), file("${sp}_${strain_set}_${MAF}_Genotype_Matrix.tsv"), val(MAF), emit: sim_geno
        tuple val(sp), val(strain_set), val(strains), val(MAF), file("renamed_chroms.vcf.gz"), file("renamed_chroms.vcf.gz.tbi"), emit: renamed_chrom_vcf_to_ld


    """
    bcftools annotate --rename-chrs ${num_chroms} ${vcf} |\\
    bcftools view -s `echo ${strains} | tr -d '\\n'` |\\
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
    sed 's/.\\/./NA/g' > ${sp}_${strain_set}_${MAF}_Genotype_Matrix.tsv
    """
}
// Temp version that just creates expected output files
process prepare_repeated_simulation_files_temp{
    executor 'local'
    input:
        tuple val(sp), val(strain_set), val(strains), file(vcf), file(index), file(plink_dir), file(num_chroms), val(MAF)

    output:
        tuple val(sp), val(strain_set), val(strains), file("TO_SIMS.bed"), file("TO_SIMS.bim"), file("TO_SIMS.fam"), file("TO_SIMS.map"), file("TO_SIMS.nosex"), file("TO_SIMS.ped"), file("TO_SIMS.log"), file("${sp}_${strain_set}_${MAF}_Genotype_Matrix.tsv"), val(MAF), emit: sim_geno
        //tuple val(sp), val(strain_set), val(strains), val(MAF), file("renamed_chroms.vcf.gz"), file("renamed_chroms.vcf.gz.tbi"), emit: renamed_chrom_vcf_to_ld


    """
    cp plink_files/TO_SIMS* .
    cp plink_files/*.tsv .
    """
    
}
/*
------------ Decomposition per chromosome
*/

process chrom_eigen_variants_sims_repeated  {

    tag { CHROM }

    cpus 6
    time '5m'
    memory 5.GB
    container = 'andersenlab/nemascan:20220407173056db3227'

    memory params.eigen_mem

    input:
        tuple val(CHROM), val(sp), val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(sex), file(ped), file(log), file(geno), val(MAF), file(get_genomatrix_eigen)

    output:
        tuple val(sp), val(strain_set), val(strains), val(MAF), file(bed), file(bim), file(fam), file(map), file(sex), file(ped), file(log), file(geno), emit: sim_geno_meta
        tuple val(sp), val(strain_set), val(strains), val(MAF), file("${CHROM}_${sp}_${strain_set}_${MAF}_independent_snvs.csv"), emit: sim_geno_eigen_join


    """
        cat ${geno} |\\
        awk -v chrom="${CHROM}" '{if(\$1 == "CHROM" || \$1 == chrom) print}' > ${CHROM}_gm.tsv
        Rscript --vanilla ${get_genomatrix_eigen} ${CHROM}_gm.tsv ${CHROM}
        mv ${CHROM}_independent_snvs.csv ${CHROM}_${sp}_${strain_set}_${MAF}_independent_snvs.csv
    """

}
// Temp version that just creates expected output files

process chrom_eigen_variants_sims_repeated_temp  {

    tag { CHROM }
    executor 'local'


    input:
        tuple val(CHROM), val(sp), val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(sex), file(ped), file(log), file(geno), val(MAF), file(get_genomatrix_eigen)

    output:
        tuple val(sp), val(strain_set), val(strains), val(MAF), file(bed), file(bim), file(fam), file(map), file(sex), file(ped), file(log), file(geno), emit: sim_geno_meta
        tuple val(sp), val(strain_set), val(strains), val(MAF), file("${CHROM}_${sp}_${strain_set}_${MAF}_independent_snvs.csv"), emit: sim_geno_eigen_join


    """
        cat > ${CHROM}_${sp}_${strain_set}_${MAF}_independent_snvs.csv
    """

}
/*
------------ Sum independent tests for all chromosomes
*/

process collect_eigen_variants_sims_repeated {

    //executor 'local'

    publishDir "${params.out}/Genotype_Matrix", mode: 'copy'
    container = 'andersenlab/nemascan:20220407173056db3227'
    cpus 1
    time '5m'
    memory 5.GB



    input:
        tuple val(sp), val(strain_set), val(strains), val(MAF), file(tests), file(bed), file(bim), file(fam), file(map), file(sex), file(ped), file(log), file(geno)


    output:
        tuple val(sp), val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(sex), file(ped), file(log), file(geno), val(MAF), file("${sp}_${strain_set}_${MAF}_total_independent_tests.txt")

    """
        cat *independent_snvs.csv |\\
        grep -v inde |\\
        awk '{s+=\$1}END{print s}' > ${sp}_${strain_set}_${MAF}_total_independent_tests.txt
    """

}
// Temp version that just creates expected output files
process collect_eigen_variants_sims_repeated_temp {

    executor 'local'


    input:
        tuple val(sp), val(strain_set), val(strains), val(MAF), file(tests), file(bed), file(bim), file(fam), file(map), file(sex), file(ped), file(log), file(geno)


    output:
        tuple val(sp), val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(sex), file(ped), file(log), file(geno), val(MAF), file("${sp}_${strain_set}_${MAF}_total_independent_tests.txt")

    """
    cat > ${sp}_${strain_set}_${MAF}_total_independent_tests.txt
    """

}

/*  
------------  Simulate effects for orthogroup variants
*/

process simulate_orthogroup_effects {
    label 'causal_ogs'
    executor 'local'
    conda '/home/rjm6024/.conda/envs/vcf_stats_1.0'


    input:
        tuple val(sp), val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(MAF), file(n_indep_tests), val(og1), val(og2), val(og3), val(og4), val(og5), val(SIMREP), file(create_causal_qtls), file(master_snps_dir)

    output:
        tuple val(sp), val(strain_set), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(SIMREP), val(MAF), file(n_indep_tests), val(og1), val(og2), val(og3), val(og4), val(og5), file("${sp}_${strain_set}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${SIMREP}_causal_og_vars.txt"), emit: pheno_inputs
        //tuple val(sp), val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(MAF), file(n_indep_tests), val(og1), val(og2), val(og3), val(og4), val(og5), val(SIMREP), file(master_snps_dir), file("${sp}_${strain_set}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${SIMREP}_causal_og_vars.txt")


    """
        python ${create_causal_qtls} ${og1} ${og2} ${og3} ${og4} ${og5} ${bim} ${master_snps_dir} ${sp}
        cat causal_og_vars.txt > ${sp}_${strain_set}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${SIMREP}_causal_og_vars.txt
    """
}

process simulate_map_phenotypes {

    label 'sim_map_phenos'
    tag {"${SIMREP} - ${H2} - ${MAF}"}

    //errorStrategy 'retry'

    publishDir "${params.out}/Simulations/${og1}_${og2}_${og3}_${og4}_${og5}_${sp}/Mappings", pattern: "*fastGWA", overwrite: true
    publishDir "${params.out}/Simulations/${og1}_${og2}_${og3}_${og4}_${og5}_${sp}/Mappings", pattern: "*loco.mlma", overwrite: true
    publishDir "${params.out}/Simulations/${og1}_${og2}_${og3}_${og4}_${og5}_${sp}/Phenotypes", pattern: "*.phen", overwrite: true
    publishDir "${params.out}/Simulations/${og1}_${og2}_${og3}_${og4}_${og5}_${sp}/Phenotypes", pattern: "*.par", overwrite: true

    cpus 5
    time '20m'
    memory 10.GB



    input:
        tuple val(sp), val(strain_set), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(SIMREP), val(MAF), file(n_indep_tests), val(og1), val(og2), val(og3), val(og4), val(og5), file(loci), val(H2)

    output:
        tuple file("TO_SIMS_${SIMREP}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}.bed"), file("TO_SIMS_${SIMREP}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}.bim"), file("TO_SIMS_${SIMREP}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}.fam"), file("TO_SIMS_${SIMREP}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}.map"), file("TO_SIMS_${SIMREP}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}.nosex"), file("TO_SIMS_${SIMREP}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}.ped"), file("TO_SIMS_${SIMREP}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}.log"), val(SIMREP), file(loci), file("${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_sims.phen"), file("${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_sims.par"), emit: sim_phen_output
        tuple file("${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_lmm-exact_inbred.fastGWA"), file("${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_lmm-exact.loco.mlma"), file("${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_lmm-exact.log"), file("${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_lmm-exact_inbred.log"), emit: sim_GCTA_mapping_results
        
        path "${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_lmm-exact_inbred.fastGWA", emit: lmm_exact_inbred_analyze_sims
        path "${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_lmm-exact_inbred_pca.fastGWA", emit: lmm_exact_inbred_pca_analyze_sims
        path "${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_lmm-exact.loco.mlma", emit: lmm_exact_loco_analyze_sims
        path "${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_lmm-exact_pca.loco.mlma", emit: lmm_exact_loco_pca_analyze_sims
        path "${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_sims.phen", emit: simphen_analyze_sims
        path "${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_sims.par", emit: simgen_analyze_sims
        tuple val(sp), val(strain_set), val(SIMREP), val(H2), file(loci), file(gm), file(n_indep_tests), val(MAF),val(og1), val(og2), val(og3), val(og4), val(og5), file("${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_lmm-exact_inbred.fastGWA"), file("${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_lmm-exact_inbred_pca.fastGWA"),file("${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_lmm-exact.loco.mlma"), file("${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_lmm-exact_pca.loco.mlma"), file("${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_sims.phen"), emit: gcta_intervals

    """
    gcta64 --bfile TO_SIMS \\
         --simu-qt \\
         --simu-causal-loci ${loci} \\
         --simu-hsq ${H2} \\
         --simu-rep 1 \\
         --thread-num 5 \\
         --out ${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_sims
    plink --bfile TO_SIMS \\
        --make-bed \\
        --snps-only \\
        --biallelic-only \\
        --maf ${MAF} \\
        --set-missing-var-ids @:# \\
        --geno \\
        --recode \\
        --out TO_SIMS_${SIMREP}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set} \\
        --allow-extra-chr \\
        --pheno ${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_sims.phen
    gcta64 --bfile TO_SIMS_${SIMREP}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set} \\
            --autosome --maf ${MAF} --make-grm \\
            --out TO_SIMS_${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_gcta_grm \\
            --thread-num 5
    gcta64 --bfile TO_SIMS_${SIMREP}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set} \\
            --autosome --maf ${MAF} --make-grm-inbred \\
            --out TO_SIMS_${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_gcta_grm_inbred \\
            --thread-num 5
    gcta64 --grm TO_SIMS_${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_gcta_grm_inbred \\
            --pheno ${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_sims.phen \\
            --reml --out check_vp \\
            --thread-num 5
    vp=`grep Vp check_vp.hsq | head -1 | cut -f2`
    if (( \$(echo "0.00001 > \$vp" |bc -l) ));
      then
        awk '{print \$1, \$2, \$3*1000}' ${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_sims.phen > temp.phen;
        rm ${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_sims.phen
        mv temp.phen ${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_sims.phen
    fi

    gcta64 --grm TO_SIMS_${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_gcta_grm \\
           --make-bK-sparse ${params.sparse_cut} \\
           --out ${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_sparse_grm \\
           --thread-num 5
    gcta64 --grm TO_SIMS_${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_gcta_grm \\
           --pca 1 \\
           --out ${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_sparse_grm \\
           --thread-num 5
    gcta64 --mlma-loco \\
           --bfile TO_SIMS_${SIMREP}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set} \\
           --grm ${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_sparse_grm \\
           --out ${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_lmm-exact \\
           --pheno ${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_sims.phen \\
           --maf ${MAF} \\
           --thread-num 5
    gcta64 --mlma-loco \\
           --bfile TO_SIMS_${SIMREP}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set} \\
           --grm ${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_sparse_grm \\
           --qcovar ${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_sparse_grm.eigenvec \\
           --out ${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_lmm-exact_pca \\
           --pheno ${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_sims.phen \\
           --maf ${MAF} \\
           --thread-num 5

    gcta64 --grm TO_SIMS_${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_gcta_grm_inbred \\
          --make-bK-sparse ${params.sparse_cut} \\
          --out ${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_sparse_grm_inbred \\
          --thread-num 5
    gcta64 --grm TO_SIMS_${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_gcta_grm_inbred \\
          --pca 1 \\
          --out ${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_sparse_grm_inbred \\
          --thread-num 5
    gcta64 --fastGWA-lmm-exact \\
          --grm-sparse ${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_sparse_grm_inbred \\
          --bfile TO_SIMS_${SIMREP}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set} \\
          --out ${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_lmm-exact_inbred \\
          --pheno ${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_sims.phen \\
          --maf ${MAF} \\
          --thread-num 5
    gcta64 --fastGWA-lmm-exact \\
          --grm-sparse ${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_sparse_grm_inbred \\
          --qcovar ${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_sparse_grm_inbred.eigenvec \\
          --bfile TO_SIMS_${SIMREP}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set} \\
          --out ${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_lmm-exact_inbred_pca \\
          --pheno ${SIMREP}_${H2}_${MAF}_${og1}_${og2}_${og3}_${og4}_${og5}_${sp}_${strain_set}_sims.phen \\
          --maf ${MAF} \\
          --thread-num 5
    """
}

process get_gcta_intervals_repeated {
    label 'get_gcta_intervals'

    tag {"${SIMREP} - ${H2}"}

    publishDir "${params.out}/Simulations/${og1}_${og2}_${og3}_${og4}_${og5}_${sp}/Mappings", mode: 'copy', pattern: "*processed_LMM-EXACT-INBRED_mapping.tsv"
    publishDir "${params.out}/Simulations/${og1}_${og2}_${og3}_${og4}_${og5}_${sp}/Mappings", mode: 'copy', pattern: "*processed_LMM-EXACT-INBRED_PCA_mapping.tsv"  
    publishDir "${params.out}/Simulations/${og1}_${og2}_${og3}_${og4}_${og5}_${sp}/Mappings", mode: 'copy', pattern: "*processed_LMM-EXACT-LOCO_mapping.tsv"
    publishDir "${params.out}/Simulations/${og1}_${og2}_${og3}_${og4}_${og5}_${sp}/Mappings", mode: 'copy', pattern: "*processed_LMM-EXACT-LOCO_PCA_mapping.tsv"    
    publishDir "${params.out}/Simulations/${og1}_${og2}_${og3}_${og4}_${og5}_${sp}/Mappings", mode: 'copy', pattern: "*qtl_region.tsv"

    // memory '70 GB'

    input:
    tuple val(sp), val(strain_set), val(NQTL), val(SIMREP), val(H2), file(loci), file(gm), val(effect_range), file(n_indep_tests), val(MAF),val(og1), val(og2), val(og3), val(og4), val(og5), file(lmmexact_inbred), file(lmmexact_inbred_pca), file(lmmexact_loco), file(lmmexact_loco_pca), \
    file(phenotypes), val(THRESHOLD), val(QTL_GROUP_SIZE), val(QTL_CI_SIZE), file(find_gcta_intervals_repeated), file(find_gcta_intervals_loco_repeated)

    output:
    //tuple val(strain_set), val(strains), val(NQTL), val(SIMREP), val(H2), file(loci), file(gm), val(effect_range), file(n_indep_tests), file(phenotypes), val(THRESHOLD), file("*processed_LMM-EXACT-INBRED_mapping.tsv"), emit: processed_gcta_inbred
    //tuple val(strain_set), val(strains), val(MAF), val(NQTL), val(SIMREP), val(H2), val(effect_range), file("*LMM-EXACT-INBRED_qtl_region.tsv"), emit: gcta_qtl_to_ld_inbred    
    //tuple val(strain_set), val(strains), val(NQTL), val(SIMREP), val(H2), file(loci), file(gm), val(effect_range), file(n_indep_tests), file(phenotypes), val(THRESHOLD), file("*processed_LMM-EXACT-INBRED_PCA_mapping.tsv"), emit: processed_gcta_inbred_pca
    //tuple val(strain_set), val(strains), val(MAF), val(NQTL), val(SIMREP), val(H2), val(effect_range), file("*LMM-EXACT-INBRED_PCA_qtl_region.tsv"), emit: gcta_qtl_to_ld_inbred_pca
    //tuple val(strain_set), val(strains), val(NQTL), val(SIMREP), val(H2), file(loci), file(gm), val(effect_range), file(n_indep_tests), file(phenotypes), val(THRESHOLD), file("*processed_LMM-EXACT-LOCO_mapping.tsv"), emit: processed_gcta_loco
    //tuple val(strain_set), val(strains), val(MAF), val(NQTL), val(SIMREP), val(H2), val(effect_range), file("*LMM-EXACT-LOCO_qtl_region.tsv"), emit: gcta_qtl_to_ld_loco
    //tuple val(strain_set), val(strains), val(NQTL), val(SIMREP), val(H2), file(loci), file(gm), val(effect_range), file(n_indep_tests), file(phenotypes), val(THRESHOLD), file("*processed_LMM-EXACT-LOCO_PCA_mapping.tsv"), emit: processed_gcta_loco_pca
    //tuple val(strain_set), val(strains), val(MAF), val(NQTL), val(SIMREP), val(H2), val(effect_range), file("*LMM-EXACT-LOCO_PCA_qtl_region.tsv"), emit: gcta_qtl_to_ld_loco_pca
    //tuple val(strain_set), val(strains), val(MAF), val(NQTL), val(SIMREP), val(H2), val(effect_range), file(loci), file(phenotypes), emit: simulated_phenotypes


    """
        Rscript --vanilla ${find_gcta_intervals_repeated} ${gm} ${phenotypes} ${lmmexact_inbred} ${n_indep_tests} ${sp} ${SIMREP} ${QTL_GROUP_SIZE} ${QTL_CI_SIZE} ${H2} ${params.maf} ${THRESHOLD} ${strain_set} ${MAF} ${effect_range} LMM-EXACT-INBRED
        Rscript --vanilla ${find_gcta_intervals_repeated} ${gm} ${phenotypes} ${lmmexact_inbred_pca} ${n_indep_tests} ${sp} ${SIMREP} ${QTL_GROUP_SIZE} ${QTL_CI_SIZE} ${H2} ${params.maf} ${THRESHOLD} ${strain_set} ${MAF} ${effect_range} LMM-EXACT-INBRED_PCA
        Rscript --vanilla ${find_gcta_intervals_loco_repeated} ${gm} ${phenotypes} ${lmmexact_loco} ${n_indep_tests} ${sp} ${SIMREP} ${QTL_GROUP_SIZE} ${QTL_CI_SIZE} ${H2} ${params.maf} ${THRESHOLD} ${strain_set} ${MAF} ${effect_range} LMM-EXACT-LOCO
        Rscript --vanilla ${find_gcta_intervals_loco_repeated} ${gm} ${phenotypes} ${lmmexact_loco_pca} ${n_indep_tests} ${sp} ${SIMREP} ${QTL_GROUP_SIZE} ${QTL_CI_SIZE} ${H2} ${params.maf} ${THRESHOLD} ${strain_set} ${MAF} ${effect_range} LMM-EXACT-LOCO_PCA
    """
}
