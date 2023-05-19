

process prepare_repeated_simulation_files {

    cpus 4
    memory 30.GB
    time '30m'

    input:
        tuple val(sp), val(strain_set), val(strains), file(vcf), file(index), file(num_chroms), val(MAF)

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
    tuple val(sp), val(strain_set), val(strains), file(vcf), file(index), file(num_chroms), val(MAF)

    output:
        tuple val(sp), val(strain_set), val(strains), file("TO_SIMS.bed"), file("TO_SIMS.bim"), file("TO_SIMS.fam"), file("TO_SIMS.map"), file("TO_SIMS.nosex"), file("TO_SIMS.ped"), file("TO_SIMS.log"), file("${sp}_${strain_set}_${MAF}_Genotype_Matrix.tsv"), val(MAF), emit: sim_geno
        tuple val(sp), val(strain_set), val(strains), val(MAF), file("renamed_chroms.vcf.gz"), file("renamed_chroms.vcf.gz.tbi"), emit: renamed_chrom_vcf_to_ld


    """
    cat > TO_SIMS.bed
    cat > TO_SIMS.bim
    cat > TO_SIMS.fam
    cat > TO_SIMS.map
    cat > TO_SIMS.nosex
    cat > TO_SIMS.ped
    cat > TO_SIMS.log 
    cat > ${sp}_${strain_set}_${MAF}_Genotype_Matrix.tsv
    cat > renamed_chroms.vcf.gz
    cat > renamed_chroms.vcf.gz.tbi
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
    conda '/home/rjm6024/.conda/envs/vcf_stats_1.0'

    input:
        tuple val(sp), val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(MAF), file(n_indep_tests), val(og1), val(og2), val(og3), val(og4), val(og5), val(SIMREP), file(create_causal_qtls), file(master_snps_dir)

    output:
        tuple val(sp), val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(MAF), file(n_indep_tests), val(og1), val(og2), val(og3), val(og4), val(og5), val(SIMREP), file(master_snps_dir), file("${sp}_${strain_set}_${MAF}_${SIMREP}_causal_og_vars.txt")


    """
        python ${create_causal_qtls} ${og1} ${og2} ${og3} ${og4} ${og5} ${bim} ${master_snps_dir} ${sp}
        cat causal_og_vars.txt > ${sp}_${strain_set}_${MAF}_${SIMREP}_causal_og_vars.txt
    """
}