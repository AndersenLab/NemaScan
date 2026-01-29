/*
==============================================
~ > *                                    * < ~
~ ~ > *                                * < ~ ~
~ ~ ~ > *  DOWNLOAD VCF FROM CENDR   * < ~ ~ ~
~ ~ > *                                * < ~ ~
~ > *                                    * < ~
==============================================
*/

process pull_vcf {

    tag {"PULLING VCF FROM CeNDR"}
    label "pull_vcf"

    output:
        path "*hard-filter.isotype.vcf.gz", emit: hard_vcf 
        path "*hard-filter.isotype.vcf.gz.tbi", emit: hard_vcf_index 
        path "*impute.isotype.vcf.gz", emit: impute_vcf 
        path "*impute.isotype.vcf.gz.tbi", emit: impute_vcf_index 
        path "*.strain-annotation*.tsv", emit: ann_vcf

    script:
        """
            wget https://storage.googleapis.com/caendr-site-public-bucket/dataset_release/${params.species}/${params.vcf}/variation/WI.${params.vcf}.small.hard-filter.isotype.vcf.gz
            tabix -p vcf WI.${params.vcf}.small.hard-filter.isotype.vcf.gz

            wget https://storage.googleapis.com/caendr-site-public-bucket/dataset_release/${params.species}/${params.vcf}/variation/WI.${params.vcf}.impute.isotype.vcf.gz
            tabix -p vcf WI.${params.vcf}.impute.isotype.vcf.gz

            wget https://storage.googleapis.com/caendr-site-public-bucket/dataset_release/${params.species}/${params.vcf}/variation/WI.${params.vcf}.strain-annotation.tsv

        """
    
    stub:
        """
        touch vcf WI.${params.vcf}.small.hard-filter.isotype.vcf.gz
        touch vcf WI.${params.vcf}.small.hard-filter.isotype.vcf.gz.tbi
        touch vcf WI.${params.vcf}.impute.isotype.vcf.gz
        touch vcf WI.${params.vcf}.impute.isotype.vcf.gz.tbi
        touch vcf WI.${params.vcf}.strain-annotation.tsv
        """
}


/*
==============================================
~ > *                                    * < ~
~ ~ > *                                * < ~ ~
~ ~ ~ > *  UPDATE ANNOTATION INPUTS  * < ~ ~ ~
~ ~ > *                                * < ~ ~
~ > *                                    * < ~
==============================================
*/

process update_annotations {

    label "xs"
    label "update_annotations"
    label "R"

    publishDir "${save_dir}", mode: 'copy'

    input:
        tuple val(gtf_to_refflat), val(save_dir), file(update_annotations)

    output:
        tuple file("*canonical_geneset.gtf.gz"), file("c_${params.species}_${params.wbb}_refFlat.txt")

    script:
        """
        Rscript --vanilla ${update_annotations} ${params.wbb} ${params.species} ${gtf_to_refflat}
        """

    stub:
        """
        touch canonical_geneset.gtf.gz
        touch c_${params.species}_${params.wbb}_refFlat.txt
        """

}   

/*
===============================================================
~ > *                                                     * < ~
~ ~ > *                                                 * < ~ ~
~ ~ ~ > *  FIX STRAIN NAMES TO MATCH THOSE ON CAENDR  * < ~ ~ ~
~ ~ > *                                                 * < ~ ~
~ > *                                                     * < ~
===============================================================
*/

/*
THIS WILL NEED TO BE UPDATED TO HANDLE OTHER SPECIES
*/


process fix_strain_names_bulk {

    label "md"
    label "fix_strain_names_bulk"
    label "R"

    tag {"BULK TRAIT"}

    publishDir "${params.out}/Phenotypes", mode: 'copy', pattern: "*pr_*.tsv"
    publishDir "${params.out}/Phenotypes", mode: 'copy', pattern: "strain_issues.txt"

    input:
        tuple file(phenotypes), file(isotype_lookup), file(fix_script), val(run_fix)

    output:
        path "pr_*.tsv", emit: fixed_strain_phenotypes 
        path "Phenotyped_Strains.txt", emit: phenotyped_strains_to_analyze 
        path "strain_issues.txt", emit: strain_issues

    script:
        """
        awk '{
            print \$0;
            split(\$2,A,"|");
            for (I in A) {
                printf "%s\\t%s\\t%s\\n", A[I], A[I], \$3;
            }
        }' ${isotype_lookup} > expanded_isotype_lookup.tsv

        # for now, don't fix isotypes for non elegans
        Rscript --vanilla ${fix_script} ${phenotypes} $run_fix expanded_isotype_lookup.tsv

        # check to make sure there are more than 40 strains for a mapping.
        if [[ \$(wc -l <Phenotyped_Strains.txt) -le 40 ]]
        then
            echo "ERROR: Please provide at least 40 strains for a GWAS mapping."
            exit 1
        fi
        """
}

process fix_strain_names_alt {

    label "xs"
    label "fix_strain_names_alt"
    label "R"

    publishDir "${params.out}/Phenotypes", mode: 'copy', pattern: "*.txt"

    input:
        tuple file(phenotypes), file(isotype_lookup), file(fix_script), val(run_fix)

    output:
        path "Phenotyped_Strains.txt", emit: phenotyped_strains_to_analyze 
        file("strain_issues.txt")

    script:
        """
        Rscript --vanilla ${fix_script} ${phenotypes} $run_fix $isotype_lookup
        """
}


/*
===========================================================
~ > *                                                 * < ~
~ ~ > *                                             * < ~ ~
~ ~ ~ > *  CONVERT THE VCF TO A GENOTYPE MATRIX   * < ~ ~ ~
~ ~ > *                                             * < ~ ~
~ > *                                                 * < ~
===========================================================
*/

process vcf_to_geno_matrix {

    label "vcf_to_geno_matrix"
    machineType "n1-standard-2"
    // machineType { 
    //     if (task.attempt == 1) return 'n1-standard-2'
    //     else if (task.attempt == 2) return 'n1-standard-4' 
    //     else return 'n1-standard-6'
    // }

    errorStrategy 'retry'
    maxRetries 3

    publishDir "${params.out}/Genotype_Matrix", mode: 'copy'

    input:
        tuple file(vcf), file(index), file(strains)

    output:
        file("Genotype_Matrix.tsv") 

    """
        bcftools view -S ${strains} -Ou ${vcf} -e 'CHROM="MtDNA"' |\\
        bcftools filter -i N_MISSING=0 -Oz --threads 5 -o Phenotyped_Strain_VCF.vcf.gz
        tabix -p vcf Phenotyped_Strain_VCF.vcf.gz
        # LD pruning
        plink --vcf Phenotyped_Strain_VCF.vcf.gz \\
              --threads 5 \\
              --snps-only \\
              --biallelic-only \\
              --maf ${params.maf} \\
              --set-missing-var-ids @:# \\
              --indep-pairwise 50 10 0.8 \\
              --geno \\
              --allow-extra-chr
        awk -F":" '\$1=\$1' OFS="\\t" plink.prune.in | \\
        sort -k1,1d -k2,2n > markers.txt
        bcftools query -l Phenotyped_Strain_VCF.vcf.gz |\\
        sort > sorted_samples.txt 
        bcftools view -v snps \\
        -S sorted_samples.txt \\
        -R markers.txt -Ou \\
        Phenotyped_Strain_VCF.vcf.gz |\\
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
        sed 's/.\\/./NA/g' > Genotype_Matrix.tsv
    """

}



/*
============================================================
~ > *                                                  * < ~
~ ~ > *                                              * < ~ ~
~ ~ ~ > *  EIGEN DECOMPOSITION OF GENOTYPE MATRIX  * < ~ ~ ~
~ ~ > *                                              * < ~ ~
~ > *                                                  * < ~
============================================================
*/


/*
------------ Decomposition per chromosome
*/

process chrom_eigen_variants {

    tag { CHROM }
    label "chrom_eigen_variants"
    label "R"
    machineType "n1-standard-2"
    // machineType { 
    //     if (task.attempt == 1) return 'n1-standard-2'
    //     else if (task.attempt == 2) return 'n1-standard-4' 
    //     else return 'n1-standard-6'
    // }

    errorStrategy 'retry'
    maxRetries 3

    input:
        tuple val(CHROM), file(genotypes), file(get_genomatrix_eigen)


    output:
        file("${CHROM}_independent_snvs.csv")


    """
        cat ${genotypes} |\\
        awk -v chrom="${CHROM}" '{if(\$1 == "CHROM" || \$1 == chrom) print}' > ${CHROM}_gm.tsv
        # Get independent tests from eigen values for looser p-value correction
        Rscript --vanilla ${get_genomatrix_eigen} ${CHROM}_gm.tsv ${CHROM}
    """

}

/*
------------ Sum independent tests for all chromosomes
*/

process collect_eigen_variants {

    label "local"

    publishDir "${params.out}/Genotype_Matrix", mode: 'copy'

    input:
        file(chrom_tests) //from sig_snps_geno_matrix.collect()

    output:
        file("total_independent_tests.txt") //into independent_tests

    """
        cat *independent_snvs.csv |\\
        grep -v inde |\\
        awk '{s+=\$1}END{print s}' > total_independent_tests.txt
    """

}