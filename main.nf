#! usr/bin/env nextflow

date = new Date().format( 'yyyyMMdd' )

/*
~ ~ ~ > * Parameters: common to all analyses
*/
params.traitfile = null
params.vcf 		 = null
params.help 	 = null
params.e_mem 	 = "100"
params.eigen_mem = params.e_mem + " GB"
params.cendr_v   = "20180527"
params.fix_names = null

/*
~ ~ ~ > * Parameters: for simulations
*/


if (params.simulate){
    /*
    ~ ~ ~ > * number of qtl
    */

    nqtl = Channel.fromPath("${params.simulate_nqtl}")
              .splitCsv()

    /*
    ~ ~ ~ > * heritability
    */
    sim_h2 = Channel.fromPath("${params.simulate_h2}")
                  .splitCsv()

    /*
    ~ ~ ~ > * effect size
    */
    simulate_eff = Channel.fromPath("${params.simulate_eff}")
                  .splitCsv()

    /*
    ~ ~ ~ > * minor allele frequency
    */
    sim_maf = Channel.fromPath("${params.simulate_maf}")
                  .splitCsv()

    /*
    ~ ~ ~ > * strain set
    */

    File pop_file = new File(params.simulate_strains);

    sim_strains = Channel.from(pop_file.collect { it.tokenize( ' ' ) })
                 .map { SM, STRAINS -> [SM, STRAINS] }


}


/*
~ ~ ~ > * Parameters: for burden mapping
*/
params.refflat   = "${params.data_dir}/annotations/c_${params.species}_${params.wbb}_refFlat.txt"
params.freqUpper = 0.05
params.minburden = 2


/*
~ ~ ~ > * Parameters: for GCTA mapping
*/

Channel
	.from("${params.refflat}")
	.set{genes_to_burden}


/*
~ ~ ~ > * OUTPUT DIRECTORY
*/

params.out = "Analysis_Results-${date}"


if (params.help) {
	log.info '''


O~~~     O~~                                   O~~ ~~
O~ O~~   O~~                                 O~~    O~~
O~~ O~~  O~~   O~~    O~~~ O~~ O~~    O~~     O~~         O~~~   O~~    O~~ O~~
O~~  O~~ O~~ O~   O~~  O~~  O~  O~~ O~~  O~~    O~~     O~~    O~~  O~~  O~~  O~~
O~~   O~ O~~O~~~~~ O~~ O~~  O~  O~~O~~   O~~       O~~ O~~    O~~   O~~  O~~  O~~
O~~    O~ ~~O~         O~~  O~  O~~O~~   O~~ O~~    O~~ O~~   O~~   O~~  O~~  O~~
O~~      O~~  O~~~~   O~~~  O~  O~~  O~~ O~~~  O~~ ~~     O~~~  O~~ O~~~O~~~  O~~



    '''
    log.info "----------------------------------------------------------------"
    log.info "                      USAGE                                     "
    log.info "----------------------------------------------------------------"
    log.info ""
    log.info "nextflow main.nf --traitdir=test_bulk --p3d=TRUE --sthresh=BF # download VCF from CeNDR"
    log.info ""
    log.info "Profiles available:"
    log.info "mappings              Profile                Perform GWA mappings with a provided trait file"
    log.info "simulations           Profile                Perform phenotype simulations with GCTA"
    log.info "----------------------------------------------------------------"
    log.info "             -profile annotations USAGE"
    log.info "----------------------------------------------------------------"
    log.info "----------------------------------------------------------------"
    log.info "nextflow main.nf --vcf input_data/elegans/genotypes/WI.20180527.impute.vcf.gz -profile annotations --species elegans"
    log.info "----------------------------------------------------------------"
    log.info "Mandatory arguments:"
    log.info "--wb_build               String                Wormbase build number, must be greater than WS270"
    log.info "--species                String                What species to download information for (elegans, briggsae, or tropicalis)"
    log.info "----------------------------------------------------------------"
    log.info "             -profile mappings USAGE"
    log.info "----------------------------------------------------------------"
    log.info "----------------------------------------------------------------"
    log.info "nextflow main.nf --vcf input_data/elegans/genotypes/WI.20180527.impute.vcf.gz --traitfile input_data/elegans/phenotypes/PC1.tsv -profile mappings --p3d TRUE"
    log.info "----------------------------------------------------------------"
    log.info "----------------------------------------------------------------"
    log.info "Mandatory arguments:"
    log.info "--traitfile              String                Name of file that contains phenotypes. File should be tab-delimited with the columns: strain trait1 trait2 ..."
    log.info "--vcf                    String                Name of VCF to extract variants from. There should also be a tabix-generated index file with the same name in the directory that contains the VCF. If none is provided, the pipeline will download the latest VCF from CeNDR"
    log.info "Optional arguments:"
    log.info "--maf                    String                Minimum minor allele frequency to use for single-marker mapping (Default: 0.05)"
    log.info "--lmm                    String                Perform GCTA mapping with --fastGWA-lmm algorithm (Default: RUN, option to not run is null)"
    log.info "--lmm-exact              String                Perform GCTA mapping with --fastGWA-lmm-exact algorithm (Default: RUN, option to not run is null)"
    log.info "--sparse_cut             String                Any off-diagonal value in the genetic relatedness matrix greater than this is set to 0 (Default: 0.05)"
    log.info "----------------------------------------------------------------"
    log.info "             -profile simulations USAGE"
    log.info "----------------------------------------------------------------"
    log.info "----------------------------------------------------------------"
    log.info "nextflow main.nf --vcf input_data/elegans/genotypes/WI.20180527.impute.vcf.gz -profile simulations"
    log.info "----------------------------------------------------------------"
    log.info "----------------------------------------------------------------"
    log.info "Mandatory arguments:"
    log.info "--simulate_nqtl          File.                 A CSV file with the number of QTL to simulate per phenotype, one value per line (Default is located: input_data/all_species/simulate_nqtl.csv)"
    log.info "--simulate_h2            File                  A CSV file with phenotype heritability, one value per line (Default is located: input_data/all_species/simulate_h2.csv)"
    log.info "Optional arguments:"
    log.info "--simulate_reps          String                The number of replicates to simulate per number of QTL and heritability (Default: 2)"
    log.info "--simulate_maf           File                  A CSV file where each line is a minor allele frequency threshold to test for simulations (Default: input_data/all_species/simulate_maf.csv)"
    log.info "--simulate_eff           File                  A CSV file where each line is an effect size to test for simulations (Default: input_data/all_species/simulate_effect_sizes.csv)"
    log.info "--simulate_strains       File                  A TSV file with two columns: the first is a name for the strain set and the second is a comma-separated strain list without spaces (Default: input_data/all_species/simulate_strains.csv)"
    log.info "--simulate_qtlloc        File                  A BED file with three columns: chromosome name (numeric 1-6), start postion, end postion. The genomic range specified is where markers will be pulled from to simulate QTL (Default: null [which defaults to using the whole genome to randomly simulate a QTL])"
    log.info "----------------------------------------------------------------"
    log.info "----------------------------------------------------------------"
   	log.info "Optional arguments (General):"
   	log.info "--out                    String                Name of folder that will contain the results"
    log.info "--e_mem                  String                Value that corresponds to the amount of memory to allocate for eigen decomposition of chromosomes (DEFAULT = 100)"
    log.info "--cendr_v                String                CeNDR release (DEFAULT = 20180527)"
    log.info "Optional arguments (Marker):"
    log.info "--sthresh                String                Significance threshold for QTL - Options: BF - for bonferroni correction, EIGEN - for SNV eigen value correction, or another number e.g. 4"
    log.info "--group_qtl              Integer               If two QTL are less than this distance from each other, combine the QTL into one, (DEFAULT = 1000)"
    log.info "--ci_size                Integer               Number of SNVs to the left and right of the peak marker used to define the QTL confidence interval, (DEFAULT = 150)"
    log.info "Optional arguments (Burden):"
    log.info "--freqUpper              Float                 Maximum allele frequency for a variant to be considered for burden mapping, (DEFAULT = 0.05)"
    log.info "--minburden              Interger              Minimum number of strains to have a variant for the variant to be considered for burden mapping, (DEFAULT = 2)"
    log.info "--genes                  String                refFlat file format that contains start and stop genomic coordinates for genes of interest, (DEFAULT = bin/gene_ref_flat.Rda)"
    log.info ""
    log.info "--------------------------------------------------------"
    log.info "Information describing the stucture of the input files can be located in input_files/README.txt"
    log.info ""
    log.info ""
    log.info "Flags:"
    log.info "--help                                      Display this message"
    log.info ""
    log.info "--------------------------------------------------------"
    log.info ""
    log.info " Required software packages to be in users path"
    log.info "BCFtools               v1.9"
    log.info "plink                  v1.9"
    log.info "R-cegwas2              Found on GitHub"
    log.info "R-tidyverse            v1.2.1"
    log.info "R-correlateR           Found on GitHub"
    log.info "R-rrBLUP               v4.6"
    log.info "R-sommer               v3.5"
    log.info "R-RSpectra             v0.13-1"
    log.info "R-ggbeeswarm           v0.6.0"
    log.info "--------------------------------------------------------"
    exit 1
} else {
	log.info '''


O~~~     O~~                                   O~~ ~~
O~ O~~   O~~                                 O~~    O~~
O~~ O~~  O~~   O~~    O~~~ O~~ O~~    O~~     O~~         O~~~   O~~    O~~ O~~
O~~  O~~ O~~ O~   O~~  O~~  O~  O~~ O~~  O~~    O~~     O~~    O~~  O~~  O~~  O~~
O~~   O~ O~~O~~~~~ O~~ O~~  O~  O~~O~~   O~~       O~~ O~~    O~~   O~~  O~~  O~~
O~~    O~ ~~O~         O~~  O~  O~~O~~   O~~ O~~    O~~ O~~   O~~   O~~  O~~  O~~
O~~      O~~  O~~~~   O~~~  O~  O~~  O~~ O~~~  O~~ ~~     O~~~  O~~ O~~~O~~~  O~~



    '''
log.info ""
log.info "Trait File                              = ${params.maps}"
log.info "VCF                                     = ${params.vcf}"
log.info "Significance Threshold                  = ${params.sthresh}"
log.info "Result Directory                        = ${params.out}"
log.info "Eigen Memory allocation                 = ${params.eigen_mem}"
log.info ""
}

/*
~ ~ ~ > * COMBINE VCF AND VCF INDEX INTO A CHANNEL
*/

if (params.vcf) {

    vcf = Channel.fromPath("${params.vcf}")

    vcf_index = Channel.fromPath("${params.vcf}" + ".tbi")

    vcf
        .spread(vcf_index)
        .into{vcf_to_names;
              vcf_to_simulations;
              vcf_to_whole_gcta_grm;
              vcf_to_whole_gcta_map;
              vcf_to_whole_genome;
              vcf_to_extract_ld;
              vcf_to_fine_map;
              vcf_to_burden;
              vcf_to_query_vcf}

} else {

    process pull_vcf {

        tag {"PULLING VCF FROM CeNDR"}
        executor 'local'

        output:
            file("*.vcf.gz") into dl_vcf
            file("*.vcf.gz.tbi") into dl_vcf_index

        """
            wget https://storage.googleapis.com/elegansvariation.org/releases/${params.cendr_v}/variation/WI.${params.cendr_v}.impute.vcf.gz
            tabix -p vcf WI.${params.cendr_v}.impute.vcf.gz
        """
    }

    dl_vcf
        .spread(dl_vcf_index)
        .into{vcf_to_names;
              vcf_to_simulations;
              vcf_to_whole_gcta_grm;
              vcf_to_whole_gcta_map;
              vcf_to_whole_genome;
              vcf_to_extract_ld;
              vcf_to_fine_map;
              vcf_to_burden;
              vcf_to_query_vcf}

}

/*
~ ~ ~ > * INITIATE MAPPING QTL GROUPING PARAMETER
*/

Channel
    .from("${params.ci_size}")
    .into{qtl_ci_size;
          qtl_ci_size_sim;
          qtl_ci_size_sim_gcta}

/*
~ ~ ~ > * INITIATE MAPPING QTL CONFIDENCE INTERVAL SIZE PARAMETER
*/

Channel
    .from("${params.group_qtl}")
    .into{qtl_snv_grouping_maps;
          qtl_snv_grouping_sims;
          qtl_snv_grouping_sims_gcta}

/*
~ ~ ~ > * INITIATE MAPPING METHOD CHANNEL
*/

Channel
    .from("${params.p3d}")
    .into{p3d_full;
          p3d_fine;
          p3d_full_sim;
          p3d_fine_sim}

/*
~ ~ ~ > * INITIATE THRESHOLD CHANNEL
*/

Channel
    .from("${params.sthresh}")
    .into{sig_threshold_full;
          sig_threshold_fine;
          sig_threshold_full_sim;
          sig_threshold_fine_sim;
          sig_threshold_gcta_sim}

/*
~ ~ ~ > * INITIATE PHENOTYPE CHANNEL
*/

Channel
    .fromPath("${params.traitfile}")
    .set{ traits_to_strainlist }

/*
~ ~ ~ > * INITIATE CHROMOSOME RENAME CHANNEL
*/

Channel
    .fromPath("${params.numeric_chrom}")
    .into{ rename_chroms_gcta;
           rename_chroms_sims;
           rename_chroms_sims_ld }


/*
==============================================
~ > *                                    * < ~
~ ~ > *                                * < ~ ~
~ ~ ~ > *  UPDATE ANNOTATION INPUTS  * < ~ ~ ~
~ ~ > *                                * < ~ ~
~ > *                                    * < ~
==============================================
*/

if (params.annotate) {

save_dir = "${params.input_data}/${params.species}/annotations"

Channel
    .fromPath("${params.script_loc}")
    .set{path_to_converter}

    process update_annotations {

    executor 'local'

    publishDir "${save_dir}", mode: 'copy'

    input:
        val(gtf_to_refflat) from path_to_converter
        val(save_dir)

    output:
        set file("*canonical_geneset.gtf.gz"), file("c_${params.species}_${params.wb_build}_refFlat.txt") into updated_annotations

    when:
        params.annotate

    """
        Rscript --vanilla `which update_annotations.R` ${params.wb_build} ${params.species} ${gtf_to_refflat}
    """

    }
}

/*
==============================================================
~ > *                                                    * < ~
~ ~ > *                                                * < ~ ~
~ ~ ~ > *  FIX STRAIN NAMES TO MATCH THOSE ON CENDR  * < ~ ~ ~
~ ~ > *                                                * < ~ ~
~ > *                                                    * < ~
==============================================================
*/

/*
THIS WILL NEED TO BE UPDATED TO HANDLE OTHER SPECIES
*/


if (params.maps) {

    process fix_strain_names_bulk {

    executor 'local'

    input:
        file(phenotypes) from traits_to_strainlist

    output:
        file("pr_*.tsv") into fixed_strain_phenotypes
        file("Phenotyped_Strains.txt") into phenotyped_strains_to_analyze

    when:
        params.maps

    """
        Rscript --vanilla `which Fix_Isotype_names_bulk.R` ${phenotypes} ${params.fix_names}
    """

}

fixed_strain_phenotypes
    .flatten()
    .map { file -> tuple(file.baseName.replaceAll(/pr_/,""), file) }
    .into{ traits_to_map;
           traits_to_burden;
           traits_to_gcta_grm;
           traits_to_gcta_map }

phenotyped_strains_to_analyze
    .into{strain_list_genome;
          strain_list_emma;
          strain_list_finemap;
          strain_list_gcta_grml;
          strain_list_gcta_map}


} else if (params.simulate) {

    simulation_replicates = Channel.from(1..params.simulate_reps)

/*
    process simulate_strain_names {

        executor 'local'

        input:
            set file(vcf), file(index), val(strain_set) from vcf_to_names

        output:
            file("sorted_samples.txt") into phenotyped_strains_to_analyze

        when:
            params.simulate

        """
        bcftools view -s ${strains} ${vcf} |\\
        bcftools query -l |\\
        sort > sorted_samples.txt
        """
    }
*/


sim_strains
    .into{strain_list_genome;
          strain_list_emma;
          strain_list_simulate;
          strain_list_finemap;
          strain_list_gcta_grml;
          strain_list_gcta_map}

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

if (params.maps) {

    process vcf_to_geno_matrix {

        publishDir "${params.out}/Genotype_Matrix", mode: 'copy'

        cpus 1

        input:
            set file(vcf), file(index) from vcf_to_whole_genome
            file(strains) from strain_list_genome

        output:
            file("Genotype_Matrix.tsv") into geno_matrix

        when:
            params.maps

        """

            bcftools view -S ${strains} ${vcf} |\\
            bcftools filter -i N_MISSING=0 -Oz -o Phenotyped_Strain_VCF.vcf.gz

            tabix -p vcf Phenotyped_Strain_VCF.vcf.gz

            plink --vcf Phenotyped_Strain_VCF.vcf.gz \\
                --snps-only \\
                --biallelic-only \\
                --maf 0.05 \\
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
            -R markers.txt \\
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

    geno_matrix
        .into{eigen_gm;
              mapping_gm}


    /*
    ============================================================
    ~ > *                                                  * < ~
    ~ ~ > *                                              * < ~ ~
    ~ ~ ~ > *  EIGEN DECOMPOSITION OF GENOTYPE MATRIX  * < ~ ~ ~
    ~ ~ > *                                              * < ~ ~
    ~ > *                                                  * < ~
    ============================================================
    */

    CONTIG_LIST = ["I", "II", "III", "IV", "V", "X"]
    contigs = Channel.from(CONTIG_LIST)

    /*
    ------------ Decomposition per chromosome
    */

    process chrom_eigen_variants {

        tag { CHROM }

        cpus 6
        memory params.eigen_mem

        input:
            file(genotypes) from eigen_gm
            each CHROM from contigs

        output:
            file("${CHROM}_independent_snvs.csv") into sig_snps_geno_matrix
            file(genotypes) into concat_geno

        when:
            params.maps


        """
            cat Genotype_Matrix.tsv |\\
            awk -v chrom="${CHROM}" '{if(\$1 == "CHROM" || \$1 == chrom) print}' > ${CHROM}_gm.tsv
            Rscript --vanilla `which Get_GenoMatrix_Eigen.R` ${CHROM}_gm.tsv ${CHROM}
        """

    }

    /*
    ------------ Sum independent tests for all chromosomes
    */

    process collect_eigen_variants {

        publishDir "${params.out}/Genotype_Matrix", mode: 'copy'

        cpus 1

        input:
            file(chrom_tests) from sig_snps_geno_matrix.collect()

        output:
            file("total_independent_tests.txt") into independent_tests

        when:
            params.maps

        """
            cat *independent_snvs.csv |\\
            grep -v inde |\\
            awk '{s+=\$1}END{print s}' > total_independent_tests.txt
        """

    }


independent_tests
    .spread(mapping_gm)
    .spread(traits_to_map)
    .spread(p3d_full)
    .spread(sig_threshold_full)
    .spread(qtl_snv_grouping_maps)
    .spread(qtl_ci_size)
    .into{mapping_data_emma;
          mapping_data_gcta}

}



/*
======================================
~ > *                            * < ~
~ ~ > *                        * < ~ ~
~ ~ ~ > *  GWAS SIMULATIONS  * < ~ ~ ~
~ ~ > *                        * < ~ ~
~ > *                            * < ~
======================================
*/


if(params.maps){
    strain_list_emma
        .spread(traits_to_gcta_grm)
        .spread(vcf_to_whole_gcta_grm)
        .spread(rename_chroms_gcta)
        .into{gcta_prep_inputs;
              print_inputs}

} else if (params.simulate){
    strain_list_simulate
    .spread(vcf_to_simulations)
    .spread(rename_chroms_sims)
    .spread(sim_maf)
    .into{simulation_prep_inputs;
          print_inputs}
}


if(params.simulate){

    process prepare_simulation_files {

        cpus 4

        input:
            set val(strain_set), val(strains), file(vcf), file(index), file(num_chroms), val(MAF) from simulation_prep_inputs

        output:
            set val(strain_set), val(strains), file("TO_SIMS.bed"), file("TO_SIMS.bim"), file("TO_SIMS.fam"), file("TO_SIMS.map"), file("TO_SIMS.nosex"), file("TO_SIMS.ped"), file("TO_SIMS.log"), file("${strain_set}_${MAF}_Genotype_Matrix.tsv"), val(MAF) into sim_geno
            set val(strain_set), val(strains), val(MAF), file("renamed_chroms.vcf.gz"), file("renamed_chroms.vcf.gz.tbi") into renamed_chrom_vcf_to_ld

        when:
            params.simulate

        """

        bcftools annotate --rename-chrs rename_chromosomes ${vcf} |\\
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
        sed 's/.\\/./NA/g' > ${strain_set}_${MAF}_Genotype_Matrix.tsv

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

    CONTIG_LIST = ["1", "2", "3", "4", "5", "6"]
    contigs = Channel.from(CONTIG_LIST)

    /*
    ------------ Decomposition per chromosome
    */

    process chrom_eigen_variants_sims {

        tag { CHROM }

        cpus 6
        memory params.eigen_mem

        input:
            set val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(sex), file(ped), file(log), file(geno), val(MAF) from sim_geno
            each CHROM from contigs

        output:
            set val(strain_set), val(strains), val(MAF), file(bed), file(bim), file(fam), file(map), file(sex), file(ped), file(log), file(geno) into sim_geno_meta
            set val(strain_set), val(strains), val(MAF), file("${CHROM}_${strain_set}_${MAF}_independent_snvs.csv") into sim_geno_eigen_join

        when:
            params.simulate


        """
            cat ${geno} |\\
            awk -v chrom="${CHROM}" '{if(\$1 == "CHROM" || \$1 == chrom) print}' > ${CHROM}_gm.tsv
            Rscript --vanilla `which Get_GenoMatrix_Eigen.R` ${CHROM}_gm.tsv ${CHROM}

            mv ${CHROM}_independent_snvs.csv ${CHROM}_${strain_set}_${MAF}_independent_snvs.csv
        """

    }

    /*
    ------------ Sum independent tests for all chromosomes
    */

    sim_geno_eigen_join.groupTuple(by:[0,1,2]).join(sim_geno_meta, by:[0,1,2]).into{combine_independent_tests;combine_independent_tests_print}

    process collect_eigen_variants_sims {

        executor 'local'

        publishDir "${params.out}/Genotype_Matrix", mode: 'copy'

        cpus 1

        input:
            set val(strain_set), val(strains), val(MAF), file(tests), file(bed), file(bim), file(fam), file(map), file(sex), file(ped), file(log), file(geno) from combine_independent_tests


        output:
            set val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(sex), file(ped), file(log), file(geno), val(MAF), file("${strain_set}_${MAF}_total_independent_tests.txt") into sim_inputs

        when:
            params.simulate

        """
            cat *independent_snvs.csv |\\
            grep -v inde |\\
            awk '{s+=\$1}END{print s}' > ${strain_set}_${MAF}_total_independent_tests.txt
        """

    }


    sim_inputs
        .spread(nqtl)
        .into{sim_nqtl_inputs_loc;
             sim_nqtl_inputs}

if(params.simulate_qtlloc){

    qtl_locations = Channel.fromPath("${params.simulate_qtlloc}")

    sim_nqtl_inputs_loc
        .spread(qtl_locations)
        .spread(simulate_eff)
        .set{sim_nqtl_inputs_qtl_locations}

    process simulate_effects_loc {

        tag {NQTL}

        cpus 4

        input:
            set val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(MAF), file(n_indep_tests), val(NQTL), file(qtl_loc_bed), val(effect_range) from sim_nqtl_inputs_qtl_locations
            each SIMREP from simulation_replicates

        output:
            set val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(MAF), file(n_indep_tests), val(NQTL), val(SIMREP), val(effect_range), file("causal.variants.sim.${NQTL}.${SIMREP}.txt") into sim_phen_inputs

        when:
            params.simulate
            params.simulate_qtlloc

        """

         Rscript --vanilla `which create_causal_QTLs.R` ${bim} ${NQTL} ${effect_range} ${qtl_loc_bed}

         mv causal.variants.sim.${NQTL}.txt causal.variants.sim.${NQTL}.${SIMREP}.txt

        """
    }

} else {

        sim_nqtl_inputs
            .spread(simulate_eff)
            .set{sim_nqtl_inputs_with_effects}

        process simulate_effects_genome {

        tag {NQTL}

        cpus 4

        input:
            set val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(MAF), file(n_indep_tests), val(NQTL), val(effect_range) from sim_nqtl_inputs_with_effects
            each SIMREP from simulation_replicates

        output:
            set val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(MAF), file(n_indep_tests), val(NQTL), val(SIMREP), val(effect_range), file("causal.variants.sim.${NQTL}.${SIMREP}.txt") into sim_phen_inputs

        when:
            params.simulate

        """

         Rscript --vanilla `which create_causal_QTLs.R` ${bim} ${NQTL} ${effect_range}

         mv causal.variants.sim.${NQTL}.txt causal.variants.sim.${NQTL}.${SIMREP}.txt

        """
    }

}



    sim_phen_inputs
        .spread(sim_h2)
        .set{sim_phen_h2_input}

/* SIMULATIONS INCLUDING LOCO+INBRED and LMM-EXACT NON-INBRED
    process simulate_map_phenotypes {

        tag {"${NQTL} - ${SIMREP} - ${H2} - ${MAF}"}

        publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Mappings", mode: 'copy', pattern: "*fastGWA", overwrite: true
        publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Mappings", mode: 'copy', pattern: "*loco.mlma", overwrite: true
        publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Phenotypes", mode: 'copy', pattern: "*.phen", overwrite: true
        publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Phenotypes", mode: 'copy', pattern: "*.par", overwrite: true

        cpus 4

        errorStrategy 'ignore'

        input:
            set val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(MAF), file(n_indep_tests), val(NQTL), val(SIMREP), val(effect_range), file(loci), val(H2) from sim_phen_h2_input

        output:
            set file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.bed"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.bim"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.fam"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.map"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.nosex"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.ped"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.log"), val(NQTL), val(SIMREP), file(loci), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.par") into sim_phen_output
            set file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact.fastGWA"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred.fastGWA"),file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact.loco.mlma"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred.loco.mlma"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact.log"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred.log") into sim_GCTA_mapping_results

            set val(strain_set), val(strains), val(NQTL), val(SIMREP), val(H2), file(loci), file(gm), val(MAF), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.par"), val(effect_range), file(n_indep_tests) into sim_phen_to_emma
            file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact.fastGWA") into lmm_exact_analyze_sims
            file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred.fastGWA") into lmm_exact_inbred_analyze_sims
            file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact.loco.mlma") into lmm_exact_loco_analyze_sims
            file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred.loco.mlma") into lmm_exact_loco_inbred_analyze_sims
            file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen") into simphen_analyze_sims
            file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.par") into simgen_analyze_sims
            set val(strain_set), val(strains), val(NQTL), val(SIMREP), val(H2), file(loci), file(gm), val(effect_range), file(n_indep_tests), val(MAF), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact.fastGWA"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred.fastGWA"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact.loco.mlma"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred.loco.mlma"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen") into gcta_intervals

        when:
            params.simulate

        """

        gcta64 --bfile TO_SIMS \\
             --simu-qt \\
             --simu-causal-loci ${loci} \\
             --simu-hsq ${H2} \\
             --simu-rep 1 \\
             --out ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims

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
            --pheno ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen

        gcta64 --bfile TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set} --autosome --maf ${MAF} --make-grm --out TO_SIMS_${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_gcta_grm --thread-num 10
        gcta64 --bfile TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set} --autosome --maf ${MAF} --make-grm-inbred --out TO_SIMS_${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_gcta_grm_inbred --thread-num 10

        gcta64 --grm TO_SIMS_${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_gcta_grm_inbred --pheno ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen --reml --out check_vp --thread-num 10
        vp=`grep Vp check_vp.hsq | head -1 | cut -f2`
        if (( \$(echo "0.00001 > \$vp" |bc -l) ));
        then
        awk '{print \$1, \$2, \$3*1000}' ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen > temp.phen;
        rm ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen
        mv temp.phen ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen
        fi



        gcta64 --grm TO_SIMS_${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_gcta_grm --make-bK-sparse ${params.sparse_cut} --out ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sparse_grm
        """gcta64 --fastGWA-lmm-exact \\
            --grm-sparse ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sparse_grm \\
            --bfile TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set} \\
            --out ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact \\
            --pheno ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen \\
            --maf ${MAF}"""
        gcta64 --mlma-loco \\
            --bfile TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set} \\
            --grm ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sparse_grm \\
            --out ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact \\
            --pheno ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen \\
            --maf ${MAF}


        gcta64 --grm TO_SIMS_${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_gcta_grm_inbred --make-bK-sparse ${params.sparse_cut} --out ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sparse_grm_inbred
        gcta64 --fastGWA-lmm-exact \\
            --grm-sparse ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sparse_grm_inbred \\
            --bfile TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set} \\
            --out ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred \\
            --pheno ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen \\
            --maf ${MAF}

        """gcta64 --mlma-loco \\
            --bfile TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set} \\
            --grm ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sparse_grm_inbred \\
            --out ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred \\
            --pheno ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen \\
            --maf ${MAF}"""

        """
    }

    */

    process simulate_map_phenotypes {

        tag {"${NQTL} - ${SIMREP} - ${H2} - ${MAF}"}

        publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Mappings", pattern: "*fastGWA", overwrite: true
        publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Mappings", pattern: "*loco.mlma", overwrite: true
        publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Phenotypes", pattern: "*.phen", overwrite: true
        publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Phenotypes", pattern: "*.par", overwrite: true

        cpus 4

        errorStrategy 'ignore'

        input:
            set val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(MAF), file(n_indep_tests), val(NQTL), val(SIMREP), val(effect_range), file(loci), val(H2) from sim_phen_h2_input

        output:
            set file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.bed"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.bim"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.fam"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.map"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.nosex"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.ped"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.log"), val(NQTL), val(SIMREP), file(loci), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.par") into sim_phen_output
            set file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred.fastGWA"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact.loco.mlma"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact.log"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred.log") into sim_GCTA_mapping_results
            file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred.fastGWA") into lmm_exact_inbred_analyze_sims
            file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact.loco.mlma") into lmm_exact_loco_analyze_sims
            file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen") into simphen_analyze_sims
            file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.par") into simgen_analyze_sims
            set val(strain_set), val(strains), val(NQTL), val(SIMREP), val(H2), file(loci), file(gm), val(effect_range), file(n_indep_tests), val(MAF), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred.fastGWA"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact.loco.mlma"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen") into gcta_intervals

        when:
            params.simulate

        """

        gcta64 --bfile TO_SIMS \\
             --simu-qt \\
             --simu-causal-loci ${loci} \\
             --simu-hsq ${H2} \\
             --simu-rep 1 \\
             --out ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims

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
            --pheno ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen

        gcta64 --bfile TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set} --autosome --maf ${MAF} --make-grm --out TO_SIMS_${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_gcta_grm --thread-num 10
        gcta64 --bfile TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set} --autosome --maf ${MAF} --make-grm-inbred --out TO_SIMS_${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_gcta_grm_inbred --thread-num 10

        gcta64 --grm TO_SIMS_${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_gcta_grm_inbred --pheno ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen --reml --out check_vp --thread-num 10
        vp=`grep Vp check_vp.hsq | head -1 | cut -f2`
        if (( \$(echo "0.00001 > \$vp" |bc -l) ));
        then
        awk '{print \$1, \$2, \$3*1000}' ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen > temp.phen;
        rm ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen
        mv temp.phen ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen
        fi



        gcta64 --grm TO_SIMS_${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_gcta_grm --make-bK-sparse ${params.sparse_cut} --out ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sparse_grm
        gcta64 --mlma-loco \\
            --bfile TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set} \\
            --grm ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sparse_grm \\
            --out ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact \\
            --pheno ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen \\
            --maf ${MAF}

        gcta64 --grm TO_SIMS_${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_gcta_grm_inbred --make-bK-sparse ${params.sparse_cut} --out ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sparse_grm_inbred
        gcta64 --fastGWA-lmm-exact \\
            --grm-sparse ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sparse_grm_inbred \\
            --bfile TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set} \\
            --out ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred \\
            --pheno ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen \\
            --maf ${MAF}

        """
    }

/*
    sim_phen_to_emma
        .spread(qtl_snv_grouping_sims)
        .spread(qtl_ci_size_sim)
        .spread(p3d_full_sim)
        .spread(sig_threshold_full_sim)
        .into{sim_emma_inputs;
              sim_emma_fine_inputs}
*/


/*
------------ For simulations, i want to separate mapping and defining threshold, so we can look at various threshold after the long mapping step
------------ further optimization can be made by making kinship matrix outside of mapping process
*/

/*
    process sim_emmma_maps {

        memory '48 GB'

        publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Mappings", mode: 'copy', pattern: "*processed_mapping.tsv"

        input:
        set val(strain_set), val(strains), val(NQTL), val(SIMREP), val(H2), file(loci), file(gm), val(MAF), file(pheno), file(sim_params), val(effect_range), file(n_indep_tests), val(QTL_GROUP_SIZE), val(QTL_CI_SIZE), val(P3D), val(THRESHOLD) from sim_emma_inputs

        output:
        set val(NQTL), val(SIMREP), val(H2), file("*raw_mapping.tsv"), file("*processed_mapping.tsv") into pr_sim_emma_maps
        set val(strain_set), val(strains), val(MAF), val(NQTL), val(SIMREP), val(H2), val(effect_range), file("*_emma_qtl_region.tsv") into emmma_qtl_to_ld

        """

        Rscript --vanilla `which Run_Sims_EMMA.R` ${gm} ${pheno} ${task.cpus} ${P3D} ${NQTL} ${SIMREP} ${QTL_GROUP_SIZE} ${QTL_CI_SIZE} ${H2} ${params.maf} ${THRESHOLD} ${strain_set} ${MAF} ${effect_range}

        """
    }
*/
    gcta_intervals
        .spread(sig_threshold_gcta_sim)
        .spread(qtl_snv_grouping_sims_gcta)
        .spread(qtl_ci_size_sim_gcta)
        .into{find_gcta_intervals}

    process get_gcta_intervals {

    publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Mappings", mode: 'copy', pattern: "*processed_aggregate_mapping.tsv"
    publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Mappings", mode: 'copy', pattern: "*processed_LMM-EXACT-INBRED_mapping.tsv"
    publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Mappings", mode: 'copy', pattern: "*processed_LMM-EXACT-LOCO_mapping.tsv"
    publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Mappings", mode: 'copy', pattern: "*qtl_region.tsv"
    publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Mappings", mode: 'copy', pattern: "*LD.tsv"

    memory '48 GB'

    errorStrategy 'ignore'

    input:
    set val(strain_set), val(strains), val(NQTL), val(SIMREP), val(H2), file(loci), file(gm), val(effect_range), file(n_indep_tests), val(MAF), file(lmmexact_inbred), file(lmmexact_loco), file(phenotypes), val(THRESHOLD), val(QTL_GROUP_SIZE), val(QTL_CI_SIZE) from find_gcta_intervals

    output:
    set val(strain_set), val(strains), val(NQTL), val(SIMREP), val(H2), file(loci), file(gm), val(effect_range), file(n_indep_tests), file(phenotypes), val(THRESHOLD), file("*processed_aggregate_mapping.tsv"), file("*processed_LMM-EXACT-INBRED_mapping.tsv"), file("*processed_LMM-EXACT-LOCO_mapping.tsv") into processed_gcta
    set val(strain_set), val(strains), val(MAF), val(NQTL), val(SIMREP), val(H2), val(effect_range), file("*aggregate_qtl_region.tsv") into gcta_qtl_to_ld
    set val(strain_set), val(strains), val(MAF), val(NQTL), val(SIMREP), val(H2), val(effect_range), file(loci), file(phenotypes) into simulated_phenotypes

    """
    Rscript --vanilla `which Aggregate_Mappings.R` ${lmmexact_loco} ${lmmexact_inbred}
    Rscript --vanilla `which Find_Aggregate_Intervals.R` ${gm} ${phenotypes} temp.aggregate.mapping.tsv ${n_indep_tests} ${NQTL} ${SIMREP} ${QTL_GROUP_SIZE} ${QTL_CI_SIZE} ${H2} ${params.maf} ${THRESHOLD} ${strain_set} ${MAF} ${effect_range} aggregate
    Rscript --vanilla `which Find_GCTA_Intervals.R` ${gm} ${phenotypes} ${lmmexact_inbred} ${n_indep_tests} ${NQTL} ${SIMREP} ${QTL_GROUP_SIZE} ${QTL_CI_SIZE} ${H2} ${params.maf} ${THRESHOLD} ${strain_set} ${MAF} ${effect_range} LMM-EXACT-INBRED
    Rscript --vanilla `which Find_GCTA_Intervals_LOCO.R` ${gm} ${phenotypes} ${lmmexact_loco} ${n_indep_tests} ${NQTL} ${SIMREP} ${QTL_GROUP_SIZE} ${QTL_CI_SIZE} ${H2} ${params.maf} ${THRESHOLD} ${strain_set} ${MAF} ${effect_range} LMM-EXACT-LOCO

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
}

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

/*
------------ FOR SOME DUMB REASON THE WHEN DIRECTIVE IS NOT WORKING AS IT SHOULD, SO I HAVE WRAPPED THESE PROCESSES IN A CONDITIONAL
*/

if (params.maps) {

    process prepare_gcta_files {

        cpus 4

        input:
            set file(strains), val(TRAIT), file(traits), file(vcf), file(index), file(num_chroms) from gcta_prep_inputs

        output:
            set val(TRAIT), file("plink_formated_trats.tsv"), file("${TRAIT}.bed"), file("${TRAIT}.bim"), file("${TRAIT}.fam"), file("${TRAIT}.map"), file("${TRAIT}.nosex"), file("${TRAIT}.ped"), file("${TRAIT}.log") into gcta_grm_inputs

        when:
            params.maps == "RUN"

        """

        bcftools annotate --rename-chrs rename_chromosomes ${vcf} |\\
        bcftools view -S ${strains} |\\
        bcftools filter -i N_MISSING=0 -Oz -o renamed_chroms.vcf.gz

        tabix -p vcf renamed_chroms.vcf.gz

        plink --vcf renamed_chroms.vcf.gz \\
        --snps-only \\
        --biallelic-only \\
        --maf 0.05 \\
        --set-missing-var-ids @:# \\
        --indep-pairwise 50 10 0.8 \\
        --geno \\
        --not-chr MtDNA \\
        --allow-extra-chr

        tail -n +2 ${traits} | awk 'BEGIN {OFS="\\t"}; {print \$1, \$1, \$2}' > plink_formated_trats.tsv

        plink --vcf renamed_chroms.vcf.gz \\
        --make-bed \\
        --snps-only \\
        --biallelic-only \\
        --maf 0.05 \\
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

        cpus 4

        input:
            set val(TRAIT), file(traits), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log) from gcta_grm_inputs

        output:
            set val(TRAIT), file(traits), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file("${TRAIT}_gcta_grm.grm.bin"), file("${TRAIT}_gcta_grm.grm.id"), file("${TRAIT}_gcta_grm.grm.N.bin"), file("${TRAIT}_heritability.hsq"), file("${TRAIT}_heritability.log"), file("${TRAIT}_gcta_grm_inbred.grm.bin"), file("${TRAIT}_gcta_grm_inbred.grm.id"), file("${TRAIT}_gcta_grm_inbred.grm.N.bin"), file("${TRAIT}_heritability_inbred.hsq"), file("${TRAIT}_heritability_inbred.log") into gcta_mapping_inputs

        when:
            params.maps

        """

        gcta64 --bfile ${TRAIT} --autosome --maf 0.05 --make-grm --out ${TRAIT}_gcta_grm --thread-num 10
        gcta64 --bfile ${TRAIT} --autosome --maf 0.05 --make-grm-inbred --out ${TRAIT}_gcta_grm_inbred --thread-num 10

        gcta64 --grm ${TRAIT}_gcta_grm --pheno plink_formated_trats.tsv --reml --out ${TRAIT}_heritability --thread-num 10
        gcta64 --grm ${TRAIT}_gcta_grm_inbred --pheno plink_formated_trats.tsv --reml --out ${TRAIT}_heritability_inbred --thread-num 10

        """
    }

    gcta_mapping_inputs
        .into{gcta_lmm_exact;
              gcta_lmm;
              gcta_mlma_loco}

    process gcta_lmm_exact_mapping {

        cpus 4

        publishDir "${params.out}/Mapping/Raw", pattern: "*fastGWA", overwrite: true
        publishDir "${params.out}/Mapping/Raw", pattern: "*loco.mlma", overwrite: true

        errorStrategy 'ignore'

        input:
        set val(TRAIT), file(traits), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(grm_bin), file(grm_id), file(grm_nbin), file(h2), file(h2log), file(grm_bin_inbred), file(grm_id_inbred), file(grm_nbin_inbred), file(h2_inbred), file(h2log_inbred) from gcta_lmm_exact

        output:
        set file("${TRAIT}_lmm-exact_inbred.fastGWA"), file("${TRAIT}_lmm-exact.loco.mlma") into lmm_exact_output


        when:
          params.lmm_exact

        """

        gcta64 --grm ${TRAIT}_gcta_grm --make-bK-sparse ${params.sparse_cut} --out ${TRAIT}_sparse_grm
        gcta64 --mlma-loco \\
            --grm ${TRAIT}_sparse_grm \\
            --bfile ${TRAIT} \\
            --out ${TRAIT}_lmm-exact \\
            --pheno ${traits} \\
            --maf ${params.maf}

        gcta64 --grm ${TRAIT}_gcta_grm_inbred --make-bK-sparse ${params.sparse_cut} --out ${TRAIT}_sparse_grm_inbred
        gcta64 --fastGWA-lmm-exact \\
            --grm-sparse ${TRAIT}_sparse_grm \\
            --bfile ${TRAIT} \\
            --out ${TRAIT}_lmm-exact_inbred \\
            --pheno ${traits} \\
            --maf ${params.maf}

        """
    }


    lmm_exact_output
        .into{find_gcta_intervals}

    process gcta_intervals_maps {

    publishDir "${params.out}/Mapping/Processed", mode: 'copy', pattern: "*AGGREGATE_mapping.tsv"
    publishDir "${params.out}/Mapping/Processed", mode: 'copy', pattern: "*AGGREGATE_qtl_region.tsv"

    memory '48 GB'

    input:
    set file(lmmexact_inbred), file(lmmexact_loco) from find_gcta_intervals
    set file(tests), file(geno), val(TRAIT), file(pheno), val(P3D), val(sig_thresh), val(qtl_grouping_size), val(qtl_ci_size) from mapping_data_gcta

    output:
    set file(geno), file(pheno), file("*AGGREGATE_mapping.tsv"), file("*AGGREGATE_qtl_region.tsv") into processed_gcta

    """
    Rscript --vanilla `which Aggregate_Mappings.R` ${lmmexact_loco} ${lmmexact_inbred}
    Rscript --vanilla `which Find_Aggregate_Intervals_Maps.R` ${geno} ${pheno} temp.aggregate.mapping.tsv ${tests} ${qtl_grouping_size} ${qtl_ci_size} ${sig_thresh} ${TRAIT}_AGGREGATE

    """
}


  processed_gcta
      .into{plotting_data}

  process generate_plots {


    publishDir "${params.out}/Plots/LDPlots", mode: 'copy', pattern: "*_LD.plot.png"
    publishDir "${params.out}/Plots/EffectPlots", mode: 'copy', pattern: "*_effect.plot.png"
    publishDir "${params.out}/Plots/ManhattanPlots", mode: 'copy', pattern: "*_manhattan.plot.png"

    input:
    set file(geno), file(pheno), file(aggregate_mapping), file(aggregate_regions) from plotting_data

    output:
    file("*.png") into plots

    """

    Rscript --vanilla `which pipeline.plotting.mod.R` ${aggregate_mapping} `which sweep_summary.tsv`

    """
}




    /*
    set file("*LMM_EXACT_INBRED_qtl_region.tsv"), file("*LMM_EXACT_LOCO_qtl_region.tsv"), file("*AGGREGATE_qtl_region.tsv") into gcta_qtl_to_ld
    Rscript --vanilla `which Find_GCTA_Intervals_Maps.R` ${geno} ${pheno} ${lmmexact_inbred} ${tests} ${qtl_grouping_size} ${qtl_ci_size} ${sig_thresh} ${TRAIT}_LMM_EXACT_INBRED
    Rscript --vanilla `which Find_GCTA_Intervals_Maps_LOCO.R` ${geno} ${pheno} ${lmmexact_loco} ${tests} ${qtl_grouping_size} ${qtl_ci_size} ${sig_thresh} ${TRAIT}_LMM_EXACT_LOCO
    publishDir "${params.out}/Mapping/QTL_Regions", mode: 'copy', pattern: "*LMM_EXACT_INBRED_qtl_region.tsv"
    publishDir "${params.out}/Mapping/QTL_Regions", mode: 'copy', pattern: "*LMM_EXACT_LOCO_qtl_region.tsv"
    ^^ when we're ready for fine mapping ^^
    */


    /*
    ------------ EMMA
    */
/*
    process rrblup_maps {

        cpus 4

        tag { TRAIT }

        publishDir "${params.out}/Mapping/EMMA/Data", mode: 'copy', pattern: "*processed_mapping.tsv"

        input:
        set file("independent_snvs.csv"), file(geno), val(TRAIT), file(pheno), val(P3D), val(sig_thresh), val(qtl_grouping_size), val(qtl_ci_size) from mapping_data_emma

        output:
        file("*raw_mapping.tsv") into raw_map
        set val(TRAIT), file(geno), file(pheno) into processed_map_to_ld
        file("*processed_mapping.tsv") into processed_map_to_summary_plot
        set val(TRAIT), file("*processed_mapping.tsv") into pr_maps_trait

        """

        tests=`cat independent_snvs.csv | grep -v inde`

        Rscript --vanilla `which Run_Mappings.R` ${geno} ${pheno} ${task.cpus} ${P3D} \$tests ${sig_thresh} ${qtl_grouping_size} ${qtl_ci_size}

        if [ -e Rplots.pdf ]; then
        rm Rplots.pdf
        fi

        """
    }
*/
}
