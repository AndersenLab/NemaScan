#! usr/bin/env nextflow

if( !nextflow.version.matches('20.0+') ) {
    println "This workflow requires Nextflow version 20.0 or greater -- You are running version $nextflow.version"
    println "On QUEST, you can use `module load python/anaconda3.6; source activate /projects/b1059/software/conda_envs/nf20_env`"
    exit 1
}

nextflow.preview.dsl=2

date = new Date().format( 'yyyyMMdd' )

/*
~ ~ ~ > * Parameters: common to all analyses
*/
//params.traitfile   = null
//params.vcf         = null
params.help        = null
if(params.simulate) {
    params.e_mem   = "100"
} else {
    params.e_mem   = "10" // I noticed it was crashing with 100 gb for mappings... maybe too much allocation?
}
params.eigen_mem   = params.e_mem + " GB"
//params.R_libpath   = "/projects/b1059/software/R_lib_3.6.0"
params.out         = "Analysis_Results-${date}"
params.debug       = null
params.species     = "elegans"
params.wbb         = "WS276"
params.data_dir    = "input_data/${params.species}"
params.numeric_chrom = "input_data/all_species/rename_chromosomes"
params.sparse_cut  = 0.01
params.group_qtl   = 1000
params.ci_size     = 150
params.sthresh     = "BF"
params.p3d         = "TRUE"
params.maf         = 0.05

if(params.debug) {
    println """

        *** Using debug mode ***

    """
    // debug for now with small vcf
    params.vcf = "330_TEST.vcf.gz"
    vcf = Channel.fromPath("${workflow.projectDir}/test_data/330_TEST.vcf.gz")
    vcf_index = Channel.fromPath("${workflow.projectDir}/test_data/330_TEST.vcf.gz.tbi")
    params.traitfile = "${workflow.projectDir}/test_data/example_trait.tsv"
} else { // does this work with gcp config? which takes preference?
    vcf = Channel.fromPath("/projects/b1059/analysis/WI-${params.vcf}/isotype_only/WI.${params.vcf}.hard-filter.isotype.vcf.gz")
    vcf_index = Channel.fromPath("/projects/b1059/analysis/WI-${params.vcf}/isotype_only/WI.${params.vcf}.hard-filter.isotype.vcf.gz.tbi")
    impute_vcf = Channel.fromPath("/projects/b1059/analysis/WI-${params.vcf}/imputed/WI.${params.vcf}.impute.isotype.vcf.gz")
}


/*
~ ~ ~ > * Parameters: for burden mapping
*/
params.refflat   = "${params.data_dir}/annotations/c_${params.species}_${params.wbb}_refFlat.txt"
params.freqUpper = 0.05
params.minburden = 2


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
log.info "Trait File                              = ${params.traitfile}"
log.info "VCF                                     = ${params.vcf}"
log.info "Significance Threshold                  = ${params.sthresh}"
log.info "Result Directory                        = ${params.out}"
log.info ""
}

// add more params to show user ^^ how does this work with different profiles?


/*
~ ~ ~ > * WORKFLOW
*/
workflow {

    // for mapping
    if(params.maps) {

        // Fix strain names
        Channel.fromPath("${params.traitfile}")
            .combine(Channel.fromPath("${workflow.projectDir}/${params.data_dir}/isotypes/strain_isotype_lookup.tsv")) | fix_strain_names_bulk
        traits_to_map = fix_strain_names_bulk.out.fixed_strain_phenotypes
                .flatten()
                .map { file -> tuple(file.baseName.replaceAll(/pr_/,""), file) }

        // Genotype matrix
        pheno_strains = fix_strain_names_bulk.out.phenotyped_strains_to_analyze

        vcf.spread(vcf_index)
                .combine(pheno_strains) | vcf_to_geno_matrix

        // EIGEN
        contigs = Channel.from(["I", "II", "III", "IV", "V", "X"])
        contigs.combine(vcf_to_geno_matrix.out) | chrom_eigen_variants
        chrom_eigen_variants.out.collect() | collect_eigen_variants

        // GWAS mapping
        pheno_strains
            .spread(traits_to_map)
            .spread(vcf.spread(vcf_index))
            .spread(Channel.fromPath("${params.numeric_chrom}")) | prepare_gcta_files | gcta_grm | gcta_lmm_exact_mapping

        // process GWAS mapping
        traits_to_map
            .spread(collect_eigen_variants.out)
            .spread(vcf_to_geno_matrix.out)
            .combine(Channel.from("${params.p3d}"))
            .combine(Channel.from("${params.sthresh}"))
            .combine(Channel.from("${params.group_qtl}"))
            .combine(Channel.from("${params.ci_size}"))
            .join(gcta_lmm_exact_mapping.out) | gcta_intervals_maps

        // plot
        gcta_intervals_maps.out.maps_to_plot
            .combine(Channel.fromPath("${workflow.projectDir}/bin/sweep_summary.tsv")) | generate_plots 

        // LD b/w regions
        //gcta_intervals_maps.out.maps_to_plot | LD_between_regions

        // summarize all peaks
        peaks = gcta_intervals_maps.out.qtl_peaks
            .collectFile(keepHeader: true, name: "QTL_peaks.tsv", storeDir: "${params.out}/Mapping/Processed")

        // prep LD files
        peaks
            .splitCsv(sep: '\t', skip: 1)
            .join(gcta_intervals_maps.out.maps_to_plot, by: 2)
            .spread(vcf.spread(vcf_index))
            .spread(pheno_strains) //| prep_ld_files

        // divergent regions and haplotypes
        peaks
            .combine(Channel.fromPath("${workflow.projectDir}/${params.data_dir}/isotypes/divergent_bins.bed"))
            .combine(Channel.fromPath("${workflow.projectDir}/${params.data_dir}/isotypes/divergent_df_isotype.bed"))
            .combine(Channel.fromPath("${workflow.projectDir}/${params.data_dir}/isotypes/haplotype_df_isotype.bed"))
            .combine(Channel.fromPath("${workflow.projectDir}/${params.data_dir}/isotypes/div_isotype_list.txt")) | divergent_and_haplotype

        // generate main html report
        peaks
            .spread(traits_to_map)
            .combine(divergent_and_haplotype.out.div_done) | html_report_main

    } else if(params.annotate) {

        // what does annotate do?? just this one process?
        save_dir = "${params.input_data}/${params.species}/annotations"

        Channel.fromPath("${params.script_loc}")
            .combine(save_dir) | update_annotations

    } else if(params.simulate) {

        // vcf = Channel.fromPath("input_data/elegans/genotypes/WI.20180527.impute.vcf.gz")
        // vcf_index = Channel.fromPath("input_data/elegans/genotypes/WI.20180527.impute.vcf.gz.tbi")

        // for simulations
        File pop_file = new File(params.simulate_strains);

        Channel.from(pop_file.collect { it.tokenize( ' ' ) })
            .map { SM, STRAINS -> [SM, STRAINS] }
            .spread(vcf.spread(vcf_index))
            .spread(Channel.fromPath("${params.numeric_chrom}"))
            .spread(Channel.fromPath("${params.simulate_maf}").splitCsv()) | prepare_simulation_files

        // eigen
        contigs = Channel.from(["1", "2", "3", "4", "5", "6"])
        contigs.combine(prepare_simulation_files.out.sim_geno) | chrom_eigen_variants_sims

        chrom_eigen_variants_sims.out.sim_geno_eigen_join
            .groupTuple(by:[0,1,2]).
            join(chrom_eigen_variants_sims.out.sim_geno_meta, by:[0,1,2]) | collect_eigen_variants_sims

        // simulate qtl or genome
        if(params.simulate_qtlloc){

            collect_eigen_variants_sims.out
                .spread(Channel.fromPath("${params.simulate_nqtl}").splitCsv())
                .spread(Channel.fromPath("${params.simulate_qtlloc}"))
                .spread(Channel.fromPath("${params.simulate_eff}").splitCsv())
                .combine(Channel.from(1..params.simulate_reps)) | simulate_effects_loc

            sim_phen_inputs = simulate_effects_loc.out

        } else {

            collect_eigen_variants_sims.out
                .spread(Channel.fromPath("${params.simulate_nqtl}").splitCsv())
                .spread(Channel.fromPath("${params.simulate_eff}").splitCsv())
                .combine(Channel.from(1..params.simulate_reps)) | simulate_effects_genome

            sim_phen_inputs = simulate_effects_genome.out

        }

        sim_phen_inputs
            .spread(Channel.fromPath("${params.simulate_h2}").splitCsv()) | simulate_map_phenotypes

        // simulation mappings
        simulate_map_phenotypes.out.gcta_intervals
            .spread(Channel.from("${params.sthresh}"))
            .spread(Channel.from("${params.group_qtl}"))
            .spread(Channel.from("${params.ci_size}")) | get_gcta_intervals
    }

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

    //executor 'local'

    publishDir "${save_dir}", mode: 'copy'

    input:
        tuple val(gtf_to_refflat), val(save_dir)

    output:
        tuple file("*canonical_geneset.gtf.gz"), file("c_${params.species}_${params.wb_build}_refFlat.txt")

    """
        # add R_libpath to .libPaths() into the R script, create a copy into the NF working directory 
        echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/update_annotations.R > update_annotations.R
        Rscript --vanilla update_annotations.R ${params.wb_build} ${params.species} ${gtf_to_refflat}
    """

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


process fix_strain_names_bulk {

    //executor 'local'

    tag {"BULK TRAIT"}

    publishDir "${params.out}/Phenotypes", mode: 'copy', pattern: "*pr_*.tsv"
    publishDir "${params.out}/Phenotypes", mode: 'copy', pattern: "strain_issues.txt"

    input:
        tuple file(phenotypes), file(isotype_lookup)

    output:
        path "pr_*.tsv", emit: fixed_strain_phenotypes 
        path "Phenotyped_Strains.txt", emit: phenotyped_strains_to_analyze 
        file("strain_issues.txt")

    """
        # add R_libpath to .libPaths() into the R script, create a copy into the NF working directory 
        echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - `which Fix_Isotype_names_bulk.R` > Fix_Isotype_names_bulk.R 

        Rscript --vanilla Fix_Isotype_names_bulk.R ${phenotypes} fix $isotype_lookup
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

    //executor 'local'

    publishDir "${params.out}/Genotype_Matrix", mode: 'copy'

    cpus 1

    input:
        tuple file(vcf), file(index), file(strains)

    output:
        file("Genotype_Matrix.tsv") 

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

    cpus 6
    memory params.eigen_mem

    input:
        tuple val(CHROM), file(genotypes)


    output:
        file("${CHROM}_independent_snvs.csv")


    """
        cat Genotype_Matrix.tsv |\\
        awk -v chrom="${CHROM}" '{if(\$1 == "CHROM" || \$1 == chrom) print}' > ${CHROM}_gm.tsv

        echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - `which Get_GenoMatrix_Eigen.R` > Get_GenoMatrix_Eigen.R
        Rscript --vanilla Get_GenoMatrix_Eigen.R ${CHROM}_gm.tsv ${CHROM}
    """

}

/*
------------ Sum independent tests for all chromosomes
*/

process collect_eigen_variants {

    //executor 'local'

    publishDir "${params.out}/Genotype_Matrix", mode: 'copy'

    cpus 1

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

    cpus 4

    input:
        tuple file(strains), val(TRAIT), file(traits), file(vcf), file(index), file(num_chroms)

    output:
        tuple val(TRAIT), file("plink_formated_trats.tsv"), file("${TRAIT}.bed"), file("${TRAIT}.bim"), file("${TRAIT}.fam"), file("${TRAIT}.map"), file("${TRAIT}.nosex"), file("${TRAIT}.ped"), file("${TRAIT}.log")

    """

    bcftools annotate --rename-chrs ${num_chroms} ${vcf} |\\
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
        tuple val(TRAIT), file(traits), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log)

    output:
        tuple val(TRAIT), file(traits), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file("${TRAIT}_gcta_grm.grm.bin"), file("${TRAIT}_gcta_grm.grm.id"), file("${TRAIT}_gcta_grm.grm.N.bin"), file("${TRAIT}_heritability.hsq"), file("${TRAIT}_heritability.log"), file("${TRAIT}_gcta_grm_inbred.grm.bin"), file("${TRAIT}_gcta_grm_inbred.grm.id"), file("${TRAIT}_gcta_grm_inbred.grm.N.bin"), file("${TRAIT}_heritability_inbred.hsq"), file("${TRAIT}_heritability_inbred.log")

    when:
        params.maps

    """

    gcta64 --bfile ${TRAIT} --autosome --maf 0.05 --make-grm --out ${TRAIT}_gcta_grm --thread-num 10
    gcta64 --bfile ${TRAIT} --autosome --maf 0.05 --make-grm-inbred --out ${TRAIT}_gcta_grm_inbred --thread-num 10

    gcta64 --grm ${TRAIT}_gcta_grm --pheno plink_formated_trats.tsv --reml --out ${TRAIT}_heritability --thread-num 10
    gcta64 --grm ${TRAIT}_gcta_grm_inbred --pheno plink_formated_trats.tsv --reml --out ${TRAIT}_heritability_inbred --thread-num 10

    """
}


process gcta_lmm_exact_mapping {

    cpus 4

    publishDir "${params.out}/Mapping/Raw", pattern: "*fastGWA", overwrite: true
    publishDir "${params.out}/Mapping/Raw", pattern: "*loco.mlma", overwrite: true

    errorStrategy 'ignore'

    input:
    tuple val(TRAIT), file(traits), file(bed), file(bim), file(fam), file(map), \
    file(nosex), file(ped), file(log), file(grm_bin), file(grm_id), file(grm_nbin), \
    file(h2), file(h2log), file(grm_bin_inbred), file(grm_id_inbred), file(grm_nbin_inbred), \
    file(h2_inbred), file(h2log_inbred)

    output:
    tuple val(TRAIT), file("${TRAIT}_lmm-exact_inbred.fastGWA"), file("${TRAIT}_lmm-exact.loco.mlma")


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



process gcta_intervals_maps {

    publishDir "${params.out}/Mapping/Processed", mode: 'copy', pattern: "*AGGREGATE_mapping.tsv"
    publishDir "${params.out}/Mapping/Processed", mode: 'copy', pattern: "*AGGREGATE_qtl_region.tsv" //would be nice to put all these files per trait into one file

    memory '48 GB'

    input:
        tuple val(TRAIT), file(pheno), file(tests), file(geno), val(P3D), val(sig_thresh), \
        val(qtl_grouping_size), val(qtl_ci_size), file(lmmexact_inbred), file(lmmexact_loco)

    output:
        tuple file(geno), file(pheno), val(TRAIT), file("*AGGREGATE_mapping.tsv"), emit: maps_to_plot
        path "*AGGREGATE_qtl_region.tsv", emit: qtl_peaks

    """
    echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - `which Aggregate_Mappings.R` > Aggregate_Mappings.R
    Rscript --vanilla Aggregate_Mappings.R ${lmmexact_loco} ${lmmexact_inbred}

    echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - `which Find_Aggregate_Intervals_Maps.R` > Find_Aggregate_Intervals_Maps.R
    Rscript --vanilla Find_Aggregate_Intervals_Maps.R ${geno} ${pheno} temp.aggregate.mapping.tsv ${tests} ${qtl_grouping_size} ${qtl_ci_size} ${sig_thresh} ${TRAIT}_AGGREGATE

    """
}


process generate_plots {


    publishDir "${params.out}/Plots/LDPlots", mode: 'copy', pattern: "*_LD.plot.png"
    publishDir "${params.out}/Plots/EffectPlots", mode: 'copy', pattern: "*_effect.plot.png"
    publishDir "${params.out}/Plots/ManhattanPlots", mode: 'copy', pattern: "*_manhattan.plot.png"

    input:
        tuple file(geno), file(pheno), val(TRAIT), file(aggregate_mapping), file(sweep_summary)

    output:
        file("*.png")

    """
    echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - `which pipeline.plotting.mod.R` > pipeline.plotting.mod.R
    Rscript --vanilla pipeline.plotting.mod.R ${aggregate_mapping} ${sweep_summary}

    """
}


process LD_between_regions {

  publishDir "${params.out}/Mapping/Processed", mode: 'copy', pattern: "*LD_between_QTL_regions.tsv"

  input:
        tuple file(geno), file(pheno), val(TRAIT), file(aggregate_mapping)

  output:
        tuple val(TRAIT), path("*LD_between_QTL_regions.tsv") optional true
        val TRAIT, emit: linkage_done

  """
    echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - `which LD_between_regions.R` > LD_between_regions.R 
    Rscript --vanilla LD_between_regions.R ${geno} ${aggregate_mapping} ${TRAIT}
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

    input:
        tuple val(TRAIT), val(CHROM), val(marker), val(start_pos), val(peak_pos), val(end_pos), \
        val(peak_id), val(h2), file(geno), file(pheno), file(aggregate_mapping), file(vcf), file(index), file(phenotype)

    output:
        tuple val(TRAIT), val(CHROM), val(marker), val(start_pos), val(peak_pos), val(end_pos), \
        val(peak_id), val(h2), file(geno), file(pheno), file(aggregate_mapping), file(vcf), file(index), \
        file(strains), path("*ROI_Genotype_Matrix.tsv"), path("*LD.tsv") 

    """
        echo "HELLO"
        cat ${aggregate_mapping} |\\
        awk '\$0 !~ "\\tNA\\t" {print}' |\\
        awk '!seen[\$2,\$5,\$12,\$13,\$14]++' |\\
        awk 'NR>1{print \$2, \$5, \$12, \$13, \$14}' OFS="\\t" > ${TRAIT}_QTL_peaks.tsv
        filename='${TRAIT}_QTL_peaks.tsv'
        echo Start
        while read p; do 
            chromosome=`echo \$p | cut -f1 -d' '`
            trait=`echo \$p | cut -f2 -d' '`
            start_pos=`echo \$p | cut -f3 -d' '`
            peak_pos=`echo \$p | cut -f4 -d' '`
            end_pos=`echo \$p | cut -f5 -d' '`
        
        cat ${phenotype} | awk '\$0 !~ "strain" {print}' | cut -f1 > phenotyped_samples.txt
        bcftools view --regions \$chromosome:\$start_pos-\$end_pos ${vcf} \
        -S phenotyped_samples.txt |\\
        bcftools filter -i N_MISSING=0 |\\
        awk '\$0 !~ "#" {print \$1":"\$2}' > \$trait.\$chromosome.\$start_pos.\$end_pos.txt
        bcftools view --regions \$chromosome:\$start_pos-\$end_pos ${vcf} \
        -S phenotyped_samples.txt |\\
        bcftools filter -i N_MISSING=0 -Oz -o finemap.vcf.gz
        plink --vcf finemap.vcf.gz \\
            --snps-only \\
            --maf 0.05 \\
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
        sed 's/ */\\t/g' \$trait.\$chromosome:\$start_pos-\$end_pos.QTL.ld |\\
        cut -f2-10 |\\
        sed 's/^23/X/g' | sed 's/\\t23\\t/\\tX\\t/g' > \$trait.\$chromosome.\$start_pos.\$end_pos.LD.tsv
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




/*
------ Slice out the QTL region for plotting divergent region and haplotype data.
*/


process divergent_and_haplotype {

  //executor 'local'

  publishDir "${params.out}/Divergent_and_haplotype", mode: 'copy'


  input:
    tuple file("QTL_peaks.tsv"), file("divergent_bins"), file(divergent_df_isotype), file(haplotype_df_isotype), file(div_isotype_list)

  output:
    tuple file("all_QTL_bins.bed"), file("all_QTL_div.bed"), file("haplotype_in_QTL_region.txt"), file("div_isotype_list.txt") //, emit: div_hap_table
    val true, emit: div_done


  """
  awk NR\\>1 QTL_peaks.tsv | awk -v OFS='\t' '{print \$1,\$4,\$6}' > QTL_region.bed

  bedtools intersect -wa -a ${divergent_bins} -b QTL_region.bed | sort -k1,1 -k2,2n | uniq > all_QTL_bins.bed

  bedtools intersect -a ${divergent_df_isotype} -b QTL_region.bed | sort -k1,1 -k2,2n | uniq > all_QTL_div.bed

  bedtools intersect -a ${haplotype_df_isotype} -b QTL_region.bed -wo | sort -k1,1 -k2,2n | uniq > haplotype_in_QTL_region.txt

  cp ${div_isotype_list} . 

  """

}

// generate trait-specific html reports
process html_report_main {

  //executor 'local'
  errorStrategy 'ignore'

  tag {TRAIT}
  memory '16 GB'
  

  publishDir "${params.out}/Reports", mode: 'copy'


  input:
    tuple file("QTL_peaks.tsv"), val(TRAIT), file(pheno), val(div_done)

  output:
    tuple file("NemaScan_Report_*.Rmd"), file("NemaScan_Report_*.html")


  """
    cat `which NemaScan_Report_main.Rmd` | sed "s/TRAIT_NAME_HOLDER/${TRAIT}/g" > NemaScan_Report_${TRAIT}_main.Rmd 
    cat `which NemaScan_Report_region_template.Rmd` > NemaScan_Report_region_template.Rmd 

    echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" > .Rprofile

    # probably need to change root dir...
    Rscript -e "rmarkdown::render('NemaScan_Report_${TRAIT}_main.Rmd', knit_root_dir='${workflow.launchDir}/${params.out}')"

  """
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


process prepare_simulation_files {

    cpus 4

    input:
        tuple val(strain_set), val(strains), file(vcf), file(index), file(num_chroms), val(MAF)

    output:
        tuple val(strain_set), val(strains), file("TO_SIMS.bed"), file("TO_SIMS.bim"), file("TO_SIMS.fam"), file("TO_SIMS.map"), file("TO_SIMS.nosex"), file("TO_SIMS.ped"), file("TO_SIMS.log"), file("${strain_set}_${MAF}_Genotype_Matrix.tsv"), val(MAF), emit: sim_geno
        tuple val(strain_set), val(strains), val(MAF), file("renamed_chroms.vcf.gz"), file("renamed_chroms.vcf.gz.tbi"), emit: renamed_chrom_vcf_to_ld


    """

    bcftools annotate --rename-chrs num_chroms ${vcf} |\\
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
------------ Decomposition per chromosome
*/

process chrom_eigen_variants_sims {

    tag { CHROM }

    cpus 6
    memory params.eigen_mem

    input:
        tuple val(CHROM), val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(sex), file(ped), file(log), file(geno), val(MAF)

    output:
        tuple val(strain_set), val(strains), val(MAF), file(bed), file(bim), file(fam), file(map), file(sex), file(ped), file(log), file(geno), emit: sim_geno_meta
        tuple val(strain_set), val(strains), val(MAF), file("${CHROM}_${strain_set}_${MAF}_independent_snvs.csv"), emit: sim_geno_eigen_join


    """
        cat ${geno} |\\
        awk -v chrom="${CHROM}" '{if(\$1 == "CHROM" || \$1 == chrom) print}' > ${CHROM}_gm.tsv

        # add R_libpath to .libPaths() into the R script, create a copy into the NF working directory 
        echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/Get_GenoMatrix_Eigen.R > Get_GenoMatrix_Eigen.R
        Rscript --vanilla Get_GenoMatrix_Eigen.R ${CHROM}_gm.tsv ${CHROM}

        mv ${CHROM}_independent_snvs.csv ${CHROM}_${strain_set}_${MAF}_independent_snvs.csv
    """

}

/*
------------ Sum independent tests for all chromosomes
*/

process collect_eigen_variants_sims {

    //executor 'local'

    publishDir "${params.out}/Genotype_Matrix", mode: 'copy'

    cpus 1

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

    cpus 4

    input:
        tuple val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(MAF), file(n_indep_tests), val(NQTL), file(qtl_loc_bed), val(effect_range), val(SIMREP)

    output:
        tuple val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(MAF), file(n_indep_tests), val(NQTL), val(SIMREP), val(effect_range), file("causal.variants.sim.${NQTL}.${SIMREP}.txt")

    """
        # add R_libpath to .libPaths() into the R script, create a copy into the NF working directory 
        echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/create_causal_QTLs.R > create_causal_QTLs.R
        Rscript --vanilla create_causal_QTLs.R ${bim} ${NQTL} ${effect_range} ${qtl_loc_bed}

        mv causal.variants.sim.${NQTL}.txt causal.variants.sim.${NQTL}.${SIMREP}.txt

    """
}


process simulate_effects_genome {

    tag {NQTL}

    cpus 4

    input:
        tuple val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(MAF), file(n_indep_tests), val(NQTL), val(effect_range), val(SIMREP)

    output:
        tuple val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(MAF), file(n_indep_tests), val(NQTL), val(SIMREP), val(effect_range), file("causal.variants.sim.${NQTL}.${SIMREP}.txt")


    """
        # add R_libpath to .libPaths() into the R script, create a copy into the NF working directory 
        echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/create_causal_QTLs.R > create_causal_QTLs.R
        Rscript --vanilla create_causal_QTLs.R ${bim} ${NQTL} ${effect_range}

        mv causal.variants.sim.${NQTL}.txt causal.variants.sim.${NQTL}.${SIMREP}.txt

    """
}


process simulate_map_phenotypes {

    tag {"${NQTL} - ${SIMREP} - ${H2} - ${MAF}"}

    publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Mappings", pattern: "*fastGWA", overwrite: true
    publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Mappings", pattern: "*loco.mlma", overwrite: true
    publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Phenotypes", pattern: "*.phen", overwrite: true
    publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Phenotypes", pattern: "*.par", overwrite: true

    cpus 4

    errorStrategy 'ignore'

    input:
        tuple val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(MAF), file(n_indep_tests), val(NQTL), val(SIMREP), val(effect_range), file(loci), val(H2)

    output:
        tuple file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.bed"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.bim"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.fam"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.map"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.nosex"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.ped"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.log"), val(NQTL), val(SIMREP), file(loci), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.par"), emit: sim_phen_output
        tuple file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred.fastGWA"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact.loco.mlma"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact.log"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred.log"), emit: sim_GCTA_mapping_results
        path "${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred.fastGWA", emit: lmm_exact_inbred_analyze_sims
        path "${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact.loco.mlma", emit: lmm_exact_loco_analyze_sims
        path "${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen", emit: simphen_analyze_sims
        path "${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.par", emit: simgen_analyze_sims
        tuple val(strain_set), val(strains), val(NQTL), val(SIMREP), val(H2), file(loci), file(gm), val(effect_range), file(n_indep_tests), val(MAF), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred.fastGWA"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact.loco.mlma"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen"), emit: gcta_intervals

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


process get_gcta_intervals {

    publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Mappings", mode: 'copy', pattern: "*processed_aggregate_mapping.tsv"
    publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Mappings", mode: 'copy', pattern: "*processed_LMM-EXACT-INBRED_mapping.tsv"
    publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Mappings", mode: 'copy', pattern: "*processed_LMM-EXACT-LOCO_mapping.tsv"
    publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Mappings", mode: 'copy', pattern: "*qtl_region.tsv"

    memory '48 GB'

    input:
    tuple val(strain_set), val(strains), val(NQTL), val(SIMREP), val(H2), file(loci), file(gm), val(effect_range), file(n_indep_tests), val(MAF), file(lmmexact_inbred), file(lmmexact_loco), file(phenotypes), val(THRESHOLD), val(QTL_GROUP_SIZE), val(QTL_CI_SIZE)

    output:
    tuple val(strain_set), val(strains), val(NQTL), val(SIMREP), val(H2), file(loci), file(gm), val(effect_range), file(n_indep_tests), file(phenotypes), val(THRESHOLD), file("*processed_aggregate_mapping.tsv"), file("*processed_LMM-EXACT-INBRED_mapping.tsv"), file("*processed_LMM-EXACT-LOCO_mapping.tsv"), emit: processed_gcta
    tuple val(strain_set), val(strains), val(MAF), val(NQTL), val(SIMREP), val(H2), val(effect_range), file("*aggregate_qtl_region.tsv"), emit: gcta_qtl_to_ld
    tuple val(strain_set), val(strains), val(MAF), val(NQTL), val(SIMREP), val(H2), val(effect_range), file(loci), file(phenotypes), emit: simulated_phenotypes

    """
        # add R_libpath to .libPaths() into the R script, create a copy into the NF working directory 
        echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/Aggregate_Mappings.R > Aggregate_Mappings.R
        Rscript --vanilla Aggregate_Mappings.R ${lmmexact_loco} ${lmmexact_inbred}

        echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/Find_Aggregate_Intervals.R > Find_Aggregate_Intervals.R
        Rscript --vanilla Find_Aggregate_Intervals.R ${gm} ${phenotypes} temp.aggregate.mapping.tsv ${n_indep_tests} ${NQTL} ${SIMREP} ${QTL_GROUP_SIZE} ${QTL_CI_SIZE} ${H2} ${params.maf} ${THRESHOLD} ${strain_set} ${MAF} ${effect_range} aggregate
        
        echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/Find_GCTA_Intervals.R > Find_GCTA_Intervals.R
        Rscript --vanilla Find_GCTA_Intervals.R ${gm} ${phenotypes} ${lmmexact_inbred} ${n_indep_tests} ${NQTL} ${SIMREP} ${QTL_GROUP_SIZE} ${QTL_CI_SIZE} ${H2} ${params.maf} ${THRESHOLD} ${strain_set} ${MAF} ${effect_range} LMM-EXACT-INBRED
        
        echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/Find_GCTA_Intervals_LOCO.R > Find_GCTA_Intervals_LOCO.R
        Rscript --vanilla Find_GCTA_Intervals_LOCO.R ${gm} ${phenotypes} ${lmmexact_loco} ${n_indep_tests} ${NQTL} ${SIMREP} ${QTL_GROUP_SIZE} ${QTL_CI_SIZE} ${H2} ${params.maf} ${THRESHOLD} ${strain_set} ${MAF} ${effect_range} LMM-EXACT-LOCO

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

/*
=====================================
~ > *                           * < ~
~ ~ > *                       * < ~ ~
~ ~ ~ > *  GENERATE REPORT  * < ~ ~ ~
~ ~ > *                       * < ~ ~
~ > *                           * < ~
=====================================
*/

workflow.onComplete {

    summary = """

    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]

    { Parameters }
    ---------------------------
    Phenotype File                          = ${params.traitfile}
    VCF                                     = ${params.vcf}

    Significance Threshold                  = ${params.sthresh}
    P3D                                     = ${params.p3d}
    Max AF for Burden Mapping               = ${params.freqUpper}
    Min Strains with Variant for Burden     = ${params.minburden}
    Threshold for grouping QTL              = ${params.group_qtl}
    Number of SNVs to define CI             = ${params.ci_size}
    Eigen Memory allocation                 = ${params.eigen_mem}
    Path to R libraries.                    = ${params.R_libpath}

    Mapping                                 = ${params.maps}
    Simulation                              = ${params.simulate}
    Simulate QTL effects                    = ${params.simulate_qtlloc}
    Annotation                              = ${params.annotate}
    Result Directory                        = ${params.out}

    """

    println summary

    def outlog = new File("${params.out}/log.txt")
    outlog.newWriter().withWriter {
        outlog << param_summary
        outlog << summary
    }
}

    // mail summary
    //if (params.email) {
    //    ['mail', '-s', 'cegwas2-nf', params.email].execute() << summary
    //}




