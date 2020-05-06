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
params.fix_names = "fix"

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
}


/*
~ ~ ~ > * Parameters: for burden mapping 
*/
params.refflat   = "${params.data_dir}/annotations/${params.species}_${params.wbb}_refFlat.txt"
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
    log.info "             -profile mappings USAGE"
    log.info "----------------------------------------------------------------" 
    log.info "----------------------------------------------------------------"   
    log.info "nextflow main.nf --vcf input_data/elegans/genotypes/WI.20180527.impute.vcf.gz --traitfile input_data/elegans/phenotypes/PC1.tsv -profile mappings --p3d TRUE" 
    log.info "----------------------------------------------------------------" 
    log.info "----------------------------------------------------------------" 
    log.info "Mandatory arguments:"
    log.info "--traitfile              String                Name of file that contains phenotypes. File should be tab-delimited with the columns: strain trait1 trait2 ..."
    log.info "--vcf                    String                Name of VCF to extract variants from. There should also be a tabix-generated index file with the same name in the directory that contains the VCF. If none is provided, the pipeline will download the latest VCF from CeNDR"
    log.info "--p3d                    BOOLEAN               Set to FALSE for EMMA algortith, TRUE for EMMAx"
    log.info "Optional arguments:"
    log.info "--maf                    String                Minimum minor allele frequency to use for single-marker mapping (Default: 0.05)"
    log.info "--lmm                    String                Perform GCTA mapping with --fastGWA-lmm algorithm (Default: RUN, option to not run is null)"
    log.info "--lmm-exact              String                Perform GCTA mapping with --fastGWA-lmm-exact algorithm (Default: RUN, option to not run is null)"
    log.info "--sparse_cut             String                Any off-diagonal value in the genetic relatedness matrix greater than this is set to 0 (Default: 0.05)"
    log.info "----------------------------------------------------------------"
    log.info "             -profile simulations USAGE"
    log.info "----------------------------------------------------------------"    
    log.info "Mandatory arguments:"
    log.info "--simulate_nqtl          String                A CSV file with the number of QTL to simulate per phenotype, one value per line (Default is located: input_data/all_species/simulate_nqtl.csv)"
    log.info "--simulate_h2            String                A CSV file with phenotype heritability, one value per line (Default is located: input_data/all_species/simulate_h2.csv)"
    log.info "Optional arguments:"
    log.info "--simulate_reps          String                The number of replicates to simulate per number of QTL and heritability (Default: 2)"
    log.info "--simulate_maf           String                Minimum minor allele frequency to use for single-marker mapping (Default: 0.05)"
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
log.info "VCF                                     = ${params.simulate}"
log.info "CeNDR Release                           = ${params.refflat}"
log.info "P3D                                     = ${params.p3d}"
log.info "Significance Threshold                  = ${params.sthresh}"
log.info "Max AF for Burden Mapping               = ${params.freqUpper}"
log.info "Min Strains with Variant for Burden     = ${params.minburden}"
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
          qtl_ci_size_sim}

/*
~ ~ ~ > * INITIATE MAPPING QTL CONFIDENCE INTERVAL SIZE PARAMETER
*/

Channel
    .from("${params.group_qtl}")
    .into{qtl_snv_grouping_maps;
          qtl_snv_grouping_sims}

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
          sig_threshold_fine_sim}

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
           rename_chroms_sims }

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

    process simulate_strain_names {

        executor 'local'

        input:
            set file(vcf), file(index) from vcf_to_names

        output:
            file("sorted_samples.txt") into phenotyped_strains_to_analyze

        when:
            params.simulate

        """
        bcftools query -l ${vcf} |\\
        sort > sorted_samples.txt 
        """
    }


phenotyped_strains_to_analyze
    .into{strain_list_genome;
          strain_list_emma;
          strain_list_simulate;
          strain_list_finemap;
          strain_list_gcta_grml;
          strain_list_gcta_map}

}


/*
===================================================================
~ > *                                                         * < ~
~ ~ > *                                                     * < ~ ~
~ ~ ~ > *  CONVERT THE VCF TO A GENOTYPE MATRIX FOR EMMA  * < ~ ~ ~
~ ~ > *                                                     * < ~ ~
~ > *                                                         * < ~
===================================================================
*/

process vcf_to_geno_matrix {

    executor 'local'

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

    executor 'local'

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


if(params.maps){
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
        .into{gcta_prep_inputs;
              print_inputs}

} else if (params.simulate){
    strain_list_simulate
    .spread(vcf_to_simulations)
    .into{simulation_prep_inputs;
          print_inputs}
}


if(params.simulate){

    process prepare_simulation_files {

        cpus 4

        input:
            file(num_chroms) from rename_chroms_sims
            set file(strains), file(vcf), file(index) from simulation_prep_inputs

        output:
            set file("TO_SIMS.bed"), file("TO_SIMS.bim"), file("TO_SIMS.fam"), file("TO_SIMS.map"), file("TO_SIMS.nosex"), file("TO_SIMS.ped"), file("TO_SIMS.log") into sim_inputs
            file("Genotype_Matrix.tsv") into sim_geno_matrix

        when:
            params.simulate

        """

        bcftools annotate --rename-chrs rename_chromosomes ${vcf} |\\
        bcftools view -S ${strains} |\\
        bcftools filter -i N_MISSING=0 -Oz -o renamed_chroms.vcf.gz

        tabix -p vcf renamed_chroms.vcf.gz

        plink --vcf renamed_chroms.vcf.gz \\
        --snps-only \\
        --biallelic-only \\
        --maf ${params.simulate_maf} \\
        --set-missing-var-ids @:# \\
        --indep-pairwise 50 10 0.8 \\
        --geno \\
        --not-chr MtDNA \\
        --allow-extra-chr

        plink --vcf renamed_chroms.vcf.gz \\
        --make-bed \\
        --snps-only \\
        --biallelic-only \\
        --maf ${params.simulate_maf} \\
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
        sed 's/.\\/./NA/g' > Genotype_Matrix.tsv

        """
    }

    sim_inputs
        .spread(nqtl)
        .set{sim_nqtl_inputs}

    process simulate_effects {

        tag {NQTL}

        cpus 4

        input:
            set file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), val(NQTL) from sim_nqtl_inputs
            each SIMREP from simulation_replicates

        output:
            set file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), val(NQTL), val(SIMREP), file("causal.variants.sim.${NQTL}.${SIMREP}.txt") into sim_phen_inputs

        when:
            params.simulate

        """

         Rscript --vanilla `which create_causal_QTLs.R` ${bim} ${NQTL}

         mv causal.variants.sim.${NQTL}.txt causal.variants.sim.${NQTL}.${SIMREP}.txt

        """
    }

    sim_phen_inputs
        .spread(sim_h2)
        .set{sim_phen_h2_input}

    process simulate_map_phenotypes {

        tag {"${NQTL} - ${SIMREP} - ${H2}"}

        publishDir "${params.out}/Simulations/${NQTL}/Mappings", mode: 'copy', pattern: "*fastGWA"
        publishDir "${params.out}/Simulations/${NQTL}/Phenotypes", mode: 'copy', pattern: "${NQTL}_${SIMREP}_${H2}_sims.phen"
        publishDir "${params.out}/Simulations/${NQTL}/Phenotypes", mode: 'copy', pattern: "${NQTL}_${SIMREP}_${H2}_sims.par"

        cpus 4

        input:
            set file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), val(NQTL), val(SIMREP), file(loci), val(H2) from sim_phen_h2_input

        output:
            set file("TO_SIMS_${NQTL}_${SIMREP}.bed"), file("TO_SIMS_${NQTL}_${SIMREP}.bim"), file("TO_SIMS_${NQTL}_${SIMREP}.fam"), file("TO_SIMS_${NQTL}_${SIMREP}.map"), file("TO_SIMS_${NQTL}_${SIMREP}.nosex"), file("TO_SIMS_${NQTL}_${SIMREP}.ped"), file("TO_SIMS_${NQTL}_${SIMREP}.log"), val(NQTL), val(SIMREP), file(loci), file("${NQTL}_${SIMREP}_${H2}_sims.phen"), file("${NQTL}_${SIMREP}_${H2}_sims.par") into sim_phen_output
            set file("${NQTL}_${SIMREP}_${H2}_lmm-exact.fastGWA"), file("${NQTL}_${SIMREP}_${H2}_lmm-exact_inbred.fastGWA"), file("${NQTL}_${SIMREP}_${H2}_lmm-exact.log"), file("${NQTL}_${SIMREP}_${H2}_lmm-exact_inbred.log") into sim_GCTA_mapping_results
            set val(NQTL), val(SIMREP), val(H2), file(loci), file("${NQTL}_${SIMREP}_${H2}_sims.phen"), file("${NQTL}_${SIMREP}_${H2}_sims.par") into sim_phen_to_emma

        when:
            params.simulate

        """

        gcta64 --bfile TO_SIMS \\
             --simu-qt \\
             --simu-causal-loci ${loci} \\
             --simu-hsq ${H2} \\
             --simu-rep 1 \\
             --out ${NQTL}_${SIMREP}_${H2}_sims

        plink --bfile TO_SIMS \\
            --make-bed \\
            --snps-only \\
            --biallelic-only \\
            --maf ${params.simulate_maf} \\
            --set-missing-var-ids @:# \\
            --geno \\
            --recode \\
            --out TO_SIMS_${NQTL}_${SIMREP} \\
            --allow-extra-chr \\
            --pheno ${NQTL}_${SIMREP}_${H2}_sims.phen

        gcta64 --bfile TO_SIMS_${NQTL}_${SIMREP} --autosome --maf ${params.simulate_maf} --make-grm --out TO_SIMS_${NQTL}_${SIMREP}_${H2}_gcta_grm --thread-num 10
        gcta64 --bfile TO_SIMS_${NQTL}_${SIMREP} --autosome --maf ${params.simulate_maf} --make-grm-inbred --out TO_SIMS_${NQTL}_${SIMREP}_${H2}_gcta_grm_inbred --thread-num 10


        gcta64 --grm TO_SIMS_${NQTL}_${SIMREP}_${H2}_gcta_grm_inbred --pheno ${NQTL}_${SIMREP}_${H2}_sims.phen --reml --out check_vp --thread-num 10

        vp=`grep Vp check_vp.hsq | head -1 | cut -f2`

        if (( \$(bc <<< "\$vp==0") > 0 )); 
        then
        awk '{print \$1, \$2, \$3*1000}' ${NQTL}_${SIMREP}_${H2}_sims.phen > temp.phen;
        rm ${NQTL}_${SIMREP}_${H2}_sims.phen
        mv temp.phen ${NQTL}_${SIMREP}_${H2}_sims.phen
        fi

        gcta64 --grm TO_SIMS_${NQTL}_${SIMREP}_${H2}_gcta_grm --make-bK-sparse ${params.sparse_cut} --out ${NQTL}_${SIMREP}_${H2}_sparse_grm

        gcta64 --fastGWA-lmm-exact \\
            --grm-sparse ${NQTL}_${SIMREP}_${H2}_sparse_grm \\
            --bfile TO_SIMS_${NQTL}_${SIMREP} \\
            --out ${NQTL}_${SIMREP}_${H2}_lmm-exact \\
            --pheno ${NQTL}_${SIMREP}_${H2}_sims.phen \\
            --maf ${params.simulate_maf}

        gcta64 --grm TO_SIMS_${NQTL}_${SIMREP}_${H2}_gcta_grm_inbred --make-bK-sparse ${params.sparse_cut} --out ${NQTL}_${SIMREP}_${H2}_sparse_grm_inbred

        gcta64 --fastGWA-lmm-exact \\
            --grm-sparse ${NQTL}_${SIMREP}_${H2}_sparse_grm_inbred \\
            --bfile TO_SIMS_${NQTL}_${SIMREP} \\
            --out ${NQTL}_${SIMREP}_${H2}_lmm-exact_inbred \\
            --pheno ${NQTL}_${SIMREP}_${H2}_sims.phen \\
            --maf ${params.simulate_maf}

        """
    }

    sim_phen_to_emma
        .spread(sim_geno_matrix)
        .spread(qtl_snv_grouping_sims)
        .spread(qtl_ci_size_sim)
        .spread(p3d_full_sim)
        .spread(sig_threshold_full_sim)
        .into{sim_emma_inputs;
              sim_emma_fine_inputs}

/*
------------ For simulations, i want to separate mapping and defining threshold, so we can look at various threshold after the long mapping step
------------ further optimization can be made by making kinship matrix outside of mapping process
*/

    process sim_emmma_maps {

        cpus 4

        publishDir "${params.out}/Simulations/${NQTL}/Mappings", mode: 'copy', pattern: "*processed_mapping.tsv"

        input:
        set val(NQTL), val(SIMREP), val(H2), file(loci), file(pheno), file(sim_params), file(geno), val(QTL_GROUP_SIZE), val(QTL_CI_SIZE), val(P3D), val(THRESHOLD) from sim_emma_inputs

        output:
        set val(NQTL), val(SIMREP), val(H2), file("*raw_mapping.tsv"), file("*processed_mapping.tsv") into pr_sim_emma_maps

        """

        Rscript --vanilla `which Run_Sims_EMMA_SJW.R` ${geno} ${pheno} ${task.cpus} ${P3D} ${NQTL} ${SIMREP} ${QTL_GROUP_SIZE} ${QTL_CI_SIZE} ${H2} ${params.maf} ${THRESHOLD}
        
        """
    }

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
            file(num_chroms) from rename_chroms_gcta
            set file(strains), val(TRAIT), file(traits), file(vcf), file(index) from gcta_prep_inputs

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

        publishDir "${params.out}/Mapping/lmm_exact/Data", mode: 'copy', pattern: "*_lmm-exact.fastGWA"
        publishDir "${params.out}/Mapping/lmm_exact/Data", mode: 'copy', pattern: "*_lmm-exact_inbred.fastGWA"

        input:
        set val(TRAIT), file(traits), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(grm_bin), file(grm_id), file(grm_nbin), file(h2), file(h2log), file(grm_bin_inbred), file(grm_id_inbred), file(grm_nbin_inbred), file(h2_inbred), file(h2log_inbred) from gcta_lmm_exact

        output:
        set val(TRAIT), file(traits), file("${TRAIT}_lmm-exact.fastGWA"), file("${TRAIT}_lmm-exact_inbred.fastGWA") into lmm_exact_output

        when:
          params.lmm_exact

        """

        gcta64 --grm ${TRAIT}_gcta_grm --make-bK-sparse ${params.sparse_cut} --out ${TRAIT}_sparse_grm

        gcta64 --fastGWA-lmm-exact \\
            --grm-sparse ${TRAIT}_sparse_grm \\
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

    process gcta_lmm_mapping {

        cpus 4

        publishDir "${params.out}/Mapping/lmm/Data", mode: 'copy', pattern: "*_lmm.fastGWA"
        publishDir "${params.out}/Mapping/lmm/Data", mode: 'copy', pattern: "*_lmm_inbred.fastGWA"

        input:
        set val(TRAIT), file(traits), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(grm_bin), file(grm_id), file(grm_nbin), file(h2), file(h2log), file(grm_bin_inbred), file(grm_id_inbred), file(grm_nbin_inbred), file(h2_inbred), file(h2log_inbred) from gcta_lmm

        output:
        set val(TRAIT), file(traits), file("${TRAIT}_lmm.fastGWA"), file("${TRAIT}_lmm_inbred.fastGWA") into lmm_output

        when:
          params.lmm    

        """

        gcta64 --grm ${TRAIT}_gcta_grm --make-bK-sparse ${params.sparse_cut} --out ${TRAIT}_sparse_grm

        gcta64 --fastGWA-lmm \\
            --grm-sparse ${TRAIT}_sparse_grm \\
            --bfile ${TRAIT} \\
            --out ${TRAIT}_lmm \\
            --pheno ${traits} \\
            --maf ${params.maf}

        gcta64 --grm ${TRAIT}_gcta_grm_inbred --make-bK-sparse ${params.sparse_cut} --out ${TRAIT}_sparse_grm_inbred

        gcta64 --fastGWA-lmm \\
            --grm-sparse ${TRAIT}_sparse_grm \\
            --bfile ${TRAIT} \\
            --out ${TRAIT}_lmm_inbred \\
            --pheno ${traits} \\
            --maf ${params.maf}

        """
    }

    /*
    ------------ EMMA
    */

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

}




