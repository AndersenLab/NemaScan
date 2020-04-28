#! usr/bin/env nextflow

date = new Date().format( 'yyyyMMdd' )

/*
~ ~ ~ > * Parameters: common to all analyses 
*/
params.traitfile = null
params.vcf 		 = null
params.sthresh   = null
params.help 	 = null
params.e_mem 	 = "100"
params.eigen_mem = params.e_mem + " GB"
params.cendr_v   = "20180527"
params.fix_names = "fix"


/*
~ ~ ~ > * Parameters: for burden mapping 
*/
params.refflat   = "${params.data_dir}/${params.species}_${params.wbb}_refFlat.txt"
params.freqUpper = 0.05
params.minburden = 2

/*
~ ~ ~ > * Parameters: for EMMA mapping 
*/
params.p3d 		 = null
params.group_qtl = 1000
params.ci_size   = 150

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
    log.info "nextflow main.nf --traitfile=test_bulk --vcf=bin/WI.20180527.impute.vcf.gz --p3d=TRUE --sthresh=EIGEN # run all traits from a single file"
    log.info "nextflow main.nf --traitdir=test_bulk --p3d=TRUE --sthresh=BF # download VCF from CeNDR"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--traitfile              String                Name of file that contains phenotypes. File should be tab-delimited with the columns: strain trait1 trait2 ..."
    log.info "--vcf                    String                Name of VCF to extract variants from. There should also be a tabix-generated index file with the same name in the directory that contains the VCF. If none is provided, the pipeline will download the latest VCF from CeNDR"
    log.info "--p3d                    BOOLEAN               Set to FALSE for EMMA algortith, TRUE for EMMAx"
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
log.info "Phenotype Directory                     = ${params.traitdir}"
log.info "VCF                                     = ${params.vcf}"
log.info "CeNDR Release                           = ${params.cendr_v}"
log.info "P3D                                     = ${params.p3d}"
log.info "Significance Threshold                  = ${params.sthresh}"
log.info "Max AF for Burden Mapping               = ${params.freqUpper}"
log.info "Min Strains with Variant for Burden     = ${params.minburden}"
log.info "Significance Threshold                  = ${params.sthresh}"
log.info "Gene File                               = ${params.genes}"
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
        .into{vcf_to_whole_gcta_grm;
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
        .into{vcf_to_whole_gcta_grm;
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
    .set{qtl_ci_size}

/*
~ ~ ~ > * INITIATE MAPPING QTL CONFIDENCE INTERVAL SIZE PARAMETER
*/

Channel
    .from("${params.group_qtl}")
    .set{qtl_snv_groupinng}

/*
~ ~ ~ > * INITIATE MAPPING METHOD CHANNEL
*/

Channel
    .from("${params.p3d}")
    .into{p3d_full;
          p3d_fine}

/*
~ ~ ~ > * INITIATE THRESHOLD CHANNEL
*/

Channel
    .from("${params.sthresh}")
    .into{sig_threshold_full;
          sig_threshold_fine}

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
    .set{ rename_chroms }

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

    executor 'local'

    input:
        file(phenotypes) from traits_to_strainlist

    output:
        file("pr_*.tsv") into fixed_strain_phenotypes
        file("Phenotyped_Strains.txt") into phenotyped_strains_to_analyze

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
          strain_list_finemap;
          strain_list_gcta_grml;
          strain_list_gcta_map}


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
    .spread(qtl_snv_groupinng)
    .spread(qtl_ci_size)
    .into{mapping_data_emma;
          mapping_data_gcta}

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

strain_list_gcta_grml
    .spread(traits_to_gcta_grm)
    .spread(vcf_to_whole_gcta_grm)
    .set{gcta_prep_inputs}


process prepare_gcta_files {

    cpus 4

    input:
    file(num_chroms) from rename_chroms
    set file(strains), val(TRAIT), file(traits), file(vcf), file(index) from gcta_prep_inputs

    output:
    set val(TRAIT), file(traits), file("${TRAIT}.bed"), file("${TRAIT}.bim"), file("${TRAIT}.fam"), file("${TRAIT}.map"), file("${TRAIT}.nosex"), file("${TRAIT}.ped"), file("${TRAIT}.log") into gcta_grm_inputs

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


    """

    gcta64 --bfile ${TRAIT} --autosome --maf 0.05 --make-grm --out ${TRAIT}_gcta_grm --thread-num 10

    """
}

