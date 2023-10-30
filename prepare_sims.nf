if( !nextflow.version.matches('>20.0') ) {
    println "This workflow requires Nextflow version 20.0 or greater -- You are running version $nextflow.version"
    println "On QUEST, you can use `module load python/anaconda3.6; source activate /projects/b1059/software/conda_envs/nf20_env`"
    exit 1
}

nextflow.preview.dsl=2
// nextflow.enable.dsl=2


date = new Date().format( 'yyyyMMdd' )

params.out = "SimFiles-${date}"

params.maf = 0.05

params.cores = 4

params.ld = false

rename_key = Channel.fromPath("input_data/all_species/rename_chromosomes")
maf_file = Channel.fromPath("input_data/all_species/simulate_maf.csv").splitCsv()

workflow{
    //load the strain set file 
    File pop_file = new File("/projects/b1059/projects/Ryan/ortholog_sims/20230614_pre_sim_ogs/data/08.18.23_all_outliers.txt") ;
    
    //Array used to attach the correct vcf to strain sets 
    sp_ids = [["c_elegans", "/projects/b1059/data/c_elegans/WI/variation/20220216/vcf/WI.20220216.hard-filter.isotype.vcf.gz", "/projects/b1059/data/c_elegans/WI/variation/20220216/vcf/WI.20220216.hard-filter.isotype.vcf.gz.tbi"],
                ["c_briggsae", "/projects/b1059/data/c_briggsae/WI/variation/20210803/vcf/WI.20210803.hard-filter.isotype.vcf.gz", "/projects/b1059/data/c_briggsae/WI/variation/20210803/vcf/WI.20210803.hard-filter.isotype.vcf.gz.tbi"],
                ["c_tropicalis", "/projects/b1059/data/c_tropicalis/WI/variation/20210901/vcf/WI.20210901.hard-filter.isotype.vcf.gz", "/projects/b1059/data/c_tropicalis/WI/variation/20210901/vcf/WI.20210901.hard-filter.isotype.vcf.gz.tbi"]]
if(params.ld == false) {

Channel.from(pop_file.collect { it.tokenize( ' ' ) })
          .map {SP, SM, STRAINS -> [SP, SM, STRAINS] }
          .combine(Channel.from(sp_ids), by:[0])
          .combine(rename_key)
          .combine(maf_file)
          .map {
            tuple ->[tuple[0], //extract species ID
                    (tuple[1]), //extract sample name
                    (tuple[2]), //extract sample list
                    file(tuple[3]), // convert path to file obj for vcf
                    file(tuple[4]), // index
                    file(tuple[5]), // rename key 
                    tuple[6]] // MAF
       } | prepare_simulation_files


}

if(params.ld == true) {

Channel.from(pop_file.collect { it.tokenize( ' ' ) })
          .map {SP, SM, STRAINS -> [SP, SM, STRAINS] }
          .combine(Channel.from(sp_ids), by:[0])
          .combine(rename_key)
          .combine(maf_file)
          .map {
            tuple ->[tuple[0], //extract species ID
                    (tuple[1]), //extract sample name
                    (tuple[2]), //extract sample list
                    file(tuple[3]), // convert path to file obj for vcf
                    file(tuple[4]), // index
                    file(tuple[5]), // rename key 
                    tuple[6]] // MAF
       } | prepare_simulation_files_ld


}

}
process prepare_simulation_files {
    container 'andersenlab/nemascan:20220407173056db3227'
    cpus 4
    memory 30.GB
    time '30m'
    
    publishDir "${params.out}/${sp}/${strain_set}", pattern: "*vcf.gz", overwrite: true
    publishDir "${params.out}/${sp}/${strain_set}", pattern: "*vcf.gz.tbi", overwrite: true
    publishDir "${params.out}/${sp}/${strain_set}", pattern: "*bim", overwrite: true
    publishDir "${params.out}/${sp}/${strain_set}", pattern: "*log", overwrite: true
    publishDir "${params.out}/${sp}/${strain_set}", pattern: "*frq", overwrite: true


    input:
        tuple val(sp), val(strain_set), val(strains), file(vcf), file(index), file(num_chroms), val(MAF)

    output:
        tuple val(strain_set), val(strains), file("${strain_set}_${MAF}.bed"), file("${strain_set}_${MAF}.bim"), file("${strain_set}_${MAF}.fam"), file("${strain_set}_${MAF}.map"), file("${strain_set}_${MAF}.nosex"), file("${strain_set}_${MAF}.ped"), file("${strain_set}_${MAF}.log"), file("${strain_set}_${MAF}_Genotype_Matrix.tsv"),  val(MAF), file("${strain_set}_${MAF}.frq"),  emit: sim_geno
        tuple val(strain_set), val(strains), val(MAF), file("renamed_chroms.vcf.gz"), file("renamed_chroms.vcf.gz.tbi"), emit: renamed_chrom_vcf_to_ld // This output is no longer used. Should be removed in the future.


    """
    bcftools annotate --rename-chrs ${num_chroms} ${vcf} |\\
    bcftools view -s `echo ${strains} | tr -d '\\n'` |\\
    bcftools filter -i N_MISSING=0 -Oz -o renamed_chroms.vcf.gz
    tabix -p vcf renamed_chroms.vcf.gz

   
    plink --vcf renamed_chroms.vcf.gz \\
    --make-bed \\
    --snps-only \\
    --biallelic-only \\
    --maf ${MAF} \\
    --set-missing-var-ids @:# \\
    --geno \\
    --recode \\
    --out ${strain_set}_${MAF} \\
    --allow-extra-chr

    plink --vcf renamed_chroms.vcf.gz \\
    --snps-only \\
    --biallelic-only \\
    --maf ${MAF} \\
    --set-missing-var-ids @:# \\
    --geno \\
    --freq \\
    --out ${strain_set}_${MAF}

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


process prepare_simulation_files_ld {
    container 'andersenlab/nemascan:20220407173056db3227'
    cpus 4
    memory 30.GB
    time '30m'
    
    publishDir "${params.out}/${sp}/${strain_set}", pattern: "*vcf.gz", overwrite: true
    publishDir "${params.out}/${sp}/${strain_set}", pattern: "*vcf.gz.tbi", overwrite: true
    publishDir "${params.out}/${sp}/${strain_set}", pattern: "*bim", overwrite: true
    publishDir "${params.out}/${sp}/${strain_set}", pattern: "*log", overwrite: true
    publishDir "${params.out}/${sp}/${strain_set}", pattern: "*frq", overwrite: true


    input:
        tuple val(sp), val(strain_set), val(strains), file(vcf), file(index), file(num_chroms), val(MAF)

    output:
        tuple val(strain_set), val(strains), file("${strain_set}_${MAF}.bed"), file("${strain_set}_${MAF}.bim"), file("${strain_set}_${MAF}.fam"), file("${strain_set}_${MAF}.map"), file("${strain_set}_${MAF}.nosex"), file("${strain_set}_${MAF}.ped"), file("${strain_set}_${MAF}.log"), file("${strain_set}_${MAF}_Genotype_Matrix.tsv"),  val(MAF), file("${strain_set}_${MAF}.frq"),  emit: sim_geno
        tuple val(strain_set), val(strains), val(MAF), file("renamed_chroms.vcf.gz"), file("renamed_chroms.vcf.gz.tbi"), emit: renamed_chrom_vcf_to_ld // This output is no longer used. Should be removed in the future.


    """
    bcftools annotate --rename-chrs ${num_chroms} ${vcf} |\\
    bcftools view -s `echo ${strains} | tr -d '\\n'` |\\
    bcftools filter -i N_MISSING=0 -Oz -o renamed_chroms.vcf.gz
    tabix -p vcf renamed_chroms.vcf.gz

   
    plink --vcf renamed_chroms.vcf.gz \\
    --make-bed \\
    --snps-only \\
    --biallelic-only \\
    --maf ${MAF} \\
    --set-missing-var-ids @:# \\
    --indep-pairwise 50 10 0.8 \\
    --geno \\
    --recode \\
    --out ${strain_set}_${MAF} \\
    --allow-extra-chr

    plink --vcf renamed_chroms.vcf.gz \\
    --make-bed \\
    --snps-only \\
    --biallelic-only \\
    --maf ${MAF} \\
    --set-missing-var-ids @:# \\
    --extract ${strain_set}_${MAF}.prune.in \\
    --geno \\
    --freq \\
    --allow-extra-chr \\
    --out ${strain_set}_${MAF}

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