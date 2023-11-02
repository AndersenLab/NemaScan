#load quest version of bcftools 
module load bcftools

#rename chromosomes
rename_file=input_data/all_species/rename_chromosomes
vcf=input_data/c_elegans/genotypes/c_elegans.test.vcf.gz

bcftools annotate -O z --rename-chrs ${rename_file} ${vcf} > input_data/c_elegans/genotypes/c_elegans.test.rename.vcf.gz
