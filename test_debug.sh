#!/bin/bash
#SBATCH -J nemascan_test       
#SBATCH -A b1042               
#SBATCH -p genomicsguestA               
#SBATCH -t 12:00:00            
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=16G               
#SBATCH -o nemascan.test.slurm.txt

# script to test all profiles of nemascan to make sure everything is working for a new release

# first: run nextflow run main.nf --debug to make sure the classic debug mode works well
# second: make sure to commit all changes
# third: make new directory (i.e. 20211020_test_debug) and cd into this directory then run `sbatch ../test_debug.sh`


module purge
module load python/anaconda3.6
source activate nf20_env
module add singularity
cp ../test_template.sh .

# step 0: print out git commit
git status > git.log
git log | head -n 6 >> git.log

# step 1: test c_elegans mapping profile
sbatch test_template.sh 'test_mapping_elegans' '--vcf 20220216 --traitfile ../../input_data/c_elegans/phenotypes/test_pheno.tsv --out test_mapping_elegans -resume -N kathrynevans2015@u.northwestern.edu'

# step 1: test c_elegans mapping profile no mediation
sbatch test_template.sh 'test_mapping_elegans_nomediate' '--vcf 20220216 --traitfile ../../input_data/c_elegans/phenotypes/test_pheno.tsv --out test_mapping_elegans_nomediate --mediation false -resume -N kathrynevans2015@u.northwestern.edu'

# step 1: test c_elegans mapping profile no finemap
sbatch test_template.sh 'test_mapping_elegans_nofinemap' '--vcf 20220216 --traitfile ../../input_data/c_elegans/phenotypes/test_pheno.tsv --out test_mapping_elegans_nofinemap --finemap false -resume -N kathrynevans2015@u.northwestern.edu'

# step 2: test c_elegans genomatrix pofile
sbatch test_template.sh 'test_genomatrix_elegans' '--vcf 20220216 --strains ../../input_data/c_elegans/phenotypes/strain_file.tsv -profile genomatrix --out test_genomatrix_elegans -N kathrynevans2015@u.northwestern.edu'

# step 2: test using own VCF
sbatch test_template.sh 'test_userVCF' '--vcf /projects/b1059/data/c_elegans/WI/variation/20220216/vcf/WI.20220216.small.hard-filter.isotype.vcf.gz --traitfile ../../input_data/c_elegans/phenotypes/test_pheno.tsv --out test_userVCF -N kathrynevans2015@u.northwestern.edu'

# step 4: test c_elegans simulations profile
sbatch test_template.sh 'test_simulations_elegans' '--vcf 20220216 -profile simulations --out test_simulations_elegans --simulate_nqtl ../../input_data/all_species/simulate_nqtl.csv \
--simulate_h2 ../../input_data/all_species/simulate_h2.csv --simulate_eff ../../input_data/all_species/simulate_effect_sizes.csv --simulate_strains ../../input_data/all_species/simulate_strains.tsv \
--simulate_maf ../../input_data/all_species/simulate_maf.csv -N kathrynevans2015@u.northwestern.edu'

# test sthresh
sbatch test_template.sh 'test_sthresh' '--vcf 20220216 --traitfile ../../input_data/c_elegans/phenotypes/test_pheno.tsv --out test_sthresh --mediation false --sthresh 4.5 -resume -N kathrynevans2015@u.northwestern.edu'

# test MAF 
sbatch test_template.sh 'test_maf' '--vcf 20220216 --traitfile ../../input_data/c_elegans/phenotypes/test_pheno.tsv --out test_maf --finemap false --mediation false --maf 0.1 -resume -N kathrynevans2015@u.northwestern.edu'

# test giving the wrong vcf 
sbatch test_template.sh 'error_vcf' '--vcf 20220215 --traitfile ../../input_data/c_elegans/phenotypes/test_pheno.tsv --out error_vcf --finemap false --mediation false -resume -N kathrynevans2015@u.northwestern.edu'

# test 20210121 vcf 
sbatch test_template.sh 'old_vcf' '--vcf 20210121 --traitfile ../../input_data/c_elegans/phenotypes/test_pheno.tsv --out old_vcf --finemap false --mediation false -resume -N kathrynevans2015@u.northwestern.edu'

# test download vcf from cendr
sbatch test_template.sh 'cendr_vcf' '--traitfile ../../input_data/c_elegans/phenotypes/test_pheno.tsv --out cendr_vcf --finemap false --mediation false -resume -N kathrynevans2015@u.northwestern.edu'

# make sure the mito qtl works
sbatch test_template.sh 'mito_test' '--traitfile ../../input_data/c_elegans/phenotypes/mito_test.tsv --out mito_test -resume -N kathrynevans2015@u.northwestern.edu'



# skip the local, gcp, and annotation profiles - test gcp next on gcp itself

# step 5: test c_tropicalis mapping
# send text when finishes because this is the last step
sbatch test_template.sh 'test_mapping_tropicalis' '--vcf 20210901 --traitfile ../../input_data/c_tropicalis/phenotypes/test_pheno.tsv --out test_mapping_tropicalis --sthresh EIGEN --species c_tropicalis -N kathrynevans2015@u.northwestern.edu'

# step 6: test c_briggsae mapping - not ready yet
sbatch test_template.sh 'test_mapping_briggsae' '--vcf 20210803 --traitfile ../../input_data/c_briggsae/phenotypes/test_pheno.tsv --out test_mapping_briggsae --species c_briggsae --sthresh EIGEN'

# when all these have finished and results look good... you can remove the entire directory