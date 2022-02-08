#!/bin/bash
#SBATCH -J nemascan_test       
#SBATCH -A b1042               
#SBATCH -p genomicsguestA               
#SBATCH -t 12:00:00            
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=16G               
#SBATCH -o nemascan.test.slurm.out


# script to test all profiles of nemascan to make sure everything is working for a new release

# first: make sure to commit all changes
# second: make new directory (i.e. 20211020_test_debug) and cd into this directory then run `sbatch ../test_debug.sh`

module purge
module load python/anaconda3.6
source activate nf20_env
module add singularity

# step 0: print out git commit
git status > git.log
git log | head -n 6 >> git.log

# step 1: test c_elegans mapping profile
nextflow run ../main.nf --vcf 20210121 --traitfile ../input_data/c_elegans/phenotypes/test_pheno.tsv --out test_mapping_elegans -resume -N 2088410137@vtext.com

# step 1: test c_elegans mapping profile no mediation
nextflow run ../main.nf --vcf 20210121 --traitfile ../input_data/c_elegans/phenotypes/test_pheno.tsv --out test_mapping_elegans_nomediate --mediation false -resume -N 2088410137@vtext.com

# step 1: test c_elegans mapping profile no finemap
nextflow run ../main.nf --vcf 20210121 --traitfile ../input_data/c_elegans/phenotypes/test_pheno.tsv --out test_mapping_elegans_nofinemap --finemap false -resume -N 2088410137@vtext.com

# step 2: test c_elegans genomatrix pofile
nextflow run ../main.nf --vcf 20210121 --strains ../input_data/c_elegans/phenotypes/strain_file.tsv -profile genomatrix --out test_genomatrix_elegans -N 2088410137@vtext.com

# step 2: test using own VCF
nextflow run ../main.nf --vcf '/projects/b1059/data/c_elegans/WI/variation/20210121/vcf/WI.20210121.small.hard-filter.isotype.vcf.gz' --traitfile ../input_data/c_elegans/phenotypes/test_pheno.tsv --out test_userVCF -N 2088410137@vtext.com

# step 3: test c_elegans mappings_docker profile
nextflow run ../main.nf --vcf 20210121 --traitfile ../input_data/c_elegans/phenotypes/test_pheno.tsv -profile mappings_docker --out test_mdocker_elegans -resume -N 2088410137@vtext.com

# step 4: test c_elegans simulations profile
nextflow run ../main.nf --vcf 20210121 -profile simulations --out test_simulations_elegans --simulate_nqtl ../input_data/all_species/simulate_nqtl.csv \
--simulate_h2 ../input_data/all_species/simulate_h2.csv --simulate_eff ../input_data/all_species/simulate_effect_sizes.csv --simulate_strains ../input_data/all_species/simulate_strains.tsv \
--simulate_maf ../input_data/all_species/simulate_maf.csv -N 2088410137@vtext.com

# test sthresh
nextflow run ../main.nf --vcf 20210121 --traitfile ../input_data/c_elegans/phenotypes/test_pheno.tsv --out test_sthresh --mediation false --sthresh 4.5 -resume -N 2088410137@vtext.com

# test MAF 
nextflow run ../main.nf --vcf 20210121 --traitfile ../input_data/c_elegans/phenotypes/test_pheno.tsv --out test_maf --finemap false --mediation false --maf 0.1 -resume -N 2088410137@vtext.com



# skip the local, gcp, and annotation profiles - test gcp next on gcp itself

# step 5: test c_tropicalis mapping
# send text when finishes because this is the last step
nextflow run ../main.nf --vcf 20210901 --traitfile ../input_data/c_tropicalis/phenotypes/test_pheno.tsv --out test_mapping_tropicalis --sthresh EIGEN --species c_tropicalis -N 2088410137@vtext.com

# step 6: test c_briggsae mapping - not ready yet
nextflow run ../main.nf --vcf 20210803 --traitfile ../input_data/c_briggsae/phenotypes/test_pheno.tsv --out test_mapping_briggsae --species c_briggsae --sthresh EIGEN

# when all these have finished and results look good... you can remove the entire directory