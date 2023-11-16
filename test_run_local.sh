#!/bin/bash
#SBATCH --account=b1042 ## Required: your allocation/account name, i.e. eXXXX, pXXXX or bXXXX
#SBATCH --partition=genomics ## Required: (buyin, short, normal, long, gengpu, genhimem, etc)
#SBATCH --time=04:00:00 ## Required: How long will the job need to run (remember different partitions h$
#SBATCH --nodes=1 ## how many computers/nodes do you need (no default)
#SBATCH --ntasks-per-node=4 ## how many cpus or processors do you need on per computer/node (default va$
#SBATCH --mem=30G ## how much RAM do you need per computer/node (this affects your FairShare score so b$
#SBATCH --job-name=test_local ## When you run squeue -u  this is how you can identify the job
#SBATCH --output=%j-%x.log ## standard out and standard error goes to this file



module load singularity
module load mamba
eval "$(conda shell.bash hook)"
conda activate /projects/b1059/software/conda_envs/nf20_env

nextflow run main.nf -profile simulations \
--species c_elegans \
--vcf /projects/b1059/projects/Ryan/NemaScan_Updates/NemaScan/input_data/c_elegans/genotypes/c_elegans.test.rename.vcf.gz \
--simulate_strains /projects/b1059/projects/Ryan/NemaScan_Updates/NemaScan/ce_test.vcf_input.strains.txt \
--simulate_nqtl /projects/b1059/projects/Ryan/Caenorhabditis_GWAS/best_panel_subsample/20231017_CE_96.192_Alloutlier/nqtl.csv \
--simulate_h2 /projects/b1059/projects/Ryan/Caenorhabditis_GWAS/best_panel_subsample/20231017_CE_96.192_Alloutlier/h2.csv \
--simulate_eff /projects/b1059/projects/Ryan/Caenorhabditis_GWAS/best_panel_subsample/20231017_CE_96.192_Alloutlier/effect_sizes.csv \
--simulate_reps 5 \
--out local_execution_test
