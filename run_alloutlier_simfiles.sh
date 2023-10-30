#!/bin/bash
#SBATCH --account=b1042 ## Required: your allocation/account name, i.e. eXXXX, pXXXX or bXXXX
#SBATCH --partition=genomicsguestA ## Required: (buyin, short, normal, long, gengpu, genhimem, etc)
#SBATCH --time=3:00:00 ## Required: How long will the job need to run (remember different partitions h$
#SBATCH --nodes=1 ## how many computers/nodes do you need (no default)
#SBATCH --ntasks-per-node=1 ## how many cpus or processors do you need on per computer/node (default va$
#SBATCH --mem=20G ## how much RAM do you need per computer/node (this affects your FairShare score so b$
#SBATCH --job-name=run_simfiles ## When you run squeue -u  this is how you can identify the job
#SBATCH --output=20230824_sim_files_0.01.log ## standard out and standard error goes to this file

module load singularity

nextflow run prepare_sims.nf --ld true --out 20230824_sim_files_ld_0.01 