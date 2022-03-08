#!/bin/bash
#SBATCH -J nemascan_test       
#SBATCH -A b1042               
#SBATCH -p genomicsguestA               
#SBATCH -t 2:00:00            
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=16G

nextflow run ../main.nf ${1}