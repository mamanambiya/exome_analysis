#!/bin/sh
#SBATCH --account=cbio
#SBATCH --partition=ada
#SBATCH --nodes=1 --ntasks=5 --mem=7000
#SBATCH --time=72:00:00
#SBATCH --job-name="exome_analysis_nextflow"
#SBATCH --mail-user=mbymam001@myuct.ac.za
#SBATCH --mail-type=BEGIN,END,FAIL
#SBTACH -o /scratch/mamana/exome_aibst/LOG/exome_analysis.nextflow.out

/home/mamana/exome_analysis/HPC/exome_analysis_nextflow.sh
