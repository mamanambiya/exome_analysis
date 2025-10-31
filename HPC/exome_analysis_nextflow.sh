#!/usr/bin/env bash

#SBATCH --nodes=1 --ntasks=5 --mem=7000
#SBATCH --time=48:00:00
#SBATCH --job-name="exome1"
#SBATCH --mail-user=mbymam001@myuct.ac.za
#SBATCH --mail-type=BEGIN,END,FAIL
#SBTACH -o /scratch/mamana/exome_aibst/LOG/.nextflow.out

TODAY=$(date +"%Y_%m_%d__%H_%M_%S")
HOMEDIR="${HOME}/exome_aibst"
OUTDIR="/scratch/mamana/exome_aibst"
NAME="exome_analysis_nextflow"

cd ${OUTDIR}
nextflow -log exome_analysis.nextflow.log \
    run ~/exome_aibst/main.nf \
    -c ~/exome_aibst/HPC/exome_analysis_nextflow.config \
    -profile singularity,slurm -resume

# nextflow ~/exome_aibst/main.nf -c ~/exome_aibst/ILIFU/exome_analysis_nextflow.config -profile singularity,slurm -resume

# nextflow -log ${OUTDIR}/LOG/exome_analysis.nextflow.log \
#     run ${HOMEDIR}/main.nf \
#     -c ${HOMEDIR}/ILIFU/${NAME}.config \
#     -w ${OUTDIR}/work \
#     -resume \
#     -with-trace ${OUTDIR}/LOG/${NAME}.txt \
#     -with-timeline ${OUTDIR}/LOG/${NAME}.html \
#     -profile slurm,Singularity
# #    -with-dag ${OUTDIR}/LOG/${NAME}.dot \

#sbatch /users/mamana/exome_aibst/ILIFU/exome_analysis_nextflow.sh
