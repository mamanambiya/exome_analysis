#!/bin/bash
#PBS -q UCTlong
#PBS -l nodes=1:ppn=1:series600
#PBS -j oe
#PBS -o /scratch/mamana/exome_aibst/LOG/exome_analysis.out
#PBS -M mbymam001@myuct.ac.za
#PBS -m ae

TODAY=$(date +"%Y_%m_%d__%H_%M_%S")
HOMEDIR="${HOME}/exome_aibst"
OUTDIR="/scratch/mamana/exome_aibst"
NAME="exome_analysis_nextflow"

. $HOME/.bashrc

## Nextflow script here
##cd ${HOMEDIR}
cd ${OUTDIR}
nextflow -log ${OUTDIR}/LOG/exome_analysis.nextflow.log \
    run ${HOMEDIR}/exome_analysis.nf \
    -c ${HOMEDIR}/HPC/${NAME}.config \
    -w ${OUTDIR}/work \
    -resume \
    -with-trace ${OUTDIR}/LOG/${NAME}.txt \
    -with-timeline ${OUTDIR}/LOG/${NAME}.html \
    -profile slurm,Singularity
#    -with-dag ${OUTDIR}/LOG/${NAME}.dot \
