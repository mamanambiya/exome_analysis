#!/bin/bash

TODAY=$(date +"%Y_%m_%d__%H_%M_%S")
HOMEDIR="${HOME}/exome_analysis"
OUTDIR="/scratch/mamana/exome_aibst"
NAME="exome_analysis_nextflow"

. $HOME/.bashrc

## Activate environment
#source activate ngs


## Nextflowscript here
#cd ${HOMEDIR}
cd ${OUTDIR}
nextflow -log ${OUTDIR}/LOG/exome_analysis.nextflow_test.log \
    run ${HOMEDIR}/exome_analysis_test.nf \
    -c ${HOMEDIR}/HPC/${NAME}_test.config \
    -w ${OUTDIR}/work \
    -with-dag ${OUTDIR}/LOG/${NAME}.png \
    -resume \
    -profile slurm,Singularity
