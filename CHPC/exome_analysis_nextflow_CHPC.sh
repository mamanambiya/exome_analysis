#!/bin/bash
#PBS -q UCTlong
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -o /mnt/lustre/users/mmbiyavanga/exome_aibst/LOG/exome_analysis.out
#PBS -M mbymam001@myuct.ac.za
#PBS -m ae


TODAY=$(date +"%Y_%m_%d__%H_%M_%S")
HOMEDIR="/home/mmbiyavanga/exome_aibst"
OUTDIR="/mnt/lustre/users/mmbiyavanga/exome_aibst"
NAME="exome_analysis_nextflow_CHPC"

. $HOME/.bashrc

## Activate environment
source activate ngs


## Nextflowscript here
cd ${HOMEDIR}
#cd ${OUTDIR}
nextflow -log ${OUTDIR}/LOG/exome_analysis.nextflow.log \
    run ${HOMEDIR}/exome_analysis.nf \
    -c ${HOMEDIR}/CHPC/${NAME}.config \
    -w ${OUTDIR}/work \
    -resume \
    -with-dag ${OUTDIR}/LOG/${NAME}.png \
    -with-trace ${OUTDIR}/LOG/${NAME}.txt \
    -with-timeline ${OUTDIR}/LOG/${NAME}.html
