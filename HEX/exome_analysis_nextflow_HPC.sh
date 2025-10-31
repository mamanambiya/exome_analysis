#!/bin/bash
#PBS -q UCTlong
#PBS -l nodes=1:ppn=1:series600
#PBS -j oe
#PBS -o /researchdata/fhgfs/mamana/exome_aibst/LOG/exome_analysis.out
#PBS -M mbymam001@myuct.ac.za
#PBS -m ae


TODAY=$(date +"%Y_%m_%d__%H_%M_%S")
HOMEDIR="${HOME}/exome_aibst"
OUTDIR="/researchdata/fhgfs/mamana/exome_aibst"
NAME="exome_analysis_nextflow_HPC"

. $HOME/.bashrc

## Activate environment
source activate ngs


## Nextflowscript here
#cd ${HOMEDIR}
cd ${OUTDIR}
nextflow -log ${OUTDIR}/LOG/exome_analysis.nextflow.log \
    run ${HOMEDIR}/exome_analysis.nf \
    -c ${HOMEDIR}/HPC/${NAME}.config \
    -w ${OUTDIR}/work \
    -resume \
    -with-trace ${OUTDIR}/LOG/${NAME}.txt \
    -with-timeline ${OUTDIR}/LOG/${NAME}.html \
    -profile pbs
#    -with-dag ${OUTDIR}/LOG/${NAME}.dot \

"UniqueID", "POS", "REF", "ALT", "AC", "Gene", "KG_AF", "KG_AFR_AF", "KG_EUR_AF", "KG_AMR_AF", "KG_EAS_AF", "gnomAD_AF", "gnom\
AD_AFR_AF", "gnomAD_FIN_AF", "ExAC_AF", "ExAC_AFR_AF", "AGVP_AF", "SAHGP_AF", "Effect", "CLNDN" "CDS" "GWASCAT_TRAIT" "GWASCAT_P_VALUE" "GWASCAT_PUBMED_ID