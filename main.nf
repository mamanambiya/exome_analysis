"""
Author: Mamana M.
Affiliation: University of Cape Town
Aim: A simple Snakemake workflow to download and process 1000 Genomes project data.
Date: Mon June 8 14:03:11 CET 2015
Run: snakemake -s down_1000g.snakefile --configfile config.yaml --timestamp -p -j 1000 --cluster "qsub {params.cluster}" --js jobscript.sh

Latest modification:
  - TODO
"""



CHRMS = params.chromosomes.split(',')
println "Project : $workflow.projectDir"
// println "Git info: $workflow.repository - $workflow.revision [$workflow.commitId]"
println "Cmd line: $workflow.commandLine"
println "Chromosomes used: ${CHRMS.join(',')}"



// All POP ..
def POPS_ALL = []
pop_to_dataset = [:]
pop_dataset_sample = []
params.POPS.each { entry ->
    POPS_ALL.addAll(entry.value.split(','))
    entry.value.split(',').each{ pop ->
        dataset = entry.key
        pop_to_dataset[pop] = dataset
        if (dataset in params.POP_samples.keySet()){
            pop_dataset_sample << [ pop, dataset, file(params.POP_samples[dataset])]
        }
    }
}


// All POP groups
def dataset_all = [:]
params.POPS.each { entry->
    dataset_all[entry.key] = entry.value.split(',')
}

def dataset_all1 = []
params.GROUP_POPS.each { entry->
    dataset_all1.add([entry.key, entry.value.split(',')])
}
dataset_all_C = Channel.from(dataset_all.keySet())
//dataset_all_KEY = params.GROUP_POPS.keySet().toList()
//dataset_all_VAL = params.GROUP_POPS.values().toList()

//GROUP_POPS_ANALYSIS = params.GROUP_POPS_ANALYSIS.split(',')

dataset_pops = [:]
params.POPS.each { entry ->
    dataset_pops[entry.key] = entry.value.split(',')
}
dataset_files_annot = params.dataset_files_annot.split(',')
dataset_files = []
dataset_files_all = []
dataset_files_ = [:]
dataset_files_only = [] // Full dataset, AIBSST only for now
dataset_files_chr = []
dataset_files_pop = []
params.dataset_full_files.each { dataset ->
    dataset_files_only << [dataset.key, file(params.dataset_full_files[dataset.key])]
}
params.dataset_files.each{ dataset ->
    if (!(dataset.key in dataset_files_.keySet())){
        dataset_files_[dataset.key] = []
    }
    CHRMS.each { chrm ->
        dataset_chrm = file(sprintf(dataset.value, chrm))
        if(!file(dataset_chrm).exists()){
            System.err.println "File ${dataset_chrm} not found. Please check your config file."
            exit 1
        }
        else {
            if (dataset.key in dataset_files_annot){
                dataset_files << [dataset.key, chrm, dataset_chrm]
                //dataset_files_all << [dataset.key, chrm, file(params.dataset_full_files[dataset.key])]
            }
            dataset_files_[dataset.key] << dataset_chrm
            dataset_files_chr << dataset.key+"=="+dataset_chrm
            params.dataset_pops.each{ dataset_pops ->
                if(dataset_pops.key == dataset.key){
                    dataset_pops.value.split(',').each { pop ->
                        params.dataset_samples.each { dataset_samples ->
                            if(dataset_samples.key == dataset.key) {
                                if(!file(dataset_samples.value).exists()){
                                    System.err.println "File ${dataset_samples.value} not found. Please check your config file."
                                    exit 1
                                }
                                else {
                                    dataset_files_pop << [dataset.key, pop, chrm, file(dataset_chrm), file(dataset_samples.value)]
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
//dataset_files = Channel.from(dataset_files)
dataset_files_chrs = Channel.from(dataset_files_all)
dataset_files_only_cha = Channel.from(dataset_files_only)

"""
Step: Split dataset per population
"""
process split_dataset_pop {
    tag "split_dataset_${dataset}_${pop}_${chrm}"
    label "medmem"
//    publishDir "${params.work_dir}/data/MERGED_POP/FRQ/${dataset}", mode:'symlink'
    input:
        set dataset, pop, chrm, file(dataset_chrm_vcf), file(dataset_sample) from dataset_files_pop
    output:
        set pop, chrm, file(pop_chrm_vcf), file(pop_sample_update), dataset into split_dataset_pop,split_dataset_pop_combine
    script:
        pop_chrm_vcf = "${pop}_${dataset}_chr${chrm}.vcf.gz"
        pop_sample = "${pop}.sample"
        pop_sample_update = "${pop}_update.sample"
        """
        grep ${pop} ${dataset_sample} | awk '{ \$2=\$2"\\t${dataset}"; print \$0}' > ${pop_sample_update}
        grep ${pop} ${dataset_sample} | awk '{print \$1}' > ${pop_sample}
        ## Keep only samples for population and Recalculate AC, AN, AF
        tabix ${dataset_chrm_vcf}
        bcftools view \
            --samples-file ${pop_sample} \
            ${dataset_chrm_vcf} | \
        bcftools annotate \
            -x INFO,^FORMAT \
            --set-id +'%CHROM\\_%POS\\_%REF\\_%ALT' | \
        bcftools +fill-tags | \
            bgzip -c > ${pop}.tmp3.vcf.gz
        bcftools +fixref \
            ${pop}.tmp3.vcf.gz \
            -Oz -o ${pop_chrm_vcf} -- \
            -f ${params.ref_genome} -m flip -d
        """
}

dataset_pop_combine = split_dataset_pop_combine.groupTuple()
process concat_dataset_pop {
    tag "concat_dataset_${pop}_${dataset}"
    label "bigmem"
    input:
        set pop, chrms, vcfs, samples, datasets from dataset_pop_combine
    output:
        set pop, file(vcf_out), sample, dataset into concat_dataset_pop,concat_dataset_pop_groups
    script:
        dataset = datasets[0]
        sample = samples[0]
        vcf_out = "${pop}_${dataset}.vcf.gz"
        """
        bcftools concat \
            ${vcfs.join(' ')} \
            -Oz -o ${pop}.tmp1.vcf.gz
        ## Recalculate AC, AN, AF
        bcftools +fill-tags ${pop}.tmp1.vcf.gz -Oz -o ${pop}.tmp2.vcf.gz
        bcftools sort ${pop}.tmp2.vcf.gz -Oz -o ${vcf_out}
        bcftools index --tbi -f ${vcf_out}
        rm ${pop}.tmp*.vcf.gz
        """
}


dataset_pop_groups = [:]
concat_dataset_pop_groups_list = concat_dataset_pop_groups.toSortedList().val
params.dataset_groups.each{ dataset_group ->
    group = dataset_group.key
    pops = dataset_group.value.split(',')
    pops.each { pop ->
        concat_dataset_pop_groups_list.each { pop_, pop_vcf, pop_sample, dataset ->
            if (pop == pop_) {
                if (!(group in dataset_pop_groups)){
                    dataset_pop_groups[group] = [group, file(pop_vcf), file(pop_sample)]
                }
                else {
                    dataset_pop_groups[group][1] = dataset_pop_groups[group][1] + " " + pop_vcf
                    dataset_pop_groups[group][2] = dataset_pop_groups[group][2] + " " + pop_sample
                }
            }
        }

    }
}


"""
Step: Merge populations into group
"""
process merge_pop_groups_1 {
    tag "merge_pop_groups_${group}"
    label "extrabig"
    input:
        set group, pop_vcfs, pop_samples from dataset_pop_groups.values()
    output:
        set group, file(vcf_out), file(sample_out) into merge_pop_groups_1
    script:
        vcf_out = "${group}.tmp2.vcf.gz"
        sample_out = "${group}.sample"
        """
        bcftools merge \
            ${pop_vcfs} \
            -Oz -o ${vcf_out}
        cat ${pop_samples} > ${sample_out}
        """
}

"""
Step: Merge populations into group - Recalculate AC, AF, MAF
"""
process merge_pop_groups_2 {
    tag "merge_pop_groups_${group}"
    label "extrabig"
    input:
    set group, file(group_vcf), file(group_sample) from merge_pop_groups_1

    output:
    set group, file(vcf_out), file(group_sample) into merge_pop_groups_2

    script:
    vcf_out = "${group}.tmp2.vcf.gz"
    """
    ## Recalculate AC, AN, AF
    bcftools +fill-tags ${group_vcf} -Oz -o ${vcf_out}
    """
}

"""
Step: Merge populations into group - Sort VCF
"""
process merge_pop_groups_3 {
    tag "merge_pop_groups_${group}"
    label "extrabig"
    input:
    set group, file(group_vcf), file(group_sample) from merge_pop_groups_2

    output:
    set group, file(vcf_out), file(group_sample) into merge_pop_groups_3

    script:
    vcf_out = "${group}.tmp3.vcf.gz"
    """
    bcftools sort ${group_vcf} -Oz -o ${vcf_out}
    """
}

"""
Step: Merge populations into group - Set IDs to POS_REF_ALT
"""
process merge_pop_groups {
    tag "merge_pop_groups_${group}"
    label "extrabig"
    input:
    set group, file(group_vcf), file(group_sample) from merge_pop_groups_3

    output:
    set group, file(group_vcf), file(group_sample) into merge_pop_groups
    set group, file(vcf_out_plink), file(group_sample) into merge_pop_groups_plink

    script:
    vcf_out_plink = "${group}_plink.vcf.gz"
    """
    bcftools annotate \
        -x INFO,^FORMAT \
        --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' ${group_vcf} -Oz -o ${vcf_out_plink}
    """
}

//
// """
// Step: Convert from VCF to plink for PCA analysis
// """
// process group_vcf_to_plink {
//     tag "group_vcf_to_plink_${group}"
//     label "extrabig"
//     input:
//         set group, file(group_vcf), file(group_sample) from merge_pop_groups_plink
//     output:
//         set group, file(plink_bed), file(plink_bim), file(plink_fam), file(group_sample) into group_vcf_to_plink
//     script:
//         plink_bed = "${group}_pruned.bed"
//         plink_bim = "${group}_pruned.bim"
//         plink_fam = "${group}_pruned.fam"
//         """
//         nblines=\$(zcat ${group_vcf} | grep -v '^#' | wc -l)
//         if (( \$nblines > 0 ))
//         then
//             plink2 \
//                 --vcf ${group_vcf} \
//                 --indep-pairwise 50 5 0.5 \
//                 --allow-no-sex \
//                 --make-bed \
//                 --snps-only --max-alleles 2  \
//                 --out ${group}
//             plink2 \
//                 --vcf ${group_vcf} \
//                 --extract ${group}.prune.in \
//                 --make-bed --snps-only --max-alleles 2 \
//                 --out ${group}_pruned
//             rm -rf ${group}.{bed,bim,fam}
//         else
//             touch ${plink_bed}
//             touch ${plink_bim}
//             touch ${plink_fam}
//         fi
//         """
// }
//
//
// """
// Step: Run Smartpca PCA analysis for group
// """
// process smartpca_group {
//     tag "smartpca_group_${group}"
//     label "extrabig"
// //    publishDir "${params.work_dir}/data/${dataset}/ALL/VCF", mode: 'symlink'
//     input:
//         set group, file(group_bed), file(group_bim), file(group_fam), file(group_sample) from group_vcf_to_plink
//     output:
//         set group, file(group_evec), file(group_eval), file(group_grmjunk), file(group_sample) into smartpca_group
//     script:
//         group_evec = "${group}.evec"
//         group_eval = "${group}.eval"
//         group_grmjunk = "${group}.evec_grmjunk"
//         """
//         nblines=\$(cat ${group_bim} | grep -v '^#' | wc -l)
//         if (( \$nblines > 0 ))
//         then
//             plink2 \
//                 --bfile ${group_bed.baseName} \
//                 --allow-no-sex \
//                 --recode \
//                 --out ${group}
//             ## Create parameter file for smartpca
//             echo -e \
//             "genotypename:    ${group}.ped
//             snpname:         ${group}.map
//             indivname:       ${group_fam}
//             evecoutname:     ${group_evec}
//             evaloutname:     ${group_eval}
//             altnormstyle:    NO
//             numoutevec:      5
//             numoutlieriter:  2
//             familynames:     NO
//             grmoutname:      ${group_grmjunk}" > ${group}.EIGENSTRAT.par
//             ## Run smartpca
//             smartpca \
//                     -p ${group}.EIGENSTRAT.par \
//                     > ${group}.EIGENSTRAT.log
//             rm -rf ${group}.{ped,map}
//         else
//             touch ${group_eval}
//         fi
//         """
// }
//
//
// """
// Step: Add group to evec file from smartpca
// """
// process update_evec {
//     tag "update_evec_${group}"
//     label "medmem"
// //    publishDir "${params.work_dir}/data/${dataset}/ALL/VCF", mode: 'symlink'
//     input:
//         set group, file(group_evec), file(group_eval), file(group_grmjunk), file(group_sample) from smartpca_group
//     output:
//         set group, file(group_evec_update), file(group_eval), file(group_grmjunk), file(group_sample) into update_evec
//     script:
//         group_evec_update = "${file(group_evec).baseName}_update.evec"
//         evec_file = group_evec
//         evec_out = group_evec_update
//         annot_file = group_sample
//         template "update_evec.py"
// }
//
//
// """
// Step: Plot PCA analysis for group
// """
// process plot_pca_group {
//     tag "plot_pca_group_${group}"
//     label "medmem"
//     publishDir "${params.work_dir}/REPORTS/PCA", mode:'copy', pattern: "*tiff"
//     input:
//         set group, file(group_evec), file(group_eval), file(group_grmjunk), file(group_sample) from update_evec
//     output:
//         set group, file(group_evec), file(group_sample), file("${output_pdf}*tiff") into plot_pca_group
//     script:
//         output_pdf = "${group}"
//         input_evec = group_evec
//         template "plot_pca.R"
// }
//

'''
Step: Prepare all data to be used
'''
process dataPreparation {
    tag "dataPrep"
    label "small"
    publishDir "${params.work_dir}/PGX_DATA", mode:'copy'

    output:
        file("pgxClinicalrsID.txt") into clinVarData_all
    script:
        """
        # The clinicalAnnotations.csv was downloaded from PharmGKB on 03 July 2017
        # This was used to find rsIDs for genes in this study

        # Select rsIDs from file and use this list to search annotated vcf
        # awk -F',' 'NR>1 {gsub(/"/, "", \$1); print \$1}' ${params.clinicalVariants_db} > pgxClinicalrsID.txt
        awk 'NR>1 {print \$1"\\t"\$2"\\t"\$3}' ${params.clinicalVariants_db} > pgxClinicalrsID.txt
        """
}



//
//"""
//Step 1.2: Annotate original VCFs AIBST using snpEff
//"""
//process annotate_dataset_snpeff_orig{
//    tag "snpeff_${file(vcf_file.baseName).baseName}"
//    label "bigmem"
//    publishDir "${params.work_dir}/data/${dataset}/VCF/", mode:'symlink'
//    input:
//        set val(dataset), file(vcf_file) from dataset_files_only_cha
//    output:
//        set val(dataset), file("${vcf_out}.gz") into annotate_dataset_snpeff_org
//        set val(dataset), file("${base}_snpeff.html"), file("${base}_snpeff.csv") into annotate_dataset_snpeff_org_extra
//    script:
//        base = file(vcf_file.baseName).baseName
//        vcf_out = "${base}_snpeff.vcf"
//        """
//        snpEff \
//            ${params.snpEff_human_db} -lof \
//            -stats ${base}_snpeff.html \
//            -csvStats ${base}_snpeff.csv \
//            -dataDir ${params.snpEff_database} \
//            ${vcf_file} -v > ${vcf_out}
//        bgzip -f ${vcf_out}
//        """
//}

//
//"""
//Step 1: Split VCF AIBST
//"""
//process split_vcf{
//    tag "split_vcf_${dataset}_${chrm}"
//    label "bigmem"
//    publishDir "${params.work_dir}/data/${dataset}/VCF/CHRS", mode:'copy'
//    input:
//        set val(dataset), val(chrm), file(vcf_file) from dataset_files_chrs
//    output:
//        set val(dataset), val(chrm), file(vcf_out), file("${vcf_out}.tbi") into split_vcf,split_vcf_1,split_vcf_2
//    script:
//        vcf_file = params.dataset_full_files[dataset]
//        vcf_out = "${file(file(vcf_file).baseName).baseName}_chr${chrm}.vcf.gz"
//        """
//        vcftools \
//            --gzvcf ${vcf_file} \
//            --chr ${chrm} \
//            --recode --stdout \
//            | bgzip -c > ${vcf_out}
//        bcftools index --tbi ${vcf_out}
//        """
//}


"""
Step 2: Filter Biallelic sites only in VCFs AIBST
"""
process biall_dataset {
    tag "biall_dataset_${file(vcf_file.baseName).baseName}"
    label "small"
    publishDir "${params.work_dir}/data/${dataset}/VCF/CHRS", mode: 'symlink'
    input:
        set val(dataset), val(chrm), file(vcf_file) from dataset_files
    output:
        set val(dataset), val(chrm), file(vcf_out), file("${vcf_out}.tbi") into biall_dataset,biall_dataset_1,biall_dataset_2
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_biall.vcf.gz"
        """
        tabix ${vcf_file}
        bcftools \
            view -m2 -M2 -v snps \
            ${vcf_file} \
            -Oz -o ${vcf_out}
        bcftools index --tbi -f ${vcf_out}
        """
}

"""
Step 3: Phase VCFs AIBST using eagle
"""
process phase_dataset {
    tag "phase_dataset_${file(vcf_file.baseName).baseName}"
    label "small"
    publishDir "${params.work_dir}/data/${dataset}/VCF/CHRS", mode: 'copy'
    input:
        set val(dataset), val(chrm), file(vcf_file), file(vcf_file_tbi) from biall_dataset_1
    output:
        set val(dataset), val(chrm), file("${vcf_out}.vcf.gz"), file("${vcf_out}.vcf.gz.tbi") into phase_dataset,phase_dataset_1
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_phased"
        """
        eagle \
            --vcf=${vcf_file} \
            --geneticMapFile=${params.genetic_map} \
            --chrom=${chrm} \
            --vcfOutFormat=z \
            --outPrefix=${vcf_out} > /dev/null 2>&1
        bcftools index --tbi -f ${vcf_out}.vcf.gz
        """
}


"""
Step 4: Annotate VCFs AIBST using snpEff
"""
process annotate_dataset_snpeff{
    tag "snpeff_${file(vcf_file.baseName).baseName}"
    label "extrabig"
    publishDir "${params.work_dir}/data/${dataset}/ANN/CHRS", mode:'symlink'
    input:
        set val(dataset), val(chrm), file(vcf_file) from dataset_files
    output:
        set val(dataset), val(chrm), file("${vcf_out}.gz"), file("${vcf_out}.gz.tbi") into annotate_dataset_snpeff
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_snpeff.vcf"
        """
        snpEff -Xmx${task.memory.toGiga()}g \
            ${params.snpEff_human_db} -lof \
            -stats ${dataset}_snpeff.html \
            -csvStats ${dataset}_snpeff.csv \
            -dataDir ${params.snpEff_database} \
            ${vcf_file} -v > ${vcf_out}
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        """
}
annotate_dataset_snpeff.into{ annotate_dataset_snpeff; annotate_dataset_snpeff_sub}

"""
Step 5: Annotate dbSNP IDs using snpsift
"""
annotate_dataset_snpeff.into { annotate_dataset_snpeff; annotate_dataset_snpeff_1 }
process annotate_dataset_dbsnp{
    tag "dbSNP_${file(vcf_file.baseName).baseName}"
    label "extrabig"
    publishDir "${params.work_dir}/data/${dataset}/ANN/CHRS", mode:'symlink'
    input:
        set val(dataset), val(chrm), file(vcf_file), file(vcf_file_tbi) from annotate_dataset_snpeff_1
    output:
        set val(dataset), val(chrm), file("${vcf_out}.gz"), file("${vcf_out}.gz.tbi") into annotate_dataset_dbsnp,annotate_dataset_dbsnp_6
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_dbsnp.vcf"
        """
        snpsift \
            annotate \
            ${params.snpEff_dbsnp_vcf} \
            -dataDir ${params.snpEff_database} \
            ${vcf_file} > ${vcf_out} -v
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        """
}
annotate_dataset_dbsnp.into{ annotate_dataset_dbsnp; annotate_dataset_dbsnp_sub}


"""
Step 6: Annotate GWAS Catalogue using snpsift
"""
annotate_dataset_dbsnp.into {annotate_dataset_dbsnp; annotate_dataset_dbsnp_1}
process annotate_dataset_gwascat{
    tag "gwascat_${file(vcf_file.baseName).baseName}"
    label "extrabig"
    publishDir "${params.work_dir}/data/${dataset}/ANN", mode:'symlink'

    input:
        set val(dataset), val(chrm), file(vcf_file), file(vcf_file_tbi) from annotate_dataset_dbsnp_1
    output:
        set val(dataset), val(chrm), file("${vcf_out}.gz"), file("${vcf_out}.gz.tbi") into annotate_dataset_gwascat
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_gwascat.vcf"
        """
        snpsift \
            gwasCat \
            -db ${params.gwascat_b37} \
            ${vcf_file} \
            -v \
            -dataDir ${params.snpEff_database} \
            > ${vcf_out}
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        """
}
annotate_dataset_gwascat.into{ annotate_dataset_gwascat; annotate_dataset_gwascat_sub}


'''
Step 7: Annotate with snpEff using clinvar
'''
annotate_dataset_gwascat.into { annotate_dataset_gwascat; annotate_dataset_gwascat_1}
process annotate_dataset_clinvar {
    tag "clinvar_${file(vcf_file.baseName).baseName}"
    label "extrabig"
    time = { 6.hour * task.attempt }
    publishDir "${params.work_dir}/data/${dataset}/ANN/CHRS", mode:'symlink'
    input:
        set val(dataset), val(chrm), file(vcf_file), file(vcf_file_tbi) from annotate_dataset_gwascat_1
    output:
        set val(dataset), val(chrm), file("${vcf_out}.gz"), file("${vcf_out}.gz.tbi") into annotate_dataset_clinvar
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_clinvar.vcf"
        """
        snpsift \
            annotate \
            ${params.clinvar} \
            ${vcf_file} \
            -v \
            -dataDir ${params.snpEff_database} \
            > ${vcf_out}
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        """
}


'''
Step 8: Annotate with snpEff using cosmic
'''
annotate_dataset_clinvar.into { annotate_dataset_clinvar; annotate_dataset_clinvar_2}
process annotate_cosmic_baylor {
    tag "cosmic_${file(vcf_file.baseName).baseName}"
    label "extrabig"
    time { 6.hour * task.attempt }
    publishDir "${params.work_dir}/data/${dataset}/ANN/CHRS", mode:'symlink'
    input:
        set val(dataset), val(chrm), file(vcf_file), file(vcf_file_tbi) from annotate_dataset_clinvar_2
    output:
        set val(dataset), val(chrm), file("${vcf_out}.gz"), file("${vcf_out}.gz.tbi") into annotate_dataset_cosmic
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_cosmic.vcf"
        """
        snpsift \
            annotate \
            ${params.cosmic} \
            ${vcf_file} \
            -v \
            -dataDir ${params.snpEff_database} \
            > ${vcf_out}
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        """
}
annotate_dataset_cosmic.into {annotate_dataset_cosmic; annotate_dataset_cosmic_sub}


'''
Step 9.1: combine frequency annotation files AGVP, SAHGP, TRYPANOGEN, gnomAD, ExAC per chromosome
'''
mafs_annotations = params.mafs_annotations
mafs_annotations_data = []
mafs_annotations.each{ dataset ->
    CHRMS.each{ chrm ->
        mafs_annotations_data << [dataset.key, chrm, file(dataset.value)]
    }
}
mafs_annotations_data_cha = Channel.from(mafs_annotations_data)

process mafs_annot_chrs {
    tag "mafs_annot_${mafs_dataset}_${chrm}"
    label "bigmem"
    publishDir "${params.reference_dir}/pop_mafs/${mafs_dataset}", mode:'copy'
    input:
        set val(mafs_dataset), val(chrm), file(mafs_file) from mafs_annotations_data_cha
    output:
        set val(mafs_dataset), val(chrm), file(outTSV) into mafs_annot_dataset
    script:
        outVCF = ""
        inVCF = ""
        inTSV = mafs_file
        outTSV = "${mafs_dataset}_chr${chrm}_mafs.tsv"
        template "annotateVCFwithTSV.py"
}


'''
Step 9.2: Annotate whole merge database with mafs from AGVP, SAHGP, TRYPANOGEN, gnomAD, ExAC
'''
mafs_annot_dataset.into { mafs_annot_dataset; mafs_annot_dataset_1 }
annotate_dataset_cosmic.into { annotate_dataset_cosmic; annotate_dataset_cosmic_1}
mafs_annot_dataset_list = mafs_annot_dataset_1.toSortedList().val
annotate_dataset_cosmic_list = annotate_dataset_cosmic_1.toSortedList().val
annotate_dataset_cosmic_2 = [:]
annotate_dataset_cosmic_list.each{ dataset, chrm, vcf_file, vcf_file_tbi ->
    annotate_dataset_cosmic_2[chrm] = [dataset, chrm, file(vcf_file)]
    annot = []
    datasets = []
    mafs_annot_dataset_list.each{ mafs_dataset, chr, mafs_file ->
        if (chrm == chr){
            annot << file(mafs_file)
            if( mafs_dataset.size() >= 6 ) {
                mafs_dataset = mafs_dataset[0..5]
            }
            datasets << mafs_dataset
        }
    }
    annotate_dataset_cosmic_2[chrm] << annot.join(';')
    annotate_dataset_cosmic_2[chrm] << datasets.join('-')
}
annotate_dataset_cosmic_2_cha = Channel.from(annotate_dataset_cosmic_2.values())
process annotate_mafs_dataset {
    tag "mafs_${file(vcf_file.baseName).baseName}"
    label "bigmem"
    input:
        set val(dataset), val(chrm), file(vcf_file), val(mafs_files), val(mafs_dataset) from annotate_dataset_cosmic_2_cha
    output:
        set val(dataset), val(chrm), file(outVCF) into annotate_mafs_dataset
    script:
        outVCF = "${file(vcf_file.baseName).baseName}_mafs-${mafs_dataset}.vcf"
        inVCF = vcf_file
        inTSV = mafs_files
        outTSV = ''
        template "annotateVCFwithTSV.py"
}

process annotate_mafs_dataset_gz {
    tag "mafs_${file(vcf_in.baseName).baseName}"
    label "bigmem"
    publishDir "${params.work_dir}/data/${dataset}/ALL/ANN/CHRS", mode:'symlink'
    input:
        set val(dataset), val(chrm), file(vcf_in) from annotate_mafs_dataset
    output:
        set val(dataset), val(chrm), file("${vcf_in}.gz") into annotate_mafs_dataset_gz,annotate_mafs_dataset_gz_1
    script:
        println chrm
        """
        bgzip -f ${vcf_in}
        """
}

//annotate_mafs_dataset.into {annotate_mafs_dataset; annotate_mafs_dataset_sub}

////
////'''
////Step 9.2.1: Annotate VCF for Ancestral Allele (AA) using in-house python script
////'''
////annotate_mafs_dataset.into {annotate_mafs_dataset; annotate_mafs_dataset_2}
////process add_ANC_to_VCF_merge {
////    tag "AA_${chrm}_${file(vcf_file.baseName).baseName}"
////    memory { 10.GB * task.attempt }
////    time { 8.hour * task.attempt }
////    publishDir "${params.work_dir}/VCF_ANN/MERGED", mode:'symlink'
////    input:
////        set val(dataset),  val(chrm), file(vcf_file), file(vcf_file_tbi) from annotate_mafs_dataset_2
////    output:
////        set val(dataset), val(chrm), file("${vcf_out}.gz"), file("${vcf_out}.gz.tbi") into vcf_anc_dataset
////    script:
////        vcf_out = "${file(vcf_file.baseName).baseName}_anc.vcf"
////        """
////        gunzip -c ${vcf_file} > ${vcf_file.baseName}
////        ~/.conda/envs/ngs_py27/bin/python2.7 ${params.homedir}/templates/add-ANC-to-vcf_new.py \
////            --in ${vcf_file.baseName} \
////            --genomedata ${params.genomedata_path} \
////            --out ${vcf_out}
////        bgzip -f ${vcf_out}
////        bcftools index --tbi -f ${vcf_out}.gz
////        rm -rf ${vcf_file.baseName}
////        """
////}
//
//
'''
Step 9.3: Reduced VCFs to PGX variants only
'''
//annotate_mafs_dataset_9_3 = vcf_anc_dataset_2.combine(gene_regions_slop_cha)
gene_regions_slop_cha = Channel.fromPath(params.gene_regions_slop)
annotate_mafs_dataset_9_3 = annotate_mafs_dataset_gz.combine(gene_regions_slop_cha)
process subset_pgx {
    tag "subset_pgx_${file(vcf_file.baseName).baseName}"
    label "bigmem"
    publishDir "${params.work_dir}/data/${dataset}/PGX_ONLY/ANN/CHRS", mode:'symlink'
    input:
        set val(dataset), val(chrm), file(vcf_file), file(gene_regions_slop) from annotate_mafs_dataset_9_3
    output:
        set val(dataset), val(chrm), file("${vcf_out}.gz") into subset_pgx
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_pgx.vcf"
        """
        tabix ${vcf_file}
        bcftools view ${vcf_file} \
            --regions-file ${gene_regions_slop} | \
        bgzip -c > ${file(vcf_file.baseName).baseName}.tmp.vcf.gz
        bcftools sort ${file(vcf_file.baseName).baseName}.tmp.vcf.gz -Oz -o ${vcf_out}.gz
        rm ${file(vcf_file.baseName).baseName}.tmp.vcf.gz
        ## TODO remove INTRON variants
        """
}


'''
Step : Split sample files by population
'''
mafs_annotations_data_cha = Channel.from(mafs_annotations_data)
process split_POP_samples {
    tag "split_${pop_sample}"
    label "small"
    publishDir "${params.work_dir}/data/samples/", mode:'copy'
    input:
        set val(POP), val(dataset), file(dataset_sample) from pop_dataset_sample
    output:
        set val(POP), file(pop_sample), val(dataset) into split_POP_samples
    script:
        pop_sample = "${POP}.sample"
        """
        awk '/${POP}/ {print \$0"\t${dataset}"}' ${dataset_sample} > ${pop_sample}
        """
}
split_POP_samples.into {split_POP_samples; split_POP_samples_sub}



'''
Step 4.0: Merge population sample files
'''
split_POP_samples.into { split_POP_samples; split_POP_samples_1 }
pop_samples_per_dataset = [:]
split_POP_samples_1.toSortedList().val.each { pop, pop_sample, dataset ->
    if(!(dataset in pop_samples_per_dataset.keySet())){
        pop_samples_per_dataset[dataset] = []
    }
    pop_samples_per_dataset[dataset] << file(pop_sample)
}
process merge_pop_sample{
    tag "merge_sample_${dataset}"
    label "small"
    publishDir "${params.work_dir}/data/samples/", mode:'symlink'
    input:
        val(dataset) from pop_samples_per_dataset.keySet()
    output:
        set val(dataset), file(dataset_sample), file("${dataset}.populations.txt"), file("${dataset}.superPopulations.txt") into merge_pop_sample
    script:
        dataset_sample = "${dataset}.sample"
        """
        cat ${pop_samples_per_dataset[dataset].join(' ')} > ${dataset_sample}
        # Make files with unique (super)population
        awk '{print \$2}' ${dataset_sample} | sort | uniq > ${dataset}.populations.txt
        awk '{print \$3}' ${dataset_sample} | sort | uniq > ${dataset}.superPopulations.txt
        """
}
merge_pop_sample.into{ merge_pop_sample; merge_pop_sample_sub}

//
//"""
//Step 10a: Combine chromosome VCFs AIBST into one
//"""
////vcf_anc_dataset.into{ vcf_anc_dataset; vcf_anc_dataset_1}
////annotate_mafs_dataset.into {annotate_mafs_dataset; annotate_mafs_dataset_3}
//vcf_anc_dataset_list = [:]
//merge_pop_sample.into{ merge_pop_sample; merge_pop_sample_1 }
//merge_pop_sample_1_list = merge_pop_sample_1.toSortedList().val
//annotate_mafs_dataset_3_list = annotate_mafs_dataset_gz_1.toSortedList().val
////vcf_anc_dataset_1.toSortedList().val.each { dataset, chrm, vcf_file, vcf_file_tbi ->
//annotate_mafs_dataset_3_list.each { dataset, chrm, vcf_file ->
//    if (!(dataset in vcf_anc_dataset_list.keySet())) {
//        merge_pop_sample_1_list.each{ dataset1, dataset_sample, dataset_populations, dataset_superPopulations ->
//            if (dataset1 == dataset){
//                vcf_anc_dataset_list[dataset] = [dataset, file(dataset_sample)]
//            }
//        }
//        vcf_anc_dataset_list[dataset][2] = ''
//    }
//    vcf_anc_dataset_list[dataset][2] += ' '+vcf_file
//}
//vcf_anc_dataset_cha = Channel.from(vcf_anc_dataset_list.values())
//process concat_dataset {
//    tag "concat_dataset_${dataset}"
//    label "small"
//    publishDir "${params.work_dir}/data/${dataset}/ALL/VCF", mode: 'symlink'
//    input:
//        set val(dataset), file(dataset_sample), val(dataset_vcfs) from vcf_anc_dataset_cha
//    output:
//        set val(dataset), file(dataset_sample), file(vcf_out), file("${vcf_out}.tbi") into concat_dataset_all,concat_dataset_all_140
//    script:
//        vcf_out = "${dataset}_fully_annot.vcf.gz"
//        """
//        bcftools concat \
//            ${dataset_vcfs} \
//            -Oz -o ${dataset}.tmp.vcf.gz
//        ## Recalculate AC, AN, AF
//        bcftools +fill-tags  ${dataset}.tmp.vcf.gz -Oz -o ${dataset}.tmp1.vcf.gz
//        bcftools sort ${dataset}.tmp1.vcf.gz -Oz -o ${vcf_out}
//        bcftools index --tbi -f ${vcf_out}
//        rm ${dataset}.tmp*.vcf.gz
//        """
//}


////
////"""
////Step 10.3: Annotate VCF AIBST using snpEff
////"""
////concat_dataset_all.into{ concat_dataset_all; concat_dataset_10_3}
////process annotate_dataset_all_snpeff{
////    tag "snpeff_${file(vcf_file.baseName).baseName}"
////    label "bigmem"
////    publishDir "${params.work_dir}/data/${dataset}/ALL/VCF", mode:'symlink'
////    input:
////        set val(dataset), file(dataset_sample), file(vcf_file), file(vcf_tbi) from concat_dataset_10_3
////    output:
////        set val(dataset), file(dataset_sample), file("${vcf_out}.gz"), file("${vcf_out}.gz.tbi"), file("${vcf_out}_snpeff.html"), file("${vcf_out}_snpeff.csv") into annotate_dataset_all_snpeff
////    script:
////        vcf_out = "${file(vcf_file.baseName).baseName}_snpeff.vcf"
////        """
////        snpEff \
////            ${params.snpEff_human_db} -lof \
////            -stats ${vcf_out}_snpeff.html \
////            -csvStats ${vcf_out}_snpeff.csv \
////            -dataDir ${params.snpEff_database} \
////            ${vcf_file} -v > ${vcf_out}
////        bgzip -f ${vcf_out}
////        bcftools index --tbi -f ${vcf_out}.gz
////        """
////}
////
////"""
////Annotate for dbNSFP
////"""
////annotate_dataset_all_snpeff.into{ annotate_dataset_all_snpeff; annotate_dataset_all_snpeff_10_3}
////process annotate_dataset_all_dbnsfp{
////    tag "snpeff_${file(vcf_file.baseName).baseName}"
////    label "bigmem"
////    publishDir "${params.work_dir}/data/${dataset}/ALL/VCF", mode:'symlink'
////    input:
////        set val(dataset), file(dataset_sample), file(vcf_file), file(vcf_tbi), file(vcf_snpeff), file(vcf_snpeff) from annotate_dataset_all_snpeff_10_3
////    output:
////        set val(dataset), file(dataset_sample), file("${vcf_out}.gz"), file("${vcf_out}.gz.tbi") into annotate_dataset_all_dbnsfp
////    script:
////        vcf_out = "${file(vcf_file.baseName).baseName}_dbnsfp.vcf"
////        """
////        snpsift dbnsfp \
////            -f Ancestral_allele,SIFT_score,SIFT_converted_rankscore,SIFT_pred,Uniprot_acc_Polyphen2,Uniprot_id_Polyphen2,Uniprot_aapos_Polyphen2,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred,LRT_score,LRT_converted_rankscore,LRT_pred,LRT_Omega,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MutationTaster_model,MutationTaster_AAE,MutationAssessor_UniprotID,MutationAssessor_variant,MutationAssessor_score,MutationAssessor_score_rankscore,MutationAssessor_pred,FATHMM_score,FATHMM_converted_rankscore,FATHMM_pred,PROVEAN_score,PROVEAN_converted_rankscore,PROVEAN_pred,Transcript_id_VEST3,Transcript_var_VEST3,VEST3_score,VEST3_rankscore,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,MetaLR_score,MetaLR_rankscore,MetaLR_pred,Reliability_index,M-CAP_score,M-CAP_rankscore,M-CAP_pred,REVEL_score,REVEL_rankscore,MutPred_score,MutPred_rankscore,MutPred_protID,MutPred_AAchange,MutPred_Top5features,CADD_raw,CADD_raw_rankscore,CADD_phred,DANN_score,DANN_rankscore,fathmm-MKL_coding_score,fathmm-MKL_coding_rankscore,fathmm-MKL_coding_pred,fathmm-MKL_coding_group,Eigen_coding_or_noncoding,Eigen-raw,Eigen-phred,Eigen-PC-raw,Eigen-PC-phred,Eigen-PC-raw_rankscore,GenoCanyon_score,GenoCanyon_score_rankscore,integrated_fitCons_score,integrated_fitCons_score_rankscore,integrated_confidence_value,GM12878_fitCons_score,GM12878_fitCons_score_rankscore,GM12878_confidence_value,H1-hESC_fitCons_score,H1-hESC_fitCons_score_rankscore,H1-hESC_confidence_value,HUVEC_fitCons_score,HUVEC_fitCons_score_rankscore,HUVEC_confidence_value,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phyloP20way_mammalian,phyloP20way_mammalian_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,phastCons20way_mammalian,phastCons20way_mammalian_rankscore,SiPhy_29way_pi,SiPhy_29way_logOdds,SiPhy_29way_logOdds_rankscore  \
////            -v -db ${params.dbnsfp_db} \
////            ${vcf_file} -v > ${vcf_out}
////        bgzip -f ${vcf_out}
////        bcftools index --tbi -f ${vcf_out}.gz
////        """
////}
//
//"""
//Step 10b: Combine chromosome VCFs AIBST into one for PGX variant dataset
//"""
//subset_pgx.into{ subset_pgx; subset_pgx_1 }
//vcf_anc_dataset_pgx_list = [:]
//subset_pgx_1.toSortedList().val.each { dataset, chrm, vcf_file, vcf_file_tbi ->
//    if (!(dataset in vcf_anc_dataset_pgx_list.keySet())) {
//        merge_pop_sample_1_list.each{ dataset1, dataset_sample, dataset_populations, dataset_superPopulations ->
//            if (dataset1 == dataset){
//                vcf_anc_dataset_pgx_list[dataset] = [dataset, file(dataset_sample)]
//            }
//        }
//        vcf_anc_dataset_pgx_list[dataset][2] = ''
//    }
//    vcf_anc_dataset_pgx_list[dataset][2] += ' '+vcf_file
//}
//vcf_anc_dataset_pgx = Channel.from(vcf_anc_dataset_pgx_list.values())
//process concat_dataset_pgx {
//    tag "concat_dataset_pgx_${dataset}"
//    label "small"
//    publishDir "${params.work_dir}/data/${dataset}/PGX_ONLY/VCF", mode: 'symlink'
//    input:
//        set val(dataset), file(dataset_sample), val(dataset_vcfs) from vcf_anc_dataset_pgx
//    output:
//        set val(dataset), file(dataset_sample), file(vcf_out), file("${vcf_out}.tbi") into concat_dataset_pgx,concat_dataset_pgx_6,concat_dataset_pgx_10c,concat_dataset_pgx_gene_stats,concat_dataset_all_14_0
//    script:
//        vcf_out = "${dataset}_annot_pgx.vcf.gz"
//        """
//        bcftools concat \
//            ${dataset_vcfs} \
//            -Oz -o ${dataset}.tmp.vcf.gz
//        zcat ${dataset}.tmp.vcf.gz | \
//        snpsift extractFields - -e "." CHROM POS ANN[*].GENE | \
//        grep -E "UGT|HLA" | \
//        awk '{print \$1"\\t"\$2}' > ${dataset}.ignore.tsv
//        bcftools view \
//            --targets-file ^${dataset}.ignore.tsv \
//            ${dataset}.tmp.vcf.gz \
//            -Oz -o ${dataset}.tmp1.vcf.gz
//        ## Recalculate AC, AN, AF
//        bcftools +fill-tags \
//            ${dataset}.tmp1.vcf.gz \
//            -Oz -o  ${dataset}.tmp2.vcf.gz
//        bcftools sort \
//            ${dataset}.tmp2.vcf.gz \
//            -Oz -o ${vcf_out}
//        bcftools index --tbi -f ${vcf_out}
//        rm ${dataset}.tmp*.vcf.gz
//        """
//}
//
//
//
//"""
//Step 10.3: Annotate PGx VCF AIBST using snpEff
//"""
//concat_dataset_pgx.into{ concat_dataset_pgx; concat_dataset_pgx_1}
//process annotate_dataset_pgx_snpeff{
//    tag "snpeff_${file(vcf_file.baseName).baseName}"
//    label "bigmem"
//    publishDir "${params.work_dir}/data/${dataset}/PGX_ONLY/VCF", mode:'symlink'
//    input:
//        set val(dataset), file(dataset_sample), file(vcf_file), file(vcf_tbi) from concat_dataset_pgx_1
//    output:
//        set val(dataset), file(dataset_sample), file("${vcf_out}.gz"), file("${vcf_out}.gz.tbi"), file("${vcf_out}_snpeff.html") into annotate_dataset_pgx_snpeff
//    script:
//        vcf_out = "${file(vcf_file.baseName).baseName}_snpeff.vcf"
//        """
//        snpEff \
//            ${params.snpEff_human_db} -lof \
//            -stats ${vcf_out}_snpeff.html \
//            -csvStats ${vcf_out}_snpeff.csv \
//            -dataDir ${params.snpEff_database} \
//            ${vcf_file} -v > ${vcf_out}
//        bgzip -f ${vcf_out}
//        bcftools index --tbi -f ${vcf_out}.gz
//        """
//}
//
//
//"""
//Annotate for dbNSFP
//"""
//concat_dataset_pgx.into{ concat_dataset_pgx; concat_dataset_pgx_2}
//process annotate_dataset_pgx_dbnsfp{
//    tag "snpeff_${file(vcf_file.baseName).baseName}"
//    label "bigmem"
//    publishDir "${params.work_dir}/data/${dataset}/ALL/VCF", mode:'symlink'
//    input:
//        set val(dataset), file(dataset_sample), file(vcf_file), file(vcf_tbi) from concat_dataset_pgx_2
//    output:
//        set val(dataset), file(dataset_sample), file("${vcf_out}.gz"), file("${vcf_out}.gz.tbi") into annotate_dataset_pgx_dbnsfp
//    script:
//        vcf_out = "${file(vcf_file.baseName).baseName}_dbnsfp.vcf"
//        """
//        snpsift dbnsfp \
//            -f Ancestral_allele,SIFT_score,SIFT_converted_rankscore,SIFT_pred,Uniprot_acc_Polyphen2,Uniprot_id_Polyphen2,Uniprot_aapos_Polyphen2,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred,LRT_score,LRT_converted_rankscore,LRT_pred,LRT_Omega,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MutationTaster_model,MutationTaster_AAE,MutationAssessor_UniprotID,MutationAssessor_variant,MutationAssessor_score,MutationAssessor_score_rankscore,MutationAssessor_pred,FATHMM_score,FATHMM_converted_rankscore,FATHMM_pred,PROVEAN_score,PROVEAN_converted_rankscore,PROVEAN_pred,Transcript_id_VEST3,Transcript_var_VEST3,VEST3_score,VEST3_rankscore,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,MetaLR_score,MetaLR_rankscore,MetaLR_pred,Reliability_index,M-CAP_score,M-CAP_rankscore,M-CAP_pred,REVEL_score,REVEL_rankscore,MutPred_score,MutPred_rankscore,MutPred_protID,MutPred_AAchange,MutPred_Top5features,CADD_raw,CADD_raw_rankscore,CADD_phred,DANN_score,DANN_rankscore,fathmm-MKL_coding_score,fathmm-MKL_coding_rankscore,fathmm-MKL_coding_pred,fathmm-MKL_coding_group,Eigen_coding_or_noncoding,Eigen-raw,Eigen-phred,Eigen-PC-raw,Eigen-PC-phred,Eigen-PC-raw_rankscore,GenoCanyon_score,GenoCanyon_score_rankscore,integrated_fitCons_score,integrated_fitCons_score_rankscore,integrated_confidence_value,GM12878_fitCons_score,GM12878_fitCons_score_rankscore,GM12878_confidence_value,H1-hESC_fitCons_score,H1-hESC_fitCons_score_rankscore,H1-hESC_confidence_value,HUVEC_fitCons_score,HUVEC_fitCons_score_rankscore,HUVEC_confidence_value,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phyloP20way_mammalian,phyloP20way_mammalian_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,phastCons20way_mammalian,phastCons20way_mammalian_rankscore,SiPhy_29way_pi,SiPhy_29way_logOdds,SiPhy_29way_logOdds_rankscore  \
//            -v -db ${params.dbnsfp_db} \
//            ${vcf_file} -v > ${vcf_out}
//        bgzip -f ${vcf_out}
//        bcftools index --tbi -f ${vcf_out}.gz
//        """
//}
//
////// ANALYSIS START HERE
//
//////// Generate dataflow for use in processes
//concat_dataset_pgx.into{ concat_dataset_pgx; concat_dataset_pgx_2}
//group_vcf_with_pop_sample = [] // For POP, dataset, Full GROUP_VCF, GROUP_SAMPLE, dataset_populations, dataset_superPopulations
//group_pop_data = [:] // contains annotated vcfs of all
//concat_dataset_pgx_2.toSortedList().val.each { dataset, dataset_sample, dataset_vcf, dataset_vcf_tbi ->
//    if ( dataset in dataset_files_annot) {
//        data = file(dataset_sample).readLines().collect { it.split()[1] }.unique()
//        data.each { POP ->
//            group_vcf_with_pop_sample << [POP, dataset, dataset_vcf, dataset_sample]
//        }
//        // All dataset vcf and sample files
//        if ( !(dataset in group_pop_data.keySet()) ){
//            group_pop_data[dataset] = [dataset_vcf, dataset_sample]
//        }
//    }
//}
//
//// TODO Start HERE
//'''
//Step 11b: Analysis of LoF
//'''
//concat_dataset_pgx.into{ concat_dataset_pgx; concat_dataset_pgx_1}
//process lof_analysis_pgx{
//    tag "lof_analysis_${dataset}"
//    label "small"
//    publishDir "${params.work_dir}/data/PGX_ONLY/LOF/${dataset}", mode:'symlink'
//    input:
//        set val(dataset), file(dataset_sample), file(dataset_vcf), file(dataset_vcf_tbi) from concat_dataset_pgx_1
//    output:
//        set val(dataset), file(dataset_vcf), file("${dataset_out}.lof${lof_cutoff1}.vcf.gz"), file("${dataset_out}.lof${lof_cutoff1}.vcf.gz.tbi"), file("${dataset_out}.lof${lof_cutoff1}.INFO"), file("${dataset_out}.lof${lof_cutoff1}.LoFConsequencesIDacGENE"), file("${dataset_out}.lof${lof_cutoff1}.LofPerSample.012"), file("${dataset_out}.lof${lof_cutoff1}.LofPerSample.012.indv"), file(dataset_sample) into lof_analysis_pgx
//    script:
//        dataset_out = "${file(dataset_vcf.baseName).baseName}"
//        lof_cutoff = params.lof_cutoff
//        lof_cutoff1 = lof_cutoff.replaceAll('\\.', '')
//        """
//        ## Filter VCf for LOF > 1
//        zcat ${dataset_vcf} | snpsift filter "LOF[*].PERC >= ${lof_cutoff}" | bgzip -c > ${dataset_out}.lof.tmp.vcf.gz
//        ## Remove multi-allelic variants # --max-alleles 2
//        ## Remove singleton and doubleton # --mac 2
//        ## Remove variants with a Lof_flag
//        ## TODO Curate variants with >30% allele freq: remove rs11356919 (AF 69%) due to lack of evidence
//        vcftools --gzvcf ${dataset_out}.lof.tmp.vcf.gz \
//           --mac 2 --recode-INFO-all --recode --stdout | bgzip -c > ${dataset_out}.lof1.vcf.gz
//        bcftools view -m2 -M2 -v snps ${dataset_out}.lof1.vcf.gz | bcftools +fill-tags | \
//        grep -v 'SINGLE_EXON|NAGNAG_SITE|PHYLOCSF_WEAK|PHYLOCSF_UNLIKELY_ORF|PHYLOCSF_TOO_SHORT' | \
//            bgzip -c > ${dataset_out}.lof${lof_cutoff1}.vcf.gz
//        # CALCULATE WHICH GENES HARBOUR LoF MUTATIONS and HOW MANY
//        # Extract information from vcf info file
//        zcat ${dataset_out}.lof${lof_cutoff1}.vcf.gz | snpsift extractFields - -e "." \
//            CHROM POS REF ALT AC ANN[0].GENE MAF KG_AF KG_AFR_AF KG_EUR_AF KG_AMR_AF KG_EAS_AF gnomAD_AF gnomAD_AFR_AF gnomAD_FIN_AF ExAC_AF ExAC_AFR_AF AGVP_AF SAHGP_AF "ANN[0].EFFECT" CLNDN CDS GWASCAT_TRAIT GWASCAT_P_VALUE GWASCAT_PUBMED_ID\
//            > ${dataset_out}.lof${lof_cutoff1}.INFO
//        bcftools index --tbi -f ${dataset_out}.lof${lof_cutoff1}.vcf.gz
//        # Create a unique ID for each of the LoF variants
//        awk 'NR>1 {print \$1":"\$2":"\$4}' ${dataset_out}.lof${lof_cutoff1}.INFO > ${dataset_out}.lof${lof_cutoff1}.UniqueID
//        # Get allele counts for each of the genes
//        awk 'NR>1 {print \$0}' ${dataset_out}.lof${lof_cutoff1}.INFO | cut -c5- > ${dataset_out}.lof${lof_cutoff1}.AC
//        # Get the relevant gene names for each of the LoF variants
//        #awk -F' ' 'NR>1 {print \$NF}' ${dataset_out}.lof${lof_cutoff1}.INFO > ${dataset_out}.lof${lof_cutoff1}.Gene
//        # Combine the Unique ID, AC/AF and Gene files
//        paste ${dataset_out}.lof${lof_cutoff1}.UniqueID ${dataset_out}.lof${lof_cutoff1}.AC > ${dataset_out}.lof${lof_cutoff1}.LoFConsequencesIDacGENE
//        #rm ${dataset_out}.lof${lof_cutoff1}.UniqueID ${dataset_out}.lof${lof_cutoff1}.AC ${dataset_out}.lof${lof_cutoff1}.Gene
//        # CALCULATE NUMBER OF LoF / SAMPLE
//        # Calculate the number of non-ref alleles per sample
//        vcftools --gzvcf ${dataset_vcf} --012 --out ${dataset_out}.lof${lof_cutoff1}.LofPerSample
//        """
//}
//
//
//'''
//Step 11c: Plot analysis of LoF
//'''
//lof_analysis_pgx.into{ lof_analysis_pgx; lof_analysis_pgx_1 }
//process plot_lof_analysis_pgx{
//    tag "plot_lof_${dataset}"
//    label "small"
//    publishDir "${params.work_dir}/REPORTS/${dataset}/PGX_ONLY/LOF/", mode:'copy'
//    input:
//        set val(dataset), file(pgx_dataset_vcf), file(dataset_vcf), file(dataset_vcf_tbi), file(dataset_INFO), file(dataset_LoFConsequencesIDacGENE), file(dataset_LofPerSample_012), file(dataset_LofPerSample_012_indv), file(dataset_sample) from lof_analysis_pgx_1
//    output:
//        set val(dataset), file(dataset_pgxLoFPerGeneCount), file(dataset_pgxLoFPerCombinedAF), file(dataset_allLofPop), file(dataset_LoFConsequencesIDacGENE) into plot_lof_analysis_pgx
//    script:
//        dataset_allLofPop = "${dataset}_allLofPop.csv"
//        dataset_pgxLoFPerGeneCount = "${dataset}_pgxLoFPerGeneCount.tiff"
//        dataset_pgxLoFPerCombinedAF = "${dataset}_pgxLoFPerCombinedAF.tiff"
//        template "step11b_lof_analysis.R"
//}
//
//
//'''
//Step 12: Population Frequency analysis
//'''
////merge_pop_vcf_all.into {merge_pop_vcf_all; merge_pop_vcf__lof_analysis}
////annot_vep_merge_pop_all.into { annot_vep_merge_pop_all; annot_vep_merge_pop__popFrq }
//concat_dataset_pgx.into{ concat_dataset_pgx; concat_dataset_pgx_2}
////merge_pop_sample__cha.into { merge_pop_sample__cha; merge_pop_sample__popFrq_cha }
//process popFrq_1 {
//    tag "popFreq_1_${dataset}"
//    label "small"
//    publishDir "${params.work_dir}/data/MERGED_POP/FRQ/${dataset}", mode:'symlink'
//
//    input:
//        set val(dataset), file(dataset_sample), file(dataset_vcf), file(dataset_vcf_tbi) from concat_dataset_pgx_2
//    output:
//        set val(dataset), file("${dataset}.biall.recode.vcf.gz"), file(dataset_sample), file("${dataset}.populations.txt"), file("${dataset}.superPopulations.txt") into popFrq_1_all
//        set val(dataset), file(dataset_frq) into dataset_pgx_summary,dataset_pgx_summary_12_7
//    script:
//        dataset_frq = "${dataset}.frq"
//        """
//        # Make files with unique (super)population
//        awk '{print \$2}' ${dataset_sample} | sort | uniq > ${dataset}.populations.txt
//        awk '{print \$3}' ${dataset_sample} | sort | uniq > ${dataset}.superPopulations.txt
//        # Calculate bi-allelic SNPs
//        bcftools view \
//            -m2 -M2 -v snps ${dataset_vcf} | \
//        bcftools +fill-tags | \
//        bcftools annotate \
//            --set-id +'%CHROM\\_%POS\\_%REF\\_%ALT' | \
//        bgzip -c > ${dataset}.biall.recode.vcf.gz
//        echo -e "CHRM:POS\\tAIBST_AF" > ${dataset_frq}
//        bcftools query \
//            -f '%CHROM:%POS\\t%INFO/MAF\\n' \
//            ${dataset}.biall.recode.vcf.gz >> ${dataset_frq}
//        """
//}
//
//popFrq_1_all.into { popFrq_1_all; popFrq_1__popFrq_2 }
//pop_samples_per_group_1 = [] // For POP, GROUP_POP, Biall GROUP_VCF, GROUP_SAMPLE, GROUP_POP_populations, GROUP_POP_superPopulations
//popFrq_1_list = popFrq_1__popFrq_2.toSortedList()
//popFrq_1_list.val.each { dataset, dataset_biall_vcf, dataset_sample, dataset_populations, dataset_superPopulations ->
//    data = file(dataset_populations).readLines().toSorted()
//    data.each { POP ->
//        pop_samples_per_group_1 << [POP, dataset, dataset_biall_vcf, dataset_sample, dataset_populations, dataset_superPopulations ]
//    }
//}
//pop_samples_per_group_1_chan = Channel.from(pop_samples_per_group_1)
//
//'''
//Step 12.1: Compute frequencies per population
//           - Create a file with each sample in a given population and
//           - Calculate allele frequencies of bi-allelic SNPS for individual populations based on the minor allele
//'''
//process popFrq_2 {
//    tag "popFreq_2_${POP}"
//    label "small"
//    publishDir "${params.work_dir}/data/MERGED_POP/FRQ/${dataset}", mode:'symlink'
//
//    input:
//        set val(POP), val(dataset), file(dataset_biall_vcf), file(dataset_sample), file(dataset_populations), file(dataset_superPopulations) from pop_samples_per_group_1_chan
//    output:
//        set val(POP), val(dataset), file("${POP}.samples"), file("${POP}.vcf.gz"), file("${POP}.Ind.frq"), file(report) into popFrq_2_all
//    script:
//        report = "${POP}_count.csv"
//        """
//        grep ${POP} ${dataset_sample} | cut -f1 > ${POP}.samples
//        ## Keep only samples for population and Recalculate AC, AN, AF
//        bcftools view \
//            --samples-file ${POP}.samples \
//            ${dataset_biall_vcf} | \
//        bcftools +fill-tags | \
//        bcftools annotate \
//            --set-id +'%CHROM\\_%POS\\_%REF\\_%ALT' | \
//            bgzip -c > ${POP}_tmp.vcf.gz
//        ## Compute frequency
//        echo "rsID\tREF\tALT\t${POP}_Alt_FREQ" > ${POP}.Ind.frq
//        bcftools view \
//            --min-ac 1:minor \
//            ${POP}_tmp.vcf.gz | \
//            bgzip -c > ${POP}.vcf.gz
//        bcftools query \
//            -f '%ID\\t%REF\\t%ALT\\t%INFO/MAF\\n' \
//            ${POP}_tmp.vcf.gz >> ${POP}.Ind.frq
//        vcftools \
//            --gzvcf ${dataset_biall_vcf} \
//            --singletons \
//            --out ${dataset}
//        grep ${POP} ${dataset_sample} | cut -f1 > ${POP}.samples
//        while read SAMPLE;
//            do grep \$SAMPLE ${dataset}.singletons;
//        done < ${POP}.samples >> ${POP}.singleton.count
//        echo -e "${POP};\$( zcat ${POP}.vcf.gz | grep -v "^CHR" | wc -l );\$( cat ${POP}.singleton.count | wc -l )" > ${report}
//        """
//}
//
//
//
//'''
//Step 12.2: Compute global frequencies for dataset
//            - Calculate the global allele frequency of bi-allelic SNPs
//'''
//popFrq_1_all.into { popFrq_1_all; popFrq_1__popFrq_3 }
//process popFrq_3 {
//    tag "popFreq_3_${dataset}"
//    label "small"
//    publishDir "${params.work_dir}/data/MERGED_POP/FRQ/${dataset}", mode:'symlink'
//
//    input:
//        set val(dataset), file(dataset_biall_vcf), file(dataset_sample), file(dataset_populations), file(dataset_superPopulations) from popFrq_1__popFrq_3
//    output:
//        set val(dataset), file("${dataset}.biall.recode.vcf.gz"), file(dataset_sample), file("${dataset}.biall.INFO.rsID") into popFrq_3_all
//    script:
//        """
//        # Create file with allele count/number/freq INFO fields from vcf
//        zcat ${dataset_biall_vcf} \
//            | snpsift extractFields - -e "." CHROM POS ID REF ALT AC AN AF "ANN[0].GENE" MAF KG_AFR_AF ExAC_AF ExAC_AFR_AF gnomAD_AF gnomAD_AFR_AF AGVP_AF "ANN[0].EFFECT" \
//            > ${dataset}.biall.INFO
//        echo "UniqueID\tCHROM\tPOS\tID\tREF\tALT\tAC\tAN\tAF\tGene\tAIBST_Alt_FREQ\tKG_AFR_Alt_FREQ\tExAC_Alt_FRE\tExAC_AFR_Alt_FRE\tgnomAD_Alt_FRE\tgnomAD_AFR_Alt_FREQ\tAGVP_Alt_FREQ\tEFFECT" > ${dataset}.biall.INFO.rsID
//        awk 'NR>1 {print \$1":"\$2":"\$5"\t"\$0}' ${dataset}.biall.INFO >> ${dataset}.biall.INFO.rsID
//        #vcftools --gzvcf ${dataset_biall_vcf} \
//        #    --get-INFO AC --get-INFO AN --get-INFO AF --get-INFO MAF \
//        #    --out ${dataset}.biall
//        # Create file with rsIDs for vcf (double hash saves header)
//        #zcat ${dataset}.biall.vcf.gz | grep -v '##' | awk '{print \$3}' > ${dataset}.biall.rsID
//        # Combine INFO and rsID files
//        #paste ${dataset}.biall.INFO ${dataset}.biall.rsID > ${dataset}.biall.INFO.rsID
//        """
//}
//
//
//'''
//Step 12.3: Compute global frequencies for dataset
//'''
//popFrq_3_all.into{popFrq_3_all; popFrq_3__2}
//popFrq_2_all.into{popFrq_2_all; popFrq_2__2}
//Ind_group = [:]
//Ind_group1 = [:]
//temp = [:]
//temp1 = [:]
//temp2 = [:]
//popFrq_2__2.toSortedList().val.each { POP, dataset, POP_sample, POP_vcf, POP_Ind_frq, POP_count ->
//    if ( !(dataset in Ind_group.keySet()) ){
//        Ind_group[dataset] = [dataset]
//        temp[dataset] = []
//        temp1[dataset] = []
//    }
//    temp[dataset] << POP+":"+POP_Ind_frq
//    temp1[dataset] << POP
//    if ( !(dataset in Ind_group1.keySet()) ){
//        Ind_group1[dataset] = [dataset]
//        temp2[dataset] = []
//    }
//    temp2[dataset] << POP+":"+POP_count
//}
//popFrq_3__2.toSortedList().val.each{ dataset, dataset_biall_vcf, GROUP_sample, GROUP_biall_info_rsid ->
//    Ind_group[dataset] << GROUP_biall_info_rsid
//    Ind_group[dataset] << temp[dataset].toSorted().collect{ it.split(':')[1] }.join(' ')
//    Ind_group[dataset] << temp1[dataset].toSorted().join(' ')
//
//    Ind_group1[dataset] << temp2[dataset].toSorted().collect{ it.split(':')[1] }.join(' ')
//    Ind_group1[dataset] << temp1[dataset].toSorted().join(' ')
//}
//
//
//
//process popFrq_3_1 {
//    tag "popFreq_3_1_${dataset}"
//    label "small"
//    publishDir "${params.work_dir}/data/${dataset}/ALL/VCF", mode:'symlink'
//    input:
//        val(dataset) from Ind_group.keySet()
//    output:
//        set val(dataset), file("${dataset}_dataset.csv"), file("${dataset}_datasetAnnotated.csv") into popFrq_3_1_all
//    script:
//        GROUP_biall_info_rsid = Ind_group[dataset][1]
//        POP_Ind_frq_files = file(Ind_group[dataset][2])
//        POPs = Ind_group[dataset][3]
//        template "create_datasetAnnotated.R"
//}
//
//process popFrq_combine_count {
//    tag "popFrq_combine_count_${dataset}"
//    label "small"
//    publishDir "${params.work_dir}/REPORTS/${dataset}/PGX_ONLY/FRQ", mode:'symlink'
//    input:
//        set val(dataset), val(POP_counts), POPs from Ind_group1.values()
//    output:
//        set val(dataset), file(report) into popFrq_combine_count
//    script:
//        report  = "${dataset}_polymophism_count.csv"
//        """
//        cat ${POP_counts} > ${report}
//        """
//}
//
//popFrq_3_1_all.into{ popFrq_3_1_all; popFrq_3_1_1}
////// Add to group data array that contains group:vcf,sapmle,annotatedDatapopulation groups
//popFrq_3_1_1.toSortedList().val.each{ dataset, dataset_dataset, dataset_datasetAnnotated ->
//    group_pop_data[dataset] << dataset_datasetAnnotated
//}
//
//
//'''
//Step 12.4: Plot total polymosphism per population
//'''
//process plot_PolymorphicPerPopulation{
//    tag "plot_poly_${dataset}"
//    label "small"
//    publishDir "${params.work_dir}/REPORTS/${dataset}/PGX_ONLY/FRQ/", mode:'copy'
//    input:
//        set val(dataset), file(poly_count) from popFrq_combine_count
//    output:
//        set val(dataset), file(dataset_pgxPolymorphicPerPopulation) into plot_PolymorphicPerPopulation
//    script:
//        dataset_sample = group_pop_data[dataset][1]
//        dataset_pgxPolymorphicPerPopulation = "${dataset}_pgxPolymorphicPerPopulation.tiff"
//        template "step12_4_polymorphic_per_population.R"
//}
//
//plot_PolymorphicPerPopulation.into{ plot_PolymorphicPerPopulation; plot_PolymorphicPerPopulation_sub}
//
//
//'''
//Step 12.5: List singleton markers
//'''
//popFrq_1_all.into { popFrq_1_all; popFrq_1__popFrq_4 }
//process popFrq_4 {
//    tag "popFreq_4_${dataset}"
//    label "small"
//    publishDir "${params.work_dir}/REPORTS/${dataset}/PGX_ONLY/FRQ", mode:'symlink'
//    input:
//        set val(dataset), file(dataset_biall_vcf), file(dataset_sample), file(dataset_populations), file(dataset_superPopulations) from popFrq_1__popFrq_4
//    output:
//        set val(dataset), file(dataset_sample), file("${dataset}.singletons.per.sample") into popFrq_4_all
//    script:
//        """
//        # List singleton markers
//        vcftools --gzvcf ${dataset_biall_vcf} --singletons --out ${dataset}
//        # Number of singletons per sample
//        awk '{print \$1" "\$2}' ${dataset_sample} > ${dataset}.txt
//        echo "singletons" > ${dataset}.singleton.count
//        while read SAMPLE; do
//            grep \$SAMPLE ${dataset}.singletons | wc -l;
//        done < ${dataset}.txt >> ${dataset}.singleton.count
//        echo "sample pop\n \$(cat ${dataset}.txt)" > ${dataset}.txt
//        paste ${dataset}.txt ${dataset}.singleton.count > ${dataset}.singletons.per.sample
//        rm ${dataset}.singleton.count
//        """
//}
//popFrq_4_all.into{popFrq_4_all; popFrq_4__sub}
//
//
//'''
//Step 12.6: Plot average singletons per population
//'''
//popFrq_4_all.into{ popFrq_4_all; popFrq_4_1 }
//process plot_pgxAverageSingletonsPerPopulation{
//    tag "plot_singl_${dataset}"
//    label "small"
////    publishDir "${params.work_dir}/data/PGX_ONLY/FRQ/${dataset}", mode:'copy'
//    publishDir "${params.work_dir}/REPORTS/${dataset}/PGX_ONLY/FRQ/", mode:'copy'
//    input:
//        set val(dataset), file(dataset_sample), file(dataset_singletons_per_sample) from popFrq_4_1
//    output:
//        set val(dataset), file(dataset_pgxAverageSingletonsPerPopulation) into plot_pgxAverageSingletonsPerPopulation
//    script:
//        dataset_pgxAverageSingletonsPerPopulation = "${dataset}_pgxAverageSingletonsPerPopulation.tiff"
//        template "step12_6_average_singletons_per_population.R"
//}
//
//// """
//// Step 12.6: Allele frequency comparison
//// """
//// Channel
////     .from(mafs_annotations_data)
////     .groupTuple()
////     .set{frequency_comparison_data}
////
//// // process combine_frq_dataset{
//// //     tag "combine_frq_${dataset}"
//// //     label "medmem"
//// //     input:
//// //         set dataset, chrm, dataset_frqs from frequency_comparison_data
//// //     output:
//// //         set dataset, file(dataset_combined_frq) into combine_frq_dataset
//// //     script:
//// //         dataset_combined_frq = "${dataset}.frq"
//// //         """
//// //         head -n 1 ${dataset_frqs[0]} | sed 's/[\\.#]//g' > ${dataset_combined_frq}
//// //         tail -q -n +2 ${dataset_frqs.join(' ')} >> ${dataset_combined_frq}
//// //         """
//// // }
//
//"""
//Step 12.7:
//"""
//mafs_annotations_afr_data = []
//mafs_annotations.each{ dataset ->
//        mafs_annotations_afr_data << [dataset.key, file(dataset.value)]
//}
//mafs_annotations_afr_data_cha = Channel.from(mafs_annotations_afr_data)
//process filter_frq_dataset{
//    tag "filter_frq_${dataset}"
//    label "medmem"
//    input:
//        set dataset, file(frq_file), dataset_AIBST, file(AIBST_frq) from mafs_annotations_afr_data_cha.combine(dataset_pgx_summary_12_7)
//    output:
//        set val("${dataset_AIBST}_${dataset}"), file(frq_output) into filter_frq_dataset
//    script:
//        frq_output = "${dataset_AIBST}.${dataset}.pgx.frq"
//        frq_to_filter = AIBST_frq
//        dataset_to_filter = dataset_AIBST
//        dataset_frq_to_filter = dataset
//        template "filter_frq.py"
//}
//
//
//'''
//Step 12.8: Plot Allele frequency comparison
//'''
//process plot_AlleleFreqComparison{
//    tag "plot_AlleleFreqComparison_${datasets}"
//    label "small"
//    publishDir "${params.work_dir}/REPORTS/${dataset1}/PGX_ONLY/FRQ/", mode:'copy'
//    input:
//        set val(datasets), file(datasets_frq) from filter_frq_dataset
//    output:
//        set val(datasets), file(ExternalAlleleFreqComparison_tiff) into plot_AlleleFreqComparison
//    script:
//        ExternalAlleleFreqComparison_tiff = "${datasets}_ExternalAlleleFreqComparison.tiff"
//        dataset1 = datasets.split('_')[0]
//        dataset2 = datasets.split('_')[1]
//        template "step12_8_plot_allele_frequency_comparison.R"
//}
//
//'''
//Step 13: FST analysis
//'''
//group_vcf_with_pop_sample_cha = Channel.from( group_vcf_with_pop_sample )
//
//// Building dataflow for dataset, annotated VCF, sample file
//concat_dataset_pgx.into { concat_dataset_pgx; concat_dataset_pgx_3 }
//pop_samples_per_group = []
//pop2_per_group = []
//split_POP_samples.into{split_POP_samples; split_POP_samples_2}
//split_POP_samples_list = split_POP_samples_2.toSortedList().val
//pop2_per_group_cha = concat_dataset_pgx_3.flatMap { dataset, dataset_sample, dataset_vcf, dataset_vcf_tbi ->
//    if ( dataset in dataset_files_annot){
//        pop_samples_per_group << [dataset, dataset_vcf, dataset_sample]
//        tmp = []
//        for (data1 in split_POP_samples_list){
//            POP1 = data1[0]
//            SAMPLE1 = data1[1]
//            for (data2 in split_POP_samples_list) {
//                POP2 = data2[0]
//                SAMPLE2 = data2[1]
//                if ( POP1 != POP2 && !([POP2, POP1] in tmp) ){
//                    pop2_per_group << [ dataset, [POP1, POP2], file(SAMPLE1), file(SAMPLE2), dataset_vcf, dataset_sample ]
//                    tmp << [POP1, POP2]
//                }
//            }
//        }
//    }
//    return pop2_per_group
//}
//
//
///// Have to do FST analysis for all variants and compare with PGx
//process fst_analysis {
//    //echo true
//    tag "fst_${POPS[0]}_${POPS[1]}"
//    label "small"
////    publishDir "${params.work_dir}/data/PGX_ONLY/FST/${dataset}", mode:'symlink'
//    publishDir "${params.work_dir}/REPORTS/${dataset}/PGX_ONLY/FST", mode:'copy', pattern:'*.tiff'
//    input:
//        set val(dataset), val(POPS), file(pop1_sample), file(pop2_sample), file(dataset_vcf), file(dataset_sample) from pop2_per_group_cha
//    output:
//        set val(dataset), val(POPS), file("${fst_basename}.weir.fst"), file(fst_log) into fst_analysis_all
//    script:
//        fst_basename = "${dataset}.${POPS[0]}_${POPS[1]}"
//        fst_log = "${dataset}.${POPS[0]}_${POPS[1]}.fst.log"
//        """
//        ## Filter for synonymous variants
//        zcat ${dataset_vcf} | snpsift filter "ANN[*].EFFECT has 'synonymous_variant'" | bcftools view --min-ac 2 | bgzip -c > ${dataset}.fst.syn.vcf.gz
//        ## Compute Fst for pair of populations, here the log file will be used to extract the weighted mean fst
//        vcftools \
//            --gzvcf ${dataset}.fst.syn.vcf.gz \
//            --weir-fst-pop ${pop1_sample} \
//            --weir-fst-pop ${pop2_sample} \
//            --out ${fst_basename} \
//            2> ${fst_log}
//        """
//}
//
//
//'''
//Step 13.1: Combine FST results
//'''
//fst_analysis_all.into { fst_analysis_all; fst_analysis__combine_fst }
//combine_fst_data = [:]
//combine_weir_fst_data = [:]
//fst_analysis_list = fst_analysis__combine_fst.toSortedList().val
//fst_analysis_list.each{ dataset, pop2, weir_fst, fst_log ->
//    if ( dataset in combine_fst_data.keySet() ){
//        pop = []
//        pop << pop2.join('__')
//        pop << weir_fst
//        pop << fst_log
//        combine_fst_data[dataset] << pop.join('__')
//    }
//    else {
//        pop = []
//        pop << pop2.join('__')
//        pop << weir_fst
//        pop << fst_log
//        combine_fst_data[dataset] = [pop.join('__')]
//    }
//}
//
//'''
//Step 13.2
//'''
//// Combine pop fst
//combine_fst_data_cha = Channel.from( combine_fst_data.keySet() )
//process combine_fst_analysis {
//    tag "comb_fst_${dataset}"
//    label "small"
//    publishDir "${params.work_dir}/data/PGX_ONLY/FST/${dataset}", mode:'symlink'
//    input:
//        val(dataset) from combine_fst_data_cha
//    output:
//        set val(dataset), file(fst_out) into combine_fst_analysis,combine_fst_analysis_1,combine_fst_analysis_2
//    script:
//        datas = combine_fst_data[dataset].join(' ')
//        fst_out = "${dataset}.weighted.fst.estimates.txt"
//        """
//        for data in ${datas};do
//            ## Split string into array
//            data_array=(\${data//__/ })
//            pop1=\${data_array[0]}
//            pop2=\${data_array[1]}
//            log_=\${data_array[3]}
//            # Get mean Fst estimates for each of the comparisons
//            weir_fst=`grep "Weir and Cockerham weighted Fst estimate" \$log_ | awk -F":" '{print \$NF}'`
//            echo \"\$pop1 \$pop2 \$weir_fst\" >> ${fst_out}
//        done
//        """
//}
//
//
//'''
//Step 13.3
//'''
//// combine SNP weir fst
//combine_weir_fst_data_cha = Channel.from( combine_fst_data.keySet() )
//process combine_weir_fst_analysis {
//    tag "comb_weir_fst_${dataset}"
//    label "small"
//    publishDir "${params.work_dir}/REPORTS/${dataset}/PGX_ONLY/FST", mode:'symlink', pattern: '*csv'
//
//    input:
//        val(dataset) from combine_weir_fst_data_cha
//    output:
//        set val(dataset), file("${fst_out}.csv"), file("${fst_out}_HighDiff.csv"), file("${fst_out}_HighDiff_pos.csv"), file("${fst_out}_HighDiff.vcf.gz"), file("${fst_out}_HighDiff.ann.csv") into combine_weir_fst_analysis
//    script:
//        datas = combine_fst_data[dataset].join(' ')
//        dataset_vcf = file(group_pop_data[dataset][0])
//        fst_out = "${dataset}.weighted_weir-fst_estimates"
//        """
//        # Select columns where at least one of the Fst comparisons is highly differentiated (i.e. >0.5)
//        python2.7 ${params.homedir}/templates/combine_weir_fst.py  \
//            --fst_input '${datas}' \
//            --fst_output ${fst_out}
//        vcftools \
//            --gzvcf ${dataset_vcf} \
//            --positions ${"${fst_out}_HighDiff_pos.csv"} \
//            --recode-INFO-all --recode --stdout \
//            | gzip -c > ${"${fst_out}_HighDiff.vcf.gz"}
//        zcat ${fst_out}_HighDiff.vcf.gz \
//            | snpsift extractFields - -e "." CHROM POS ID REF ALT AC AF MAF "ANN[0].GENE" AGVP_AF KG_AFR_AF ExAC_AFR_AF gnomAD_AFR_AF "ANN[0].EFFECT" CLNDN CDS GWASCAT_TRAIT GWASCAT_P_VALUE GWASCAT_PUBMED_ID\
//            | awk '{print \$1":"\$2"\\t"\$0}' > ${fst_out}_HighDiff.ann.csv
//        """
//}
//
//'''
//Step 13.4: Plot total FST matrix
//'''
//process plot_fst_matrix{
//    tag "plot_fst_matrix_${dataset}"
//    label "small"
//    publishDir "${params.work_dir}/REPORTS/${dataset}/PGX_ONLY/FST/", mode:'copy'
//    input:
//        set val(dataset), file(weighted_fst_estimates) from combine_fst_analysis_1
//    output:
//        set val(dataset), file(weighted_fst_estimates_png) into plot_fst_matrix
//    script:
//        weighted_fst_estimates_png = "${dataset}_weighted_fst_estimates_png.tiff"
//        template "step13_4_plot_fst_matrix.R"
//}
//
//
//'''
//Step 13.5: Filter dataset VCF for highly differentiated SNPs
//'''
//combine_weir_fst_analysis.into { combine_weir_fst_analysis; combine_weir_fst_analysis_1 }
//process plot_highDiff {
//    tag "plot_highDiff_${dataset}"
//    publishDir "${params.work_dir}/REPORTS/${dataset}/PGX_ONLY/FST/", mode:'copy'
//    input:
//        set val(dataset), file(fst_csv), file(fst_HighDiff_csv), file(fst_HighDiff_pos),  file(fst_HighDiff_vcf),  file(fst_HighDiff_ann) from combine_weir_fst_analysis_1
//    output:
//        set val(dataset), file(HighGene_out), file(HighFst_plot) into vcf_highDiff_weir_fsf
//    script:
//        dataset_vcf = file(group_pop_data[dataset][0])
//        datasetAnnotated = file(group_pop_data[dataset][2])
//        HighGene_out = "${dataset}.highFstMaxGenes.csv"
//        HighFst_plot = "${dataset}.pgxHighDiffPGxPlot.tiff"
//        template "step13_5_highly_diff_variants.R"
//}
//
//
//
//'''
//Step 13.2: transform FST results to matrix
//'''
//process matrix_combine_fst_analysis {
//    tag "matrix_comb_fst_${dataset}"
//    label "small"
//    publishDir "${params.work_dir}/data/PGX_ONLY/FST/${dataset}", mode:'copy'
//    publishDir "${params.work_dir}/REPORTS/${dataset}/PGX_ONLY/FST/", mode:'copy'
//
//    input:
//        set val(dataset), file(fst_in) from combine_fst_analysis_2
//    output:
//        set val(dataset), file(fst_out) into matrix_combine_fst_analysis
//    script:
//        fst_out = "${dataset}.weighted.fst.estimates.matrix.csv"
//        """
//        python2.7 ${params.homedir}/templates/convert_fst_result_to_matrix.py \
//            --fts_2by2_input ${fst_in} \
//            --fst_matrix_output ${fst_out}
//        """
//}
//
//
//"""
//14.0: Generate clinical vcf file
//"""
//clinVarData_all.into { clinVarData_all; clinVarData_all_2 }
//process clinical_dataset{
//    tag "clinical_${dataset}"
//    label "bigmem"
//    publishDir "${params.work_dir}/data/${dataset}/ALL/VCF", mode:'symlink', pattern: "*.vcf.gz"
//    publishDir "${params.work_dir}/REPORTS/${dataset}/CLINICAL", mode:'symlink', pattern: "*.INFO"
//    input:
//        set val(dataset), file(dataset_sample), file(vcf_file), file(vcf_tbi) from concat_dataset_all_14_0
//        file(pgxClinicalrsID) from clinVarData_all_2
//    output:
//        set val(dataset), file(dataset_sample), file("${dataset}.clinical.vcf.gz"), file("${dataset}.clinical.vcf.gz.tbi"), file(pgxClinicalrsID) into clinical_dataset
//        set val(dataset), file(dataset_sample), file("${dataset}.clinical.vcf.gz"), file("${dataset}.clinical.INFO"), file("${dataset}.ClinPerSample.012"), file("${dataset}.ClinPerSample.012.indv"), file("${dataset}.ClinPerSample.012.pos"), file(pgxClinicalrsID) into clinical_dataset_012,clinical_dataset_012_143,clinical_dataset_012_144
//        file("${dataset}.clinical.INFO") into clinical_dataset_info
//    script:
//        """
//        # Extract clinical sites for population
//        bcftools \
//           view \
//           --targets-file ${pgxClinicalrsID} \
//           ${vcf_file} \
//           -Oz -o ${dataset}.clinical_.vcf.gz
//        ## Recalculate AC, AN, AF
//        bcftools +fill-tags ${dataset}.clinical_.vcf.gz -Oz -o ${dataset}.clinical.vcf.gz
//        bcftools index --tbi -f ${dataset}.clinical.vcf.gz
//        vcftools --gzvcf ${dataset}.clinical.vcf.gz \
//            --012 --out ${dataset}.ClinPerSample
//        zcat ${dataset}.clinical.vcf.gz | snpsift extractFields - -e "." \
//            CHROM POS ID REF ALT AC ANN[0].GENE MAF KG_AF KG_AFR_AF KG_EUR_AF KG_AMR_AF KG_EAS_AF gnomAD_AF gnomAD_AFR_AF gnomAD_FIN_AF ExAC_AF ExAC_AFR_AF AGVP_AF SAHGP_AF "ANN[0].EFFECT" CLNDN CDS GWASCAT_TRAIT GWASCAT_P_VALUE GWASCAT_PUBMED_ID\
//            > ${dataset}.clinical.INFO
//        rm -f ${dataset}.clinical_.vcf.gz
//        """
//}
//
//
//'''
//Step 14.1: Clinical variants analysis. Generate VCF of clinical variants only per population
//'''
//group_vcf_with_pop_sample_cha.into { group_vcf_with_pop_sample_cha; group_vcf_with_pop_sample_cha__clinVar }
//clinVarData_all.into { clinVarData_all; clinVarData__clinVar }
//process clinVar {
//    tag "clinVar_${POP}"
//    label "small"
//    publishDir "${params.work_dir}/data/PGX_ONLY/CLINICAL/${GROUP_POP}", mode:'symlink'
//
//    input:
//        set val(POP), val(GROUP_POP), file(GROUP_POP_vcf), file(GROUP_POP_sample) from group_vcf_with_pop_sample_cha__clinVar
//        file(pgxClinicalrsID) from clinVarData__clinVar
//    output:
//        set val(POP), val(GROUP_POP), file("${POP}.clinical.vcf.gz"), file("${POP}.clinical.ind.frq"), file("${POP}.rsIndexClin.txt") into clinVar_all
//    script:
//        """
//        grep ${POP} ${GROUP_POP_sample} | cut -f1 > ${POP}.samples
//        # Select rsIDs from file and use this list to search annotated vcf
//        # Extract sites for population
//        vcftools --gzvcf ${GROUP_POP_vcf} \
//            --keep ${POP}.samples \
//            --recode-INFO-all \
//            --recode --stdout | \
//            bgzip -c > ${POP}.tmp.clinical.vcf.gz
//        bcftools \
//           view \
//           --targets-file ${pgxClinicalrsID} \
//           ${POP}.tmp.clinical.vcf.gz \
//           -Oz -o ${POP}.clinical.vcf.gz
//        # Generate sites frequency for population
//        vcftools --gzvcf ${POP}.clinical.vcf.gz \
//            --freq2 \
//            --out ${POP}.clinical.ind
//        # Create an index from vcf file that can be used to map back rsIDs as freq only outputs coordinates
//        zcat ${POP}.clinical.vcf.gz | awk '!/^(\$|#)/ {print \$3"\t"\$1":"\$2}' - > ${POP}.rsIndexClin.txt
//        rm ${POP}.tmp.clinical.vcf.gz
//        """
//}
//
//
/////
//////'''
//////Step 14.2: Clinical variants analysis. Calculate the number of non-ref alleles per sample for the whole dataset (AIBST)
//////'''
//////clinVarData_all.into { clinVarData_all; clinVarData__clinVar1 }
//////process clinVar1 {
//////    tag "clinVar1_${dataset}"
//////    label "small"
//////    publishDir "${params.work_dir}/data/PGX_ONLY/CLINICAL/${dataset}", mode:'symlink'
//////    input:
//////        val(dataset) from dataset_files_annot
//////        file(pgxClinicalrsID) from clinVarData__clinVar1
//////    output:
//////        set val(dataset), file("${dataset}.clinical.vcf.gz"), file("${dataset}.ClinPerSample.012"), file("${dataset}.ClinPerSample.012.indv"), file("${dataset}.ClinPerSample.012.pos"), file("pgxClinicalrsID.txt") into clinVar1_all
//////    script:
//////        dataset_vcf = file(group_pop_data[dataset][0])
//////        """
//////        # Calculate the number of non-ref alleles per sample
//////        vcftools --gzvcf ${dataset_vcf} \
//////            --snps ${pgxClinicalrsID} \
//////            --recode-INFO-all --recode --stdout | \
//////            bgzip -c > ${dataset}.clinical.vcf.gz
//////        vcftools --gzvcf ${dataset}.clinical.vcf.gz \
//////            --012 --out ${dataset}.ClinPerSample
//////        """
//////}
//
//'''
//Step 14.3: Plot Clinical variants analysis. (AIBST)
//'''
//process clinVar_plot {
//    tag "clinVar_plot_${dataset}"
//    label "small"
//    publishDir "${params.work_dir}/REPORTS/${dataset}/PGX_ONLY/CLINICAL", mode:'copy'
//    input:
//        set val(dataset), file(dataset_sample_file), file(dataset_clinical_vcf), file(dataset_clinical_INFO), file(ClinPerSample_012_file), file(ClinPerSample_012_indv_file), file(ClinPerSample_012_pos_file), file(pgxClinicalrsID_txt) from clinical_dataset_012_143
//    output:
//        set val(dataset), file(pgxClinPerSamplePerPopulation_tiff_file), file(clinSampleTable_csv_file), file(clinPopTable_csv_file) into clinVar_plot
//    script:
//        dataset_allLofPop_file = "${params.work_dir}/data/PGX_ONLY/LOF/${dataset}/${dataset}_allLofPop.csv"
//        clinSampleTable_csv_file = "${file(dataset_clinical_vcf.baseName).baseName}_clinSampleTable.csv"
//        clinPopTable_csv_file = "${file(dataset_clinical_vcf.baseName).baseName}_clinPopTable.csv"
//        pgxClinPerSamplePerPopulation_tiff_file = "${file(dataset_clinical_vcf.baseName).baseName}_pgxClinPerSamplePerPopulation.tiff"
//        template "step14_3_clinical_variants.R"
//}
//
//
//'''
//Step 14.4: Plot Clinical variants per population (AIBST pops, 1KG, gnomAD)
//'''
//process clinVar_plot_pop {
//    tag "clinVar_pop_plot_${dataset}"
//    label "small"
//    publishDir "${params.work_dir}/REPORTS/${dataset}/PGX_ONLY/CLINICAL", mode:'copy'
//    input:
//        set val(dataset), file(dataset_sample_file), file(dataset_clinical_vcf), file(dataset_clinical_INFO), file(ClinPerSample_012_file), file(ClinPerSample_012_indv_file), file(ClinPerSample_012_pos_file), file(pgxClinicalrsID_txt) from clinical_dataset_012_144
//    output:
//        set val(dataset), file(tiff_file), file(tsv_file) into clinVar_plot_pop
//    script:
//        datasetAnnotated = group_pop_data[dataset][2]
//        pgxClinicalLevel1 = dataset_clinical_INFO
//        tiff_file = "${file(dataset_clinical_vcf.baseName).baseName}_pgxClinicalEvidencePopulation.tiff"
//        tsv_file = "${file(dataset_clinical_vcf.baseName).baseName}_pgxClinicalEvidencePopulation.tsv"
//        template "step14_4_clinical_variants_pop.R"
//}
//
//
//
//
//'''
//Step 15.1: Rare Global Common Variants
//'''
//concat_dataset_pgx.into { concat_dataset_pgx; concat_dataset_pgx_4 }
//process rareGlobCommon {
//    tag "rareGlobCommon_${dataset}"
//    label "small"
//    publishDir "${params.work_dir}/data/PGX_ONLY/RARE_GLOB_COM/${dataset}", mode:'symlink'
//
//    input:
//        set val(dataset), file(dataset_sample), file(dataset_vcf), file(dataset_vcf_tbi) from concat_dataset_pgx_4
//    output:
//        set val(dataset), file(dataset_sample), file("${dataset}.rare.vcf.gz"), file("${dataset}.rare.IDs"), file("${dataset}.rare.pos") into rareGlobCommon
//    script:
//        """
//        # Make a vcf with all rare variants (0.05%) in the global dataset
//        vcftools --gzvcf ${dataset_vcf} \
//            --max-maf 0.01 --recode-INFO-all \
//            --recode --stdout | \
//            bgzip -c > ${dataset}.rare.vcf.gz
//        # Get variant IDs for all rare variants
//        zcat ${dataset}.rare.vcf.gz | grep -v ^# | awk '{print \$3}' > ${dataset}.rare.IDs
//        # Set ID to CHRM_POS_REF_ALT
//        bcftools query \
//            -f '%CHROM\\t%POS\\n' \
//            ${dataset}.rare.vcf.gz \
//            > ${dataset}.rare.pos
//        """
//
//}
//
//
//''''
//Step 15.2: Rare Global Common Variants per population
//'''
//group_vcf_with_pop_sample_cha.into { group_vcf_with_pop_sample_cha; group_vcf_with_pop_sample_1 }
//rareGlobCommon.into { rareGlobCommon; rareGlobCommon_1 }
//rareGlobCommon_list = rareGlobCommon_1.toSortedList().val
//
//rareGlobCommon_pop_all = group_vcf_with_pop_sample_1.flatMap{ POP, dataset, dataset_POP_vcf, dataset_POP_sample ->
//    group_vcf_with_pop_list = []
//    rareGlobCommon_list.each{ dataset_, dataset_sample, dataset_rare_vcf, dataset_rare_ids, dataset_rare_pos  ->
//        if (dataset == dataset_){
//            group_vcf_with_pop_list << [POP, dataset, dataset_rare_vcf, dataset_sample, dataset_rare_ids, dataset_rare_pos]
//        }
//    }
//    return group_vcf_with_pop_list
//}
//
//process rareGlobCommon_pop {
//    tag "rareGlobCommon_pop_${POP}"
//    label "small"
//    publishDir "${params.work_dir}/data/PGX_ONLY/RARE_GLOB_COM/${dataset}", mode:'symlink'
//
//    input:
//        set val(POP), val(dataset), file(dataset_rare_vcf), file(dataset_sample), file(dataset_rare_ids), file(dataset_rare_pos) from rareGlobCommon_pop_all
//    output:
//        set val(POP), val(dataset), file("${POP}.pop.diff.vcf.gz") into rareGlobCommon_pop
//    script:
//        """
//        grep ${POP} ${dataset_sample} | cut -f1 > ${POP}.samples
//        # Make a vcf for each population for variants that occur at 5% in that population
//        # And are rare in global dataset
//        vcftools --gzvcf ${dataset_rare_vcf} \
//            --keep ${POP}.samples \
//            --snps ${dataset_rare_ids} \
//            --maf 0.01 \
//            --recode-INFO-all \
//            --recode --stdout | \
//        bgzip -c > ${POP}.pop.diff.vcf.gz
//        """
//}
//
//
//
//''''
//Step 15.3: Rare Global Common Variants per population
//'''
//rareGlobCommon_pop.into { rareGlobCommon_pop; rareGlobCommon_pop_1 }
//rareGlobCommon_pop_combine_all = [:]
//rareGlobCommon_pop_combine_pop_all = [:]
//rareGlobCommon_pop_1.toSortedList().val.each{ POP, dataset, POP_diff_vcf ->
//    if ( !(dataset in rareGlobCommon_pop_combine_all.keySet()) ){
//        rareGlobCommon_pop_combine_all[dataset] = []
//        rareGlobCommon_pop_combine_pop_all[dataset] = []
//    }
//    // Put all population VCFs in list for group_pop
//    rareGlobCommon_pop_combine_all[dataset] << POP_diff_vcf
//    rareGlobCommon_pop_combine_pop_all[dataset] << [POP_diff_vcf, POP].join('__')
//}
//rareGlobCommon_pop_combine_cha = Channel.from( rareGlobCommon_pop_combine_pop_all.keySet() )
//
//process rareGlobCommon_pop_combine {
//    tag "rareGlobCommon_pop_combine_${dataset}"
//    label "small"
//    publishDir "${params.work_dir}/data/PGX_ONLY/RARE_GLOB_COM/${dataset}", mode:'symlink'
//
//    input:
//        val(dataset) from rareGlobCommon_pop_combine_cha
//    output:
//        set val(dataset), file(dataset_diff_snps_vcf), file(dataset_diff_snps_genotypes), file("${dataset}.population.pop.diff.snps.genotypes"), file("${dataset}.rsIndexPopDiff.txt"), file("${dataset}.genes.pop.diff.snps.genotypes"), file("${dataset}.population.gene.pop.diff.snps.genotypes") into rareGlobCommon_pop_combine_cha_all
//        file("*pop.diff.snps.Ind.frq") into rareGlobCommon_pop_combine_cha_all_1
//    script:
//        dataset_diff_snps_vcf = "${dataset}.diff.snps.vcf.gz"
//        dataset_diff_snps_genotypes = "${dataset}.diff.snps.genotypes"
//        """
//        # Extract SNPs that are differentiated
//        data_array=\"${rareGlobCommon_pop_combine_pop_all[dataset].join(' ')}\"
//        for POP_VCF in \$data_array;do
//            POP_VCF_array=(\${POP_VCF//__/ })
//            vcf=\${POP_VCF_array[0]}
//            pop=\${POP_VCF_array[1]}
//            zcat \$vcf | awk '!/^(\$|#)/ {print \$0" ""'"\$pop"'"}' >> ${dataset}.diff.snps.genotypes
//            # Get variant IDs for all differentiated variants
//            zcat \$vcf | grep -v ^# | awk '{print \$3}' >> ${dataset}.diff.snps.rsIDs
//            zcat \$vcf | awk '!/^(\$|#)/ {print \$3" "\$1":"\$2":"\$5}' >> ${dataset}.rsIndexPopDiff.txt
//            #Extract gene name and variant type
//            zcat \$vcf | snpsift extractFields - -e "." ANN[0].GENE ANN[0].EFFECT | \
//                awk 'NR>1 {print \$0}' >> ${dataset}.temp.genes.pop.diff.snps.genotypes
//        done
//        # Generate new vcf
//        vcftools --gzvcf ${group_pop_data[dataset][0]} \
//            --recode --snps ${dataset}.diff.snps.rsIDs \
//            --recode-INFO-all --stdout | bgzip -c > ${dataset}.diff.snps.vcf.gz
//        for POP_VCF in \$data_array;do
//            POP_VCF_array=(\${POP_VCF//__/ })
//            vcf=\${POP_VCF_array[0]}
//            pop=\${POP_VCF_array[1]}
//            grep \$pop ${group_pop_data[dataset][1]} | cut -f1 > \$pop.samples
//            vcftools --gzvcf ${dataset}.diff.snps.vcf.gz \
//                --keep \$pop.samples \
//                --freq2 \
//                --out \$pop.pop.diff.snps.Ind
//        done
//        zcat ${dataset}.diff.snps.vcf.gz | awk '!/^(\$|#)/ {print \$3" "\$1":"\$2":"\$5}' > ${dataset}.rsIndexPopDiff.txt
//        # Extract populations
//        awk '{print \$NF}' ${dataset}.diff.snps.genotypes > ${dataset}.temp.population.pop.diff.snps.genotypes
//        # Extract rsID only
//        awk '{print \$3}' ${dataset}.diff.snps.genotypes > ${dataset}.temp.rs.pop.diff.snps.genotypes
//        # Population and rsID
//        paste ${dataset}.temp.population.pop.diff.snps.genotypes ${dataset}.temp.rs.pop.diff.snps.genotypes > ${dataset}.population.pop.diff.snps.genotypes
//        # Extract population only
//        awk '{print \$1}' ${dataset}.population.pop.diff.snps.genotypes > ${dataset}.population.only.pop.diff.snps.genotypes
//        #Extract gene name and variant type
//        paste ${dataset}.temp.genes.pop.diff.snps.genotypes ${dataset}.rsIndexPopDiff.txt > ${dataset}.genes.pop.diff.snps.genotypes
//        paste ${dataset}.population.only.pop.diff.snps.genotypes ${dataset}.genes.pop.diff.snps.genotypes > ${dataset}.population.gene.pop.diff.snps.genotypes
//        """
//}
//
//
//''''
//Step 16.1: Consequences
//'''
//concat_dataset_pgx.into { concat_dataset_pgx; concat_dataset_pgx_5 }
//process csq {
//   tag "csq_${dataset}"
//   label "small"
//   publishDir "${params.work_dir}/data/PGX_ONLY/CSQ/${dataset}", mode:'symlink'
//   publishDir "${params.work_dir}/REPORTS/${dataset}/PGX_ONLY/CSQ", mode:'copy', pattern: "*_summaryVariantsPerGene.tsv"
//   input:
//       set val(dataset), file(dataset_sample), file(dataset_vcf), file(dataset_vcf_tbi), file(so_term_file) from concat_dataset_pgx_5.combine([file(params.so_term_file)])
//   output:
//       set val(dataset), file(dataset_sample), file("${dataset}.singletons.vcf.gz"), file("${dataset}_SO.terms.MAF.summary"), file("${dataset}.singletons.uniqueID"), file("${dataset}.nonsing-0.01MAF.uniqueID"), file("${dataset}.0.01-0.05MAF.uniqueID"), file("${dataset}.greater.0.05.uniqueID") into csq
//   script:
//       """
//       # Create different vcfs for different frequency classes
//       # Singletons
//       vcftools --gzvcf ${dataset_vcf} \
//           --max-mac 1 --recode-INFO-all \
//           --recode --stdout | \
//           bgzip -c > ${dataset}.singletons.vcf.gz
//       # Non-singletons to MAF 0.01
//       vcftools --gzvcf ${dataset_vcf} \
//           --recode --mac 2 --max-maf 0.01 \
//           --recode-INFO-all --stdout | \
//           bgzip -c > ${dataset}.nonsing-0.01MAF.vcf.gz
//       # MAF 0.01-0.05
//       vcftools --gzvcf ${dataset_vcf} \
//           --recode --maf 0.01 --max-maf 0.05 \
//           --recode-INFO-all --stdout | \
//           bgzip -c > ${dataset}.0.01-0.05MAF.vcf.gz
//       # MAF > 0.05
//       vcftools --gzvcf ${dataset_vcf} \
//           --recode --maf 0.05 \
//           --recode-INFO-all --stdout | \
//           bgzip -c > ${dataset}.greater.0.05MAF.vcf.gz
//
//       # Count number of appearances of each of the Sequence Ontology (SO) terms
//       # so_term_file was modified from:
//       # https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
//       while read SO; do
//           zcat ${dataset_vcf} | grep \$SO | wc -l
//       done < ${so_term_file} > ${dataset}_allCount
//       while read SO; do
//           zcat ${dataset}.singletons.vcf.gz | grep \$SO | wc -l
//       done < ${so_term_file} > ${dataset}_singletonsCount
//       while read SO; do
//           zcat ${dataset}.nonsing-0.01MAF.vcf.gz | grep \$SO | wc -l
//       done < ${so_term_file} > ${dataset}_nonsing-0.01MAFCount
//       while read SO; do
//           zcat ${dataset}.0.01-0.05MAF.vcf.gz | grep \$SO | wc -l
//       done < ${so_term_file} > ${dataset}_0.01-0.05MAFCount
//       while read SO; do
//           zcat ${dataset}.greater.0.05MAF.vcf.gz | grep \$SO | wc -l
//       done < ${so_term_file} > ${dataset}_greater0.05MAFCount
//       echo "terms 	countAll	countSingletons	countnonSing-0.01MAF	count0.01-0.05MAF	countGreater0.05" \
//           > ${dataset}_SO.terms.MAF.summary
//       paste ${so_term_file} ${dataset}_allCount ${dataset}_singletonsCount  ${dataset}_nonsing-0.01MAFCount ${dataset}_0.01-0.05MAFCount ${dataset}_greater0.05MAFCount \
//           >> ${dataset}_SO.terms.MAF.summary
//
//       # Create UniqueID for each of the vcf frequency classes
//       zcat ${dataset}.singletons.recode.vcf.gz | grep -v ^# | awk '{print \$1":"\$2":"\$5}' \
//           >  ${dataset}.singletons.uniqueID
//       zcat ${dataset}.nonsing-0.01MAF.recode.vcf.gz | grep -v ^# | awk '{print \$1":"\$2":"\$5}' \
//           >  ${dataset}.nonsing-0.01MAF.uniqueID
//       zcat ${dataset}.0.01-0.05MAF.recode.vcf.gz | grep -v ^# | awk '{print \$1":"\$2":"\$5}' \
//           >  ${dataset}.0.01-0.05MAF.uniqueID
//       zcat ${dataset}.greater.0.05MAF.recode.vcf.gz | grep -v ^# | awk '{print \$1":"\$2":"\$5}' \
//           >  ${dataset}.greater.0.05.uniqueID
//
//       """
//}
//
//Channel
//        .fromPath( params.final_all_pharmacogenes_details )
//        .splitCsv()
//        .set{ gene_list }
//
//
//process gene_vcfs {
//   tag "gene_vcfs_${dataset}_${gene}"
//   label "small"
//   input:
//       set val(dataset), file(dataset_sample), file(dataset_vcf), file(dataset_vcf_tbi), gene from concat_dataset_pgx_gene_stats.combine( gene_list )
//   output:
//       set dataset, gene into gene_vcfs
//   script:
//        """
//        # Generate vcf for each of the genes in the study independently
//        snpsift filter "ANN[*].GENE = '${gene}'" ${dataset_vcf} -v | bgzip -c > ${dataset}_${gene}.vcf.gz
//        """
//}
//
//''''
//Step 16.2: Plot consequences
//'''
//csq.into{ csq; csq_1 }
//process plot_csq {
//   tag "plot_csq_${dataset}"
//   label "small"
////    publishDir "${params.work_dir}/data/PGX_ONLY/CSQ/${dataset}", mode:'symlink'
//   publishDir "${params.work_dir}/REPORTS/${dataset}/PGX_ONLY/CSQ", mode:'copy'
//   input:
//       set val(dataset), file(dataset_sample), file(dataset_singletons_vcf_gz), file(dataset_SO_terms_MAF_summary), file(dataset_singletons_uniqueID), file(dataset_nonsing_0_01MAF_uniqueID), file(dataset_0_01_0_05MAF_uniqueID), file(dataset_greater_0_05_uniqueID) from csq_1
//   output:
//       set val(dataset), file(dataset_pgxFunctionalClassesCounts_tiff), file(dataset_pgxFunctionalClassesByFrequency_tiff) into plot_csq
//   script:
//       dataset_SO_terms_MAF_summary = dataset_SO_terms_MAF_summary
//       dataset_pgxFunctionalClassesCounts_tiff = "${dataset}_pgxFunctionalClassesCounts.tiff"
//       dataset_pgxFunctionalClassesByFrequency_tiff = "${dataset}_pgxFunctionalClassesByFrequency.tiff"
//       template "step16_2_consequence.R"
//}
//
//
//''''
//Step 17.1: RareMissense
//'''
////csq_gene_vcf.into { csq_gene_vcf; csq_gene_vcf_1 }
////process RareMissense {
////   tag "RareMissense_${dataset}"
////   label "small"
////   publishDir "${params.work_dir}/REPORTS/${dataset}/PGX_ONLY/RareMissense/${dataset}", mode:'copy'
////   input:
////       set val(dataset), file(gene_vcfs_url) from csq_gene_vcf_1
////   output:
////       set val(dataset), file("${dataset}_rare.missense.gene.txt") into RareMissense
////   script:
////       """
////       echo "rare_missense" > ${gene_vcfs_url}/${dataset}_temp.rare.missense.gene.txt
////       # Analyse each of the genes in the study independently, appending counts
////       while read gene;
////       do
////           vcftools --gzvcf ${gene_vcfs_url}/${dataset}_\$gene.vcf.gz --max-maf 0.01 --recode-INFO-all --recode --stdout | \
////           grep missense_variant | wc -l >> ${gene_vcfs_url}/${dataset}_temp.rare.missense.gene.txt
////       done < ${params.final_all_pharmacogenes_details}
////       paste ${gene_vcfs_url}/${dataset}_symbol.txt ${gene_vcfs_url}/${dataset}_temp.rare.missense.gene.txt > ${dataset}_rare.missense.gene.txt
////       #rm ${gene_vcfs_url}/${dataset}_temp.rare.missense.gene.txt
////       """
////}
//
//
////
////''''
////Step 18.1: VIP Mask
////'''
////concat_dataset_all.into{ concat_dataset_all; concat_dataset_all_18}
////process vipmask{
//// tag "vipmask_${dataset}"
//// cpus { 2 * task.attempt }
//// memory { 2.GB * task.cpus }
//// publishDir "${params.work_dir}/REPORTS/${dataset}/PGX_ONLY/VIPMASK/${dataset}", mode:'copy'
//// input:
////     set val(dataset), file(dataset_sample), file(dataset_vcf), file(dataset_vcf_tbi) from concat_dataset_all_18
//// output:
////     set val(dataset), file("${dataset_out}_maskGeneSummary.tsv") into vipmask
//// script:
////     dataset_out = "${file(dataset_vcf.baseName).baseName}"
////     """
////     ###sed s/^chr//g 20141020.strict_mask.whole_genome.bed > 20141020.strict_mask.whole_genome_final.bed
////
////     # Make vcf of masked/unmasked regions
////     vcftools --gzvcf ${dataset_vcf} \
////        --bed ${params.strict_mask_whole_genome_final} \
////        --recode --recode-INFO-all --stdout \
////        | bgzip -c > ${dataset_out}_accessible.vcf.gz
////     vcftools --gzvcf ${dataset_vcf} \
////        --exclude-bed ${params.strict_mask_whole_genome_final} \
////        --recode --recode-INFO-all --stdout \
////        | bgzip -c > ${dataset_out}_inaccessible.vcf.gz
////
////     # Look at segmental duplications
////     # Make two vcf files, one within segmental dups, another outside of them
////     vcftools --gzvcf ${dataset_vcf} \
////        --bed ${params.GRCh37GenomicSuperDup} \
////        --recode --recode-INFO-all --stdout \
////        | bgzip -c > ${dataset_out}_segmental.vcf.gz
////     vcftools --gzvcf ${dataset_vcf} \
////        --exclude-bed ${params.GRCh37GenomicSuperDup} \
////        --recode --recode-INFO-all --stdout \
////        | bgzip -c > ${dataset_out}_nonsegmental.vcf.gz
////
////
////     # Make a vcf of SNPs found in both segmental duplications and inaccessible regions
////     vcftools --gzvcf ${dataset_out}_segmental.vcf.gz \
////         --exclude-bed ${params.strict_mask_whole_genome_final} \
////         --stdout --recode --recode-INFO-all \
////         | bgzip -c > ${dataset_out}.segmental.inaccessible.vcf.gz
////
////
////     # Calculate the number of variants per gene
////     echo total_variants > ${dataset_out}_totalVariantsPerGene.txt
////     while read gene;
////     do
////         zcat ${params.work_dir}/data/PGX_ONLY/CSQ/${dataset}/GENE/${dataset}_\$gene.vcf.gz \
////             | grep -v ^#  | wc -l >> ${dataset_out}_totalVariantsPerGene.txt
////     done < ${params.final_all_pharmacogenes_details}
////
////
////     # Look at number of variants in inaccessible regions
////     # Calculate the number of variants per gene
////     echo inaccessible_variants > ${dataset_out}_inaccessiblePerGene.txt
////     while read geneIn;
////     do
////         zcat  ${dataset_out}_inaccessible.vcf.gz \
////             | grep -v ^# | grep -w \$geneIn | wc -l >> ${dataset_out}_inaccessiblePerGene.txt
////     done < ${params.final_all_pharmacogenes_details}
////
////
////     # Look at number of variants in segmental duplications
////     # Calculate the number of variants per gene
////     echo segmental_variants > ${dataset_out}_segmentalPerGene.txt
////     while read geneSeg;
////     do
////         zcat ${dataset_out}_segmental.vcf.gz \
////             | grep -v ^# | grep -w \$geneSeg | wc -l >> ${dataset_out}_segmentalPerGene.txt
////     done < ${params.final_all_pharmacogenes_details}
////
////
////     # Look at number of variants in inaccessible and segmental duplications
////     # Calculate the number of variants per gene
////     echo inaccess_segmental_variants > ${dataset_out}_inaccessSegmentalPerGene.txt
////     while read geneBoth;
////     do
////         zcat ${dataset_out}_inaccessible.vcf.gz \
////             |grep -v ^# | grep -w \$geneBoth | wc -l >> ${dataset_out}_inaccessSegmentalPerGene.txt
////     done < ${params.final_all_pharmacogenes_details}
////
////
////     # Add heading to VIP genes gene symbol
////     echo gene > ${dataset_out}_symbol.txt
////     awk '{print \$1}' ${params.final_all_pharmacogenes_details} >> ${dataset_out}_symbol.txt
////
////     # Calculate total length of each gene that was extracted
////     echo geneBedLengthTotal > ${dataset_out}_geneBedLengthTotal.txt
////     while read geneID
////     do
////         grep -w \$geneID ${params.gene_regions_slop} | \
////         sort -k1,1 -k2,2n | bedtools merge -i stdin | \
////         awk -F'\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}' \
////         >> ${dataset_out}_geneBedLengthTotal.txt
////     done < ${params.final_all_pharmacogenes_details}
////
////     # Calculate transcript length of each gene that was extracted
////     echo geneBedLengthTranscript > ${dataset_out}_geneBedLengthTranscript.txt
////     while read geneID
////     do
////         grep -w \$geneID ${params.gene_regions_slop} | \
////         sort -k1,1 -k2,2n | bedtools merge -i stdin | \
////         awk -F'\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}' \
////         >> ${dataset_out}_geneBedLengthTranscript.txt
////     done < ${params.final_all_pharmacogenes_details}
////
////     # Combine all data together
////     paste ${dataset_out}_symbol.txt ${dataset_out}_totalVariantsPerGene.txt ${dataset_out}_inaccessiblePerGene.txt ${dataset_out}_segmentalPerGene.txt ${dataset_out}_inaccessSegmentalPerGene.txt ${dataset_out}_geneBedLengthTranscript.txt ${dataset_out}_geneBedLengthTotal.txt > ${dataset_out}_maskGeneSummary.tsv
////
////     """
////}
//
//
//"""
//SUMMARY
//"""
//
//'''
//All SNVs
//'''
//annotate_dataset_snpeff_org.into{ annotate_dataset_snpeff_org; annotate_dataset_snpeff_org_variant_summary}
//process variant_summary {
//   tag "variant_summary_${dataset}"
//   label "small"
//   echo true
//   publishDir "${params.work_dir}/REPORTS/${dataset}/SUMMARY/ALL", mode:'copy', pattern: "*csv"
//   input:
//       set dataset, file(vcf_file) from annotate_dataset_snpeff_org_variant_summary
//   output:
//       set dataset, file(report) into variant_summary_report
//       set dataset, file("${vcf_file}"), file("${base}_coding.vcf.gz"), file("${base}.singletons") into variant_summary
//   script:
//       base = file(vcf_file.baseName).baseName
//       report = "${base}_all_variants_report.csv"
//       """
//       ## Get total number positions of variants (Exome postions with SNV)
//       echo -e "Total positions ; \$(zcat ${vcf_file} | grep -v "^#" | wc -l)" > ${report}
//
//       ## Total rare variants
//       echo -e "Rare variants ; \$(vcftools --gzvcf ${vcf_file} --max-maf 0.01 --recode-INFO-all --recode --stdout | grep -v "^#" | wc -l)" >> ${report}
//
//       ## Total biallelic variants
//       echo -e "Biallelic positions ; \$(bcftools view -m2 -M2 ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//
//       ## Total biallelic SNV variants
//       echo -e "Biallelic SNV positions ; \$(bcftools view -m2 -M2 -v snps ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//
//       ## Total biallelic INDELS variants
//       echo -e "Biallelic INDELS positions ; \$(bcftools view -m2 -M2 -v indels ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//
//       ## Total multiallelic variants
//       echo -e "Multiallelic positions ; \$(bcftools view -m3 ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//
//       ## Total multiallelic SNV variants
//       echo -e "Multiallelic SNV positions ; \$(bcftools view -m3 -v snps ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//
//       ## Total multiallelic INDELS variants
//       echo -e "Multiallelic INDELS positions ; \$(bcftools view -m3 -v indels ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//
//       ## Total Singleton positions
//       vcftools \
//           --gzvcf ${vcf_file} \
//           --singletons \
//           --out ${base}
//       echo -e "Total Singletons ; \$(grep -v "^CHR" ${base}.singletons| wc -l)" >> ${report}
//
//       ## Nonsynonymous/missense
//       echo -e "Nonsynonymous SNV ; \$(zcat ${vcf_file} | snpsift filter "(ANN[*].EFFECT has 'missense_variant')" | grep -v "^#" | wc -l)" >> ${report}
//       ## Synonymous
//       echo -e "Synonymous SNV ; \$(zcat ${vcf_file} | snpsift filter "(ANN[*].EFFECT has 'synonymous_variant')" | grep -v "^#" | wc -l)" >> ${report}
//       ## Stop gained
//       echo -e "Stop gained SNV ; \$(zcat ${vcf_file} | snpsift filter "(ANN[*].EFFECT has 'stop_gained')" | grep -v "^#" | wc -l)" >> ${report}
//       ## Stop lost
//       echo -e "Stop lost SNV ; \$(zcat ${vcf_file} | snpsift filter "(ANN[*].EFFECT has 'stop_lost')" | grep -v "^#" | wc -l)" >> ${report}
//       ## Splicing
//       echo -e "Splicing SNV ; \$(zcat ${vcf_file} | snpsift filter "(ANN[*].EFFECT has 'splice_site_region')" | grep -v "^#" | wc -l)" >> ${report}
//       ## LOF
//       echo -e "LOF SNV ; \$(zcat ${vcf_file} | snpsift filter "(LOF[*].PERC >= 0.5)" | grep -v "^#" | wc -l)" >> ${report}
//
//       ## Generate VCF of protein conding variants
//       zcat ${vcf_file} | \
//       snpsift filter " \
//           (ANN[*].EFFECT has 'synonymous_variant') || \
//           (ANN[*].EFFECT has 'missense_variant') || \
//           (ANN[*].EFFECT has 'stop_gained') || \
//           (ANN[*].EFFECT has 'frameshift_variant') || \
//           (ANN[*].EFFECT has 'stop_lost') || \
//           (ANN[*].EFFECT has 'inframe_insertion') || \
//           (ANN[*].EFFECT has 'inframe_deletion') || \
//           (ANN[*].EFFECT has 'coding_sequence_variant') " | \
//       bgzip -c > ${base}_coding.vcf.gz
//       """
//}
//
////
////'''
////All PGx SNVs
////'''
////annotate_dataset_snpeff_org.into{ annotate_dataset_snpeff_org; annotate_dataset_snpeff_org_variant_summary}
////process variant_summary {
////   tag "variant_summary_${dataset}"
////   label "small"
////   echo true
////   publishDir "${params.work_dir}/REPORTS/${dataset}/SUMMARY/ALL", mode:'copy', pattern: "*csv"
////   input:
////       set dataset, file(vcf_file) from annotate_dataset_snpeff_org_variant_summary
////   output:
////       set dataset, file(report) into variant_summary_report
////       set dataset, file("${vcf_file}"), file("${base}_coding.vcf.gz"), file("${base}.singletons") into variant_summary
////   script:
////       base = file(vcf_file.baseName).baseName
////       report = "${base}_all_variants_report.csv"
////       """
////       zcat ${vcf_file} | \
////       snpsift filter " \
////           (ANN[*].EFFECT has 'synonymous_variant') || \
////           (ANN[*].EFFECT has 'missense_variant') || \
////           (ANN[*].EFFECT has 'stop_gained') || \
////           (ANN[*].EFFECT has 'frameshift_variant') || \
////           (ANN[*].EFFECT has 'stop_lost') || \
////           (ANN[*].EFFECT has 'inframe_insertion') || \
////           (ANN[*].EFFECT has 'inframe_deletion') || \
////           (ANN[*].EFFECT has 'coding_sequence_variant') " | \
////       bgzip -c > ${base}_coding.vcf.gz
////        params.gene_regions_slop
////       ## Get total number positions of variants (Exome postions with SNV)
////       echo -e "Total positions ; \$(zcat ${vcf_file} | grep -v "^#" | wc -l)" > ${report}
////
////       ## Total rare variants
////       echo -e "Rare variants ; \$(vcftools --gzvcf ${vcf_file} --max-maf 0.01 --recode-INFO-all --recode --stdout | grep -v "^#" | wc -l)" >> ${report}
////
////       ## Total biallelic variants
////       echo -e "Biallelic positions ; \$(bcftools view -m2 -M2 ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
////
////       ## Total biallelic SNV variants
////       echo -e "Biallelic SNV positions ; \$(bcftools view -m2 -M2 -v snps ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
////
////       ## Total biallelic INDELS variants
////       echo -e "Biallelic INDELS positions ; \$(bcftools view -m2 -M2 -v indels ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
////
////       ## Total multiallelic variants
////       echo -e "Multiallelic positions ; \$(bcftools view -m3 ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
////
////       ## Total multiallelic SNV variants
////       echo -e "Multiallelic SNV positions ; \$(bcftools view -m3 -v snps ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
////
////       ## Total multiallelic INDELS variants
////       echo -e "Multiallelic INDELS positions ; \$(bcftools view -m3 -v indels ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
////
////       ## Total Singleton positions
////       vcftools \
////           --gzvcf ${vcf_file} \
////           --singletons \
////           --out ${base}
////       echo -e "Total Singletons ; \$(grep -v "^CHR" ${base}.singletons| wc -l)" >> ${report}
////
////       ## Nonsynonymous/missense
////       echo -e "Nonsynonymous SNV ; \$(zcat ${vcf_file} | snpsift filter "(ANN[*].EFFECT has 'missense_variant')" | grep -v "^#" | wc -l)" >> ${report}
////       ## Synonymous
////       echo -e "Synonymous SNV ; \$(zcat ${vcf_file} | snpsift filter "(ANN[*].EFFECT has 'synonymous_variant')" | grep -v "^#" | wc -l)" >> ${report}
////       ## Stop gained
////       echo -e "Stop gained SNV ; \$(zcat ${vcf_file} | snpsift filter "(ANN[*].EFFECT has 'stop_gained')" | grep -v "^#" | wc -l)" >> ${report}
////       ## Stop lost
////       echo -e "Stop lost SNV ; \$(zcat ${vcf_file} | snpsift filter "(ANN[*].EFFECT has 'stop_lost')" | grep -v "^#" | wc -l)" >> ${report}
////       ## Splicing
////       echo -e "Splicing SNV ; \$(zcat ${vcf_file} | snpsift filter "(ANN[*].EFFECT has 'splice_site_region')" | grep -v "^#" | wc -l)" >> ${report}
////       ## LOF
////       echo -e "LOF SNV ; \$(zcat ${vcf_file} | snpsift filter "(LOF[*].PERC >= 0.5)" | grep -v "^#" | wc -l)" >> ${report}
////       """
////}
//
//
//'''
//All coding SNVs
//'''
//variant_summary.into{ variant_summary; variant_summary_1}
//process variant_summary_coding {
//   tag "Summary_coding_${dataset}"
//   label "small"
//   publishDir "${params.work_dir}/REPORTS/${dataset}/SUMMARY/CODING", mode:'copy', pattern: "*csv"
//   input:
//       set dataset, file(vcf_file), file(vcf_file_coding), file(dataset_singletons) from variant_summary_1
//   output:
//       set dataset, file(report) into variant_summary_report_coding,variant_summary_report_coding_1
//   script:
//       base = file(vcf_file.baseName).baseName
//       report = "${base}_all_coding_variants_report.csv"
//       """
//       ## Get number of protein coding variants (Exome SNV)
//       echo -e "Region;All Coding" > ${report}
//       echo -e "Total Exome SNV ; \$(bcftools view --min-ac 1 ${base}_coding.vcf.gz | grep -v "^#" | wc -l)" >> ${report}
//       ## Total rare variants
//       echo -e "Rare variants ; \$(vcftools --gzvcf ${base}_coding.vcf.gz --max-maf 0.01 --recode-INFO-all --recode --stdout | grep -v "^#" | wc -l)" >> ${report}
//       echo -e "Exome Singletons ; \$( vcftools --gzvcf ${base}_coding.vcf.gz --singletons --stdout | grep -v "^CHR" | wc -l)" >> ${report}
//       ## Singletons
//       ## Nonsynonymous/missense
//       echo -e "Nonsynonymous SNV ; \$(zcat ${base}_coding.vcf.gz | snpsift filter "(ANN[*].EFFECT has 'missense_variant')" | grep -v "^#" | wc -l)" >> ${report}
//       ## Synonymous
//       echo -e "Synonymous SNV ; \$(zcat ${base}_coding.vcf.gz | snpsift filter "(ANN[*].EFFECT has 'synonymous_variant')" | grep -v "^#" | wc -l)" >> ${report}
//       ## Stop gained
//       echo -e "Stop gained SNV ; \$(zcat ${base}_coding.vcf.gz | snpsift filter "(ANN[*].EFFECT has 'stop_gained')" | grep -v "^#" | wc -l)" >> ${report}
//       ## Stop lost
//       echo -e "Stop lost SNV ; \$(zcat ${base}_coding.vcf.gz | snpsift filter "(ANN[*].EFFECT has 'stop_lost')" | grep -v "^#" | wc -l)" >> ${report}
//       ## Splicing
//       echo -e "Splicing SNV ; \$(zcat ${base}_coding.vcf.gz | snpsift filter "(ANN[*].EFFECT has 'splice_site_region')" | grep -v "^#" | wc -l)" >> ${report}
//       ## LOF
//       echo -e "LOF SNV ; \$(zcat ${base}_coding.vcf.gz | snpsift filter "(LOF[*].PERC >= 0.5)" | grep -v "^#" | wc -l)" >> ${report}
//       """
//}
//
//
//
//"""
//"""
//annotate_dataset_snpeff_org.into{ annotate_dataset_snpeff_org; annotate_dataset_snpeff_org_1 }
//process stats_dataset_snpeff_orig{
//   tag "stats_${file(vcf_file.baseName).baseName}"
//   label "small"
//   publishDir "${params.work_dir}/data/${dataset}/VCF/", mode:'copy'
//   input:
//       set val(dataset), file(vcf_file) from annotate_dataset_snpeff_org_1
//   output:
//       set dataset, file("${base}.samples") into stats_dataset_snpeff_org
//       set dataset, file("${base}.sites") into dataset_sites
//   script:
//       base = file(vcf_file.baseName).baseName
//       """
//       bcftools query -l ${vcf_file} > ${base}.samples
//       bcftools query -f'%CHROM:%POS\\t%POS\\t%REF\\t%ALT\\n' ${vcf_file} > ${base}.sites
//       """
//}
//
//samples = stats_dataset_snpeff_org.splitCsv(elem:1).collect { it[1][0] }
//
//
//'''
//All SNV files per ind
//'''
//variant_summary.into{ variant_summary; variant_summary_2}
//process variant_summary_ind_file {
//   tag "vcf_ind_${sample}"
//   label "small"
//   input:
//       set dataset, file(vcf_file), file(vcf_file_coding), file(dataset_singletons) from variant_summary_2
//       each sample from samples
//   output:
//       set dataset, sample, file("${base}_${sample}.vcf.gz"), file("${base}_${sample}_coding.vcf.gz"), file("${sample}.singletons") into variant_summary_ind_file
//   script:
//       base = file(vcf_file.baseName).baseName
//       """
//       ## Generate vcf for single sample and only heterozygote variants
//       bcftools view -s ${sample} ${vcf_file} | bcftools view -i '(GT="het" & GT!~"\\.") || (GT="AA" & GT!~"\\.")' | bgzip -c > ${base}_${sample}.vcf.gz
//       ## Generate VCF of protein conding variants
//       bcftools view -s ${sample} ${vcf_file_coding} | bcftools view -i '(GT="het" & GT!~"\\.") || (GT="AA" & GT!~"\\.")' | bgzip -c > ${base}_${sample}_coding.vcf.gz
//       awk '/${sample}/ {print \$1"\t"\$2}' ${dataset_singletons} > ${sample}.singletons
//       """
//}
//
//
//'''
//All SNVs summary per ind
//'''
//variant_summary_ind_file.into{ variant_summary_ind_file; variant_summary_ind_file_1}
//process variant_summary_ind {
//   tag "summary_ind_${sample}"
//   label "small"
//   input:
//       set dataset, sample, file(vcf_file), file(vcf_file_coding), file(singletons) from variant_summary_ind_file_1
//   output:
//       set dataset, sample, file(report) into variant_summary_ind
//   script:
//       base = file(vcf_file.baseName).baseName
//       report = "${base}_${sample}_all_variants_report.csv"
//       """
//       ## Get total number positions of variants (Exome postions with SNV)
//       echo -e "Total positions ; \$(zcat ${vcf_file} | grep -v "^#" | wc -l)" > ${report}
//       ## Total rare variants
//       echo -e "Rare variants ; \$(vcftools --gzvcf ${vcf_file} --max-maf 0.01 --recode-INFO-all --recode --stdout | grep -v "^#" | wc -l)" >> ${report}
//       ## Total biallelic variants
//       echo -e "Biallelic positions ; \$(bcftools view -m2 -M2 ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//       ## Total biallelic SNV variants
//       echo -e "Biallelic SNV positions ; \$(bcftools view -m2 -M2 -v snps ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//       ## Total biallelic INDELS variants
//       echo -e "Biallelic INDELS positions ; \$(bcftools view -m2 -M2 -v indels ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//       ## Total multiallelic variants
//       echo -e "Multiallelic positions ; \$(bcftools view -m3 ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//       ## Total multiallelic SNV variants
//       echo -e "Multiallelic SNV positions ; \$(bcftools view -m3 -v snps ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//       ## Total multiallelic INDELS variants
//       echo -e "Multiallelic INDELS positions ; \$(bcftools view -m3 -v indels ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//       ## Total Singleton positions
//       echo -e "Total Singletons ; \$(grep -v "^CHR" ${singletons} | wc -l)" >> ${report}
//       """
//}
//
//
//'''
//All coding SNVs summary per ind
//'''
//variant_summary_ind_file.into{ variant_summary_ind_file; variant_summary_ind_file_2 }
//process coding_variant_summary_ind {
//   tag "summary_coding_ind_${sample}"
//   label "small"
//   input:
//       set dataset, sample, file(vcf_file), file(vcf_file_coding), file(singletons) from variant_summary_ind_file_2
//   output:
//       set dataset, sample, file(report) into coding_variant_summary_ind
//   script:
//       base = file(vcf_file.baseName).baseName
//       report = "${base}_${sample}_conding_variants_report.csv"
//       """
//       bcftools index --tbi -f ${vcf_file_coding}
//       ## Get number of protein coding variants (Exome SNV)
//       echo -e "Total Exome SNV ; \$(bcftools view --min-ac 1 ${vcf_file_coding} | grep -v "^#" | wc -l)" >> ${report}
//       ## Total rare variants
//       echo -e "Rare variants ; \$(vcftools --gzvcf ${vcf_file_coding} --max-maf 0.01 --recode-INFO-all --recode --stdout | grep -v "^#" | wc -l)" >> ${report}
//       ## Singletons
//       echo -e "Exome Singletons ; \$(bcftools view -i 'AC=1' ${vcf_file_coding} --regions-file ${singletons} | grep -v "^#" | wc -l)" >> ${report}
//       ## Nonsynonymous/missense
//       echo -e "Nonsynonymous SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'missense_variant')" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//       ## Synonymous
//       echo -e "Synonymous SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'synonymous_variant')" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//       ## Stop gained
//       echo -e "Stop gained SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'stop_gained')" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//       ## Stop lost
//       echo -e "Stop lost SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'stop_lost')" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//       ## Splicing
//       echo -e "Splicing SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'ssplice_site_region')" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//       ## LOF
//       echo -e "LOF SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(LOF[*].PERC >= 0.5)" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//       """
//}
//
//"""
//Combined per ind coding SNVs
//"""
//coding_variant_summary_ind.into{ coding_variant_summary_ind; coding_variant_summary_ind_1}
//coding_variant_summary_ind_1_list = [:]
//coding_variant_summary_ind_1.toSortedList().val.each{ dataset, sample, report ->
//   if(!(dataset in coding_variant_summary_ind_1_list.keySet())){
//       coding_variant_summary_ind_1_list[dataset] = [dataset, sample+":"+report]
//   }
//   else{
//       coding_variant_summary_ind_1_list[dataset][1] += ',' + sample+":"+report
//   }
//
//}
//process variant_summary_ind_comb {
//   tag "variant_summary_${dataset}"
//   label "small"
//   publishDir "${params.work_dir}/REPORTS/${dataset}/SUMMARY/CODING/IND", mode: 'copy', pattern: "*csv"
//   input:
//       set dataset, reports_list from coding_variant_summary_ind_1_list.values()
//   output:
//       set dataset, file(combine_report) into variant_summary_ind_comb,variant_summary_ind_comb_1
//       set dataset, file(combine_report_summary) into variant_summary_ind_comb_summary,variant_summary_ind_comb_1_summary
//       set dataset, file(combine_report_accum) into variant_summary_ind_comb_accum
//   script:
//       title = "All Coding"
//       combine_report = "${dataset}_coding_variants_report_inds.csv"
//       combine_report_summary = "${dataset}_coding_variants_report_inds_summary.csv"
//       combine_report_accum = "${dataset}_coding_variants_report_inds_accum.csv"
//       template "combine_report_samples.py"
//}
//
//
////// PRIVATE VARIANTS
//
//"""
//Generate private variants list
//"""
//dataset_sites.into{ dataset_sites; dataset_sites_1 }
//dataset_sites_1_list = dataset_sites_1.toSortedList().val
//get_private_variants_data = []
//dataset_sites_1_list.each { dataset, sites ->
//   get_private_variants_data << [dataset, file(sites), mafs_annotations.values().join(';')]
//}
//
//process get_private_variants {
//   tag "private_${dataset}"
//   label "bigmem"
//   publishDir "${params.work_dir}/REPORTS/${dataset}/SUMMARY", mode: 'copy', pattern: "*csv"
//   input:
//       set dataset, file(mainTSV), otherTSVs from get_private_variants_data
//   output:
//       set dataset, file(private_list_file) into get_private_variants,get_private_variants_1,get_private_variants_2,get_private_variants_3
//   script:
//       private_list_file = "${dataset}_private_sites.csv"
//       template "intersect_variants.py"
//}
//
//'''
//All private SNVs
//'''
//variant_summary.into{ variant_summary; variant_summary_1}
//get_private_variants.into{get_private_variants; get_private_variants_1}
//variant_summary_private_data = []
//variant_summary_adme_data = []
//variant_summary_known_data = []
//variant_summary_1_list = variant_summary_1.toSortedList().val
//get_private_variants_list = get_private_variants.toSortedList().val
//variant_summary_1_list.each{ dataset, vcf_file, vcf_file_coding, singletons ->
//   get_private_variants_list.each{ dataset_, private_list_file ->
//       if(dataset_ == dataset){
//           // for private variants
//           variant_summary_private_data << [dataset, file(vcf_file), file(vcf_file_coding), file(singletons), file(private_list_file)]
//           // for known variants
//           variant_summary_known_data << [dataset, file(vcf_file), file(vcf_file_coding), file(singletons), file(private_list_file)]
//       }
//   }
//   // For adme variants
//   variant_summary_adme_data << [dataset, file(vcf_file), file(vcf_file_coding), file(singletons), file(params.gene_regions_slop)]
//}
//process variant_summary_private {
//   tag "summary_private_${dataset}"
//   label "bigmem"
//   publishDir "${params.work_dir}/REPORTS/${dataset}/SUMMARY/ALL/PRIVATE", mode:'copy', pattern: "*csv"
//   input:
//       set dataset, file(vcf_file), file(vcf_file_coding), file(singletons), file(private_list_file) from variant_summary_private_data
//   output:
//       set dataset, file(report) into variant_summary_private_report,variant_summary_private_report_1
//       set dataset, file("${base}_private.vcf.gz"), file("${base}_coding_private.vcf.gz") into variant_summary_private
//   script:
//       base = file(vcf_file.baseName).baseName
//       report = "${base}_all_private_variants_report.csv"
//       """
//       ## Subset vcf for private variants only
//       bcftools index --tbi ${vcf_file}
//       bcftools \
//           view \
//           --targets-file ${private_list_file} \
//           ${vcf_file} \
//           -Oz -o  ${base}_private.vcf.gz
//       bcftools index --tbi ${base}_private.vcf.gz
//       ## Get total number positions of variants (Exome postions with SNV)
//       echo -e "Region;All Private" > ${report}
//       echo -e "Total positions ; \$(zcat ${base}_private.vcf.gz | grep -v "^#" | wc -l)" >> ${report}
//       ## Total rare variants
//       echo -e "Rare variants ; \$(vcftools --gzvcf ${base}_private.vcf.gz --max-maf 0.01 --recode-INFO-all --recode --stdout | grep -v "^#" | wc -l)" >> ${report}
//       ## Total biallelic variants
//       echo -e "Biallelic positions ; \$(bcftools view -m2 -M2 ${base}_private.vcf.gz | grep -v "^#" | wc -l )" >> ${report}
//       ## Total biallelic SNV variants
//       echo -e "Biallelic SNV positions ; \$(bcftools view -m2 -M2 -v snps ${base}_private.vcf.gz | grep -v "^#" | wc -l )" >> ${report}
//       ## Total biallelic INDELS variants
//       echo -e "Biallelic INDELS positions ; \$(bcftools view -m2 -M2 -v indels ${base}_private.vcf.gz | grep -v "^#" | wc -l )" >> ${report}
//       ## Total multiallelic variants
//       echo -e "Multiallelic positions ; \$(bcftools view -m3 ${base}_private.vcf.gz | grep -v "^#" | wc -l )" >> ${report}
//       ## Total multiallelic SNV variants
//       echo -e "Multiallelic SNV positions ; \$(bcftools view -m3 -v snps ${base}_private.vcf.gz | grep -v "^#" | wc -l )" >> ${report}
//       ## Total multiallelic INDELS variants
//       echo -e "Multiallelic INDELS positions ; \$(bcftools view -m3 -v indels ${base}_private.vcf.gz | grep -v "^#" | wc -l )" >> ${report}
//       ## Total Singleton positions
//       echo -e "Total Singletons ; \$(bcftools view -i 'AC=1' ${base}_private.vcf.gz | grep -v "^#" | wc -l)" >> ${report}
//       ## Generate VCF of protein conding variants
//       bcftools index --tbi ${vcf_file_coding}
//       bcftools \
//           view \
//           --targets-file ${private_list_file} ${vcf_file_coding} \
//           -Oz -o ${base}_coding_private.vcf.gz
//       """
//}
//
//
//'''
//All private coding SNVs
//'''
//variant_summary_private.into{ variant_summary_private; variant_summary_private_1}
//process variant_summary_private_coding {
//   tag "Summary_private_coding_${dataset}"
//   label "bigmem"
//   publishDir "${params.work_dir}/REPORTS/${dataset}/SUMMARY/CODING/PRIVATE", mode:'copy', pattern: "*csv"
//   input:
//       set dataset, file(vcf_file), file(vcf_file_coding) from variant_summary_private_1
//   output:
//       set dataset, file(report) into variant_summary_report_private_coding,variant_summary_report_private_coding_1
//   script:
//       base = file(vcf_file.baseName).baseName
//       report = "${base}_all_coding_private_variants_report.csv"
//       """
//       bcftools index --tbi ${vcf_file_coding}
//       echo -e "Region;Private Coding" > ${report}
//       ## Get number of protein coding variants (Exome SNV)
//       echo -e "Total Exome SNV ; \$(bcftools view --min-ac 1 ${vcf_file_coding} | grep -v "^#" | wc -l)" >> ${report}
//       ## Total rare variants
//       echo -e "Rare variants ; \$(vcftools --gzvcf ${vcf_file_coding} --max-maf 0.01 --recode-INFO-all --recode --stdout | grep -v "^#" | wc -l)" >> ${report}
//       echo -e "Exome Singletons ; \$(bcftools view -i 'AC=1' ${vcf_file_coding} | grep -v "^#" | wc -l)" >> ${report}
//       ## Singletons
//       ## Nonsynonymous/missense
//       echo -e "Nonsynonymous SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'missense_variant')" | grep -v "^#" | wc -l)" >> ${report}
//       ## Synonymous
//       echo -e "Synonymous SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'synonymous_variant')" | grep -v "^#" | wc -l)" >> ${report}
//       ## Stop gained
//       echo -e "Stop gained SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'stop_gained')" | grep -v "^#" | wc -l)" >> ${report}
//       ## Stop lost
//       echo -e "Stop lost SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'stop_lost')" | grep -v "^#" | wc -l)" >> ${report}
//       ## Splicing
//       echo -e "Splicing SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'splice_site_region')" | grep -v "^#" | wc -l)" >> ${report}
//       ## LOF
//       echo -e "LOF SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(LOF[*].PERC >= 0.5)" | grep -v "^#" | wc -l)" >> ${report}
//       """
//}
//
//
//'''
//All private SNV files per ind
//'''
//variant_summary_ind_file.into{ variant_summary_ind_file; variant_summary_ind_file_3}
//variant_summary_ind_file_3_data = []
//variant_summary_ind_file_adme_data = []
//variant_summary_ind_file_3_list = variant_summary_ind_file_3.toSortedList().val
//variant_summary_ind_file_3_list.each{ dataset, sample, vcf_file, vcf_file_coding, singletons ->
//   get_private_variants_list.each{ dataset_, private_list_file ->
//       if(dataset_ == dataset){
//           variant_summary_ind_file_3_data << [dataset, sample, file(vcf_file), file(vcf_file_coding), file(singletons), file(private_list_file)]
//       }
//   }
//   //// For ADME gene variants per individual
//   variant_summary_ind_file_adme_data << [dataset, sample, file(vcf_file), file(vcf_file_coding), file(singletons), file(params.gene_regions_slop)]
//}
//process variant_private_summary_ind_file {
//   tag "vcf_ind_${sample}"
//   label "small"
//   input:
//       set dataset, sample, file(vcf_file), file(vcf_file_coding), file(singletons), file(private_list_file) from variant_summary_ind_file_3_data
//   output:
//       set dataset, sample, file("${base}_private.vcf.gz"), file("${base}_coding_private.vcf.gz"), file(singletons) into variant_private_summary_ind_file
//   script:
//       base = file(vcf_file.baseName).baseName
//       """
//       ## Generate vcf for single sample and only heterozygote variants
//       bcftools view \
//           --targets-file ${private_list_file} \
//           ${vcf_file} |\
//           bcftools view \
//               -i '(GT="het" & GT!~"\\.") || (GT="AA" & GT!~"\\.")' |\
//           bgzip -c > ${base}_private.vcf.gz
//       ## Generate VCF of protein conding variants
//       bcftools view \
//           --targets-file ${private_list_file} \
//           ${vcf_file_coding} |\
//           bcftools view -i '(GT="het" & GT!~"\\.") || (GT="AA" & GT!~"\\.")' |\
//           bgzip -c > ${base}_coding_private.vcf.gz
//
//       """
//}
//
//
//'''
//All private SNVs summary per ind
//'''
//variant_private_summary_ind_file.into{ variant_private_summary_ind_file; variant_private_summary_ind_file_1}
//process variant_private_summary_ind {
//   tag "summary_private_ind_${sample}"
//   label "small"
//   input:
//       set dataset, sample, file(vcf_file), file(vcf_file_coding), file(singletons) from variant_private_summary_ind_file_1
//   output:
//       set dataset, sample, file(report) into variant_private_summary_ind
//   script:
//       base = file(vcf_file.baseName).baseName
//       report = "${base}_${sample}_all_variants_private_report.csv"
//       """
//       ## Get total number positions of variants (Exome postions with SNV)
//       echo -e "Region;Private per Ind" > ${report}
//       echo -e "Total positions ; \$(zcat ${vcf_file} | grep -v "^#" | wc -l)" > ${report}
//       ## Total rare variants
//       echo -e "Rare variants ; \$(vcftools --gzvcf ${vcf_file} --max-maf 0.01 --recode-INFO-all --recode --stdout | grep -v "^#" | wc -l)" >> ${report}
//       ## Total biallelic variants
//       echo -e "Biallelic positions ; \$(bcftools view -m2 -M2 ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//       ## Total biallelic SNV variants
//       echo -e "Biallelic SNV positions ; \$(bcftools view -m2 -M2 -v snps ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//       ## Total biallelic INDELS variants
//       echo -e "Biallelic INDELS positions ; \$(bcftools view -m2 -M2 -v indels ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//       ## Total multiallelic variants
//       echo -e "Multiallelic positions ; \$(bcftools view -m3 ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//       ## Total multiallelic SNV variants
//       echo -e "Multiallelic SNV positions ; \$(bcftools view -m3 -v snps ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//       ## Total multiallelic INDELS variants
//       echo -e "Multiallelic INDELS positions ; \$(bcftools view -m3 -v indels ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//       ## Total Singleton positions
//       echo -e "Total Singletons ; \$(wc -l ${singletons})" >> ${report}
//       """
//}
//
//
//'''
//All coding private SNVs summary per ind
//'''
//variant_private_summary_ind_file.into{ variant_private_summary_ind_file; variant_private_summary_ind_file_2}
//process coding_private_variant_summary_ind {
//   tag "summary_coding_private_ind_${sample}"
//   label "small"
//   maxForks 10
//   input:
//       set dataset, sample, file(vcf_file), file(vcf_file_coding), file(singletons) from variant_private_summary_ind_file_2
//   output:
//       set dataset, sample, file(report) into coding_private_variant_summary_ind
//   script:
//       base = file(vcf_file.baseName).baseName
//       report = "${base}_${sample}_conding_private_variants_report.csv"
//       """
//       ## Get number of protein coding variants (Exome SNV)
//       echo -e "Region;Private per Ind" > ${report}
//       echo -e "Total Exome SNV ; \$(bcftools view --min-ac 1 ${vcf_file_coding} | grep -v "^#" | wc -l)" >> ${report}
//       ## Total rare variants
//       echo -e "Rare variants ; \$(vcftools --gzvcf ${vcf_file_coding} --max-maf 0.01 --recode-INFO-all --recode --stdout | grep -v "^#" | wc -l)" >> ${report}
//       ## Singletons
//       echo -e "Exome Singletons ; \$(bcftools view -i 'AC=1' ${vcf_file_coding} | grep -v "^#" | wc -l)" >> ${report}
//       ## Nonsynonymous/missense
//       echo -e "Nonsynonymous SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'missense_variant')" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//       ## Synonymous
//       echo -e "Synonymous SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'synonymous_variant')" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//       ## Stop gained
//       echo -e "Stop gained SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'stop_gained')" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//       ## Stop lost
//       echo -e "Stop lost SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'stop_lost')" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//       ## Splicing
//       echo -e "Splicing SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'ssplice_site_region')" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//       ## LOF
//       echo -e "LOF SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(LOF[*].PERC >= 0.5)" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//       """
//}
//
//"""
//Combined per ind coding private SNVs
//"""
//coding_private_variant_summary_ind.into{ coding_private_variant_summary_ind; coding_private_variant_summary_ind_1}
//coding_private_variant_summary_ind_1_list = [:]
//coding_private_variant_summary_ind_1.toSortedList().val.each{ dataset, sample, report ->
//    if(!(dataset in coding_private_variant_summary_ind_1_list.keySet())){
//        coding_private_variant_summary_ind_1_list[dataset] = [dataset, sample+":"+report]
//    }
//    else{
//        coding_private_variant_summary_ind_1_list[dataset][1] += ',' + sample+":"+report
//    }
//
//}
//process variant_private_summary_ind_comb {
//    tag "variant_private_summary_${dataset}"
//    label "small"
//    publishDir "${params.work_dir}/REPORTS/${dataset}/SUMMARY/CODING/IND", mode: 'copy', pattern: "*csv"
//    input:
//        set dataset, reports_list from coding_private_variant_summary_ind_1_list.values()
//    output:
//        set dataset, file(combine_report) into variant_private_summary_ind_comb,variant_private_summary_ind_comb_1
//        set dataset, file(combine_report_summary) into variant_private_summary_ind_comb_summary,variant_private_summary_ind_comb_1_summary
//        set dataset, file(combine_report_accum) into variant_private_summary_ind_comb_accum
//    script:
//        title = "Private per Ind"
//        combine_report = "${dataset}_coding_variants_private_report_inds.csv"
//        combine_report_summary = "${dataset}_coding_variants_private_report_inds_summary.csv"
//        combine_report_accum = "${dataset}_coding_variants_private_report_inds_accum.csv"
//        template "combine_report_samples.py"
//}
//
//
/////// KNOWN VARIANTS
//
//process variant_summary_known {
//    tag "summary_known_${dataset}"
//    label "bigmem"
//    publishDir "${params.work_dir}/REPORTS/${dataset}/SUMMARY/ALL/KNOWN", mode:'copy', pattern: "*csv"
//    input:
//        set dataset, file(vcf_file), file(vcf_file_coding), file(singletons), file(private_list_file) from variant_summary_known_data
//    output:
//        set dataset, file(report) into variant_summary_known_report,variant_summary_known_report_1
//        set dataset, file("${base}_known.vcf.gz"), file("${base}_coding_known.vcf.gz") into variant_summary_known
//    script:
//        base = file(vcf_file.baseName).baseName
//        report = "${base}_all_known_variants_report.csv"
//        """
//        ## Subset vcf for known variants only
//        bcftools index --tbi ${vcf_file}
//        bcftools \
//            view \
//            --targets-file ^${private_list_file} \
//            ${vcf_file} \
//            -Oz -o  ${base}_known.vcf.gz
//        bcftools index --tbi ${base}_known.vcf.gz
//        ## Get total number positions of variants (Exome postions with SNV)
//        echo -e "Region;All known" > ${report}
//        echo -e "Total positions ; \$(zcat ${base}_known.vcf.gz | grep -v "^#" | wc -l)" >> ${report}
//        ## Total rare variants
//        echo -e "Rare variants ; \$(vcftools --gzvcf ${base}_known.vcf.gz --max-maf 0.01 --recode-INFO-all --recode --stdout | grep -v "^#" | wc -l)" >> ${report}
//        ## Total biallelic variants
//        echo -e "Biallelic positions ; \$(bcftools view -m2 -M2 ${base}_known.vcf.gz | grep -v "^#" | wc -l )" >> ${report}
//        ## Total biallelic SNV variants
//        echo -e "Biallelic SNV positions ; \$(bcftools view -m2 -M2 -v snps ${base}_known.vcf.gz | grep -v "^#" | wc -l )" >> ${report}
//        ## Total biallelic INDELS variants
//        echo -e "Biallelic INDELS positions ; \$(bcftools view -m2 -M2 -v indels ${base}_known.vcf.gz | grep -v "^#" | wc -l )" >> ${report}
//        ## Total multiallelic variants
//        echo -e "Multiallelic positions ; \$(bcftools view -m3 ${base}_known.vcf.gz | grep -v "^#" | wc -l )" >> ${report}
//        ## Total multiallelic SNV variants
//        echo -e "Multiallelic SNV positions ; \$(bcftools view -m3 -v snps ${base}_known.vcf.gz | grep -v "^#" | wc -l )" >> ${report}
//        ## Total multiallelic INDELS variants
//        echo -e "Multiallelic INDELS positions ; \$(bcftools view -m3 -v indels ${base}_known.vcf.gz | grep -v "^#" | wc -l )" >> ${report}
//        ## Total Singleton positions
//        echo -e "Total Singletons ; \$(bcftools view -i 'AC=1' ${base}_known.vcf.gz | grep -v "^#" | wc -l)" >> ${report}
//        ## Generate VCF of protein conding variants
//        bcftools index --tbi ${vcf_file_coding}
//        bcftools \
//            view \
//            --targets-file ^${private_list_file} ${vcf_file_coding} \
//            -Oz -o ${base}_coding_known.vcf.gz
//        """
//}
//
//
//'''
//All known coding SNVs
//'''
//variant_summary_known.into{ variant_summary_known; variant_summary_known_1}
//process variant_summary_known_coding {
//    tag "Summary_known_coding_${dataset}"
//    label "bigmem"
//    publishDir "${params.work_dir}/REPORTS/${dataset}/SUMMARY/CODING/KNOWN", mode:'copy', pattern: "*csv"
//    input:
//        set dataset, file(vcf_file), file(vcf_file_coding) from variant_summary_known_1
//    output:
//        set dataset, file(report) into variant_summary_report_known_coding,variant_summary_report_known_coding_1
//    script:
//        base = file(vcf_file.baseName).baseName
//        report = "${base}_all_coding_known_variants_report.csv"
//        """
//        bcftools index --tbi ${vcf_file_coding}
//        echo -e "Region;known Coding" > ${report}
//        ## Get number of protein coding variants (Exome SNV)
//        echo -e "Total Exome SNV ; \$(bcftools view --min-ac 1 ${vcf_file_coding} | grep -v "^#" | wc -l)" >> ${report}
//        echo -e "Exome Singletons ; \$(bcftools view -i 'AC=1' ${vcf_file_coding} | grep -v "^#" | wc -l)" >> ${report}
//        ## Singletons
//        ## Nonsynonymous/missense
//        echo -e "Nonsynonymous SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'missense_variant')" | grep -v "^#" | wc -l)" >> ${report}
//        ## Synonymous
//        echo -e "Synonymous SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'synonymous_variant')" | grep -v "^#" | wc -l)" >> ${report}
//        ## Stop gained
//        echo -e "Stop gained SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'stop_gained')" | grep -v "^#" | wc -l)" >> ${report}
//        ## Stop lost
//        echo -e "Stop lost SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'stop_lost')" | grep -v "^#" | wc -l)" >> ${report}
//        ## Splicing
//        echo -e "Splicing SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'splice_site_region')" | grep -v "^#" | wc -l)" >> ${report}
//        ## LOF
//        echo -e "LOF SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(LOF[*].PERC >= 0.5)" | grep -v "^#" | wc -l)" >> ${report}
//        """
//}
//
//
//'''
//All known SNV files per ind
//'''
//process variant_known_summary_ind_file {
//    tag "vcf_ind_${sample}"
//    label "small"
//    input:
//        set dataset, sample, file(vcf_file), file(vcf_file_coding), file(singletons), file(private_list_file) from variant_summary_ind_file_3_data
//    output:
//        set dataset, sample, file("${base}_known.vcf.gz"), file("${base}_coding_known.vcf.gz"), file(singletons) into variant_known_summary_ind_file
//    script:
//        base = file(vcf_file.baseName).baseName
//        """
//        ## Generate vcf for single sample and only heterozygote variants
//        bcftools view \
//            --targets-file ^${private_list_file} \
//            ${vcf_file} |\
//            bcftools view \
//                -i '(GT="het" & GT!~"\\.") || (GT="AA" & GT!~"\\.")' |\
//            bgzip -c > ${base}_known.vcf.gz
//        ## Generate VCF of protein conding variants
//        bcftools view \
//            --targets-file ^${private_list_file} \
//            ${vcf_file_coding} |\
//            bcftools view -i '(GT="het" & GT!~"\\.") || (GT="AA" & GT!~"\\.")' |\
//            bgzip -c > ${base}_coding_known.vcf.gz
//
//        """
//}
//
//
//'''
//All known SNVs summary per ind
//'''
//variant_known_summary_ind_file.into{ variant_known_summary_ind_file; variant_known_summary_ind_file_1}
//process variant_known_summary_ind {
//    tag "summary_known_ind_${sample}"
//    label "small"
//    input:
//        set dataset, sample, file(vcf_file), file(vcf_file_coding), file(singletons) from variant_known_summary_ind_file_1
//    output:
//        set dataset, sample, file(report) into variant_known_summary_ind
//    script:
//        base = file(vcf_file.baseName).baseName
//        report = "${base}_${sample}_all_variants_known_report.csv"
//        """
//        ## Get total number positions of variants (Exome postions with SNV)
//        echo -e "Region;known per Ind" > ${report}
//        echo -e "Total positions ; \$(zcat ${vcf_file} | grep -v "^#" | wc -l)" > ${report}
//        ## Total rare variants
//        echo -e "Rare variants ; \$(vcftools --gzvcf ${vcf_file} --max-maf 0.01 --recode-INFO-all --recode --stdout | grep -v "^#" | wc -l)" >> ${report}
//        ## Total biallelic variants
//        echo -e "Biallelic positions ; \$(bcftools view -m2 -M2 ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//        ## Total biallelic SNV variants
//        echo -e "Biallelic SNV positions ; \$(bcftools view -m2 -M2 -v snps ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//        ## Total biallelic INDELS variants
//        echo -e "Biallelic INDELS positions ; \$(bcftools view -m2 -M2 -v indels ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//        ## Total multiallelic variants
//        echo -e "Multiallelic positions ; \$(bcftools view -m3 ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//        ## Total multiallelic SNV variants
//        echo -e "Multiallelic SNV positions ; \$(bcftools view -m3 -v snps ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//        ## Total multiallelic INDELS variants
//        echo -e "Multiallelic INDELS positions ; \$(bcftools view -m3 -v indels ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//        ## Total Singleton positions
//        echo -e "Total Singletons ; \$(wc -l ${singletons})" >> ${report}
//        """
//}
//
//
//'''
//All coding known SNVs summary per ind
//'''
//variant_known_summary_ind_file.into{ variant_known_summary_ind_file; variant_known_summary_ind_file_2}
//process coding_known_variant_summary_ind {
//    tag "summary_coding_known_ind_${sample}"
//    label "small"
//    maxForks 10
//    input:
//        set dataset, sample, file(vcf_file), file(vcf_file_coding), file(singletons) from variant_known_summary_ind_file_2
//    output:
//        set dataset, sample, file(report) into coding_known_variant_summary_ind
//    script:
//        base = file(vcf_file.baseName).baseName
//        report = "${base}_${sample}_conding_known_variants_report.csv"
//        """
//        ## Get number of protein coding variants (Exome SNV)
//        echo -e "Region;known per Ind" > ${report}
//        echo -e "Total Exome SNV ; \$(bcftools view --min-ac 1 ${vcf_file_coding} | grep -v "^#" | wc -l)" >> ${report}
//        ## Total rare variants
//        echo -e "Rare variants ; \$(vcftools --gzvcf ${vcf_file_coding} --max-maf 0.01 --recode-INFO-all --recode --stdout | grep -v "^#" | wc -l)" >> ${report}
//        ## Singletons
//        echo -e "Exome Singletons ; \$(bcftools view -i 'AC=1' ${vcf_file_coding} | grep -v "^#" | wc -l)" >> ${report}
//        ## Nonsynonymous/missense
//        echo -e "Nonsynonymous SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'missense_variant')" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//        ## Synonymous
//        echo -e "Synonymous SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'synonymous_variant')" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//        ## Stop gained
//        echo -e "Stop gained SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'stop_gained')" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//        ## Stop lost
//        echo -e "Stop lost SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'stop_lost')" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//        ## Splicing
//        echo -e "Splicing SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'ssplice_site_region')" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//        ## LOF
//        echo -e "LOF SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(LOF[*].PERC >= 0.5)" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//        """
//}
//
//"""
//Combined per ind coding known SNVs
//"""
//coding_known_variant_summary_ind.into{ coding_known_variant_summary_ind; coding_known_variant_summary_ind_1}
//coding_known_variant_summary_ind_1_list = [:]
//coding_known_variant_summary_ind_1.toSortedList().val.each{ dataset, sample, report ->
//    if(!(dataset in coding_known_variant_summary_ind_1_list.keySet())){
//        coding_known_variant_summary_ind_1_list[dataset] = [dataset, sample+":"+report]
//    }
//    else{
//        coding_known_variant_summary_ind_1_list[dataset][1] += ',' + sample+":"+report
//    }
//
//}
//process variant_known_summary_ind_comb {
//    tag "variant_known_summary_${dataset}"
//    label "small"
//    publishDir "${params.work_dir}/REPORTS/${dataset}/SUMMARY/CODING/KNOWN", mode: 'copy', pattern: "*csv"
//    input:
//        set dataset, reports_list from coding_known_variant_summary_ind_1_list.values()
//    output:
//        set dataset, file(combine_report) into variant_known_summary_ind_comb,variant_known_summary_ind_comb_1
//        set dataset, file(combine_report_summary) into variant_known_summary_ind_comb_summary,variant_known_summary_ind_comb_1_summary
//        set dataset, file(combine_report_accum) into variant_known_summary_ind_comb_accum
//    script:
//        title = "known per Ind"
//        combine_report = "${dataset}_coding_variants_known_report_inds.csv"
//        combine_report_summary = "${dataset}_coding_variants_known_report_inds_summary.csv"
//        combine_report_accum = "${dataset}_coding_variants_known_report_inds_accum.csv"
//        template "combine_report_samples.py"
//}
//
//
//
/////// ADME GENE VARIANTS
//
//process variant_summary_adme {
//    tag "summary_adme_${dataset}"
//    label "bigmem"
//    publishDir "${params.work_dir}/REPORTS/${dataset}/SUMMARY/ALL/ADME", mode:'copy', pattern: "*csv"
//    input:
//        set dataset, file(vcf_file), file(vcf_file_coding), file(singletons), file(adme_list_file) from variant_summary_adme_data
//    output:
//        set dataset, file(report) into variant_summary_adme_report,variant_summary_adme_report_1
//        set dataset, file("${base}_adme.vcf.gz"), file("${base}_coding_adme.vcf.gz") into variant_summary_adme,variant_summary_adme_2
//    script:
//        base = file(vcf_file.baseName).baseName
//        report = "${base}_all_adme_variants_report.csv"
//        """
//        ## Subset vcf for adme variants only
//        bcftools index --tbi ${vcf_file}
//        bcftools \
//            view \
//            --targets-file ${adme_list_file} \
//            ${vcf_file} \
//            -Oz -o  ${base}_adme.vcf.gz
//        bcftools index --tbi ${base}_adme.vcf.gz
//        echo -e "Region;All ADME" > ${report}
//        ## Get total number positions of variants (Exome postions with SNV)
//        echo -e "Total positions ; \$(zcat ${base}_adme.vcf.gz | grep -v "^#" | wc -l)" > ${report}
//        ## Total rare variants
//        echo -e "Rare variants ; \$(vcftools --gzvcf ${base}_adme.vcf.gz --max-maf 0.01 --recode-INFO-all --recode --stdout | grep -v "^#" | wc -l)" >> ${report}
//        ## Total biallelic variants
//        echo -e "Biallelic positions ; \$(bcftools view -m2 -M2 ${base}_adme.vcf.gz | grep -v "^#" | wc -l )" >> ${report}
//        ## Total biallelic SNV variants
//        echo -e "Biallelic SNV positions ; \$(bcftools view -m2 -M2 -v snps ${base}_adme.vcf.gz | grep -v "^#" | wc -l )" >> ${report}
//        ## Total biallelic INDELS variants
//        echo -e "Biallelic INDELS positions ; \$(bcftools view -m2 -M2 -v indels ${base}_adme.vcf.gz | grep -v "^#" | wc -l )" >> ${report}
//        ## Total multiallelic variants
//        echo -e "Multiallelic positions ; \$(bcftools view -m3 ${base}_adme.vcf.gz | grep -v "^#" | wc -l )" >> ${report}
//        ## Total multiallelic SNV variants
//        echo -e "Multiallelic SNV positions ; \$(bcftools view -m3 -v snps ${base}_adme.vcf.gz | grep -v "^#" | wc -l )" >> ${report}
//        ## Total multiallelic INDELS variants
//        echo -e "Multiallelic INDELS positions ; \$(bcftools view -m3 -v indels ${base}_adme.vcf.gz | grep -v "^#" | wc -l )" >> ${report}
//        ## Total Singleton positions
//        echo -e "Total Singletons ; \$(bcftools view -i 'AC=1' ${base}_adme.vcf.gz | grep -v "^#" | wc -l)" >> ${report}
//        ## Generate VCF of protein conding variants
//        bcftools index --tbi ${vcf_file_coding}
//        bcftools \
//            view \
//            --targets-file ${adme_list_file} ${vcf_file_coding} \
//            -Oz -o ${base}_coding_adme.vcf.gz
//        """
//}
//
//
//
//'''
//All adme coding SNVs
//'''
//variant_summary_adme.into{ variant_summary_adme; variant_summary_adme_1}
//process variant_summary_adme_coding {
//    tag "Summary_adme_coding_${dataset}"
//    label "bigmem"
//    publishDir "${params.work_dir}/REPORTS/${dataset}/SUMMARY/CODING/ADME", mode:'copy', pattern: "*csv"
//    input:
//        set dataset, file(vcf_file), file(vcf_file_coding) from variant_summary_adme_1
//    output:
//        set dataset, file(report) into variant_summary_report_adme_coding,variant_summary_report_adme_coding_1
//    script:
//        base = file(vcf_file.baseName).baseName
//        report = "${base}_all_coding_adme_variants_report.csv"
//        """
//        bcftools index --tbi ${vcf_file_coding}
//        echo -e "Region;Coding ADME" > ${report}
//        ## Get number of protein coding variants (Exome SNV)
//        echo -e "Total Exome SNV ; \$(bcftools view --min-ac 1 ${vcf_file_coding} | grep -v "^#" | wc -l)" >> ${report}
//        ## Total rare variants
//        echo -e "Rare variants ; \$(vcftools --gzvcf ${vcf_file_coding} --max-maf 0.01 --recode-INFO-all --recode --stdout | grep -v "^#" | wc -l)" >> ${report}
//        ## Singletons
//        echo -e "Exome Singletons ; \$(bcftools view -i 'AC=1' ${vcf_file_coding} | grep -v "^#" | wc -l)" >> ${report}
//        ## Nonsynonymous/missense
//        echo -e "Nonsynonymous SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'missense_variant')" | grep -v "^#" | wc -l)" >> ${report}
//        ## Synonymous
//        echo -e "Synonymous SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'synonymous_variant')" | grep -v "^#" | wc -l)" >> ${report}
//        ## Stop gained
//        echo -e "Stop gained SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'stop_gained')" | grep -v "^#" | wc -l)" >> ${report}
//        ## Stop lost
//        echo -e "Stop lost SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'stop_lost')" | grep -v "^#" | wc -l)" >> ${report}
//        ## Splicing
//        echo -e "Splicing SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'splice_site_region')" | grep -v "^#" | wc -l)" >> ${report}
//        ## LOF
//        echo -e "LOF SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(LOF[*].PERC >= 0.5)" | grep -v "^#" | wc -l)" >> ${report}
//        """
//}
//
//
//'''
//All coding adme private SNVs
//'''
//process variant_summary_coding_adme_private {
//    tag "coding_adme_private_${dataset}"
//    label "bigmem"
//    publishDir "${params.work_dir}/REPORTS/${dataset}/SUMMARY/CODING/ADME", mode:'copy', pattern: "*csv"
//    input:
//        set dataset, file(vcf_file), file(vcf_file_adme), dataset1, file(private_list_file) from variant_summary_adme_2.combine(get_private_variants_1)
//    output:
//        set dataset, file(report) into variant_summary_report_coding_adme_private
//    script:
//        base = file(vcf_file_adme.baseName).baseName
//        ap_vcf_file = "${base}_private.vcf.gz"
//        report = "${base}_all_coding_adme_private_variants_report.csv"
//        """
//        ## Subset vcf for private variants only from coding adme
//        bcftools index --tbi ${vcf_file_adme}
//        bcftools \
//           view \
//           --targets-file ${private_list_file} \
//           ${vcf_file_adme} \
//           -Oz -o  ${base}_private.vcf.gz
//        bcftools index --tbi ${base}_private.vcf.gz
//        echo -e "Region;Coding ADME Private" > ${report}
//        ## Get number of protein coding variants (Exome SNV)
//        echo -e "Total Exome SNV ; \$(bcftools view --min-ac 1 ${ap_vcf_file} | grep -v "^#" | wc -l)" >> ${report}
//        ## Total rare variants
//        echo -e "Rare variants ; \$(vcftools --gzvcf ${ap_vcf_file} --max-maf 0.01 --recode-INFO-all --recode --stdout | grep -v "^#" | wc -l)" >> ${report}
//        ## Singletons
//        echo -e "Exome Singletons ; \$(bcftools view -i 'AC=1' ${ap_vcf_file} | grep -v "^#" | wc -l)" >> ${report}
//        ## Nonsynonymous/missense
//        echo -e "Nonsynonymous SNV ; \$(zcat ${ap_vcf_file} | snpsift filter "(ANN[*].EFFECT has 'missense_variant')" | grep -v "^#" | wc -l)" >> ${report}
//        ## Synonymous
//        echo -e "Synonymous SNV ; \$(zcat ${ap_vcf_file} | snpsift filter "(ANN[*].EFFECT has 'synonymous_variant')" | grep -v "^#" | wc -l)" >> ${report}
//        ## Stop gained
//        echo -e "Stop gained SNV ; \$(zcat ${ap_vcf_file} | snpsift filter "(ANN[*].EFFECT has 'stop_gained')" | grep -v "^#" | wc -l)" >> ${report}
//        ## Stop lost
//        echo -e "Stop lost SNV ; \$(zcat ${ap_vcf_file} | snpsift filter "(ANN[*].EFFECT has 'stop_lost')" | grep -v "^#" | wc -l)" >> ${report}
//        ## Splicing
//        echo -e "Splicing SNV ; \$(zcat ${ap_vcf_file} | snpsift filter "(ANN[*].EFFECT has 'splice_site_region')" | grep -v "^#" | wc -l)" >> ${report}
//        ## LOF
//        echo -e "LOF SNV ; \$(zcat ${ap_vcf_file} | snpsift filter "(LOF[*].PERC >= 0.5)" | grep -v "^#" | wc -l)" >> ${report}
//        """
//}
//
//
//'''
//All adme SNV files per ind
//'''
//process variant_adme_summary_ind_file {
//    tag "vcf_ind_${sample}"
//    label "small"
//    input:
//        set dataset, sample, file(vcf_file), file(vcf_file_coding), file(singletons), file(adme_list_file) from variant_summary_ind_file_adme_data
//    output:
//        set dataset, sample, file("${base}_adme.vcf.gz"), file("${base}_coding_adme.vcf.gz"), file(singletons) into variant_adme_summary_ind_file
//    script:
//        base = file(vcf_file.baseName).baseName
//        """
//        ## Generate vcf for single sample and only heterozygote variants
//        bcftools index --tbi ${vcf_file}
//        bcftools view \
//            --targets-file ${adme_list_file} \
//            ${vcf_file} |\
//            bcftools view \
//                -i '(GT="het" & GT!~"\\.") || (GT="AA" & GT!~"\\.")' |\
//            bgzip -c > ${base}_adme.vcf.gz
//        ## Generate VCF of protein conding variants
//        bcftools view \
//            --targets-file ${adme_list_file} \
//            ${vcf_file_coding} |\
//            bcftools view -i '(GT="het" & GT!~"\\.") || (GT="AA" & GT!~"\\.")' |\
//            bgzip -c > ${base}_coding_adme.vcf.gz
//
//        """
//}
//
//
//
//
//'''
//All adme SNVs summary per ind
//'''
//variant_adme_summary_ind_file.into{ variant_adme_summary_ind_file; variant_adme_summary_ind_file_1}
//process variant_adme_summary_ind {
//    tag "summary_adme_ind_${sample}"
//    label "small"
//    input:
//        set dataset, sample, file(vcf_file), file(vcf_file_coding), file(singletons) from variant_adme_summary_ind_file_1
//    output:
//        set dataset, sample, file(report) into variant_adme_summary_ind
//    script:
//        base = file(vcf_file.baseName).baseName
//        report = "${base}_${sample}_all_variants_adme_report.csv"
//        """
//        ## Get total number positions of variants (Exome postions with SNV)
//        echo -e "Region;ADME per Ind" > ${report}
//        echo -e "Total positions ; \$(zcat ${vcf_file} | grep -v "^#" | wc -l)" >> ${report}
//        ## Total rare variants
//        echo -e "Rare variants ; \$(vcftools --gzvcf ${vcf_file} --max-maf 0.01 --recode-INFO-all --recode --stdout | grep -v "^#" | wc -l)" >> ${report}
//        ## Total biallelic variants
//        echo -e "Biallelic positions ; \$(bcftools view -m2 -M2 ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//        ## Total biallelic SNV variants
//        echo -e "Biallelic SNV positions ; \$(bcftools view -m2 -M2 -v snps ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//        ## Total biallelic INDELS variants
//        echo -e "Biallelic INDELS positions ; \$(bcftools view -m2 -M2 -v indels ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//        ## Total multiallelic variants
//        echo -e "Multiallelic positions ; \$(bcftools view -m3 ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//        ## Total multiallelic SNV variants
//        echo -e "Multiallelic SNV positions ; \$(bcftools view -m3 -v snps ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//        ## Total multiallelic INDELS variants
//        echo -e "Multiallelic INDELS positions ; \$(bcftools view -m3 -v indels ${vcf_file} | grep -v "^#" | wc -l )" >> ${report}
//        ## Total Singleton positions
//        echo -e "Total Singletons ; \$(wc -l ${singletons})" >> ${report}
//        """
//}
//
//
//
//'''
//All coding adme SNVs summary per ind
//'''
//variant_adme_summary_ind_file.into{ variant_adme_summary_ind_file; variant_adme_summary_ind_file_2; variant_adme_summary_ind_file_3}
//process coding_adme_variant_summary_ind {
//    tag "summary_coding_adme_ind_${sample}"
//    label "small"
//    maxForks 10
//    input:
//        set dataset, sample, file(vcf_file), file(vcf_file_coding), file(singletons) from variant_adme_summary_ind_file_2
//    output:
//        set dataset, sample, file(report) into coding_adme_variant_summary_ind
//    script:
//        base = file(vcf_file.baseName).baseName
//        report = "${base}_${sample}_coding_adme_variants_report.csv"
//        """
//        ## Get number of protein coding variants (Exome SNV)
//        echo -e "Region;ADME per Ind" > ${report}
//        echo -e "Total Exome SNV ; \$(bcftools view --min-ac 1 ${vcf_file_coding} | grep -v "^#" | wc -l)" >> ${report}
//        ## Total rare variants
//        echo -e "Rare variants ; \$(vcftools --gzvcf ${vcf_file_coding} --max-maf 0.01 --recode-INFO-all --recode --stdout | grep -v "^#" | wc -l)" >> ${report}
//        ## Singletons
//        echo -e "Exome Singletons ; \$(bcftools view -i 'AC=1' ${vcf_file_coding} | grep -v "^#" | wc -l)" >> ${report}
//        ## Nonsynonymous/missense
//        echo -e "Nonsynonymous SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'missense_variant')" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//        ## Synonymous
//        echo -e "Synonymous SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'synonymous_variant')" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//        ## Stop gained
//        echo -e "Stop gained SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'stop_gained')" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//        ## Stop lost
//        echo -e "Stop lost SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'stop_lost')" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//        ## Splicing
//        echo -e "Splicing SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'ssplice_site_region')" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//        ## LOF
//        echo -e "LOF SNV ; \$(zcat ${vcf_file_coding} | snpsift filter "(LOF[*].PERC >= 0.5)" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//        """
//}
//
//
//"""
//Combined per ind coding adme SNVs
//"""
//coding_adme_variant_summary_ind.into{ coding_adme_variant_summary_ind; coding_adme_variant_summary_ind_1}
//coding_adme_variant_summary_ind_1_list = [:]
//coding_adme_variant_summary_ind_1.toSortedList().val.each{ dataset, sample, report ->
//    if(!(dataset in coding_adme_variant_summary_ind_1_list.keySet())){
//        coding_adme_variant_summary_ind_1_list[dataset] = [dataset, sample+":"+report]
//    }
//    else{
//        coding_adme_variant_summary_ind_1_list[dataset][1] += ',' + sample+":"+report
//    }
//
//}
//process variant_adme_summary_ind_comb {
//    tag "variant_adme_summary_${dataset}"
//    label "small"
//    publishDir "${params.work_dir}/REPORTS/${dataset}/SUMMARY/CODING/ADME/IND", mode: 'copy', pattern: "*csv"
//    input:
//        set dataset, reports_list from coding_adme_variant_summary_ind_1_list.values()
//    output:
//        set dataset, file(combine_report) into variant_adme_summary_ind_comb,variant_adme_summary_ind_comb_1
//        set dataset, file(combine_report_summary) into variant_adme_summary_ind_comb_summary,variant_adme_summary_ind_comb_1_summary
//        set dataset, file(combine_report_accum) into variant_adme_summary_ind_comb_accum
//    script:
//        title = "Coding ADME per Ind"
//        combine_report = "${dataset}_coding_variants_adme_report_inds.csv"
//        combine_report_summary = "${dataset}_coding_variants_adme_report_inds_summary.csv"
//        combine_report_accum = "${dataset}_coding_variants_adme_report_inds_accum.csv"
//        template "combine_report_samples.py"
//}
//
//
//'''
//All coding adme private SNVs summary per ind
//'''
//process coding_adme_private_variant_summary_ind {
//    tag "summary_coding_adme_private_ind_${sample}"
//    label "small"
//    maxForks 10
//    input:
//        set dataset, sample, file(vcf_file), file(vcf_file_coding), file(singletons), dataset1, file(private_list_file) from variant_adme_summary_ind_file_3.combine(get_private_variants_3)
//    output:
//        set dataset, sample, file(report) into coding_adme_private_variant_summary_ind,coding_adme_private_variant_summary_ind_1
//    script:
//        base = file(vcf_file_coding.baseName).baseName
//        report = "${base}_${sample}_conding_adme_private_variants_report.csv"
//        ap_vcf_file_coding = "${base}_private.vcf.gz"
//        """
//        ## Subset vcf for private variants only from coding adme
//        bcftools index --tbi ${vcf_file_coding}
//        bcftools \
//           view \
//           --targets-file ${private_list_file} \
//           ${vcf_file_coding} \
//           -Oz -o  ${ap_vcf_file_coding}
//        ## Get number of protein coding variants (Exome SNV)
//        echo -e "Region;ADME per Ind" > ${report}
//        echo -e "Total Exome SNV ; \$(bcftools view --min-ac 1 ${ap_vcf_file_coding} | grep -v "^#" | wc -l)" >> ${report}
//        ## Total rare variants
//        echo -e "Rare variants ; \$(vcftools --gzvcf ${ap_vcf_file_coding} --max-maf 0.01 --recode-INFO-all --recode --stdout | grep -v "^#" | wc -l)" >> ${report}
//        ## Singletons
//        echo -e "Exome Singletons ; \$(bcftools view -i 'AC=1' ${ap_vcf_file_coding} | grep -v "^#" | wc -l)" >> ${report}
//        ## Nonsynonymous/missense
//        echo -e "Nonsynonymous SNV ; \$(zcat ${ap_vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'missense_variant')" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//        ## Synonymous
//        echo -e "Synonymous SNV ; \$(zcat ${ap_vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'synonymous_variant')" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//        ## Stop gained
//        echo -e "Stop gained SNV ; \$(zcat ${ap_vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'stop_gained')" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//        ## Stop lost
//        echo -e "Stop lost SNV ; \$(zcat ${ap_vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'stop_lost')" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//        ## Splicing
//        echo -e "Splicing SNV ; \$(zcat ${ap_vcf_file_coding} | snpsift filter "(ANN[*].EFFECT has 'ssplice_site_region')" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//        ## LOF
//        echo -e "LOF SNV ; \$(zcat ${ap_vcf_file_coding} | snpsift filter "(LOF[*].PERC >= 0.5)" | bcftools view --min-ac 1 | grep -v "^#" | wc -l)" >> ${report}
//        """
//}
//
//
//"""
//Combined per ind coding adme SNVs
//"""
//coding_adme_private_variant_summary_ind_1_list = [:]
//coding_adme_private_variant_summary_ind_1.toSortedList().val.each{ dataset, sample, report ->
//    if(!(dataset in coding_adme_private_variant_summary_ind_1_list.keySet())){
//        coding_adme_private_variant_summary_ind_1_list[dataset] = [dataset, sample+":"+report]
//    }
//    else{
//        coding_adme_private_variant_summary_ind_1_list[dataset][1] += ',' + sample+":"+report
//    }
//
//}
//process variant_adme_private_summary_ind_comb {
//    tag "variant_adme_private_summary_${dataset}"
//    label "small"
//    publishDir "${params.work_dir}/REPORTS/${dataset}/SUMMARY/CODING/ADME/IND", mode: 'copy', pattern: "*csv"
//    input:
//        set dataset, reports_list from coding_adme_private_variant_summary_ind_1_list.values()
//    output:
//        set dataset, file(combine_report) into variant_adme_private_summary_ind_comb,variant_adme_private_summary_ind_comb_1
//        set dataset, file(combine_report_summary) into variant_adme_private_summary_ind_comb_summary,variant_adme_private_summary_ind_comb_1_summary
//        set dataset, file(combine_report_accum) into variant_adme_private_summary_ind_comb_accum
//    script:
//        title = "Coding ADME Private per Ind"
//        combine_report = "${dataset}_coding_variants_adme_private_report_inds.csv"
//        combine_report_summary = "${dataset}_coding_variants_adme_private_report_inds_summary.csv"
//        combine_report_accum = "${dataset}_coding_variants_adme_private_report_inds_accum.csv"
//        template "combine_report_samples.py"
//}
//
//
//"""
//Final report
//"""
//variant_summary_report_coding_1
//        .join(variant_summary_ind_comb_1_summary, by: 0)
//        .join(variant_summary_report_private_coding_1, by: 0)
//        .join(variant_private_summary_ind_comb_1_summary, by: 0)
//        .join(variant_adme_summary_ind_comb_1_summary, by: 0)
//        .map { it -> [ it[0], it[1..-1].join(",") ] }
//        .set{summary}
//
//process final_summary {
//    tag "final_summary_${dataset}"
//    label "small"
//    publishDir "${params.work_dir}/REPORTS/${dataset}/SUMMARY", mode: 'copy', pattern: "*csv"
//    input:
//        set dataset, file(reports) from summary
//    output:
//        set dataset, file(combine_report) into final_summary
//    script:
//        combine_report = "${dataset}_variants_report.csv"
//        template "final_summary.py"
//}
//
//
//"""
//Step:
//"""
//
////"""
////Combine accummulative reports
////"""
//////variant_summary_ind_comb_accum
//////        .join(variant_summary_ind_comb_1_summary, by: 0)
//////        .join(variant_known_summary_ind_comb_accum, by: 0)
//////        .join(variant_adme_summary_ind_comb_accum, by: 0)
//////        .map { it -> [ it[0], it[1..-1].join(",") ] }
//////        .set{summary_accum}.view()
////
////
////
////////'''
////////Merge genes.pop.diff.snps.genotypes and population.pop.diff.snps.genotypes
////////'''
////////rareGlobCommon2_all.into{rareGlobCommon2_all; rareGlobCommon2_all_1}
////////process rareGlobCommon2_1 {
////////    echo true
////////    tag "rareGlobCommon2_1_${GROUP_POP}"
////////    cpus { 2 * task.attempt }
////////    memory { 2.GB * task.cpus }
////////    publishDir "${params.work_dir}/data/MERGED_POP/POPDIFF/${GROUP_POP}", mode:'symlink'
////////
////////    input:
////////        set val(GROUP_POP), file(GROUP_POP_diff_snps_vcf), val(GROUP_POP_diff_snps_genotypes), val(GROUP_POP_population_pop_diff_snps_genotypes), file(GROUP_POP_rsIndexPopDiff_txt), val(GROUP_POP_genes_pop_diff_snps_genotypes) from rareGlobCommon2_all_1
////////    output:
////////        set val(GROUP_POP), file("${GROUP_POP}.population.gene.pop.diff.snps.genotypes") into rareGlobCommon2_1
////////    script:
////////        pop_diff = GROUP_POP_population_pop_diff_snps_genotypes.text
////////        println file(pop_diff).readLines()
////////}
//////
//////
////////''''
////////Step 8.3: Compute allele frequencies for each population for rare/common sites
////////'''
////////rareGlobCommon2_all.into { rareGlobCommon2_all; rareGlobCommon2_pop }
////////rareGlobCommon2_pop_list = [] // For POP, GROUP_POP,  GROUP_POP_sample, GROUP_diff_snps_vcf
////////rareGlobCommon2_pop.toSortedList().val.each { GROUP_POP, GROUP_diff_snps_vcf, GROUP_POP_genotype, GROUP_POP_population_genotypes, GROUP_POP_rsIndexPopDiff, GROUP_POP_genes_genotypes, GROUP_POP_pop_only_genotypes ->
////////    data = file(group_pop_data[GROUP_POP][1]).readLines().collect { it.split()[1]}.unique()
////////    data.each { POP ->
////////        rareGlobCommon2_pop_list << [ POP, GROUP_POP, group_pop_data[GROUP_POP][1], GROUP_diff_snps_vcf ]
////////    }
////////}
////////rareGlobCommon2_pop_chan = Channel.from( rareGlobCommon2_pop_list )
////////
////////process rareGlobCommon3 {
////////    echo true
////////    tag "rareGlobCommon3_${POP}"
////////    cpus { 2 * task.attempt }
////////    memory { 2.GB * task.cpus }
////////    publishDir "${params.work_dir}/data/MERGED_POP/POPDIFF/${GROUP_POP}", mode:'symlink'
////////
////////    input:
////////        set val(POP), val(GROUP_POP), file(GROUP_POP_sample), file(GROUP_diff_snps_vcf) from rareGlobCommon2_pop_chan
////////    output:
////////        set val(POP), val(GROUP_POP), file("${POP}.diff.snps.ind.frq") into rareGlobCommon3_all
////////    script:
////////        """
////////        grep ${POP} ${GROUP_POP_sample} | cut -f1 > ${POP}.samples
////////        # Extract sites for population
////////        # Generate sites frequency for population
////////        vcftools --gzvcf ${GROUP_diff_snps_vcf} \
////////            --keep ${POP}.samples \
////////            --recode-INFO-all \
////////            --recode --stdout | \
////////        vcftools --vcf - \
////////            --freq2 \
////////            --out ${POP}.diff.snps.ind
////////        """
////////}
////////
////////
////////''''
////////Step 8.4: Plot allele frequencies for each population for rare/common sites
////////'''
////////rareGlobCommon2_all.into{rareGlobCommon2_all; rareGlobCommon2_4}
////////process plot_rareGlobCommon {
////////    tag "rareGlobCommon4_${GROUP_POP}"
////////    cpus { 2 * task.attempt }
////////    memory { 2.GB * task.cpus }
////////    publishDir "${params.work_dir}/data/MERGED_POP/POPDIFF/${GROUP_POP}", mode:'symlink'
////////
////////    input:
////////        set val(GROUP_POP), file(GROUP_POP_diff_snps_vcf_gz), file(GROUP_POP_diff_snps_genotypes), file(GROUP_POP_population_pop_diff_snps_genotypes), file(GROUP_POP_rsIndexPopDiff_txt), file(GROUP_POP_genes_pop_diff_snps_genotypes), file(GROUP_POP_population_gene_pop_diff_snps_genotypes) from rareGlobCommon2_4
////////    output:
////////        set val(GROUP_POP), file("${GROUP_POP}.highDiffPlot.tiff") into plot_rareGlobCommon
////////    script:
////////        cadd_annotations = file(params.cadd_annotations)
////////        GROUP_POP_sample = file(group_pop_data[GROUP_POP][1])
////////        GROUP_POP_datasetAnnotated = file(group_pop_data[GROUP_POP][2])
////////        GROUP_POP_tiff = "${GROUP_POP}.highDiffPlot.tiff"
////////        GROUP_POP_csv = "${GROUP_POP}.highDiffPlot.csv"
////////        println GROUP_POP_diff_snps_genotypes.readLines()
//////////        if (GROUP_POP_diff_snps_genotypes.size() == 0){
//////////            println GROUP_POP_diff_snps_genotypes
//////////        }
////////        template "step8_populationFreq.R"
////////}
//////
//////
//////
////////////workflow.onComplete {
////////////    def subject = 'My pipeline execution'
////////////    def recipient = 'mypandos@gmail.com'
////////////
////////////    ['mail', '-s', subject, recipient].execute() << """
////////////
////////////    Pipeline execution summary
////////////    ---------------------------
////////////    Completed at: ${workflow.complete}
////////////    Duration    : ${workflow.duration}
////////////    Success     : ${workflow.success}
////////////    workDir     : ${workflow.workDir}
////////////    exit status : ${workflow.exitStatus}
////////////    Error report: ${workflow.errorReport ?: '-'}
////////////    """
////////////}
////////
/////////* TODO
////////masterChannel = Channel.from(params.files).println()
////////[condition1, [condition1.bam, input.bam]]
////////[condition2, [condition2.bam, input.bam]]
////////
////////the map must return the new object, hence
////////Channel
////////        .from(params.files)
////////        .map { condition, list ->
////////    def files = list.collect{ file(it)}
////////    return tuple(condition, files)
////////}
////////*/
