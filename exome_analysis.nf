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


// All POP
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
dataset_all_KEY = params.GROUP_POPS.keySet().toList()
dataset_all_VAL = params.GROUP_POPS.values().toList()

GROUP_POPS_ANALYSIS = params.GROUP_POPS_ANALYSIS.split(',')

dataset_pops = [:]
params.POPS.each { entry ->
    dataset_pops[entry.key] = entry.value.split(',')
}
dataset_files_annot = params.dataset_files_annot.split(';')
dataset_files = []
dataset_files_all = []
dataset_files_ = [:]
dataset_files_only = []
params.dataset_files.each{ dataset ->
    dataset_files_only << [dataset.key, file(params.dataset_full_files[dataset.key])]
    CHRMS.each { chrm ->
        if (dataset.key in dataset_files_annot){
            dataset_files << [dataset.key, chrm, file(sprintf(dataset.value, chrm))]
            dataset_files_all << [dataset.key, chrm, file(params.dataset_full_files[dataset.key])]
        }
        if (!(dataset.key in dataset_files_.keySet())){
            dataset_files_[dataset.key] = []
        }
        dataset_files_[dataset.key] << file(sprintf(dataset.value, chrm))
    }
}
dataset_files = Channel.from(dataset_files)
dataset_files_chrs = Channel.from(dataset_files_all)
dataset_files_only_cha = Channel.from(dataset_files_only)

'''
Step: Prepare all data to be used
'''
process dataPreparation {
    tag "dataPrep"
    cpus { 1 * task.attempt }
    memory { 2.GB * task.cpus }
    publishDir "${params.work_dir}/PGX_DATA", overwrite: true, mode:'symlink'

    output:
        file("pgxClinicalrsID.txt") into clinVarData_all
    script:
        """
        # The clinicalAnnotations.csv was downloaded from PharmGKB on 03 July 2017
        # This was used to find rsIDs for genes in this study

        # Select rsIDs from file and use this list to search annotated vcf
        awk -F',' 'NR>1 {gsub(/"/, "", \$1); print \$1}' ${params.clinicalVariants_db} > pgxClinicalrsID.txt
        """
}


"""
Step 1: Split VCF AIBST
"""
process split_vcf{
    tag "split_vcf_${dataset}_${chrm}"
    memory { 2.GB * task.attempt }
    cpus { 1 * task.attempt }
    publishDir "${params.work_dir}/data/${dataset}/VCF/CHRS", overwrite: true, mode:'copy'
    input:
        set val(dataset), val(chrm), file(vcf_file) from dataset_files_chrs
    output:
        set val(dataset), val(chrm), file(vcf_out), file("${vcf_out}.tbi") into split_vcf
    script:
        vcf_file = params.dataset_full_files[dataset]
        vcf_out = "${file(file(vcf_file).baseName).baseName}_chr${chrm}.vcf.gz"
        """
        vcftools \
            --gzvcf ${vcf_file} \
            --chr ${chrm} \
            --recode --stdout \
            | bgzip -c > ${vcf_out}
        bcftools index --tbi ${vcf_out}
        """
}
split_vcf.into{ split_vcf; split_vcf_sub}


"""
Step 1.2: Annotate original VCFs AIBST using snpEff
"""
process annotate_dataset_snpeff_orig{
    tag "snpeff_${file(vcf_file.baseName).baseName}"
    memory { 20.GB * task.attempt }
    cpus { 10 * task.attempt }
    publishDir "${params.work_dir}/data/${dataset}/VCF/", overwrite: true, mode:'symlink'
    input:
        set val(dataset), file(vcf_file) from dataset_files_only_cha
    output:
        set val(dataset), file("${vcf_out}.gz"), file("${base}_snpeff.html") into annotate_dataset_snpeff_org
    script:
        base = file(vcf_file.baseName).baseName
        vcf_out = "${base}_snpeff.vcf"
        """
        ${params.snpEff} \
            ${params.snpEff_human_db} -lof \
            -stats ${base}_snpeff.html \
            -csvStats ${base}_snpeff.csv \
            -dataDir ${params.snpEff_database} \
            -c ${params.snpeff_config} \
            ${vcf_file} -v > ${vcf_out}
        bgzip -f ${vcf_out}
        """
}

annotate_dataset_snpeff_org.into{ annotate_dataset_snpeff_org; annotate_dataset_snpeff_org_1 }
"""
Step 1.3: Annotate original VCFs AIBST using snpEff
"""
process stats_dataset_snpeff_orig{
    tag "stats_${file(vcf_file.baseName).baseName}"
    memory { 20.GB * task.attempt }
    cpus { 10 * task.attempt }
    publishDir "${params.work_dir}/data/${dataset}/VCF/", overwrite: true, mode:'symlink'
    input:
        set val(dataset), file(vcf_file), file(vcf_html) from annotate_dataset_snpeff_org_1
    output:
        set val(dataset), file(vcf_file) into stats_dataset_snpeff_org
    script:
        base = file(vcf_file.baseName).baseName
        """
        bcftools query -l ${vcf_file} > ${base}.samples
        """
}

"""
Step 2: Filter Biallelic sites only in VCFs AIBST
"""
split_vcf.into{ split_vcf; split_vcf_1}
process biall_dataset {
    tag "biall_dataset_${file(vcf_file.baseName).baseName}"
    memory { 20.GB * task.attempt }
    cpus { 10 * task.attempt }
    publishDir "${params.work_dir}/data/${dataset}/VCF/CHRS", overwrite: true, mode: 'symlink'
    input:
        set val(dataset), val(chrm), file(vcf_file), file(vcf_file_tbi) from split_vcf_1
    output:
        set val(dataset), val(chrm), file(vcf_out), file("${vcf_out}.tbi") into biall_dataset
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_biall.vcf.gz"
        """
        bcftools \
            view -m2 -M2 -v snps \
            ${vcf_file} \
            -Oz -o ${vcf_out}
        bcftools index --tbi -f ${vcf_out}
        """
}
biall_dataset.into{ biall_dataset; biall_dataset_sub}


"""
Step 3: Phase VCFs AIBST using eagle
"""
biall_dataset.into{ biall_dataset; biall_dataset_1}
process phase_dataset {
    tag "phase_dataset_${file(vcf_file.baseName).baseName}"
    memory { 20.GB * task.attempt }
    cpus { 10 * task.attempt }
    publishDir "${params.work_dir}/data/${dataset}/VCF/CHRS", overwrite: true, mode: 'copy'
    input:
        set val(dataset), val(chrm), file(vcf_file), file(vcf_file_tbi) from biall_dataset_1
    output:
        set val(dataset), val(chrm), file("${vcf_out}.vcf.gz"), file("${vcf_out}.vcf.gz.tbi") into phase_dataset
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_phased"
        """
        eagle \
            --vcf=${vcf_file} \
            --geneticMapFile=${params.genetic_map} \
            --chrom=${chrm} \
            --genoErrProb 0.003 --pbwtOnly \
            --vcfOutFormat=z \
            --outPrefix=${vcf_out} 2>&1 | tee ${vcf_out}.log
        bcftools index --tbi -f ${vcf_out}.vcf.gz
        """
}
phase_dataset.into{ phase_dataset; phase_dataset_sub}


"""
Step 4: Annotate VCFs AIBST using snpEff
"""
phase_dataset.into{ phase_dataset; phase_dataset_1 }
process annotate_dataset_snpeff{
    tag "snpeff_${file(vcf_file.baseName).baseName}"
    memory { 20.GB * task.attempt }
    cpus { 10 * task.attempt }
    //maxForks 10
    publishDir "${params.work_dir}/data/${dataset}/ANN/CHRS", overwrite: true, mode:'symlink'
    input:
        set val(dataset), val(chrm), file(vcf_file), file(vcf_filetbi) from phase_dataset_1
    output:
        set val(dataset), val(chrm), file("${vcf_out}.gz"), file("${vcf_out}.gz.tbi") into annotate_dataset_snpeff
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_snpeff.vcf"
        """
        ${params.snpEff} \
            ${params.snpEff_human_db} -lof \
            -stats ${dataset}_snpeff.html \
            -csvStats ${dataset}_snpeff.csv \
            -dataDir ${params.snpEff_database} \
            -c ${params.snpeff_config} \
            ${vcf_file} -v > ${vcf_out}
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        """
}
annotate_dataset_snpeff.into{ annotate_dataset_snpeff; annotate_dataset_snpeff_sub}

"""
Step 5: Annotate dbSNP IDs using snpSift
"""
annotate_dataset_snpeff.into { annotate_dataset_snpeff; annotate_dataset_snpeff_1 }
process annotate_dataset_dbsnp{
    tag "dbSNP_${file(vcf_file.baseName).baseName}"
    memory { 20.GB * task.attempt }
    cpus { 10 * task.attempt }
    //maxForks 10
    publishDir "${params.work_dir}/data/${dataset}/ANN/CHRS", overwrite: true, mode:'symlink'
    input:
        set val(dataset), val(chrm), file(vcf_file), file(vcf_file_tbi) from annotate_dataset_snpeff_1
    output:
        set val(dataset), val(chrm), file("${vcf_out}.gz"), file("${vcf_out}.gz.tbi") into annotate_dataset_dbsnp
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_dbsnp.vcf"
        """
        ${params.snpSift} \
            annotate \
            ${params.snpEff_dbsnp_vcf} \
            -c ${params.snpeff_config} \
            -dataDir ${params.snpEff_database} \
            ${vcf_file} > ${vcf_out} -v
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        """
}
annotate_dataset_dbsnp.into{ annotate_dataset_dbsnp; annotate_dataset_dbsnp_sub}


"""
Step 6: Annotate GWAS Catalogue using snpSift
"""
annotate_dataset_dbsnp.into {annotate_dataset_dbsnp; annotate_dataset_dbsnp_1}
process annotate_dataset_gwascat{
    tag "gwascat_${file(vcf_file.baseName).baseName}"
    cpus { 10 * task.attempt }
    memory { 20.GB * task.attempt }
    //maxForks 10
    publishDir "${params.work_dir}/data/${dataset}/ANN", overwrite: true, mode:'symlink'

    input:
        set val(dataset), val(chrm), file(vcf_file), file(vcf_file_tbi) from annotate_dataset_dbsnp_1
    output:
        set val(dataset), val(chrm), file("${vcf_out}.gz"), file("${vcf_out}.gz.tbi") into annotate_dataset_gwascat
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_gwascat.vcf"
        """
        ${params.snpSift} \
            gwasCat \
            -db ${params.gwascat_b37} \
            ${vcf_file} \
            -v -c ${params.snpeff_config} \
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
    memory { 20.GB * task.attempt }
    cpus { 10  * task.attempt }
    //maxForks 10
    time = { 6.hour * task.attempt }
    publishDir "${params.work_dir}/data/${dataset}/ANN/CHRS", overwrite: true, mode:'symlink'
    input:
        set val(dataset), val(chrm), file(vcf_file), file(vcf_file_tbi) from annotate_dataset_gwascat_1
    output:
        set val(dataset), val(chrm), file("${vcf_out}.gz"), file("${vcf_out}.gz.tbi") into annotate_dataset_clinvar
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_clinvar.vcf"
        """
        ${params.snpSift} \
            annotate \
            ${params.clinvar} \
            ${vcf_file} \
            -v -c ${params.snpeff_config} \
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
    memory { 20.GB * task.attempt }
    cpus { 10  * task.attempt }
    time { 6.hour * task.attempt }
    //maxForks 10
    publishDir "${params.work_dir}/data/${dataset}/ANN/CHRS", overwrite: true, mode:'symlink'
    input:
        set val(dataset), val(chrm), file(vcf_file), file(vcf_file_tbi) from annotate_dataset_clinvar_2
    output:
        set val(dataset), val(chrm), file("${vcf_out}.gz"), file("${vcf_out}.gz.tbi") into annotate_dataset_cosmic
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_cosmic.vcf"
        """
        ${params.snpSift} \
            annotate \
            ${params.cosmic} \
            ${vcf_file} \
            -v -c ${params.snpeff_config} \
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
    memory { 20.GB * task.attempt }
    cpus { 10  * task.attempt }
    time { 3.hour * task.attempt }
    //maxForks 10
    publishDir "${params.reference_dir}/pop_mafs/${mafs_dataset}", overwrite: true, mode:'copy'
    input:
        set val(mafs_dataset), val(chrm), file(mafs_file) from mafs_annotations_data_cha
    output:
        set val(mafs_dataset), val(chrm), file(tsv_out) into mafs_annot_dataset
    script:
        tsv_out = "${mafs_dataset}_chr${chrm}_mafs.tsv"
        """
        ${params.homedir}/templates/annotateVCFwithTSV.py \
            --inTSV ${mafs_file} \
            --chrm ${chrm} \
            --outTSV ${tsv_out}
        """
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
    memory { 20.GB * task.attempt }
    time { 5.hour * task.attempt }
    cpus { 10  * task.attempt }
    maxRetries 5
    //maxForks 10
    publishDir "${params.work_dir}/data/${dataset}/ALL/ANN/CHRS", overwrite: true, mode:'symlink'
    input:
        set val(dataset), val(chrm), file(vcf_file), val(mafs_files), val(mafs_dataset) from annotate_dataset_cosmic_2_cha
    output:
        set val(dataset),  val(chrm), file("${vcf_out}.gz"), file("${vcf_out}.gz.tbi") into annotate_mafs_dataset
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_mafs-${mafs_dataset}.vcf"
        """
        gunzip -c ${vcf_file} > ${vcf_file.baseName}
        ${params.homedir}/templates/annotateVCFwithTSV.py \
            --inVCF ${vcf_file.baseName} \
            --inTSV \'${mafs_files}\' \
            --outVCF ${vcf_out}
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        rm -f ${vcf_file.baseName}
        """
}
annotate_mafs_dataset.into {annotate_mafs_dataset; annotate_mafs_dataset_sub}


'''
Step 9.2.1: Annotate VCF for Ancestral Allele (AA) using in-house python script
'''
annotate_mafs_dataset.into {annotate_mafs_dataset; annotate_mafs_dataset_2}
process add_ANC_to_VCF_merge {
    //echo true
    tag "AA_${chrm}_${file(vcf_file.baseName).baseName}"
    memory { 10.GB * task.attempt }
    time { 8.hour * task.attempt }
    publishDir "${params.work_dir}/VCF_ANN/MERGED", overwrite: true, mode:'symlink'
    input:
        set val(dataset),  val(chrm), file(vcf_file), file(vcf_file_tbi) from annotate_mafs_dataset_2
    output:
        set val(dataset), val(chrm), file("${vcf_out}.gz"), file("${vcf_out}.gz.tbi") into vcf_anc_dataset
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_anc.vcf"
        """
        gunzip -c ${vcf_file} > ${vcf_file.baseName}
        ~/.conda/envs/ngs_py27/bin/python2.7 ${params.homedir}/templates/add-ANC-to-vcf_new.py \
            --in ${vcf_file.baseName} \
            --genomedata ${params.genomedata_path} \
            --out ${vcf_out}
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        rm -rf ${vcf_file.baseName}
        """
}


'''
Step 9.3: Reduced VCFs to PGX variants only
'''
vcf_anc_dataset.into {vcf_anc_dataset; vcf_anc_dataset_2}
gene_regions_slop_cha = Channel.fromPath(params.gene_regions_slop)
annotate_mafs_dataset_9_3 = vcf_anc_dataset_2.combine(gene_regions_slop_cha)
process subset_pgx {
    tag "subset_pgx_${file(vcf_file.baseName).baseName}"
    memory { 10.GB * task.attempt }
    cpus { 5 * task.attempt }
    publishDir "${params.work_dir}/data/${dataset}/PGX_ONLY/ANN/CHRS", overwrite: true, mode:'symlink'
    input:
        set val(dataset),  val(chrm), file(vcf_file), file(vcf_file_tbi), file(gene_regions_slop) from annotate_mafs_dataset_9_3
    output:
        set val(dataset), val(chrm), file("${vcf_out}.gz"), file("${vcf_out}.gz.tbi") into subset_pgx
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_pgx.vcf"
        """
        bcftools view ${vcf_file} \
            --regions-file ${gene_regions_slop} \
            -Oz -o ${file(vcf_file.baseName).baseName}.tmp.vcf.gz
        bcftools sort ${file(vcf_file.baseName).baseName}.tmp.vcf.gz -Oz -o ${vcf_out}.gz
        bcftools index --tbi -f ${vcf_out}.gz
        rm ${file(vcf_file.baseName).baseName}.tmp.vcf.gz
        """
}
subset_pgx.into {subset_pgx; subset_pgx_sub}


'''
Step : Split sample files by population
'''
mafs_annotations_data_cha = Channel.from(mafs_annotations_data)
process split_POP_samples {
    tag "split_${pop_sample}"
    memory { 1.GB * task.attempt }
    cpus { 1 * task.attempt }
    time { 1.hour * task.attempt }
    publishDir "${params.work_dir}/data/samples/", overwrite: true, mode:'copy'
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
    cpus { 1 * task.attempt }
    memory { 2.GB * task.cpus }
    publishDir "${params.work_dir}/data/samples/", overwrite: true, mode:'symlink'
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


"""
Step 10a: Combine chromosome VCFs AIBST into one
"""
vcf_anc_dataset.into{ vcf_anc_dataset; vcf_anc_dataset_1}
vcf_anc_dataset_list = [:]
merge_pop_sample.into{ merge_pop_sample; merge_pop_sample_1 }
merge_pop_sample_1_list = merge_pop_sample_1.toSortedList().val
vcf_anc_dataset_1.toSortedList().val.each { dataset, chrm, vcf_file, vcf_file_tbi ->
    if (!(dataset in vcf_anc_dataset_list.keySet())) {
        merge_pop_sample_1_list.each{ dataset1, dataset_sample, dataset_populations, dataset_superPopulations ->
            if (dataset1 == dataset){
                vcf_anc_dataset_list[dataset] = [dataset, file(dataset_sample)]
            }
        }
        vcf_anc_dataset_list[dataset][2] = ''
    }
    vcf_anc_dataset_list[dataset][2] += ' '+vcf_file
}
vcf_anc_dataset_cha = Channel.from(vcf_anc_dataset_list.values())
process concat_dataset {
    tag "concat_dataset_${dataset}"
    cpus { 10 * task.attempt }
    memory { 20.GB * task.attempt }
    publishDir "${params.work_dir}/data/${dataset}/ALL/VCF", overwrite: true, mode: 'symlink'
    input:
        set val(dataset), file(dataset_sample), val(dataset_vcfs) from vcf_anc_dataset_cha
    output:
        set val(dataset), file(dataset_sample), file(vcf_out), file("${vcf_out}.tbi") into concat_dataset_all
    script:
        vcf_out = "${dataset}_annot.vcf.gz"
        """
        bcftools concat \
            ${dataset_vcfs} \
            -Oz -o ${dataset}.tmp.vcf.gz
        ## Recalculate AC, AN, AF
        bcftools +fill-tags  ${dataset}.tmp.vcf.gz -Oz -o  ${dataset}.tmp1.vcf.gz
        #${params.snpSift} sort ${dataset}.tmp1.vcf.gz | bgzip -c > ${vcf_out}
        bcftools sort ${dataset}.tmp1.vcf.gz -Oz -o ${vcf_out}
        bcftools index --tbi -f ${vcf_out}
        rm ${dataset}.tmp*.vcf.gz
        """
}
concat_dataset_all.into{ concat_dataset_all; concat_dataset_sub}


"""
Step 10.3: Annotate PGx VCF AIBST using snpEff
"""
concat_dataset_all.into{ concat_dataset_all; concat_dataset_10_3}
process annotate_dataset_all_snpeff{
    tag "snpeff_${file(vcf_file.baseName).baseName}"
    memory { 20.GB * task.attempt }
    cpus { 10 * task.attempt }
    //maxForks 10
    publishDir "${params.work_dir}/data/${dataset}/ALL/VCF", overwrite: true, mode:'symlink'
    input:
        set val(dataset), file(dataset_sample), file(vcf_file), file(vcf_tbi) from concat_dataset_10_3
    output:
        set val(dataset), file(dataset_sample), file("${vcf_out}.gz"), file("${vcf_out}.gz.tbi"), file("${vcf_out}_snpeff.html"), file("${vcf_out}_snpeff.csv") into annotate_dataset_all_snpeff
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_snpeff.vcf"
        """
        ${params.snpEff} \
            ${params.snpEff_human_db} -lof \
            -stats ${vcf_out}_snpeff.html \
            -csvStats ${vcf_out}_snpeff.csv \
            -dataDir ${params.snpEff_database} \
            -c ${params.snpeff_config} \
            ${vcf_file} -v > ${vcf_out}
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        """
}

"""
Step 10b: Combine chromosome VCFs AIBST into one for PGX variant dataset
"""
subset_pgx.into{ subset_pgx; subset_pgx_1 }
vcf_anc_dataset_pgx_list = [:]
subset_pgx_1.toSortedList().val.each { dataset, chrm, vcf_file, vcf_file_tbi ->
    if (!(dataset in vcf_anc_dataset_pgx_list.keySet())) {
        merge_pop_sample_1_list.each{ dataset1, dataset_sample, dataset_populations, dataset_superPopulations ->
            if (dataset1 == dataset){
                vcf_anc_dataset_pgx_list[dataset] = [dataset, file(dataset_sample)]
            }
        }
        vcf_anc_dataset_pgx_list[dataset][2] = ''
    }
    vcf_anc_dataset_pgx_list[dataset][2] += ' '+vcf_file
}
vcf_anc_dataset_pgx = Channel.from(vcf_anc_dataset_pgx_list.values())
process concat_dataset_pgx {
    tag "concat_dataset_pgx_${dataset}"
    cpus { 10 * task.attempt }
    memory { 20.GB * task.attempt }
    publishDir "${params.work_dir}/data/${dataset}/PGX_ONLY/VCF", overwrite: true, mode: 'symlink'
    input:
        set val(dataset), file(dataset_sample), val(dataset_vcfs) from vcf_anc_dataset_pgx
    output:
        set val(dataset), file(dataset_sample), file(vcf_out), file("${vcf_out}.tbi") into concat_dataset_pgx
    script:
        vcf_out = "${dataset}_annot_pgx.vcf.gz"
        """
        bcftools concat \
            ${dataset_vcfs} \
            -Oz -o ${dataset}.tmp.vcf.gz
        ## Recalculate AC, AN, AF
        bcftools +fill-tags  ${dataset}.tmp.vcf.gz -Oz -o  ${dataset}.tmp1.vcf.gz
        #${params.snpSift} sort ${dataset}.tmp1.vcf.gz | bgzip -c > ${vcf_out}
        bcftools sort ${dataset}.tmp1.vcf.gz -Oz -o ${vcf_out}
        bcftools index --tbi -f ${vcf_out}
        rm ${dataset}.tmp*.vcf.gz
        """
}


"""
Step 10.3: Annotate PGx VCF AIBST using snpEff
"""
concat_dataset_pgx.into{ concat_dataset_pgx; concat_dataset_pgx_1}
process annotate_dataset_pgx_snpeff{
    tag "snpeff_${file(vcf_file.baseName).baseName}"
    memory { 20.GB * task.attempt }
    cpus { 10 * task.attempt }
    //maxForks 10
    publishDir "${params.work_dir}/data/${dataset}/PGX_ONLY/VCF", overwrite: true, mode:'symlink'
    input:
        set val(dataset), file(dataset_sample), file(vcf_file), file(vcf_tbi) from concat_dataset_pgx_1
    output:
        set val(dataset), file(dataset_sample), file("${vcf_out}.gz"), file("${vcf_out}.gz.tbi") into annotate_dataset_pgx_snpeff
    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_snpeff.vcf"
        """
        ${params.snpEff} \
            ${params.snpEff_human_db} -lof \
            -stats ${vcf_out}_snpeff.html \
            -csvStats ${vcf_out}_snpeff.csv \
            -dataDir ${params.snpEff_database} \
            -c ${params.snpeff_config} \
            ${vcf_file} -v > ${vcf_out}
        bgzip -f ${vcf_out}
        bcftools index --tbi -f ${vcf_out}.gz
        """
}

//// ANALYSIS START HERE

////// Generate dataflow for use in processes
concat_dataset_pgx.into{ concat_dataset_pgx; concat_dataset_pgx_2}
group_vcf_with_pop_sample = [] // For POP, dataset, Full GROUP_VCF, GROUP_SAMPLE, dataset_populations, dataset_superPopulations
group_pop_data = [:] // contains annotated vcfs of all
concat_dataset_pgx_2.toSortedList().val.each { dataset, dataset_sample, dataset_vcf, dataset_vcf_tbi ->
    if ( dataset in dataset_files_annot) {
        data = file(dataset_sample).readLines().collect { it.split()[1] }.unique()
        data.each { POP ->
            group_vcf_with_pop_sample << [POP, dataset, dataset_vcf, dataset_sample]
        }
        // All dataset vcf and sample files
        if ( !(dataset in group_pop_data.keySet()) ){
            group_pop_data[dataset] = [dataset_vcf, dataset_sample]
        }
    }
}

// TODO Start HERE
'''
Step 11b: Analysis of LoF
'''
concat_dataset_pgx.into{ concat_dataset_pgx; concat_dataset_pgx_1}
process lof_analysis_pgx{
    tag "lof_analysis_${dataset}"
    cpus { 2 * task.attempt }
    memory { 4.GB * task.cpus }
    publishDir "${params.work_dir}/data/PGX_ONLY/LOF/${dataset}", overwrite: true, mode:'symlink'
    // publishDir "${params.work_dir}/REPORTS/PGX_ONLY/LOF/${dataset}", overwrite: true, mode:'symlink'
    input:
        set val(dataset), file(dataset_sample), file(dataset_vcf), file(dataset_vcf_tbi) from concat_dataset_pgx_1
    output:
        set val(dataset), file("${dataset_out}.lof${lof_cutoff1}.vcf.gz"), file("${dataset_out}.lof${lof_cutoff1}.vcf.gz.tbi"), file("${dataset_out}.lof${lof_cutoff1}.INFO"), file("${dataset_out}.lof${lof_cutoff1}.LoFConsequencesIDacGENE"), file("${dataset_out}.lof${lof_cutoff1}.LofPerSample.012"), file("${dataset_out}.lof${lof_cutoff1}.LofPerSample.012.indv"), file(dataset_sample) into lof_analysis_pgx
    script:
        dataset_out = "${file(dataset_vcf.baseName).baseName}"
        lof_cutoff = '0.4'
        lof_cutoff1 = lof_cutoff.replaceAll('\\.', '')
        """
        ## Filter VCf for LOF > 1
        zcat ${dataset_vcf} | ${params.snpSift} filter "LOF[*].PERC >= ${lof_cutoff}" | bgzip -c > ${dataset_out}.lof.tmp.vcf.gz
        ## Remove multi-allelic variants
        ## Remove singleton and doubleton
        ## Remove variants with a Lof_flag
        ## TODO Curate variants with >30% allele freq: remove rs11356919 (AF 69%) due to lack of evidence
        vcftools --gzvcf ${dataset_out}.lof.tmp.vcf.gz \
           --max-alleles 2 --recode-INFO-all --recode --stdout | \
        grep -v 'SINGLE_EXON|NAGNAG_SITE|PHYLOCSF_WEAK|PHYLOCSF_UNLIKELY_ORF|PHYLOCSF_TOO_SHORT' | \
            bgzip -c > ${dataset_out}.lof${lof_cutoff1}.vcf.gz
        # CALCULATE WHICH GENES HARBOUR LoF MUTATIONS and HOW MANY
        # Extract information from vcf info file
        zcat ${dataset_out}.lof${lof_cutoff1}.vcf.gz | ${params.snpSift} extractFields - -e "." \
            CHROM POS REF ALT AC ANN[0].GENE KG_AF KG_AFR_AF KG_EUR_AF KG_AMR_AF KG_EAS_AF gnomAD_AF gnomAD_AFR_AF gnomAD_FIN_AF ExAC_AF ExAC_AFR_AF AGVP_AF SAHGP_AF "ANN[0].EFFECT" CLNDN CDS GWASCAT_TRAIT GWASCAT_P_VALUE GWASCAT_PUBMED_ID\
            > ${dataset_out}.lof${lof_cutoff1}.INFO
        bcftools index --tbi -f ${dataset_out}.lof${lof_cutoff1}.vcf.gz
        # Create a unique ID for each of the LoF variants
        awk 'NR>1 {print \$1":"\$2":"\$4}' ${dataset_out}.lof${lof_cutoff1}.INFO > ${dataset_out}.lof${lof_cutoff1}.UniqueID
        # Get allele counts for each of the genes
        awk 'NR>1 {print \$0}' ${dataset_out}.lof${lof_cutoff1}.INFO | cut -c5- > ${dataset_out}.lof${lof_cutoff1}.AC
        # Get the relevant gene names for each of the LoF variants
        #awk -F' ' 'NR>1 {print \$NF}' ${dataset_out}.lof${lof_cutoff1}.INFO > ${dataset_out}.lof${lof_cutoff1}.Gene
        # Combine the Unique ID, AC/AF and Gene files
        paste ${dataset_out}.lof${lof_cutoff1}.UniqueID ${dataset_out}.lof${lof_cutoff1}.AC > ${dataset_out}.lof${lof_cutoff1}.LoFConsequencesIDacGENE
        #rm ${dataset_out}.lof${lof_cutoff1}.UniqueID ${dataset_out}.lof${lof_cutoff1}.AC ${dataset_out}.lof${lof_cutoff1}.Gene
        # CALCULATE NUMBER OF LoF / SAMPLE
        # Calculate the number of non-ref alleles per sample
        vcftools --gzvcf ${dataset_vcf} --012 --out ${dataset_out}.lof${lof_cutoff1}.LofPerSample
        """
}


'''
Step 11c: Plot analysis of LoF
'''
lof_analysis_pgx.into{ lof_analysis_pgx; lof_analysis_pgx_1 }
process plot_lof_analysis_pgx{
    tag "plot_lof_${dataset}"
    cpus { 2 * task.attempt }
    memory { 2.GB * task.cpus }
    publishDir "${params.work_dir}/data/PGX_ONLY/LOF/${dataset}", overwrite: true, mode:'copy'
    publishDir "${params.work_dir}/REPORTS/${dataset}/PGX_ONLY/LOF/", overwrite: true, mode:'copy', pattern: '*tiff'
    input:
        set val(dataset), file(dataset_vcf), file(dataset_vcf_tbi), file(dataset_INFO), file(dataset_LoFConsequencesIDacGENE), file(dataset_LofPerSample_012), file(dataset_LofPerSample_012_indv), file(dataset_sample) from lof_analysis_pgx_1
    output:
        set val(dataset), file(dataset_pgxLoFPerGeneCount), file(dataset_pgxLoFPerCombinedAF), file(dataset_allLofPop) into plot_lof_analysis_pgx
    script:
        dataset_allLofPop = "${dataset}_allLofPop.csv"
        dataset_pgxLoFPerGeneCount = "${dataset}_pgxLoFPerGeneCount.tiff"
        dataset_pgxLoFPerCombinedAF = "${dataset}_pgxLoFPerCombinedAF.tiff"
        template "step11b_lof_analysis.R"
}


'''
Step 12: Population Frequency analysis
'''
//merge_pop_vcf_all.into {merge_pop_vcf_all; merge_pop_vcf__lof_analysis}
//annot_vep_merge_pop_all.into { annot_vep_merge_pop_all; annot_vep_merge_pop__popFrq }
concat_dataset_pgx.into{ concat_dataset_pgx; concat_dataset_pgx_2}
//merge_pop_sample__cha.into { merge_pop_sample__cha; merge_pop_sample__popFrq_cha }
process popFrq_1 {
    tag "popFreq_1_${dataset}"
    cpus { 2 * task.attempt }
    memory { 2.GB * task.cpus }
    publishDir "${params.work_dir}/data/MERGED_POP/FRQ/${dataset}", overwrite: true, mode:'symlink'

    input:
        set val(dataset), file(dataset_sample), file(dataset_vcf), file(dataset_vcf_tbi) from concat_dataset_pgx_2
    output:
        set val(dataset), file("${dataset}.biall.recode.vcf.gz"), file(dataset_sample), file("${dataset}.populations.txt"), file("${dataset}.superPopulations.txt") into popFrq_1_all
    script:
        """
        # Make files with unique (super)population
        awk '{print \$2}' ${dataset_sample} | sort | uniq > ${dataset}.populations.txt
        awk '{print \$3}' ${dataset_sample} | sort | uniq > ${dataset}.superPopulations.txt
        # Calculate bi-allelic SNPs
        vcftools --gzvcf ${dataset_vcf} \
            --min-alleles 2 --max-alleles 2 \
            --recode-INFO-all --recode --out ${dataset}.biall
        bgzip -f ${dataset}.biall.recode.vcf
        """
}
popFrq_1_all.into { popFrq_1_all; popFrq_1__sub }


popFrq_1_all.into { popFrq_1_all; popFrq_1__popFrq_2 }
pop_samples_per_group_1 = [] // For POP, GROUP_POP, Biall GROUP_VCF, GROUP_SAMPLE, GROUP_POP_populations, GROUP_POP_superPopulations
popFrq_1_list = popFrq_1__popFrq_2.toSortedList()
popFrq_1_list.val.each { dataset, dataset_biall_vcf, dataset_sample, dataset_populations, dataset_superPopulations ->
    data = file(dataset_populations).readLines().toSorted()
    data.each { POP ->
        pop_samples_per_group_1 << [POP, dataset, dataset_biall_vcf, dataset_sample, dataset_populations, dataset_superPopulations ]
    }
}
pop_samples_per_group_1_chan = Channel.from(pop_samples_per_group_1)

'''
Step 12.1: Compute frequencies per population
           - Create a file with each sample in a given population and
           - Calculate allele frequencies of bi-allelic SNPS for individual populations based on the derived allele
'''
process popFrq_2 {
    tag "popFreq_2_${POP}"
    cpus { 2 * task.attempt }
    memory { 2.GB * task.cpus }
    publishDir "${params.work_dir}/data/MERGED_POP/FRQ/${GROUP_POP}", overwrite: true, mode:'symlink'

    input:
        set val(POP), val(GROUP_POP), file(GROUP_POP_biall_vcf), file(GROUP_POP_sample), file(GROUP_POP_populations), file(GROUP_POP_superPopulations) from pop_samples_per_group_1_chan
    output:
        set val(POP), val(GROUP_POP), file("${POP}.samples"), file("${POP}.vcf.gz"), file("${POP}.Ind.frq") into popFrq_2_all
    script:
        """
        grep ${POP} ${GROUP_POP_sample} | cut -f1 > ${POP}.samples
        vcftools \
            --gzvcf ${GROUP_POP_biall_vcf} \
            --keep ${POP}.samples \
            --recode-INFO-all --recode --stdout | \
        bgzip -c > ${POP}.tmp.vcf.gz
        ## Recalculate AC, AN, AF
        bcftools +fill-tags  ${POP}.tmp.vcf.gz -Oz -o  ${POP}.vcf.gz
        ## Compute frequency
        ##vcftools --gzvcf ${POP}.vcf.gz --freq2 --derived --out ${POP}.Ind
        plink2 --vcf ${POP}.vcf.gz --freq --set-missing-var-ids @:# --out ${POP} || true
        echo "rsID\tREF\tALT\t${POP}_Alt_FREQ" > ${POP}.Ind.frq
        awk -F" " 'NR>1 {print \$2"\t"\$3"\t"\$4"\t"\$5}' ${POP}.frq >> ${POP}.Ind.frq
        rm -f ${POP}.tmp.vcf.gz
        """
}


'''
Step 12.2: Compute global frequencies for dataset
            - Calculate the global allele frequency of bi-allelic SNPs
'''
popFrq_1_all.into { popFrq_1_all; popFrq_1__popFrq_3 }
process popFrq_3 {
    tag "popFreq_3_${dataset}"
    cpus { 2 * task.attempt }
    memory { 2.GB * task.cpus }
    publishDir "${params.work_dir}/data/MERGED_POP/FRQ/${dataset}", overwrite: true, mode:'symlink'

    input:
        set val(dataset), file(dataset_biall_vcf), file(dataset_sample), file(dataset_populations), file(dataset_superPopulations) from popFrq_1__popFrq_3
    output:
        set val(dataset), file("${dataset}.biall.recode.vcf.gz"), file(dataset_sample), file("${dataset}.biall.INFO.rsID") into popFrq_3_all
    script:
        """
        # Calculate the global allele frequency of bi-allelic SNPs
        #vcftools --gzvcf ${dataset_biall_vcf} --freq2 --derived --out ${dataset}.all
        # Create file with allele count/number/freq INFO fields from vcf
        zcat ${dataset_biall_vcf} \
            | ${params.snpSift} extractFields - -e "." CHROM POS ID REF ALT AC AN AF MAF "ANN[0].GENE" KG_AF KG_AFR_AF ExAC_AF ExAC_AFR_AF gnomAD_AF gnomAD_AFR_AF AGVP_AF "ANN[0].EFFECT" \
            > ${dataset}.biall.INFO
        echo "UniqueID\tCHROM\tPOS\tID\tREF\tALT\tAC\tAN\tAF\tMAF\tGene\tKG_Alt_FRE\tKG_AFR_Alt_FREQ\tExAC_Alt_FRE\tExAC_AFR_Alt_FRE\tgnomAD_Alt_FRE\tgnomAD_AFR_Alt_FREQ\tAGVP_Alt_FREQ\tEFFECT" > ${dataset}.biall.INFO.rsID
        awk 'NR>1 {print \$1":"\$2":"\$5"\t"\$0}' ${dataset}.biall.INFO >> ${dataset}.biall.INFO.rsID
        #vcftools --gzvcf ${dataset_biall_vcf} \
        #    --get-INFO AC --get-INFO AN --get-INFO AF --get-INFO MAF \
        #    --out ${dataset}.biall
        # Create file with rsIDs for vcf (double hash saves header)
        #zcat ${dataset}.biall.vcf.gz | grep -v '##' | awk '{print \$3}' > ${dataset}.biall.rsID
        # Combine INFO and rsID files
        #paste ${dataset}.biall.INFO ${dataset}.biall.rsID > ${dataset}.biall.INFO.rsID
        """
}


'''
Step 12.3: Compute global frequencies for dataset
'''
popFrq_3_all.into{popFrq_3_all; popFrq_3__2}
popFrq_2_all.into{popFrq_2_all; popFrq_2__2}
Ind_group = [:]
temp = [:]
temp1 = [:]
popFrq_2__2.toSortedList().val.each { POP, dataset, POP_sample, POP_vcf, POP_Ind_frq ->
    if ( !(dataset in Ind_group.keySet()) ){
        Ind_group[dataset] = [dataset]
        temp[dataset] = []
        temp1[dataset] = []
    }
    temp[dataset] << POP+":"+POP_Ind_frq
    temp1[dataset] << POP
}
popFrq_3__2.toSortedList().val.each{ dataset, dataset_biall_vcf, GROUP_sample, GROUP_biall_info_rsid ->
    Ind_group[dataset] << GROUP_biall_info_rsid
    Ind_group[dataset] << temp[dataset].toSorted().collect{ it.split(':')[1] }.join(' ')
    Ind_group[dataset] << temp1[dataset].toSorted().join(' ')
}

process popFrq_3_1 {
    //echo true
    tag "popFreq_3_1_${dataset}"
    cpus { 2 * task.attempt }
    memory { 2.GB * task.cpus }
    publishDir "${params.work_dir}/data/MERGED_POP/FRQ/${dataset}", overwrite: true, mode:'symlink'
    input:
        val(dataset) from Ind_group.keySet()
    output:
        set val(dataset), file("${dataset}_dataset.csv"), file("${dataset}_datasetAnnotated.csv") into popFrq_3_1_all
    script:
        GROUP_biall_info_rsid = Ind_group[dataset][1]
        POP_Ind_frq_files = file(Ind_group[dataset][2])
        POPs = Ind_group[dataset][3]
        template "create_datasetAnnotated.R"
}


popFrq_3_1_all.into{ popFrq_3_1_all; popFrq_3_1_1}
//// Add to group data array that contains group:vcf,sapmle,annotatedDatapopulation groups
popFrq_3_1_1.toSortedList().val.each{ dataset, dataset_dataset, dataset_datasetAnnotated ->
    group_pop_data[dataset] << dataset_datasetAnnotated
}


'''
Step 12.4: Plot total polymosphism per population
'''
popFrq_3_1_all.into{ popFrq_3_1_all; popFrq_3_1_2 }
process plot_PolymorphicPerPopulation{
    tag "plot_poly_${dataset}"
    cpus { 2 * task.attempt }
    memory { 2.GB * task.cpus }
    publishDir "${params.work_dir}/data/PGX_ONLY/FRQ/${dataset}", overwrite: true, mode:'copy'
    publishDir "${params.work_dir}/REPORTS/${dataset}/PGX_ONLY/FRQ/", overwrite: true, mode:'copy'
    input:
        set val(dataset), file(dataset_dataset), file(dataset_datasetAnnotated) from popFrq_3_1_2
    output:
        set val(dataset), file(dataset_pgxPolymorphicPerPopulation) into plot_PolymorphicPerPopulation
    script:
        dataset_sample = group_pop_data[dataset][1]
        dataset_pgxPolymorphicPerPopulation = "${dataset}_pgxPolymorphicPerPopulation.tiff"
        template "step12_4_polymorphic_per_population.R"
}
plot_PolymorphicPerPopulation.into{ plot_PolymorphicPerPopulation; plot_PolymorphicPerPopulation_sub}


'''
Step 12.5: List singleton markers
'''
popFrq_1_all.into { popFrq_1_all; popFrq_1__popFrq_4 }
process popFrq_4 {
    tag "popFreq_4_${dataset}"
    cpus { 2 * task.attempt }
    memory { 2.GB * task.cpus }
    publishDir "${params.work_dir}/data/MERGED_POP/FRQ/${dataset}", overwrite: true, mode:'symlink'
    input:
        set val(dataset), file(dataset_biall_vcf), file(dataset_sample), file(dataset_populations), file(dataset_superPopulations) from popFrq_1__popFrq_4
    output:
        set val(dataset), file(dataset_sample), file("${dataset}.singletons.per.sample") into popFrq_4_all
    script:
        """
        # List singleton markers
        vcftools --gzvcf ${dataset_biall_vcf} --singletons --out ${dataset}
        # Number of singletons per sample
        awk '{print \$1" "\$2}' ${dataset_sample} > ${dataset}.txt
        echo "singletons" > ${dataset}.singleton.count
        while read SAMPLE; do
            grep \$SAMPLE ${dataset}.singletons | wc -l;
        done < ${dataset}.txt >> ${dataset}.singleton.count
        echo "sample pop\n \$(cat ${dataset}.txt)" > ${dataset}.txt
        paste ${dataset}.txt ${dataset}.singleton.count > ${dataset}.singletons.per.sample
        rm ${dataset}.singleton.count
        """
}
popFrq_4_all.into{popFrq_4_all; popFrq_4__sub}



'''
Step 12.6: Plot average singletons per population
'''
popFrq_4_all.into{ popFrq_4_all; popFrq_4_1 }
process plot_pgxAverageSingletonsPerPopulation{
    tag "plot_singl_${dataset}"
    cpus { 2 * task.attempt }
    memory { 2.GB * task.cpus }
    publishDir "${params.work_dir}/data/PGX_ONLY/FRQ/${dataset}", overwrite: true, mode:'copy'
    publishDir "${params.work_dir}/REPORTS/${dataset}/PGX_ONLY/FRQ/", overwrite: true, mode:'copy'
    input:
        set val(dataset), file(dataset_sample), file(dataset_singletons_per_sample) from popFrq_4_1
    output:
        set val(dataset), file(dataset_pgxAverageSingletonsPerPopulation) into plot_pgxAverageSingletonsPerPopulation
    script:
        dataset_pgxAverageSingletonsPerPopulation = "${dataset}_pgxAverageSingletonsPerPopulation.tiff"
        template "step12_6_average_singletons_per_population.R"
}


'''
Step 13: FST analysis
'''
group_vcf_with_pop_sample_cha = Channel.from( group_vcf_with_pop_sample )

// Building dataflow for dataset, annotated VCF, sample file
concat_dataset_pgx.into { concat_dataset_pgx; concat_dataset_pgx_3 }
pop_samples_per_group = []
pop2_per_group = []
split_POP_samples.into{split_POP_samples; split_POP_samples_2}
split_POP_samples_list = split_POP_samples_2.toSortedList().val
pop2_per_group_cha = concat_dataset_pgx_3.flatMap { dataset, dataset_sample, dataset_vcf, dataset_vcf_tbi ->
    if ( dataset in dataset_files_annot){
        pop_samples_per_group << [dataset, dataset_vcf, dataset_sample]
        tmp = []
        for (data1 in split_POP_samples_list){
            POP1 = data1[0]
            SAMPLE1 = data1[1]
            for (data2 in split_POP_samples_list) {
                POP2 = data2[0]
                SAMPLE2 = data2[1]
                if ( POP1 != POP2 && !([POP2, POP1] in tmp) ){
                    pop2_per_group << [ dataset, [POP1, POP2], file(SAMPLE1), file(SAMPLE2), dataset_vcf, dataset_sample ]
                    tmp << [POP1, POP2]
                }
            }
        }
    }
    return pop2_per_group
}


/// Have to do FST analysis for all variants and compare with PGx
process fst_analysis {
    //echo true
    tag "fst_${POPS[0]}_${POPS[1]}"
    cpus { 2 * task.attempt }
    memory { 2.GB * task.cpus }
    publishDir "${params.work_dir}/data/PGX_ONLY/FST/${dataset}", overwrite: true, mode:'symlink'
    //publishDir "${params.work_dir}/REPORT/PGX_ONLY/FST/${dataset}", overwrite: true, mode:'copy', pattern:'*.tiff'
    input:
        set val(dataset), val(POPS), file(pop1_sample), file(pop2_sample), file(dataset_vcf), file(dataset_sample) from pop2_per_group_cha
    output:
        set val(dataset), val(POPS), file("${fst_basename}.weir.fst"), file(fst_log) into fst_analysis_all
    script:
        fst_basename = "${dataset}.${POPS[0]}_${POPS[1]}"
        fst_log = "${dataset}.${POPS[0]}_${POPS[1]}.fst.log"
        """
        ## Filter for synonymous variants
        zcat ${dataset_vcf} | ${params.snpSift} filter "ANN[*].EFFECT has 'synonymous_variant'" | bgzip -c > ${dataset}.fst.syn.vcf.gz
        ## Compute Fst for pair of populations, here the log file will be used to extract the weighted mean fst
        vcftools \
            --gzvcf ${dataset}.fst.syn.vcf.gz \
            --weir-fst-pop ${pop1_sample} \
            --weir-fst-pop ${pop2_sample} \
            --out ${fst_basename} \
            2> ${fst_log}
        """
}
fst_analysis_all.into{fst_analysis_all; fst_analysis__sub}
fst_analysis__sub.subscribe {
      //println "|-- Finished for ${it[0]}, ${it[-1]}"
}


'''
Step 13.1: Combine FST results
'''
fst_analysis_all.into { fst_analysis_all; fst_analysis__combine_fst }
combine_fst_data = [:]
combine_weir_fst_data = [:]
fst_analysis_list = fst_analysis__combine_fst.toSortedList().val
fst_analysis_list.each{ dataset, pop2, weir_fst, fst_log ->
    if ( dataset in combine_fst_data.keySet() ){
        pop = []
        pop << pop2.join('__')
        pop << weir_fst
        pop << fst_log
        combine_fst_data[dataset] << pop.join('__')
    }
    else {
        pop = []
        pop << pop2.join('__')
        pop << weir_fst
        pop << fst_log
        combine_fst_data[dataset] = [pop.join('__')]
    }
}

// Combine pop fst
combine_fst_data_cha = Channel.from( combine_fst_data.keySet() )
process combine_fst_analysis {
    //echo true
    tag "comb_fst_${dataset}"
    cpus { 1 * task.attempt }
    memory { 1.GB * task.cpus }
    publishDir "${params.work_dir}/data/PGX_ONLY/FST/${dataset}", overwrite: true, mode:'symlink'

    input:
        val(dataset) from combine_fst_data_cha
    output:
        set val(dataset), file(fst_out) into combine_fst_analysis
    script:
        datas = combine_fst_data[dataset].join(' ')
        fst_out = "${dataset}.weighted.fst.estimates.txt"
        """
        for data in ${datas};do
            ## Split string into array
            data_array=(\${data//__/ })
            pop1=\${data_array[0]}
            pop2=\${data_array[1]}
            log_=\${data_array[3]}
            # Get mean Fst estimates for each of the comparisons
            weir_fst=`grep "Weir and Cockerham weighted Fst estimate" \$log_ | awk -F":" '{print \$NF}'`
            echo \"\$pop1 \$pop2 \$weir_fst\" >> ${fst_out}
        done
        """
}


// combine SNP weir fst
combine_weir_fst_data_cha = Channel.from( combine_fst_data.keySet() )
process combine_weir_fst_analysis {
    tag "comb_weir_fst_${dataset}"
    cpus { 1 * task.attempt }
    memory { 1.GB * task.cpus }
    publishDir "${params.work_dir}/data/PGX_ONLY/FST/${dataset}", overwrite: true, mode:'symlink'

    input:
        val(dataset) from combine_weir_fst_data_cha
    output:
        set val(dataset), file("${fst_out}.csv"), file("${fst_out}_HighDiff.csv"), file("${fst_out}_HighDiff_pos.csv"), file("${fst_out}_HighDiff.vcf.gz"), file("${fst_out}_HighDiff.ann.csv") into combine_weir_fst_analysis
    script:
        datas = combine_fst_data[dataset].join(' ')
        dataset_vcf = file(group_pop_data[dataset][0])
        fst_out = "${dataset}.weighted_weir-fst_estimates"
        """
        python2.7 ${params.homedir}/templates/combine_weir_fst.py  \
            --fst_input '${datas}' \
            --fst_output ${fst_out}
        vcftools \
            --gzvcf ${dataset_vcf} \
            --positions ${"${fst_out}_HighDiff_pos.csv"} \
            --recode-INFO-all --recode --stdout \
            | gzip -c > ${"${fst_out}_HighDiff.vcf.gz"}
        zcat ${fst_out}_HighDiff.vcf.gz \
            | ${params.snpSift} extractFields - -e "." CHROM POS ID REF ALT AC AF MAF "ANN[0].GENE" KG_AF KG_AFR_AF ExAC_AF ExAC_AFR_AF gnomAD_AF gnomAD_AFR_AF SAHGP_AF "ANN[0].EFFECT" CLNDN CDS GWASCAT_TRAIT GWASCAT_P_VALUE GWASCAT_PUBMED_ID\
            > ${fst_out}_HighDiff.ann.csv
        """
}
combine_weir_fst_analysis.into{combine_weir_fst_analysis; combine_weir_fst_analysis__sub}

//
//'''
//Filter dataset VCF for highly differentiated SNPs
//'''
//combine_weir_fst_analysis.into { combine_weir_fst_analysis; combine_weir_fst_analysis_1 }
//process vcf_highDiff_weir_fsf {
//    tag "vcf_highDiff_weir_fst_${dataset}"
//    cpus { 1 * task.attempt }
//    memory { 1.GB * task.cpus }
//    publishDir "${params.work_dir}/data/PGX_ONLY/FST/${dataset}", overwrite: true, mode:'symlink'
//
//    input:
//        set val(dataset), file(fst_csv), file(fst_HighDiff_csv) from combine_weir_fst_analysis_1
//    output:
//        set val(dataset), file("${dataset}.csv"), file("${dataset}_HighDiff.csv") into vcf_highDiff_weir_fsf
//    script:
//        dataset_vcf = file(group_pop_data[dataset][0])
//        //println dataset_vcf
////        fst_out = "${dataset}.weighted_weir-fst_estimates"
//        """
//        """
//}
//vcf_highDiff_weir_fsf.into{vcf_highDiff_weir_fsf; vcf_highDiff_weir_fsf_sub}


'''
Step 13.2: transform FST results to matrix
'''
combine_fst_analysis.into{ combine_fst_analysis; combine_fst_analysis_1}
process matrix_combine_fst_analysis {
    tag "matrix_comb_fst_${dataset}"
    cpus { 1 * task.attempt }
    memory { 1.GB * task.cpus }
    publishDir "${params.work_dir}/data/PGX_ONLY/FST/${dataset}", overwrite: true, mode:'copy'
    publishDir "${params.work_dir}/REPORTS/${dataset}/PGX_ONLY/FST/", overwrite: true, mode:'copy'

    input:
        set val(dataset), file(fst_in) from combine_fst_analysis_1
    output:
        set val(dataset), file(fst_out) into matrix_combine_fst_analysis
    script:
        fst_out = "${dataset}.weighted.fst.estimates.matrix.csv"
        """
        python2.7 ${params.homedir}/templates/convert_fst_result_to_matrix.py \
            --fts_2by2_input ${fst_in} \
            --fst_matrix_output ${fst_out}
        """
}
matrix_combine_fst_analysis.into{matrix_combine_fst_analysis; matrix_combine_fst_analysis__sub}


'''
Step 14.1: Clinical variants analysis. Generate VCF of clinical variants only per population
'''
group_vcf_with_pop_sample_cha.into { group_vcf_with_pop_sample_cha; group_vcf_with_pop_sample_cha__clinVar }
clinVarData_all.into { clinVarData_all; clinVarData__clinVar }
process clinVar {
    //echo true
    tag "clinVar_${POP}"
    cpus { 2 * task.attempt }
    memory { 2.GB * task.cpus }
    publishDir "${params.work_dir}/data/PGX_ONLY/CLINICAL/${GROUP_POP}", overwrite: true, mode:'symlink'

    input:
        set val(POP), val(GROUP_POP), file(GROUP_POP_vcf), file(GROUP_POP_sample) from group_vcf_with_pop_sample_cha__clinVar
        file(pgxClinicalrsID) from clinVarData__clinVar
    output:
        set val(POP), val(GROUP_POP), file("${POP}.clinical.vcf.gz"), file("${POP}.clinical.ind.frq"), file("${POP}.rsIndexClin.txt") into clinVar_all
    script:
        """
        grep ${POP} ${GROUP_POP_sample} | cut -f1 > ${POP}.samples
        # Select rsIDs from file and use this list to search annotated vcf
        # Extract sites for population
        vcftools --gzvcf ${GROUP_POP_vcf} \
            --snps ${pgxClinicalrsID} \
            --keep ${POP}.samples \
            --recode-INFO-all \
            --recode --stdout | \
            bgzip -c > ${POP}.clinical.vcf.gz
        # Generate sites frequency for population
        vcftools --gzvcf ${POP}.clinical.vcf.gz \
            --freq2 --derived \
            --out ${POP}.clinical.ind
        # Create an index from vcf file that can be used to map back rsIDs as freq only outputs coordinates
        zcat ${POP}.clinical.vcf.gz | awk '!/^(\$|#)/ {print \$3"\t"\$1":"\$2}' - > ${POP}.rsIndexClin.txt
        """
}


'''
Step 14.2: Clinical variants analysis. Calculate the number of non-ref alleles per sample for the whole dataset (AIBST)
'''
clinVarData_all.into { clinVarData_all; clinVarData__clinVar1 }
process clinVar1 {
    tag "clinVar1_${dataset}"
    cpus { 2 * task.attempt }
    memory { 2.GB * task.cpus }
    publishDir "${params.work_dir}/data/PGX_ONLY/CLINICAL/${dataset}", overwrite: true, mode:'symlink'
    input:
        val(dataset) from dataset_files_annot
        file(pgxClinicalrsID) from clinVarData__clinVar1
    output:
        set val(dataset), file("${dataset}.clinical.vcf.gz"), file("${dataset}.ClinPerSample.012"), file("${dataset}.ClinPerSample.012.indv"), file("${dataset}.ClinPerSample.012.pos"), file("pgxClinicalrsID.txt") into clinVar1_all
    script:
        dataset_vcf = file(group_pop_data[dataset][0])
        """
        # Calculate the number of non-ref alleles per sample
        vcftools --gzvcf ${dataset_vcf} \
            --snps ${pgxClinicalrsID} \
            --recode-INFO-all --recode --stdout | \
            bgzip -c > ${dataset}.clinical.vcf.gz
        vcftools --gzvcf ${dataset}.clinical.vcf.gz \
            --012 --out ${dataset}.ClinPerSample
        """
}

'''
Step 14.3: Plot Clinical variants analysis. (AIBST)
'''
clinVar1_all.into { clinVar1_all; clinVar1_all_1 }
plot_lof_analysis_pgx.into {plot_lof_analysis_pgx; plot_lof_analysis_pgx_14_3}
process clinVar_plot {
    tag "clinVar_plot_${dataset}"
    cpus { 2 * task.attempt }
    memory { 2.GB * task.cpus }
    publishDir "${params.work_dir}/data/PGX_ONLY/CLINICAL/${dataset}", overwrite: true, mode:'symlink'
    publishDir "${params.work_dir}/REPORTS/${dataset}/PGX_ONLY/CLINICAL", overwrite: true, mode:'copy'
    input:
        set val(dataset), file(dataset_clinical_vcf), file(ClinPerSample_012_file), file(ClinPerSample_012_indv_file), file(ClinPerSample_012_pos_file), file(pgxClinicalrsID_txt) from clinVar1_all_1
    output:
        set val(dataset), file(pgxClinPerSamplePerPopulation_tiff_file), file(clinSampleTable_csv_file), file(clinPopTable_csv_file) into clinVar_plot
    script:
        dataset_allLofPop_file = "${params.work_dir}/data/PGX_ONLY/LOF/${dataset}/${dataset}_allLofPop.csv"
        dataset_sample_file = "${params.work_dir}/data/samples/${dataset}.sample"
        clinSampleTable_csv_file = "${file(dataset_clinical_vcf.baseName).baseName}_clinSampleTable.csv"
        clinPopTable_csv_file = "${file(dataset_clinical_vcf.baseName).baseName}_clinPopTable.csv"
        pgxClinPerSamplePerPopulation_tiff_file = "${file(dataset_clinical_vcf.baseName).baseName}_pgxClinPerSamplePerPopulation.tiff"
        template "step14_3_clinical_variants.R"
}


'''
Step 14.4: Plot Clinical variants per population (AIBST pops, 1KG, gnomAD)
'''
clinVar1_all.into { clinVar1_all; clinVar1_all_2 }
clinVar_plot.into {clinVar_plot; clinVar_plot_1}
process clinVar_plot_pop {
    tag "clinVar_pop_plot_${dataset}"
    cpus { 2 * task.attempt }
    memory { 2.GB * task.cpus }
    publishDir "${params.work_dir}/data/PGX_ONLY/CLINICAL/${dataset}", overwrite: true, mode:'copy'
    publishDir "${params.work_dir}/REPORTS/${dataset}/PGX_ONLY/CLINICAL", overwrite: true, mode:'copy', pattern: '*tiff'
    input:
        set val(dataset), file(dataset_clinical_vcf), file(ClinPerSample_012_file), file(ClinPerSample_012_indv_file), file(ClinPerSample_012_pos_file), file(pgxClinicalrsID_txt) from clinVar1_all_2
    output:
        set val(dataset), file(tiff_file) into clinVar_plot_pop
    script:
        datasetAnnotated = group_pop_data[dataset][2]
        dataset_sample_file = "${params.work_dir}/data/samples/${dataset}.sample"
        pgxClinicalLevel1 = params.pgxClinicalLevel1
        tiff_file = "${file(dataset_clinical_vcf.baseName).baseName}_pgxClinicalEvidencePopulation.tiff"
        template "step14_4_clinical_variants_pop.R"
}


'''
Step 15.1: Rare Global Common Variants
'''
concat_dataset_pgx.into { concat_dataset_pgx; concat_dataset_pgx_4 }
process rareGlobCommon {
    tag "rareGlobCommon_${dataset}"
    cpus { 2 * task.attempt }
    memory { 2.GB * task.cpus }
    publishDir "${params.work_dir}/data/PGX_ONLY/RARE_GLOB_COM/${dataset}", overwrite: true, mode:'symlink'

    input:
        set val(dataset), file(dataset_sample), file(dataset_vcf), file(dataset_vcf_tbi) from concat_dataset_pgx_4
    output:
        set val(dataset), file(dataset_sample), file("${dataset}.rare.vcf.gz"), file("${dataset}.rare.IDs") into rareGlobCommon
    script:
        """
        # Make a vcf with all rare variants (0.05%) in the global dataset
        vcftools --gzvcf ${dataset_vcf} \
            --max-maf 0.005 --recode-INFO-all \
            --recode --stdout | \
            bgzip -c > ${dataset}.rare.vcf.gz
        # Get variant IDs for all rare variants
        zcat  ${dataset}.rare.vcf.gz | grep -v ^# | awk '{print \$3}' > ${dataset}.rare.IDs
        """

}


''''
Step 15.2: Rare Global Common Variants per population
'''
group_vcf_with_pop_sample_cha.into { group_vcf_with_pop_sample_cha; group_vcf_with_pop_sample_1 }
rareGlobCommon.into { rareGlobCommon; rareGlobCommon_1 }
rareGlobCommon_list = rareGlobCommon_1.toSortedList().val

rareGlobCommon_pop_all = group_vcf_with_pop_sample_1.flatMap{ POP, dataset, dataset_POP_vcf, dataset_POP_sample ->
    group_vcf_with_pop_list = []
    rareGlobCommon_list.each{ dataset_, dataset_sample, dataset_rare_vcf, dataset_rare_ids  ->
        if (dataset == dataset_){
            group_vcf_with_pop_list << [POP, dataset, dataset_rare_vcf, dataset_sample, dataset_rare_ids]
        }
    }
    return group_vcf_with_pop_list
}

process rareGlobCommon_pop {
    tag "rareGlobCommon_pop_${POP}"
    cpus { 2 * task.attempt }
    memory { 2.GB * task.cpus }
    publishDir "${params.work_dir}/data/PGX_ONLY/RARE_GLOB_COM/${dataset}", overwrite: true, mode:'symlink'

    input:
        set val(POP), val(dataset), file(dataset_rare_vcf), file(dataset_sample), file(dataset_rare_ids) from rareGlobCommon_pop_all
    output:
        set val(POP), val(dataset), file("${POP}.pop.diff.vcf.gz") into rareGlobCommon_pop
    script:
        """
        grep ${POP} ${dataset_sample} | cut -f1 > ${POP}.samples
        # Make a vcf for each population for variants that occur at 5% in that population
        # And are rare in global dataset
        vcftools --gzvcf ${dataset_rare_vcf} \
            --keep ${POP}.samples \
            --snps ${dataset_rare_ids} \
            --maf 0.05 \
            --recode-INFO-all \
            --recode --stdout | \
        bgzip -c > ${POP}.pop.diff.vcf.gz
        """
}


''''
Step 15.3: Rare Global Common Variants per population
'''
rareGlobCommon_pop.into { rareGlobCommon_pop; rareGlobCommon_pop_1 }
rareGlobCommon_pop_combine_all = [:]
rareGlobCommon_pop_combine_pop_all = [:]
rareGlobCommon_pop_1.toSortedList().val.each{ POP, dataset, POP_diff_vcf ->
    if ( !(dataset in rareGlobCommon_pop_combine_all.keySet()) ){
        rareGlobCommon_pop_combine_all[dataset] = []
        rareGlobCommon_pop_combine_pop_all[dataset] = []
    }
    // Put all population VCFs in list for group_pop
    rareGlobCommon_pop_combine_all[dataset] << POP_diff_vcf
    rareGlobCommon_pop_combine_pop_all[dataset] << [POP_diff_vcf, POP].join('__')
}
rareGlobCommon_pop_combine_cha = Channel.from( rareGlobCommon_pop_combine_pop_all.keySet() )

process rareGlobCommon_pop_combine {

    tag "rareGlobCommon_pop_combine_${dataset}"
    cpus { 2 * task.attempt }
    memory { 2.GB * task.cpus }
    publishDir "${params.work_dir}/data/PGX_ONLY/RARE_GLOB_COM/${dataset}", overwrite: true, mode:'symlink'

    input:
        val(dataset) from rareGlobCommon_pop_combine_cha
    output:
        set val(dataset), file(dataset_diff_snps_vcf), file(dataset_diff_snps_genotypes), file("${dataset}.population.pop.diff.snps.genotypes"), file("${dataset}.rsIndexPopDiff.txt"), file("${dataset}.genes.pop.diff.snps.genotypes"), file("${dataset}.population.gene.pop.diff.snps.genotypes") into rareGlobCommon_pop_combine_cha_all
        file("*pop.diff.snps.Ind.frq") into rareGlobCommon_pop_combine_cha_all_1
    script:
        dataset_diff_snps_vcf = "${dataset}.diff.snps.vcf.gz"
        dataset_diff_snps_genotypes = "${dataset}.diff.snps.genotypes"
        """
        # Extract SNPs that are differentiated
        data_array=\"${rareGlobCommon_pop_combine_pop_all[dataset].join(' ')}\"
        #touch ${dataset}.diff.snps.genotypes
        #touch ${dataset}.diff.snps.rsIDs
        for POP_VCF in \$data_array;do
            POP_VCF_array=(\${POP_VCF//__/ })
            vcf=\${POP_VCF_array[0]}
            pop=\${POP_VCF_array[1]}
            zcat \$vcf | awk '!/^(\$|#)/ {print \$0" ""'"\$pop"'"}' >> ${dataset}.diff.snps.genotypes
            # Get variant IDs for all differentiated variants
            zcat \$vcf | grep -v ^# | awk '{print \$3}' >> ${dataset}.diff.snps.rsIDs
            zcat \$vcf | awk '!/^(\$|#)/ {print \$3" "\$1":"\$2":"\$5}' >> ${dataset}.rsIndexPopDiff.txt
            #Extract gene name and variant type
            zcat \$vcf | ${params.snpSift} extractFields - -e "." ANN[0].GENE ANN[0].EFFECT | \
                awk 'NR>1 {print \$0}' >> ${dataset}.temp.genes.pop.diff.snps.genotypes
        done
        # Generate new vcf
        vcftools --gzvcf ${group_pop_data[dataset][0]} \
            --recode --snps ${dataset}.diff.snps.rsIDs \
            --recode-INFO-all --stdout | bgzip -c > ${dataset}.diff.snps.vcf.gz
        for POP_VCF in \$data_array;do
            POP_VCF_array=(\${POP_VCF//__/ })
            vcf=\${POP_VCF_array[0]}
            pop=\${POP_VCF_array[1]}
            grep \$pop ${group_pop_data[dataset][1]} | cut -f1 > \$pop.samples
            vcftools --gzvcf ${dataset}.diff.snps.vcf.gz \
                --keep \$pop.samples \
                --freq2 --derived \
                --out \$pop.pop.diff.snps.Ind
        done
        zcat ${dataset}.diff.snps.vcf.gz | awk '!/^(\$|#)/ {print \$3" "\$1":"\$2":"\$5}' > ${dataset}.rsIndexPopDiff.txt
        # Extract populations
        awk '{print \$NF}' ${dataset}.diff.snps.genotypes > ${dataset}.temp.population.pop.diff.snps.genotypes
        # Extract rsID only
        awk '{print \$3}' ${dataset}.diff.snps.genotypes > ${dataset}.temp.rs.pop.diff.snps.genotypes
        # Population and rsID
        paste ${dataset}.temp.population.pop.diff.snps.genotypes ${dataset}.temp.rs.pop.diff.snps.genotypes > ${dataset}.population.pop.diff.snps.genotypes
        # Extract population only
        awk '{print \$1}' ${dataset}.population.pop.diff.snps.genotypes > ${dataset}.population.only.pop.diff.snps.genotypes
        #Extract gene name and variant type
        paste ${dataset}.temp.genes.pop.diff.snps.genotypes ${dataset}.rsIndexPopDiff.txt > ${dataset}.genes.pop.diff.snps.genotypes
        paste ${dataset}.population.only.pop.diff.snps.genotypes ${dataset}.genes.pop.diff.snps.genotypes > ${dataset}.population.gene.pop.diff.snps.genotypes
        """
}


''''
Step 16.1: Consequences
'''
concat_dataset_pgx.into { concat_dataset_pgx; concat_dataset_pgx_5 }
process csq {
    tag "csq_${dataset}"
    cpus { 2 * task.attempt }
    memory { 2.GB * task.cpus }
    publishDir "${params.work_dir}/data/PGX_ONLY/CSQ/${dataset}", overwrite: true, mode:'symlink'
    input:
        set val(dataset), file(dataset_sample), file(dataset_vcf), file(dataset_vcf_tbi) from concat_dataset_pgx_5
    output:
        set val(dataset), file(dataset_sample), file("${dataset}.singletons.vcf.gz"), file("${dataset}_SO.terms.MAF.summary"), file("${dataset}.singletons.uniqueID"), file("${dataset}.nonsing-0.01MAF.uniqueID"), file("${dataset}.0.01-0.05MAF.uniqueID"), file("${dataset}.greater.0.05.uniqueID") into csq
        set val(dataset), file("GENE") into csq_gene_vcf
    script:
        """
        # Create different vcfs for different frequency classes
        # Singletons
        vcftools --gzvcf ${dataset_vcf} \
            --max-mac 1 --recode-INFO-all \
            --recode --stdout | \
            bgzip -c > ${dataset}.singletons.vcf.gz
        # Non-singletons to MAF 0.01
        vcftools --gzvcf ${dataset_vcf} \
            --recode --mac 2 --max-maf 0.01 \
            --recode-INFO-all --stdout | \
            bgzip -c > ${dataset}.nonsing-0.01MAF.vcf.gz
        # MAF 0.01-0.05
        vcftools --gzvcf ${dataset_vcf} \
            --recode --maf 0.01 --max-maf 0.05 \
            --recode-INFO-all --stdout | \
            bgzip -c > ${dataset}.0.01-0.05MAF.vcf.gz
        # MAF > 0.05
        vcftools --gzvcf ${dataset_vcf} \
            --recode --maf 0.05 \
            --recode-INFO-all --stdout | \
            bgzip -c > ${dataset}.greater.0.05MAF.vcf.gz

        # Count number of appearances of each of the Sequence Ontology (SO) terms
        # so_term_file was modified from:
        # http://www.ensembl.org/info/genome/variation/predicted_data.html
        while read SO; do
            zcat ${dataset_vcf} | grep \$SO | wc -l
        done < ${params.so_term_file} > ${dataset}_allCount
        while read SO; do
            zcat ${dataset}.singletons.vcf.gz | grep \$SO | wc -l
        done < ${params.so_term_file} > ${dataset}_singletonsCount
        while read SO; do
            zcat ${dataset}.nonsing-0.01MAF.vcf.gz | grep \$SO | wc -l
        done < ${params.so_term_file} > ${dataset}_nonsing-0.01MAFCount
        while read SO; do
            zcat ${dataset}.0.01-0.05MAF.vcf.gz | grep \$SO | wc -l
        done < ${params.so_term_file} > ${dataset}_0.01-0.05MAFCount
        while read SO; do
            zcat ${dataset}.greater.0.05MAF.vcf.gz | grep \$SO | wc -l
        done < ${params.so_term_file} > ${dataset}_greater0.05MAFCount
        echo "terms 	countAll	countSingletons	countnonSing-0.01MAF	count0.01-0.05MAF	countGreater0.05" \
            > ${dataset}_SO.terms.MAF.summary
        paste ${params.so_term_file} ${dataset}_allCount ${dataset}_singletonsCount  ${dataset}_nonsing-0.01MAFCount ${dataset}_0.01-0.05MAFCount ${dataset}_greater0.05MAFCount \
            >> ${dataset}_SO.terms.MAF.summary

        # Create UniqueID for each of the vcf frequency classes
        zcat ${dataset}.singletons.recode.vcf.gz | grep -v ^# | awk '{print \$1":"\$2":"\$5}' \
            >  ${dataset}.singletons.uniqueID
        zcat ${dataset}.nonsing-0.01MAF.recode.vcf.gz | grep -v ^# | awk '{print \$1":"\$2":"\$5}' \
            >  ${dataset}.nonsing-0.01MAF.uniqueID
        zcat ${dataset}.0.01-0.05MAF.recode.vcf.gz | grep -v ^# | awk '{print \$1":"\$2":"\$5}' \
            >  ${dataset}.0.01-0.05MAF.uniqueID
        zcat ${dataset}.greater.0.05MAF.recode.vcf.gz | grep -v ^# | awk '{print \$1":"\$2":"\$5}' \
            >  ${dataset}.greater.0.05.uniqueID

        # Create individual files for each coding SO term for all genes
        mkdir -p GENE

        while read SO; do
            echo \$SO > GENE/${dataset}_\$SO.gene.txt
        done < ${params.so_coding_terms}

        # Add heading to VIP genes gene symbol
        echo gene > GENE/${dataset}_symbol.txt
        cat ${params.final_all_pharmacogenes_details} >> GENE/${dataset}_symbol.txt

        # Analyse each of the genes in the study independently, appending counts
        while read gene;
        do
            SnpSift filter "ANN[*].GENE = '\$gene'" ${dataset_vcf} -v | bgzip -c > GENE/${dataset}_\$gene.vcf.gz
            while read term
            do
                zcat GENE/${dataset}_\$gene.vcf | grep \$term | wc -l >> GENE/${dataset}_\$term.gene.txt
            done < ${params.so_coding_terms}
        done < ${params.final_all_pharmacogenes_details}

        paste GENE/${dataset}_symbol.txt GENE/*gene* > GENE/${dataset}_summaryVariantsPerGene
        #rm GENE/*txt
        """
}


''''
Step 16.2: Plot consequences
'''
csq.into{ csq; csq_1 }
process plot_csq {
    tag "plot_csq_${dataset}"
    cpus { 2 * task.attempt }
    memory { 10.GB * task.cpus }
    publishDir "${params.work_dir}/data/PGX_ONLY/CSQ/${dataset}", overwrite: true, mode:'symlink'
    publishDir "${params.work_dir}/REPORTS/${dataset}/PGX_ONLY/CSQ", overwrite: true, mode:'symlink'
    input:
        set val(dataset), file(dataset_sample), file(dataset_singletons_vcf_gz), file(dataset_SO_terms_MAF_summary), file(dataset_singletons_uniqueID), file(dataset_nonsing_0_01MAF_uniqueID), file(dataset_0_01_0_05MAF_uniqueID), file(dataset_greater_0_05_uniqueID) from csq_1
    output:
        set val(dataset), file(dataset_pgxFunctionalClassesCounts_tiff), file(dataset_pgxFunctionalClassesByFrequency_tiff) into plot_csq
    script:
        dataset_SO_terms_MAF_summary = dataset_SO_terms_MAF_summary
        dataset_pgxFunctionalClassesCounts_tiff = "${dataset}_pgxFunctionalClassesCounts.tiff"
        dataset_pgxFunctionalClassesByFrequency_tiff = "${dataset}_pgxFunctionalClassesByFrequency.tiff"
        template "step16_2_consequence.R"
}


''''
Step 17.1: RareMissense
'''
csq_gene_vcf.into { csq_gene_vcf; csq_gene_vcf_1 }
process RareMissense {
    tag "RareMissense_${dataset}"
    cpus { 2 * task.attempt }
    memory { 2.GB * task.cpus }
    publishDir "${params.work_dir}/data/PGX_ONLY/RareMissense/${dataset}", overwrite: true, mode:'symlink'
    input:
        set val(dataset), file(gene_vcfs_url) from csq_gene_vcf_1
    output:
        set val(dataset), file("${dataset}_rare.missense.gene.txt") into RareMissense
    script:
        """
        echo "rare_missense" > ${gene_vcfs_url}/${dataset}_temp.rare.missense.gene.txt
        # Analyse each of the genes in the study independently, appending counts
        while read gene;
        do
            vcftools --gzvcf ${gene_vcfs_url}/${dataset}_\$gene.vcf.gz --max-maf 0.005 --recode-INFO-all --recode --stdout | \
            grep missense_variant | wc -l >> ${gene_vcfs_url}/${dataset}_temp.rare.missense.gene.txt
        done < ${params.final_all_pharmacogenes_details}
        paste ${gene_vcfs_url}/${dataset}_symbol.txt ${gene_vcfs_url}/${dataset}_temp.rare.missense.gene.txt > ${dataset}_rare.missense.gene.txt
        #rm ${gene_vcfs_url}/${dataset}_temp.rare.missense.gene.txt
        """
}


//
// ''''
// Step 18.1: VIP Mask
// '''
// concat_dataset_all.into{ concat_dataset_all; concat_dataset_all_18}
// process vipmask{
//     tag "vipmask_${dataset}"
//     cpus { 2 * task.attempt }
//     memory { 2.GB * task.cpus }
//     publishDir "${params.work_dir}/data/PGX_ONLY/VIPMASK/${dataset}", overwrite: true, mode:'symlink'
//     input:
//         set val(dataset), file(dataset_sample), file(dataset_vcf), file(dataset_vcf_tbi) from concat_dataset_all_18
//     output:
//         set val(dataset), file("${dataset_out}_maskGeneSummary.txt") into vipmask
//     script:
//         dataset_out = "${file(dataset_vcf.baseName).baseName}"
//         """
//         ###sed s/^chr//g 20141020.strict_mask.whole_genome.bed > 20141020.strict_mask.whole_genome_final.bed
//
//         # Make vcf of masked/unmasked regions
//         vcftools --gzvcf ${dataset_vcf} \
//            --bed ${params.strict_mask_whole_genome_final} \
//            --recode --recode-INFO-all --stdout \
//            | bgzip -c > ${dataset_out}_accessible.vcf.gz
//         vcftools --gzvcf ${dataset_vcf} \
//            --exclude-bed ${params.strict_mask_whole_genome_final} \
//            --recode --recode-INFO-all --stdout \
//            | bgzip -c > ${dataset_out}_inaccessible.vcf.gz
//
//         # Look at segmental duplications
//         # Make two vcf files, one within segmental dups, another outside of them
//         vcftools --gzvcf ${dataset_vcf} \
//            --bed ${params.GRCh37GenomicSuperDup} \
//            --recode --recode-INFO-all --stdout \
//            | bgzip -c > ${dataset_out}_segmental.vcf.gz
//         vcftools --gzvcf ${dataset_vcf} \
//            --exclude-bed ${params.GRCh37GenomicSuperDup} \
//            --recode --recode-INFO-all --stdout \
//            | bgzip -c > ${dataset_out}_nonsegmental.vcf.gz
//
//
//         # Make a vcf of SNPs found in both segmental duplications and inaccessible regions
//         vcftools --gzvcf ${dataset_out}_segmental.vcf.gz \
//             --exclude-bed ${params.strict_mask_whole_genome_final} \
//             --stdout --recode --recode-INFO-all \
//             | bgzip -c > ${dataset_out}.segmental.inaccessible.vcf.gz
//
//
//         # Calculate the number of variants per gene
//         echo total_variants > ${dataset_out}_totalVariantsPerGene.txt
//         while read gene;
//         do
//             zcat ${params.work_dir}/data/PGX_ONLY/CSQ/${dataset}/GENE/${dataset}_\$gene.vcf.gz \
//                 | grep -v ^#  | wc -l >> ${dataset_out}_totalVariantsPerGene.txt
//         done < ${params.final_all_pharmacogenes_details}
//
//
//         # Look at number of variants in inaccessible regions
//         # Calculate the number of variants per gene
//         echo inaccessible_variants > ${dataset_out}_inaccessiblePerGene.txt
//         while read geneIn;
//         do
//             zcat  ${dataset_out}_inaccessible.vcf.gz \
//                 | grep -v ^# | grep -w \$geneIn | wc -l >> ${dataset_out}_inaccessiblePerGene.txt
//         done < ${params.final_all_pharmacogenes_details}
//
//
//         # Look at number of variants in segmental duplications
//         # Calculate the number of variants per gene
//         echo segmental_variants > ${dataset_out}_segmentalPerGene.txt
//         while read geneSeg;
//         do
//             zcat ${dataset_out}_segmental.vcf.gz \
//                 | grep -v ^# | grep -w \$geneSeg | wc -l >> ${dataset_out}_segmentalPerGene.txt
//         done < ${params.final_all_pharmacogenes_details}
//
//
//         # Look at number of variants in inaccessible and segmental duplications
//         # Calculate the number of variants per gene
//         echo inaccess_segmental_variants > ${dataset_out}_inaccessSegmentalPerGene.txt
//         while read geneBoth;
//         do
//             zcat ${dataset_out}_inaccessible.vcf.gz \
//                 |grep -v ^# | grep -w \$geneBoth | wc -l >> ${dataset_out}_inaccessSegmentalPerGene.txt
//         done < ${params.final_all_pharmacogenes_details}
//
//
//         # Add heading to VIP genes gene symbol
//         echo gene > ${dataset_out}_symbol.txt
//         awk '{print \$1}' ${params.final_all_pharmacogenes_details} >> ${dataset_out}_symbol.txt
//
//         # Calculate total length of each gene that was extracted
//         echo geneBedLengthTotal > ${dataset_out}_geneBedLengthTotal.txt
//         while read geneID
//         do
//             grep -w \$geneID ${params.gene_regions_slop25} | \
//             sort -k1,1 -k2,2n | bedtools merge -i stdin | \
//             awk -F'\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}' \
//             >> ${dataset_out}_geneBedLengthTotal.txt
//         done < ${params.final_all_pharmacogenes_details}
//
//         # Calculate transcript length of each gene that was extracted
//         echo geneBedLengthTranscript > ${dataset_out}_geneBedLengthTranscript.txt
//         while read geneID
//         do
//             grep -w \$geneID ${params.gene_regions} | \
//             sort -k1,1 -k2,2n | bedtools merge -i stdin | \
//             awk -F'\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}' \
//             >> ${dataset_out}_geneBedLengthTranscript.txt
//         done < ${params.final_all_pharmacogenes_details}
//
//         # Combine all data together
//         paste ${dataset_out}_symbol.txt ${dataset_out}_totalVariantsPerGene.txt ${dataset_out}_inaccessiblePerGene.txt ${dataset_out}_segmentalPerGene.txt ${dataset_out}_inaccessSegmentalPerGene.txt ${dataset_out}_geneBedLengthTranscript.txt ${dataset_out}_geneBedLengthTotal.txt > ${dataset_out}_maskGeneSummary.txt
//
//         """
// }









//////'''
//////Merge genes.pop.diff.snps.genotypes and population.pop.diff.snps.genotypes
//////'''
//////rareGlobCommon2_all.into{rareGlobCommon2_all; rareGlobCommon2_all_1}
//////process rareGlobCommon2_1 {
//////    echo true
//////    tag "rareGlobCommon2_1_${GROUP_POP}"
//////    cpus { 2 * task.attempt }
//////    memory { 2.GB * task.cpus }
//////    publishDir "${params.work_dir}/data/MERGED_POP/POPDIFF/${GROUP_POP}", overwrite: true, mode:'symlink'
//////
//////    input:
//////        set val(GROUP_POP), file(GROUP_POP_diff_snps_vcf), val(GROUP_POP_diff_snps_genotypes), val(GROUP_POP_population_pop_diff_snps_genotypes), file(GROUP_POP_rsIndexPopDiff_txt), val(GROUP_POP_genes_pop_diff_snps_genotypes) from rareGlobCommon2_all_1
//////    output:
//////        set val(GROUP_POP), file("${GROUP_POP}.population.gene.pop.diff.snps.genotypes") into rareGlobCommon2_1
//////    script:
//////        pop_diff = GROUP_POP_population_pop_diff_snps_genotypes.text
//////        println file(pop_diff).readLines()
//////}
////
////
////''''
////Step 8.3: Compute allele frequencies for each population for rare/common sites
////'''
////rareGlobCommon2_all.into { rareGlobCommon2_all; rareGlobCommon2_pop }
////rareGlobCommon2_pop_list = [] // For POP, GROUP_POP,  GROUP_POP_sample, GROUP_diff_snps_vcf
////rareGlobCommon2_pop.toSortedList().val.each { GROUP_POP, GROUP_diff_snps_vcf, GROUP_POP_genotype, GROUP_POP_population_genotypes, GROUP_POP_rsIndexPopDiff, GROUP_POP_genes_genotypes, GROUP_POP_pop_only_genotypes ->
////    data = file(group_pop_data[GROUP_POP][1]).readLines().collect { it.split()[1]}.unique()
////    data.each { POP ->
////        rareGlobCommon2_pop_list << [ POP, GROUP_POP, group_pop_data[GROUP_POP][1], GROUP_diff_snps_vcf ]
////    }
////}
////rareGlobCommon2_pop_chan = Channel.from( rareGlobCommon2_pop_list )
////
////process rareGlobCommon3 {
////    echo true
////    tag "rareGlobCommon3_${POP}"
////    cpus { 2 * task.attempt }
////    memory { 2.GB * task.cpus }
////    publishDir "${params.work_dir}/data/MERGED_POP/POPDIFF/${GROUP_POP}", overwrite: true, mode:'symlink'
////
////    input:
////        set val(POP), val(GROUP_POP), file(GROUP_POP_sample), file(GROUP_diff_snps_vcf) from rareGlobCommon2_pop_chan
////    output:
////        set val(POP), val(GROUP_POP), file("${POP}.diff.snps.ind.frq") into rareGlobCommon3_all
////    script:
////        """
////        grep ${POP} ${GROUP_POP_sample} | cut -f1 > ${POP}.samples
////        # Extract sites for population
////        # Generate sites frequency for population
////        vcftools --gzvcf ${GROUP_diff_snps_vcf} \
////            --keep ${POP}.samples \
////            --recode-INFO-all \
////            --recode --stdout | \
////        vcftools --vcf - \
////            --freq2 \
////            --out ${POP}.diff.snps.ind
////        """
////}
////rareGlobCommon3_all.into{ rareGlobCommon3_all; rareGlobCommon3__sub }
////rareGlobCommon3__sub.subscribe {
////      //println "|-- Finished for ${it[-1]}"
////}
////
////''''
////Step 8.4: Plot allele frequencies for each population for rare/common sites
////'''
////rareGlobCommon2_all.into{rareGlobCommon2_all; rareGlobCommon2_4}
////process rareGlobCommon4 {
////    echo true
////    tag "rareGlobCommon4_${GROUP_POP}"
////    cpus { 2 * task.attempt }
////    memory { 2.GB * task.cpus }
////    publishDir "${params.work_dir}/data/MERGED_POP/POPDIFF/${GROUP_POP}", overwrite: true, mode:'symlink'
////
////    input:
////        set val(GROUP_POP), file(GROUP_POP_diff_snps_vcf_gz), file(GROUP_POP_diff_snps_genotypes), file(GROUP_POP_population_pop_diff_snps_genotypes), file(GROUP_POP_rsIndexPopDiff_txt), file(GROUP_POP_genes_pop_diff_snps_genotypes), file(GROUP_POP_population_gene_pop_diff_snps_genotypes) from rareGlobCommon2_4
////    output:
////        set val(GROUP_POP), file("${GROUP_POP}.highDiffPlot.tiff") into rareGlobCommon4_all
////    script:
////        cadd_annotations = file(params.cadd_annotations)
////        GROUP_POP_sample = file(group_pop_data[GROUP_POP][1])
////        GROUP_POP_datasetAnnotated = file(group_pop_data[GROUP_POP][2])
////        GROUP_POP_tiff = "${GROUP_POP}.highDiffPlot.tiff"
////        GROUP_POP_csv = "${GROUP_POP}.highDiffPlot.csv"
////        println GROUP_POP_diff_snps_genotypes.readLines()
//////        if (GROUP_POP_diff_snps_genotypes.size() == 0){
//////            println GROUP_POP_diff_snps_genotypes
//////        }
////        template "step8_populationFreq.R"
////}
//
//

//////workflow.onComplete {
//////    def subject = 'My pipeline execution'
//////    def recipient = 'mypandos@gmail.com'
//////
//////    ['mail', '-s', subject, recipient].execute() << """
//////
//////    Pipeline execution summary
//////    ---------------------------
//////    Completed at: ${workflow.complete}
//////    Duration    : ${workflow.duration}
//////    Success     : ${workflow.success}
//////    workDir     : ${workflow.workDir}
//////    exit status : ${workflow.exitStatus}
//////    Error report: ${workflow.errorReport ?: '-'}
//////    """
//////}
//
///* TODO
//masterChannel = Channel.from(params.files).println()
//[condition1, [condition1.bam, input.bam]]
//[condition2, [condition2.bam, input.bam]]
//
//the map must return the new object, hence
//Channel
//        .from(params.files)
//        .map { condition, list ->
//    def files = list.collect{ file(it)}
//    return tuple(condition, files)
//}
//*/
