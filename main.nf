#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { dataset_qc_alt; get_vcf_site ; tabix_data_proc_pop; sites_only } from './module/vcf_qc'
include { fst_analysis; combine_fst_weir_estimates; combine_fst_weir_estimates as combine_fst_weir_estimates_cutoff; combine_fst_mean_weir_estimates; combine_fst_mean_weir_estimates as combine_fst_mean_weir_estimates_1; combine_fst_mean_weir_estimates as combine_fst_mean_weir_estimates_2; get_fst_weir_estimates; get_fst_weir_estimates_cutoff; combine_weir_fst_analysis; plot_fst_matrix ; generate_fst_matrix } from './module/fst'
include { check_files; singleton_pop; split_dataset_vcf_pop; split_dataset_vcf_pop_tabix; concat_pop; group_pops_data; merge_pop_groups; cat_chrm_groups; combine_csv } from './module/exome_analysis'
include { annotate_snpeff; annotate_dbsnp; annotate_clinvar; annotate_cosmic; combine_mafs_dataset; fix_mafs_headers; annotate_mafs; fill_tags; get_pop; concat_mafs; concat_vcf; annotate_vcf_from_vcf; concat_chunks_vcf } from './module/annotate_vcf.nf'
include { get_sites as get_sites1; get_sites as get_sites2; get_gene_vcf; get_gene_vcf1 as get_novel_vcf; get_gene_vcf as get_pgx_pop; get_gene_vcf_simple; get_gene_vcf1; vcf_to_plink; generate_pheno_fam_sample; generate_pheno_fam_sample1;  get_protein_coding_variants; get_vcf_ind; merge_groups; generate_chunks_vcf; split_vcf_to_chunk;} from './module/subset_vcf.nf'
include { novel_in_existing_database; get_novel_sites; novel_count; novel_count_combine ;combine_novel_count as combine_novel_count3;combine_novel_count as combine_novel_count5;combine_novel_count as combine_novel_count4;  combine_novel_count as combine_novel_count1; combine_novel_count as combine_novel_count2; count_variants; combine_counts; generate_circos_dataframe; generate_circos_plot; singleton_dataset_chrm} from './module/novel_vcf'
include { get_freq; get_freq2; get_freq_3; get_freq_3_1; cat_freq; cat_sites; data_for_upset ; plot_upset_R; fisher_test_plink; fisher_test_plink_pops; filter_fisher_test; split_id_freq ; split_id_freq_1; split_id_freq_2; get_sites_simple; get_map_vcf; get_map; get_freq_3_2; get_freq_3_3; get_freq_from_vcf_sites; add_clinvar_to_freq } from './module/freq_vcf'
include { hdv; hdv_dataset; combine_hdvs; combine_hdvs as combine_hdvs_all; combine_hdvs as combine_hdvs_base; combine_hdvs as combine_hdvs_base_all ; annotate_hdvs; } from './module/hdv_vcf'
include { csq; plot_csq  } from './module/csq'
include { pruning_vcf; vcf_to_plink1; smartpca_dataset; plot_pca_group; update_evec; updated_ped } from './module/structure'

workflow data_proc {
    main:
        // println "Project : $workflow.projectDir"
        // println "Cmd line: $workflow.commandLine"
        // println "Chromosomes used: ${params.chromosomes.join(',')}"

        // // Check all files
        check_files([params.dbsnp_vcf, "${params.dbsnp_vcf}.tbi", params.clinvar, "${params.clinvar}.tbi", params.cosmic, "${params.cosmic}.tbi"])

        dataset_pop_files = []
        dataset_chrm_files = []
        // params.dataset_files.each{ dataset_name, dataset_vcf, dataset_sample, dataset_pops ->
        //     params.chromosomes.each { chrm ->
        //         dataset_vcf_chrm = file(sprintf(dataset_vcf, chrm)) //TODO use string replacement https://www.nextflow.io/docs/latest/script.html#string-replacement
        //         check_files([dataset_vcf_chrm, dataset_sample])
        //         // only on AIBST
        //         // if(dataset_name == 'AIBST'){
        //             dataset_chrm_files << [dataset_name, file(dataset_vcf_chrm), file(dataset_sample), chrm]
        //         // }
        //     }
        // }

        datasets_map = []
        params.dataset_files.each { dataset_name, dataset_vcf, dataset_sample, dataset_pops ->
            vcfs = file(dataset_vcf)
            if (vcfs instanceof List){
                vcfs.each{ vcf_ ->
                    if( file(vcf_).getExtension() == "gz" ){
                        vcf_idx = "${vcf_}.tbi"
                    }
                    if( file(vcf_).getExtension() == "bcf" ){
                        vcf_idx = "${vcf_}.csi"
                    }
                    check_files([ vcf_, vcf_idx, dataset_sample ])
                    datasets_map << [ dataset_name, file(vcf_), file(vcf_idx), file(dataset_sample) ]
                }
            }
            else{
                if( file(dataset_vcf).getExtension() == "gz" ){
                    vcf_idx = "${dataset_vcf}.tbi"
                }
                if( file(dataset_vcf).getExtension() == "bcf" ){
                    vcf_idx = "${dataset_vcf}.csi"
                }
                check_files([ dataset_vcf, vcf_idx, dataset_sample ])
                datasets_map << [ dataset_name, file(dataset_vcf), file(vcf_idx), file(dataset_sample) ]
            }
        }

        datasets_map_data = Channel.from(datasets_map)

        get_map(datasets_map_data)
        generate_chunks_vcf(get_map.out.map{ dataset, vcf, vcf_idx, sample_file, map_file -> [ dataset, file(vcf), file(vcf_idx), file(sample_file), file(map_file), '', params.chunk_size ] })
        chunks_datas = generate_chunks_vcf.out.flatMap{ dataset, vcf, vcf_idx, sample_file, chunk_file ->
            datas = []
            chunks = file(chunk_file).text.split()
            chunks.each{ chunk_data ->
                data = chunk_data.split(',')
                chrm = data[0]
                start = data[1]
                end = data[2]
                datas << [dataset, chrm, start, end, file(vcf), file(vcf_idx), file(sample_file)]
            }
            return datas
        }
        split_vcf_to_chunk(chunks_datas)

        // // Data QC
        dataset_qc_alt(split_vcf_to_chunk.out)

        // // Annotate
        //////annotate_snpeff(annotate_mafs.out.combine([[params.snpeff_human_db, params.snpeff_database]]))
        check_files([params.dbsnp_vcf, "${params.dbsnp_vcf}.tbi", params.clinvar, "${params.clinvar}.tbi", params.cosmic, "${params.cosmic}.tbi"])
        annotate_snpeff(dataset_qc_alt.out.combine([[params.snpeff_human_db, params.snpeff_database]]))
        annotate_dbsnp(annotate_snpeff.out.combine([[file(params.dbsnp_vcf), file("${params.dbsnp_vcf}.tbi")]]))
        annotate_clinvar(annotate_dbsnp.out.combine([[file(params.clinvar), file("${params.clinvar}.tbi")]]))
        annotate_cosmic(annotate_clinvar.out.combine([[file(params.cosmic), file("${params.cosmic}.tbi")]]))
        process_mafs_annot()
        vcf_annot_data_ = []
        vcf_annot_data = annotate_cosmic.out.combine(process_mafs_annot.out.mafs, by:3)
        annotate_vcf_from_vcf(vcf_annot_data)

        // get_vcf_site(annotate_cosmic.out)
        // // Step 1: Annotate
        // // vcfs_mafs = get_vcf_site.out.map{ it1 -> combine_mafs(it1, params.mafs_annotations)}
        // vcfs_mafs = get_vcf_site.out.map{ it -> 
        //     combine_mafs_datas = []
        //     maf_datasets = []
        //     maf_files = []
        //     dataset_name = it[0]
        //     dataset_vcf = it[1]
        //     dataset_sample = it[2]
        //     dataset_sites = it[3]
        //     chrm = it[4]
        //     params.mafs_annotations.each{ maf_dataset, maf_file ->
        //         maf_file = file(sprintf(maf_file, chrm))
        //         if(!file(maf_file).exists()) exit 1, "File ${maf_file} not found. Please check your config file."
        //         if( maf_dataset.size() >= 6 ) {
        //             maf_dataset = maf_dataset[0..5]
        //         }
        //         maf_datasets << maf_dataset
        //         maf_files << maf_file
        //     }
        //     combine_mafs_datas = [dataset_name, file(dataset_vcf), file(dataset_sample), file(dataset_sites), chrm, maf_datasets.join('-'), maf_files.join(';')]
        //     return combine_mafs_datas
        //     }
        // combine_mafs_dataset(vcfs_mafs)
        // annotate_mafs(combine_mafs_dataset.out)

        // VCF all variants annotated
        fill_tags(annotate_vcf_from_vcf.out)

        // Concat chunks into chromosome
        concat_chunks_vcf(fill_tags.out.groupTuple(by:[0,3]))

        singleton_dataset_chrm(concat_chunks_vcf.out.map{ dataset, vcf, sample, chrm -> [ dataset, file(vcf), chrm]})

        // Get Protein Coding regions
        // get_protein_coding_variants(fill_tags.out[0])

        // Get PGX variants dataset 
        //// get_gene_vcf(fill_tags.out[0].combine(Channel.from(params.gene_lists).map{ name, bed -> [name, file(bed)]}))
        get_gene_vcf(concat_chunks_vcf.out.combine(Channel.from(params.gene_lists).map{ name, bed -> [name, file(bed)]}))
        concat_vcf(get_gene_vcf.out.groupTuple())
        sites_only(concat_vcf.out)
        // concat_vcf.out.map{ dataset, vcf, sample -> [ 'ALL', dataset, file(vcf), file(sample) ]}.groupTuple().map{ group, datasets, vcfs, samples -> [ group, vcfs]}
        get_freq_3_2(concat_vcf.out.map{ dataset, vcf, sample -> [ dataset, file(vcf), file(sample), '' ]})
        split_id_freq_1(get_freq_3_2.out)
        
        // Group datasets
        merge_datasets_group(sites_only.out)
        // VCF for each sample for all, coding regions and Pgx
        // get_sample_vcfs(fill_tags.out[0]) 

    emit:
        annotate_cosmic
        qc_dataset_chrm = dataset_qc_alt.out // before annotations
        // fill_tags_dataset_chrm = fill_tags.out[0] // With sample file
        fill_tags_dataset_chrm = concat_chunks_vcf.out // With sample file
        singl_dataset_chrm = singleton_dataset_chrm.out
        //// fill_tags_dataset = fill_tags.out[1] // Without sample file
        // exome_dataset_chrm = get_protein_coding_variants.out
        pgx_dataset_chrm = get_gene_vcf.out
        pgx_dataset = concat_vcf.out
        all_datasets_annotated = merge_datasets_group.out.all_dataset_annotation
        // Per sample VCF
        // samples_vcf = get_sample_vcfs.out.samples_vcf
        // samples_pgx_vcf = get_sample_vcfs.out.samples_pgx_vcf
        // samples_exome_vcf = get_sample_vcfs.out.samples_exome_vcf
        // datasets_chrm
        datasets_map_data
}


def get_dataset_pop(){
    add_pop_groups_datas = []
    params.dataset_files.each{ dataset_name, dataset_vcf_, dataset_sample_, dataset_pops ->
        add_pop_groups_datas << [dataset_name, dataset_name]
        dataset_pops.split(',').each{ pop ->
            add_pop_groups_datas << [pop, dataset_name]
        }
    }
    return add_pop_groups_datas
}

workflow merge_datasets_group{
    take: data
    main:
        merge_groups(data.map{ dataset, vcf, sample -> [ 'ALL', file(vcf) ]}.groupTuple())
        get_freq_from_vcf_sites(merge_groups.out.map{ dataset, vcf -> [ dataset, file(vcf), '', '' ]})
        check_files([ params.clinvar1 ])
        add_clinvar_to_freq(get_freq_from_vcf_sites.out.map{ dataset, vcf_file, sample_file, chrm, freq_file -> [ dataset, freq_file, file(params.clinvar1), file(params.pharmgkb_file), file(params.gwas_file) ] })
    emit:
        all_dataset_annotation = add_clinvar_to_freq.out
        // all_dataset_annotation = get_freq_from_vcf_sites.out
}
workflow data_proc_pop {
    take: data
    main:
        /// PGx variants 
        gene_lists = Channel.from(params.gene_lists)
            .map{ name, bed -> [name, file(bed)]}
        /// Per population
        dataset_dataset_pops = data
            .flatMap{ dataset, dataset_vcf, dataset_sample, chrm ->
                add_pop_groups_datas = []
                params.dataset_files.each{ dataset_name, dataset_vcf_, dataset_sample_, dataset_pops ->
                    dataset_pops.split(',').each{ pop ->
                        if(dataset == dataset_name){
                            add_pop_groups_datas << [pop, dataset, file(dataset_vcf), file(dataset_sample), chrm]
                        }
                    }
                }
                return add_pop_groups_datas
            }
        split_dataset_vcf_pop(dataset_dataset_pops)
        pop_data = split_dataset_vcf_pop.out[0]
            .combine(gene_lists)
            .map{ pop, dataset, pop_vcf, pop_sample, chrm, pgx_name, pgx_bed -> [ pop, file(pop_vcf), file(pop_sample), chrm, pgx_name, file(pgx_bed) ] } //TODO add dataset here

        get_pgx_pop(pop_data)
        pop_pgx_data_chrm = get_pgx_pop.out.combine(get_dataset_pop(), by:0).map{ pop, pop_vcf, pop_sample, chrm, dataset -> [ pop, dataset, pop_vcf, pop_sample, chrm ] }

        /// Concatenate vcf by population 
        concat_pop(pop_pgx_data_chrm.groupTuple())  
        /// Get PGX variants by population 
        pop_pgx_data = concat_pop.out
            .map{ pop, dataset, pop_vcf, pop_sample -> [ pop, file(pop_vcf), file(pop_sample), dataset ] }

    emit:
        pop_exome_chrm = split_dataset_vcf_pop.out
        pop_exome = concat_pop.out
        pop_pgx = pop_pgx_data
        pop_pgx_chrm = pop_pgx_data_chrm
}

workflow process_mafs_annot {
    main:
        datas = []
        params.mafs_annotations.each { dataset, vcf ->
            vcfs = file(vcf)
            if (vcfs instanceof List){
                vcfs.each{ vcf_ ->
                    if( file(vcf_).getExtension() == "gz" ){
                        vcf_idx = "${vcf_}.tbi"
                    }
                    if( file(vcf_).getExtension() == "bcf" ){
                        vcf_idx = "${vcf_}.csi"
                    }
                    check_files([ vcf_, vcf_idx ])
                    datas << [ dataset, file(vcf_), file(vcf_idx) ]
                }
            }
            else{
                if( file(vcf).getExtension() == "gz" ){
                    vcf_idx = "${vcf}.tbi"
                }
                if( file(vcf).getExtension() == "bcf" ){
                    vcf_idx = "${vcf}.csi"
                }
                check_files([ vcf, vcf_idx ])
                datas << [ dataset, file(vcf), file(vcf_idx) ]
            }
        }
        mafs_ = Channel.from(datas)
        get_vcf_chrom(mafs_)
    emit:
        mafs = get_vcf_chrom.out.vcf_chrm
}

workflow get_vcf_chrom {
    take:
        data
    main:
        get_map_vcf( data )
        vcf_chrm = get_map_vcf.out
            .flatMap{ dataset, vcf, vcf_idx, dataset_map ->
                chrm_annot_data = []
                chrms = file(dataset_map).text.split()
                chrms.each{ chrm ->
                    chrm_annot_data << [ dataset, vcf, vcf_idx, chrm ] //chrm.toInteger() 
                }
                return chrm_annot_data
            }
    emit:
        vcf_chrm
}

workflow pca{
    take: 
        data
    main:
        pruning_vcf(data.map{ group_3, chrm_3, group_vcf_3, group_sample_3 -> [ group_3, file(group_vcf_3), file(group_sample_3) ] })
        groups_data = pruning_vcf.out.groupTuple(by:[0])
        cat_chrm_groups(groups_data)
        vcf_to_plink1(cat_chrm_groups.out)
        updated_ped(vcf_to_plink1.out)
        smartpca_dataset(updated_ped.out)
        update_evec(smartpca_dataset.out)
        plot_pca_group(update_evec.out)
    emit:
        data
}

workflow pca_pgx{
    take: 
        data
    main:
        pruning_vcf(data)
        vcf_to_plink1(pruning_vcf.out)
        updated_ped(vcf_to_plink1.out)
        smartpca_dataset(updated_ped.out)
        update_evec(smartpca_dataset.out)
        plot_pca_group(update_evec.out)
    emit:
        data
}

workflow hdv_all_chrm_dataset {
    take: data
    main:
        // Get allele frequency
        get_freq(data)
        split_id_freq(get_freq.out)
        cat_freq(split_id_freq.out.map{ dataset, vcf, sample, frq1, frq2, frq3, chrm -> [ dataset, file(frq2), chrm] }.groupTuple())
        // Get HDV - whole genome
        // hdv_dataset(cat_freq.out)
        
    emit:
        data
        // hdvs = hdv_dataset.out
}

workflow hdv_pgx_dataset {
    take: data
    main:
        /* Get allele frequency */
        get_freq_3_1(data.map{ def dataset, def vcf, def sample, def chrm, def group, def pvalue_file -> [ dataset, file(vcf), file(sample), chrm, file(pvalue_file)] })
        split_id_freq(get_freq_3_1.out)
        /* HDV pgx variants */
        fold = 2
        ac = 3
        group_treshold = 1
        test_pops = 'AGVP_AF,KG_AF,KG_AFR_AF,gnomAD_AF,gnomAD_AFR_AF,TOPMED_AF'
        hdv_dataset(split_id_freq.out.map{ def dataset, def freq_file, def pvalue_file, def chrm -> [ dataset, file(freq_file), file(pvalue_file), chrm, fold, test_pops, ac, group_treshold, params.fst_cutoff ] })
        combine_hdvs(hdv_dataset.out.groupTuple(by:[0]).map{ dataset, hdvs -> [ dataset, hdvs, "dataset_${fold}fold_ac${ac}_g${group_treshold}" ]})
    emit:
        get_freq_3_1
        hdvs = hdv_dataset.out
}

workflow combine_hdvs_pops {
    take: data
    main:
        // data.groupTuple(by:[0]).map{ pop, hdvs -> [ pop, hdvs.join(' ') ] }
        // fold = 2
        // ac = 3
        // group_treshold = 1

        // Combine all
        combine_hdvs(
            data.map{ pop, hdv_all, hdv_base, fold, ac, group_treshold -> 
                if(file(hdv_all).countLines() > 1 && ac == 3 && fold == 2){
                    return [ pop, hdv_all, fold, ac, group_treshold ]}
            }.groupTuple(by:[0,2,3,4]).map{ pop, hdvs_all, fold, ac, group_treshold -> [ pop, hdvs_all, fold, ac, group_treshold, "pop_${fold}fold_ac${ac}_g${group_treshold}_all" ]})
        combine_hdvs_all(
            combine_hdvs.out.map{ pop, hdv, fold, ac, group_treshold -> [ 'ALL', pop, hdv, fold, ac, group_treshold ] }.groupTuple(by:[0,3,4,5]).map{ dataset, pops, hdvs, fold, ac, group_treshold  -> [ dataset, hdvs, fold, ac, group_treshold, "All_hdvs_pops_${fold}fold_ac${ac}_g${group_treshold}_all" ]})

        // Combine base 
        combine_hdvs_base(
            data.map{ pop, hdv_all, hdv_base, fold, ac, group_treshold -> 
                if(file(hdv_base).countLines() > 1 && ac == 3 && fold == 2){
                    return [ pop, hdv_base, fold, ac, group_treshold ]}
            }
            .groupTuple(by:[0,2,3,4]).map{ pop, hdvs_base, fold, ac, group_treshold -> [ pop, hdvs_base, fold, ac, group_treshold, "pop_${fold}fold_ac${ac}_g${group_treshold}_base" ]})
        combine_hdvs_base_all(
            combine_hdvs_base.out.map{ pop, hdv, fold, ac, group_treshold -> [ 'ALL', pop, hdv, fold, ac, group_treshold ] }.groupTuple(by:[0,3,4,5]).map{ dataset, pops, hdvs, fold, ac, group_treshold -> [ dataset, hdvs, fold, ac, group_treshold, "All_hdvs_pops_${fold}fold_ac${ac}_g${group_treshold}_base" ]})
        hdvs_pops = combine_hdvs.out.concat( combine_hdvs_all.out, combine_hdvs_base.out, combine_hdvs_base_all.out)
    emit:
        // data
        hdvs_pops
}

workflow annotate_hdvs_all {
    take: data
    main:
        annotate_hdvs(data.map{ dataset, dataset_hdv, fold, ac, group_treshold, group, group_annot -> [ dataset, file(dataset_hdv), file(group_annot) ] })

    emit:
        annotate_hdvs.out
}

workflow hdv_pgx_pop {
    take: data
    main:
        /* Get allele frequency */
        get_freq_3_3(data)
        split_id_freq_2(get_freq_3_3.out)
        /* HDV pgx variants */
        fold = 2
        ac = 3
        group_treshold = 1
        fold1 = Channel.of(2)
        ac1 = Channel.of(3)
        // group_treshold1 = Channel.of(1,2,3,4,5,6,7,8,9,10)
        group_treshold1 = Channel.of(1,2,3,4,5)
        fold_ac_group_treshold = fold1.combine(ac1).combine(group_treshold1)
        base_pops = 'KNK,KNL,NGY,SAV,TZA,TZB,ZWD,ZWS'
        // base_pops = params.FST_TEST
        test_pops = params.FST_TEST
        // test_pops = 'KNK,KNL,KNM,KNP,NGH,NGI,NGY,SAV,TZA,TZB,ZWD,ZWS,ACB,ASW,LWK,ESN,GWD,MSL,YRI,Zulu'
        hdv_pop_data = split_id_freq_2.out.combine(fold_ac_group_treshold)
            .map{ dataset, freq_file, chrm, pvalue_file, fst_file, fold, ac, group_treshold ->  [ dataset, file(freq_file), file(pvalue_file), file(fst_file), chrm, fold, test_pops, base_pops, ac, group_treshold, params.fst_cutoff ]}
        hdv_dataset(hdv_pop_data)
        // hdv_dataset.out.map{ pop, hdv_all, hdv_base, fold, ac, group_treshold -> 
        //         if(file(hdv_all).countLines() > 1 && fold == 2 && ac == 3 && group_treshold == 2){
        //             return [ pop, hdv_all, fold, ac, group_treshold ]}
        //     }
        //     .groupTuple(by:[0,2,3,4]).map{ pop, hdvs_all, fold, ac, group_treshold -> [ pop, hdvs_all, fold, ac, group_treshold, "pop_${fold}fold_ac${ac}_g${group_treshold}_all" ]}
        combine_hdvs_pops(hdv_dataset.out)
        // hdvs_pops = combine_hdvs.out.map{ dataset, top_hdv, site_hdv, all_hdv, gene_hdv, clndn_hdv -> [ 'AIBST', dataset+'=='+file(site_hdv) ]}.groupTuple().map{ dataset, sites -> [dataset, sites, "hdv_pgx_pop_${fold}fold_ac${ac}_g${group_treshold}"] }
        // all_hdv_pops = combine_hdvs.out.map{ dataset, top_hdv, site_hdv, all_hdv, gene_hdv, clndn_hdv -> [ 'AIBST', dataset+'=='+file(all_hdv) ]}.groupTuple().map{ dataset, all_hdvs -> [dataset, all_hdvs.join(','), "hdv_pgx_pop_${fold}fold_ac${ac}_g${group_treshold}"] }
        // generate_circos_dataframe(hdvs_pops)
        // generate_circos_plot(generate_circos_dataframe.out)

    emit:
        // data
        hdvs_pops = combine_hdvs_pops.out.hdvs_pops
}


//// Variant counts

workflow count_pgx_pop{
    take: data
    main:
        count_data = data.map{ dataset, chrm, pop, pop_vcf, pop_sample, dataset_vcf, dataset_singl_sites, dataset_singls -> [ pop, file(pop_vcf), dataset, file(dataset_singl_sites), chrm ] }
        count_dataset(count_data)
        // novel_dataset(count_data)
        // //// Get sites
        // get_sites = novel_dataset.out.novel_in_dbs
        //     .map{ pop, pop_vcf, not_other_dbs, novel_sites, chrm -> [ pop, pop_vcf, novel_sites, chrm ]}
        // get_sites_simple(get_sites)
        // //// Circos plot
        // generate_circos_dataframe(get_sites_simple.out.map{ dataset, vcf, site, chrm -> [ 'AIBST', dataset+'=='+file(site) ]}.groupTuple().map{ dataset, sites -> [dataset, sites, 'novel_pgx_pop'] })
        // generate_circos_plot(generate_circos_dataframe.out)
        // //// Upset plot
        // cat_sites(get_sites_simple.out.map{ pop, vcf, sites, chrm -> [ pop, file(sites), chrm] }.groupTuple())
        // upset_data = cat_sites.out
        //     .map{ pop, sites -> [ 'AIBST', pop, file(sites)] }
        //     .groupTuple(by:[0])
        //     .map{ dataset, pops, sites -> [ dataset, sites ]}
        // data_for_upset(upset_data)
        // plot_upset_R(data_for_upset.out)
    emit:
        counts = count_dataset
}

workflow count_all_dataset{
    take: data
    main:
        count_data = data.map{ dataset, chrm, dataset_vcf, dataset_sample, dataset_vcf_, dataset_singl_sites, dataset_singls -> [ dataset, file(dataset_vcf), dataset, file(dataset_singl_sites), chrm ] }
        count_dataset(count_data)
        novel_dataset(count_data)
    emit:
        counts = count_dataset
}

workflow count_pgx_dataset {
    take: data
    main:
        count_data = data.map{ dataset, chrm, dataset_vcf, dataset_sample, dataset_vcf_, dataset_singl_sites, dataset_singls -> [ dataset, file(dataset_vcf), dataset, file(dataset_singl_sites), chrm ] }
        count_dataset(count_data)
        novel_dataset(count_data)
        novel_dataset_vcfs = novel_dataset.out.novel_dataset_vcf.groupTuple(by:0)
        concat_vcf(novel_dataset_vcfs)
    emit:
        count_dataset
        novel_vcf = novel_dataset.out.novel_dataset_vcf
}

workflow count_exome_dataset{
    take: data
    main:
        count_data = data.map{ dataset, chrm, dataset_vcf, dataset_sample, dataset_vcf_, dataset_singl_sites, dataset_singls -> [ dataset, file(dataset_vcf), dataset, file(dataset_singl_sites), chrm ] }
        count_dataset(count_data)
        novel_dataset(count_data)
    emit:
        count_dataset
}


workflow count_all_ind {
    take: data
    main:
        count_data = data.map{ dataset, chrm, sample_vcf, sample_id, dataset_vcf, dataset_singl_sites, dataset_singls -> [ sample_id, file(sample_vcf), dataset, file(dataset_singl_sites), chrm ] }
        count_dataset(count_data)
        novel_dataset(count_data)
    emit:
        data
        // count_dataset
}

workflow count_all_ind_pgx {
    take: data
    main:
        count_data = data.map{ dataset, chrm, sample_vcf, sample_id, dataset_vcf, dataset_singl_sites, dataset_singls -> [ sample_id, file(sample_vcf), dataset, file(dataset_singl_sites), chrm ] }
        count_dataset(count_data)
        novel_dataset(count_data)
    emit:
        data
        // count_dataset
}

workflow count_all_ind_exome {
    take: data
    main:
        count_data = data.map{ dataset, chrm, sample_vcf, sample_id, dataset_vcf, dataset_singl_sites, dataset_singls -> [ sample_id, file(sample_vcf), dataset, file(dataset_singl_sites), chrm ] }
        count_dataset(count_data)
        novel_dataset(count_data)
    emit:
        data
        // count_dataset
}

workflow get_sample_vcfs {
    take: data
    main:
        get_pop(data)
        samples_data = get_pop.out
            .flatMap{ dataset, dataset_vcf, pop_info, chrm ->
                datas = []
                if(dataset == 'AIBST'){
                    file(pop_info)
                        .readLines()
                        .each{ sample_id ->
                            datas << [ sample_id, file(dataset_vcf), dataset, chrm]
                    }
                }
                return datas
            }

        // All variants
        get_vcf_ind(samples_data)

        // Get Coding regions
        get_protein_coding_variants(get_vcf_ind.out)
        samples_exome = get_protein_coding_variants.out
                            .combine(get_vcf_ind.out, by:[0,3])
                            .map{ dataset_id, chrm, sample_vcf, sample_file, dataset_vcf, sample_id -> [ sample_id, file(sample_vcf), dataset_id, chrm ] }

        // Get Pgx
        ind_pgx_data = get_vcf_ind.out
                        .map{ dataset, vcf, sample_id, chrm -> [ sample_id, file(vcf), dataset, chrm ] }
                        .combine(Channel.from(params.gene_lists).map{ name, bed -> [name, file(bed)]})
        get_gene_vcf(ind_pgx_data)
        samples_pgx = get_gene_vcf.out
                        .combine(get_vcf_ind.out.map{ dataset, vcf, sample_id, chrm -> [ sample_id, file(vcf), dataset, chrm ] }, by:[0,3])
                        .map{ sample_id, chrm, sample_vcf, sample_file, dataset_vcf, dataset -> [ dataset, file(sample_vcf), sample_id, chrm ] }

    emit:
        samples_vcf = get_vcf_ind.out
        samples_exome_vcf = samples_exome
        samples_pgx_vcf = samples_pgx
}

workflow novel_in_existing_database_wf {
    take: data
    main:
        all_chrm_dataset_data = data
            .flatMap{ dataset1, vcf1, sample_id1, chrm1 ->
                datas = [] 
                if(dataset1 == 'AIBST'){
                    datas << [ dataset1, file(vcf1), chrm1 ]
                    return datas
                }
                if(sample_id1 == 'AIBST'){
                    datas << [ dataset1, file(vcf1), chrm1 ]
                    return datas
                }
            }
        novel_in_existing_database( all_chrm_dataset_data )
    emit:
        novel_in_dbs = novel_in_existing_database.out
}

workflow combine_novel_count1_wf{
    take: data
    main:
        combine_data = data
            .flatMap{ dataset1, vcf1, count1, csv1, chrm -> 
                datas1 = []
                datas1 << [ dataset1, file(count1) ] 
                return datas1
            }.groupTuple()
            .flatMap{ dataset2, counts2 -> 
                datas2 = []
                datas2 << [ dataset2, counts2, 'false', 'all' ]
                return datas2
            }
        combine_novel_count1(combine_data)
        combine_counts_data = combine_novel_count1.out
            .groupTuple(by:[3])
            .map{ datasets, count_totals, count_summaries, group -> [ group, count_totals, true, group ]}
        // Combine totals
        combine_combine_counts(combine_counts_data)
    emit:
        novel_counts = combine_novel_count1.out
}

workflow count_dataset {
    take: data
    main:
        chrm_dataset_data = data
            .flatMap{ dataset, vcf, sample_id, singletons, chrm ->
                datas = []
                // For datasets
                if(dataset == 'AIBST'){
                    datas << [ dataset, file(vcf), sample_id, file(singletons), chrm ]
                    return datas
                }
                // Individuals
                if(sample_id == 'AIBST'){
                    datas << [ dataset, file(vcf), sample_id, file(singletons), chrm ]
                    return datas
                }
            }
        count_variants(chrm_dataset_data).view()
        dataset_counts = count_variants.out
            .flatMap{ dataset, csv, sample_id, chrm ->
                datas1 = [] 
                datas1 << [ dataset, dataset, file(csv), chrm ]
                return datas1
            }
            .groupTuple()
            .flatMap{ dataset, sample_ids, count_csvs, chrms ->
                datas2 = [] 
                datas2 << [ dataset, count_csvs, true, 'all' ]
                return datas2
            }
        combine_counts(dataset_counts) // Variant counts in AIBST
        // combine_counts_data = combine_counts.out
        //     .groupTuple(by:[3])
        //     .map{ datasets, count_totals, count_summaries, group -> [ group, count_totals, true, group ]}
        // // Combine totals
        // combine_combine_counts(combine_counts_data)
    emit:
        data
        // combine_counts
        // novel_in_existing_database.out        // combine_novel_count1.out
}

workflow novel_dataset {
    take: data
    main:
        count_data = data.map{ dataset, vcf, sample_id, singletons, chrm -> [ dataset, file(vcf), sample_id, chrm ] }
        // Not in other DBs
        novel_in_existing_database_wf( count_data )
        combine_novel_count1_wf(novel_in_existing_database_wf.out.novel_in_dbs) // AIBST novel not in other dataset counts
        combine_novel_count1_wf.out.novel_counts
        novel_data = novel_in_existing_database_wf.out.novel_in_dbs
            .map{ dataset, vcf, not_in_dbs, novel_sites, chrm ->
                if (!(file(novel_sites).isEmpty())) {
                    return [dataset, file(vcf), '', chrm, 'novel', file(novel_sites)]
                }
            }
        // Novel sites count
        get_novel_vcf(novel_data)
        count_data = get_novel_vcf.out
            .map{ dataset, dataset_novel_vcf, label, chrm -> [ dataset, file(dataset_novel_vcf), dataset, label, chrm ] }
            .combine( data, by:[0,4] )
            .map{ dataset, chrm, dataset_novel_vcf, dataset_id1, label, dataset_vcf, dataset_id2, dataset_singletons -> [ dataset, file(dataset_novel_vcf), dataset_id2, file(dataset_singletons), chrm ] }
        count_variants(count_data)
        dataset_novel_counts = count_variants.out
            .flatMap{ dataset2, csv2, sample_id2, chrm2 ->
                datas = []
                // For datasets
                if(dataset2 == 'AIBST'){
                    datas << [ dataset2, dataset2, file(csv2), chrm2 ]
                    return datas
                }
                // Individuals
                if(sample_id2 == 'AIBST'){
                    datas << [ dataset2, sample_id2, file(csv2), chrm2 ]
                    return datas
                }
            }
            .groupTuple(by:[0,1])
            .flatMap{ dataset3, sample_id, count_csvs3, chrms3 ->
                datas1 = [] 
                datas1 << [ dataset3, count_csvs3, true, sample_id]
                return datas1
            }
        combine_counts( dataset_novel_counts ) // Novel variant counts in AIBST
        combine_counts_data = combine_counts.out
            .groupTuple(by:[3])
            .map{ datasets, count_totals, count_summaries, group -> [ group, count_totals, true, group ]}
        // Combine totals
        combine_combine_counts(combine_counts_data)
    emit:
        data
        novel_in_dbs = novel_in_existing_database_wf.out.novel_in_dbs
        novel_pop_vcf = count_variants.out
        novel_dataset_vcf = get_novel_vcf.out
        // novel_in_existing_database.out        // combine_novel_count1.out

}

workflow combine_combine_counts {
    take: data
    main:
        combine_counts( data )//
    emit:
        data
}
// workflow count_all_chrm_dataset {
//     take: data
//     main:
//         all_chrm_dataset_data = data
//             .flatMap{ dataset, vcf, sample, chrm ->
//                 datas = [] 
//                 if(dataset == 'AIBST'){
//                     datas << [ dataset, file(vcf), chrm ]
//                 }
//                 return datas
//             }
//         count_variants(all_chrm_dataset_data)
//         dataset_counts = count_variants.out
//             .flatMap{ dataset, csv, chrm ->
//                 datas = [] 
//                 datas << [ 'AIBST', dataset, file(csv), chrm ]
//                 return datas
//             }
//             .groupTuple()
//             .flatMap{ dataset, pops, count_csvs, chrms ->
//                 datas = [] 
//                 datas << [ dataset, count_csvs, true, 'all' ]
//                 return datas
//             }
//         combine_counts(dataset_counts) // Novel variant counts in AIBST
//     emit:
//         data
//         // novel_in_existing_database.out        // combine_novel_count1.out
// }

workflow novel_all_chrm_pop {
    take: data
    main:
        // Novel per population
        novel_in_existing_database(data.map{ pop, dataset, vcf, sample -> [ pop, file(vcf), dataset] })
        get_freq(novel_in_existing_database.out.map{ pop, vcf, count, csv, chrm -> [ pop, file(vcf), file(count), chrm ]})
        get_novel_sites(get_freq.out.map{ dataset, vcf, sample, sites, mafs, frq, chrm -> [ dataset, file(vcf), file(mafs), '' ] })
        get_gene_vcf(get_novel_sites.out.map{ dataset, vcf, sites, chrm -> [ dataset, file(vcf), '', chrm, '', file(sites)] })
        novel_count(get_gene_vcf.out.map{ dataset, vcf, sample, chrm -> [ dataset, file(vcf), chrm ] })
        novel_count_combine(novel_count.out.map{ dataset2, csv2, chrm2 -> [ 'AIBST', dataset2, file(csv2), chrm2] }.groupTuple().map{ dataset, pops, count_csvs, chrms -> [ dataset, pops, count_csvs, chrms, 'all']})
    emit:
        data
        // novel_count_combine.out
}

workflow novel_pgx_dataset {
    take: data
    main:
        novel_in_existing_database(data.map{ dataset, vcf, csv, chrm -> [ dataset, file(vcf), chrm] }) // Soft novel sites, not in dbsnp
        combine_novel_count5(novel_in_existing_database.out.map{ dataset, vcf, count, csv, chrm -> [ dataset, file(count)] }.groupTuple().map{ dataset, counts -> [ dataset, counts, 'false', 'pgx']})
        get_freq(novel_in_existing_database.out.map{ pop, vcf, count, csv, chrm -> [ pop, file(vcf), file(count), chrm ]})
        get_novel_sites(get_freq.out.map{ dataset, vcf, sample, sites, mafs, frq, chrm -> [ dataset, file(vcf), file(mafs), '' ] })
        get_gene_vcf(get_novel_sites.out.map{ dataset, vcf, sites, chrm -> [ dataset, file(vcf), '', chrm, '', file(sites)] })
        novel_count(get_gene_vcf.out.map{ dataset, vcf, sample, chrm -> [ dataset, file(vcf), chrm ] })
        novel_count_combine(novel_count.out.map{ dataset2, csv2, chrm2 -> [ 'AIBST', dataset2, file(csv2), chrm2] }.groupTuple().map{ dataset, pops, count_csvs, chrms -> [ dataset, pops, count_csvs, chrms, 'pgx']})
    emit:
        data
        // combine_novel_count5.out
}

// workflow novel_pgx_pop {
//     take: data
//     main:
//         novel_in_existing_database(data.map{ pop, vcf, sample, chrm -> [ pop, file(vcf), pop] }) // Soft novel sites, not in dbsnp
//         combine_novel_count3(novel_in_existing_database.out.map{ pop, vcf, count, csv, chrm -> [ "AIBST", file(count)] }.groupTuple().map{ dataset, csvs -> [ dataset, csvs, 'false', 'pop_pgx'] })
//         get_freq2(novel_in_existing_database.out.map{ pop, vcf, count, csv, chrm -> [ pop, file(vcf), file(count), chrm ]})
//         get_novel_sites(get_freq2.out.map{ dataset, vcf, sample, sites, mafs, frq, chrm -> [ dataset, file(vcf), file(mafs), '' ] })
//         get_gene_vcf(get_novel_sites.out.map{ dataset, vcf, sites, chrm -> [ dataset, file(vcf), '', chrm, '', file(sites)] })
//         novel_count(get_gene_vcf.out.map{ dataset, vcf, sample, chrm -> [ dataset, file(vcf), chrm ] })
//         novel_count_combine(novel_count.out.map{ dataset2, csv2, chrm2 -> [ 'AIBST', dataset2, file(csv2), chrm2] }.groupTuple().map{ dataset, pops, count_csvs, chrms -> [ dataset, pops, count_csvs, chrms, 'pop_pgx']})
//         // get_freq2.out.map{ dataset, vcf, sample, sites, mafs, frq, chrm -> [ 'AIBST', dataset+'=='+file(sites) ]}
//         generate_circos_dataframe(get_freq2.out.map{ dataset, vcf, sample, sites, mafs, frq, chrm -> [ 'AIBST', dataset+'=='+file(sites) ]}.groupTuple().map{ dataset, sites -> [dataset, sites, 'novel_pgx_pop'] })
//         // generate_circos_plot(generate_circos_dataframe.out)
//      emit:
//         data
//         // combine_novel_count3.out
//         // combine_novel_count4.out
// }

// workflow stats_pgx {
//     take: data
//     main:
//         get_gene_vcf1(data.combine(params.drug_disease_variants).map{ dataset, vcf, sample, list_name, list_file -> [ dataset, file(vcf), file(sample), '', list_name, file(list_file) ]})
//         get_freq(get_gene_vcf1.out)
//         // get_freq.out
//         // get_freq1.out[0]
//         // cat_freq(get_freq.out.map{ dataset, vcf, sample, sites, frq, frq1, chrm -> [ dataset, file(frq), chrm] }.groupTuple())

//         // freq_for_upset(cat_freq.out)
//         // plot_upset_R(freq_for_upset.out)
//     emit:
//         data
        
// }

workflow csq_pgx {
    take: data
    main:
        csq(data)
        plot_csq(csq.out)
    emit:
        data
        
}

workflow group_pops_pgx{
    take: data
    main:
        groups_data_pgx = data
            .flatMap{ def pop, def pop_vcf, def pop_sample, def chrm ->
                datas = []
                params.dataset_groups.each{ def group, def group_pops ->
                    pops = group_pops.split(',')
                    if(pop in pops){
                        datas << [ group, pop, file(pop_vcf), file(pop_sample) ]
                    }
                }
                return datas
            }
            .groupTuple()
            .map{ def group, def pops, def pop_vcfs, def pop_samples ->
                if(pops.size() > 1){ [ group, '', pop_vcfs, pop_samples ] }
            }
        merge_pop_groups(groups_data_pgx)
    emit:
        groups_pgx = merge_pop_groups.out
}

workflow group_pops {
    take: data
    main:
        pop_data = data.groupTuple(by:[0,1])
            .map{ def pop, def dataset, def pop_vcfs, def pop_samples, def chrms -> [ pop, pop_vcfs, pop_samples ]}
        cat_chrm_groups(pop_data)
    emit:
        pop_data = cat_chrm_groups.out
}

workflow group_pops_by_chrm{
    take: data
    main:
        groups_data_chrm = data
            .flatMap{ def pop, def dataset, def pop_vcf, def pop_sample, def chrm ->
                datas = []
                params.dataset_groups.each{ def group, def group_pops ->
                    pops = group_pops.split(',')
                    if(pop in pops){
                        datas << [ chrm, group, pop, file(pop_vcf), file(pop_sample) ]
                    }
                }
                return datas
            }
            .groupTuple(by:[0,1])
            .flatMap{ def chrm, def group, def pops, def pop_vcfs, def pop_samples ->
                datas1 = []
                if(pops.size() > 1){ 
                    datas1 << [ group, chrm, pop_vcfs, pop_samples ]
                    return datas1
                }
            }
        merge_pop_groups(groups_data_chrm)
    emit:
        groups_chrm = merge_pop_groups.out
}

workflow fisher_group_pgx{
    take: data
    main:
        fisher_group_pgx_cha = data
            .flatMap{ def group_5, def chrm_5, def group_vcf_5, def group_sample_5 ->
                datas = []
                fisher_test_group = params.FISHER_TEST.split(',')
                fisher_test_group.each{ def group_f ->
                    if(group_f == group_5){
                        datas << [ group_5, file(group_vcf_5), file(group_sample_5), chrm_5 ]
                    }
                }
                return datas
            }
        vcf_to_plink(fisher_group_pgx_cha)
        generate_pheno_fam_sample(vcf_to_plink.out)
        fisher_test_plink(generate_pheno_fam_sample.out.map{ def dataset, def ped, def bim, def fam, def pheno, def chrm -> [dataset, file(ped), file(bim), file(fam), file(pheno), chrm, 'pgx']})
    emit:
        pvalues = fisher_test_plink.out
}

workflow fisher_group{
    take: data
    main:
        fisher_group_cha = data
            .flatMap{ def group_4, def chrm_4, def group_vcf_4, def group_sample_4 ->
                datas = []
                fisher_test_group = params.FISHER_TEST.split(',')
                fisher_test_group.each{ def group_f ->
                    if(group_f == group_4){
                        datas << [ group_4, file(group_vcf_4), file(group_sample_4), chrm_4 ]
                    }
                }
                return datas
            }
        vcf_to_plink(fisher_group_cha)
        generate_pheno_fam_sample(vcf_to_plink.out)
        fisher_test_plink(generate_pheno_fam_sample.out.map{ def dataset, def ped, def bim, def fam, def pheno, def chrm -> [dataset, file(ped), file(bim), file(fam), file(pheno), chrm, '']})
        combine_csv(fisher_test_plink.out.groupTuple(by:[0]))
    emit:
        // data
        pvalues = fisher_test_plink.out
}

workflow fst {
    take: data
    main:
        // pops = params.FST_TEST.split(',')
        // data_2pops = data.combine(data, by:[0,1,2,3])
        
        // pop2_per_group_cha = data_2pops.flatMap { dataset, var, dataset_vcf, dataset_sample, pop1, pop1_vcf, pop1_sample, pop2, pop2_vcf, pop2_sample  ->
        //             pop2_per_group = []
        //             tmp = []
        //             if ( pop1 in pops && pop2 in pops && pop1 != pop2 && !([pop2, pop1] in tmp) ){
        //             //     if(pop1 in [ 'KNK', 'KNL', 'KNM', 'KNP', 'NGH', 'NGI', 'NGY', 'SAV', 'TZA', 'TZB', 'ZWD', 'ZWS' ] && pop2 in [ 'KNK', 'KNL', 'KNM', 'KNP', 'NGH', 'NGI', 'NGY', 'SAV', 'TZA', 'TZB', 'ZWD', 'ZWS' ]){
        //                 pop2_per_group << [ dataset, pop1, file(pop1_vcf), file(pop1_sample), pop2, file(pop2_vcf), file(pop2_sample)]
        //                 tmp << [pop1, pop2]
        //                 }
        //             // }
        //     return pop2_per_group
        // }}
        fst_analysis(data)
        /* Everything, use for matrix and mean weir */
        get_fst_weir_estimates(fst_analysis.out)
        combine_fst_weir_estimates(get_fst_weir_estimates.out.groupTuple(by:[0]).map{ def dataset, def pops, def weir_fsts, def mean_fsts -> [ dataset, weir_fsts ]})
        /* FST > 0.5 */
        get_fst_weir_estimates_cutoff(fst_analysis.out.map{ def dataset, def pops, def fst_weir_estimate, def fst_log -> [ dataset, pops, fst_weir_estimate, fst_log, params.fst_cutoff]})
        combine_fst_weir_estimates_cutoff(get_fst_weir_estimates_cutoff.out.groupTuple(by:[0]).map{ def dataset, def pops, def weir_fsts, def mean_fsts -> [ dataset, weir_fsts ]})
        combine_fst_mean_weir_estimates_1(
            get_fst_weir_estimates.out.map{ def dataset, def pops, def weir_fst, def mean_fst -> [ pops.split('__')[0], weir_fst]}.mix(get_fst_weir_estimates.out.map{ def dataset, def pops, def weir_fst, def mean_fst -> [ pops.split('__')[1], weir_fst]}).groupTuple(by:[0])
        )
        combine_fst_mean_weir_estimates(get_fst_weir_estimates.out.groupTuple(by:[0]).map{ def dataset, def pops, def weir_fst, def mean_fst -> [ dataset, mean_fst ]})
        /// Plot matrix
        generate_fst_matrix(combine_fst_mean_weir_estimates.out)
        // plot_fst_matrix(generate_fst_matrix.out)
        // combine_weir_fst_analysis(get_fst_weir_estimates.out.map{ dataset, pops, weir_fst, mean_fst -> [ dataset, weir_fst ]}.groupTuple(by:[0]).map{ dataset, weir_fsts -> [ dataset, weir_fsts.join(','), '0.25']})
    emit:
        data
        fst_2pops_all = combine_fst_weir_estimates.out
        fst_2pops = combine_fst_mean_weir_estimates_1.out
}

workflow group_2pops {
    take: data
    main:
        // pops = params.FST_TEST.split(',')
        test_pops = Channel.from(params.FST_TEST.split(','))
        tmp1 = []
        
        my_pops = test_pops.combine(test_pops)
                .flatMap{ def pops_2 ->
                my_pops_ = []
                    if(pops_2[0] != pops_2[1]){
                        if ([pops_2[1], pops_2[0]] !in tmp1){
                            tmp1 << [pops_2[0], pops_2[1]]
                            my_pops_ << ["${pops_2[0]}_${pops_2[1]}"]
                        }
                    }
                    return my_pops_
                }
        data_2pops = data.combine(data)
        pop2_per_group_cha = data_2pops.map { def pop1, def pop1_vcf, def pop1_sample, def pop2, def pop2_vcf, def pop2_sample  -> [ "${pop1}_${pop2}", 'ALL', [file(pop1_vcf), file(pop2_vcf)], [file(pop1_sample), file(pop2_sample)] ]}.combine(my_pops, by:0)

        merge_pop_groups(pop2_per_group_cha)
    emit:
        pops2_data = pop2_per_group_cha
        pops_vcf = merge_pop_groups.out
}

workflow group_2pops_pgx {
    take: data
    main:
        // pops = params.FST_TEST.split(',')
        test_pops = Channel.from(params.FST_TEST.split(','))
        tmp1 = []
        
        my_pops = test_pops.combine(test_pops)
                .flatMap{ def pops_2 ->
                my_pops_ = []
                    if(pops_2[0] != pops_2[1]){
                        if ([pops_2[1], pops_2[0]] !in tmp1){
                            tmp1 << [pops_2[0], pops_2[1]]
                            my_pops_ << ["${pops_2[0]}_${pops_2[1]}"]
                        }
                    }
                    return my_pops_
                }
        data_2pops = data.combine(data)
        pop2_per_group_cha = data_2pops.map { def pop1, def pop1_vcf, def pop1_sample, def dataset1, def pop2, def pop2_vcf, def pop2_sample, def dataset2  -> [ "${pop1}_${pop2}", 'ALL', [file(pop1_vcf), file(pop2_vcf)], [file(pop1_sample), file(pop2_sample)] ]}.combine(my_pops, by:0)

        merge_pop_groups(pop2_per_group_cha)
    emit:
        pops2_data = pop2_per_group_cha
        pops_vcf = merge_pop_groups.out
}

workflow fisher_2pops_pgx{
    take: data
    main:
        vcf_to_plink(data.map{ def pops, def dataset, def pops_vcf, def pops_samples -> [  pops, pops_vcf, pops_samples, dataset ]})
        generate_pheno_fam_sample1(vcf_to_plink.out)
        fisher_test_plink_pops(generate_pheno_fam_sample1.out.map{ def dataset, def ped, def bim, def fam, def pheno, def chrm -> [dataset, file(ped), file(bim), file(fam), file(pheno), chrm, 'pgx']})
        filter_fisher_test(fisher_test_plink_pops.out)
        combine_csv(filter_fisher_test.out.groupTuple(by:[2]).map{ def pops, def fisher_files, def dataset ->  [ dataset, fisher_files, pops ] })
    emit:
        data
        pvalues = combine_csv.out
}

workflow fisher_2pops{
    take: data
    main:
        vcf_to_plink(data.map{ def pops, def dataset, def pops_vcf, def pops_samples -> [  pops, pops_vcf, pops_samples, dataset ]})
        generate_pheno_fam_sample1(vcf_to_plink.out)
        fisher_test_plink_pops(generate_pheno_fam_sample1.out.map{ def dataset, def ped, def bim, def fam, def pheno, def chrm -> [dataset, file(ped), file(bim), file(fam), file(pheno), chrm, 'pgx']})
        filter_fisher_test(fisher_test_plink_pops.out)
        combine_csv(filter_fisher_test.out.groupTuple(by:[2]).map{ def pops, def fisher_files, def dataset ->  [ dataset, fisher_files, pops ] })
    emit:
        data
        pvalues = combine_csv.out
}

// workflow {
//     main:
//         /// 0. Data processing
//         data_proc()  // Dataset: all chrm and pgx
//     emit:
//         data_proc.out.datasets_map_data
// }   

workflow {
    main:
        /// 0. Data processing
        data_proc()  // Dataset: all chrm and pgx
        // data_proc_pop(data_proc.out.fill_tags_dataset_chrm)// Population: all chrm and pgx
        data_proc_pop(data_proc.out.fill_tags_dataset_chrm)// Population: all chrm and pgx

        /// Group All variants data by chromosome
        group_pops_by_chrm(data_proc_pop.out.pop_exome_chrm)

        // Pop all chrms
        group_pops(data_proc_pop.out.pop_exome_chrm)
        group_2pops(group_pops.out.pop_data)


        /// Groups PGX data
        group_pops_pgx(data_proc_pop.out.pop_pgx)
        group_2pops_pgx(data_proc_pop.out.pop_pgx)

        /* Fisher test */
        //// fisher_group(group_pops_by_chrm.out.groups_chrm)
        //// fisher_group(group_pops_pgx.out.groups_pgx)
        fisher_group_pgx(group_pops_pgx.out.groups_pgx)
        fisher_2pops_pgx(group_2pops_pgx.out.pops_vcf)
        // fisher_2pops(group_2pops.out.pops_vcf)
        

        /* Fst for all populations  */
        // fst(group_2pops_pgx.out.pops2_data)
        fst(group_2pops.out.pops2_data) //*
        // fst(group_pops_pgx.out.combine(data_proc_pop.out.pop_pgx.map{ pop, pop_vcf, pop_sample, dataset -> [dataset, pop, file(pop_vcf), file(pop_sample)]}, by:0))
        // TODO

        /* PCA Analysis */
        // pca(group_pops_by_chrm.out.groups_chrm)


        //# 1. HDV
        //## 1.1. Dataset
        //### 1.1.1 All variants -- Not needed
        // hdv_all_chrm_dataset_aibst = data_proc.out.fill_tags_dataset_chrm
        //     .flatMap{ dataset, vcf, sample, chrm ->
        //         datas = [] 
        //         if(dataset == 'AIBST'){
        //             datas << [ dataset, file(vcf), file(sample), chrm ]
        //         }
        //         return datas
        //     }
        // hdv_all_chrm_dataset(hdv_all_chrm_dataset_aibst)
        
        // /* 1.1.2 PGx variants  */
        // hdv_pgx_dataset_cha = data_proc.out.pgx_dataset_chrm
        //     .combine(fisher_group.out.pvalues.map{ dataset, fisher_file, chrm -> [ dataset, file(fisher_file), '', chrm ]}, by:[3])
        //     .flatMap{ chrm, dataset, vcf, sample, group, pvalue_file, it ->
        //         datas = [] 
        //         if(dataset == 'AIBST'){
        //             datas << [ dataset, file(vcf), file(sample), chrm, group, file(pvalue_file) ]
        //         }
        //         return datas
        //     }
        // hdv_pgx_dataset(hdv_pgx_dataset_cha)
        // // hdv_pgx_dataset.out.get_freq_3


        /* 1.2. Population  */
        hdv_pgx_pop_cha = data_proc_pop.out.pop_pgx
            .map{ def pop, def pop_vcf, def pop_sample, def dataset -> [ pop, file(pop_vcf), dataset, '' ]}
            .combine(fisher_2pops_pgx.out.pvalues)
            .combine(fst.out.fst_2pops, by:[0])
            .flatMap{ def pop, def pop_vcf_chrm, def dataset, def chrm, def group, def pvalue_file, def fst_file  ->
                datas = [] 
                pops = params.FST_TEST.split(',')
                if(pop in pops){
                    datas << [ pop, file(pop_vcf_chrm), group, chrm, file(pvalue_file), file(fst_file) ]
                }
                return datas
            }
        hdv_pgx_pop(hdv_pgx_pop_cha)
        // annotate_hdvs_all(hdv_pgx_pop.out.hdvs_pops.combine(data_proc.out.all_datasets_annotated))





//         /* COUNTS */
//         /* 2. Novel */

//         /* 2.1. All chromosomes / Genome --- per chromosome */
//         // count_all_dataset(data_proc.out.fill_tags_dataset_chrm.combine(data_proc.out.singl_dataset_chrm, by:[0,3]))
//         // count_all_ind(data_proc.out.samples_vcf.combine(data_proc.out.singl_dataset_chrm, by:[0,3]))

//         /* 2.2. Pgx variants */
//         // count_pgx_dataset(data_proc.out.pgx_dataset_chrm.combine(data_proc.out.singl_dataset_chrm, by:[0,3]))
//         // count_all_ind_pgx(data_proc.out.samples_pgx_vcf.combine(data_proc.out.singl_dataset_chrm, by:[0,3]))

        // ## 2.3 Coding protein variants
        // count_exome_dataset(data_proc.out.exome_dataset_chrm.combine(data_proc.out.singl_dataset_chrm, by:[0,3]))
        // count_all_ind_exome(data_proc.out.samples_exome_vcf.combine(data_proc.out.singl_dataset_chrm, by:[0,3]))

        //### 2.1.2. Population
        // count_pop_chrm = data_proc_pop.out.pop_pgx_chrm
        //     .map{ pop, dataset, pop_vcf, pop_sample, chrm ->
        //         if(dataset == 'AIBST'){
        //             datas = [ dataset, pop, file(pop_vcf), chrm, file(pop_sample) ]
        //         }
        //         return datas
        //     }
        //     .combine(data_proc.out.singl_dataset_chrm, by:[0,3])
        // count_pgx_pop(count_pop_chrm)


//         /* Drug */ 
//         // stats_pgx(data_proc.out.pgx_dataset)
        


        /* CSQ */
        csq_data = data_proc.out.pgx_dataset
            .combine([file(params.so_term_file)])
            .flatMap{ dataset3, vcf3, sample3, so_term ->
                datas3 = [] 
                if(dataset3 == 'AIBST'){
                    datas3 << [ dataset3, file(vcf3), file(sample3), file(so_term) ]
                }
                return datas3
            }
        // csq_pgx(csq_data)


        println params.outDir
}



    


    // ////////////////////////////////////// VCF annotation with SnpEff
    // // Step 1: Annotate
    // check_files([params.dbsnp_vcf, "${params.dbsnp_vcf}.tbi", params.clinvar, "${params.clinvar}.tbi", params.cosmic, "${params.cosmic}.tbi"])
    // vcfs_mafs = get_gene_vcf.out.map{ it1 -> 
    //         combine_mafs(it1, params.mafs_annotations)
    //     }
    // combine_mafs_dataset(vcfs_mafs)
    // annotate_mafs (combine_mafs_dataset.out)
    // annotate_snpeff(annotate_mafs.out.combine([[params.snpeff_human_db, params.snpeff_database]]))
    // annotate_dbsnp(annotate_snpeff.out.combine([[file(params.dbsnp_vcf), file("${params.dbsnp_vcf}.tbi")]]))
    // annotate_clinvar(annotate_dbsnp.out.combine([[file(params.clinvar), file("${params.clinvar}.tbi")]]))
    // annotate_cosmic(annotate_clinvar.out.combine([[file(params.cosmic), file("${params.cosmic}.tbi")]]))
    
    // // Step 2: Split by population
    // // dataset_dataset_pops = annotate_cosmic.out.map{ it2 ->
    // //         add_pop_groups(it2, params.dataset_files)
    // //     }
    // dataset_dataset_pops = annotate_cosmic.out
    //     .flatMap{ dataset, dataset_vcf, dataset_sample, chrm ->
    //          add_pop_groups_datas = []
    //          params.dataset_files.each{ dataset_name, dataset_vcf_, dataset_sample_, dataset_pops ->
    //             dataset_pops.split(',').each{ pop ->
    //                 if(dataset == dataset_name){
    //                     add_pop_groups_datas << [pop, dataset, file(dataset_vcf), file(dataset_sample), chrm]
    //                 }
    //             }
    //         }
    //         return add_pop_groups_datas
    //     }
    // split_dataset_vcf_pop(dataset_dataset_pops) 
    // // Step 3: Concatenate vcf by population 
    // concat_dataset_pop(split_dataset_vcf_pop.out.groupTuple())
    // // Step 4: Group populations
    // // group_pops_data_cha = concat_dataset_pop.out[0]
    // //     .flatMap{ it -> 
    // //         group_pops_data(it, params.dataset_groups) }
    // //     .groupTuple()
    // //     .map{ group, pop_vcfs, pop_samples ->
    // //         if(pop_vcfs.size() > 1){ [group, pop_vcfs, pop_samples] }}
    // // merge_pop_groups(group_pops_data_cha)
    // // Step 5: Get singleton
    // // singleton_pop(concat_dataset_pop.out[0])
    // // Step 6: Novel per population
    // // novel_in_existing_database(concat_dataset_pop.out[1]) // Soft novel sites, not in dbsnp
    // novel_in_existing_database(concat_dataset_pop.out) // Soft novel sites, not in dbsnp
    // get_sites(novel_in_existing_database.out) | get_novel_sites // Generate list of novel sites, mot in other datasbases
    // // Remove empty list
    // novel_sites_noempty = get_novel_sites.out
    //     .map{ pop, vcf, bed_list ->
    //         if(bed_list.size() > 1){ [pop, file(vcf), file(bed_list)] }}
    // get_gene_vcf_simple(novel_sites_noempty) | novel_count  // Generate novel vcf and counts
    // group_novel_count = novel_count.out
    //     .map{ pop, csv ->
    //         [ 'AIBST', pop, file(csv)] }
    //     .groupTuple()
    // novel_count_combine(group_novel_count)
// }

