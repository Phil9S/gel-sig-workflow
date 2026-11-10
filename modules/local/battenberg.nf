// Battenberg clonal and subclonal copy number for a tumour/normal WGS pair.
// Reference paths in params.genomes are relative to the staged reference bundle.
// Large intermediates stay in the work directory.
process BATTENBERG {
    tag "${meta.id}"
    label 'process_high'
    container "${params.container_battenberg}"

    input:
    tuple val(meta),
          path(tumour_bam, stageAs: 'tumour/*'),
          path(tumour_bai, stageAs: 'tumour/*'),
          path(normal_bam, stageAs: 'normal/*'),
          path(normal_bai, stageAs: 'normal/*'),
          path(ref_dir,    stageAs: 'battenberg_ref')

    output:
    tuple val(meta), path("${meta.tumour_name}_copynumber.txt"),                                  emit: copynumber
    tuple val(meta), path("${meta.tumour_name}_purity_ploidy.txt"), path("${meta.tumour_name}_rho_and_psi.txt"), emit: purity_ploidy
    tuple val(meta), path("${meta.tumour_name}_copynumber_extended.txt"),                         emit: copynumber_extended
    tuple val(meta), path("*.png"),                                                               emit: plots, optional: true
    tuple val(meta), path("${meta.tumour_name}{_segment_masking_details.txt,_refit_suggestions.txt,_copynumber_solution_status.txt,.BAFsegmented.txt,.logRsegmented.txt,_chrX_copynumber.txt}"), emit: aux, optional: true
    tuple val(meta), path("${meta.id}.battenberg.log"),                                           emit: log

    script:
    def g = params.genomes[meta.build]
    // optional reference arguments: omitted when unset for this genome build
    def ref_args = [
        impute_info:           g.battenberg_impute_info,
        g1000_loci_prefix:     g.battenberg_g1000_loci_prefix,
        g1000_alleles_prefix:  g.battenberg_g1000_alleles_prefix,
        gc_prefix:             g.battenberg_gc_prefix,
        rt_prefix:             g.battenberg_rt_prefix,
        problem_loci:          g.battenberg_problem_loci,
        chrom_coord_file:      g.battenberg_chrom_coord_file,
        beagle_jar:            g.battenberg_beagle_jar,
        beagle_ref_template:   g.battenberg_beagle_ref_template,
        beagle_plink_template: g.battenberg_beagle_plink_template
    ].findAll { k, v -> v != null && v.toString() != '' }
     .collect { k, v -> "--${k} '${v}'" }
     .join(' ')
    // keep the Beagle JVM heap inside the task allocation
    def beagle_mem = [[params.battenberg_beagle_max_mem.toString().toInteger(), task.memory.toGiga() - 4].min(), 1].max()
    """
    # Battenberg/htslib expect <bam>.bai next to each BAM
    [ -e "${tumour_bam}.bai" ] || ln -s "\$(basename "${tumour_bai}")" "${tumour_bam}.bai"
    [ -e "${normal_bam}.bai" ] || ln -s "\$(basename "${normal_bai}")" "${normal_bam}.bai"

    run_battenberg.R \\
        --samplename '${meta.tumour_name}' \\
        --normalname '${meta.normal_name}' \\
        --tumour_bam '${tumour_bam}' \\
        --normal_bam '${normal_bam}' \\
        --sex '${meta.sex}' \\
        --genomebuild '${meta.build}' \\
        --ref_dir battenberg_ref \\
        ${ref_args} \\
        --use_beagle ${params.battenberg_use_beagle} \\
        --beagle_max_mem ${beagle_mem} \\
        --beagle_window ${params.battenberg_beagle_window} \\
        --beagle_overlap ${params.battenberg_beagle_overlap} \\
        --cpus ${task.cpus} \\
        --min_ploidy ${params.battenberg_min_ploidy} \\
        --max_ploidy ${params.battenberg_max_ploidy} \\
        --min_rho ${params.battenberg_min_rho} \\
        --max_rho ${params.battenberg_max_rho} \\
        --min_goodness ${params.battenberg_min_goodness} \\
        --platform_gamma ${params.battenberg_platform_gamma} \\
        --phasing_gamma ${params.battenberg_phasing_gamma} \\
        --segmentation_gamma ${params.battenberg_segmentation_gamma} \\
        --segmentation_kmin ${params.battenberg_segmentation_kmin} \\
        --phasing_kmin ${params.battenberg_phasing_kmin} \\
        --clonality_dist_metric ${params.battenberg_clonality_dist_metric} \\
        --ascat_dist_metric ${params.battenberg_ascat_dist_metric} \\
        --uninformative_baf_threshold ${params.battenberg_uninformative_baf_threshold} \\
        --min_normal_depth ${params.battenberg_min_normal_depth} \\
        --min_base_qual ${params.battenberg_min_base_qual} \\
        --min_map_qual ${params.battenberg_min_map_qual} \\
        --calc_seg_baf_option ${params.battenberg_calc_seg_baf_option} \\
        --max_allowed_state ${params.battenberg_max_allowed_state} \\
        --enhanced_grid_search ${params.battenberg_enhanced_grid_search} \\
        2>&1 | tee ${meta.id}.battenberg.log
    """

    stub:
    """
    printf 'chr\\tstartpos\\tendpos\\tnMaj1_A\\tnMin1_A\\tfrac1_A\\tnMaj2_A\\tnMin2_A\\tfrac2_A\\n1\\t1\\t249000000\\t1\\t1\\t1\\tNA\\tNA\\tNA\\n' > ${meta.tumour_name}_copynumber.txt
    cp ${meta.tumour_name}_copynumber.txt ${meta.tumour_name}_copynumber_extended.txt
    printf 'purity\\tploidy\\tpsi\\n0.7\\t2.1\\t2.05\\n' > ${meta.tumour_name}_purity_ploidy.txt
    printf 'rho\\tpsi\\tploidy\\tdistance\\tis.best\\nASCAT\\t0.7\\t2.05\\t2.1\\tNA\\tNA\\nFRAC_GENOME\\t0.7\\t2.05\\t2.1\\t0.05\\tTRUE\\nREF_SEG\\t0.7\\t2.05\\t2.1\\t0.06\\tFALSE\\n' > ${meta.tumour_name}_rho_and_psi.txt
    touch ${meta.id}.battenberg.log
    """
}
