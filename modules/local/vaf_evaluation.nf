// Evaluate a Battenberg fit with somatic SNV VAFs on clonal segments
process VAF_EVALUATION {
    tag "${meta.id}"
    label 'process_low'
    container "${params.container_r}"

    input:
    tuple val(meta), path(copynumber), path(purity_ploidy), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.vaf_eval.metrics.tsv"), emit: metrics
    tuple val(meta), path("${meta.id}.vaf_eval.states.tsv"),  emit: states
    tuple val(meta), path("${meta.id}.vaf_eval.snvs.tsv.gz"), emit: snvs
    tuple val(meta), path("${meta.id}.vaf_eval.pdf"),         emit: plots

    script:
    """
    evaluate_vaf.R \\
        --sample '${meta.id}' \\
        --tumour_name '${meta.tumour_name}' \\
        --normal_name '${meta.normal_name}' \\
        --copynumber '${copynumber}' \\
        --purity_ploidy '${purity_ploidy}' \\
        --vcf '${vcf}' \\
        --min_depth ${params.vaf_min_depth} \\
        --min_alt_reads ${params.vaf_min_alt_reads} \\
        --pass_only ${params.vaf_pass_only} \\
        --autosomes_only ${params.vaf_autosomes_only} \\
        --max_total_cn ${params.vaf_max_total_cn} \\
        --min_snvs ${params.vaf_min_snvs} \\
        --min_snvs_per_state ${params.vaf_min_snvs_per_state} \\
        --ccf_tolerance ${params.vaf_ccf_tolerance} \\
        --max_exceed_fraction ${params.vaf_max_exceed_fraction} \\
        --out_prefix '${meta.id}.vaf_eval'
    """

    stub:
    """
    printf 'sample\\tqc_flag\\n${meta.id}\\tPASS\\n' > ${meta.id}.vaf_eval.metrics.tsv
    printf 'sample\\tstate\\n' > ${meta.id}.vaf_eval.states.tsv
    printf 'sample\\n' | gzip > ${meta.id}.vaf_eval.snvs.tsv.gz
    touch ${meta.id}.vaf_eval.pdf
    """
}
