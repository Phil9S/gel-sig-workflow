// One row per sample joining whichever Battenberg, VAF, ASCAT and sWGS tables were produced
process CN_SUMMARY {
    label 'process_single'
    container "${params.container_r}"

    input:
    path(tables, stageAs: 'tables/*')

    output:
    path("cn_pipeline_summary.tsv"), emit: summary

    script:
    """
    cn_summary.R --tables_dir tables --out cn_pipeline_summary.tsv
    """

    stub:
    """
    printf 'sample\\n' > cn_pipeline_summary.tsv
    """
}
