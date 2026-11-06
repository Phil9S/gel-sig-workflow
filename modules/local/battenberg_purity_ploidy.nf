// One-row table of the Battenberg purity, ploidy and psi for a sample
process BATTENBERG_PURITY_PLOIDY {
    tag "${meta.id}"
    label 'process_single'
    container "${params.container_r}"

    input:
    tuple val(meta), path(purity_ploidy), path(rho_and_psi)

    output:
    tuple val(meta), path("${meta.id}.battenberg_purity_ploidy.tsv"), emit: tsv

    script:
    """
    extract_purity_ploidy.R \\
        --sample '${meta.id}' \\
        --tumour_name '${meta.tumour_name}' \\
        --normal_name '${meta.normal_name}' \\
        --cancer '${meta.cancer}' \\
        --build '${meta.build}' \\
        --sex '${meta.sex}' \\
        --purity_ploidy '${purity_ploidy}' \\
        --rho_psi '${rho_and_psi}' \\
        --out '${meta.id}.battenberg_purity_ploidy.tsv'
    """

    stub:
    """
    printf 'sample\\ttumour_name\\tnormal_name\\tcancer\\tbuild\\tsex\\tpurity\\tploidy\\tpsi\\tgoodness_of_fit\\tsolution\\n${meta.id}\\t${meta.tumour_name}\\t${meta.normal_name}\\t${meta.cancer}\\t${meta.build}\\t${meta.sex}\\t0.7\\t2.1\\t2.05\\t0.05\\tFRAC_GENOME\\n' > ${meta.id}.battenberg_purity_ploidy.tsv
    """
}
