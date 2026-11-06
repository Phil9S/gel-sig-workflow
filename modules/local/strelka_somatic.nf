// Strelka2 somatic calling, used when a sample has no snv_vcf
process STRELKA_SOMATIC {
    tag "${meta.id}"
    label 'process_high'
    container "${params.container_strelka}"

    input:
    tuple val(meta),
          path(tumour_bam, stageAs: 'tumour/*'),
          path(tumour_bai, stageAs: 'tumour/*'),
          path(normal_bam, stageAs: 'normal/*'),
          path(normal_bai, stageAs: 'normal/*'),
          path(fasta),
          path(fasta_fai),
          path(call_regions),
          path(call_regions_tbi)

    output:
    tuple val(meta), path("${meta.id}.strelka.somatic_snvs.vcf.gz"), path("${meta.id}.strelka.somatic_snvs.vcf.gz.tbi"),     emit: snvs
    tuple val(meta), path("${meta.id}.strelka.somatic_indels.vcf.gz"), path("${meta.id}.strelka.somatic_indels.vcf.gz.tbi"), emit: indels

    script:
    def regions = call_regions ? "--callRegions ${call_regions}" : ''
    def mem_gb = [task.memory.toGiga() - 2, 2].max()
    """
    [ -e "${tumour_bam}.bai" ] || ln -s "\$(basename "${tumour_bai}")" "${tumour_bam}.bai"
    [ -e "${normal_bam}.bai" ] || ln -s "\$(basename "${normal_bai}")" "${normal_bam}.bai"
    # Strelka expects the index at <fasta>.fai
    [ -e "${fasta}.fai" ] || ln -s "${fasta_fai}" "${fasta}.fai"

    configureStrelkaSomaticWorkflow.py \\
        --tumorBam ${tumour_bam} \\
        --normalBam ${normal_bam} \\
        --referenceFasta ${fasta} \\
        ${regions} \\
        --runDir strelka

    python strelka/runWorkflow.py -m local -j ${task.cpus} -g ${mem_gb}

    mv strelka/results/variants/somatic.snvs.vcf.gz       ${meta.id}.strelka.somatic_snvs.vcf.gz
    mv strelka/results/variants/somatic.snvs.vcf.gz.tbi   ${meta.id}.strelka.somatic_snvs.vcf.gz.tbi
    mv strelka/results/variants/somatic.indels.vcf.gz     ${meta.id}.strelka.somatic_indels.vcf.gz
    mv strelka/results/variants/somatic.indels.vcf.gz.tbi ${meta.id}.strelka.somatic_indels.vcf.gz.tbi
    """

    stub:
    """
    touch ${meta.id}.strelka.somatic_snvs.vcf.gz ${meta.id}.strelka.somatic_snvs.vcf.gz.tbi
    touch ${meta.id}.strelka.somatic_indels.vcf.gz ${meta.id}.strelka.somatic_indels.vcf.gz.tbi
    """
}
