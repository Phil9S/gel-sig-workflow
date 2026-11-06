// swgs-absolutecn v2.0.0 stage 2 with Battenberg ploidy/purity, running the release's own
// scripts from params.swgs_scripts_dir in the v2.0.0 container.
// https://github.com/Phil9S/swgs-absolutecn/releases/tag/v2.0.0

// Downsampling and QDNAseq share a task because a sample kept at full depth is a symlink
// to the input BAM, which is only mounted in the task that staged it.
process SWGS_DOWNSAMPLE_QDNASEQ {
    tag "${meta.id}"
    label 'process_medium'
    container "${params.container_swgs}"

    input:
    tuple val(meta), path(bam, stageAs: 'input/*'), path(bai, stageAs: 'input/*')

    output:
    tuple val(meta), path("${params.swgs_project}_${meta.id}_fitted_QC.tsv"),                              emit: qc
    tuple val(meta), path("${params.swgs_project}_${meta.id}_${params.swgs_bin}kb_relSmoothedCN.rds"),    emit: rds
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"),                                   emit: bam

    script:
    def scripts = params.swgs_scripts_dir
    def qc = "${params.swgs_project}_${meta.id}_fitted_QC.tsv"
    """
    [ -e "${bam}.bai" ] || ln -s "\$(basename "${bai}")" "${bam}.bai"

    # processPrecomputed.R / splitInput.R: stage-2 QC row with downsample_depth
    swgs_precomputed_qc.R \\
        --sample '${meta.id}' \\
        --patient '${meta.id}' \\
        --bam ${bam} \\
        --bin ${params.swgs_bin} \\
        --genome ${meta.build} \\
        --ploidy ${meta.ploidy} \\
        --purity ${meta.purity} \\
        --smooth ${meta.smooth} \\
        --qc_out ${qc} \\
        --sheet_out sample_sheet.tsv \\
        --bamcheck_out ${meta.id}_bam.ok

    # placeholder for the stage-1 RDS (not read when ploidy/purity are precomputed)
    echo '${meta.id}' > precomputed_relSmoothedCN.rds

    Rscript ${scripts}/downsampleBams.R \\
        --bam ${bam} \\
        --meta ${qc} \\
        --rds precomputed_relSmoothedCN.rds \\
        --bamOut ${meta.id}.bam \\
        --sampleName '${meta.id}' \\
        --prplpu TRUE \\
        --fileType BAM

    Rscript ${scripts}/qdnaseqMod.R \\
        --bam ${meta.id}.bam \\
        --sampleName '${meta.id}' \\
        --meta sample_sheet.tsv \\
        --bin ${params.swgs_bin} \\
        --genome ${meta.build} \\
        --rdsOut ${params.swgs_project}_${meta.id}_${params.swgs_bin}kb_relSmoothedCN.rds \\
        --autoSmooth ${params.swgs_auto_smooth ? 'TRUE' : 'FALSE'} \\
        --smoothThreshold ${params.swgs_smooth_threshold} \\
        --use_seed ${params.swgs_use_seed ? 'TRUE' : 'FALSE'} \\
        --seed_val ${params.swgs_seed} \\
        --pairedEnd ${params.swgs_paired_end ? 'TRUE' : 'FALSE'}
    """

    stub:
    """
    printf 'SAMPLE_ID\\tploidy\\tpurity\\tdownsample_depth\\tuse\\n${meta.id}\\t${meta.ploidy}\\t${meta.purity}\\t1\\tTRUE\\n' > ${params.swgs_project}_${meta.id}_fitted_QC.tsv
    touch ${params.swgs_project}_${meta.id}_${params.swgs_bin}kb_relSmoothedCN.rds ${meta.id}.bam ${meta.id}.bam.bai
    """
}

process SWGS_REL_TO_ABS {
    tag "${meta.id}"
    label 'process_low'
    container "${params.container_swgs}"

    input:
    tuple val(meta), path(rds), path(qc)

    output:
    tuple val(meta),
          path("${params.swgs_project}_${meta.id}_${params.swgs_bin}kb_ds_abs_fits.tsv"),
          path("${params.swgs_project}_${meta.id}_${params.swgs_bin}kb_ds_absCopyNumber.rds"),
          path("${params.swgs_project}_${meta.id}_${params.swgs_bin}kb_ds_absCopyNumber_segTable.tsv"), emit: results
    tuple val(meta), path("plots/${meta.id}.png"),                                                     emit: plot

    script:
    def prefix = "${params.swgs_project}_${meta.id}_${params.swgs_bin}kb"
    """
    mkdir -p plots
    Rscript ${params.swgs_scripts_dir}/qdnaseqRelToAbs.R \\
        --rds ${rds} \\
        --meta ${qc} \\
        --genome ${meta.build} \\
        --sampleName '${meta.id}' \\
        --tsvOut ${prefix}_ds_abs_fits.tsv \\
        --rdsOut ${prefix}_ds_absCopyNumber.rds \\
        --tabOut ${prefix}_ds_absCopyNumber_segTable.tsv \\
        --plotOut plots/${meta.id}.png
    """

    stub:
    def prefix = "${params.swgs_project}_${meta.id}_${params.swgs_bin}kb"
    """
    mkdir -p plots && touch plots/${meta.id}.png ${prefix}_ds_absCopyNumber.rds
    printf 'SAMPLE_ID\\tsegments.post\\n${meta.id}\\t1\\n' > ${prefix}_ds_abs_fits.tsv
    printf 'chromosome\\tstart\\tend\\tsegVal\\tsample\\n' > ${prefix}_ds_absCopyNumber_segTable.tsv
    """
}

process SWGS_COMBINE {
    label 'process_low'
    container "${params.container_swgs}"

    input:
    path(fits, stageAs: 'fits/*')
    path(rds, stageAs: 'rds/*')
    path(segtables, stageAs: 'seg/*')

    output:
    path("${params.swgs_project}_${params.swgs_bin}kb_ds_abs_fits.tsv"),                emit: fits
    path("${params.swgs_project}_${params.swgs_bin}kb_ds_absCopyNumber.rds"),           emit: rds
    path("${params.swgs_project}_${params.swgs_bin}kb_ds_absCopyNumber_segTable.tsv"),  emit: segtable

    script:
    def prefix = "${params.swgs_project}_${params.swgs_bin}kb"
    """
    Rscript ${params.swgs_scripts_dir}/combineResults.R \\
        --tsv "\$(ls fits/*.tsv | paste -sd, -)" \\
        --rds "\$(ls rds/*.rds | paste -sd, -)" \\
        --tab "\$(ls seg/*.tsv | paste -sd, -)" \\
        --tsvOut ${prefix}_ds_abs_fits.tsv \\
        --rdsOut ${prefix}_ds_absCopyNumber.rds \\
        --tabOut ${prefix}_ds_absCopyNumber_segTable.tsv
    """

    stub:
    def prefix = "${params.swgs_project}_${params.swgs_bin}kb"
    """
    touch ${prefix}_ds_absCopyNumber.rds
    printf 'SAMPLE_ID\\n' > ${prefix}_ds_abs_fits.tsv
    printf 'chromosome\\tstart\\tend\\tsegVal\\tsample\\n' > ${prefix}_ds_absCopyNumber_segTable.tsv
    """
}
