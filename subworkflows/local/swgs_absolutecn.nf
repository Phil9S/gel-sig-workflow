// Downsampled shallow WGS absolute copy number with swgs-absolutecn v2.0.0 stage 2,
// using the Battenberg purity/ploidy. One bin size per run.
include { SWGS_DOWNSAMPLE_QDNASEQ; SWGS_REL_TO_ABS; SWGS_COMBINE } from '../../modules/local/swgs'

workflow SWGS_ABSOLUTECN {

    take:
    ch_tumour       // [ meta (with purity, ploidy, smooth), tumour_bam, tumour_bai ]

    main:
    def fitting_dir = "${params.outdir}/sWGS_fitting/${params.swgs_project}_${params.swgs_bin}kb"

    SWGS_DOWNSAMPLE_QDNASEQ(ch_tumour)

    SWGS_REL_TO_ABS(SWGS_DOWNSAMPLE_QDNASEQ.out.rds.join(SWGS_DOWNSAMPLE_QDNASEQ.out.qc))

    SWGS_COMBINE(
        SWGS_REL_TO_ABS.out.results.map { _meta, fits, _rds, _seg -> fits }.collect(),
        SWGS_REL_TO_ABS.out.results.map { _meta, _fits, rds, _seg -> rds }.collect(),
        SWGS_REL_TO_ABS.out.results.map { _meta, _fits, _rds, seg -> seg }.collect()
    )

    // cohort QC table written by processPrecomputed.R in swgs-absolutecn
    SWGS_DOWNSAMPLE_QDNASEQ.out.qc
        .map { _meta, tsv -> tsv }
        .collectFile(name: "${params.swgs_project}_${params.swgs_bin}kb_combined_QC_predownsample.tsv",
                     keepHeader: true, skip: 1, sort: { f -> f.name },
                     storeDir: "${fitting_dir}/absolute_PRE_down_sampling")

    emit:
    abs_cn       = SWGS_REL_TO_ABS.out.results
    combined_rds = SWGS_COMBINE.out.rds
    fits_table   = SWGS_COMBINE.out.fits
}
