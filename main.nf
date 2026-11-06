#!/usr/bin/env nextflow
// Tumour/normal WGS copy number pipeline: Battenberg, then SNV VAF evaluation,
// downsampled sWGS and SNP6 ASCAT at the Battenberg purity/ploidy, then a cohort summary.
// Samplesheet columns are described in subworkflows/local/utils.nf.

include { validateParameters; readSamplesheet; readSwgsSmooth; addPurityPloidy } from './subworkflows/local/utils'
include { BATTENBERG                } from './modules/local/battenberg'
include { BATTENBERG_PURITY_PLOIDY  } from './modules/local/battenberg_purity_ploidy'
include { SNV_VAF_QC                } from './subworkflows/local/snv_vaf_qc'
include { SWGS_ABSOLUTECN           } from './subworkflows/local/swgs_absolutecn'
include { ASCAT_SNP6                } from './subworkflows/local/ascat_snp6'
include { CN_SUMMARY                } from './modules/local/cn_summary'

workflow {

    validateParameters()

    log.info """
     gel-sig-workflow copy number pipeline  v${workflow.manifest.version}
     samples        :${params.samples}
     outdir         : ${params.outdir}
     VAF evaluation : ${params.skip_vaf_evaluation ? 'skipped' : (params.call_snvs ? 'on (Strelka2 when no snv_vcf)' : 'on (supplied VCFs only)')}
     sWGS           : ${params.skip_swgs ? 'skipped' : "swgs-absolutecn v2.0.0, ${params.swgs_bin} kb bins"}
     ASCAT (SNP6)   : ${params.skip_ascat ? 'skipped' : 'on'}
     profile        : ${workflow.profile}
    """.stripIndent()

    // [ meta, tumour_bam, tumour_bai, normal_bam, normal_bai, snv_vcf | [] ]
    ch_samples = channel.fromList(readSamplesheet(params.samples))

    ch_bams = ch_samples.map { meta, tbam, tbai, nbam, nbai, _vcf -> tuple(meta, tbam, tbai, nbam, nbai) }

    // Battenberg
    ch_battenberg_in = ch_bams.map { meta, tbam, tbai, nbam, nbai ->
        tuple(meta, tbam, tbai, nbam, nbai, file(params.genomes[meta.build].battenberg_ref_dir, checkIfExists: true))
    }
    BATTENBERG(ch_battenberg_in)

    // Purity and ploidy per sample
    BATTENBERG_PURITY_PLOIDY(BATTENBERG.out.purity_ploidy)

    // [ meta, meta + purity/ploidy/psi ] for samples with a usable Battenberg fit
    ch_fit = BATTENBERG_PURITY_PLOIDY.out.tsv
        .map { meta, tsv -> tuple(meta, addPurityPloidy(meta, tsv)) }
        .filter { _meta, meta_fit -> meta_fit != null }

    // [ meta, purity_ploidy.tsv ] restricted to the same samples
    ch_purity_ok = BATTENBERG_PURITY_PLOIDY.out.tsv
        .join(ch_fit)
        .map { meta, tsv, _meta_fit -> tuple(meta, tsv) }

    // [ meta_fit, tumour_bam, tumour_bai, normal_bam, normal_bai ]
    ch_bams_fit = ch_bams
        .join(ch_fit)
        .map { _meta, tbam, tbai, nbam, nbai, meta_fit -> tuple(meta_fit, tbam, tbai, nbam, nbai) }

    ch_tables = BATTENBERG_PURITY_PLOIDY.out.tsv
        .map { _meta, tsv -> tsv }
        .collectFile(name: 'battenberg_purity_ploidy.tsv', keepHeader: true, skip: 1, sort: { f -> f.name }, storeDir: "${params.outdir}/summary")

    // Check Battenberg clonal segments against somatic SNV VAFs
    if (!params.skip_vaf_evaluation) {
        SNV_VAF_QC(ch_samples, BATTENBERG.out.copynumber, ch_purity_ok)
        ch_tables = ch_tables.mix(SNV_VAF_QC.out.metrics_table)
    }

    // Downsampled sWGS absolute copy number
    if (!params.skip_swgs) {
        // smoothing is added here, not in the shared meta, so changing it does not re-run Battenberg
        def swgs_smooth = readSwgsSmooth(params.samples)
        SWGS_ABSOLUTECN(ch_bams_fit.map { meta, tbam, tbai, _nbam, _nbai -> tuple(meta + [smooth: swgs_smooth[meta.id]], tbam, tbai) })
        ch_tables = ch_tables.mix(SWGS_ABSOLUTECN.out.fits_table)
    }

    // ASCAT at SNP6 loci
    if (!params.skip_ascat) {
        ASCAT_SNP6(ch_bams_fit)
        ch_tables = ch_tables.mix(ASCAT_SNP6.out.metrics_table)
    }

    // Cohort summary
    CN_SUMMARY(ch_tables.collect())
}
