// Evaluate Battenberg clonal segments with somatic SNV VAFs, from a supplied VCF or Strelka2 calls.
include { STRELKA_SOMATIC } from '../../modules/local/strelka_somatic'
include { VAF_EVALUATION  } from '../../modules/local/vaf_evaluation'

workflow SNV_VAF_QC {

    take:
    ch_samples      // [ meta, tumour_bam, tumour_bai, normal_bam, normal_bai, snv_vcf | [] ]
    ch_copynumber   // [ meta, <tumour>_copynumber.txt ]
    ch_purity       // [ meta, <sample>.battenberg_purity_ploidy.tsv ]

    main:
    ch_branched = ch_samples.branch { _meta, _tbam, _tbai, _nbam, _nbai, vcf ->
        supplied: vcf
        call:     params.call_snvs
        skip:     true
    }

    ch_branched.skip.subscribe { row ->
        log.warn("Sample '${row[0].id}': no snv_vcf supplied and --call_snvs false - VAF evaluation skipped")
    }

    ch_strelka_in = ch_branched.call.map { meta, tbam, tbai, nbam, nbai, _vcf ->
        def g = params.genomes[meta.build]
        def fasta = file(g.fasta, checkIfExists: true)
        def fai = file(g.fasta_fai ?: "${g.fasta}.fai", checkIfExists: true)
        def regions = g.strelka_call_regions ? file(g.strelka_call_regions, checkIfExists: true) : []
        def regions_tbi = g.strelka_call_regions ? file("${g.strelka_call_regions}.tbi", checkIfExists: true) : []
        tuple(meta, tbam, tbai, nbam, nbai, fasta, fai, regions, regions_tbi)
    }
    STRELKA_SOMATIC(ch_strelka_in)

    ch_vcf = ch_branched.supplied
        .map { meta, _tbam, _tbai, _nbam, _nbai, vcf -> tuple(meta, vcf) }
        .mix(STRELKA_SOMATIC.out.snvs.map { meta, vcf, _tbi -> tuple(meta, vcf) })

    VAF_EVALUATION(ch_copynumber.join(ch_purity).join(ch_vcf))

    ch_metrics_table = VAF_EVALUATION.out.metrics
        .map { _meta, tsv -> tsv }
        .collectFile(name: 'vaf_evaluation_metrics.tsv', keepHeader: true, skip: 1, sort: { f -> f.name }, storeDir: "${params.outdir}/summary")

    emit:
    vcf           = ch_vcf
    metrics       = VAF_EVALUATION.out.metrics
    metrics_table = ch_metrics_table
}
