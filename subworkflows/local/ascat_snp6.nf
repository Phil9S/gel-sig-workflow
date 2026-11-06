// ASCAT at SNP6 loci with purity/psi fixed to the Battenberg fit.
// alleleCounter runs once per unique BAM so a normal shared by several tumours is counted once.
include { ASCAT_SNP6_REFERENCE; ALLELECOUNTER_SNP6; ASCAT_SNP6_RUN } from '../../modules/local/ascat_snp6'

workflow ASCAT_SNP6 {

    take:
    ch_samples      // [ meta (with purity, ploidy, psi), tumour_bam, tumour_bai, normal_bam, normal_bai ]

    main:
    // SNP6 reference per genome build present in the samplesheet
    ch_ref_input = ch_samples
        .map { meta, _tbam, _tbai, _nbam, _nbai -> meta.build }
        .unique()
        .map { build ->
            def g = params.genomes[build]
            tuple(build,
                  file(g.snp6_positions, checkIfExists: true),
                  g.ascat_g1000_dir ? file(g.ascat_g1000_dir, checkIfExists: true) : [])
        }
    ASCAT_SNP6_REFERENCE(ch_ref_input)
    ch_ref = ASCAT_SNP6_REFERENCE.out.ref       // [ build, ref_dir ]

    // One allele counting job per unique BAM and build
    ch_requests = ch_samples.flatMap { meta, tbam, tbai, nbam, nbai ->
        [
            [ "${meta.build}:${tbam}".toString(), meta, 'tumour', meta.tumour_name, tbam, tbai ],
            [ "${meta.build}:${nbam}".toString(), meta, 'normal', meta.normal_name, nbam, nbai ]
        ]
    }
    ch_count_jobs = ch_requests
        .unique { request -> request[0] }
        .map { key, meta, _role, name, bam, bai -> tuple(meta.build, [key: key, name: name, build: meta.build], bam, bai) }
        .combine(ch_ref, by: 0)
        .map { _build, acmeta, bam, bai, ref -> tuple(acmeta, bam, bai, ref) }
    ALLELECOUNTER_SNP6(ch_count_jobs)

    // Counts back to samples: [ meta, role, counts ]
    ch_counts = ch_requests
        .map { key, meta, role, _name, _bam, _bai -> tuple(key, meta, role) }
        .combine(ALLELECOUNTER_SNP6.out.counts.map { acmeta, counts -> tuple(acmeta.key, counts) }, by: 0)
        .map { _key, meta, role, counts -> tuple(meta, role, counts) }

    ch_tumour_counts = ch_counts.filter { _meta, role, _counts -> role == 'tumour' }.map { meta, _role, counts -> tuple(meta, counts) }
    ch_normal_counts = ch_counts.filter { _meta, role, _counts -> role == 'normal' }.map { meta, _role, counts -> tuple(meta, counts) }

    ch_ascat_input = ch_tumour_counts
        .join(ch_normal_counts)
        .map { meta, tcounts, ncounts -> tuple(meta.build, meta, tcounts, ncounts) }
        .combine(ch_ref, by: 0)
        .map { _build, meta, tcounts, ncounts, ref -> tuple(meta, tcounts, ncounts, ref) }
    ASCAT_SNP6_RUN(ch_ascat_input)

    ch_metrics_table = ASCAT_SNP6_RUN.out.metrics
        .map { _meta, tsv -> tsv }
        .collectFile(name: 'ascat_snp6_metrics.tsv', keepHeader: true, skip: 1, sort: { f -> f.name }, storeDir: "${params.outdir}/summary")

    emit:
    segments      = ASCAT_SNP6_RUN.out.segments
    metrics       = ASCAT_SNP6_RUN.out.metrics
    metrics_table = ch_metrics_table
}
