// ASCAT at SNP6 loci using alleleCounter on tumour/normal WGS BAMs

// SNP6-restricted ASCAT references, once per build.
// Without g1000_dir the 1000 Genomes references in the container under /ref/<build>/ are used.
process ASCAT_SNP6_REFERENCE {
    tag "${build}"
    label 'process_single'
    container "${params.container_ascat}"

    input:
    tuple val(build), path(positions), path(g1000_dir, stageAs: 'g1000_ref')

    output:
    tuple val(build), path("ascat_snp6_${build}"), emit: ref

    script:
    def g = params.genomes[build]
    def base = g1000_dir ? 'g1000_ref/' : "/ref/${build}/"
    def rt = g.ascat_rt_file ? "--rt ${base}${g.ascat_rt_file}" : ''
    """
    prepare_snp6_reference.sh \\
        --positions ${positions} \\
        --loci-prefix ${base}${g.ascat_g1000_loci_prefix} \\
        --alleles-prefix ${base}${g.ascat_g1000_alleles_prefix} \\
        --gc ${base}${g.ascat_gc_file} \\
        ${rt} \\
        --out ascat_snp6_${build}
    """

    stub:
    """
    mkdir -p ascat_snp6_${build}
    for c in \$(seq 1 22) X; do printf '1\\t1\\n' > ascat_snp6_${build}/SNP6_loci_chr\$c.txt; done
    touch ascat_snp6_${build}/GC_SNP6.txt
    """
}

// alleleCounter at SNP6 loci, once per unique BAM
process ALLELECOUNTER_SNP6 {
    tag "${acmeta.name}"
    label 'process_medium'
    container "${params.container_allelecounter}"

    input:
    tuple val(acmeta), path(bam, stageAs: 'input/*'), path(bai, stageAs: 'input/*'), path(ref)

    output:
    tuple val(acmeta), path("alleleFrequencies_chr*.txt"), emit: counts

    script:
    def chr_prefix = params.genomes[acmeta.build].chr_prefix ? 'true' : 'false'
    """
    [ -e "${bam}.bai" ] || ln -s "\$(basename "${bai}")" "${bam}.bai"

    run_allelecounter.sh \\
        --bam ${bam} \\
        --ref-dir ${ref} \\
        --out-prefix alleleFrequencies_chr \\
        --chr-prefix ${chr_prefix} \\
        --min-base-qual ${params.ascat_min_base_qual} \\
        --min-map-qual ${params.ascat_min_map_qual} \\
        --cpus ${task.cpus}
    """

    stub:
    """
    for c in \$(seq 1 22) X; do touch alleleFrequencies_chr\$c.txt; done
    """
}

// ASCAT fit with purity and psi fixed to the Battenberg solution
process ASCAT_SNP6_RUN {
    tag "${meta.id}"
    label 'process_medium'
    container "${params.container_ascat}"

    input:
    tuple val(meta), path(tumour_counts, stageAs: 'tumour_counts/*'), path(normal_counts, stageAs: 'normal_counts/*'), path(ref)

    output:
    tuple val(meta), path("${meta.id}.ascat.segments.txt"),     emit: segments
    tuple val(meta), path("${meta.id}.ascat.segments_raw.txt"), emit: segments_raw
    tuple val(meta), path("${meta.id}.ascat.metrics.tsv"),      emit: metrics
    tuple val(meta), path("${meta.id}.ascat.RData"),            emit: rdata
    tuple val(meta), path("${meta.id}.ascat_unconstrained.segments*.txt"), emit: unconstrained, optional: true
    tuple val(meta), path("*.png"),                             emit: plots, optional: true

    script:
    """
    run_ascat_snp6.R \\
        --sample '${meta.id}' \\
        --tumour_name '${meta.tumour_name}' \\
        --normal_name '${meta.normal_name}' \\
        --sex ${meta.sex} \\
        --build ${meta.build} \\
        --purity ${meta.purity} \\
        --psi ${meta.psi} \\
        --ploidy ${meta.ploidy} \\
        --tumour_counts_prefix tumour_counts/alleleFrequencies_chr \\
        --normal_counts_prefix normal_counts/alleleFrequencies_chr \\
        --alleles_prefix ${ref}/SNP6_alleles_chr \\
        --gc_file ${ref}/GC_SNP6.txt \\
        --rt_file ${ref}/RT_SNP6.txt \\
        --gamma ${params.ascat_gamma} \\
        --penalty ${params.ascat_penalty} \\
        --min_counts ${params.ascat_min_counts} \\
        --seed ${params.ascat_seed} \\
        --run_unconstrained ${params.ascat_run_unconstrained}
    """

    stub:
    """
    printf 'sample\\tchr\\tstartpos\\tendpos\\tnMajor\\tnMinor\\n' > ${meta.id}.ascat.segments.txt
    printf 'sample\\tchr\\tstartpos\\tendpos\\tnMajor\\tnMinor\\tnAraw\\tnBraw\\n' > ${meta.id}.ascat.segments_raw.txt
    printf 'sample\\tbattenberg_purity\\tascat_fixed_goodness_of_fit\\n${meta.id}\\t${meta.purity}\\t99\\n' > ${meta.id}.ascat.metrics.tsv
    touch ${meta.id}.ascat.RData
    """
}
