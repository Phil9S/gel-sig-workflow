// Parameter validation and samplesheet parsing

// True when a samplesheet value is present, i.e. not blank or NA
def isSet(value) {
    if (value == null) {
        return false
    }
    def v = value.toString().trim()
    return !(v in ['', 'NA', 'NaN', 'na', 'null', 'NULL', 'None', '.'])
}

def parseBool(value, fallback) {
    if (!isSet(value)) {
        return fallback
    }
    return value.toString().trim().toUpperCase() in ['TRUE', 'T', 'YES', 'Y', '1']
}

def toNumber(value) {
    if (!isSet(value)) {
        return null
    }
    try {
        return value.toString().trim().toDouble()
    }
    catch (Exception e) {
        return null
    }
}

def normaliseSex(value, String sample) {
    def s = value?.toString()?.trim()?.toUpperCase()
    if (s in ['XX', 'F', 'FEMALE']) {
        return 'XX'
    }
    if (s in ['XY', 'M', 'MALE']) {
        return 'XY'
    }
    error("Sample '${sample}': sex must be XX or XY (female/male also accepted), found '${value}'")
}

// Find the BAM index as <bam>.bai or <name>.bai
def findBamIndex(bam, String sample) {
    def candidates = [
        file("${bam}.bai"),
        file(bam.toString().replaceAll(/\.bam$/, '.bai'))
    ]
    def index = candidates.find { f -> f.exists() }
    if (!index) {
        error("Sample '${sample}': no BAM index found for ${bam} (looked for ${candidates.join(' and ')})")
    }
    return index
}

// swgs-absolutecn v2.0.0 runs one bin size in kb per run, limited to its QDNAseqmod bin annotations
def validateSwgsBin() {
    def allowed = [1, 5, 15, 30, 50, 100, 500, 1000]
    def bin = params.swgs_bin.toString().trim()
    if (!(bin ==~ /\d+/) || !(bin.toInteger() in allowed)) {
        error("--swgs_bin: '${params.swgs_bin}' is not supported; give a single bin size from ${allowed.join(', ')} (swgs-absolutecn v2.0.0 runs one bin size per run)")
    }
}

def validateParameters() {
    def modes = ['copy', 'copyNoFollow', 'link', 'move', 'rellink', 'symlink']
    if (!(params.publish_dir_mode in modes)) {
        error("--publish_dir_mode must be one of ${modes.join(', ')}")
    }
    if (!params.skip_swgs) {
        validateSwgsBin()
    }
    if (params.battenberg_beagle_max_mem.toString().toDouble() <= 0) {
        error("--battenberg_beagle_max_mem must be > 0")
    }
}

// Parse and validate config/samples.tsv into [ meta, tumour_bam, tumour_bai, normal_bam, normal_bai, snv_vcf | [] ].
// Optional columns are snv_vcf (any somatic VCF with AD or AF; Strelka2 is run when absent)
// and swgs_smooth (extra CBS smoothing for noisy sWGS profiles). Other extra columns are ignored.
def readSamplesheet(String path) {
    def sheet = file(path, checkIfExists: true)
    def rows = sheet.splitCsv(header: true, sep: '\t', strip: true)
    if (!rows) {
        error("Samplesheet '${path}' contains no samples")
    }

    def required = ['sample', 'tumour_bam', 'normal_bam', 'tumour_name', 'normal_name', 'sex', 'build', 'cancer']
    def missing = required.findAll { col -> !rows[0].containsKey(col) }
    if (missing) {
        error("Samplesheet '${path}' is missing required column(s): ${missing.join(', ')}")
    }

    def ids = rows.collect { row -> row.sample }
    def duplicated = ids.findAll { id -> ids.count(id) > 1 }.unique()
    if (duplicated) {
        error("Samplesheet '${path}' contains duplicated sample IDs: ${duplicated.join(', ')}")
    }

    def needs_snv_reference = !params.skip_vaf_evaluation && params.call_snvs

    return rows.collect { row ->
        def id = row.sample
        required.each { col ->
            if (!isSet(row[col])) {
                error("Sample '${id}': column '${col}' is empty")
            }
        }

        def build = row.build
        if (!params.genomes.containsKey(build)) {
            error("Sample '${id}': build '${build}' is not configured in params.genomes (available: ${params.genomes.keySet().join(', ')})")
        }
        def genome = params.genomes[build]
        if (!genome.battenberg_ref_dir) {
            error("Sample '${id}': params.genomes.${build}.battenberg_ref_dir is not set (see conf/genomes.config)")
        }
        if (!params.skip_ascat && !genome.snp6_positions) {
            error("Sample '${id}': params.genomes.${build}.snp6_positions is not set; set it or use --skip_ascat")
        }
        def has_vcf = isSet(row.snv_vcf)
        if (needs_snv_reference && !has_vcf && !genome.fasta) {
            error("Sample '${id}': no snv_vcf given and params.genomes.${build}.fasta is not set for Strelka2 calling")
        }

        // keep meta to sample identity only: every key is part of every task hash
        def meta = [
            id          : id,
            tumour_name : row.tumour_name,
            normal_name : row.normal_name,
            sex         : normaliseSex(row.sex, id),
            build       : build,
            cancer      : row.cancer
        ]
        def tumour_bam = file(row.tumour_bam, checkIfExists: true)
        def normal_bam = file(row.normal_bam, checkIfExists: true)
        def snv_vcf = has_vcf ? file(row.snv_vcf, checkIfExists: true) : []

        return [meta, tumour_bam, findBamIndex(tumour_bam, id), normal_bam, findBamIndex(normal_bam, id), snv_vcf]
    }
}

// Per-sample sWGS smoothing flag from the swgs_smooth column, defaulting to --swgs_smooth
def readSwgsSmooth(String path) {
    def rows = file(path, checkIfExists: true).splitCsv(header: true, sep: '\t', strip: true)
    return rows.collectEntries { row -> [(row.sample): parseBool(row.swgs_smooth, params.swgs_smooth)] }
}

// Add Battenberg purity, ploidy and psi to meta.
// Returns null when there is no usable fit so the sample is dropped downstream.
def addPurityPloidy(Map meta, tsv) {
    def rows = tsv.splitCsv(header: true, sep: '\t')
    def row = rows ? rows[0] : [:]
    def purity = toNumber(row.purity)
    def ploidy = toNumber(row.ploidy)
    def psi = toNumber(row.psi)
    if (purity == null || ploidy == null) {
        log.warn("Sample '${meta.id}': no Battenberg purity/ploidy available - skipping sWGS and ASCAT for this sample")
        return null
    }
    return meta + [purity: purity, ploidy: ploidy, psi: psi != null ? psi : ploidy]
}
