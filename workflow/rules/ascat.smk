rule ascat:
    input:
        unpack(get_bams)
    output:
        fit="{sample}_fit.tsv"
    singularity:
        "library://phil9s/gel-pl/ascat_v3_1_2:latest"
    resources:
        mem_mb=8000,
        threads=8,
        time="02:00:00"
    script:
        "scripts/test_sing.R"

