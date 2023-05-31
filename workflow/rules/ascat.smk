rule ascat:
    input:
        unpack(get_sample_bams)
    output:
        out_dir+"fitting/{sample}/{sample}_fit.tsv"
    params:
        sample=get_sample_params,
        outloc=out_dir
    singularity:
        "library://phil9s/gel-pl/ascat_v3_1_2:latest"
    resources:
        mem_mb=8000,
        threads=23,
        time="01:00:00"
    shell:
        """
        Rscript --vanilla workflow/scripts/run_ascat.R {wildcards.sample} \
            {input} {params[sample][tumour_name]} {params[sample][normal_name]} \
            {params[sample][sex]} {params[sample][build]} {output} {params[outloc]} {resources.threads} 
        """

