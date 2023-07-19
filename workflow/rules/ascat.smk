rule ascat:
    input:
        unpack(get_sample_bams)
    output:
        out_dir+"fitting/{sample}/{sample}.segments_raw.txt"
    params:
        sample=get_sample_params,
        outloc=out_dir
    singularity:
        image_base_url+"ascat_v3_1_2:latest"
    resources:
        mem_mb=8000,
        threads=23,
        time="01:00:00" # slurm
        #time=60 # lsf
    threads: 23
    shell:
        """
        Rscript --vanilla workflow/scripts/run_ascat.R {wildcards.sample} \
            {input} {params[sample][tumour_name]} {params[sample][normal_name]} \
            {params[sample][sex]} {params[sample][build]} {output} {params[outloc]} {resources.threads}

        ln -T -f {params[outloc]}fitting/{wildcards.sample}/{params[sample][tumour_name]}.segments_raw.txt \
    {params[outloc]}fitting/{wildcards.sample}/{wildcards.sample}.segments_raw.txt 
        """

