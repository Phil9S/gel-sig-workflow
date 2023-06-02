rule combine:
    input:
        expand(out_dir+"fitting/{sample}/{sample}.segments_raw.txt",sample=samples.index)
    output:
        out_dir+"copy_number_signatures/copy_number_segments.tsv"
    singularity:
        image_base_url+"cinsignaturequantification_v1_1_2:latest"
    resources:
        mem_mb=8000,
        threads=1,
        time="01:00:00"
    shell:
        """
        Rscript --vanilla workflow/scripts/combine_ascat.R {output} {input}
        """
