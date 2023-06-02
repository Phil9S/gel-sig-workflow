rule get_signature_activities:
    input:
        out_dir+"copy_number_signatures/copy_number_segments.tsv"
    output:
        out_dir+"copy_number_signatures/cn_signatures.rds"
    params:
        genome=genome,
        method=method
    singularity:
        image_base_url+"cinsignaturequantification_v1_1_2:latest"
    resources:
        mem_mb=4000,
        threads=1,
        time="01:00:00"
    shell:
        """
        Rscript --vanilla workflow/scripts/run_signature_quantification.R {input} {output} {params.genome} {params.method} {resources.threads}
        """
