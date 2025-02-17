rule get_signature_activities_cancer:
    input:
        out_dir+"copy_number_signatures/{cancer}/copy_number_segments.tsv"
    output:
        out_dir+"copy_number_signatures/{cancer}/cn_signatures.rds"
    params:
        genome=genome,
        method=method
    singularity:
        image_base_url+"cinsignaturequantification:latest"
    resources:
        mem_mb=4000,
        threads=1,
        #time="01:00:00" # slurm
        time=60 # lsf
    threads: 1
    shell:
        """
        Rscript --vanilla workflow/scripts/run_signature_quantification.R {input} {output} {params.genome} {params.method} {resources.threads}
        """
