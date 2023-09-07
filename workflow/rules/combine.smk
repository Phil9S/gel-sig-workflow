rule combine:
    input:
        expand(out_dir+"copy_number_signatures/{cancer}/copy_number_segments.tsv",cancer=cancersUniq)
    output:
        out_dir+"copy_number_signatures/copy_number_segments.tsv"
    singularity:
        image_base_url+"cinsignaturequantification_v1_1_2:latest"
    resources:
        mem_mb=8000,
        threads=1,
        #time="01:00:00" # slurm
        time=60 # lsf
    threads: 1
    shell:
        """
        Rscript --vanilla workflow/scripts/combine_ascat.R {output} {input}
        """
