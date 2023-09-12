def get_cancer_sets(cancersUniq):
    sdict = {}
    for cu in cancersUniq:
        samples_sub = samples.loc[samples["cancer"].isin([cu]), "sample"]
        sdict[cu] = samples_sub.tolist()
    return sdict

cancer_subsets = get_cancer_sets(cancersUniq)
#print(cancer_subsets)

rule combine_cancer:
    input:
        lambda wildcards: \
            [out_dir+"fitting/{{cancer}}/{0}/{0}.segments_raw.txt".format(samplesub) \
                for samplesub in cancer_subsets[wildcards.cancer]\
            ]
    output:
        out_dir+"copy_number_signatures/{cancer}/copy_number_segments.tsv"
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
        Rscript --vanilla workflow/scripts/combine_ascat_cancer.R {output} {input}
        """
