rule combine:
    input:
        expand(out_dir+"fitting/{sample}/{sample}_fit.tsv",sample=samples.index)
    output:
        out_dir+"test_out.txt"
    shell:
        """
        touch {output}
        """
