rule combine:
    input:
        expand("{sample}_fit.tsv",sample=samples.index)
    output:
        "test_out.txt"
    shell:
        """
        touch test_out.txt
        """
