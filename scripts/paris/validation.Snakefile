
rule all:
    input:
        "data/qc/self-corr_objf.png"

rule aggregate:
    input:
        "data/output/signature_of_10-genes"
    output:
        "data/output/signatures_z10.tsv"
    shell:
        "cat {input}/signature_* | sort | uniq > {output}"

rule correlations_genes_samples:
    input:
        task="correlations_genes-samples.py",
        ranks="data/inter/paris+ranks.h5an.gz",
        sign="data/output/signatures_z10.tsv"
    output:
        gcorr="data/inter/gcorr.h5an.gz",
        scorr="data/inter/scorr.h5an.gz",
        # plot_gcorr="data/qc/genes_corr.png",
        plot_selfcorr="data/qc/self-corr_objf.png",
    shell:
        # "python3 {input.task} {input.ranks} 10 {output.gcorr} {output.scorr} {output.plot_gcorr} {output.plot_selfcorr} {input.sign}"
        "python3 {input.task} {input.ranks} 10 {output.gcorr} {output.scorr} tmp.png {output.plot_selfcorr} {input.sign}"

