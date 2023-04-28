
rule all:
    input:
        "data/qc/self-corr_objf.png",
        "data/output/sgcorr.csv"

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
        gccorr="data/inter/gccorr.h5an.gz",
        # sccorr="data/inter/sccorr.h5an.gz",
        # sgcorr="data/inter/sgcorr.h5an.gz",
        # plot_gcorr="data/qc/genes_corr.png",
        # plot_selfcorr="data/qc/self-corr_objf.png",
    shell:
        # "python3 {input.task} {input.ranks} 10 {output.gccorr} {output.sccorr} {output.plot_gcorr} {output.plot_selfcorr} {input.sign}"
        # "python3 {input.task} {input.ranks} 10 {output.gccorr} {output.sccorr} {output.sgcorr} tmp.png {output.plot_selfcorr} {input.sign}"
        "python3 {input.task} {input.ranks} 10 {output.gccorr} {input.sign}"

rule correlations_signatures:
    input:
        task="correlations_signatures.py",
        ranks="data/inter/paris+ranks.h5an.gz",
        gccorr="data/inter/gccorr.h5an.gz",
        sign="data/output/signatures_z10.tsv"
    output:
        sccorr="data/inter/sccorr.h5an.gz",
    shell:
        "python3 {input.task} {input.ranks} 10 {input.gccorr} {output.sccorr} {input.sign}"

rule qc_selfcorr_score:
    input:
        task="qc_selfcorr-score.py",
        sccorr="data/inter/sccorr.h5an.gz",
        sign="data/output/signatures_z10.tsv"
    output:
        plot_selfcorr="data/qc/self-corr_objf.png",
    shell:
        "python3 {input.task} {input.sccorr} 10 {output.plot_selfcorr} {input.sign}"

rule correlations_signatures_genes:
    input:
        task="correlations_signatures-genes.py",
        ranks="data/inter/paris+ranks.h5an.gz",
        sccorr="data/inter/sccorr.h5an.gz",
    output:
        sgcorr="data/inter/sgcorr.h5an.gz",
        sgcorr_csv="data/output/sgcorr.csv",
    shell:
        "python3 {input.task} {input.ranks} 10 {input.sccorr} {output.sgcorr} {output.sgcorr_csv}"
