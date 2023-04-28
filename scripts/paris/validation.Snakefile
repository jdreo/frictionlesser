
rule all:
    input:
        "data/qc/self-corr_objf.png",
        "data/output/correlations_signatures-genes.csv",
        "data/output/scores_signatures-samples.csv"

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
    shell:
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
        sgcorr_csv="data/output/correlations_signatures-genes.csv",
    shell:
        "python3 {input.task} {input.ranks} 10 {input.sccorr} {output.sgcorr} {output.sgcorr_csv}"

rule scores_signatures_samples:
    input:
        task="scores_signatures-samples.py",
        sign="data/output/signatures_z10.tsv",
    output:
        ssscores="data/output/scores_signatures-samples.csv"
    shell:
        "python3 {input.task} 10 {output.ssscores} {input.sign}"
