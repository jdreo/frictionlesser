
# FIXME how to do that from wildcards?
configfile: "config_validation.yaml"

rule all:
    input:
        "data/qc/self-corr_objf.png",
        "data/output/correlations_signatures-genes.csv",
        "data/output/scores_signatures-samples.csv",
        "data/qc/gcpval_clustering.png",
        "data/qc/gcpval_membership.csv",

rule aggregate:
    input:
        "data/output/signature_of_{params.size}-genes",
    output:
        "data/output/signatures_z{params.size}.tsv",
    shell:
        "cat {input}/signature_* | sort | uniq > {output}"

rule correlations_genes_samples:
    params:
        size=lambda wildcards: config["size"]
    input:
        task="correlations_genes-samples.py",
        ranks="data/inter/paris+ranks.h5an.gz",
        sign="data/output/signatures_z{params.size}.tsv",
    output:
        gccorr="data/inter/gccorr.h5an.gz",
    shell:
        "python3 {input.task} {input.ranks} {params.size} {output.gccorr} {input.sign}"

rule correlations_signatures:
    params:
        size=lambda wildcards: config["size"]
    input:
        task="correlations_signatures.py",
        ranks="data/inter/paris+ranks.h5an.gz",
        gccorr="data/inter/gccorr.h5an.gz",
        sign="data/output/signatures_z{params.size}.tsv",
    output:
        sccorr="data/inter/sccorr.h5an.gz",
    shell:
        "python3 {input.task} {input.ranks} {params.size} {input.gccorr} {output.sccorr} {input.sign}"

rule qc_selfcorr_score:
    params:
        size=lambda wildcards: config["size"]
    input:
        task="qc_selfcorr-score.py",
        sccorr="data/inter/sccorr.h5an.gz",
        sign="data/output/signatures_z{params.size}.tsv",
    output:
        plot_selfcorr="data/qc/self-corr_objf.png",
    shell:
        "python3 {input.task} {input.sccorr} {params.size} {output.plot_selfcorr} {input.sign}"

rule correlations_signatures_genes:
    params:
        size=lambda wildcards: config["size"]
    input:
        task="correlations_signatures-genes.py",
        ranks="data/inter/paris+ranks.h5an.gz",
        sccorr="data/inter/sccorr.h5an.gz",
    output:
        sgcorr="data/inter/sgcorr.h5an.gz",
        sgcorr_csv="data/output/correlations_signatures-genes.csv",
    shell:
        "python3 {input.task} {input.ranks} {params.size} {input.sccorr} {output.sgcorr} {output.sgcorr_csv}"

rule scores_signatures_samples:
    params:
        size=lambda wildcards: config["size"]
    input:
        task="scores_signatures-samples.py",
        sign="data/output/signatures_z{params.size}.tsv",
    output:
        ssscores="data/output/scores_signatures-samples.csv",
    shell:
        "python3 {input.task} {params.size} {output.ssscores} {input.sign}"

rule qc_genes_corr_clustering:
    input:
        task="qc_genes-corr_clustering.py",
        gccorr="data/inter/gccorr.h5an.gz",
    output:
        plot_cluster="data/qc/gcpval_clustering.png",
        gcpval_members="data/qc/gcpval_membership.csv",
    shell:
        "python3 {input.task} {input.gccorr} {output.plot_cluster} {output.gcpval_members}"

rule qc_observed_genes:
    params:
        size=lambda wildcards: config["size"]
    input:
        task="qc_observed-genes.py",
        ranks="data/inter/paris+ranks.h5an.gz",
        sign="data/output/signatures_z{params.size}.tsv",
    output:
        plot="qc_observed-genes"
    shell:
        "python3 {params.size} {input.ranks} {output.plot} {input.sign}"
