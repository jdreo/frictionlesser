
rule all:
    input:
        "clustermap_correlations.png",
        "clustermap_pvalues.png",
        "data/qc/observed_genes.txt",
        "data/qc/observed_genes_distribution.png"

rule aggregate:
    input:
        "data/output/signature_of_10-genes"
    output:
        "data/output/signatures_z10.tsv"
    shell:
        "cat {input}/signature_* > {output}"

rule qc_observed_genes:
    input:
        task="qc-observed-genes.py",
        ranks="data/inter/ranks.tsv",
        A="data/output/emil_all.csv",
        B="data/output/signatures_z10.tsv"
    output:
        text="data/qc/observed_genes.txt",
        plot="data/qc/observed_genes_distribution.png"
    shell:
        "python3 {input.task} 10 {input.ranks} {output.plot} {input.A} {input.B} > {output.text}"

rule genes_standardized_scores:
    input:
        task="genes-standardized-over-samples__sign-to-npy.py",
        ranks="data/inter/ranks.tsv",
        A="data/output/emil_all.csv",
        B="data/output/signatures_z10.tsv"
    output:
        scores="data/inter/genes-standardized-scores.npy",
        indices="data/inter/genes-standardized-scores_genes.csv"
    shell:
        "python3 {input.task} {input.ranks} 10 {input.A} {input.B} {output.scores} {output.indices}"

rule signatures_correlations:
    input:
        task="signatures-correlations__sign-to-npy.py",
        ranks="data/inter/ranks.tsv",
        scores="data/inter/genes-standardized-scores.npy",
        indices="data/inter/genes-standardized-scores_genes.csv",
        A="data/output/emil_all.csv",
        B="data/output/signatures_z10.tsv"
    output:
        corr="data/inter/signatures-average-correlations.npy",
        orig="data/inter/signatures-average-correlations_origins.csv"
    shell:
        "python3 {input.task} {input.ranks} {input.scores} {input.indices} 10 {output.corr} {output.orig} {input.A} {input.B}"

rule signatures_pvalues:
    input:
        task="signatures-pvalues.py",
        ranks="data/inter/ranks.tsv",
        corr="data/inter/signatures-average-correlations.npy",
    output:
        "data/inter/signatures-pvalues.npy",
    shell:
        "python3 {input.task} {input.ranks} {input.corr} {output}"

rule clustering_corr:
    input:
        task="signatures-clustering.py",
        corr="data/inter/signatures-average-correlations.npy",
        pval="data/inter/signatures-pvalues.npy",
        orig="data/inter/signatures-average-correlations_origins.npy"
    output:
        plot="clustermap_correlations.png",
        membership="data/inter/cluster_membership.csv"
    shell:
        "python3 {input.task} {input.corr} {input.orig} {input.pval} 0.01 1.5 {output.plot} {output.membership}"

