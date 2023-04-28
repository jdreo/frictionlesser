
rule all:
    input:
        "data/qc/clustermap_correlations.png",
        # "data/qc/clustermap_pvalues.png",
        "data/qc/observed_genes.txt",
        "data/qc/observed_genes_distribution.png",
        "data/qc/selfcorr-score_local.png",
        "data/qc/jaccard-distances.png"

rule aggregate:
    input:
        "data/output/signature_of_10-genes"
    output:
        "data/output/signatures_z10.tsv"
    shell:
        "cat {input}/signature_* | sort | uniq > {output}"

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
# Two signatures sets.
    input:
        task="genes-standardized-over-samples__sign-to-npy.py",
        ranks="data/inter/ranks.tsv",
        A="data/output/emil_all.csv",
        B="data/output/signatures_z10.tsv"
    output:
        scores="data/inter/genes-standardized-scores.npy",
        indices="data/inter/genes-standardized-scores_genes.csv"
    shell:
        "python3 {input.task} {input.ranks} 10 {output.scores} {output.indices} {input.A} {input.B}"

rule genes_standardized_scores_local:
# Single signature set.
    input:
        task="genes-standardized-over-samples__sign-to-npy.py",
        ranks="data/inter/ranks.tsv",
        signs="data/output/signatures_z10.tsv"
    output:
        scores="data/inter/genes-standardized-scores_local.npy",
        indices="data/inter/genes-standardized-scores_genes_local.csv"
    shell:
        "python3 {input.task} {input.ranks} 10 {output.scores} {output.indices} {input.signs}"

rule qc_selfcorr_score_local:
    input:
        task="qc_selfcorr-score.py",
        ranks="data/inter/ranks.tsv",
        scores="data/inter/genes-standardized-scores_local.npy",
        indices="data/inter/genes-standardized-scores_genes_local.csv",
        signs="data/output/signatures_z10.tsv"
    output:
        corr="data/inter/signatures-average-selfcorrelations_local.npy",
        plot="data/qc/selfcorr-score_local.png"
    shell:
        "python3 {input.task} {input.ranks} {input.scores} {input.indices} 10 {output.corr} {output.plot} {input.signs}"

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
        plot="data/qc/clustermap_correlations.png",
        membership="data/inter/cluster_membership.csv"
    shell:
        "python3 {input.task} {input.corr} {input.orig} {input.pval} 0.01 1.5 {output.plot} {output.membership}"

rule qc_jaccard:
    input:
        task="signatures_jaccard-distances_plot.py",
        A="data/output/emil_all.csv",
        B="data/output/signatures_z10.tsv"
    output:
        "data/qc/jaccard-distances.png"
    shell:
        "python3 {input.task} 10 0.5 {output} {input.A} {input.B}"
