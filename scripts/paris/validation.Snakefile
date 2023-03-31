
SIZE=10

rule all:
    input:
        "clustermap.png"

# rule aggregate:
#     input:
#         "data/output/signature_of_{SIZE}-genes/"
#     output:
#         "data/output/signatures_z{SIZE}.tsv"
#     shell:
#         "cat data/output/signature_of_{SIZE}-genes/signature_* > {output}"

# rule ranksgenes_correlations_npy:
#     input:
#         ranks="data/inter/ranks.tsv",
#         A="data/output/emil_all.csv",
#         B="data/output/signatures_z{SIZE}.tsv"
#     output:
#         "data/inter/genes-ranks-correlations.npy"
#     shell:
#         "python3 genes-correlations__sign-to-npy.py {input.ranks} {SIZE} {input.A} {input.B}"

# rule npycorrsign_signavcorr_npy:
#     input:
#         ranks="data/inter/ranks.tsv",
#         corrs="data/inter/genes-ranks-correlations.npy",
#         A="data/output/emil_all.csv",
#         B="data/output/signatures_z{SIZE}.tsv"
#     output:
#         "data/inter/signatures-average-correlations.npy",
#         "data/inter/signatures-average-correlations_full.npy"
#     shell:
#         "python3 average-correlations__rkcorr-to-npy.py {input.ranks} {input.corrs} {SIZE} {input.A} {input.B}"

rule genes_standardized_scores:
    input:
        task="genes-standardized-over-samples__sign-to-npy.py",
        ranks="data/inter/ranks.tsv",
        A="data/output/emil_all.csv",
        B="data/output/signatures_z{SIZE}.tsv"
    output:
        scores="data/inter/genes-standardized-scores.npy",
        indices="data/inter/genes-standardized-scores_genes.csv"
    shell:
        "python3 {input.task} {input.ranks} {SIZE} {input.A} {input.B} {output.scores} {output.indices}"

rule signatures_correlations:
    input:
        task="signatures-sum-standarized-scores__sign-to-npy.py",
        scores="data/inter/genes-standardized-scores.npy",
        indices="data/inter/genes-standardized-scores_genes.csv",
        A="data/output/emil_all.csv",
        B="data/output/signatures_z{SIZE}.tsv"
    output:
        corr="data/inter/signatures-average-correlations.npy",
        orig="data/inter/signatures-average-correlations_origins.npy"
    shell:
        "python3 {input.task} {input.scores} {input.indices} {SIZE} {input.A} {input.B} {output.corr} {output.orig}"
