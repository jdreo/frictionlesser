
SIZE=10

rule all:
    input:
        "data/inter/signatures-average-correlations.npy",
        "data/inter/signatures-average-correlations_full.npy"

rule aggregate:
    input:
        "data/output/signature_of_{SIZE}-genes/"
    output:
        "data/output/signatures_z{SIZE}.tsv"
    shell:
        "cat data/output/signature_of_{SIZE}-genes/signature_* > {output}"

rule ranksgenes_correlations_npy:
    input:
        ranks="data/inter/ranks.tsv",
        A="data/output/emil_all.csv",
        B="data/output/signatures_z{SIZE}.tsv"
    output:
        "data/inter/genes-ranks-correlations.npy"
    shell:
        "python3 genes-correlations__sign-to-npy.py {input.ranks} {SIZE} {input.A} {input.B}"

rule npycorrsign_signavcorr_npy:
    input:
        ranks="data/inter/ranks.tsv",
        corrs="data/inter/genes-ranks-correlations.npy",
        A="data/output/emil_all.csv",
        B="data/output/signatures_z{SIZE}.tsv"
    output:
        "data/inter/signatures-average-correlations.npy",
        "data/inter/signatures-average-correlations_full.npy"
    shell:
        "python3 average-correlations__rkcorr-to-npy.py {input.ranks} {input.corrs} {SIZE} {input.A} {input.B}"
