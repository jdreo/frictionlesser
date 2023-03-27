
SIZE=10

rule all:
    input:
        "signatures-average-correlations.npy",
        "signatures-average-correlations_full.npy"

rule ranksgenes_correlations_npy:
    input:
        ranks="ranks.tsv",
        A="signatures/emil.csv",
        B="signatures/johann.tsv"
    output:
        "genes-ranks-correlations.npy"
    shell:
        "python3 genes-correlations__sign-to-npy.py {input.ranks} {SIZE} {input.A} {input.B}"

rule npycorrsign_signavcorr_npy:
    input:
        ranks="ranks.tsv",
        corrs="genes-ranks-correlations.npy",
        A="signatures/emil.csv",
        B="signatures/johann.tsv"
    output:
        "signatures-average-correlations.npy",
        "signatures-average-correlations_full.npy"
    shell:
        "python3 average-correlations__rkcorr-to-npy.py {input.ranks} {input.corrs} {SIZE} {input.A} {input.B}"
