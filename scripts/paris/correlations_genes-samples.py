import sys
import numpy as np
import anndata as ad
from sortedcollections import OrderedSet
import seaborn

import signatures

if __name__ == "__main__":

    assert(len(sys.argv) >= 6)
    franks = sys.argv[1]
    size = int(sys.argv[2])
    if size == 0:
        size = None
    fout_gcorr = sys.argv[3]
    fout_scorr = sys.argv[4]
    fplot = sys.argv[5]
    fsignatures = sys.argv[6:]

    print("Load annotated ranks data from: ", franks, file=sys.stderr, flush=True)
    adata = ad.read(franks)
    print(adata, file=sys.stderr, flush=True)
    ncells = adata.shape[0]
    ngenes_all = adata.shape[1]
    print("Loaded", ncells,"cells and", ngenes_all,"genes", file=sys.stderr, flush=True)

    print("Load signatures from: ", fsignatures, file=sys.stderr, flush=True)
    # genome is the genes actually used in at least one signature.
    genome = signatures.load_genome(fsignatures, filter_size = size)
    ngenes = len(genome)
    print("Found",ngenes,"unique genes in signatures of size",size, file=sys.stderr, flush=True)
    assert(len(genome) > 0)

    loaded_signatures = {}
    genomes = {}
    breaks = {}
    unique_signatures = OrderedSet()
    loaded_scores = {}
    unique_scores = {}
    for fsign in fsignatures:
        s,g,b,S = signatures.load([fsign], size, with_scores=True)
        loaded_signatures[fsign] = s
        genomes[fsign] = g
        breaks[fsign] = b
        loaded_scores[fsign] = S
        unique_signatures |= s
        for sk in S:
            if __debug__:
                if sk in unique_scores:
                    assert(unique_scores[sk] == S[sk])
            unique_scores[sk] = S[sk]
    nsignatures = len(unique_signatures)
    print("Found",nsignatures,"unique signatures", file=sys.stderr, flush=True)

    print("Gather samples...", file=sys.stderr, flush=True)
    samples = OrderedSet()
    for s in adata.obs["sample"]:
        samples.add(s)
    print("Found",len(samples),"samples", file=sys.stderr, flush=True)
    nsamples = len(samples)
    assert(len(samples) > 0)

    print("Compute genes correlations...", file=sys.stderr, flush=True)
    corr_cellgene = ad.AnnData(
        np.zeros((ncells,len(genome))),
        adata.obs, {"gene":np.asarray(genome)},
        dtype=adata.layers["ranks"].dtype )

    ranks = adata.layers["ranks"]
    # print(ranks.shape, file=sys.stderr, flush=True)
    assert(ranks.shape == adata.X.shape)
    for i,s in enumerate(samples):
        cells = adata.obs["sample"] == s
        for j,g in enumerate(genome):
            cranks = ranks[ cells, adata.var["id"] == g ]
            nc = cranks.shape[0]
            smean = np.mean(cranks)
            sdev  = np.std (cranks)
            if sdev != 0:
                corr_cellgene[cells,j] = np.clip(
                        (cranks-smean)/sdev * 1/np.sqrt(nc-1) * 1/nsamples,
                    -1, 1)
            else:
                corr_cellgene[cells,j] = np.zeros(cranks.shape)

    print("Save genes correlations...", file=sys.stderr, flush=True)
    print(corr_cellgene, file=sys.stderr, flush=True)
    corr_cellgene.write(fout_gcorr, compression = "gzip")

    print("Compute signatures correlations...", file=sys.stderr, flush=True)
    corr_cellsign = ad.AnnData(
        np.zeros((ncells,nsignatures)),
        adata.obs,
        {"signature":np.asarray(unique_signatures)},
        dtype=corr_cellgene.X.dtype )

    for signature in unique_signatures:
        sgenes = np.zeros(ngenes, dtype=bool)
        for g in signature:
            # Logical "and" on locations of each genes.
            sgenes |= (corr_cellgene.var["gene"] == g)

        # Correlations for those genes.
        # Use the data array to avoid useless filtering.
        scorrs = corr_cellgene.X[:, sgenes]

        # Keepdims avoid collapsing dimension size == 1.
        csum = np.sum(scorrs, axis=1, keepdims=True)

        # Sum across genes.
        corr_cellsign[:, corr_cellsign.var["signature"] == signature] = csum

    # varp is a dictionary holding (var × var) arrays.
    # We want the dot product of all columns.
    corr_cellsign.varp["correlations"] = corr_cellsign.X.T @ corr_cellsign.X

    print(corr_cellsign, file=sys.stderr, flush=True)
    print(corr_cellsign.varp["correlations"], file=sys.stderr, flush=True)

    print("Save signatures correlations...", file=sys.stderr, flush=True)
    # AnnData cannot implicitely convert signature's OrderedSet to strings,
    # so we do it manually.
    ssigns = []
    for oset in corr_cellsign.var["signature"]:
        ssigns.append(" ".join(oset)) # Space-separated.
    corr_cellsign.var["signature"] = ssigns
    print(corr_cellsign, file=sys.stderr, flush=True)
    corr_cellsign.write(fout_scorr, compression = "gzip")

    print("Plot self-correlations of signatures against scores...", file=sys.stderr, flush=True)
    selfcorr = np.diag(corr_cellsign.X)
    if __debug__:
        for c in selfcorr:
            assert(-1 <= c <= 1)

    objf = [unique_scores[frozenset(s.split())] for s in corr_cellsign.var["signature"]]
    if __debug__:
        for s in objf:
            assert(s is not None)
            assert(s > 0)

    plot = seaborn.scatterplot(x=objf, y=selfcorr, s=5)
    plot.set(xlabel="objective function", ylabel="self-correlation",
             title="Relationship between\nobjective function and self-correlation of genes, for signatures in\n{}".format(fsignatures))
    fig = plot.get_figure()
    fig.savefig(fplot, dpi=600)

    print("Done", file=sys.stderr, flush=True)
