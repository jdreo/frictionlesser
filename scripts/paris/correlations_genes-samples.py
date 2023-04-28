import sys
import numpy as np
import anndata as ad
from sortedcollections import OrderedSet
import seaborn

import signatures

if __name__ == "__main__":

    assert(len(sys.argv) >= 8)
    franks = sys.argv[1]
    size = int(sys.argv[2])
    if size == 0:
        size = None
    fout_gcorr = sys.argv[3]
    fout_scorr = sys.argv[4]
    fplot_gcorr = sys.argv[5]
    fplot_selfcorr = sys.argv[6]
    fsignatures = sys.argv[7:]

    # Precision for floating-point computation asserts.
    epsilon = 1e-6

    ###########################################################################
    # LOAD
    ###########################################################################

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
    origins = {}
    for fsign in fsignatures:
        s,g,b,S = signatures.load([fsign], size, with_scores=True)
        loaded_signatures[fsign] = s
        genomes[fsign] = g
        breaks[fsign] = b
        loaded_scores[fsign] = S
        for sign in s:
            # origins.setdefault(" ".join(sorted(sign)),set()).add(fsign)
            origins.setdefault(sign,set()).add(fsign)
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

    ###########################################################################
    # GENES CORRELATIONS
    ###########################################################################

    print("Compute cell-gene z-scores sample fractions...", file=sys.stderr, flush=True)
    corr_cellgene = ad.AnnData(
        np.zeros((ncells,len(genome))),
        adata.obs, {"gene":np.asarray(genome)},
        dtype=adata.layers["ranks"].dtype )

    ranks = adata.layers["ranks"]
    # print(ranks.shape, file=sys.stderr, flush=True)
    assert(ranks.shape == adata.X.shape)
    for i,s in enumerate(samples):
        cells = (adata.obs["sample"] == s)
        assert(np.sum(cells) > 1)
        for j,g in enumerate(genome): # FIXME bitand merge truth tables for all genes to avoid this loop
            cranks = ranks[ cells, adata.var["id"] == g ]
            # print("cells ranks:",cranks, file=sys.stderr, flush=True)
            assert(len(cranks.shape) == 1) # only one gene should match.
            sn = cranks.shape[0]
            assert(sn > 0)
            # print("sn:",sn, file=sys.stderr, flush=True)
            smean = np.mean(cranks)
            assert(smean == (sn+1)/2)
            # sdev  = np.std (cranks) # sample mean
            sdev = np.sqrt(np.sum(np.power(cranks-smean,2)))
            if sdev != 0:
                # Standardized z-score across cells,
                # and pre-division for preparing the incoming dot product.
                corr_cellgene[cells,j] = (cranks-smean)/sdev * 1/np.sqrt(nsamples)

                # print("corr:",corr_cellgene[cells,j].X, file=sys.stderr, flush=True)

                # FIXME why the sqrt?
                # ssp2 = np.sqrt(np.sum(np.power(corr_cellgene[cells,j].X,2)))
                ssp2 = np.sum(np.power(corr_cellgene[cells,j].X, 2)) * nsamples
                # print("ssp2:", ssp2, file=sys.stderr, flush=True)
                assert(np.abs(1-ssp2) <= epsilon)
            else:
                corr_cellgene[cells,j] = np.zeros(cranks.shape)
                assert(np.sum(np.power(corr_cellgene[cells,j].X, 2)) * nsamples == 0)

    # Assert sum of scores for each gene.
    for g in genome:
        cgene = corr_cellgene[ :, corr_cellgene.var["gene"] == g ]
        # print(np.sum(cgene.X, axis=1), file=sys.stderr, flush=True)
        # assert(all([0 <= s <= 1 for s in np.sqrt(np.sum(np.power(cgene.X,2), axis=1))]))
        assert(all(0 <= s <= 1 for s in np.sum(np.power(cgene.X,2), axis=1)))

    print("Compute pairwise average correlation matrix for genes over samples...", file=sys.stderr, flush=True)
    corr_cellgene.varp["correlations"] = corr_cellgene.X.T @ corr_cellgene.X
    # print(np.diag(corr_cellgene.varp["correlations"]), file=sys.stderr, flush=True)

    print("Save cell-gene z-score and gene correlations...", file=sys.stderr, flush=True)
    print(corr_cellgene, file=sys.stderr, flush=True)
    corr_cellgene.write(fout_gcorr, compression = "gzip")

    ###########################################################################
    # PLOT
    ###########################################################################



    ###########################################################################
    # SIGNATURES CORRELATIONS
    ###########################################################################

    # for s in unique_signatures:
    #     print("Uniq Sign:",s)
    # for k in origins:
    #     print("Origin:",k,"=>",origins[k])
    var = {"signature":np.asarray([frozenset(s) for s in unique_signatures])}
    aorig = []
    for s in var["signature"]:
        aorig.append(origins[s])
    var["origin"] = np.asarray(aorig)

    print("Compute cell-signatures average correlations over samples...", file=sys.stderr, flush=True)
    corr_cellsign = ad.AnnData(
        np.zeros((ncells,nsignatures)),
        adata.obs,
        var,
        dtype=corr_cellgene.X.dtype )

    for signature in unique_signatures:
        # FIXME annotate with origin file
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
        sngenes = len([i for i in sgenes if i])
        if size:
            assert(sngenes == size)
        corr_cellsign[:, corr_cellsign.var["signature"] == signature] = csum / sngenes

    print("Compute pairwise correlation matrix for signatures...", file=sys.stderr, flush=True)
    # varp is a dictionary holding (var × var) arrays.
    # We want the dot product of all columns.
    smean = np.mean(corr_cellsign.X, axis=1, keepdims=True)
    corr_cellsign.layers["zscore"] = (corr_cellsign.X - smean) / np.sqrt(np.sum(np.power(corr_cellsign.X - smean,2), axis=1, keepdims=True))
    corr_cellsign.varp["correlations"] = corr_cellsign.X.T @ corr_cellsign.X

    print(corr_cellsign, file=sys.stderr, flush=True)
    # print(corr_cellsign.varp["correlations"], file=sys.stderr, flush=True)

    print("Save cell-signatures correlations...", file=sys.stderr, flush=True)
    # AnnData cannot implicitely convert signature's OrderedSet to strings,
    # so we do it manually.
    ssigns = []
    for oset in corr_cellsign.var["signature"]:
        ssigns.append(" ".join(oset)) # Space-separated.
    corr_cellsign.var["signature"] = ssigns
    origs = []
    for orig in corr_cellsign.var["origin"]:
        origs.append(" ".join(orig))
    corr_cellsign.var["origin"] = origs
    print(corr_cellsign, file=sys.stderr, flush=True)
    corr_cellsign.write(fout_scorr, compression = "gzip")

    ###########################################################################
    # PLOT
    ###########################################################################

    print("Plot self-correlations of signatures against scores...", file=sys.stderr, flush=True)
    selfcorr = np.abs(np.diag(corr_cellsign.varp["correlations"]))
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
    fig.savefig(fplot_selfcorr, dpi=600)

    print("Done", file=sys.stderr, flush=True)
