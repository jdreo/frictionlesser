import sys
import numpy as np
import anndata as ad

import signatures

if __name__ == "__main__":

    assert(len(sys.argv) >= 6)
    franks = sys.argv[1]
    size = int(sys.argv[2])
    if size == 0:
        size = None
    fgccorr = sys.argv[3]
    fout_sccorr = sys.argv[4]
    fsignatures = sys.argv[5:]

    # Precision for floating-point computation asserts.
    epsilon = 1e-6

    print("Load signatures from: ", fsignatures, file=sys.stderr, flush=True)
    unique_signatures, unique_scores, origins = signatures.load_unique_signatures(fsignatures, size)
    nsignatures = len(unique_signatures)
    print("Found",nsignatures,"unique signatures", file=sys.stderr, flush=True)

    print("Load annotated ranks data from: ", franks, file=sys.stderr, flush=True)
    cells_allgenes = ad.read(franks)
    print(cells_allgenes, file=sys.stderr, flush=True)
    ncells = cells_allgenes.shape[0]
    ngenes_all = cells_allgenes.shape[1]

    print("Load z-scores of genes in signatures: ", fgccorr, file=sys.stderr, flush=True)
    cells_sgenes = ad.read(fgccorr)
    print(cells_sgenes, file=sys.stderr, flush=True)
    ngenes = cells_sgenes.shape[1]

    ###########################################################################
    # SIGNATURES CORRELATIONS
    # (Average correlation across genes in each signatures)
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
    cells_signs = ad.AnnData(
        np.zeros((ncells,nsignatures)),
        cells_allgenes.obs,
        var,
        dtype=cells_allgenes.X.dtype )
    cells_signs.varp["av.correlations"] = np.zeros((nsignatures,nsignatures))

    # cells_signs structure:
    #
    #   ← signatures →
    # ┌─────────────────┐
    # │ var:"signature" │
    # │     "origin"    │
    # └─────────────────┘
    # ┏━━━━━━━━━━━━━━━━━┓  ┌────────────────┐
    # ┃ X:              ┃  │ obs:           │  ↑
    # ┃sum.zscore/ngenes┃  │ "sample"       │ cells
    # ┃                 ┃  │                │  ↓
    # ┗━━━━━━━━━━━━━━━━━┛  └────────────────┘
    # ┌─────────────────┐
    # │ varp:           │
    # │ "correlations"  │  ↑
    # │                 │ signatures
    # │                 │  ↓
    # │                 │
    # └─────────────────┘
    #   ← signatures →
    #

    def genes_of(s):
        s_genes = np.zeros(ngenes, dtype=bool)
        for g in s:
            # Logical "and" on locations of each genes.
            s_genes |= (cells_sgenes.var["id"] == g)
        return s_genes

    print("Compute pairwise correlation matrix for signatures...", file=sys.stderr, flush=True)

    # for s0 in unique_signatures:
    #     s0_genes = genes_of(s0)
    #     s0_size = len(s0)
    #     c0_sum = np.sum(cells_sgenes.varp["correlations"][s0_genes,:], axis=0)
    #     is0 = cells_signs.var["signature"] == s0
    #     for s1 in unique_signatures:
    #         s1_genes = genes_of(s1)
    #         s1_size = len(s1)
    #         c_sum = np.sum(c0_sum[s1_genes]) # on axis=1
    #         is1 = cells_signs.var["signature"] == s0
    #         size2 = s0_size * s1_size
    #         corr = c_sum / size2
    #         assert(-1-epsilon <= corr <= 1+epsilon)
    #         cells_signs.varp["av.correlations"][is1,is0] = corr
    #         cells_signs.varp["av.correlations"][is0,is1] = corr

    # assert(all(-1-epsilon <= c <= 1+epsilon for c in np.nditer(cells_signs.varp["av.correlations"])))


    for signature in unique_signatures:
        sgenes = genes_of(signature)

        # Correlations for those genes.
        # Use the data array to avoid useless filtering.
        # zscores = cells_sgenes.X[:, sgenes]
        zscores = cells_sgenes.layers["zscore"][:, sgenes]

        # Keepdims avoid collapsing dimension.
        zsum = np.sum(zscores, axis=1, keepdims=True)
        assert(zsum.shape == (ncells,1))

        # Average z-score.
        sngenes = np.sum(sgenes)
        if size:
            assert(sngenes == size)
        cells_signs[:, cells_signs.var["signature"] == signature] = zsum / sngenes # Squared nb of genes at the end.

    # varp is a dictionary holding (var × var) arrays.
    # We want the dot product of all columns.
    # Average of z-scores is already standardized.
    cells_signs.varp["correlations"] = cells_signs.X.T @ cells_signs.X
    assert(all(-1-epsilon <= c <= 1+epsilon for c in np.nditer(cells_signs.varp["correlations"])))

    print(cells_signs, file=sys.stderr, flush=True)
    # print(cells_signs.varp["correlations"], file=sys.stderr, flush=True)

    print("Save cell-signatures correlations...", file=sys.stderr, flush=True)
    # AnnData cannot implicitely convert signature's OrderedSet to strings,
    # so we do it manually.
    ssigns = []
    for oset in cells_signs.var["signature"]:
        ssigns.append(" ".join(oset)) # Space-separated.
    cells_signs.var["signature"] = ssigns
    origs = []
    for orig in cells_signs.var["origin"]:
        origs.append(" ".join(orig))
    cells_signs.var["origin"] = origs
    print(cells_signs, file=sys.stderr, flush=True)
    cells_signs.write(fout_sccorr, compression = "gzip")

