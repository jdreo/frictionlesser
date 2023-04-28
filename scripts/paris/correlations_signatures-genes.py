import sys
import numpy as np
import anndata as ad

import signatures

if __name__ == "__main__":

    assert(len(sys.argv) == 6)
    franks = sys.argv[1]
    size = int(sys.argv[2])
    if size == 0:
        size = None
    fsccorr = sys.argv[3]
    fout_sgcorr = sys.argv[4]
    fout_sgcorr_csv = sys.argv[5]

    # Precision for floating-point computation asserts.
    epsilon = 1e-6

    print("Load annotated cells-signatures data from: ", fsccorr, file=sys.stderr, flush=True)
    cells_signs = ad.read(fsccorr)
    print(cells_signs, file=sys.stderr, flush=True)
    ncells = cells_signs.shape[0]

    print("Load annotated ranks data from: ", franks, file=sys.stderr, flush=True)
    cells_allgenes = ad.read(franks)
    print(cells_allgenes, file=sys.stderr, flush=True)
    ncells = cells_allgenes.shape[0]
    ngenes_all = cells_allgenes.shape[1]

    ###########################################################################
    # SIGNATURES-GENES CORRELATIONS
    ###########################################################################

    print("Compute signatures-genes correlations...", file=sys.stderr, flush=True)

    # print("Standardization...", file=sys.stderr, flush=True)
    # "keepdims" option not suported by sparse matrices, but is still the default behavior.
    # rmean = np.mean(cells_allgenes.X, axis=1) #, keepdims=True)
    # FIXME move the zscore layer computation in the cells_allgenes section and save.
    # cells_allgenes.layers["zscore"] = (cells_allgenes.X - rmean) / np.sqrt(np.sum(np.power(cells_allgenes.X - rmean, 2), axis=1)) #, keepdims=True))

    print("Big matrix multiplication...", file=sys.stderr, flush=True)
    signs_allgenes = ad.AnnData(
        # Correlation between signatures z-scores vector (i.e. standardized average ranks vector),
        # and genes z-scores vector (i.e. standardized ranks vector).
        #                 (sign × cells) @ (cells × genes) = (sign × genes)
        cells_signs.X.T @ cells_allgenes.layers["zscore"], #np.zeros((nsignatures,ngenes_all)),
        cells_signs.var, # Metadata about signatures
        cells_allgenes.var,         # Metadata about all genes
        dtype=cells_signs.X.dtype )

    assert(all(-1-epsilon <= c <= 1+epsilon for c in np.nditer(signs_allgenes.X)))

    # signs_allgenes structure:
    #
    #   ← all genes  →
    # ┌───────────────┐
    # │ var:"id", […] │
    # └───────────────┘
    # ┏━━━━━━━━━━━━━━━┓  ┌────────────────┐
    # ┃ X:            ┃  │ obs:           │  ↑
    # ┃ correlations  ┃  │ "signature",   │ signatures
    # ┃               ┃  │ "origin"       │  ↓
    # ┗━━━━━━━━━━━━━━━┛  └────────────────┘
    #   ← all genes  →

    print("Save signatures-genes correlations...", file=sys.stderr, flush=True)
    signs_allgenes.write(fout_sgcorr, compression="gzip")
    print(signs_allgenes, file=sys.stderr, flush=True)

    print("Export signatures-genes correlations to CSV...", file=sys.stderr, flush=True)
    # Set names.
    signs_allgenes.obs_names = signs_allgenes.obs["signature"]
    signs_allgenes.var_names = signs_allgenes.var["id"]

    sdf = signs_allgenes.to_df()
    sdf.sort_values(by=["signature"], ascending=True, inplace=True)
    sdf.to_csv(fout_sgcorr_csv)

