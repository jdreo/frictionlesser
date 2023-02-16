import sys
from scipy.sparse import load_npz
import scanpy
import pandas as pd
import anndata as ad
# from sklearn.preprocessing import normalize
# from itertools import zip_longest
# import numpy
# import scipy

def check(counts, meta, genes):
    assert(counts.shape[0] > 0)
    assert(counts.shape[1] > 0)
    assert(counts.shape[0] == len(meta))
    assert(counts.shape[1] == len(genes))
    print("OK", file=sys.stderr, flush=True)

if __name__ == "__main__":

    if len(sys.argv) == 4:

        fcounts = sys.argv[1]
        ffeatures = sys.argv[2]
        fmeta = sys.argv[3]


        print("Load data files...", file=sys.stderr, flush=True)
        # Load data.
        with open(fcounts) as fd:
            # This is a numpy sparse matrix: cells × genes.
            counts = load_npz(fcounts)

        # Those are pandas' dataframes.
        meta  = pd.read_csv(fmeta)
        genes = pd.read_csv(ffeatures)

        print("Loaded",len(genes["id"]),"genes and",len(meta["cell_types_v2"]),"cells", file=sys.stderr, flush=True)
        check(counts, meta["cell_types_v2"], genes["id"])


        print("Build annotated data...", file=sys.stderr, flush=True)
        # adata = ad.AnnData(counts, meta, genes, dtype=counts.dtype) # Seurat does not support sparse matrices.
        adata = ad.AnnData(counts.todense(), meta, genes, dtype=counts.dtype)
        check(adata.X, adata.obs_names, adata.var_names)


        # print("Save raw annotated data...", file=sys.stderr, flush=True)
        # adata.write(fcounts+".hdf5")
        # print("OK", file=sys.stderr, flush=True)

    elif len(sys.argv) == 2:
        print("Load annotated data...", file=sys.stderr, flush=True)
        adata = ad.read(sys.argv[1])
        check(adata.X, adata.obs_names, adata.var_names)

    else:
        assert(len(sys.argv) == 2 or len(sys.argv) == 4)

    print("Seurat preprocessing...", file=sys.stderr, flush=True)
    scanpy.pp.recipe_seurat(adata)
    check(adata.X, adata.obs_names, adata.var_names)
    adata.write(fcounts+".seurat.hdf5")
    adata.file.close()

    print("Done", file=sys.stderr, flush=True)
