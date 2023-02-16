import sys
from scipy.sparse import load_npz
import scanpy as sc
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
        adata = ad.AnnData(counts.tocsr(), meta, genes, dtype=counts.dtype)
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

    # https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html
    print("3k-PBMCs preprocessing...", file=sys.stderr, flush=True)

    # sc.pl.highest_expr_genes(adata, n_top=20, ) # TODO checkpoint

    # print("Basic filtering...", file=sys.stderr, flush=True)
    # sc.pp.filter_cells(adata, min_genes=200)
    # sc.pp.filter_genes(adata, min_cells=3)

    # adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    # sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True) # TODO checkpoint
    # sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt') # TODO checkpoint
    # sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts') # TODO checkpoint

    # print("Actually do the basic filtering...", file=sys.stderr, flush=True)
    # adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    # adata = adata[adata.obs.pct_counts_mt < 5, :]

    print("Total-count normalize...", file=sys.stderr, flush=True)
    # (library-size correct) the data matrix X to 10,000 reads per cell, so that counts become comparable among cells.
    sc.pp.normalize_total(adata, target_sum=1e4)

    print("Logarithmize the data...", file=sys.stderr, flush=True)
    sc.pp.log1p(adata)

    print("Identify highly-variable genes...", file=sys.stderr, flush=True)
    # sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pp.highly_variable_genes(adata, n_top_genes=10000)

    # sc.pl.highly_variable_genes(adata) # TODO checkpoint

    print("Actually do the highly-variable filtering...", file=sys.stderr, flush=True)
    adata = adata[:, adata.var.highly_variable]

    check(adata.X, adata.obs_names, adata.var_names)

    print("Write data to", fcounts+".pbmc3k.hdf5", file=sys.stderr, flush=True)
    adata.write(fcounts+".pbmc3k.hdf5")

    adata.file.close() # In case it was backed.
    print("Done", file=sys.stderr, flush=True)

