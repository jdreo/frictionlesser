import sys
import numpy
# import scanpy
# import pandas as pd
import anndata as ad

def check(counts, meta, genes):
    assert(counts.shape[0] > 0)
    assert(counts.shape[1] > 0)
    assert(counts.shape[0] == len(meta))
    assert(counts.shape[1] == len(genes))
    print("OK", file=sys.stderr, flush=True)


if __name__ == "__main__":

    print("Load annotated data...", file=sys.stderr, flush=True)
    adata = ad.read(sys.argv[1])
    check(adata.X, adata.obs_names, adata.var_names)
    print("Loaded",len(adata.var_names),"genes and",len(adata.obs_names),"cells", file=sys.stderr, flush=True)

    print("Output as CSV...", file=sys.stderr, flush=True)
    counts_full = numpy.asarray(adata.X.todense())

    # Column names are sample, not cells IDs.
    print("GENE,", ",".join(adata.obs["sample"]), sep="")
    i = 0
    for row in counts_full.T:
        print("\r",i, end="", flush=True, file=sys.stderr)
        print(adata.var["id"][i], ",".join(str(x) for x in row))
        i += 1
    assert( i == len(adata.var_names))

    print("\nDone", file=sys.stderr, flush=True)
