import anndata as ad
import signatures

if __name__ == "__main__":
    import sys

    assert(len(sys.argv) == 4)

    ranks = signatures.load_ranks_csv(sys.argv[1])
    adata = ad.read(sys.argv[2])
    adata.layers["ranks"] = ranks.T
    adata.write(sys.argv[3], compression="gzip")

    print(adata, file=sys.stderr, flush=True)
