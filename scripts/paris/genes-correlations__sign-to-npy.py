import scipy
import numpy
from sortedcollections import OrderedSet
from scipy.spatial.distance import squareform, pdist, cdist, jaccard
from fastcluster import linkage
from scipy.stats import spearmanr

from signatures import *

if __name__=="__main__":
    import sys
    from matplotlib import gridspec
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.colors import Normalize
    # import pandas
    import csv

    rankfile = sys.argv[1]

    size = int(sys.argv[2])
    if size == 0:
        size = None

    signatures_files = sys.argv[3:5] # Only two files


    genome = load_genome(signatures_files, filter_size=size)
    assert(len(genome) > 0)
    print("Genome across",len(signatures_files), "files:", len(genome), "genes", file=sys.stderr, flush=True)

    print("Load ranks...", file=sys.stderr, flush=True)
    genes_idx = {}
    ranks_l = []
    with open(rankfile) as fd:
        csv = csv.reader(fd, delimiter=" ")
        header = next(csv)
        n=0
        for row in csv:
            gene = row[0]
            assert(len(gene) > 0)
            if gene in genome:
                genes_idx[gene] = n
                n += 1
                print("\r",n,end=" ", file=sys.stderr,flush=True)
                ranks_l.append(numpy.array([float(i) for i in row[1:]]))
    assert(len(ranks_l) > 0)
    assert(len(ranks_l) == len(genes_idx))
    ngenes = len(genes_idx)
    ncells = len(ranks_l[0])
    print("\rLoaded ranks for",ngenes,"genes,",ncells,"cells", file=sys.stderr, flush=True)
    # Convert to array.
    ranks = numpy.array(ranks_l)
    assert(ranks.shape == (ngenes,ncells))

    ranked_genes = OrderedSet(genes_idx.keys())
    unloaded_genes = genome - ranked_genes
    if len(unloaded_genes) > 0:
        print("WARNING: Genes of genome not found in ranks file:",unloaded_genes, file=sys.stderr, flush=True)
        # sys.exit(1)
    # common_genome = genome & ranked_genes

    print("Compute genes correlation matrix...", file=sys.stderr, flush=True)
    genes_correlations = numpy.absolute(spearmanr(ranks, axis=1)[0]) # [0] = statistic, [1] = pvalue
    assert(genes_correlations.shape == (ngenes,ngenes))

    fname = "genes-ranks-correlations.npy"
    print("Save genes correlation matrix to `", fname,"`...", file=sys.stderr, flush=True)
    numpy.save(fname, genes_correlations)

    print("Done", file=sys.stderr, flush=True)
