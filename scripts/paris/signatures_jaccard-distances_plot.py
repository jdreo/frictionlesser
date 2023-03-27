import scipy
import numpy
from sortedcollections import OrderedSet
from scipy.spatial.distance import squareform, pdist, cdist, jaccard
from fastcluster import linkage

# geneset = set of genes.
# e.g. {"A","B"}
# signature = dictionary of all seen genes ("genome"), with value being true if the gene is in the signature.
# e.g. {"A":1, "B":1, "C": 0, "D":0}

from signatures import *

if __name__=="__main__":
    import sys
    from matplotlib import gridspec
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.colors import Normalize

    size = int(sys.argv[1])
    if size == 0:
        size = None

    sim_thresh = float(sys.argv[2])

    filenames = sys.argv[3:5] # Only two files

    genesets,genome,breaks = load(filenames, filter_size=size)
    assert(len(genesets) > 0)
    assert(len(genome) > 0)
    #print(genome)
    #for s in genesets:
    #    print(s)
    print(len(filenames), "files,", len(genesets), "unique genesets,", len(genome), "genes", flush=True)
    sizes = []
    for geneset in genesets:
        s = len(geneset)
        sizes.append(s)
    if not all(s == len(genesets[0]) for s in sizes):
        print("Statistics or sizes:",scipy.stats.describe(sizes), flush=True)

    signatures = genesets_to_signatures(genesets,genome)
    # with numpy.printoptions(threshold=numpy.inf):
        # print(signatures)
    raw_distance_matrix = self_similarity_matrix(signatures)

    # full_distance_matrix, res_order, res_linkage = compute_serial_matrix(raw_distance_matrix, "ward")
    full_distance_matrix = raw_distance_matrix

    # PLOT
    fig = plt.figure(figsize=(10,10))
    fig.tight_layout()
    gs = gridspec.GridSpec(4, 3, width_ratios=[40,1,1], height_ratios=[20,20,1,1])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[3])
    ax2r = fig.add_subplot(gs[4],sharey=ax2)
    ax2rr = fig.add_subplot(gs[5],sharey=ax2)
    ax2b = fig.add_subplot(gs[6],sharex=ax2)
    ax2bb = fig.add_subplot(gs[9],sharex=ax2)

    cmap = cm.get_cmap("Reds")
    normalizer = Normalize(0,1)
    im = cm.ScalarMappable(norm=normalizer, cmap=cmap)
    fig.colorbar(im, ax=[ax1,ax2,ax2r,ax2rr,ax2b,ax2bb])

    # Full distance matrx
    mask = numpy.tri(full_distance_matrix.shape[0], k=-1).transpose() # k>0 = above
    mat = numpy.ma.array(full_distance_matrix, mask=mask)
    pcm = colormesh(mat, ax1, cmap, normalizer)
    ax1.set_title("Jaccard similarities of all {} signatures of size {} on {} genes".format(len(genesets),len(next(iter(genesets))),len(genome)))

    # Dual distance matrix
    genesets_1,genome_1,breaks_1 = load([filenames[0]], filter_size=size)
    assert(len(genesets_1) > 0)
    assert(len(genome_1) > 0)
    signatures_1 = genesets_to_signatures(genesets_1, genome) # against common genome

    genesets_2,genome_2,breaks_2 = load([filenames[1]], filter_size=size)
    assert(len(genesets_2) > 0)
    assert(len(genome_2) > 0)
    signatures_2 = genesets_to_signatures(genesets_2, genome) # against common genome

    # print(signatures_1.shape, signatures_2.shape)

    sub_distance_matrix = 1 - cdist(signatures_1, signatures_2, "jaccard")



    # Sort columns on sum.
    # dual_distance_matrix = sub_distance_matrix[:, (sub_distance_matrix.sum(axis=0)*-1).argsort()] # rows
    # dual_distance_matrix = sub_distance_matrix[(sub_distance_matrix.sum(axis=1)).argsort(), :] # columns
    sorted_distance_matrix = sub_distance_matrix[:, (sub_distance_matrix.sum(axis=0)*-1).argsort()] # rows
    dual_distance_matrix = sorted_distance_matrix[(sorted_distance_matrix.sum(axis=1)).argsort(), :] # columns

    # argsort the columns, sorting by the last row first, then the second-to-last row, continuing up to the first row
    # dual_distance_matrix = sub_distance_matrix[numpy.lexsort(sub_distance_matrix)]
    # dual_distance_matrix = sub_distance_matrix.T[numpy.lexsort(sub_distance_matrix.T)]

    colormesh(dual_distance_matrix, ax2, cmap, normalizer)
    ax2.set_ylabel(filenames[0])
    # ax2.xaxis.tick_top()

    # Sum over rows
    # sum_rows = dual_distance_matrix.mean(1)
    sum_rows = dual_distance_matrix.sum(1)
    colormesh(sum_rows.reshape(-1,1), ax2r, cmap, normalizer)
    ax2r.set_ylabel("Sum of similarities for {} signatures in `{}`".format(len(sum_rows),filenames[0]))
    ax2r.yaxis.set_label_position("right")
    ax2r.yaxis.tick_right()

    # Sum over columns
    # sum_cols = dual_distance_matrix.mean(0)
    sum_cols = dual_distance_matrix.sum(0)
    colormesh(sum_cols.reshape(1,-1), ax2b, cmap, normalizer)
    # ax2b.set_xlabel(filenames[1])
    ax2b.set_xlabel("Sum of similarities for {} signatures in `{}`".format(len(sum_cols),filenames[1]))

    def thresh(slice):
        return any([s > sim_thresh for s in slice])

    # Threshold count over rows
    match_rows = numpy.apply_along_axis(thresh, axis=1, arr=dual_distance_matrix)
    colormesh(match_rows.reshape(-1,1), ax2rr, cmap, normalizer)
    sm = sum(match_rows)
    n = len(match_rows)
    ax2rr.set_ylabel("{:.0%} of signatures in `{}` are similar enough (>{}) to at least another one in `{}`".format(sm/n, filenames[0], sim_thresh, filenames[1]))
    ax2rr.yaxis.set_label_position("right")
    ax2rr.yaxis.tick_right()

    # Threshold count over columns
    match_cols = numpy.apply_along_axis(thresh, axis=0, arr=dual_distance_matrix)
    colormesh(match_cols.reshape(1,-1), ax2bb, cmap, normalizer)
    sm = sum(match_cols)
    n = len(match_cols)
    ax2bb.set_xlabel("{:.0%} of signatures in `{}` are similar enough (>{}) to at least another one in `{}`".format(sm/n, filenames[1], sim_thresh, filenames[0]))

    plt.show()
