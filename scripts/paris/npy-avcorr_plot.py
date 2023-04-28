import numpy

from signatures import *

if __name__=="__main__":
    import sys
    from matplotlib import gridspec
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.colors import Normalize

    fgenescorr = sys.argv[1]
    fsignavcorr = sys.argv[2]
    sim_thresh = float(sys.argv[3])

    if len(sys.argv) == 6:
        legend_A = sys.argv[4]
        legend_B = sys.argv[5]
    else:
        legend_A = "A"
        legend_B = "B"

    genes_corr_matrix = numpy.load(fgenescorr)
    print("Loaded a genes correlations matrix of shape", genes_corr_matrix.shape, file=sys.stderr, flush=True)

    distance_matrix = numpy.load(fsignavcorr)
    print("Loaded a signatures correlations matrix of shape", distance_matrix.shape, file=sys.stderr, flush=True)

    # PLOT
    fig = plt.figure(figsize=(10,10))
    # fig.tight_layout()
    gs = gridspec.GridSpec(4, 3, width_ratios=[40,1,1], height_ratios=[20,20,1,1])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[3])
    ax2r = fig.add_subplot(gs[4],sharey=ax2)
    ax2rr = fig.add_subplot(gs[5],sharey=ax2)
    ax2b = fig.add_subplot(gs[6],sharex=ax2)
    ax2bb = fig.add_subplot(gs[9],sharex=ax2)

    # # PLOT
    # fig = plt.figure(figsize=(10,10))
    # # fig.tight_layout()
    # gs = gridspec.GridSpec(3, 3, width_ratios=[40,1,1], height_ratios=[20,1,1])
    # ax2 = fig.add_subplot(gs[0])
    # ax2r = fig.add_subplot(gs[1],sharey=ax2)
    # ax2rr = fig.add_subplot(gs[2],sharey=ax2)
    # ax2b = fig.add_subplot(gs[3],sharex=ax2)
    # ax2bb = fig.add_subplot(gs[6],sharex=ax2)

    cmap = cm.get_cmap("Reds")
    normalizer = Normalize(0,1)
    im = cm.ScalarMappable(norm=normalizer, cmap=cmap)
    # fig.colorbar(im, ax=[ax2,ax2r,ax2rr,ax2b,ax2bb])
    fig.colorbar(im, ax=[ax1,ax2,ax2r,ax2rr,ax2b,ax2bb])


    # Genes ranks correlations.
    # Sort columns on sum.
    sorted_genes_corr_matrix = genes_corr_matrix[:, (genes_corr_matrix.sum(axis=0)*-1).argsort()] # rows
    dual_genes_corr_matrix = sorted_genes_corr_matrix[(sorted_genes_corr_matrix.sum(axis=1)).argsort(), :] # columns

    mask = numpy.tri(dual_genes_corr_matrix.shape[0], k=-1).transpose() # k>0 = above
    mat = numpy.ma.array(numpy.rot90(dual_genes_corr_matrix), mask=mask)
    # mat = numpy.ma.array(numpy.rot90(dual_genes_corr_matrix))
    pcm = colormesh(mat, ax1, cmap, normalizer)
    ax1.set_title("Sorted ranks correlations between all {} genes".format(genes_corr_matrix.shape[0]))


    # Sort columns on sum.
    sorted_distance_matrix = distance_matrix[:, (distance_matrix.sum(axis=0)*-1).argsort()] # rows
    dual_distance_matrix = sorted_distance_matrix[(sorted_distance_matrix.sum(axis=1)).argsort(), :] # columns

    colormesh(dual_distance_matrix, ax2, cmap, normalizer)
    ax2.set_ylabel(legend_A)
    # ax2.xaxis.tick_top()
    ax2.set_title("Sorted average ranks correlations between {} signatures".format(distance_matrix.shape))

    # Sum over rows
    sum_rows = dual_distance_matrix.mean(1)
    # sum_rows = dual_distance_matrix.sum(1)
    colormesh(sum_rows.reshape(-1,1), ax2r, cmap, normalizer)
    ax2r.set_ylabel("Average correlations for {} signatures in `{}`".format(len(sum_rows),legend_A))
    ax2r.yaxis.set_label_position("right")
    ax2r.yaxis.tick_right()

    # Sum over columns
    sum_cols = dual_distance_matrix.mean(0)
    # sum_cols = dual_distance_matrix.sum(0)
    colormesh(sum_cols.reshape(1,-1), ax2b, cmap, normalizer)
    # ax2b.set_xlabel("1")
    ax2b.set_xlabel("Average correlations for {} signatures in `{}`".format(len(sum_cols),legend_B))

    def thresh(slice):
        return any([s > sim_thresh for s in slice])

    # Threshold count over rows
    match_rows = numpy.apply_along_axis(thresh, axis=1, arr=dual_distance_matrix)
    colormesh(match_rows.reshape(-1,1), ax2rr, cmap, normalizer)
    sm = sum(match_rows)
    n = len(match_rows)
    ax2rr.set_ylabel("{:.0%} of signatures in `{}` are similar enough (>{}) to at least another one in `{}`".format(sm/n, legend_A, sim_thresh, legend_B))
    ax2rr.yaxis.set_label_position("right")
    ax2rr.yaxis.tick_right()

    # Threshold count over columns
    match_cols = numpy.apply_along_axis(thresh, axis=0, arr=dual_distance_matrix)
    colormesh(match_cols.reshape(1,-1), ax2bb, cmap, normalizer)
    sm = sum(match_cols)
    n = len(match_cols)
    ax2bb.set_xlabel("{:.0%} of signatures in `{}` are similar enough (>{}) to at least another one in `{}`".format(sm/n, legend_B, sim_thresh, legend_A))

    plt.show()

