import sys
import csv
import scipy
import seaborn
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

if __name__ == "__main__":

    assert(len(sys.argv) == 4)
    fgccorr = sys.argv[1]
    fplot = sys.argv[2]
    fout = sys.argv[3]

    print("Load annotated genes correlations data from: ", fgccorr, file=sys.stderr, flush=True)
    cells_sgenes = ad.read(fgccorr)
    print(cells_sgenes, file=sys.stderr, flush=True)

    # Data transformation.
    # data = -1 * np.log10(np.abs(correlations)+1)
    # data = np.abs(correlations)
    data = -1 * np.log10(np.abs(cells_sgenes.varp["p-values"])+1e-20)

    print("Compute the clustering...", file=sys.stderr, flush=True)
    linkage = scipy.cluster.hierarchy.complete(data)
    membership="distance"
    memberships = scipy.cluster.hierarchy.fcluster(linkage, criterion='maxclust', t = 10 )
    # memberships = scipy.cluster.hierarchy.fcluster(linkage, criterion=membership, t = tthresh)

    print("Save cluster membership...", file=sys.stderr, flush=True)
    # np.save(fout, membership)
    with open(fout, 'w') as fd:
        writer = csv.DictWriter(fd, fieldnames=["signature_id","cluster_id"])
        writer.writeheader()
        for s,c in enumerate(memberships):
            writer.writerow({"signature_id":s, "cluster_id":c})

    print("Plot...", file=sys.stderr, flush=True)
    import sys
    sys.setrecursionlimit(999999)

    palette_memberships = seaborn.color_palette("Paired",len(memberships))
    memberships_colors = [palette_memberships[i] for i in memberships]

    # First, plot the regular clustermap.
    linkage="complete"
    cmap = "viridis"
    # fig = seaborn.clustermap(data, method=linkage, row_colors=[origins, memberships_colors], col_colors=[origins,memberships_colors], cmap=cmap)#, vmin=0, vmax=5)#, vmin=0, vmax=1)#vmax=np.log10(1+1))
    fig = seaborn.clustermap(data, method=linkage, row_colors=[memberships_colors], col_colors=[memberships_colors], cmap=cmap)#, vmin=0, vmax=5)#, vmin=0, vmax=1)#vmax=np.log10(1+1))

    # # Then, apply mask based on p-values.
    # pmask = pvalues > pthresh
    # trimask = np.tril(np.ones(data.shape))
    # mask = pmask
    # values = fig.ax_heatmap.collections[0].get_array().reshape(data.shape)
    # masked_values = np.ma.array(values, mask = mask)
    # fig.ax_heatmap.collections[0].set_array(masked_values)
    # # Color of masked values.
    # fig.ax_heatmap.set_facecolor("white")

    # Origins legend.
    # handles = [Patch(facecolor=color) for color in palette_origins]
    # plt.legend(handles, names, title='Origin',
    #        bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right')
    # plt.title("Hierarchical clustering ({linkage} linkage) of average correlations across samples, between all signatures. Only correlations having p-values ≤ {pthresh} are shown. Membership on {membership} with parameter t={tthresh}".format(linkage=linkage, pthresh=pthresh, membership=membership, tthresh=tthresh))
    plt.title("Hierarchical clustering ({linkage} linkage) of average correlations across samples, between all signatures. Membership on {membership}".format(linkage=linkage, membership=membership))
    # FIXME : see 10% quantiles instead of mean (move 1/#samples to the signatures computation? and replace with other statistics).

    fig.savefig(fplot, dpi=600)
