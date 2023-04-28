import numpy
import scipy
import seaborn
import sys
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import csv

if __name__ == "__main__":

    assert(len(sys.argv) == 8)

    fcorrs = sys.argv[1]
    forigins = sys.argv[2]
    fpvalues = sys.argv[3]
    pthresh = float(sys.argv[4])
    tthresh = float(sys.argv[5])

    fplot = sys.argv[6]
    fout = sys.argv[7]

    print("Load correlations...", file=sys.stderr, flush=True)
    correlations = numpy.load(fcorrs)

    print("Load p-values...", file=sys.stderr, flush=True)
    pvalues = numpy.load(fpvalues)

    print("Load the csv file of origins...", file=sys.stderr, flush=True)
    nsigns = correlations.shape[0]
    origins = [""] * nsigns
    names = set()
    with open(forigins) as fd:
        csvreader = csv.DictReader(fd)
        for row in csvreader:
            names.add(row["origin"])
            origins[int(row["signature_index"])] = row["origin"]

    print("Make a qualitative palette from the origins...", file=sys.stderr, flush=True)
    palette_origins = seaborn.color_palette("Set2",len(names))
    for i,orig in enumerate(origins):
        for j,n in enumerate(names):
            if orig == n:
                origins[i] = palette_origins[j]
                break

    # Data transformation.
    # data = numpy.log10(numpy.abs(correlations)+1)
    data = numpy.abs(correlations)

    print("Compute the clustering...", file=sys.stderr, flush=True)
    linkage = scipy.cluster.hierarchy.complete(data)
    membership="distance"
    # memberships = scipy.cluster.hierarchy.fcluster(linkage, criterion='maxclust', t = 10 )
    memberships = scipy.cluster.hierarchy.fcluster(linkage, criterion=membership, t = tthresh)

    print("Save cluster membership...", file=sys.stderr, flush=True)
    # numpy.save(fout, membership)
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
    fig = seaborn.clustermap(data, method=linkage, row_colors=[origins, memberships_colors], col_colors=[origins,memberships_colors], cmap=cmap, vmin=0, vmax=1)#vmax=numpy.log10(1+1))

    # Then, apply mask based on p-values.
    pmask = pvalues > pthresh
    trimask = numpy.tril(numpy.ones(data.shape))
    mask = pmask
    values = fig.ax_heatmap.collections[0].get_array().reshape(data.shape)
    masked_values = numpy.ma.array(values, mask = mask)
    fig.ax_heatmap.collections[0].set_array(masked_values)
    # Color of masked values.
    fig.ax_heatmap.set_facecolor("white")

    # Origins legend.
    handles = [Patch(facecolor=color) for color in palette_origins]
    plt.legend(handles, names, title='Origin',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right')
    plt.title("Hierarchical clustering ({linkage} linkage) of average correlations across samples, between all signatures. Only correlations having p-values ≤ {pthresh} are shown. Membership on {membership} with parameter t={tthresh}".format(linkage=linkage, pthresh=pthresh, membership=membership, tthresh=tthresh))
    # FIXME : see 10% quantiles instead of mean (move 1/#samples to the signatures computation? and replace with other statistics).

    fig.savefig(fplot, dpi=600)
