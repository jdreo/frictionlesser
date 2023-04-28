import numpy
import scipy
import seaborn
import sys
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import csv

if __name__ == "__main__":

    fsignavcorrfull = sys.argv[1]
    forigins = sys.argv[2]

    sign_avcorr_full = numpy.load(fsignavcorrfull)

    # Load the csv file of origins.
    nsigns = sign_avcorr_full.shape[0]
    origins = [""] * nsigns
    names = set()
    with open(forigins) as fd:
        csvreader = csv.DictReader(fd)
        for row in csvreader:
            names.add(row["origin"])
            origins[int(row["gene_index"])] = row["origin"]

    # Make a qualitative palette from the origins.
    palette_origins = seaborn.color_palette("Set2",len(names))
    for i,orig in enumerate(origins):
        for j,n in enumerate(names):
            if orig == n:
                origins[i] = palette_origins[j]
                break

    # Compute the clustering.
    linkage = scipy.cluster.hierarchy.complete(sign_avcorr_full)
    memberships = scipy.cluster.hierarchy.fcluster(linkage, t = 10, criterion='maxclust')

    palette_memberships = seaborn.color_palette("Paired",len(memberships))
    memberships_colors = [palette_memberships[i] for i in memberships]

    # Plot.
    fig = seaborn.clustermap(sign_avcorr_full, method="complete", row_colors=[origins, memberships_colors], col_colors=[origins,memberships_colors], cmap="crest")

    # Origins legend.
    handles = [Patch(facecolor=color) for color in palette_origins]
    plt.legend(handles, names, title='Origin',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right')
    plt.title("Hierarchical clustering (complete linkage) of average correlations between all signatures.")

    fig.savefig("clustermap.png", dpi=600)
