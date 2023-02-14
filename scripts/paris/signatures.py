import scipy
import numpy
from sortedcollections import OrderedSet
from scipy.spatial.distance import squareform, pdist, jaccard
from fastcluster import linkage

# geneset = set of genes.
# e.g. {"A","B"}
# signature = dictionary of all seen genes ("genome"), with value being true if the gene is in the signature.
# e.g. {"A":1, "B":1, "C": 0, "D":0}

def load(filenames, filter_size=None):
    assert(len(filenames) > 0)

    genesets = OrderedSet()
    # current_id = 0
    genome = OrderedSet()
    #genesets_scores = {}
    breaks = []

    for filename in filenames:
        with open(filename) as fd:
            lines = fd.readlines()
            kept_lines = 0
            file_genesets = OrderedSet()
            for line in lines:
                #print(line)
                fields = line.split()
                if len(fields) == 1:
                    fields = line.split(",")
                assert(len(fields)>1)
                i = 0
                if fields[1].isdigit():
                    # Johann’ format: <score> <ngenes> <genes>…
                    score = float(fields[0])
                    n = int(fields[1])
                    i = 2
                if filter_size:
                    if len(fields[i:]) != filter_size:
                        continue
                kept_lines += 1
                genes = frozenset(fields[i:])
                #genesets_scores[genes] = score
                #assert(len(genes) == n)
                genome |= genes # Merge genes into all genes.
                # geneset = {"id":current_id, "score":score, "genes": genes}
                # geneset = {"score":score, "genes": genes}
                # if geneset not in geneset:
                    # geneset.append( geneset )
                # current_id += 1
                file_genesets |= OrderedSet([genes]) # Add a genes as a geneset if not already here.
                #print(genes)
        if kept_lines > 1:
            if filter_size:
                print("File `",filename,"` had",len(file_genesets),"unique signatures among",kept_lines,"signatures of size",filter_size, flush=True)
            else:
                print("File `",filename,"` had",len(file_genesets),"unique signatures among",kept_lines,"signatures of all size", flush=True)
        breaks.append(len(file_genesets))

        genesets |= file_genesets

    return genesets,genome,breaks


def genesets_to_signatures(genesets,genome):
    signatures = numpy.zeros( (len(genesets),len(genome)), bool )
    for s,geneset in enumerate(genesets):
        for g,gene in enumerate(genome):
            # if gene in signature["genome"]:
            if gene in geneset:
                signatures[s][g] = 1
            else:
                signatures[s][g] = 0
    return signatures


def self_similarity_matrix(signatures):
    dissim = pdist(signatures, 'jaccard') # jaccard = DISsimilarity
    return 1-squareform(dissim)


def seriation(Z,N,cur_index):
    '''
        input:
            - Z is a hierarchical tree (dendrogram)
            - N is the number of points given to the clustering process
            - cur_index is the position in the tree for the recursive traversal
        output:
            - order implied by the hierarchical tree Z

        seriation computes the order implied by a hierarchical tree (dendrogram)
    '''
    if cur_index < N:
        return [cur_index]
    else:
        left = int(Z[cur_index-N,0])
        right = int(Z[cur_index-N,1])
        return (seriation(Z,N,left) + seriation(Z,N,right))


def compute_serial_matrix(dist_mat,method="ward"):
    '''
        input:
            - dist_mat is a distance matrix
            - method = ["ward","single","average","complete"]
        output:
            - seriated_dist is the input dist_mat,
              but with re-ordered rows and columns
              according to the seriation, i.e. the
              order implied by the hierarchical tree
            - res_order is the order implied by
              the hierarhical tree
            - res_linkage is the hierarhical tree (dendrogram)

        compute_serial_matrix transforms a distance matrix into
        a sorted distance matrix according to the order implied
        by the hierarchical tree (dendrogram)
    '''
    N = len(dist_mat)
    flat_dist_mat = squareform(dist_mat)
    res_linkage = linkage(flat_dist_mat, method=method,preserve_input=True)
    res_order = seriation(res_linkage, N, N + N-2)
    seriated_dist = numpy.zeros((N,N))
    a,b = numpy.triu_indices(N,k=1)
    seriated_dist[a,b] = dist_mat[ [res_order[i] for i in a], [res_order[j] for j in b]]
    seriated_dist[b,a] = seriated_dist[a,b]

    return seriated_dist, res_order, res_linkage

def colormesh(mat, ax):
    cmap = cm.get_cmap("Blues")
    cmap.set_bad("white")

    pcm = ax.pcolormesh(mat, cmap=cmap)#, edgecolors="#eeeeee", linewidth=0.01)
    #plt.imshow(mat, cmap=cmap)

    # ax.set_aspect("equal")

    xticks=numpy.arange(0,mat.shape[1],max(1,int(mat.shape[1]/100)))
    yticks=numpy.arange(0,mat.shape[0],max(1,int(mat.shape[0]/50)))

    # Major ticks
    ax.set_xticks(xticks+0.5)
    ax.set_yticks(yticks+0.5)
    # Labels
    ax.set_xticklabels(xticks, rotation=90)
    ax.set_yticklabels(yticks, rotation=0)
    # Minor ticks
    # ax.set_xticks(xticks-0.5, minor=True)
    # ax.set_yticks(yticks-0.5, minor=True)
    # ax.tick_params(which="minor", bottom=False, left=False)
    # ax.grid(True, which="minor", axis="both", color="white", linestyle="-", linewidth=1)
    return pcm


if __name__=="__main__":
    import sys
    from matplotlib import gridspec
    import matplotlib.pyplot as plt
    from matplotlib import cm

    size = int(sys.argv[1])
    if size == 0:
        size = None

    sim_thresh = float(sys.argv[2])

    filenames = sys.argv[3:]

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
    gs = gridspec.GridSpec(4, 3, width_ratios=[20,1,1], height_ratios=[10,10,1,1])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[3])
    ax2r = fig.add_subplot(gs[4],sharey=ax2)
    ax2rr = fig.add_subplot(gs[5],sharey=ax2)
    ax2b = fig.add_subplot(gs[6],sharex=ax2)
    ax2bb = fig.add_subplot(gs[9],sharex=ax2)

    # Full distance matrx
    mask = numpy.tri(full_distance_matrix.shape[0], k=-1).transpose() # k>0 = above
    mat = numpy.ma.array(full_distance_matrix, mask=mask)
    pcm = colormesh(mat, ax1)
    fig.colorbar(pcm, ax=ax1)
    ax1.set_title("Jaccard similarities of all {} signatures of size {} on {} genes".format(len(genesets),len(next(iter(genesets))),len(genome)))

    # Dual distance matrix
    sub_distance_matrix = raw_distance_matrix[breaks[1]:,:breaks[1]]
    dual_distance_matrix = sub_distance_matrix[:, (sub_distance_matrix.sum(axis=0)*-1).argsort()]

    colormesh(dual_distance_matrix, ax2)
    ax2.set_ylabel(filenames[0])
    # ax2.xaxis.tick_top()

    # Sum over rows
    sum_rows = dual_distance_matrix.sum(1)
    colormesh(sum_rows.reshape(-1,1), ax2r)
    ax2r.set_ylabel("Sum of similarities for {} signatures in `{}`".format(len(sum_rows),filenames[0]))
    ax2r.yaxis.set_label_position("right")
    ax2r.yaxis.tick_right()

    # Sum over columns
    sum_cols = dual_distance_matrix.sum(0)
    colormesh(sum_cols.reshape(1,-1), ax2b)
    # ax2b.set_xlabel(filenames[1])
    ax2b.set_xlabel("Sum of similarities for {} signatures in `{}`".format(len(sum_cols),filenames[1]))

    def thresh(slice):
        return any([s > sim_thresh for s in slice])

    # Threshold count over rows
    match_rows = numpy.apply_along_axis(thresh, axis=1, arr=dual_distance_matrix)
    colormesh(match_rows.reshape(-1,1), ax2rr)
    sm = sum(match_rows)
    n = len(match_rows)
    ax2rr.set_ylabel("{:.0%} of signatures in `{}` are similar enough (>{}) to at least another one in `{}`".format(sm/n, filenames[0], sim_thresh, filenames[1]))
    ax2rr.yaxis.set_label_position("right")
    ax2rr.yaxis.tick_right()

    # Threshold count over columns
    match_cols = numpy.apply_along_axis(thresh, axis=0, arr=dual_distance_matrix)
    colormesh(match_cols.reshape(1,-1), ax2bb)
    sm = sum(match_cols)
    n = len(match_cols)
    ax2bb.set_xlabel("{:.0%} of signatures in `{}` are similar enough (>{}) to at least another one in `{}`".format(sm/n, filenames[1], sim_thresh, filenames[0]))

    plt.show()
