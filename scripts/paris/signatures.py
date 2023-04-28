import sys
import scipy
import numpy
from sortedcollections import OrderedSet
from scipy.spatial.distance import squareform, pdist, cdist, jaccard
from fastcluster import linkage

# geneset = set of genes.
# e.g. {"A","B"}
# signature = dictionary of all seen genes ("genome"), with value being true if the gene is in the signature.
# e.g. {"A":1, "B":1, "C": 0, "D":0}

def load_genome(filenames, filter_size=None):
    """Return the set of genes encountered in all the given signature files."""
    assert(len(filenames) > 0)

    genome = OrderedSet()

    for filename in filenames:
        with open(filename) as fd:
            lines = fd.readlines()
            for line in lines:
                #print(line)
                fields = line.strip().split()
                if len(fields) == 1:
                    fields = line.strip().split(",")
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
                genes = frozenset(fields[i:])
                genome |= genes # Merge genes into all genes.

    return genome


def load(filenames, filter_size=None):
    if filter_size == 0:
        filter_size = None
    assert(len(filenames) > 0)

    genesets = OrderedSet()
    # current_id = 0
    genome = OrderedSet()
    #genesets_scores = {}
    breaks = []

    for filename in filenames:
        kept_lines = 0
        with open(filename) as fd:
            lines = fd.readlines()
            file_genesets = OrderedSet()
            for line in lines:
                fields = line.strip().split()
                # print(len(fields),fields)
                if len(fields) == 1:
                    fields = line.strip().split(",")
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
                file_genesets |= OrderedSet([genes]) # Add a signature as a geneset if not already here.
                #print(genes)
        # if kept_lines > 1:
        if filter_size:
            print("File `",filename,"` had",len(file_genesets),"unique signatures among",kept_lines,"signatures of size",filter_size, file=sys.stderr, flush=True)
        else:
            print("File `",filename,"` had",len(file_genesets),"unique signatures among",kept_lines,"signatures of all size", file=sys.stderr, flush=True)

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

def colormesh(mat, ax, cmap, norm):
    # cmap = cm.get_cmap("Blues")
    cmap.set_bad("white")

    pcm = ax.pcolormesh(mat, cmap=cmap, norm=norm)#, edgecolors="#eeeeee", linewidth=0.01)
    #plt.imshow(mat, cmap=cmap)

    # ax.set_aspect("equal")

    xticks=numpy.arange(0,mat.shape[1],max(1,int(mat.shape[1]/40)))
    yticks=numpy.arange(0,mat.shape[0],max(1,int(mat.shape[0]/30)))

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
    from matplotlib.colors import Normalize

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

