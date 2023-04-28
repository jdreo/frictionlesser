import sys
import seaborn
import numpy as np
import anndata as ad

import signatures

if __name__ == "__main__":

    assert(len(sys.argv) >= 5)
    fsccorr = sys.argv[1]
    size = int(sys.argv[2])
    if size == 0:
        size = None
    fplot_selfcorr = sys.argv[3]
    fsignatures = sys.argv[4:]

    print("Load annotated cells-signatures data from: ", fsccorr, file=sys.stderr, flush=True)
    cells_signs = ad.read(fsccorr)
    print(cells_signs, file=sys.stderr, flush=True)
    ncells = cells_signs.shape[0]
    # ngenes = cells_signs.shape[1]

    print("Load signatures from: ", fsignatures, file=sys.stderr, flush=True)
    unique_signatures, unique_scores, origins = signatures.load_unique_signatures(fsignatures, size)
    nsignatures = len(unique_signatures)
    print("Found",nsignatures,"unique signatures", file=sys.stderr, flush=True)

    ###########################################################################
    # PLOT
    ###########################################################################

    print("Plot self-correlations of signatures against scores...", file=sys.stderr, flush=True)
    selfcorr = np.abs(np.diag(cells_signs.varp["correlations"]))
    if __debug__:
        for c in selfcorr:
            assert(-1 <= c <= 1)

    objf = [unique_scores[frozenset(s.split())] for s in cells_signs.var["signature"]]
    if __debug__:
        for s in objf:
            assert(s is not None)
            assert(s > 0)

    plot = seaborn.scatterplot(x=objf, y=selfcorr, s=5)
    plot.set(xlabel="objective function", ylabel="self-correlation",
             title="Relationship between\nobjective function and self-correlation of genes, for signatures in\n{}".format(fsignatures))
    fig = plot.get_figure()
    fig.savefig(fplot_selfcorr, dpi=600)

    # plot = seaborn.scatterplot(x=objf, y=np.abs(np.diag(cells_signs.varp["av.correlations"])), s=5)
    # plot.set(xlabel="objective function", ylabel="self-correlation",
    #          title="Relationship between\nobjective function and self-correlation of genes, for signatures in\n{}".format(fsignatures))
    # fig = plot.get_figure()
    # fig.savefig(fplot_selfcorr+"_2.png", dpi=600)

    # FIXME compute linear correlation coefficient and assert

    print("Done", file=sys.stderr, flush=True)

    # FIXME move qc plot in a separate task?




# #!/usr/bin/env python3

# import scipy
# import numpy
# import seaborn
# # import matplotlib.pyplot as plt

# from signatures import *

# if __name__=="__main__":
#     import sys
#     import csv

#     assert(len(sys.argv) >= 8)

#     franks = sys.argv[1]
#     fgenes = sys.argv[2]
#     findices = sys.argv[3]

#     filter_size = int(sys.argv[4])
#     if filter_size == 0:
#         filter_size = None

#     fout = sys.argv[5]
#     fplot = sys.argv[6]

#     fsignatures = sys.argv[7:]

#     genes_scores = numpy.load(fgenes)

#     print("Recover number of samples...", file=sys.stderr, flush=True)
#     with open(franks) as fd:
#         data = csv.reader(fd, delimiter="\t")
#         header = next(data)
#         # Indices of each sample's cells. 
#         samples = {}
#         for i,s in enumerate(header[1:]):
#             samples.setdefault(s, []).append(i)
#         assert(len(samples) > 0)
#     nsamples = len(samples)

#     print("Load indices...", file=sys.stderr, flush=True)
#     genes_indices = {}
#     with open(findices) as fd: # FIXME this should really be an AnnData
#         data = csv.reader(fd, delimiter=",")
#         header = next(data)
#         for row in data:
#             genes_indices[row[1]] = int(row[0])

#     # FIXME handle more than two origin files.
#     print("Load signatures...", file=sys.stderr, flush=True)

#     # signatures, genome, breaks = load([f_A], size)
#     signatures = OrderedSet()
#     sign_scores = {}
#     for filename in fsignatures:
#         kept_lines = 0
#         with open(filename) as fd:
#             file_genesets = OrderedSet()
#             lines = fd.readlines()
#             for line in lines:
#                 fields = line.strip().split()
#                 if len(fields) == 1:
#                     fields = line.strip().split(",")
#                 assert(len(fields)>1)
#                 i = 0
#                 if fields[1].isdigit():
#                     # Johann’ format: <score> <ngenes> <genes>…
#                     score = float(fields[0])
#                     n = int(fields[1])
#                     i = 2
#                 else:
#                     score = None
#                 if filter_size:
#                     if len(fields[i:]) != filter_size:
#                         continue
#                 kept_lines += 1
#                 genes = frozenset(fields[i:])
#                 file_genesets |= OrderedSet([genes]) # Add a signature as a geneset if not already here.
#                 sign_scores[genes] = score
#         if filter_size:
#             print("File `",filename,"` had",len(file_genesets),"unique signatures among",kept_lines,"signatures of size",filter_size, file=sys.stderr, flush=True)
#         else:
#             print("File `",filename,"` had",len(file_genesets),"unique signatures among",kept_lines,"signatures of all size", file=sys.stderr, flush=True)

#         signatures |= file_genesets

#     print("Compute sum of scores...", file=sys.stderr, flush=True)
#     ncells = genes_scores.shape[1]
#     sign_scoresums = numpy.zeros((len(signatures),ncells))
#     for i,signature in enumerate(signatures):
#         sign_gene_indices = [genes_indices[g] for g in genes_indices if g in signature]
#         # Extract genes (rows) and sum over rows.
#         sign_scoresums[i] = numpy.sum(genes_scores[sign_gene_indices], axis=0) * 1/nsamples
#         print("\r", i, "/", len(signatures), end=" ", file=sys.stderr, flush=True)

#     # print("Plot scores map...", file=sys.stderr, flush=True)
#     # ax = seaborn.heatmap(sign_scoresums)
#     # plt.title("Sum of scores across signatures, for genes on each cell.")
#     # ax.figure.savefig("signatures-sum-standardized-score.png", dpi=600)

#     print("Compute dot products...", file=sys.stderr, flush=True)
#     correlations = numpy.zeros(len(signatures))
#     for i,s1 in enumerate(signatures):
#         correlations[i] = numpy.dot(sign_scoresums[i], sign_scoresums[i])

#     print("Save self-correlations...", correlations.shape, "to `", fout, "`...", file=sys.stderr, flush=True)
#     numpy.save(fout, correlations)


#     selfcorr = []
#     objf = []
#     for i,s in enumerate(signatures):
#         selfcorr.append(correlations[i])
#         assert(-1 <= correlations[i] <= 1)
#         objf.append(sign_scores[s])
#         assert(sign_scores[s] is not None)
#         assert(sign_scores[s] > 0)

#     plot = seaborn.scatterplot(x=objf, y=selfcorr, s=5)
#     plot.set(xlabel="objective function", ylabel="self-correlation",
#              title="Relationship between\nobjective function and self-correlation of genes, for signatures in\n{}".format(fsignatures))
#     fig = plot.get_figure()
#     fig.savefig(fplot, dpi=600)

#     print("Done", file=sys.stderr, flush=True)
