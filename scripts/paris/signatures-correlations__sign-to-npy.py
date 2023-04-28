#!/usr/bin/env python3

import scipy
import numpy
# import seaborn
# import matplotlib.pyplot as plt

from signatures import *

if __name__=="__main__":
    import sys
    import csv

    assert(len(sys.argv) == 9)

    franks = sys.argv[1]
    fgenes = sys.argv[2]
    findices = sys.argv[3]

    size = int(sys.argv[4])
    if size == 0:
        size = None

    fout = sys.argv[5]
    forigins = sys.argv[6]

    f_A = sys.argv[7]
    f_B = sys.argv[8]

    genes_scores = numpy.load(fgenes)

    print("Recover number of samples...", file=sys.stderr, flush=True)
    with open(franks) as fd:
        data = csv.reader(fd, delimiter="\t")
        header = next(data)
        # Indices of each sample's cells. 
        samples = {}
        for i,s in enumerate(header[1:]):
            samples.setdefault(s, []).append(i)
        assert(len(samples) > 0)
    nsamples = len(samples)

    print("Load indices...", file=sys.stderr, flush=True)
    genes_indices = {}
    with open(findices) as fd: # FIXME this should really be an AnnData
        data = csv.reader(fd, delimiter=",")
        header = next(data)
        for row in data:
            genes_indices[row[1]] = int(row[0])

    # FIXME handle more than two origin files.
    print("Load signatures...", file=sys.stderr, flush=True)
    sign_A, genome_A, breaks_A = load([f_A], size)
    sign_B, genome_B, breaks_B = load([f_B], size)
    signatures = sign_A | sign_B

    print("Compute sum of scores...", file=sys.stderr, flush=True)
    ncells = genes_scores.shape[1]
    sign_scoresums = numpy.zeros((len(signatures),ncells))
    for i,signature in enumerate(signatures):
        sign_gene_indices = [genes_indices[g] for g in genes_indices if g in signature]
        # Extract genes (rows) and sum over rows.
        sign_scoresums[i] = numpy.sum(genes_scores[sign_gene_indices], axis=0) * 1/nsamples
        print("\r", i, "/", len(signatures), end=" ", file=sys.stderr, flush=True)

    # print("Plot scores map...", file=sys.stderr, flush=True)
    # ax = seaborn.heatmap(sign_scoresums)
    # plt.title("Sum of scores across signatures, for genes on each cell.")
    # ax.figure.savefig("signatures-sum-standardized-score.png", dpi=600)

    print("Compute dot products...", file=sys.stderr, flush=True)
    correlations = numpy.zeros((len(signatures),len(signatures)))
    signature_origins = {}
    for i,s1 in enumerate(signatures):
        # Get origins while we're here.
        if s1 in sign_A and s1 in sign_B:
            signature_origins[i] = "BOTH"
        elif s1 in sign_A:
            signature_origins[i] = f_A
        elif s1 in sign_B:
            signature_origins[i] = f_B
        else:
            assert(s1 in sign_A or s1 in sign_B)

        # FIXME compute only the triangular matrix?
        for j,s2 in enumerate(signatures):
            correlations[i,j] = numpy.dot(sign_scoresums[i], sign_scoresums[j])

    print("Save correlations...", correlations.shape, "to `", fout, "`...", file=sys.stderr, flush=True)
    numpy.save(fout, correlations)

    assert(len(signature_origins) == len(signatures))
    print("Save unique signatures origin table to `", forigins, "`...", file=sys.stderr, flush=True)

    with open(forigins, 'w') as fd:
        writer = csv.DictWriter(fd, fieldnames=["signature_index","origin"])
        writer.writeheader()
        for i in signature_origins:
            writer.writerow({"signature_index":i, "origin":signature_origins[i]})

    print("Done", file=sys.stderr, flush=True)
