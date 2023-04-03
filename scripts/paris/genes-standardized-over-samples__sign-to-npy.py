#!/usr/bin/env python3

import scipy
import numpy
from sortedcollections import OrderedSet
# import seaborn
# import matplotlib.pyplot as plt

from signatures import *

if __name__=="__main__":
    import sys
    import csv

    assert(len(sys.argv) == 7)

    rankfile = sys.argv[1]

    size = int(sys.argv[2])
    if size == 0:
        size = None

    signatures_files = sys.argv[3:5] # Only two files

    fout = sys.argv[5]
    # fgenes = fout+".genes.csv"
    fgenes = sys.argv[6]

    genome = load_genome(signatures_files, filter_size=size)
    assert(len(genome) > 0)
    print("Genome across",len(signatures_files), "files:", len(genome), "genes", file=sys.stderr, flush=True)

    print("Load ranks...", file=sys.stderr, flush=True)
    genes = []
    scores_l = []
    with open(rankfile) as fd:
        data = csv.reader(fd, delimiter="\t")
        header = next(data)
        # Indices of each sample's cells. 
        samples = {}
        for i,s in enumerate(header[1:]):
            samples.setdefault(s, []).append(i)
        assert(len(samples) > 0)
        n=0
        for row in data:
            i = 0
            gene = row[0]
            assert(len(gene) > 0)
            if gene in genome:
                # [Index of gene, gene name].
                genes.append([n,gene])
                n += 1
                print("\r", n, end=" ", file=sys.stderr, flush=True)
                ranks_row = numpy.array([float(r) for r in row[1:]])
                scores_row = numpy.zeros(ranks_row.shape)
                for s in samples:
                    # All ranks in this sample at this gene.
                    svec = ranks_row[samples[s]]
                    assert(len(samples) > 0)
                    assert(len(svec) > 0)
                    stdev = numpy.std(svec)
                    if stdev != 0:
                        # Clipped, because sometime floating-point precision comes in the way.
                        # nvec = numpy.clip((svec - numpy.mean(svec)) / stdev * 1/len(samples) * 1/numpy.sqrt(len(svec)-1), -1, 1)
                        nvec = numpy.clip((svec - numpy.mean(svec)) / stdev * 1/numpy.sqrt(len(svec)-1), -1, 1)
                    else:
                        nvec = numpy.zeros(svec.shape)
                    # Replace ranks at each sample's cells with the standardized scores.
                    scores_row[samples[s]] = nvec
                scores_l.append(scores_row)
    assert(len(scores_l) > 0)
    assert(len(scores_l) == len(genes))
    ngenes = len(genes)
    ncells = len(scores_l[0])
    print("\rComputed standardized scores for",ngenes,"genes,",ncells,"cells in", len(samples), "samples.", file=sys.stderr, flush=True)
    # Convert to array.
    print("Convert to array...", file=sys.stderr, flush=True)
    ranks = numpy.array(scores_l)
    print(ranks, file=sys.stderr, flush=True)
    assert(ranks.shape == (ngenes,ncells))

    print("Check signatures consistency...", file=sys.stderr, flush=True)
    ranked_genes = OrderedSet([g[1] for g in genes])
    unloaded_genes = genome - ranked_genes
    if len(unloaded_genes) > 0:
        print("WARNING: Genes of genome not found in ranks file:",unloaded_genes, file=sys.stderr, flush=True)
    else:
        print("OK", file=sys.stderr, flush=True)

    print("Save standardized scores matrix of shape", ranks.shape, "to `", fout, "`...", file=sys.stderr, flush=True)
    numpy.save(fout, ranks)

    print("Save genes indices to `", fgenes, "`...", file=sys.stderr, flush=True)
    with open(fgenes, 'w') as fd:
        writer = csv.DictWriter(fd, fieldnames=["gene_index","gene_name"])
        writer.writeheader()
        for g in genes:
            writer.writerow({"gene_index":g[0], "gene_name":g[1]})

    # print("Plot ranks map...", file=sys.stderr, flush=True)
    # fig = seaborn.heatmap(ranks)
    # plt.title("Standardized score of genes (rows) for each cells) column.")
    # fig.savefig("genes-standardized-over-samples.png", dpi=600)

    print("Done", file=sys.stderr, flush=True)
