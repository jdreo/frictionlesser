import sys
import csv
import operator
import seaborn
import pandas

import signatures

if __name__ == "__main__":

    assert(len(sys.argv) >= 5)

    size = int(sys.argv[1])
    if size == 0:
        size = None
    franks = sys.argv[2]
    fplot = sys.argv[3]
    fsignatures = sys.argv[4:]

    print("Load all gene names...", file=sys.stderr, flush=True)
    all_genome = set()
    with open(franks) as fd:
        data = csv.reader(fd, delimiter="\t")
        header = next(data)
        for row in data:
            all_genome.add(row[0])

    print("Load signatures...", file=sys.stderr, flush=True)
    signs, genome, breaks = signatures.load(fsignatures, size)

    print("Count genes occurences in signatures...", file=sys.stderr, flush=True)
    genes_occurences = {}
    genes = []
    for signature in signs:
        for gene in signature:
            genes_occurences.setdefault(gene, 0)
            genes_occurences[gene] += 1
            genes.append(gene)

    print("Count genes occurences in signatures...", file=sys.stderr, flush=True)

    # From small count to large count.
    genes_distrib = sorted(genes_occurences.items(), key=operator.itemgetter(1), reverse=True)
    counts = []
    for gene_count in genes_distrib:
        print(gene_count[0],gene_count[1])
        counts.append(gene_count[1])

    print(len(signs),"signatures,",len(genome),"genes in those signatures,",len(all_genome),"genes in total")

    plot = seaborn.countplot(x=genes, order=[g[0] for g in genes_distrib])
    fig = plot.get_figure()
    fig.set_size_inches(40,5)
    plot.set(xlabel="Genes", ylabel="Counts",
            title="Distribution of {ngenes}/{totalgenes} gene occurences in all signatures from\n{files}".format(
                ngenes=len(genome),
                totalgenes=len(all_genome),
                files=fsignatures))
    plot.tick_params(axis='x', rotation=45, labelsize=3)
    fig.savefig(fplot, dpi=600)


    nth_firsts = 300
    plot = seaborn.countplot(x=genes[0:nth_firsts], order=[g[0] for g in genes_distrib[0:nth_firsts]])
    fig = plot.get_figure()
    fig.set_size_inches(40,5)
    plot.set(xlabel="Genes", ylabel="Counts",
            title="Distribution of {ngenes}/{totalgenes} gene occurences in all signatures from\n{files}".format(
                ngenes=len(genome),
                totalgenes=len(all_genome),
                files=fsignatures))
    plot.tick_params(axis='x', rotation=90, labelsize=3)
    fig.savefig(fplot+".firsts.png", dpi=600)
