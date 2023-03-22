import itertools
import numpy
import sys
import signatures
import csv

if __name__ == "__main__":

    rankfile = sys.argv[1]
    fgenescorr = sys.argv[2]
    size = int(sys.argv[3])
    signatures_files = sys.argv[4:6]

    print("Load genes correlations...", file=sys.stderr, flush=True)
    genes_correlations = numpy.load(fgenescorr)
    ngenes = genes_correlations.shape[0]
    assert(genes_correlations.shape[1] == ngenes)

    print("Load signatures...", file=sys.stderr, flush=True)
    genesets_1,genome_1,breaks_1 = signatures.load([signatures_files[0]], filter_size=size)
    assert(len(genesets_1) > 0)
    assert(len(genome_1) > 0)

    genesets_2,genome_2,breaks_2 = signatures.load([signatures_files[1]], filter_size=size)
    assert(len(genesets_2) > 0)
    assert(len(genome_2) > 0)

    print("Load genes indices...", file=sys.stderr, flush=True)
    genome = genome_1 | genome_2
    genes_idx = {}
    ranks_l = []
    with open(rankfile) as fd:
        csvreader = csv.reader(fd, delimiter=" ")
        header = next(csvreader)
        n=0
        for row in csvreader:
            gene = row[0]
            assert(len(gene) > 0)
            if gene in genome:
                genes_idx[gene] = n
                n += 1
                print("\r",n,end=" ", file=sys.stderr,flush=True)
        assert(len(genes_idx) > 0)
    print("genes", file=sys.stderr, flush=True)

    print("Compute paired signatures correlations...", file=sys.stderr, flush=True)
    nsign_1 = len(genesets_1)
    nsign_2 = len(genesets_2)
    signature_correlations = numpy.zeros((nsign_1, nsign_2))
    for i1,s1 in enumerate(genesets_1):
        for i2,s2 in enumerate(genesets_2):
            g1 = [genes_idx[g] for g in s1 if g in genes_idx]
            assert(len(g1) > 0)
            g2 = [genes_idx[g] for g in s2 if g in genes_idx]
            assert(len(g2) > 0)
            subidx = itertools.product(g1,g2)
            i,j = zip(*subidx)
            submat = genes_correlations[i,j]
            assert(submat.size == len(g1)*len(g2))
            avcorr = numpy.mean(submat)
            assert(avcorr > 0)
            signature_correlations[i1,i2] = avcorr

    fname = "signatures-average-correlations.npy"
    print("Save signatures correlation matrix of shape", signature_correlations.shape, "to `", fname,"`...", file=sys.stderr, flush=True)
    numpy.save(fname, signature_correlations)

    print("Compute all signatures correlations...", file=sys.stderr, flush=True)
    genesets = genesets_1 | genesets_2
    nsign = len(genesets)
    all_signature_correlations = numpy.zeros((nsign, nsign))
    signature_origins = {}
    for i1,s1 in enumerate(genesets):
        if s1 in genesets_1 and s1 in genesets_2:
            signature_origins[i1] = "BOTH"
        elif s1 in genesets_1:
            signature_origins[i1] = signatures_files[0]
        elif s1 in genesets_2:
            signature_origins[i1] = signatures_files[1]
        else:
            assert(s1 in genesets_1 or s1 in genesets_2)

        for i2,s2 in enumerate(genesets):
            # if i2 < i1:
            #     continue
            g1 = [genes_idx[g] for g in s1 if g in genes_idx]
            assert(len(g1) > 0)
            g2 = [genes_idx[g] for g in s2 if g in genes_idx]
            assert(len(g2) > 0)
            # Cross-product of all genes indices.
            subidx = itertools.product(g1,g2)
            # Corresponding tsb-matrix indices.
            i,j = zip(*subidx)
            submat = genes_correlations[i,j]
            assert(submat.size == len(g1)*len(g2))
            # Average correlation over the sub-matrix.
            # TODO would be faster over the triangular sub-matrix.
            avcorr = numpy.mean(submat)
            assert(avcorr > 0)
            all_signature_correlations[i1,i2] = avcorr

    assert(len(signature_origins) == nsign)
    fname = "signatures-origins.csv"
    print("Save unique signatures origin table to `", fname, "`...", file=sys.stderr, flush=True)

    with open(fname, 'w') as fd:
        writer = csv.DictWriter(fd, fieldnames=["gene_index","origin"])
        writer.writeheader()
        for i in signature_origins:
            writer.writerow({"gene_index":i, "origin":signature_origins[i]})

    fname = "signatures-average-correlations_full.npy"
    print("Save full signatures correlation matrix of shape", all_signature_correlations.shape, "to `", fname,"`...", file=sys.stderr, flush=True)
    numpy.save(fname, all_signature_correlations)

    print("Done", file=sys.stderr, flush=True)
