#!/usr/bin/env python3

import scipy
import numpy
import matplotlib.pyplot as plt

from signatures import *

if __name__=="__main__":
    import sys
    import csv

    assert(len(sys.argv) == 4)

    franks = sys.argv[1]
    fcorrs = sys.argv[2]
    fout = sys.argv[3]

    print("Recover number of cells...", file=sys.stderr, flush=True)
    ncells = 0
    with open(franks) as fd: # FIXME this should really be an AnnData
        data = csv.reader(fd, delimiter="\t")
        header = next(data)
        ncells = len(header) - 1 # minus the "GENE" tag

    print("Load correlations...", file=sys.stderr, flush=True)
    correlations = numpy.load(fcorrs)

    print("Compute p-values from correlations...", file=sys.stderr, flush=True)
    dist = scipy.stats.beta( ncells/2 -1, ncells/2 - 1, loc=-1, scale=2)
    pvalues = 2 * dist.cdf(-numpy.abs(correlations))

    print("Save p-values...", correlations.shape, "to `", fout, "`...", file=sys.stderr, flush=True)
    numpy.save(fout, pvalues)
    print(pvalues, file=sys.stderr, flush=True)

    print("Done", file=sys.stderr, flush=True)

