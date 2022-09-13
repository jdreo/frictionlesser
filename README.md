
Frictionless
============

*FRIedman transCripTomic signatures Occuring iN muLtiplE Samples Simultaneously*.

Emile Zakiev, Johann Dreo, Benno Schwikowski.

This software aims at finding genetic "signatures" in RNA sequencing data.
A signature is hereby define as a set of genes that are "differentially" expressed
across a set of samples (bulk or single cells).


Build
=====

You must first download the Paradiseo framework somewhere on your computer and build it.

Then use a classical CMake workflow.

```sh
mkdir build
cd build
cmake -DPARADISEO_ROOT=~/paradiseo -DPARADISEO_BUILD=~/paradiseo/build
make
```

Usage
=====

The executable is `app/frictionlesser`, you will need to feed it with either a raw expression table or a ranked table.
Expression tables are expected to follow the following space-separated table format.

|  GENE   |  Cell_0 | … | Cell_i | … | Cell_m |
|---------|---------|---|--------|---|--------|
| gene_0  |  r_00   | … |  r_i0  | … |  r_n0  |
|    …    |    …    | … |   …    | … |   …    |
| gene_j  |  r_0j   | … |  r_ij  | … |  r_mj  |
|    …    |    …    | … |   …    | … |   …    |
| gene_n  |  r_0n   | … |  r_in  | … |  r_mn  |

NOTES:

- The first cell of the first row may be anything, but a warning will be printed if it is not "GENE".
- Cells should be named from their sample name, followed by an ID of the form "-[A-Z][0-9]{2}$",
e.g. "-A01" (following Neftel et al. convention [1]).


[1] Neftel et al., An Integrative Model of Cellular States, Plasticity, and Genetics for Glioblastoma,
    Cell. 2019 Aug 8; 178(4):835-849.e21,
    https://doi.org/10.1016/j.cell.2019.06.024
    Related datasets: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131928
