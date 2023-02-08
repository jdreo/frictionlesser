
Frictionless
============

*FRIedman transCripTomic signatures Occuring iN muLtiplE Samples Simultaneously*.

Emile Zakiev, Johann Dreo, Benno Schwikowski.

This software aims at finding genetic "signatures" in RNA sequencing data.
A signature is hereby define as a set of genes that are "differentially" expressed
across a set of samples (bulk or single cells).

The problem is modelled as a multi-modal maximization
of a partition on the multi-dimensional boolean space
(of genes).
It uses partial evaluation for the stable neighborhood
(i.e. swaping genes, so having a fixed-dimensional partition).

The search algorithm is implemented using the
[ParadisEO](https://github.com/nojhan/paradiseo) framework.


Build
=====

A self-contained [Apptainer](https://apptainer.org/) container definition file is provided.
To build the frictionlesser binary, with apptainer installed, you can just do:
```sh
apptainer build --fakeroot frictionlesser.sif frictionlesser.def
```

You can follow the commands of the `%post` section of the `frictionlesser.def`
file to see how to build locally.

NOTE: when building in Debug mode, a lot more checks are performed.
It would thus be a good idea to perform at least one run with a binary built in Debug mode,
to double check that the input data are consistent.


Usage
=====

To run the container as an executable:
```sh
apptainer run frictionlesser.sif --help
```

An example pipeline showing how to use frictionlesser is provided as a
[SnakeMake](https://snakemake.readthedocs.io/) file: `scripts/test-Neftel/Snakefile`.

The executable (raw or containerized) expects at least a ranked table of
transcriptomics expressions levels following Neftel et al. format (see below).

You can generate a ranked table from raw expression table with:
```sh
frictionlesser --exprs=<expression_table.tsv> > ranked.tsv
```

NOTE: frictionlesser always output space-separated text tables,
but it can actually read CSVs, TSVs or any such text table format using spaces,
comma or semicolon as separator.


Input data
----------

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


Actual search
-------------

### Basic search

To perform a single run, the simplest way is to indicate where to find the ranked table:
```sh
frictionlesser --ranks=<ranked_table.csv> > result.txt
```

Frictionlesser will print on the standard output the found signature,
following the format: "<score> <nb_of_genes> <gene_0> … <gene_n>".


## Caches

Frictionlesser can save pre-computed caches, to speed-up parallel computations:

- one for a given transcription levels table (option `--cache-transcriptome`),
- one for a given signature size for this transcription table (option `--cache-size`).

The option `--cache-only` allows to precompute caches without performing the search,
allowing for task decoupling in a pipeline (see example pipelines in `scripts`).

To precompute the transcriptome cache only:
```sh
frictionlesser --ranks=<ranked_table.csv> --cache-transcriptome=<trans-cache.dat> --cache-only
```

To precompute the size cache (and the transcriptome one if it is not already computed, else load it):
```sh
frictionlesser --ranks=<ranked_table.csv> --cache-transcriptome=<trans-cache.dat> --cache-size=<size-cache.dat> --cache-only
```

To then load the caches and perform a run:
```sh
frictionlesser --ranks=<ranked_table.csv> --cache-transcriptome=<trans-cache.dat> --cache-size=<size-cache.dat>
```

Examples
========

An example pipeline searching for Glioblastoma signatures is provided in the
`scripts/test-Neftel` directory.

It is self-contained: it downloads the Neftel et al. [1] data by itself and then
perform a set of searches, after having computed caches.


Sources
=======

Headers are in `include/`, sources in `src/`, executable in `app/`, examples and tests in `tests/`.

The build generates `libfrictionless` and `libR`.


References
==========

[1] Neftel et al., An Integrative Model of Cellular States, Plasticity, and Genetics for Glioblastoma,
    Cell. 2019 Aug 8; 178(4):835-849.e21,
    https://doi.org/10.1016/j.cell.2019.06.024
    Related datasets: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131928
