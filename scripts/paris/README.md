
This pipeline targets finding ovarian cancer signatures in transcriptomics data from the PARIS project.
https://www.parisproject.org/


Introduction
============

The FRICTIONLESS method searches for fixed-size signatures
(i.e. a fixed number of genes) and use the Friedman statistics across *samples*.

The pipeline is split in three parts that must be ran manually:
1. pre-processing of count matrix,
2. signature search,
3. post-processing and validation.

Most of the pipeline is written in Python.
Most of the data structures are handled with the AnnData module.

Input data should be put in `data/input`, the pipelines put intermediate data in
`data/inter` and end results in `data/output`.
Quality checks data and plots are in `data/inter/qc`.

The script `expe_build.sh <prefix_of_input_data>`
may be copied in a separate directory and run to build
up what’s needed to run the whole experiment.

The main outputs are two CSV files:
- `data/output/correlations_signatures-genes.csv`.
- `data/output/scores_signatures-samples.csv`


Pre-processing pipeline
=======================

The pipeline is handled by the `preproc.Snakefile` SnakeMake script.

Starting for the count matrix produced by the project, this step mainly filter
the count matrix, compute ranks and pre-compute caches for the search.

The main output are a space-separated text file of ranks: `ranks.tsv`
and a compressed HDF5 file with counts, ranks and meta-data: `paris+ranks.h5an.gz`,
of the form:

       ←  all genes  →
     ┌────────────────┐
     │ var:"id",[…]   │
     └────────────────┘
     ┏━━━━━━━━━━━━━━━━┓  ┌────────────────┐
     ┃ X:             ┃  │ obs:           │  ↑
     ┃ counts         ┃┓ │ "sample",[…]   │ cells
     ┃                ┃┃ │                │  ↓
     ┗━━━━━━━━━━━━━━━━┛┃ └────────────────┘
      ┃layers["ranks"] ┃ 
      ┗━━━━━━━━━━━━━━━━┛ 
        ←  all genes  →

All data are put in the intermediate data directory: `data/inter`.


Filtering
---------

The filtering pipepline uses the scanpy module.

0. Starting from the `counts.npz` matrix (cells ×_genes),
   and using `features.csv` meta-data about genes and
   `meta.csv` meta-data about cells.
1. Filter out non-cancer cells.
2. Normalize total counts,
3. Take the logarithm.
4. Keep only the 10 000 highly-variable genes.


Ranking & Cache
---------------

Ranks are computed using the `frictionlesser` software and put in both a CSV
file and an AnnData h5an file.

The pipeline also pre-computes caches:
- for the input rank matrix,
- for all the signature sizes listed in the `config_preproc.yaml` configuration file.


Signature Search
================

The search is composed of approximately 10 0000 independent runs of the
frictionlesser executable, with different random seeds.

Because it exceed the forking capacity of SnakeMake, it is handled by a simple
shell script: `submit-expe_10.sh <first_seed> <nb_of_runs>`.
Seeds starts fro_`first_seed` and are then incremented by one until `nb_of_runs`
are submitted.

This script submits `expe_run.sh` on a SLURM cluster
(with the environment being handled by the `module` software).

The main output is 10 000 files, each one being one signature of the form:
`<global_score> <nb_of_samples> <score_sample_0> … <score_sample_n> <nb_of_genes> <gene_0> … <gene_m>`.
`gene` is the actual gene label, as used in the input data.
Those files are in `data/output/signatures of_XX-genes` with XX the size of the signatures.


Post-processing & validation
============================

The pipeline is handled by the `validation.Snakefile` SnakeMake script.

The current post-processing uses hard-coded signature size of 10 genes.

The pipeline mainly computes:
1. Correlations between genes.
2. Correlations between signatures.
3. Correlations between signatures and genes.
4. Signatures-samples scores.

The main outputs are two CSV files:
- `data/output/correlations_signatures-genes.csv`.
- `data/output/scores_signatures-samples.csv`


Genes correlations
------------------

This computes the average (across samples) correlations of ranks between genes
and the corresponding p-values.
NOTE: some genes may have no rank variability,
in which case their correlation is manually set to zero.

The main output of this step is an AnnData file
with only the genes (columns) that are actually in one of the
signatures: `gccorr.h5an.gz`.

        ←  genes  →
     ┌────────────────┐
     │ var:"id",[…]   │
     └────────────────┘
     ┏━━━━━━━━━━━━━━━━┓  ┌────────────────┐
     ┃ X:             ┃  │ obs:           │  ↑
     ┃ counts         ┃┓ │ "sample",[…]   │ cells
     ┃                ┃┃ │                │  ↓
     ┗━━━━━━━━━━━━━━━━┛┃┓└────────────────┘
      ┃layers["ranks"] ┃┃
      ┗━━━━━━━━━━━━━━━━┛┃
       ┃layers["zscore"]┃
       ┗━━━━━━━━━━━━━━━━┛
        ←  genes  →
     ┌────────────────┐
     │ varp:          │┐
     │ "correlations" ││  ↑
     │                ││ genes
     │                ││  ↓
     │                ││
     └────────────────┘│
      │ "p-values"     │
      └────────────────┘
         ←  genes  →

Note that those data are actually computed for all genes, and inserted back in
the `paris+ranks.h5an.gz` file.


Signatures correlations
-----------------------

This computes the average (across samples) correlations of ranks between signatures
and the corresponding p-values.

The main output of this step is an AnnData file
with only the signatures (columns) that are actually in one of the
signatures: `sccorr.h5an.gz`.

       ← signatures →
     ┌─────────────────┐
     │ var:"signature" │
     │     "origin"    │
     └─────────────────┘
     ┏━━━━━━━━━━━━━━━━━┓  ┌────────────────┐
     ┃ X:              ┃  │ obs:           │  ↑
     ┃sum.zscore/ngenes┃  │ "sample"       │ cells
     ┃                 ┃  │                │  ↓
     ┗━━━━━━━━━━━━━━━━━┛  └────────────────┘
     ┌─────────────────┐
     │ varp:           │
     │ "correlations"  │┐  ↑
     │                 ││ signatures
     │                 ││  ↓
     │                 ││
     └─────────────────┘│
      │ "p-values"      │
      └─────────────────┘
       ← signatures →


Signatures-genes correlations
-----------------------------

This computes the average (across samples) correlations of ranks between
signatures and genes.

The main output of this step is an AnnData file
with only the signatures (columns) that are actually in one of the
signatures `sgcorr.h5an.gz`, and the corresponding CSV file:
`data/output/correlations_signatures-genes.csv`.

       ← all genes  →
     ┌───────────────┐
     │ var:"id", […] │
     └───────────────┘
     ┏━━━━━━━━━━━━━━━┓  ┌────────────────┐
     ┃ X:            ┃  │ obs:           │  ↑
     ┃ correlations  ┃  │ "signature",   │ signatures
     ┃               ┃  │ "origin"       │  ↓
     ┗━━━━━━━━━━━━━━━┛  └────────────────┘
       ← all genes  →


Signatures-samples scores
-------------------------

This step is a simple recollection of the data within signatures files,
in a single matrix, as a CSV file: `data/output/scores_signatures-samples.csv`
