#!/bin/bash

module load Python/3.8.1
module load snakemake/7.16.1

cd frictionlesser/scripts/paris/

snakemake --cores 1 --snakefile preproc.Snakefile

