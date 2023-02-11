#!/bin/sh

module load Python/3.8.1
module load snakemake

mkdir -p slurmout

sbatch --job-name "PARIS" snakemake --configfile config.yaml --jobs 15000

