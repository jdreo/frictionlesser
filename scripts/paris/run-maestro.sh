#!/bin/sh

module load Python/3.8.1
module load snakemake

sbatch --job-name "PARIS" snakemake --cluster "sbatch --mem 8000 --cpus-per-task=1 --partition common,dedicated --qos fast" --configfile config.yaml --jobs 15000

