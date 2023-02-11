#!/bin/sh

runs=5000

#SBATCH --natsks=15000
#SBATCH --cpu-per-task=1

mkdir -p logs

FRICTIONLESSER="../../release/app/frictionlesser"

for size in 10 20 30; do
	mkdir -p signature_of_${size}-genes
done


for seed in $(seq $runs); do
  for size in 10 20 30; do
    srun --quiet --job-name z${size}_s${seed} --mem 16G --cpus-per-task=1 --partition common --qos normal --output signature_of_${size}-genes/signature_${seed} --error logs/%j.log     ${FRICTIONLESSER} --ranks=ranks.tsv --cache-transcriptome=cache/trans.cache.dat --cache-size=cache/size_${size}.cache.dat --ngenes=${size} --seed=${seed} &
  done
done
# wait for all the steps to end before terminating the sbatch script
wait

