#!/bin/sh

seedbase="$1"
runs="$2"

mkdir -p logs

FRICTIONLESSER="../../release/app/frictionlesser"

size=10
#for size in 10; do
  mkdir -p signature_of_${size}-genes
#done

for r in $(seq $runs); do
    seed=$((seedbase+r))
    echo "Submission $r: z${size}_s${seed}"
    sbatch --job-name=z${size}_s${seed} --mem=8000 --partition common,dedicated --qos fast ./submit-run.sh ${size} ${seed}
done

