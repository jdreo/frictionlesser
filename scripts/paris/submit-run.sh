#!/bin/sh

size=${1}
seed=${2}

FRICTIONLESSER="../../build/app/frictionlesser"

${FRICTIONLESSER} --ranks=data/inter/ranks.tsv --cache-transcriptome=cache/trans.cache.dat --cache-size=cache/size_${size}.cache.dat --ngenes=${size} --seed=${seed} --save-sol=data/inter/frictionlesser_solutions_z${size}_s${seed}.tsv > data/output/signature_of_${size}-genes/signature_${seed} 2> data/inter/logs/frictionlesser_z${size}_s${seed}.log

