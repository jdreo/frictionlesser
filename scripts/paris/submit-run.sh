#!/bin/sh

size=${1}
seed=${2}

FRICTIONLESSER="../../release/app/frictionlesser"

${FRICTIONLESSER} --ranks=ranks.tsv --cache-transcriptome=cache/trans.cache.dat --cache-size=cache/size_${size}.cache.dat --ngenes=${size} --seed=${seed} > signature_of_${size}-genes/signature_${seed} 2> logs/z${size}_s${seed}.log

