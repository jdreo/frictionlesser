#!/bin/bash

# Run experiment from scratch on Pasteur's maestro.

# tag="v0.2_PARIS_expe_2"

if [[ -z ${1-} || $1 == --help || $1 == -h ]]; then
    printf 'Usage: %s input_data_prefix

Expects 3 data files: PATH/TO/XXX_counts.npz, PATH/TO/XXX_features.csv, PATH/TO/XXX_meta.csv
with input_data_prefix being "PATH/TO/XXX_".
For example: %s PATH/TO/XXX_
' "$0" "$0"
    exit 1
fi

if [[ ${PWD##*/} == "paris" ]]; then
    printf "You should not run the build script from within the source directory, but rather create an experiment directory, copy the build script in it and run it from there.\n"
    exit 2
fi

printf "Check input data...\n"
input_data_prefix="$1"

infiles=(counts.npz features.csv meta.csv)
for infile in ${infiles[@]}; do
    if [[ ! -f ${input_data_prefix}${infile} ]]; then
        printf "ERROR: missing data files using prefix: \"${input_data_prefix}${infile}\"\n"
        exit 3
    fi
done

printf "Build softwares...\n"
module load cmake/3.19.4 gcc/10.4.0 OpenBLAS/0.3.21

printf "Build Paradiseo...\n"
# Get Paradiseo from Pasteur's Gitlab, to get the paris tag.
git clone --branch master --single-branch --recurse-submodules https://gitlab.pasteur.fr/jdreo/paradiseo.git
cd paradiseo
#git switch -
git pull --tags
#git checkout tags/${tag}
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON ..
make
cd ../..

printf "Build Frictionlesser...\n"
# Get Frictionlesser from Pasteur's Gitlab, to get the paris tag.
git clone --branch main --single-branch --recurse-submodules https://gitlab.pasteur.fr/jdreo/frictionlesser.git
cd frictionlesser
#git switch -
git pull --tags
#git checkout tags/${tag}
mkdir -p build
cd build
cmake -DBUILD_DOCS=OFF -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Release -DUSE_LOCAL_PARADISEO=ON -DPARADISEO_ROOT=../../paradiseo -DPARADISEO_BUILD=../../paradiseo/build -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON ..
make

printf "Install scipy and scanpy...\n"
module load Python/3.8.1
module load snakemake/7.16.1
module load hdf5/1.10.6
pip3 install --user Cython
pip3 install --user scipy scanpy sortedcollections fastcluster

printf "Copy input data in working directory...\n"
cd ../../
cp ${input_data_prefix}counts.npz   frictionlesser/scripts/paris/data/input/counts.npz
cp ${input_data_prefix}features.csv frictionlesser/scripts/paris/data/input/features.csv
cp ${input_data_prefix}meta.csv     frictionlesser/scripts/paris/data/input/meta.csv

printf "Link run script...\n"
ln -s frictionlesser/scripts/paris/expe_run.sh

printf "Done.\nNow run: expe_run.sh\n"
