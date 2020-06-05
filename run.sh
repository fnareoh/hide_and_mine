#!/bin/sh -x

K=$1
TAU=$2
HASH_INPUT=$3
SENS_INPUT=$4
NAME="$(basename -- $HASH_INPUT)"

export LD_LIBRARY_PATH=/lib:/usr/lib:/usr/local/lib

mkdir -p data/results
mkdir -p data/output

./src/heuristic/heuristic $K $TAU $HASH_INPUT $SENS_INPUT

./src/ilp/build/HM_ilp $K $TAU $HASH_INPUT $SENS_INPUT

python extra/evaluate_all.py $K $TAU $HASH_INPUT $SENS_INPUT $NAME
