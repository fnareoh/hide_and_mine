#!/bin/sh -x

K=$1
TAU=$2
INPUT=$(basename -- "$3")
NAME="${INPUT%.*}_k_${K}_tau_${TAU}.${INPUT##*.}"

./ilp/build/HM_ilp $K $TAU hashmark_input/$NAME sensitive_pattern/$NAME

echo " ------------- ILP ------------- "
python utils/count_losts_ghost.py $K $TAU original_data/$INPUT sensitive_pattern/$NAME output/$NAME.output_ILP
