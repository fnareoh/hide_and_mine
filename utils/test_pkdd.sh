#!/bin/sh -x

K=$1
TAU=$2
INPUT=$(basename -- "$3")
NAME="${INPUT%.*}_k_${K}_tau_${TAU}.${INPUT##*.}"

export LD_LIBRARY_PATH=/lib:/usr/lib:/usr/local/lib
./utils/pkdd $3 sensitive_pos/$NAME $K $TAU output/$NAME.cost output/$NAME.weight $TAU sensitive_pattern/$NAME output/$NAME.knapsack_pkdd output/$NAME.output_pkdd hashmark_input/$NAME

echo " ------------- PKDD ------------- "
python utils/count_losts_ghost.py $K $TAU original_data/$INPUT sensitive_pattern/$NAME output/$NAME.output_pkdd
