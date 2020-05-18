#!/bin/sh -x

K=$1
TAU=$2
INPUT=$(basename -- "$3")
NAME="${INPUT%.*}_k_${K}_tau_${TAU}.${INPUT##*.}"

./utils/generate_input.sh $K $TAU $3

./heuristic/heuristic $K $TAU hashmark_input/$NAME sensitive_pattern/$NAME

./ilp/build/HM_ilp $K $TAU hashmark_input/$NAME sensitive_pattern/$NAME

export LD_LIBRARY_PATH=/lib:/usr/lib:/usr/local/lib
./utils/pkdd $3 sensitive_pos/$NAME $K $TAU output/$NAME.cost output/$NAME.weight $TAU sensitive_pattern/$NAME output/$NAME.result_pkdd 
