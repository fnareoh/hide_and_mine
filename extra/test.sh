#!/bin/sh -x

K=$1
TAU=$2
INPUT=$(basename -- "$3")
MAX_PATTERN=$4
NAME="${INPUT%.*}_k_${K}_tau_${TAU}_m_$MAX_PATTERN.${INPUT##*.}"

export LD_LIBRARY_PATH=/lib:/usr/lib:/usr/local/lib
./extra/pkdd/main $3 data/sensitive_pos/$NAME $K $TAU data/output/$NAME.cost data/output/$NAME.weight $TAU data/sensitive_pattern/$NAME data/output/$NAME.knapsack_pkdd data/output/$NAME.output_pkdd data/hashmark_input/$NAME

./src/heuristic/heuristic $K $TAU data/hashmark_input/$NAME data/sensitive_pattern/$NAME

./src/ilp/build/HM_ilp $K $TAU data/hashmark_input/$NAME data/sensitive_pattern/$NAME

python extra/evaluate_all.py $K $TAU data/original_data/$INPUT data/sensitive_pattern/$NAME $NAME
