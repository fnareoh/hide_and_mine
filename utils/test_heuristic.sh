#!/bin/sh -x

K=$1
TAU=$2
INPUT=$(basename -- "$3")
NAME="${INPUT%.*}_k_${K}_tau_${TAU}.${INPUT##*.}"

./heuristic/heuristic $K $TAU hashmark_input/$NAME sensitive_pattern/$NAME


echo " ------------- minimize_max_unfrequent_distance_to_tau ------------- "
python utils/count_losts_ghost.py $K $TAU original_data/$INPUT sensitive_pattern/$NAME output/$NAME.output_minimize_max_unfrequent_distance_to_tau

echo " ------------- minimize_sum_unfrequent_distance_to_tau ------------- "
python utils/count_losts_ghost.py $K $TAU original_data/$INPUT sensitive_pattern/$NAME output/$NAME.output_minimize_sum_unfrequent_distance_to_tau
