#!/bin/sh -x

K=$1
TAU=$2
INPUT=$(basename -- "$3")
NAME="${INPUT%.*}_k_${K}_tau_${TAU}.${INPUT##*.}"
X=$4
Y=$5

python utils/generate_sensitive.py $K $TAU original_data/$INPUT $X $Y sensitive_pattern/$NAME

python utils/find_sen_pos.py original_data/$INPUT sensitive_pattern/$NAME
