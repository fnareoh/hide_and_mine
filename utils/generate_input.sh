#!/bin/sh -x

K=$1
TAU=$2
INPUT=$(basename -- "$3")
NAME="${INPUT%.*}_k_${K}_tau_${TAU}.${INPUT##*.}"

./utils/maxmotif F  -l $K -u $K -p $((TAU+1)) original_data/$INPUT $TAU sensitive_pattern/$NAME

python utils/find_sen_pos.py original_data/$INPUT sensitive_pattern/$NAME

export LD_LIBRARY_PATH=/lib:/usr/lib:/usr/local/lib

./utils/pkdd_modified original_data/$INPUT sensitive_pos/$NAME sensitive_pattern/$NAME $K $TAU > hashmark_input/$NAME
