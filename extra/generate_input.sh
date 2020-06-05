#!/bin/sh -x

K=$1
TAU=$2
INPUT=$(basename -- "$3")
X=$4
Y=$5
MAX_PATTERN=$6
NAME="${INPUT%.*}_k_${K}_tau_${TAU}_m_$MAX_PATTERN.${INPUT##*.}"

python extra/generate_sensitive.py $K $TAU data/original_data/$INPUT $X $Y data/sensitive_pattern/$NAME $MAX_PATTERN
python extra/find_sen_pos.py data/original_data/$INPUT data/sensitive_pattern/$NAME

echo $X > data/parameters/$NAME
echo $Y >> data/parameters/$NAME
wc -l < data/sensitive_pattern/$NAME >> data/parameters/$NAME
wc -l < data/sensitive_pos/$NAME >> data/parameters/$NAME
