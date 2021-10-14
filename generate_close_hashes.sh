#!/bin/bash -x
generate(){
  local k=$1
  local tau=$2
  local file=$3
  local y=$4
  local S=$5 #number of sensitve pattern
  local input=$(basename -- "$file")
  local NAME="${input%.*}_k_${k}_tau_${tau}_m_$S.txt"
  k_minimal_absent=(11 12 13)
  pwd=$(pwd)

  status=1
  while [ $status -ne 0 ]; do
    #generate input
    for item in $k_minimal_absent
    do
      if [ "${k}" == "${item}" ]; then
        (shuf "${pwd}/data/sensitive_pattern/P7_reads_minimal_absent/P7_reads_k_${k}.sensitive" | head -$S ) > "${pwd}/data/sensitive_pattern/${NAME}"
       fi
    done
    touch data/sensitive_pattern/$NAME

    #ln -s "${pwd}/data/hashmark_input/${input}" "${pwd}/data/hashmark_input/${NAME}"
    python3 extra/find_sen_pos.py $file data/sensitive_pattern/$NAME
    echo "0" > data/parameters/$NAME
    echo "0" >> data/parameters/$NAME
    wc -l < data/sensitive_pattern/$NAME >> data/parameters/$NAME
    wc -l < data/sensitive_pos/$NAME >> data/parameters/$NAME

    export LD_LIBRARY_PATH=/lib:/usr/lib:/usr/local/lib
    ./extra/pkdd/main $file data/sensitive_pos/$NAME $k $tau data/output/$NAME.cost data/output/$NAME.weight $tau data/sensitive_pattern/$NAME data/output/$NAME.knapsack_pkdd data/output/$NAME.output_pkdd data/hashmark_input/$NAME

    status=$?
  done

  #compile for close hashes
  cd src/heuristic
  g++ -O3 heuristic.cpp -o heuristic -DCLOSE_HASHES
  cd ../..

  #run for close hashes
  ./src/heuristic/heuristic $k $tau data/hashmark_input/$NAME data/sensitive_pattern/$NAME
  status=$?
  if [ $status -ne 0 ]; then
    exit 1
  fi

  #evaluate
  python3 extra/evaluate_all.py $k $tau data/hashmark_input/$NAME data/sensitive_pattern/$NAME $NAME close_hashes
}

generate_list_input(){
  for i_input in ${!list_input[@]}
  do
    generate $default_k $default_tau ${list_input[$i_input]} $default_y $default_S
  done
  cat data/results/${name}_*M_k_${default_k}_tau_${default_tau}_m_${default_S}.txt.comparison > data/results/${name}_len.summary
}

generate_all() {
  generate_list_input
  for i_k in ${!list_k[@]}
  do
    generate ${list_k[$i_k]} $default_tau $default_input $default_y $default_S
  done
  for i_tau in ${!list_tau[@]}
  do
    generate $default_k ${list_tau[$i_tau]} $default_input $default_y $default_S
  done
  cat data/results/${name}_k_*_tau_${default_tau}_m_${default_S}.txt.comparison > data/results/${name}.summary
  cat data/results/${name}_k_${default_k}_tau_*_m_${default_S}.txt.comparison > data/results/${name}_k_${default_k}.summary


}

generate_close_hashes(){
#Generate limited input
python3 extra/fasta_to_hasmark_input.py data/original_data/P7_reads.fa data/hashmark_input/P7_reads.txt
head -c 1000000 data/hashmark_input/P7_reads.txt > data/hashmark_input/P7_reads_1.0M.txt
head -c 1500000 data/hashmark_input/P7_reads.txt > data/hashmark_input/P7_reads_1.5M.txt
head -c 2000000 data/hashmark_input/P7_reads.txt > data/hashmark_input/P7_reads_2.0M.txt
head -c 2500000 data/hashmark_input/P7_reads.txt > data/hashmark_input/P7_reads_2.5M.txt
default_input="data/hashmark_input/P7_reads_2.5M.txt"
name=P7_reads
list_input=("data/hashmark_input/P7_reads_1.0M.txt" "data/hashmark_input/P7_reads_1.5M.txt" "data/hashmark_input/P7_reads_2.0M.txt" "data/hashmark_input/P7_reads_2.5M.txt")
for file in $list_input
do
  grep -o '#' $file | wc -l > data/parameters/$(basename -- "$file")
done
default_k=11
list_k=(11 12 13)
default_tau=20
list_tau=(10 20 30 40)
default_S=1000
list_S=()
default_y=35
list_y=(35 35 35 35)

generate_all

figures
}

figures () {
    python3 extra/plot/plot_bar.py data/results/${name}.summary k ghosts
    #python3 extra/plot/plot_bar.py data/results/${name}.summary k distortion
    python3 extra/plot/plot_bar.py data/results/${name}_k_${default_k}.summary tau ghosts
    #python3 extra/plot/plot_bar.py data/results/${name}_k_${default_k}.summary tau distortion
}

10_generate_close_hashes(){
for i in {1..10}
do
  rm data/results/*
  generate_close_hashes
  mkdir data/results${i}
  cp -r data/results/*.summary data/results${i}/
done
}

if [[ ! -f data/original_data/P7_reads.fa ]]
then
    unzip data/original_data/P7_reads.zip -d data/original_data/
fi

#generate_close_hashes
10_generate_close_hashes

