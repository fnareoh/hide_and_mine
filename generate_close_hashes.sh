generate(){
  local k=$1
  local tau=$2
  local file=$3
  local y=$4
  local S=$5 #number of sensitve pattern
  local input=$(basename -- "$file")
  local NAME="${input%.*}_k_${k}_tau_${tau}_m_$S.txt"

  #generate input
  python extra/fasta_to_hasmark_input.py $file data/hashmark_input/$NAME
  touch data/sensitive_pattern/$NAME

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
  python extra/evaluate_all.py $k $tau data/hashmark_input/$NAME data/sensitive_pattern/$NAME $NAME close_hashes
}

generate_list_input(){
  for i_input in ${!list_input[@]}
  do
    generate $default_k $default_tau ${list_input[$i_input]} $default_y $default_S
  done
}

generate_all() {
  for i_k in ${!list_k[@]}
  do
    generate ${list_k[$i_k]} $default_tau $default_input $default_y $default_S
  done
  for i_tau in ${!list_tau[@]}
  do
    generate $default_k ${list_tau[$i_tau]} $default_input $default_y $default_S
  done
  for i_S in ${!list_S[@]}
  do
    generate $default_k $default_tau $default_input ${list_y[$i_S]} ${list_S[$i_S]}
  done
  cat data/results/${name}_k_*_tau_${default_tau}_m_${default_S}.txt.comparison > data/results/${name}.summary
  cat data/results/${name}_k_${default_k}_tau_*_m_${default_S}.txt.comparison > data/results/${name}_k_${default_k}.summary
  cat data/results/${name}_k_${default_k}_tau_${default_tau}_m_*.txt.comparison > data/results/${name}_k_${default_k}_tau_${default_tau}.summary

}

generate_close_hashes(){
default_input=data/original_data/P7_reads.fa
name=P7_reads
list_input=()
default_k=11
list_k=(9 11 13 15)
default_tau=20
list_tau=(5 10 20 30)
default_S=0
list_S=()
default_y=35
list_y=(35 35 35 35)

generate_all
}

generate_close_hashes
