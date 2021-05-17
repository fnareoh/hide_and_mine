generate(){
  local k=$1
  local tau=$2
  local file=$3
  local y=$4
  local S=$5 #number of sensitve pattern
  ./extra/generate_input.sh $k $tau $file $tau $y $S
  ./extra/test.sh $k $tau $file $S
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
    generate $default_k $default_tau $default_input $default_y ${list_S[$i_S]}
  done
  cat data/results/${name}_k_*_tau_${default_tau}_m_${default_S}.txt.comparison > data/results/${name}.summary
  cat data/results/${name}_k_${default_k}_tau_*_m_${default_S}.txt.comparison > data/results/${name}_k_${default_k}.summary
  cat data/results/${name}_k_${default_k}_tau_${default_tau}_m_*.txt.comparison > data/results/${name}_k_${default_k}_tau_${default_tau}.summary

}

generate_SYN(){
default_input=data/original_data/random20M_char.txt
name=random20M_char
list_input=("data/original_data/random5M_char.txt" "data/original_data/random10M_char.txt" "data/original_data/random15M_char.txt" "data/original_data/random20M_char.txt")
default_k=5
list_k=(3 4 5 6)
default_tau=10
list_tau=(5 10 15 20)
default_S=100
list_S=(10 100 500 1000)
default_y=100000

generate_all
generate_list_input
}

generate_OLD(){
default_input=data/original_data/olden_char.txt
name=olden_char
list_input=()
default_k=6
list_k=(3 4 5 6)
default_tau=10
list_tau=(3 5 10 15)
default_S=120
list_S=(60 120 240 320)
default_y=500

generate_all
}

generate_OLD
