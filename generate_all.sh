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

generate_TRU(){
default_input=data/original_data/truck_char.txt
name=truck_char
list_input=()
default_k=3
list_k=(2 3 4 5)
default_tau=20
list_tau=(5 10 20 30)
default_S=30
list_S=(10 30 50 70)
default_y=150

generate_all
}

generate_MSN(){
default_input=data/original_data/msnbc_char.txt
name=msnbc_char
list_input=()
default_k=8
list_k=(3 4 6 8)
default_tau=200
list_tau=(100 150 200 300)
default_S=240
list_S=(60 120 240 480)
default_y=400

generate_all
}

generate_DNA(){
default_input=data/original_data/dna_char.txt
name=dna_char
list_input=()
default_k=11
list_k=(9 11 13 15)
default_tau=20
list_tau=(5 10 20 30)
default_S=50
list_S=(30 40 50 60)
default_y=35

generate_all
}

generate_DNA
