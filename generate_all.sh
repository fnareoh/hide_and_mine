generate(){
  status=1
  while [ $status -ne 0 ]; do
    local k=$1
    local tau=$2
    local file=$3
    local y=$4
    local S=$5 #number of sensitve pattern
    ./extra/generate_input.sh $k $tau $file $tau $y $S
    ./extra/test.sh $k $tau $file $S
    status=$?
  done
}

generate_list_input(){
  for i_input in ${!list_input[@]}
  do
    generate $default_k $default_tau ${list_input[$i_input]} $default_y $default_S
  done
  cat data/results/${generic_name}*_k_${default_k}_tau_${default_tau}_m_${default_S}.txt.comparison > data/results/${generic_name}.summary
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

generate_SYN(){
default_input=data/original_data/random5M_char.txt
name=random5M_char
generic_name=random
list_input=("data/original_data/random5M_char.txt" "data/original_data/random10M_char.txt" "data/original_data/random15M_char.txt" "data/original_data/random20M_char.txt")
default_k=6
list_k=(3 4 5 6)
default_tau=7
list_tau=(5 7 10 12 15)
default_S=700
list_S=(10 100 500 650 700 1000)
default_y=100000
list_y=(100000 100000 100000 100000 100000 100000)

generate_all
generate_list_input

}

generate_OLD(){
default_input=data/original_data/olden_char.txt
name=olden_char
list_input=()
default_k=6
list_k=(3 5 6 7)
default_tau=5
list_tau=(3 5 10 15)
default_S=240
list_S=(60 120 240 480)
default_y=1000
list_y=(1000 1000 1000 1000)

generate_all
}

generate_TRU(){
default_input=data/original_data/truck_char.txt
name=truck_char
list_input=()
default_k=3
list_k=(2 3 4 5)
default_tau=5
list_tau=(3 5 7 10)
default_S=120
list_S=(40 80 120 160)
default_y=250
list_y=(250 250 250 250)

generate_all
}

generate_MSN(){
default_input=data/original_data/msnbc_char.txt
name=msnbc_char
list_input=()
default_k=8
list_k=(4 6 8 10)
default_tau=200
list_tau=(50 100 200 300)
default_S=300
list_S=(100 200 300 600)
default_y=500
list_y=(500 500 500 500)

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
list_y=(35 35 35 35)

generate_all
}

./compile.sh
generate_OLD
generate_TRU
generate_MSN
generate_DNA
generate_SYN
