generate(){
  local k=$1
  local tau=$2
  local file=$3
  local y=$4
  local S=$5 #number of sensitve pattern
  ./extra/generate_input.sh $k $tau $file $tau $y $S
  ./extra/test.sh $k $tau $file $S
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
  for i_input in ${!list_input[@]}
  do
    generate $default_k $default_tau ${list_input[$i_input]} $default_y $default_S
  done
}

default_input=data/original_data/random20M_char.txt
list_input=("data/original_data/random5M_char.txt" "data/original_data/random10M_char.txt" "data/original_data/random15M_char.txt" "data/original_data/random20M_char.txt")
default_k=5
list_k=(3 4 5 6)
default_tau=10
list_tau=(5 10 15 20)
default_S=100
list_S=(10 100 500 1000)
default_y=100000

generate_all
