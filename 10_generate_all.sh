for i in {1..10}
do
  ./generate_all.sh
  mkdir data/results${i}
  cp -r data/results/*.summary data/results${i}/
done