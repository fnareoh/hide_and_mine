#!/bin/sh -x

for NAME in olden dna msnbc truck
do
for FILE in data/results/${NAME}_char.summary
do
    F="$(basename -- $FILE)"
    python3 extra/plot/plot_bar.py $FILE k ghosts
    python3 extra/plot/plot_bar.py $FILE k distortion
done

for FILE in data/results/${NAME}*char*k_?.summary
do
   if [ -f "$FILE" ]; then
   F="$(basename -- $FILE)"
   python3 extra/plot/plot_bar.py $FILE tau ghosts
   python3 extra/plot/plot_bar.py $FILE tau distortion
   fi
done

for FILE in data/results/${NAME}*char*k_??.summary
do
   if [ -f "$FILE" ]; then
   F="$(basename -- $FILE)"
   python3 extra/plot/plot_bar.py $FILE tau ghosts
   python3 extra/plot/plot_bar.py $FILE tau distortion
   fi
done

for FILE in data/results/${NAME}*char*k*tau*.summary
do
   F="$(basename -- $FILE)"
   python3 extra/plot/plot_bar.py $FILE nb_sensitive ghosts
   python3 extra/plot/plot_bar.py $FILE nb_sensitive distortion
done
done

python3 extra/plot/plot_time.py data/results/random20M_char_k_5_tau_10.summary nb_sensitive
python3 extra/plot/plot_time.py data/results/random20M_char_k_5.summary tau
python3 extra/plot/plot_time.py data/results/random20M_char.summary k
python3 extra/plot/plot_time.py data/results/random.summary n


#for i in data/figures/*.svg;do rsvg-convert -f pdf -o ${i%.*}.pdf $i;done
#pdftk data/figures/*.pdf output figures.pdf
rm data/figures/*.pdf
