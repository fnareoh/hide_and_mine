#!/bin/sh -x
rm data/results/*.svg
rm data/figures/*.svg

for NAME in olden dna msnbc truck
do
for FILE in data/results1/${NAME}_char.summary
do
    F="$(basename -- $FILE)"
    python extra/plot/plot_avg_bar.py $F k ghosts
    python extra/plot/plot_avg_bar.py $F k distortion
done

for FILE in data/results1/${NAME}*char*k_?.summary
do
   F="$(basename -- $FILE)"
   python extra/plot/plot_avg_bar.py $F tau ghosts
   python extra/plot/plot_avg_bar.py $F tau distortion
done

for FILE in data/results1/${NAME}*char*k_??.summary
do
   F="$(basename -- $FILE)"
   python extra/plot/plot_avg_bar.py $F tau ghosts
   python extra/plot/plot_avg_bar.py $F tau distortion
done

for FILE in data/results1/${NAME}*char*k*tau*.summary
do
   F="$(basename -- $FILE)"
   python extra/plot/plot_avg_bar.py $F nb_sensitive ghosts
   python extra/plot/plot_avg_bar.py $F nb_sensitive distortion
done
done

python extra/plot/plot_time.py data/results/random20M_char_k_5_tau_10.summary nb_sensitive
python extra/plot/plot_time.py data/results/random20M_char_k_5.summary tau
python extra/plot/plot_time.py data/results/random20M_char.summary k
python extra/plot/plot_time.py data/results/random.summary n


#for i in data/figures/*.svg;do rsvg-convert -f pdf -o ${i%.*}.pdf $i;done
#pdftk data/figures/*.pdf output figures.pdf
rm data/figures/*.pdf
