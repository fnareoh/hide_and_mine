#!/bin/sh -x
rm data/results/*.svg
rm data/figures/*.svg

for NAME in olden truck msnbc dna
do
for FILE in data/results/${NAME}_char.summary
do
    python extra/plot/plot_bar.py $FILE k ghosts
    python extra/plot/plot_bar.py $FILE k distortion
done

for FILE in data/results/${NAME}*char*k_?.summary
do
   python extra/plot/plot_bar.py $FILE tau ghosts
   python extra/plot/plot_bar.py $FILE tau distortion
done

for FILE in data/results/${NAME}*char*k_??.summary
do
   python extra/plot/plot_bar.py $FILE tau ghosts
   python extra/plot/plot_bar.py $FILE tau distortion
done

for FILE in data/results/${NAME}*char*k*tau*.summary
do
   python extra/plot/plot_bar.py $FILE nb_sensitive ghosts
   python extra/plot/plot_bar.py $FILE nb_sensitive distortion
done
done

python extra/plot/plot_time.py data/results/random20M_char_k_5_tau_10.summary nb_sensitive
python extra/plot/plot_time.py data/results/random20M_char_k_5.summary tau
python extra/plot/plot_time.py data/results/random20M_char.summary k
python extra/plot/plot_time.py data/results/random.summary n

python extra/plot/plot_avg_bar.py dna_char.summary k ghosts
python extra/plot/plot_avg_bar.py dna_char.summary k distortion
python extra/plot/plot_avg_bar.py dna_char_k_11.summary tau ghosts
python extra/plot/plot_avg_bar.py dna_char_k_11.summary tau distortion
python extra/plot/plot_avg_bar.py dna_char_k_11_tau_20.summary nb_sensitive ghosts
python extra/plot/plot_avg_bar.py dna_char_k_11_tau_20.summary nb_sensitive distortion

#for i in data/figures/*.svg;do rsvg-convert -f pdf -o ${i%.*}.pdf $i;done
#pdftk data/figures/*.pdf output figures.pdf
rm data/figures/*.pdf
