#!/bin/sh -x
rm results/*.svg

for FILE in results/{dna,old,tru,msn}*char.summary
do
    python utils/plot_bar.py $FILE k ghosts
    python utils/plot_bar.py $FILE k distortion
done

for FILE in results/{dna,old,tru,msn}*char*k*.summary
do
   python utils/plot_bar.py $FILE tau ghosts
   python utils/plot_bar.py $FILE tau distortion
done

for FILE in results/{dna,old,tru,msn}*char*k*tau*.summary
do
   python utils/plot_bar.py $FILE nb_sensitive ghosts
   python utils/plot_bar.py $FILE nb_sensitive distortion
done

python utils/plot_time.py results/random20M_char_k_5_tau_10.summary nb_sensitive
python utils/plot_time.py results/random20M_k_5.summary tau
python utils/plot_time.py results/random20M.summary k
python utils/plot_time.py results/random.summary n

for i in results/*.svg;do rsvg-convert -f pdf -o ${i%.*}.pdf $i;done
pdftk results/*.pdf output figures.pdf
rm results/*.pdf
