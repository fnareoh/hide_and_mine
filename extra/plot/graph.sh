#!/bin/sh -x

for NAME in olden dna msnbc truck
do
for FILE in data/results/${NAME}_char.summary
do
    F="$(basename -- $FILE)"
    python3 extra/plot/plot_bar.py $FILE k ghosts
    python3 extra/plot/plot_bar.py $FILE k distortion
    python3 extra/plot/plot_time.py $FILE k
done

for FILE in data/results/${NAME}*char*k_?.summary
do
   if [ -f "$FILE" ]; then
   F="$(basename -- $FILE)"
   python3 extra/plot/plot_bar.py $FILE tau ghosts
   python3 extra/plot/plot_bar.py $FILE tau distortion
   python3 extra/plot/plot_time.py $FILE tau
   fi
done

for FILE in data/results/${NAME}*char*k_??.summary
do
   if [ -f "$FILE" ]; then
   F="$(basename -- $FILE)"
   python3 extra/plot/plot_bar.py $FILE tau ghosts
   python3 extra/plot/plot_bar.py $FILE tau distortion
   python3 extra/plot/plot_time.py $FILE tau
   fi
done

for FILE in data/results/${NAME}*char*k*tau*.summary
do
   F="$(basename -- $FILE)"
   python3 extra/plot/plot_bar.py $FILE nb_sensitive ghosts
   python3 extra/plot/plot_bar.py $FILE nb_sensitive distortion
   python3 extra/plot/plot_time.py $FILE nb_sensitive
done
done

syn_name=random20M_char
default_k=6
default_tau=7
python3 extra/plot/plot_time.py data/results/${syn_name}_k_${default_k}_tau_${default_tau}.summary nb_sensitive
python3 extra/plot/plot_bar.py data/results/${syn_name}_k_${default_k}_tau_${default_tau}.summary nb_sensitive ghosts
python3 extra/plot/plot_time.py data/results/${syn_name}_k_${default_k}.summary tau
python3 extra/plot/plot_bar.py data/results/${syn_name}_k_${default_k}.summary tau ghosts
python3 extra/plot/plot_time.py data/results/${syn_name}.summary k
python3 extra/plot/plot_bar.py data/results/${syn_name}.summary k ghosts
python3 extra/plot/plot_time.py data/results/random.summary n
python3 extra/plot/plot_bar.py data/results/random.summary n ghosts


#for i in data/figures/*.svg;do rsvg-convert -f pdf -o ${i%.*}.pdf $i;done
#pdftk data/figures/*.pdf output figures.pdf
rm data/figures/*.pdf
