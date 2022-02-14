#!/bin/sh -x


for NAME in olden dna msnbc truck
do
for FILE in data/results1/${NAME}_char.summary
do
    F="$(basename -- $FILE)"
    python3 extra/plot/plot_avg_bar.py $F k ghosts
    python3 extra/plot/plot_avg_bar.py $F k distortion
    python3 extra/plot/plot_avg_time.py $F k
done

for FILE in data/results1/${NAME}*char*k_?.summary
do
   F="$(basename -- $FILE)"
   python3 extra/plot/plot_avg_bar.py $F tau ghosts
   python3 extra/plot/plot_avg_bar.py $F tau distortion
    python3 extra/plot/plot_avg_time.py $F tau
done

for FILE in data/results1/${NAME}*char*k_??.summary
do
   F="$(basename -- $FILE)"
   python3 extra/plot/plot_avg_bar.py $F tau ghosts
   python3 extra/plot/plot_avg_bar.py $F tau distortion
    python3 extra/plot/plot_avg_time.py $F tau
done

for FILE in data/results1/${NAME}*char*k*tau*.summary
do
   F="$(basename -- $FILE)"
   python3 extra/plot/plot_avg_bar.py $F nb_sensitive ghosts
   python3 extra/plot/plot_avg_bar.py $F nb_sensitive distortion
   python3 extra/plot/plot_avg_time.py $F nb_sensitive
done
done

syn_name=random5M_char
default_k=6
default_tau=7
python3 extra/plot/plot_avg_time.py ${syn_name}_k_${default_k}_tau_${default_tau}.summary nb_sensitive
python3 extra/plot/plot_avg_bar.py ${syn_name}_k_${default_k}_tau_${default_tau}.summary nb_sensitive ghosts
python3 extra/plot/plot_avg_time.py ${syn_name}_k_${default_k}.summary tau
python3 extra/plot/plot_avg_bar.py ${syn_name}_k_${default_k}.summary tau ghosts
python3 extra/plot/plot_avg_time.py ${syn_name}.summary k
python3 extra/plot/plot_avg_bar.py ${syn_name}.summary k ghosts
python3 extra/plot/plot_avg_time.py random.summary n
python3 extra/plot/plot_avg_bar.py random.summary n ghosts

#for i in data/figures/*.svg;do rsvg-convert -f pdf -o ${i%.*}.pdf $i;done
#pdftk data/figures/*.pdf output figures.pdf
rm data/figures/*.pdf
