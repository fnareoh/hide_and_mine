#!/bin/sh -x


syn_name=random5M_char
default_k=6
default_tau=7
python3 extra/plot/stop_ilp_avg_time.py ${syn_name}_k_${default_k}_tau_${default_tau}.summary nb_sensitive
python3 extra/plot/stop_ilp_avg_bar.py ${syn_name}_k_${default_k}_tau_${default_tau}.summary nb_sensitive ghosts
python3 extra/plot/stop_ilp_avg_time.py ${syn_name}_k_${default_k}.summary tau
python3 extra/plot/stop_ilp_avg_bar.py ${syn_name}_k_${default_k}.summary tau ghosts
python3 extra/plot/stop_ilp_avg_time.py ${syn_name}.summary k
python3 extra/plot/stop_ilp_avg_bar.py ${syn_name}.summary k ghosts
python3 extra/plot/stop_ilp_avg_time.py random.summary n
python3 extra/plot/stop_ilp_avg_bar.py random.summary n ghosts

