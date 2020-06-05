# Hide and mine solver

This repository contains solutions to the Hide and Mine problem. There are separated in two kind, one using an ILP solver ([Gurobi](https://www.gurobi.com/)) and the others are greedy heuristics.

To use them you will need to first install Gurobi (they have free academic licenses available). Both ILP and the heuristics are written in C++ 11.

# Usage

```
./compile.sh
./run.sh 3 20 data/hashmark_input/truck_char_k_3_tau_20_m_30.txt data/sensitive_pattern/truck_char_k_3_tau_20_m_30.txt
```
