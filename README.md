# Hide and mine solver

This repository contains solutions to the Hide and Mine problem. There are separated in two kind, one using an ILP solver ([Gurobi](https://www.gurobi.com/)) and the others are greedy heuristics.

To use them you will need to first install Gurobi (they have free academic licenses available). Both ILP and the heuristics are written in C++ 11.

# Usage

```
./compile.sh
./run.sh 3 20 data/hashmark_input/truck_char_k_3_tau_20_m_30.txt data/sensitive_pattern/truck_char_k_3_tau_20_m_30.txt
```

And the results will be output to a file in `data/output` as well as a comparison on the number of ghost in `data/output`.

Both binary files `./src/heuristic/heuristic` and `./src/ilp/build/HM_ilp` take 4 parameters:
* `k`: The pattern length. All sensitive patterns must have length k.
* `tau`: The threshold that defines tau-ghosts.
* `data/hashmark_input/name.txt`:The file containing the input with hashmarks to replace (in a single line).
* `data/sensitive_pattern/name.txt`: The file containing the sensitive patterns (each sensitive pattern must be in a separate line).

# Expermients

All experiments ran  on  an  Intel  i7-3770  CPU  @  3.40GHz  with  16GB  RAM. Our source code was written in C++, and we used the Gurobi solver  v.  9.0.1  (single-thread  configuration)  to  solve ILP instances.
