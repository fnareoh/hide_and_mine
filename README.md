# Hide and mine solver

This repository contains solutions to the Hide and Mine problem. There are separated in two kind, one using an ILP solver ([Gurobi](https://www.gurobi.com/)) and the others are greedy heuristics.

To use them you will need to first install Gurobi (they have free academic licenses available). Both ILP and the heuristics are written in C++ 11.

# Requirements

For IlP: ([Gurobi](https://www.gurobi.com/)) installed with a valid license (free academic licenses are available).

To reproduce the experiments: igraph needed to install pkdd.
For the evaluation spcript: Python 3.
To plot the experiment graph:  Python 3 and the packages matplotlib and numpy.

# Install 

To compile both the ILP and Heuristic, just launch the script `compile.sh` at the root of the repository.

To reproduce the experiments, we compile a sligtly modified version of this [repository](git@github.com:fnareoh/hide_and_mine.git), (dependency on igraph):
```
cd extra/pkdd/
make
cd ../..
```
If the makefile fails you may need to change the path of the igraph library in the  `makefile`, `main.cpp` and `shorten_Z.cpp` the default path is `/usr/local/include/igraph/igraph.h`.

# Usage

Both binary files `./src/heuristic/heuristic` and `./src/ilp/build/HM_ilp` take 4 parameters:
* `k`: The pattern length. All sensitive patterns must have length k.
* `tau`: The threshold that defines tau-ghosts.
* `data/hashmark_input/name.txt`:The file containing the input with hashmarks to replace (in a single line).
* `data/sensitive_pattern/name.txt`: The file containing the sensitive patterns (each sensitive pattern must be in a separate line).
For simplicity to compute the ILP and HEU result on your file you can use the `run.sh` at the root with the following arguments

```
./run.sh k tau data/hashmark_input/name.txt data/sensitive_pattern/name.txt
```

And the results will be output to a file in `data/output` as well as a comparison on the number of ghost in `data/output`.

Example:
```
./run.sh 3 20 data/hashmark_input/truck_char_k_3_tau_20_m_30.txt data/sensitive_pattern/truck_char_k_3_tau_20_m_30.txt
```

# Experiments

All experiments ran  on  an  Intel  i7-3770  CPU  @  3.40GHz  with  16GB  RAM. Our source code was written in C++, and we used the Gurobi solver  v.  9.0.1  (single-thread  configuration)  to  solve ILP instances.

To reproduce the experiments, once you have installed pkdd you can launch the scripts `generate_all.sh` and `generate_close_hashes.sh` (no arguments needed). `10_generate_all.sh` averages all results 10 times 
