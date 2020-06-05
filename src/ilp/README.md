# Hide and Mine Integer Linear Programming

`HM_ilp.py` and `HM_ilp.cpp` are wrappers (respectively in python and C++) for
[Gurobi](https://www.gurobi.com/), for the Hide and Mine problem.

# Usage

for C++:
```
cmake -H. -Bbuild -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/HM_ilp 7 5 ../data/input_k_7_tau_5.txt ../data/input_pattern.txt
```
