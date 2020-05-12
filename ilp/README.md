# Hide and Mine Integer Linear Programming

`HM_ilp.py` and `HM_ilp.cpp` are wrappers (respectively in python and C++) for
[Gurobi](https://www.gurobi.com/), for the Hide and Mine problem.

# Usage

for Python:
`python HM_ilp.py 7 5 ../data/input_k_7_tau_5.txt ../data/input_pattern.txt`

Should tell you that there is no solutions because this instance contains some
impossible replacement.

for C++:
```
cmake -H. -Bbuild -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/HM_ilp 7 5 ../data/input_k_7_tau_5.txt ../data/input_pattern.txt
```
