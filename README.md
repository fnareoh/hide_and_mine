# Usage

```
g++ heuristic.cpp -o heuristic
./heuristic 7 5 input_k_7_tau_5.txt input_pattern.txt
```

It gives the number of ghost created by two techniques :
* `highest_frequency`: chooses the letter that follows with the highest probability (among those who don't create any forbidden pattern)
* `minimize_unfrequent`: chooses the letter that minimize the number of unfrequent k-mer added. As the frequency gets updated, if a tau-ghost is created it is considered frequent and can be added as much as needed.

Impossible replacement is the number of # where no replacement (by a letter or the empty word), is possible.
