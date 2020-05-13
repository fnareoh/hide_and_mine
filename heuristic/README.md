# Hide and Mine Heuristics

It gives the number of ghost created by two heuristics :
* `highest_frequency`: chooses the letter that follows with the highest probability (among those who don't create any forbidden pattern)
* `minimize_unfrequent_naive`: chooses the letter that minimize the total number of unfrequent k-mer added by the replacement.
* `minimize_sum_unfrequent_distance_to_tau` chooses the letter that minimize the sum of all distance to tau (measured as `1/(tau-occ(kmer))`).
* `minimize_max_unfrequent_distance_to_tau` hooses the letter that minimize the maximum of all distance to tau (measured as `1/(tau-occ(kmer))`).


In all those heuristics, as the frequency gets updated, if a tau-ghost is created it is considered frequent and can be added as much as needed.
"Impossible replacement" is the number of hashmarks where no replacement (by a letter or the empty word), is possible.

# Usage

```
g++ heuristic.cpp -o heuristic
./heuristic 7 5 ../data/input_k_7_tau_5.txt ../data/input_pattern.txt
```
