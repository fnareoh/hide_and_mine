#ifndef heuristic
#define heuristic

#include <set>
#include <stdio.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <algorithm> // std::min

using namespace std;

void add_to_H_frequent_in_Z(string Zs, unordered_map<string, int> &H, int k,
                            int tau, unordered_set<string> &deleted_from_H);
void update_H_with_new_in_Zprime(string Zs, unordered_map<string, int> &H,
                                 int k, int tau, set<char> &distinct_in_Zs,
                                 unordered_set<string> &deleted_from_H);
void assign_ghost_and_sub(
    string Zs, unordered_map<string, int> &H, int k, int tau,
    set<char> &distinct_in_Zs, set<string> &sensitive_patterns, int theta,
    string cost_filename, string weight_filename,
    vector<std::pair<int, int>> &knapsack_symbols_fallback);
#endif
