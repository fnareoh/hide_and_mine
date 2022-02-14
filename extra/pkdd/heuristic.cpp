#include "heuristic.h"
#include <climits>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <math.h>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <algorithm> // std::min

#include <limits.h>
int IMAX = INT_MAX - 1;

using namespace std;

void assign_ghost_and_sub(
    string Zs, unordered_map<string, int> &H, int k, int tau,
    set<char> &distinct_in_Zs, set<string> &sensitive_patterns, int theta,
    string cost_filename, string weight_filename,
    vector<std::pair<int, int>> &knapsack_symbols_fallback) {

  vector<vector<int>> costs_for_all_hashes;
  vector<vector<int>> weights_for_all_hashes;

  int pos = 0;
  int hashmark = 0;
  vector<char> distinct_in_Zs_vec(distinct_in_Zs.begin(), distinct_in_Zs.end());
  for (string::const_iterator it = Zs.begin(); it != Zs.end(); ++it) {
    if (*it == '#') {
      hashmark++;
      vector<int> costs;
      vector<int> weights;
      int replacement_for_hash_int = -1;

      for (int it2 = 0; it2 != distinct_in_Zs_vec.size(); ++it2) {
        string new_part;
        if (distinct_in_Zs_vec[it2] == (char)0)
          new_part =
              Zs.substr(pos - (k - 1), k - 1) + Zs.substr(pos + 1, k - 1);
        else
          new_part = Zs.substr(pos - (k - 1), k - 1) + distinct_in_Zs_vec[it2] +
                     Zs.substr(pos + 1, k - 1);

        int j = 0;

        int ghosts = 0;

        bool creates_sens = false;
        while (j + k <= new_part.length()) {
          string klen = new_part.substr(j, k);
          unordered_map<string, int>::const_iterator h_it = H.find(klen);
          if (h_it != H.end()) {
            set<string>::const_iterator it3 =
                sensitive_patterns.find(h_it->first);

            if (it3 == sensitive_patterns.end()) {
              ghosts++;

            } else {
              creates_sens = true;
            }
          } else {
            set<string>::const_iterator it3 = sensitive_patterns.find(klen);

            if (it3 != sensitive_patterns.end()) {

              creates_sens = true;
            }
          }
          j++;
        }

        if (creates_sens == false) {
          costs.push_back(-1 * ghosts);
          weights.push_back(-1);
          if (replacement_for_hash_int == -1)
            replacement_for_hash_int = it2;
        } else {
          costs.push_back(0);
          weights.push_back(IMAX);
        }
      }
      costs_for_all_hashes.push_back(costs);
      weights_for_all_hashes.push_back(weights);
      knapsack_symbols_fallback.push_back(
          make_pair(hashmark, replacement_for_hash_int));
    }

    pos++;
  }

  ofstream cfile(cost_filename.c_str(), std::ofstream::out);

  for (vector<vector<int>>::const_iterator it = costs_for_all_hashes.begin();
       it != costs_for_all_hashes.end(); ++it) {
    vector<int> cur = *it;

    for (auto &it : cur) {
      cfile << it << " ";
    }
    cfile << "\n";
  }
  cfile.close();

  ofstream wfile(weight_filename.c_str(), std::ofstream::out);

  for (vector<vector<int>>::const_iterator it = weights_for_all_hashes.begin();
       it != weights_for_all_hashes.end(); ++it) {
    vector<int> cur = *it;
    for (auto &it : cur) {
      wfile << it << " ";
    }
    wfile << "\n";
  }
  wfile.close();
}

void update_H_with_new_in_Zprime(string Zs, unordered_map<string, int> &H,
                                 int k, int tau, set<char> &distinct_in_Zs,
                                 unordered_set<string> &deleted_from_H) {

  int pos = 0;

  for (string::const_iterator it = Zs.begin(); it != Zs.end(); ++it) {
    if (*it == '#') {

      unordered_map<string, int> H_for_hash;

      int cnt = 0;

      for (set<char>::iterator it2 = distinct_in_Zs.begin();
           it2 != distinct_in_Zs.end(); ++it2) {

        unordered_map<string, int> H_for_hash_tmp;

        string new_part;

        if (*it2 == (char)0)
          new_part =
              Zs.substr(pos - (k - 1), k - 1) + Zs.substr(pos + 1, k - 1);
        else
          new_part = Zs.substr(pos - (k - 1), k - 1) + *it2 +
                     Zs.substr(pos + 1, k - 1);

        int j = 0;
        while (j + k <= new_part.length()) {
          string klen = new_part.substr(j, k);
          unordered_map<string, int>::const_iterator h_it =
              H_for_hash_tmp.find(klen);

          if (h_it == H_for_hash_tmp.end())
            H_for_hash_tmp[klen] = 1;
          else
            H_for_hash_tmp[klen]++;

          j++;
        }
        if (cnt == 0) {
          H_for_hash = H_for_hash_tmp;
        } else {
          for (auto &it5 : H_for_hash_tmp) {
            unordered_map<string, int>::const_iterator it6 =
                H_for_hash.find(it5.first);

            if (it6 == H_for_hash.end()) {
              H_for_hash[it5.first] = it5.second;
            } else {
              if (it6->second > it5.second) {
                H_for_hash[it6->first] = it5.second;
              }
            }
          }
        }
        cnt++;
      }

      for (auto &it5 : H_for_hash) {
        unordered_map<string, int>::const_iterator it6 = H.find(it5.first);
        if (it6 == H.end())
          H[it5.first] = it5.second;
        else {
          H[it6->first] = it6->second + it5.second;
        }
      }
    }

    pos++;
  }

  for (unordered_map<string, int>::iterator it = H.begin(); it != H.end();) {
    unordered_set<string>::iterator it2 = deleted_from_H.find(it->first);

    if (it->second < tau || it2 != deleted_from_H.end()) {
      it = H.erase(it);
    } else
      ++it;
  }
}

void add_to_H_frequent_in_Z(string Zs, unordered_map<string, int> &H, int k,
                            int tau, unordered_set<string> &deleted_from_H) {

  vector<string> tokens;

  stringstream check1(Zs);

  string token;

  while (getline(check1, token, '#')) {
    tokens.push_back(token);
  }

  for (vector<string>::iterator it = tokens.begin(); it != tokens.end(); ++it) {
    string cur = *it;

    int j = 0;
    while (j + k <= cur.length()) {
      string klen = cur.substr(j, k);
      if (klen.length() < k) {

        continue;
      }

      unordered_map<string, int>::const_iterator h_it = H.find(klen);

      if (h_it == H.end())
        H[klen] = 1;
      else
        H[klen]++;

      j++;
    }
  }

  for (unordered_map<string, int>::iterator it = H.begin(); it != H.end();) {

    if (it->second >= tau) {

      deleted_from_H.insert(it->first);
      it = H.erase(it);
    } else
      ++it;
  }
}
