#include "gurobi_c++.h"
using namespace std;

#include <fstream>       // std::ifstream
#include <iostream>      // cin, cout
#include <map>           //map
#include <string>        // string
#include <unordered_map> //unordered_map
#include <unordered_set> // unordered_set
#include <utility>       // pair
#include <vector>        // vector

int k = 0;
int tau = 0;

// Data structure for all the information needed
struct Input {
  unordered_map<string, int> kmers_frequency;
  unordered_set<string> forbiden_patterns;
  vector<char> alphabet;
  // map instead of unordered map because hashes of pair are not defined
  map<pair<string, string>, int> context_hashmark_index;
  vector<int> context_hashmark_count;
  int nb_hashmark;
  int nb_critical;
  map<string, int> criticals_index;
  vector<string> list_criticals;
  vector<vector<vector<int>>> alpha;
  vector<pair<int, int>> forbiden_replacements;
};

bool is_critical(string kmer, Input &input) {
  return (input.kmers_frequency[kmer] < tau &&
          input.kmers_frequency[kmer] + k * input.nb_hashmark >= tau);
}

void count_critical(int i_, int j_, string replaced, Input &input) {
  for (int it = 0; it <= replaced.size() - k; it++) {
    std::string kmer = replaced.substr(it, k);
    if (is_critical(kmer, input)) {
      if (input.forbiden_patterns.count(kmer) == 1) {
        input.forbiden_replacements.push_back(make_pair(i_, j_));
      } else if (input.criticals_index.count(kmer) == 0) {
        // First time we see this critical kmer
        input.criticals_index[kmer] = input.nb_critical;
        input.list_criticals.push_back(kmer);
        vector<vector<int>> coeff(input.context_hashmark_index.size(),
                                  vector<int>(input.alphabet.size(), 0));
        input.alpha.push_back(coeff);
        input.nb_critical++;
      }
      // account for the kmer added
      input.alpha[input.criticals_index[kmer]][i_][j_]++;
    }
  }
}

// Parse the input and get all informtion needed, occurences of kmer, context
// etc
Input parse_input(string input_file, string forbiden_pattern_file) {
  Input input;
  vector<string> hashmark;
  unordered_set<char> alphabet;
  ifstream is(input_file); // open file
  char c;
  string window = "";
  bool hash_has_back_context = true;
  while (is.get(c)) {
    if (c == '\n') {
      continue;
    }
    if (c == '#') {
      hashmark.push_back(window + "#");
      window = "";
      hash_has_back_context = false;
    } else {
      alphabet.insert(c);
      if (window.size() < k - 1)
        window += c;
      else {
        if (!hash_has_back_context) {
          hashmark.back() += window;
          hash_has_back_context = true;
        }
        if (input.kmers_frequency.count(window) == 0)
          input.kmers_frequency[window] = 0;
        input.kmers_frequency[window]++;
        window += c;        // add new character
        window.erase(0, 1); // remove first character
      }
    }
  }
  is.close(); // close file
  input.nb_hashmark = hashmark.size();
  // Parse forbiden_patterns
  ifstream is_forbiden_patterns(forbiden_pattern_file); // open file
  string pattern;
  while (getline(is_forbiden_patterns, pattern))
    input.forbiden_patterns.insert(pattern);
  input.alphabet.assign(alphabet.begin(), alphabet.end());

  int hashmark_index = 0;
  input.nb_critical = 0;

  // convert hashmark in unique context, count them
  for (string &hash : hashmark) {
    auto context = make_pair(hash.substr(0, k - 1), hash.substr(k, k - 1));
    if (input.context_hashmark_index.count(context) == 0) {
      input.context_hashmark_index[context] = hashmark_index;
      input.context_hashmark_count.push_back(0);
      hashmark_index++;
    }
    input.context_hashmark_count[input.context_hashmark_index[context]]++;
  }

  // Compute the critical kmer and the number we may add through replacement
  // (alpha)
  for (map<pair<string, string>, int>::iterator it =
           input.context_hashmark_index.begin();
       it != input.context_hashmark_index.end(); ++it) {
    auto context = it->first;
    for (int l = 0; l < input.alphabet.size(); l++) {
      string replaced = context.first + input.alphabet[l] + context.second;
      count_critical(it->second, l, replaced, input);
    }
  }

  return input;
}

int main(int argc, char *argv[]) {
  if (argc != 5)
    cout << "You must give 4 arguments: "
         << "k tau input_file forbiden_patterns" << endl;
  k = stoi(argv[1]);
  tau = stoi(argv[2]);
  string input_file = argv[3];
  string forbiden_pattern_file = argv[4];

  Input input = parse_input(input_file, forbiden_pattern_file);
  cout << "Finished parsing input to get stats!" << endl;

  try {

    // Create an environment
    GRBEnv env = GRBEnv(true);
    env.set("LogFile", "HM_ilp.log");
    env.start();

    // Create an empty model
    GRBModel model = GRBModel(env);

    // Create variables
    int nb_context = input.context_hashmark_index.size();
    GRBVar **x = new GRBVar *[nb_context];
    for (int i = 0; i < nb_context; i++)
      x[i] = model.addVars(input.alphabet.size(), GRB_INTEGER);

    GRBVar *z = model.addVars(input.nb_critical, GRB_BINARY);

    // Create objective
    GRBLinExpr obj = 0;
    for (int l = 0; l < input.nb_critical; l++)
      obj += z[l];
    // Set objective: minimize the sum of z
    model.setObjective(obj, GRB_MINIMIZE);

    // x[i][j] non negatives
    for (int i = 0; i < nb_context; i++) {
      for (int j = 0; j < input.alphabet.size(); j++) {
        model.addConstr(x[i][j] >= 0, "non negative x");
      }
    }

    // No forbiden replacement
    for (auto &r : input.forbiden_replacements) {
      model.addConstr(x[r.first][r.second] == 0, "no sensitive pattern");
    }

    GRBLinExpr lhs;
    // limit added kmer or consider it as a ghost
    for (int l = 0; l < input.nb_critical; l++) {
      int e_l = tau - 1 - input.kmers_frequency[input.list_criticals[l]];
      lhs = 0;
      for (int i = 0; i < nb_context; i++) {
        for (int j = 0; j < input.alphabet.size(); j++) {
          lhs += input.alpha[l][i][j] * x[i][j];
        }
      }
      lhs -= z[l] * k * input.nb_hashmark;
      model.addConstr(lhs <= e_l, "limit occurences or count as ghost");
    }

    // All hashmark must have a replacement
    for (int i = 0; i < nb_context; i++) {
      lhs = 0;
      for (int j = 0; j < input.alphabet.size(); j++) {
        lhs += x[i][j];
      }
      model.addConstr(lhs == input.context_hashmark_count[i],
                      "all hashmark are replaced");
    }

    // Optimize model
    model.optimize();
    int status = model.get(GRB_IntAttr_Status);

    if ((status == GRB_INF_OR_UNBD) || (status == GRB_INFEASIBLE) ||
        (status == GRB_UNBOUNDED)) {
      cout << "The model cannot be solved "
           << "because it is infeasible or unbounded" << endl;
      return 1;
    }
    if (status != GRB_OPTIMAL) {
      cout << "Optimization was stopped with status " << status << endl;
      return 1;
    }

    // Print x and z result
    /*for (int i = 0; i < nb_context; i++) {
      for (int j = 0; j < input.alphabet.size(); j++) {
        cout << "x[" << i << "][" << j << "]: " << x[i][j].get(GRB_DoubleAttr_X)
             << endl;
      }
    }*/

    int sum_z = 0;
    for (int l = 0; l < input.nb_critical; l++) {
      // cout << "z[" << l << "]: " << z[l].get(GRB_DoubleAttr_X) << endl;
      sum_z += z[l].get(GRB_DoubleAttr_X);
    }
    cout << "Number of ghost: " << sum_z << endl;

  } catch (GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch (...) {
    cout << "Exception during optimization" << endl;
  }

  return 0;
}
