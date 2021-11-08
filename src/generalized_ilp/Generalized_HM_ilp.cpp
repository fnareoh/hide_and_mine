#include "gurobi_c++.h"
using namespace std;

#include <chrono>        // chrono
#include <cmath>         // pow
#include <fstream>       // std::ifstream
#include <iostream>      // cin, cout
#include <list>          // list
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
  vector<pair<string, string>> hashmark_context;
  vector<pair<int, int>> P;
  int nb_hashmark;
  int nb_critical;
  map<string, int> criticals_index;
  vector<string> list_criticals;
  vector<tuple<int, int, int>> forbiden_replacements;
};

bool is_critical(string kmer, Input &input) {
  if (input.kmers_frequency.count(kmer) == 0)
    input.kmers_frequency[kmer] = 0;
  return (input.kmers_frequency[kmer] < tau &&
          input.kmers_frequency[kmer] + k * input.nb_hashmark >= tau);
}

// Parse the input and get all informtion needed, occurences of kmer, context
// etc
void parse_input(string input_file, string forbiden_pattern_file,
                 Input &input) {
  unordered_set<char> alphabet;
  ifstream is(input_file); // open file
  char c;
  list<int> hashmark_queue; // contains the id of the hash in the window
  string window = "";       // window with chars only
  string window_hash = "";  // window with hashes
  bool hash_has_back_context = true;
  while (is.get(c)) {
    if (c == '\n') {
      continue;
    }
    if (c == '#') {
      if (!hash_has_back_context) {
        input.hashmark_context.back().second = window;
      }
      int id = input.hashmark_context.size();
      std::list<int>::iterator it;
      for (it = hashmark_queue.begin(); it != hashmark_queue.end(); ++it) {
        input.P.push_back(make_pair(*it, id));
      }
      hashmark_queue.push_back(id);
      input.hashmark_context.push_back(make_pair(window, ""));
      window = "";
      hash_has_back_context = false;
    } else {
      alphabet.insert(c);
      if (window.size() < k - 1)
        window += c;
      else {
        if (!hash_has_back_context) {
          input.hashmark_context.back().second += window;
          hash_has_back_context = true;
        }
        window += c; // add new character
        if (input.kmers_frequency.count(window) == 0)
          input.kmers_frequency[window] = 0;
        input.kmers_frequency[window]++;
        window.erase(0, 1); // remove first character
      }
    }
    window_hash += c;
    if (window_hash.size() >= k) {
      if (window_hash[0] == '#')
        hashmark_queue.pop_front();
      window_hash.erase(0, 1); // remove first character
    }
  }
  is.close(); // close file
  if (!hash_has_back_context) {
    input.hashmark_context.back().second = window;
  }

  cout << "size kmer frequency: " << input.kmers_frequency.size() << endl;
  input.nb_hashmark = input.hashmark_context.size();
  // Parse forbiden_patterns
  ifstream is_forbiden_patterns(forbiden_pattern_file); // open file
  string pattern;
  while (getline(is_forbiden_patterns, pattern))
    input.forbiden_patterns.insert(pattern);
  input.alphabet.assign(alphabet.begin(), alphabet.end());
  input.alphabet.push_back('\0');

  input.nb_critical = 0;
}

/*
void output(vector<vector<char>> replacement, Input &input,
            std::string input_file) {
  int i_hash = 0;
  std::ifstream is(input_file); // open file
  std::ofstream os("data/output/" +
                   input_file.substr(input_file.find_last_of("/\\") + 1) +
                   ".output_ilp");
  char c;
  while (is.get(c)) {
    if (c == '\n') {
      continue;
    }
    if (c == '#') {
      if (i_hash >= input.hashmark.size())
        cout << "i_hash: " << i_hash
             << " size hashmark: " << input.hashmark.size() << std::endl;
      auto context_i = make_pair(input.hashmark[i_hash].substr(0, k - 1),
                                 input.hashmark[i_hash].substr(k, k - 1));
      int index_contex_i = input.context_hashmark_index[context_i];
      char replaced;
      if (replacement[index_contex_i].size() > 0) {
        replaced = replacement[index_contex_i].back();
        replacement[index_contex_i].pop_back();
      } else {
        replaced = '#';
      }
      if (replaced != '\0')
        os << replaced;
      i_hash++;
    } else {
      os << c;
    }
  }
  is.close(); // close file
  os.close(); // close file
}*/

/*
int count_critical(int i_, int j_, string replaced, Input &input,
                   vector<GRBLinExpr> &vect_lhs, GRBVar **x) {
  for (int it = 0; it <= (int)replaced.size() - k; it++) {
    std::string kmer = replaced.substr(it, k);
    if (is_critical(kmer, input)) {
      if (input.forbiden_patterns.count(kmer) == 1) {
        input.forbiden_replacements.push_back(make_pair(i_, j_));
        return 1;
      } else if (input.criticals_index.count(kmer) == 0) {
        // First time we see this critical kmer
        input.criticals_index[kmer] = input.nb_critical;
        input.list_criticals.push_back(kmer);
        vector<int> coeff(input.alphabet.size(), 0);
        GRBLinExpr lhs = 0;
        vect_lhs.push_back(lhs);
        input.nb_critical++;
      }
      // account for the kmer added
      vect_lhs[input.criticals_index[kmer]] += x[i_][j_];
    }
  }
  return 0;
}*/

vector<int> init_replacement(int n) {
  vector<int> J(n, 0);
  return J;
}

void next_replacement(vector<int> &J, int alphabet_size) {
  for (int i = 0; i < J.size(); i++) {
    if (J[i] == alphabet_size - 1) {
      J[i] = 0;
    } else {
      J[i] += 1;
      break;
    }
  }
}

string apply_replacement(int s, int t, vector<int> &J, Input &input) {
  int p_s, p_t;
  string total_context = input.hashmark_context[s].first;
  for (int i = s; i <= t; i++) {
    if (i == s)
      p_s = total_context.size();
    if (i == t)
      p_t = total_context.size();
    total_context +=
        input.alphabet[J[i - s]] + input.hashmark_context[i].second;
  }
  cout << total_context << endl;
  return total_context;
}

int main(int argc, char *argv[]) {
  if (argc != 5)
    cout << "You must give 4 arguments: "
         << "k tau input_file forbiden_patterns" << endl;
  k = stoi(argv[1]);
  tau = stoi(argv[2]);
  string input_file = argv[3];
  string forbiden_pattern_file = argv[4];

  Input input;
  parse_input(input_file, forbiden_pattern_file, input);
  cout << "Finished parsing input !" << endl;
  std::chrono::steady_clock sc;
  auto start = sc.now();

  int nb_context = input.hashmark_context.size();
  vector<vector<char>> replacement(nb_context);
  try {

    // Create an environment
    GRBEnv env = GRBEnv(true);
    env.set("LogFile", "data/HM_ilp.log");
    env.start();

    // Create an empty model
    GRBModel model = GRBModel(env);
    // model.set("TimeLimit", "200.0");
    model.set("Threads", "1");

    // Create variables
    GRBVar **x = new GRBVar *[nb_context];
    for (int i = 0; i < nb_context; i++)
      x[i] = model.addVars(input.alphabet.size(), GRB_INTEGER);

    GRBVar **y = new GRBVar *[input.P.size()];
    for (int i = 0; i < input.P.size(); i++) {
      int nb_to_replace = input.P[i].second - input.P[i].first + 1; // t-s+1
      int nb_replacement =
          (int)pow((float)input.alphabet.size(), nb_to_replace);
      x[i] = model.addVars(input.alphabet.size(), GRB_INTEGER);
    }

    vector<GRBLinExpr> vect_lhs;
    // count added critical kmer
    for (const pair<int, int> &p_hash : input.P) {
      int s = p_hash.first;
      int t = p_hash.second;
      int p_s, p_t;
      string total_context = input.hashmark_context[s].first;
      for (int i = s; i <= t; i++) {
        if (i == s)
          p_s = total_context.size();
        if (i == t)
          p_t = total_context.size();
        total_context += "#" + input.hashmark_context[i].second;
      }
      cout << total_context << endl;
      vector<int> J = init_replacement(t - s + 1);
      int nb_replacement = (int)pow((float)input.alphabet.size(), t - s + 1);
      for (int k = 0; k < nb_replacement; k++) {
        string S = apply_replacement(s, t, J, input);
        next_replacement(J, input.alphabet.size());
      }
    }

    /*
for (map<pair<string, string>, int>::iterator it =
         input.context_hashmark_index.begin();
     it != input.context_hashmark_index.end(); ++it) {
  auto context = it->first;
  string replaced;
  for (int j = 0; j < input.alphabet.size(); j++) {
    if (input.alphabet[j] == '\0')
      replaced = context.first + context.second;
    else
      replaced = context.first + input.alphabet[j] + context.second;
    count_critical(it->second, j, replaced, input, vect_lhs, x);
  }
}
cout << "Finished counting added critical" << endl;

GRBVar *z = model.addVars(input.nb_critical, GRB_BINARY);

// Create objective
GRBLinExpr obj = 0;
for (int l = 0; l < input.nb_critical; l++)
  obj += z[l];
// Set objective: minimize the sum of z
model.setObjective(obj, GRB_MINIMIZE);

// limit added kmer or consider it as a ghost
for (int l = 0; l < input.nb_critical; l++) {
  vect_lhs[l] -= z[l] * k * input.nb_hashmark;
  int e_l = tau - 1 - input.kmers_frequency[input.list_criticals[l]];
  model.addConstr(vect_lhs[l] <= e_l, "limit occurences or count as
ghost");
}

// x[i][j] non negatives
for (int i = 0; i < nb_context; i++) {
  for (int j = 0; j < input.alphabet.size(); j++) {
    model.addConstr(x[i][j] >= 0, "non negative x");
  }
}

// No forbiden replacement
vector<int> count_forbidden_replacement(nb_context, 0);
for (auto &r : input.forbiden_replacements) {
  model.addConstr(x[r.first][r.second] == 0, "no sensitive pattern");
  count_forbidden_replacement[r.first]++;
}

GRBLinExpr lhs;
int nb_impossible_replacement = 0;
// All hashmark must have a replacement
for (int i = 0; i < nb_context; i++) {
  if (count_forbidden_replacement[i] >= input.alphabet.size()) {
    nb_impossible_replacement += input.context_hashmark_count[i];
  } else {
    lhs = 0;
    for (int j = 0; j < input.alphabet.size(); j++) {
      lhs += x[i][j];
    }
    model.addConstr(lhs == input.context_hashmark_count[i],
                    "all hashmark are replaced");
  }
}

std::cout << "Finished building model!" << std::endl;

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
for (int i = 0; i < nb_context; i++) {
  for (int j = 0; j < input.alphabet.size(); j++) {
    for (int occ = 0; occ < x[i][j].get(GRB_DoubleAttr_X); occ++)
      replacement[i].push_back(input.alphabet[j]);
  }
}

int sum_z = 0;
for (int l = 0; l < input.nb_critical; l++) {
  sum_z += z[l].get(GRB_DoubleAttr_X);
}
cout << "Number of context with impossible_replacement: "
     << nb_impossible_replacement << endl;
cout << "Number of ghost: " << sum_z << endl;

output(replacement, input, input_file);

*/
  } catch (GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch (const std::exception &exc) {
    // catch anything thrown within try block that derives from std::exception
    std::cerr << exc.what();
  } catch (...) {
    cout << "Exception during optimization" << endl;
  }

  auto end = sc.now();
  auto time_span = static_cast<chrono::duration<double>>(end - start);
  double runtime = time_span.count();
  std::ofstream os_time("data/output/" +
                        input_file.substr(input_file.find_last_of("/\\") + 1) +
                        ".output_ilp_time");
  os_time << runtime << std::endl;
  os_time.close();
  return 0;
}
