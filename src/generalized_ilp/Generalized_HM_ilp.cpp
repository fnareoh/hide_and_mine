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
#include <algorithm>     // min max
#include <cassert>       // assert

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
  vector<tuple<int, int>> forbiden_replacements;
  vector<bool> is_impossible;
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
      hashmark_queue.push_back(id);
      for (it = hashmark_queue.begin(); it != hashmark_queue.end(); ++it) {
        input.P.push_back(make_pair(*it, id));
      }
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
  input.is_impossible = vector<bool> (input.nb_hashmark,false);
  // Parse forbiden_patterns
  ifstream is_forbiden_patterns(forbiden_pattern_file); // open file
  string pattern;
  while (getline(is_forbiden_patterns, pattern))
    input.forbiden_patterns.insert(pattern);
  input.alphabet.assign(alphabet.begin(), alphabet.end());
  // No empty string replacement for generilized HM
  //input.alphabet.push_back('\0');

  input.nb_critical = 0;
}


void output(vector<char> &replacement, Input &input,
            std::string input_file) {
  std::ifstream is(input_file); // open file
  std::ofstream os("data/output/" +
                   input_file.substr(input_file.find_last_of("/\\") + 1) +
                   ".output_ilp");
  char c;
  int i_hash=0;
  while (is.get(c)) {
    if (c == '\n') {
      continue;
    }
    if (c == '#') {
      assert(i_hash <= replacement.size());
      os << replacement[i_hash];
      i_hash++;
    } else {
      os << c;
    }
  }
  is.close(); // close file
  os.close(); // close file
}

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

string apply_replacement(int s, int t, int r, vector<int> &J, Input &input){
  //compute replacement
  string total_context = input.hashmark_context[s].first;
  for (int i = s; i <= t; i++) {
    total_context +=
        input.alphabet[J[i - s]] + input.hashmark_context[i].second;
  }
  //cout << total_context << endl;

  return total_context;
}

int count_critical(int p, int r, string replacement, int p_s, int p_t, Input &input, vector<int> &J,
  vector<GRBLinExpr> &vect_lhs, GRBVar **y){
      if (replacement.size() < k) return 0;
      //count critical due to replacement
      for (int i = max(0,p_t-k+1); i <= min((int)replacement.size()-k,p_s); i++){
        string kmer = replacement.substr(i, k);
        if (input.forbiden_patterns.count(kmer) == 1) {
          input.forbiden_replacements.push_back(make_tuple(p, r));
          return 1;
        }
        if (is_critical(kmer, input)) {
            if (input.criticals_index.count(kmer) == 0) {
            // First time we see this critical kmer
            input.criticals_index[kmer] = input.nb_critical;
            input.list_criticals.push_back(kmer);
            vector<int> coeff(input.alphabet.size(), 0);
            GRBLinExpr lhs = 0;
            vect_lhs.push_back(lhs);
            input.nb_critical++;
          }
          // account for the kmer added
          vect_lhs[input.criticals_index[kmer]] += y[p][r];
      }
      }
    return 0;
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
  try {

    // Create an environment
    GRBEnv env = GRBEnv(true);
    env.set("LogFile", "data/HM_ilp.log");
    env.start();

    // Create an empty model
    GRBModel model = GRBModel(env);
    model.set("TimeLimit", "3600.0");
    model.set("Threads", "1");

    // Create variables
    GRBVar **x = new GRBVar *[nb_context];
    GRBVar **y = new GRBVar *[input.P.size()];

    vector<GRBLinExpr> vect_sum_x(nb_context,0);
    for (int i = 0; i < nb_context; i++){
      x[i] = model.addVars(input.alphabet.size(), GRB_INTEGER);
      // x[i][j] non negatives
      vect_sum_x[i] = 0;
      for (int j = 0; j < input.alphabet.size(); j++) {
        vect_sum_x[i]+=x[i][j];
        model.addConstr(x[i][j] >= 0, "non negative x");
        model.addConstr(x[i][j] <= 1, "x less 1");
      }
    }


    vector<vector<GRBLinExpr>> vect_sum_x_y;
    for (int p = 0; p < input.P.size(); p++) {
      int nb_to_replace = input.P[p].second - input.P[p].first + 1; // t-s+1
      int nb_replacement =
          (int)pow((float)input.alphabet.size(), nb_to_replace);
      y[p] = model.addVars(nb_replacement, GRB_INTEGER);
      // y[i][j] non negatives
      int s = input.P[p].first;
      int t = input.P[p].second;
      vector<int> J = init_replacement(t - s + 1);
      vect_sum_x_y.push_back(vector<GRBLinExpr>(nb_replacement,0));
      for (int r = 0; r < nb_replacement; r++) {
        model.addConstr(y[p][r] >= 0, "non negative y");
        model.addConstr(y[p][r] <= 1, "y less 1");
        for (int i = 0; i < J.size(); i++) {
          vect_sum_x_y[p][r]+=x[s+i][J[i]];
        }
        model.addConstr(t-s+y[p][r] >= vect_sum_x_y[p][r],"replacement J made of single replacement");
        next_replacement(J, input.alphabet.size());
      }
    }


    vector<GRBLinExpr> vect_lhs;
    // count added critical kmer
    int nb_impossible_replacement = 0;
    for (int p = 0; p < input.P.size(); p++) {
      int s = input.P[p].first;
      int t = input.P[p].second;
      int p_s, p_t;
      string total_context = input.hashmark_context[s].first;
      for (int i = s; i <= t; i++) {
        if (i == s)
          p_s = total_context.size();
        if (i == t)
          p_t = total_context.size();
        total_context += "#" + input.hashmark_context[i].second;
      }
      //cout << total_context << endl;
      int impossible_replacement = 0;
      vector<int> J = init_replacement(t - s + 1);
      int nb_replacement = (int)pow((float)input.alphabet.size(), t - s + 1);
      for (int r = 0; r < nb_replacement; r++) {
        string replaced = apply_replacement(s, t, r, J, input);
        int res = count_critical(p,r,replaced, p_s, p_t,input,J,vect_lhs,y);
        if (res==1) {
          model.addConstr(y[p][r] == 0, "no sensitive pattern");
          impossible_replacement+=1;
        }
        next_replacement(J, input.alphabet.size());
      }
      if (impossible_replacement >= nb_replacement) {
        nb_impossible_replacement+=t-s+1;
        for (int i = s; i <= t; i++) {
          input.is_impossible[i]=true;
        }
      }
    }
    cout << "Finished counting added critical" << endl;

    for (int i = 0; i < nb_context; i++){
      if (!input.is_impossible[i]) {
        model.addConstr(vect_sum_x[i] == 1, "replacing # by only one letter");
      }
      else {
        //cout << "Infeasible: " << i << endl;
        model.addConstr(vect_sum_x[i] == 0, "Not replaced");
      }
    }

    GRBVar *z = model.addVars(input.nb_critical, GRB_BINARY);
    // Create objective
    GRBLinExpr obj = 0;
    for (int l = 0; l < input.nb_critical; l++){
      obj += z[l];
      int e_l = tau - 1 - input.kmers_frequency[input.list_criticals[l]];
      model.addConstr(z[l] >= 0, "Non negative z");
      model.addConstr(vect_lhs[l]-k*nb_context*z[l] <= e_l, "limit occurences or count as ghost");
    }
    // Set objective: minimize the sum of z
    model.setObjective(obj, GRB_MINIMIZE);

    std::cout << "Finished building model!" << std::endl;

    // Optimize model
    model.optimize();
    int status = model.get(GRB_IntAttr_Status);

    if (status == GRB_INFEASIBLE) {
      cout << "model is infeasible" << endl;
      model.computeIIS();
     GRBConstr* c = model.getConstrs();
     cout << "Trying to remove contraint(s): " << endl;
     for (int i = 0; i < model.get(GRB_IntAttr_NumConstrs); ++i)
     {
       if (c[i].get(GRB_IntAttr_IISConstr) == 1 and c[i].get(GRB_StringAttr_ConstrName)=="replacing # by only one letter")
       {
         cout << c[i].get(GRB_StringAttr_ConstrName) << endl;
         model.remove(c[i]);
       }
     }
     model.optimize();
     status = model.get(GRB_IntAttr_Status);
     if (status == GRB_INFEASIBLE) {
       cout << "Error: Model is still infeasible !" << endl;
       return 1;
     }
   }
    if ((status == GRB_INF_OR_UNBD) || (status == GRB_UNBOUNDED)) {
      cout << "The model cannot be solved "
      << "because it is infeasible or unbounded" << endl;

      return 1;
    }
    if (status != GRB_OPTIMAL) {
      cout << "Optimization was stopped with status " << status << endl;
      if (model.get(GRB_IntAttr_SolCount) == 0) return 1;
    }

    // Gather solution
    nb_impossible_replacement=0;
    vector<char> replacement(nb_context,'#');
    for (int i = 0; i < nb_context; i++) {
        for (int j = 0; j < input.alphabet.size(); j++) {
          if (x[i][j].get(GRB_DoubleAttr_X)==1){
            replacement[i]=input.alphabet[j];
          }
        }
        if (replacement[i]=='#') nb_impossible_replacement++;
    }

    int sum_z = 0;
    for (int l = 0; l < input.nb_critical; l++) {
      sum_z += z[l].get(GRB_DoubleAttr_X);
    }
    cout << "Number of context with impossible_replacement: "
    << nb_impossible_replacement << endl;
    cout << "Number of ghost: " << sum_z << endl;

    output(replacement, input, input_file);

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
