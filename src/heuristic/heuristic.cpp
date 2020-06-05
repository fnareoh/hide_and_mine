#include <algorithm>     //std::max
#include <chrono>        // std::chrono
#include <fstream>       // std::ifstream
#include <functional>    //std::function
#include <iostream>      // std::cin, std::cout
#include <map>           // std::map
#include <string>        // std::string
#include <unordered_set> // std::unordered_set
#include <vector>        // std::vector

int k = 0;
int tau = 0;

// Data structure for all the information I need
struct Input {
  std::map<std::string, int> k_frequency;
  std::map<std::string, int> original_k_frequency;
  std::unordered_set<std::string> forbiden_patterns;
  std::unordered_set<char> alphabet;
  std::map<std::string, std::map<char, int>> context;
  std::vector<std::string> hashmark;
};

// Print function for debugging
void print_context(std::map<std::string, std::map<char, int>> &context) {
  std::cout << "{" << std::endl;
  for (std::map<std::string, std::map<char, int>>::iterator it =
           context.begin();
       it != context.end(); ++it) {
    std::cout << "\"" << it->first << "\""
              << " : {";
    for (std::map<char, int>::iterator it2 = it->second.begin();
         it2 != it->second.end(); ++it2) {
      std::cout << " \'" << it2->first << "\': " << it2->second << ",";
    }
    std::cout << "}" << std::endl;
  }
  std::cout << "}" << std::endl;
}

// Print function for debugging
void print_k_frequency(std::map<std::string, int> k_frequency) {
  std::cout << "[" << std::endl;
  for (std::map<std::string, int>::iterator it = k_frequency.begin();
       it != k_frequency.end(); ++it) {
    if (it->second < tau) {
      std::cout << "(\"" << it->first << "\""
                << ", " << it->second << ")," << std::endl;
    }
  }
  std::cout << "]" << std::endl;
}

// Returns true if hte string s has any forbiden_pattern, else no
bool has_forbiden_pattern(Input &input, std::string s) {
  for (int i = 0; i <= s.size() - k; i++) {
    if (input.forbiden_patterns.count(s.substr(i, k)) == 1) {
      return true;
    }
  }
  return false;
}

// For the hashmark given in input "context"+"#", it chooses amng the allowed
// replacement (those that don't create forbiden pattern) the one with the
// highest frequency  in the original string, thus minimizing entropy.
std::string highest_frequency(Input &input, std::string hashmark) {
  std::map<char, int> candidates;
  if (input.context.count(hashmark.substr(0, k - 1)) == 0) {
    for (const char &c : input.alphabet) {
      candidates[c] = 0;
    }
  } else
    candidates = input.context[hashmark.substr(0, k - 1)];
  std::string max_char = "";
  int max = 0;
  for (std::map<char, int>::iterator it = candidates.begin();
       it != candidates.end(); ++it) {
    if (has_forbiden_pattern(input, hashmark.substr(0, k - 1) + it->first +
                                        hashmark.substr(k, k - 1)))
      continue;
    if (it->second >= max) {
      max = it->second;
      max_char = it->first;
    }
  }
  std::string res =
      hashmark.substr(0, k - 1) + max_char + hashmark.substr(k, k - 1);
  if (has_forbiden_pattern(input, res)) {
    // std::cout << "No solution possible!" << std::endl;
    return hashmark;
  }
  return res;
}

// Ouputs -1 if the kmer is not unfrequent, else ouputs it's number of
// occurences
int is_unfrequent(Input &input, std::string kmer) {
  if (input.k_frequency.count(kmer) == 0)
    return 0;
  else if (input.k_frequency[kmer] >= tau)
    return -1;
  else
    return input.k_frequency[kmer];
}

// returns -1 if there is a forbiden_pattern and else the number of unfrequent
// kmer in the string
float naive_unfrequent_measure(Input &input, std::string s) {
  float unfrequent = 0;
  for (int i = 0; i <= s.size() - k; i++) {
    if (input.forbiden_patterns.count(s.substr(i, k)) == 1)
      return -1;
    if (is_unfrequent(input, s.substr(i, k)) >= 0)
      unfrequent += 1;
  }
  return unfrequent;
}

// returns -1 if there is a forbiden_pattern and else the sum of inverse of
// distance to tau for all unfrequent kmer
float sum_unfrequent_distance_to_tau(Input &input, std::string s) {
  float unfrequent = 0;
  for (int i = 0; i <= s.size() - k; i++) {
    if (input.forbiden_patterns.count(s.substr(i, k)) == 1)
      return -1;
    if (is_unfrequent(input, s.substr(i, k)) >= 0)
      unfrequent += (float)1 / (tau - is_unfrequent(input, s.substr(i, k)));
  }
  return unfrequent;
}

// returns -1 if there is a forbiden_pattern and else the max of inverse of
// distance to tau for all unfrequent kmer
float max_unfrequent_distance_to_tau(Input &input, std::string s) {
  float unfrequent = 0;
  for (int i = 0; i <= s.size() - k; i++) {
    if (input.forbiden_patterns.count(s.substr(i, k)) == 1)
      return -1;
    if (is_unfrequent(input, s.substr(i, k)) >= 0)
      unfrequent = std::max(
          unfrequent, (float)1 / (tau - is_unfrequent(input, s.substr(i, k))));
  }
  return unfrequent;
}

// Given a method to compute the quantity/importance of the unfrequent added
// (the unfrequent_measurement function) it chooses the allowed replacement
// that minizes it. It is used to avoid code duplication.
std::string minimize_unfrequent(
    Input &input, std::string hashmark,
    std::function<float(Input &, std::string)> unfrequent_measurement) {
  std::map<char, int> candidates;
  // If the context has never been seen before, look at the entire alphabet

  for (const char &c : input.alphabet) {
    candidates[c] = 0;
  }

  std::string min_char = "";
  float nb_unfrequent = unfrequent_measurement(
      input, hashmark.substr(0, k - 1) + hashmark.substr(k, k - 1));
  int min = ((nb_unfrequent == -1) ? 2 * k - 1 : nb_unfrequent);
  // for each possible replacement, compute the risk of adding tau-ghost
  // (given by unfrequent_measurement)
  for (std::map<char, int>::iterator it = candidates.begin();
       it != candidates.end(); ++it) {
    nb_unfrequent =
        unfrequent_measurement(input, hashmark.substr(0, k - 1) + it->first +
                                          hashmark.substr(k, k - 1));
    if (nb_unfrequent ==
        -1) { // there is a forbiden kmer, not a possible replacement
      continue;
    } else if (nb_unfrequent <= min) {
      min = nb_unfrequent;
      min_char = it->first;
    }
  }
  nb_unfrequent = unfrequent_measurement(
      input, hashmark.substr(0, k - 1) + min_char + hashmark.substr(k, k - 1));
  if (nb_unfrequent == -1) {
    // std::cout << "No solution possible!" << std::endl;
    return hashmark; // There are no possible replacement, we just leave the
                     // hashmark
  }
  return hashmark.substr(0, k - 1) + min_char + hashmark.substr(k, k - 1);
}

// Replacement function for the naive counting of unfrequent
std::string minimize_unfrequent_naive(Input &input, std::string hashmark) {
  return minimize_unfrequent(input, hashmark, naive_unfrequent_measure);
}
// Replacement function for the sum of all distance to tau
std::string minimize_sum_unfrequent_distance_to_tau(Input &input,
                                                    std::string hashmark) {
  return minimize_unfrequent(input, hashmark, sum_unfrequent_distance_to_tau);
}
// Replacement function for the  distance to tau
std::string minimize_max_unfrequent_distance_to_tau(Input &input,
                                                    std::string hashmark) {
  return minimize_unfrequent(input, hashmark, max_unfrequent_distance_to_tau);
}

// update the frequency so far based on the replacement chosen, also report any
// ghost creation
int update_frequency_and_count_ghosts(Input &input, std::string replaced) {
  int ghost = 0;
  for (int i = 0; i <= replaced.size() - k; i++) {
    std::string kmer = replaced.substr(i, k);
    if (input.k_frequency.count(kmer) == 0) {
      input.k_frequency[kmer] = 0;
    }
    if (input.k_frequency[kmer] >= tau)
      continue;
    input.k_frequency[kmer]++;
    if (input.k_frequency[kmer] >= tau) {
      // std::cout << "ghost: " << kmer << std::endl;
      ghost++;
    }
  }
  return ghost;
}

// Given a replacement function, computes all replacement and output it to a
// file
void output(
    Input &input, std::string replacement_name,
    std::function<std::string(Input &input, std::string)> replacement_function,
    std::string input_file) {
  input.k_frequency = input.original_k_frequency;
  int nb_ghosts = 0;
  int nb_impossible_replacement = 0;
  int i = 0;
  std::chrono::steady_clock sc;
  auto start = sc.now();
  std::ifstream is(input_file); // open file
  std::ofstream os("data/output/" +
                   input_file.substr(input_file.find_last_of("/\\") + 1) +
                   ".output_" + replacement_name);
  char c;
  while (is.get(c)) {
    if (c == '\n') {
      continue;
    }
    if (c == '#') {
      std::string replaced = replacement_function(input, input.hashmark[i]);
      // std::cout << nb_ghosts << std::endl;
      // std::cout << replaced << std::endl;
      if (replaced[k - 1] == '#') {
        nb_impossible_replacement++;
        os << "#";
      } else {
        nb_ghosts += update_frequency_and_count_ghosts(input, replaced);
        if (replaced.size() == 2 * k - 1)
          os << replaced[k - 1];
      }
      i++;
    } else {
      os << c;
    }
  }
  is.close(); // close file
  os.close(); // close file

  auto end = sc.now();
  auto time_span = static_cast<std::chrono::duration<double>>(end - start);
  double runtime = time_span.count();
  std::ofstream os_time("data/output/" +
                        input_file.substr(input_file.find_last_of("/\\") + 1) +
                        ".output_" + replacement_name + "_time");
  os_time << runtime << std::endl;
  os_time.close();

  std::cout << "Number of ghosts created by " << replacement_name << ": "
            << nb_ghosts << std::endl;
  std::cout << "Impossible replacement: " << nb_impossible_replacement
            << std::endl;
}

// Parse the input and get all informtion needed, occurences of kmer, context
// etc
void parse_input(std::string input_file, std::string forbiden_pattern_file,
                 Input &input) {
  std::ifstream is(input_file); // open file
  char c;
  std::string window = "";
  bool hash_has_back_context = true;
  while (is.get(c)) {
    if (c == '\n') {
      continue;
    }
    if (c == '#') {
      input.hashmark.push_back(window + "#");
      window = "";
      hash_has_back_context = false;
    } else {
      input.alphabet.insert(c);
      if (window.size() < k - 1)
        window += c;
      else {
        if (!hash_has_back_context) {
          input.hashmark.back() += window;
          hash_has_back_context = true;
        }

        // save context of length k-1
        // create an entry with empty list if none
        if (input.context.count(window) == 0)
          input.context[window] = {};
        if (input.context[window].count(c) == 0)
          input.context[window][c] = 0;
        input.context[window][c]++; // augment the count of letter
        window += c;                // add new character
        if (input.original_k_frequency.count(window) == 0)
          input.original_k_frequency[window] = 0;
        input.original_k_frequency[window]++;
        window.erase(0, 1); // remove first character
      }
    }
  }
  is.close(); // close file
  // print_context(context);
  input.k_frequency = input.original_k_frequency;

  // Parse forbiden_patterns
  std::ifstream is_forbiden_patterns(forbiden_pattern_file); // open file
  std::string pattern;
  while (getline(is_forbiden_patterns, pattern))
    input.forbiden_patterns.insert(pattern);
}

int main(int argc, char **argv) {
  if (argc != 5)
    std::cout << "You must give 4 arguments: "
              << "k tau input_file forbiden_patterns" << std::endl;
  k = std::stoi(argv[1]);
  tau = std::stoi(argv[2]);
  std::string input_file = argv[3];
  std::string forbiden_pattern_file = argv[4];

  Input input;
  parse_input(input_file, forbiden_pattern_file, input);
  std::cout << "Finished parsing input to get stats!" << std::endl;

  // output(input, "highest_frequency", highest_frequency, input_file);

  // output(input, "minimize_unfrequent_naive", minimize_unfrequent_naive,
  //       input_file);

  output(input, "minimize_sum_unfrequent_distance_to_tau",
         minimize_sum_unfrequent_distance_to_tau, input_file);

  output(input, "minimize_max_unfrequent_distance_to_tau",
         minimize_max_unfrequent_distance_to_tau, input_file);
  return 0;
}
