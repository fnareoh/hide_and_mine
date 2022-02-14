//
//  main.cpp
//  shorten
//
//  Created by 陈慧娉 on 30/01/2019.
//  Copyright © 2019 陈慧娉. All rights reserved.
//

#include "/usr/local/include/igraph/igraph.h"
//#include "/usr/local/Cellar/igraph/0.7.1_3/include/igraph/igraph.h"
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream> // std::ifstream
#include <iostream>
#include <math.h>
#include <set>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <time.h>
#include <unordered_map>
#include <vector>

#include "functions.h"
#include "heuristic.h"
#include "shorten_Z.hpp"

using namespace std;

std::string getFileExt(const std::string &s) {

  std::size_t i = s.rfind('.', s.length());
  if (i != std::string::npos) {
    return (s.substr(i + 1, s.length() - i));
  }

  return ("");
}

int main(int argc, const char *argv[]) {

  if (argc != 12) {
    cout << "Correct use: ./main text.txt sensitive.txt k tau cost_filename "
            "weight_filename theta sensitive_patterns_filename "
            "output_knapsack_filename output_file output_hashmark_filename"
         << endl;
    return -1;
  }

  FILE *ifile;                       // string
  FILE *input_sen;                   // sensitive string file
  FILE *sensitive_patterns_filename; // sensitive patterns file
  FILE *output_knapsack_filename;
  int n;
  int k = atoi(argv[3]); // length of sensitive string
  int tau = atoi(argv[4]);

  int *C;
  char *y;
  char c;

  string cost_filename(argv[5]);
  string weight_filename(argv[6]);
  int theta = -1 * atoi(argv[7]);
  string out_knap_filename(argv[9]);
  string output_filename(argv[10]);
  string output_hashmark_filename(argv[11]);

  set<string> sensitive_patterns;

  sensitive_patterns_filename = fopen(argv[8], "r");

  if (sensitive_patterns_filename == NULL) {
    printf("File sens patterns not found\n");
    return 0;
  }

  char *sp_line;
  size_t len = 80;

  sp_line = (char *)malloc(len * sizeof(char));
  while ((getline(&sp_line, &len, sensitive_patterns_filename)) != -1) {

    string sp_line_str(sp_line);
    if (sp_line_str.length() > 1)
      sensitive_patterns.insert(
          sp_line_str.substr(0, sp_line_str.length() - 1));
  }
  fclose(sensitive_patterns_filename);
  if (sp_line)
    free(sp_line);

  // cout << "----------RESULTS---HEURISTIC------------" << endl;
  // cout<<"Sensitive patterns:\n";
  // for(auto& it : sensitive_patterns)
  //    cout<<"["<<it<<"]\n";
  string myx;

  char s[255];
  sprintf(s, "%s", argv[1]);
  ifile = fopen(s, "r");
  if (ifile == 0) {
    printf("File not found\n");
    return 0;
  }

  int i = 0;
  while (1) {
    c = fgetc(ifile);
    if (c != EOF) {

      myx += c;
    } else {
      c = '\0';
      break;
    }
  }
  char *x = new char[myx.length() + 1];
  strcpy(x, myx.c_str());

  n = strlen(x) - 1;

  C = (int *)calloc(n, sizeof(int));
  y = (char *)calloc(ceil(double(n - k + 1) / 2) * k +
                         floor((double)(n - k + 1) / 2),
                     sizeof(char));

  if (C == NULL) {
    cout << "Failed to allocate memory for y\n";
    return -1;
  }
  if (y == NULL) {
    cout << "Failed to allocate memory for y\n";
    return -1;
  }

  char t[255];
  sprintf(t, "%s", argv[2]);

  input_sen = fopen(t, "r");
  if (input_sen == 0) {
    printf("File not found\n");
    return 0;
  }

  int *b = (int *)calloc(n, sizeof(int));

  int num = 0;
  fscanf(input_sen, "%d", &num);

  for (int j = 0; j <= (num - 1); j++) {
    int ret = fscanf(input_sen, "%d", &b[j]);

    if (ret == EOF)
      break;

    C[b[j]] = 1;
  }

  if (fclose(ifile))
    printf("File close error.\n");
  if (fclose(input_sen))
    printf("File close error.\n");

  for (i = n - 1; i >= n - k + 1; i--) {
    C[i] = C[n - k];
  }

  int l = sut(x, C, k, '#', y);

  string y_temp(y);
  string y_4 = "";

  string y_2(y_temp.substr(0, y_temp.length() - 1));
  y_4 += y_2;

  char *y3 = const_cast<char *>(y_2.c_str());

  vector<string> Z;
  // shorten(k - 1, y3, Z);

  free(C);
  free(y);

  std::string Zs = y_4;
  /*for (vector<string>::const_iterator it = Z.begin(); it != Z.end() - 1; ++it)
  { Zs += *it;
  }*/

  // Start of the replacement
  int delta = 0;
  for (string::const_iterator it = Zs.begin(); it != Zs.end(); ++it) {
    if (*it == '#')
      delta++;
  }

  unordered_map<string, int> H;

  unordered_set<string> deleted_from_H;

  add_to_H_frequent_in_Z(Zs, H, k, tau,
                         deleted_from_H); // computes frequency of patterns

  set<char> distinct_in_Zs(Zs.begin(), Zs.end()); // everychar in the string
  // TODO = ATCG if DNA
  if (getFileExt(argv[1]) == "sefa") {
    distinct_in_Zs = {'A', 'T', 'C', 'G'};
  }

  distinct_in_Zs.insert((char)0);

  set<char>::iterator it = distinct_in_Zs.find('#');
  if (it != distinct_in_Zs.end())
    distinct_in_Zs.erase(it);

  chrono::steady_clock sc;
  auto start = sc.now();

  update_H_with_new_in_Zprime(Zs, H, k, tau, distinct_in_Zs, deleted_from_H);

  string Z_prime;

  if (delta > 0) {

    vector<std::pair<int, int>> knapsack_symbols_fallback;
    assign_ghost_and_sub(Zs, H, k, tau, distinct_in_Zs, sensitive_patterns,
                         theta, cost_filename, weight_filename,
                         knapsack_symbols_fallback);

    string tmp = "./extra/pkdd/mcknap_clean " + to_string(delta) + " " +
                 to_string(distinct_in_Zs.size()) + " " + to_string(1) + " " +
                 to_string(theta) + " " + cost_filename + " " +
                 weight_filename + " " + to_string(distinct_in_Zs.size()) +
                 " " + to_string(delta) + " " + out_knap_filename;
    system(tmp.c_str());

    output_knapsack_filename = fopen(argv[9], "r");

    vector<std::pair<int, int>> knapsack_symbols;
    char ch;
    if (output_knapsack_filename == NULL) {
      knapsack_symbols = knapsack_symbols_fallback;
      printf("Knapsack output file not found. This happens when knapsack is "
             "unsolvable or has trivial solution and we need to change theta. "
             "we output a random solution.\n");
    } else if (fscanf(output_knapsack_filename, "%c", &ch) == EOF) {
      knapsack_symbols = knapsack_symbols_fallback;
      printf("File is Empty. Knapsack is trivial or unsolvable, we output a "
             "random solution.\n\n");
      fclose(output_knapsack_filename);
    } else {
      output_knapsack_filename = fopen(argv[9], "r");

      char *oknap_line;
      size_t len2 = 80;
      oknap_line = (char *)malloc(len2 * sizeof(char));
      while ((getline(&oknap_line, &len2, output_knapsack_filename)) != -1) {
        string str_oknap_line(oknap_line);

        if (str_oknap_line.length() > 1) {
          char *hash_id = strtok(oknap_line, " ");
          char *replacement_for_hash = strtok(NULL, " ");

          int hash_id_int = atoi(hash_id);
          int replacement_for_hash_int = atoi(replacement_for_hash);

          knapsack_symbols.push_back(
              make_pair(hash_id_int, replacement_for_hash_int));
        } else {
          cout << "Error in knapsack output file. Each line must have two "
                  "integers\n";
          return 0;
        }
      }
      if (oknap_line)
        free(oknap_line);
      fclose(output_knapsack_filename);
    }
    sort(knapsack_symbols.begin(), knapsack_symbols.end());

    vector<char> distinct_in_Zs_vec(distinct_in_Zs.begin(),
                                    distinct_in_Zs.end());
    int hash_no = 0;
    for (string::const_iterator it = Zs.begin(); it != Zs.end(); ++it) {
      if (*it == '#') {

        pair<int, int> p = knapsack_symbols.at(hash_no);
        if (p.second == -1) {
          // cout << "int char is -1, there is an impossible replacement." <<
          // endl;
          Z_prime += "#";
          // return 0;
        } else
          Z_prime += distinct_in_Zs_vec.at(p.second);
        hash_no++;
      } else
        Z_prime += *it;
    }
  } else {
    Z_prime = Zs;
  }

  std::ofstream os(output_filename);
  os << Z_prime << endl;
  os.close();

  auto end = sc.now();
  auto time_span = static_cast<chrono::duration<double>>(end - start);
  double runtime = time_span.count();
  std::ofstream os_time(output_filename + "_time");
  os_time << runtime << endl;
  os_time.close();
  /*
  cout<<"W : "<<string(x,n)<<endl;
  cout<<"X : "<<y_4<<endl;
  cout<<"Y : "<<Zs<<endl;
  cout<<"Z : "<<Z_prime<<endl;
  */

  std::ofstream os_hash(output_hashmark_filename);
  os_hash << Zs << endl;
  os_hash.close();

  exit(0);
}
