#!/usr/bin/env python3.7

# Copyright 2020, Gurobi Optimization, LLC

# This example formulates and solves the following simple MIP model:
#  maximize
#        x +   y + 2 z
#  subject to
#        x + 2 y + 3 z <= 4
#        x +   y       >= 1
#        x, y, z binary

import gurobipy as gp
from gurobipy import GRB

from collections import defaultdict
import sys
import operator
from pprint import pprint

# Input: k tau input_file forbiden_patterns_file
k = int(sys.argv[1])
tau = int(sys.argv[2])
text_file = sys.argv[3]
pattern_file = sys.argv[4]
S = ""
forbiden_patterns = []

with open(text_file, "r") as file:
    S = file.read().replace("\n", "")

with open(pattern_file, "r") as file:
    forbiden_patterns = set(file.read().split("\n"))


def pos_sharp(T):
    # function returns the position of the first #
    # if there is one, else returns -1
    for j in range(len(T)):
        if T[j] == "#":
            return j
    return -1


i = 0
kmers_occ = defaultdict(int)
occ_context_hashmark = defaultdict(int)
nb_hashmark = 0
alphabet = set()
# Compute the ditribution of context in front and behind
while i + k < len(S):
    p = pos_sharp(S[i : i + k])
    if p == -1:
        # Get all the kmer of context
        while i + k - 1 < len(S) and S[i + k - 1] != "#":
            kmers_occ[S[i : i + k]] += 1
            alphabet.add(S[i + k - 1])
            i += 1
        if i + k - 1 < len(S):
            occ_context_hashmark[(S[i : i + k - 1], S[i + k : i + 2 * k - 1])] += 1
            nb_hashmark += 1
        i += k
    else:
        occ_context_hashmark[
            (S[i + p - k + 1 : i + p], S[i + p + 1 : i + p + k - 1])
        ] += 1
        nb_hashmark += 1
        i += p + 1
context_hashmark = list(occ_context_hashmark)
alphabet = list(alphabet) + [""]

forbiden_replacement = set()
critical_kmers = {}
nb_critical = 0
alpha = []

for i in range(len(context_hashmark)):
    context = context_hashmark[i]
    for j in range(len(alphabet)):
        letter = alphabet[j]
        replacement = context[0] + letter + context[1]
        for start in range(len(replacement) - k + 1):
            kmer = replacement[start : start + k]
            if kmer in forbiden_patterns:
                # kmer is forbiden
                forbiden_replacement.add((i, j))
            elif kmers_occ[kmer] < tau and kmers_occ[kmer] + k * nb_hashmark >= tau:
                # kmer is critical
                if kmer in critical_kmers:
                    alpha[critical_kmers[kmer]][i, j] += 1
                else:
                    critical_kmers[kmer] = nb_critical
                    nb_critical += 1
                    alpha.append(
                        dict(
                            (val, 0)
                            for val in [
                                (i, j)
                                for i in range(len(context_hashmark))
                                for j in range(len(alphabet))
                            ]
                        )
                    )
                    alpha[critical_kmers[kmer]][i, j] += 1

list_critical = [""] * nb_critical
for kmer, order in critical_kmers.items():
    list_critical[order] = kmer

try:

    # Create a new model
    m = gp.Model("Hide and Mine")

    # Create variables
    x = m.addVars(len(context_hashmark), len(alphabet), vtype="I", name="x")
    z = m.addVars(nb_critical, vtype="B", name="z")

    # Set objective
    m.setObjective(z.sum(), GRB.MINIMIZE)

    for i in range(len(context_hashmark)):
        for j in range(len(alphabet)):
            m.addConstr(x[i, j] >= 0)

    for i, j in forbiden_replacement:
        m.addConstr(x[i, j] == 0, "no sensitive pattern")

    for l in range(len(alpha)):
        e_l = tau - 1 - kmers_occ[list_critical[l]]
        m.addConstr(
            x.prod(alpha[l]) - z[l] * k * nb_hashmark <= e_l,
            "limit occurences or ghost",
        )

    for i in range(len(context_hashmark)):
        print(occ_context_hashmark[context_hashmark[i]], context_hashmark[i])
        m.addConstr(
            x.sum(i, "*") == occ_context_hashmark[context_hashmark[i]],
            "all hashmark are replaced",
        )
    # Add constraint: x + 2 y + 3 z <= 4
    # m.addConstr(x <= 4, "c0 ")

    # Add constraint: x + y >= 1
    # m.addConstr(x + y >= 1, "c1")

    # Optimize model
    m.optimize()

    for v in m.getVars():
        print("%s %g" % (v.varName, v.x))

    print("Obj: %g" % m.objVal)

except gp.GurobiError as e:
    print("Error code " + str(e.errno) + ": " + str(e))

except AttributeError:
    print("Encountered an attribute error")
