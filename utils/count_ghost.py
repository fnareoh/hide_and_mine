from collections import defaultdict
import sys
import operator
from pprint import pprint

# Input: k tau input_file forbiden_patterns_file output_file
k = int(sys.argv[1])
tau = int(sys.argv[2])
input_file = sys.argv[3]
pattern_file = sys.argv[4]
# file with the replacement, we measured the ghost with this replacement
output_file = sys.argv[5]

S = ""
forbiden_patterns = []
new_S = ""

with open(input_file, "r") as file:
    S = file.read().replace("\n", "")

with open(output_file, "r") as file:
    new_S = file.read().replace("\n", "")

with open(pattern_file, "r") as file:
    forbiden_patterns = file.read().split("\n")


def pos_sharp(T):
    # function returns the position of the first #
    # if there is one, else returns -1
    for j in range(len(T)):
        if T[j] == "#":
            return j
    return -1


i = 0
kmers_occ = defaultdict(int)
# Compute the ditribution of context in front and behind
while i + k < len(S):
    p = pos_sharp(S[i : i + k])
    if p == -1:
        # Get all the kmer of context
        while i + k - 1 < len(S) and S[i + k - 1] != "#":
            kmers_occ[S[i : i + k]] += 1
            i += 1
        i += k
    else:
        i += p + 1

i = 0
new_kmers_occ = defaultdict(int)
nb_not_replaced = 0
nb_ghosts = 0
is_working = True
# Compute the ditribution of context in front and behind
while i + k < len(new_S):
    p = pos_sharp(new_S[i : i + k])
    if p == -1:
        # Get all the kmer of context
        while i + k - 1 < len(new_S) and new_S[i + k - 1] != "#":
            if new_S[i : i + k] in forbiden_patterns:
                print(new_S[i : i + k])
                is_working = False
            elif (
                kmers_occ[new_S[i : i + k]] < tau
                and new_kmers_occ[new_S[i : i + k]] == tau - 1
            ):
                nb_ghosts += 1
            new_kmers_occ[new_S[i : i + k]] += 1
            i += 1
        nb_not_replaced += 1
        i += k
    else:
        nb_not_replaced += 1
        i += p + 1

if not is_working:
    print("FORBIDDEN PATTERN WAS FOUND!")
    print("THERE IS A BUG!")
# print the tau ghost created
print("Number not replaced: ", nb_not_replaced)
print("Number of ghosts: ", nb_ghosts)
