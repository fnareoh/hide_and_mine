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
is_working = True
nb_sensitive = 0
# Compute the ditribution of context in front and behind
while i + k < len(new_S):
    p = pos_sharp(new_S[i : i + k])
    if p == -1:
        # Get all the kmer of context
        while i + k - 1 < len(new_S) and new_S[i + k - 1] != "#":
            if new_S[i : i + k] in forbiden_patterns:
                # print(new_S[i : i + k])
                is_working = False
                nb_sensitive += 1
            else:
                new_kmers_occ[new_S[i : i + k]] += 1
            i += 1
        if i + k - 1 < len(new_S):
            nb_not_replaced += 1
        i += k
    else:
        nb_not_replaced += 1
        i += p + 1

for c, v in kmers_occ.items():
    if c not in forbiden_patterns:
        new_kmers_occ[c] += 0

nb_ghosts = 0
nb_losts = 0
distortion = 0
for c, v in new_kmers_occ.items():
    if v >= tau and kmers_occ[c] < tau:
        nb_ghosts += 1
    elif v < tau and kmers_occ[c] >= tau:
        nb_losts += 1
    distortion += (kmers_occ[c] - v) ** 2

if not is_working:
    print(nb_sensitive, "FORBIDDEN PATTERN WAS/WERE FOUND!")
    print("THERE IS A BUG!")
# print the tau ghost created
print("Number not replaced: ", nb_not_replaced)
print("Number of ghosts: ", nb_ghosts)
print("Number of losts: ", nb_losts)
print("Number of ghosts and losts: ", nb_ghosts + nb_losts)
print("Distortion: ", distortion)
