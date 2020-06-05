from collections import defaultdict
import sys
import operator
from pprint import pprint
import random

# Input: k tau input_file forbiden_patterns_file output_file
k = int(sys.argv[1])
tau = int(sys.argv[2])
input_file = sys.argv[3]
x = int(sys.argv[4])
y = int(sys.argv[5])
pattern_file_to_output = sys.argv[6]
max_pattern = int(sys.argv[7])

S = ""

with open(input_file, "r") as file:
    S = file.read().replace("\n", "")

i = 0
kmers_occ = defaultdict(int)
# Compute the ditribution of context in front and behind
while i + k < len(S):
    kmers_occ[S[i : i + k]] += 1
    i += 1

potentials_sensitives = []
for c, v in kmers_occ.items():
    if x <= v <= y:
        potentials_sensitives.append(c)
nb_sensitive = min(max_pattern, len(potentials_sensitives))
print("Outputs", nb_sensitive, "patterns.")

f = open(pattern_file_to_output, "w+")
random.shuffle(potentials_sensitives)

for i in range(nb_sensitive):
    f.write(potentials_sensitives[i] + "\n")

f.close()
