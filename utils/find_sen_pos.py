import sys

input_file = sys.argv[1]
pattern_file = sys.argv[2]

S = ""
forbiden_patterns = []

with open(input_file, "r") as file:
    S = file.read().replace("\n", "")

with open(pattern_file, "r") as file:
    forbiden_patterns = file.read().split("\n")

k = len(forbiden_patterns[0])
forbiden_patterns = set(forbiden_patterns)
occ = []
i = 0

while i + k < len(S):
    if S[i : i + k] in forbiden_patterns:
        occ.append(i)
    i += 1

outF = open("sensitive_pos/" + pattern_file.split("/")[-1], "w")
outF.write(str(len(occ)) + "\n")
for pos in occ:
    outF.write(str(pos) + "\n")
outF.close()
