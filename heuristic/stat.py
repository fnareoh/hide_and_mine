from collections import defaultdict
import sys
from pprint import pprint

# Input k and a name of input file with hashmarks
# Output the distribution of strings of length k that do not contain hashmarks
k = int(sys.argv[1])
text_file = sys.argv[2]
S = ""

with open(text_file, 'r') as file:
    S = file.read().replace('\n', '')

i = 0
d = defaultdict(int)

def pos_sharp(T):
    # function returns the position of the first #
    # if there is one, else returns -1
    for j in range(len(T)):
        if T[j] == '#':
            return j
    return -1

# Compute the ditribution
while i + k < len(S):
    p = pos_sharp(S[i:i+k])
    if p == -1:
        while (i + k < len(S) and S[i+k-1] != '#'):
            d[S[i:i+k]] += 1
            i += 1
        i+=k
    else:
        i=i+p+1

# Output for c a number, output count[c] the number of kmers that apear time c times
count =defaultdict(int)
sensitive = [0]*6
for c,v in d.items():
    count[v] +=1
    if v >= 5 :
        sensitive[5]+=1
    else:
        sensitive[v]+=1
pprint(count)
pprint(sensitive) #number of kmer that appear at least 5 times and exactlu 1,2,3,4 times
