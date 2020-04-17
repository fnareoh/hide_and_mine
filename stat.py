from collections import defaultdict
import sys
from pprint import pprint

k = int(sys.argv[1])
text_file = sys.argv[2]
S = ""

with open(text_file, 'r') as file:
    S = file.read().replace('\n', '')

i = 0
d = defaultdict(int)

def pos_sharp(T):
    for j in range(len(T)):
        if S[j] == '#':
            return j
    return -1

while i +2*k-1 < len(S):
    p = pos_sharp(S[i:i+2*k-1])
    if p== -1:
        while (i +2*k-1 < len(S) and S[i+2*k-1] != '#'):
            d[S[i:i+2*k-1]] += 1
            i +=1
        i+=k
    else:
        i=i+p+k

count =defaultdict(int)
for c,v in d.items():
    count[v] +=1
print("values: Number of context that apear key times")
pprint(count)
