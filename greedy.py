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

with open(text_file, 'r') as file:
    S = file.read().replace('\n', '')

with open(pattern_file, 'r') as file:
    forbiden_patterns = file.read().split('\n')

print(len(S))
i = 0
kmers_occ = defaultdict(int)
context_front = defaultdict(list)
context_back = defaultdict(list)
hashmark = []

def pos_sharp(T):
    # function returns the position of the first #
    # if there is one, else returns -1
    for j in range(len(T)):
        if T[j] == '#':
            return j
    return -1

# Compute the ditribution of context in front and behind
while i + 2*k-2 < len(S):
    p = pos_sharp(S[i:i+2*k-1])
    if p == -1:
        # Get all the kmer of context
        for it in range(k):
            kmers_occ[S[i+it:i+it+k]]+=1
        while (i + 2*k-2 < len(S) and S[i+2*k-2] != '#'):
            context_front[S[i:i+k-1]].append(S[i+k-1])
            context_back[S[i+k:i+2*k-1]].append(S[i+k-1])
            # context[(S[i:i+k-1],S[i+k:i+2*k-1])].append(S[i+k-1]) # context on both side
            kmers_occ[S[i+k-1:i+2*k-1]]+=1
            i += 1
        if (i + 2*k-1 < len(S)):
            hashmark.append(S[i+k-1:i+3*k-2])
        i+=2*k
    else:
        hashmark.append(S[i+p-k+1:i+p+k])
        i=i+p+k

def count_dict(d):
    return [(l, d.count(l) ) for l in set(d)]

def added_kmer(I,kmers_occ):
    # computes what are the unfrequent kmer added
    added = defaultdict(int)
    for it in range(k):
        if I[it:it+k] in forbiden_patterns:
            raise Exception('forbiden pattern')
        if (kmers_occ[I[it:it+k]] < tau):
            added[I[it:it+k]] += 1
    return added

# for all possible replacement (based on the back and front context), compute
# the number of added kmer and choose the letter that minimized the number of
# added kmer
unfrequent = defaultdict(int)
for H in hashmark:
    print((H[0:k-1],H[k:2*k-1]))
    f = count_dict(context_front[H[0:k-1]])
    b = count_dict(context_back[H[k:2*k-1]])
    # pprint(context[(H[0:k-1],H[k:2*k-1])]) # context on both side
    print(f,b)
    possible = sorted(f+b, key=lambda tup: tup[1])
    most_frequent = possible[-1][0]
    # I = H[0:k-1]+most_frequent+H[k:2*k-1]
    letter = '#'
    min_added = k
    for l in possible:
        I = H[0:k-1]+l[0]+H[k:2*k-1]
        try:
            added = added_kmer(I,kmers_occ)
            nb_added = sum(added.values())
            if nb_added <= min_added:
                min_added = nb_added
                letter = l[0]
        except Exception:
            pass
    print(letter, min_added)

    I = H[0:k-1]+letter+H[k:2*k-1]
    for it in range(k):
        if (kmers_occ[I[it:it+k]] < tau):
            # print("Unfrequent kmer: ", I[it:it+k])
            unfrequent[I[it:it+k]] += 1

# print the tau ghost created
for c,v in unfrequent.items():
    if kmers_occ[c]+v >= tau:
        print("Tau-ghost: ", c, " occ: ", kmers_occ[c]+v)
