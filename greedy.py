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
hashmark = []

def pos_sharp(T):
    # function returns the position of the first #
    # if there is one, else returns -1
    for j in range(len(T)):
        if T[j] == '#':
            return j
    return -1

# Compute the ditribution of context in front and behind
while i + k < len(S):
    p = pos_sharp(S[i:i+k])
    if p == -1:
        # Get all the kmer of context
        while (i + k-1 < len(S) and S[i+k-1] != '#'):
            context_front[S[i:i+k-1]].append(S[i+k-1])
            kmers_occ[S[i:i+k]]+=1
            i += 1
        if (i + k-1 < len(S)):
            hashmark.append(S[i:i+2*k-1])
        i+=k
    else:
        hashmark.append(S[i+p-k+1:i+p])
        i+=p+1

def count_dict(d):
    return [(l, d.count(l) ) for l in set(d)]

def added_kmer(I,kmers_occ):
    # computes what are the unfrequent kmer added
    added = defaultdict(int)
    for it in range(len(I)-k+1):
        if I[it:it+k] in forbiden_patterns:
            raise Exception('forbiden pattern')
        if (kmers_occ[I[it:it+k]] < tau):
            added[I[it:it+k]] += 1
    return added

# for all possible replacement (based on the front context), compute
# the number of added kmer and choose the letter that minimized the number of
# added kmer
unfrequent = defaultdict(int)
nb_impossible_replacement = 0
# pprint(sorted([(c,v) for c,v in kmers_occ.items() if v < tau], key=lambda tup: str(tup[1])+tup[0]))
for H in hashmark:
    # print((H[0:k-1],H[k:2*k-1]))
    f = count_dict(context_front[H[0:k-1]])
    # pprint(context[(H[0:k-1],H[k:2*k-1])]) # context on both side
    possible = sorted(f+[("",0)], key=lambda tup: tup[1])
    most_frequent = possible[-1][0]
    letter = ''
    min_added = k+1
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
    # print(letter, min_added)

    I = H[0:k-1]+letter+H[k:2*k-1]
    # this is not most frequent that does not create forbidden
    # I = H[0:k-1]+most_frequent+H[k:2*k-1]
    for it in range(len(I)-k+1):
        if (I[it:it+k] in forbiden_patterns):
            I = I[0:k-1] + "#" + I[k:2*k-1]
            break
    if I[k-1]=="#":
        nb_impossible_replacement +=1
        continue
    for it in range(len(I)-k+1):
        if (kmers_occ[I[it:it+k]] < tau):
            # print("Unfrequent kmer: ", I[it:it+k])
            kmers_occ[I[it:it+k]] +=1 # update the count of unfrequent
            unfrequent[I[it:it+k]] += 1

# print the tau ghost created
print("Number impossible replacement: ", nb_impossible_replacement )
nb_ghosts = 0
for c,v in unfrequent.items():
    if kmers_occ[c] >= tau:
        print("Tau-ghost: ", c)
        nb_ghosts +=1
print("Number of ghosts: ", nb_ghosts)
