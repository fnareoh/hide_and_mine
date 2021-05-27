# Small script to have some insights on the kmer frequency changes
from collections import defaultdict
import sys
import operator
from pprint import pprint
import csv
from os import path


# Input: k tau input_file forbiden_patterns_file output_file
k = int(sys.argv[1])
tau = int(sys.argv[2])
input_file = sys.argv[3]
pattern_file = sys.argv[4]
# base file name
file_name = sys.argv[5]
method_names = []

if len(sys.argv) == 7:
    method_names = ["minimize_sum_unfrequent_distance_to_tau", "constant", "random"]

else:
    method_names = [
        "pkdd",
        "ilp",
        # "minimize_max_unfrequent_distance_to_tau",
        "minimize_sum_unfrequent_distance_to_tau",
    ]


S = ""
forbiden_patterns = []
parameters = []

with open(input_file, "r") as file:
    S = file.read().replace("\n", "")


with open(pattern_file, "r") as file:
    forbiden_patterns = file.read().split("\n")

with open("data/parameters/" + file_name, "r") as file:
    parameters = file.read().split("\n")
    parameters.pop()
    parameters = list(map(int, parameters))
parameters[-1] -= 1  # HACK: because the first line of sens_pos is the number of occ


def pos_sharp(T):
    # function returns the position of the first #
    # if there is one, else returns -1
    for j in range(len(T)):
        if T[j] == "#":
            return j
    return -1


def evaluate(name, writer):
    print("------", name, "------")
    print("data/output/" + file_name + "_" + name)
    if not path.exists("data/output/" + file_name + ".output_" + name):
        writer.writerow([name, "/", "/", "/", "/"])
        return

    new_S = ""
    print("data/output/" + file_name + ".output_" + name)
    with open("data/output/" + file_name + ".output_" + name, "r") as file:
        new_S = file.read().replace("\n", "")

    time = 0
    with open("data/output/" + file_name + ".output_" + name + "_time", "r") as file:
        time = float(file.read())
    print("time:", time)

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
    distortion_changes = []
    for c, v in new_kmers_occ.items():
        if v >= tau and kmers_occ[c] < tau:
            nb_ghosts += 1
        elif v < tau and kmers_occ[c] >= tau:
            nb_losts += 1
        distortion += (kmers_occ[c] - v) ** 2
        distortion_changes.append((v - kmers_occ[c], c, v, kmers_occ[c]))

    if not is_working:
        print(nb_sensitive, "FORBIDDEN PATTERN WAS/WERE FOUND!")
        print("THERE IS A BUG!")
    # print the tau ghost created
    print("Number not replaced: ", nb_not_replaced)
    print("Number of ghosts: ", nb_ghosts)
    print("Number of losts: ", nb_losts)
    print("Number of ghosts and losts: ", nb_ghosts + nb_losts)
    print("Distortion: ", distortion)
    writer.writerow([name, time, nb_not_replaced, nb_ghosts, distortion])
    distortion_changes.sort()
    return distortion_changes


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

# pprint(kmers_occ)
nb_kmer_per_freq = defaultdict(int)
for k, v in kmers_occ.items():
    nb_kmer_per_freq[v] += 1
possible_tau = [3, 5, 7, 10]
for t in possible_tau:
    nb_frequent = 0
    for k, v in nb_kmer_per_freq.items():
        if k >= t:
            nb_frequent += v
    k = int(sys.argv[1])
    print(
        f"Number of infrequent for k={k}, tau={t}: {4**k-nb_frequent}, number frequent: {nb_frequent}"
    )
# pprint(nb_kmer_per_freq)
exit(0)

with open("data/results/" + file_name + ".comparison", "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["k", "tau", "x", "y", "nb_sensitive", "|S|"])
    writer.writerow([k, tau] + parameters)
    writer.writerow(
        ["method", "time", "impossible_replacement", "ghosts", "distortion"]
    )
    list_kmer_changes = []
    for name in method_names:
        list_kmer_changes.append(evaluate(name, writer))

    nb_to_print = len(list_kmer_changes[0])
    changed_most_dist_tpm = list_kmer_changes[0][-nb_to_print:]
    dict_distortion_change = dict()
    print(f"k={k} tau={tau}")
    for c, km, f_new, f_old in changed_most_dist_tpm:
        dict_distortion_change[km] = [f_old, f_new, 0, 0]
    for i in range(1, 3):
        for c, km, f_new, f_old in list_kmer_changes[i]:
            if km in dict_distortion_change:
                dict_distortion_change[km][i + 1] = f_new
    # for c, v in dict_distortion_change.items():
    #    print(f"kmer: {c} original_freq: {v[0]} ", end="")
    #    for i in range(len(method_names)):
    #        print(f"{method_names[i]}: {v[i+1]} ", end="")
    #    print()

    avg_change_tpm_frequent = 0
    nb_frequent = 0
    avg_change_tpm_unfrequent = 0
    nb_unfrequent = 0
    for c, v in dict_distortion_change.items():
        if v[1] - v[0] >= 0:
            if v[0] >= tau:
                nb_frequent += 1
                avg_change_tpm_frequent = v[1] - v[0]
            else:
                nb_unfrequent += 1
                avg_change_tpm_unfrequent = v[1] - v[0]
    avg_change_tpm_frequent = avg_change_tpm_frequent / nb_frequent
    avg_change_tpm_unfrequent = avg_change_tpm_unfrequent / nb_unfrequent
    print("Average change in frequency over all frequent:", avg_change_tpm_frequent)
    print("Average change in frequency over all unfrequent:", avg_change_tpm_unfrequent)
