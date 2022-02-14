from collections import defaultdict
import sys
import operator
from pprint import pprint
import csv
from os import path


def pos_sharp(T):
    # function returns the position of the first #
    # if there is one, else returns -1
    for j in range(len(T)):
        if T[j] == "#":
            return j
    return -1


# Input: k tau input_file forbiden_patterns_file output_file
k = int(sys.argv[1])
tau = int(sys.argv[2])
input_file = sys.argv[3]
pattern_file = sys.argv[4]
# base file name
file_name = sys.argv[5]
method_names = []

if len(sys.argv) == 7:
    method_names = [
        "minimize_sum_unfrequent_distance_to_tau",
        "ilp",
        "constant",
        "random",
    ]

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

with open(pattern_file, "r") as file:
    forbiden_patterns = file.read().split("\n")

with open("data/parameters/" + file_name, "r") as file:
    parameters = file.read().split("\n")
    parameters.pop()
    parameters = list(map(int, parameters))
parameters[-1] -= 1  # HACK: because the first line of sens_pos is the number of occ


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
                    # print(new_S[i : i + k], S[i : i + k], S[: i + k].count("#"))
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
    writer.writerow([name, time, nb_not_replaced, nb_ghosts, distortion])


with open("data/results/" + file_name + ".comparison", "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["k", "tau", "x", "y", "nb_sensitive", "|S|"])
    writer.writerow([k, tau] + parameters)
    writer.writerow(
        ["method", "time", "impossible_replacement", "ghosts", "distortion"]
    )
    for name in method_names:
        evaluate(name, writer)
