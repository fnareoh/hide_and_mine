import matplotlib
import math
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import sys
from pprint import pprint
import csv
import copy
from matplotlib import rc

plt.rcParams.update({"font.size": 16})
N = 10  # Number of time the experiments have been repeated

file_name = sys.argv[1]
variation = sys.argv[2]
metric = sys.argv[3]


def isqrt(n):
    x = n
    y = (x + 1) // 2
    while y < x:
        x = y
        y = (x + n // x) // 2
    return x


def parse_content(content):
    index = 0
    data = defaultdict(list)
    while index < len(content):
        for i in range(len(content[index])):
            if content[index][i] == "":
                break
            elif content[index][i] == "n":
                data[content[index][i]].append(content[index + 1][i])
            else:
                data[content[index][i]].append(int(content[index + 1][i]))
        index += 2
        for i in range(len(content[index])):
            if content[index][i] == "":
                break
            for j in range(3):
                if len(data[content[index][i]]) == 0:
                    data[content[index][i]] = [[] for _ in range(3)]
                if content[index][i] == "time":
                    data[content[index][i]][j].append(float(content[index + j + 1][i]))
                elif content[index][i] == "method":
                    data[content[index][i]][j].append(content[index + j + 1][i])
                else:
                    data[content[index][i]][j].append(int(content[index + j + 1][i]))
        index += 4
    # pprint(data)
    return data


def parse_10_results_files(file_name, variation, metric):
    dataset = file_name.split("/")[-1][:3].upper()
    list_content = [[] for i in range(N)]
    list_data = [[] for i in range(N)]
    for i in range(N):
        with open(f"data/results{i+1}/" + file_name, newline="") as csvfile:
            list_content[i] = list(csv.reader(csvfile))
        list_data[i] = parse_content(list_content[i])
    avg = copy.deepcopy(list_data[0])
    var = copy.deepcopy(list_data[0])
    to_avg = ["distortion", "ghosts", "time"]
    for key, val in list_data[0].items():
        if key in to_avg:
            for i in range(len(val)):
                for j in range(len(val[i])):
                    avg[key][i][j] = 0
                    var[key][i][j] = 0
                    for r in range(N):
                        avg[key][i][j] += list_data[r][key][i][j]
                    if key == "time":
                        avg[key][i][j] = round(avg[key][i][j] / N, 4)
                    else:
                        avg[key][i][j] = avg[key][i][j] // N
                    for r in range(N):
                        var[key][i][j] += (
                            list_data[r][key][i][j] - avg[key][i][j]
                        ) ** 2
                    if key == "time":
                        var[key][i][j] = math.sqrt(round(var[key][i][j] / N, 4))
                    else:
                        var[key][i][j] = isqrt(var[key][i][j] // N)
    # pprint(avg)
    return avg, var


avg, var = parse_10_results_files(file_name, variation, metric)
data = avg


def reorder(list1, list2):
    zipped_lists = zip(list1, list2)
    sorted_pairs = sorted(zipped_lists)
    tuples = zip(*sorted_pairs)
    f_list1, f_list2 = [list(tuple) for tuple in tuples]
    return f_list1, f_list2


if variation == "k":
    X = [data["k"][i] for i in range(len(data["k"]))]
    labels = [
        str(data["k"][i]) + "\n" + str(data["|S|"][i]) for i in range(len(data["k"]))
    ]
elif variation == "tau":
    X = [data["tau"][i] for i in range(len(data["k"]))]
    labels = [
        str(data["tau"][i]) + "\n" + str(data["|S|"][i])
        for i in range(len(data["tau"]))
    ]
elif variation == "n":
    X = [data["n"][i] for i in range(len(data["n"]))]
    labels = [
        str(data["n"][i]) + "\n" + str(data["|S|"][i]) for i in range(len(data["n"]))
    ]
else:
    X = [data["nb_sensitive"][i] for i in range(len(data["k"]))]
    labels = [
        str(data["nb_sensitive"][i]) + "\n" + str(data["|S|"][i])
        for i in range(len(data["nb_sensitive"]))
    ]

if metric == "ghosts":
    plot = [data["ghosts"][i] for i in range(3)]
    plot_var = [var["ghosts"][i] for i in range(3)]
else:
    plot = [data["distortion"][i] for i in range(3)]
    plot_var = [var["distortion"][i] for i in range(3)]

x = np.arange(start=0, stop=3 * len(labels), step=3)  # the label locations
width = 1.6  # the width of the bars

_, sorted_labels = reorder(X, labels)
_, plot0 = reorder(X, plot[0])
_, _varplot0 = reorder(X, plot_var[0])
_, plot1 = reorder(X, plot[1])
_, _varplot1 = reorder(X, plot_var[1])
_, plot2 = reorder(X, plot[2])
_, _varplot2 = reorder(X, plot_var[2])
varplot0 = [
    [plot0[i] - max(1e-2, plot0[i] - _varplot0[i]) for i in range(len(_varplot0))],
    _varplot0,
]
varplot1 = [
    [plot1[i] - max(1e-2, plot1[i] - _varplot1[i]) for i in range(len(_varplot1))],
    _varplot1,
]
varplot2 = [
    [plot2[i] - max(1e-2, plot2[i] - _varplot2[i]) for i in range(len(_varplot2))],
    _varplot2,
]

fig, ax = plt.subplots()
if metric == "distortion":
    fig, ax = plt.subplots(figsize=(4, 5))
    plt.xticks(fontsize=14)
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))

labels = []
if data["method"][0][0] == "minimize_sum_unfrequent_distance_to_tau":
    labels = ["HEU", "CONST", "RAND"]
else:
    labels = ["TPM", "HEU", "ILP"]

rects1 = ax.bar(
    x - width / 2,
    plot0,
    width / 2,
    label=labels[0],
    yerr=varplot0,
    ecolor="dimgray",
    capsize=4,
)
rects2 = ax.bar(
    x, plot1, width / 2, label=labels[1], yerr=varplot1, ecolor="dimgray", capsize=4
)
rects3 = ax.bar(
    x + width / 2,
    plot2,
    width / 2,
    label=labels[2],
    yerr=varplot2,
    ecolor="dimgray",
    capsize=4,
)  # h2
# rects3 = ax.bar(x + width / 2, plot[2], width / 2, label="h1") #h1
# rects4 = ax.bar(x + width, plot[3], width / 2, label="h2") #h2

# Add some text for labels, title and custom x-axis tick labels, etc.
if metric == "ghosts":
    ax.set_ylabel("Ghosts")
else:
    ax.set_ylabel("Distortion")


if variation == "k":
    ax.set_xlabel("k\n|P|")
elif variation == "tau":
    ax.set_xlabel("τ" + "\n" + "|P|", fontsize=19)
else:
    ax.set_xlabel("# sensitive patterns\n|P|")

# ax.set_title(metric.capitalize() + " vs " + variation + " - " + dataset)
ax.set_xticks(x)
ax.set_xticklabels(sorted_labels)

if metric == "ghosts":
    ax.bar_label(
        rects1,
        labels=[f"{round(v/1000,1)}K" if v > 1000 else v for v in plot0],
        fontsize=7,
    )
    ax.bar_label(
        rects2,
        labels=[f"{round(v/1000,1)}K" if v > 1000 else v for v in plot1],
        fontsize=7,
    )
    ax.bar_label(
        rects3,
        labels=[f"{round(v/1000,1)}K" if v > 1000 else v for v in plot2],
        fontsize=7,
    )

# plt.legend(loc="upper center", bbox_to_anchor=(0.5, 1.5), ncol=4)
# plt.legend().remove()
plt.legend(loc="best", prop={"size": 10})
if metric == "distortion":
    plt.legend(
        loc="upper center", bbox_to_anchor=(0.5, 1.25), ncol=3, prop={"size": 12},
    )
elif labels[0] == "HEU":
    plt.legend(
        loc="upper center", bbox_to_anchor=(0.5, 1.18), ncol=3, prop={"size": 12},
    )
plt.tight_layout()
plt.gray()
# plt.show()
fig.savefig(
    "data/figures/"
    + "avg_"
    + file_name.split("/")[-1].split(".")[0]
    + "_"
    + variation
    + "_"
    + metric
    + ".svg",
    dpi=fig.dpi,
    bbox_inches="tight",
)
