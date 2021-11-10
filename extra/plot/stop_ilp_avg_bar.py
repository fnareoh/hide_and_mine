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
TIME_STOP = 3550  # Nb seconds after which the ILP is stopped

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
    avg["stopped"] = copy.deepcopy(avg["time"])
    var = copy.deepcopy(list_data[0])
    to_avg = ["distortion", "ghosts", "time"]
    for key, val in list_data[0].items():
        if key == "|S|":
            for i in range(len(val)):
                avg[key][i] = 0
                var[key][i] = 0
                for r in range(N):
                    avg[key][i] += list_data[r][key][i]
                avg[key][i] = avg[key][i] // N
                for r in range(N):
                    var[key][i] += (list_data[r][key][i] - avg[key][i]) ** 2
                var[key][i] = isqrt(var[key][i] // N)
        if key in to_avg:
            for i in range(len(val)):
                for j in range(len(val[i])):
                    avg[key][i][j] = 0
                    var[key][i][j] = 0
                    if key == "time":
                        avg["stopped"][i][j] = 0
                    for r in range(N):
                        if key == "time" and list_data[r][key][i][j] >= TIME_STOP:
                            avg["stopped"][i][j] += 1
                        else:
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
    pprint(avg)
    return avg, var


avg, var = parse_10_results_files(file_name, variation, metric)
data = avg

labels = []
if data["method"][0][0] == "minimize_sum_unfrequent_distance_to_tau":
    labels = ["G-HEU", "FIXED", "RANDOM"]
else:
    labels = ["TPM", "ILP", "HEU"]

with open("data/avg_results/" + file_name + ".ghosts", "w") as file:
    writer = csv.writer(file)
    # writer.writerow(["dataset"] + labels)
    for i in range(len(data["k"])):
        name_dataset = (
            f"{file_name[:3]}_k_{data['k'][i]}_tau_{data['tau'][i]}_S_{data['|S|'][i]}"
        )
        line = [name_dataset]
        for j in range(len(labels)):
            line.append(data["ghosts"][j][i])
        writer.writerow(line)

with open("data/avg_results/" + file_name + ".distortion", "w") as file:
    writer = csv.writer(file)
    # writer.writerow(["dataset"] + labels)
    for i in range(len(data["k"])):
        name_dataset = (
            f"{file_name[:3]}_k_{data['k'][i]}_tau_{data['tau'][i]}_S_{data['|S|'][i]}"
        )
        line = [name_dataset]
        for j in range(len(labels)):
            line.append(data["distortion"][j][i])
        writer.writerow(line)


def reorder(list1, list2):
    zipped_lists = zip(list1, list2)
    sorted_pairs = sorted(zipped_lists)
    tuples = zip(*sorted_pairs)
    f_list1, f_list2 = [list(tuple) for tuple in tuples]
    return f_list1, f_list2


if variation == "k":
    X = [data["k"][i] for i in range(len(data["k"]))]
    x_labels = [
        str(data["k"][i]) + "\n" + str(data["|S|"][i]) for i in range(len(data["k"]))
    ]
elif variation == "tau":
    X = [data["tau"][i] for i in range(len(data["k"]))]
    x_labels = [
        str(data["tau"][i]) + "\n" + str(data["|S|"][i])
        for i in range(len(data["tau"]))
    ]
elif variation == "n":
    X = [int(data["n"][i][:-1]) for i in range(len(data["n"]))]
    x_labels = [
        str(data["n"][i]) + "\n" + str(data["|S|"][i]) for i in range(len(data["n"]))
    ]
else:
    X = [data["nb_sensitive"][i] for i in range(len(data["k"]))]
    x_labels = [
        str(data["nb_sensitive"][i]) + "\n" + str(data["|S|"][i])
        for i in range(len(data["nb_sensitive"]))
    ]

if metric == "ghosts":
    plot = [data["ghosts"][i] for i in range(3)]
    plot_var = [var["ghosts"][i] for i in range(3)]
else:
    plot = [data["distortion"][i] for i in range(3)]
    plot_var = [var["distortion"][i] for i in range(3)]

print(X, x_labels)
sorted_x, sorted_labels = reorder(X, x_labels)
print(X, sorted_labels)
_, plot0 = reorder(X, plot[0])
_, _varplot0 = reorder(X, plot_var[0])
_, plot1 = reorder(X, plot[1])
_, _varplot1 = reorder(X, plot_var[1])
_, stopped1 = reorder(X, avg["stopped"][1])
_, plot2 = reorder(X, plot[2])
_, _varplot2 = reorder(X, plot_var[2])


def remove_specific_S(list_of_list):
    if variation == "nb_sensitive":
        for i in range(len(sorted_x) - 1, -1, -1):
            if sorted_x[i] in [650, 700]:
                for l in list_of_list:
                    l.pop(i)


remove_specific_S(
    [sorted_labels, plot0, _varplot0, plot1, _varplot1, plot2, _varplot2, stopped1]
)

x = np.arange(start=0, stop=3 * len(sorted_labels), step=3)  # the label locations
width = 1.6  # the width of the bars

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

offset = max(plot1 + plot0 + plot2) * 0.08
print(stopped1)
for i in range(len(stopped1)):
    if stopped1[i] > 0:
        plt.plot(x[i], plot1[i] + offset, "kX")
        ax.annotate(
            "{}/{}".format(stopped1[i], N),
            (x[i], plot1[i] + offset),
            xytext=(0, +5),  # 3 points vertical offset
            fontsize=7,
            # backgroundcolor="w",
            textcoords="offset points",
            ha="center",
            va="bottom",
        )


# Add some text for labels, title and custom x-axis tick labels, etc.
if metric == "ghosts":
    ax.set_ylabel("Ghosts")
else:
    ax.set_ylabel("Distortion")

sym_nb_hashes = "|P|"
if data["method"][0][0] == "minimize_sum_unfrequent_distance_to_tau":
    sym_nb_hashes = "δ"
if variation == "k":
    ax.set_xlabel("k\n" + sym_nb_hashes)
elif variation == "tau":
    ax.set_xlabel("τ" + "\n" + sym_nb_hashes, fontsize=19)
elif variation == "n":
    ax.set_xlabel("n\n" + sym_nb_hashes)
else:
    ax.set_xlabel("# sensitive patterns\n" + sym_nb_hashes)

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
plt.legend(loc="upper left", prop={"size": 10})
if metric == "distortion":
    plt.legend(
        loc="upper center", bbox_to_anchor=(0.5, 1.25), ncol=3, prop={"size": 12},
    )
elif labels[0] == "G-HEU":
    plt.legend(
        loc="upper center", bbox_to_anchor=(0.5, 1.18), ncol=3, prop={"size": 12},
    )
plt.tight_layout()
plt.gray()
# plt.show()
fig.savefig(
    "data/figures/"
    + "avg_stop_"
    + file_name.split("/")[-1].split(".")[0]
    + "_"
    + variation
    + "_"
    + metric
    + ".svg",
    dpi=fig.dpi,
    bbox_inches="tight",
)
