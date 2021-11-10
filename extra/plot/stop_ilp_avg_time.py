import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import sys
from pprint import pprint
import csv
import copy
import math

from matplotlib import rc

plt.rcParams.update({"font.size": 16})

file_name = sys.argv[1]
variation = sys.argv[2]
metric = "time"
dataset = file_name.split("/")[-1][:3].upper()
content = []
N = 10  # Nb of runs over which the average is done
TIME_STOP = 3550  # Nb seconds after which the ILP is stopped


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
                    if key == "time" and avg["stopped"][i][j] != N:
                        avg[key][i][j] = round(
                            avg[key][i][j] / (N - avg["stopped"][i][j]), 4
                        )
                    else:
                        avg[key][i][j] = avg[key][i][j] // N
                    for r in range(N):
                        if not (key == "time" and list_data[r][key][i][j] >= TIME_STOP):
                            if key == "time" and i == 1 and j == 3:
                                print(list_data[r][key][i][j], avg[key][i][j])
                            var[key][i][j] += (
                                list_data[r][key][i][j] - avg[key][i][j]
                            ) ** 2
                    if key == "time" and avg["stopped"][i][j] != N:
                        var[key][i][j] = math.sqrt(
                            round(var[key][i][j] / (N - avg["stopped"][i][j]), 4)
                        )
                        if key == "time" and i == 1 and j == 3:
                            print(var[key][i][j])
                    else:
                        var[key][i][j] = isqrt(var[key][i][j] // N)
    pprint(avg)
    print(avg["time"][1])
    print(var["time"][1])
    return avg, var


avg, var = parse_10_results_files(file_name, variation, metric)
data = avg


def reorder(list1, list2):
    zipped_lists = zip(list1, list2)
    sorted_pairs = sorted(zipped_lists)
    tuples = zip(*sorted_pairs)
    f_list1, f_list2 = [list(tuple) for tuple in tuples]
    return f_list1, f_list2


plt.figure(figsize=(10, 10))
fig, ax = plt.subplots(figsize=(4, 5))

labels = []
if variation == "k":
    x = [data["k"][i] for i in range(len(data["k"]))]
    labels = [
        str(data["k"][i]) + "\n" + str(data["|S|"][i]) for i in range(len(data["k"]))
    ]
elif variation == "tau":
    x = [data["tau"][i] for i in range(len(data["k"]))]
    labels = [
        str(data["tau"][i]) + "\n" + str(data["|S|"][i])
        for i in range(len(data["tau"]))
    ]

elif variation == "n":
    x = [int(data["n"][i][:-1]) for i in range(len(data["n"]))]
    labels = [
        str(data["n"][i]) + "\n" + str(data["|S|"][i]) for i in range(len(data["n"]))
    ]
else:
    x = [data["nb_sensitive"][i] for i in range(len(data["k"]))]
    labels = [
        str(data["nb_sensitive"][i]) + "\n" + str(data["|S|"][i])
        for i in range(len(data["nb_sensitive"]))
    ]

plot = [data["time"][i] for i in range(3)]
sorted_x, sorted_labels = reorder(x, labels)
_, plot0 = reorder(x, plot[0])
_, _varplot0 = reorder(x, var["time"][0])
_, plot1 = reorder(x, plot[1])
_, stopped1 = reorder(x, avg["stopped"][1])
_, _varplot1 = reorder(x, var["time"][1])
_, plot2 = reorder(x, plot[2])
_, _varplot2 = reorder(x, var["time"][2])
varplot0 = [
    [
        max(1e-2, plot0[i] - max(1e-2, plot0[i] - _varplot0[i]))
        for i in range(len(_varplot0))
    ],
    _varplot0,
]
varplot1 = [
    [
        max(1e-2, plot1[i] - max(1e-2, plot1[i] - _varplot1[i]))
        for i in range(len(_varplot1))
    ],
    _varplot1,
]
varplot2 = [
    [
        max(1e-2, plot2[i] - max(1e-2, plot2[i] - _varplot2[i]))
        for i in range(len(_varplot2))
    ],
    _varplot2,
]


def remove_specific_S(list_of_list):
    if variation == "nb_sensitive":
        for i in range(len(sorted_x) - 1, -1, -1):
            if sorted_x[i] in [650, 700]:
                for l in list_of_list:
                    l.pop(i)


def plot_discontinuous_ILP(plot1, _varplot1):
    plt.errorbar([], [], fmt="C1s-", label="ILP", yerr=[], linewidth=2.5)
    last_pos_to_plot = 0
    ploting = True
    print("plot1:", plot1)
    for i in range(len(plot1)):
        if plot1[i] == 0 and last_pos_to_plot < i and ploting:
            plt.errorbar(
                sorted_labels[last_pos_to_plot:i],
                plot1[last_pos_to_plot:i],
                fmt="C1s-",
                yerr=_varplot1[last_pos_to_plot:i],
                linewidth=2.5,
            )
            ploting = False
        elif plot1[i] != 0 and not ploting:
            ploting = True
            last_pos_to_plot = i
    if ploting:
        plt.errorbar(
            sorted_labels[last_pos_to_plot:],
            plot1[last_pos_to_plot:],
            fmt="C1s-",
            yerr=_varplot1[last_pos_to_plot:],
            linewidth=2.5,
        )


remove_specific_S(
    [sorted_labels, plot0, _varplot0, plot1, _varplot1, plot2, _varplot2, stopped1]
)
plt.errorbar(
    sorted_labels, plot0, fmt="C0o-", label="TPM", yerr=_varplot0, linewidth=2.5
)
plot_discontinuous_ILP(plot1, _varplot1)
plt.errorbar(
    sorted_labels, plot2, fmt="C2^-", label="HEU", yerr=_varplot2, linewidth=2.5
)

# Lablels
ax.set_ylabel("Runtime (log scale) (s)")
if variation == "k":
    ax.set_xlabel("k\n|P|")
elif variation == "tau":
    ax.set_xlabel("τ" + "\n" + "|P|", fontsize=19)
elif variation == "n":
    ax.set_xlabel("n" + "\n" + "|P|")
else:
    ax.set_xlabel("# sensitive patterns\n|P|")

plot_stopped1 = []
for i in range(len(stopped1)):
    if stopped1[i] > 0:
        plot_stopped1.append(sorted_labels[i])

print(plot_stopped1)
plt.plot(plot_stopped1, [3600] * len(plot_stopped1), "kX")
for i in range(len(stopped1)):
    hoff = 0
    if i == 0:
        hoff = 10
    elif i == len(stopped1) - 1:
        hoff = -10
    if stopped1[i] > 0:
        ax.annotate(
            "{}/{}".format(stopped1[i], N),
            (sorted_labels[i], 3600),
            xytext=(hoff, -20),  # 3 points vertical offset
            fontsize=12,
            # backgroundcolor="w",
            textcoords="offset points",
            ha="center",
            va="bottom",
        )

plt.xticks(fontsize=14)
plt.yscale("log")
# ax.set_title(metric.capitalize() + " vs " + variation + " - " + dataset)
ax.legend()


# plt.legend(loc="upper center", bbox_to_anchor=(0.5, 1.5), ncol=4)
# plt.legend().remove()
plt.legend(prop={"size": 12})
plt.tight_layout()
plt.gray()
# plt.show()
fig.savefig(
    "data/figures/avg_stop_time_"
    + file_name.split("/")[-1].split(".")[0]
    + "_"
    + variation
    + "_"
    + metric
    + ".svg",
    dpi=fig.dpi,
    bbox_inches="tight",
)
