import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import sys
from pprint import pprint
import csv
from matplotlib import rc

plt.rcParams.update({"font.size": 16})

NB_METHOD = 4
TIME_STOP = 3550  # Nb seconds after which the ILP is stopped
file_name = sys.argv[1]
variation = sys.argv[2]
metric = sys.argv[3]
if variation not in ["k", "tau", "n"]:
    variation = "nb_sensitive"
if metric != "ghosts":
    metric = "distortion"
dataset = file_name.split("/")[-1][:3].upper()
content = []
with open(file_name, newline="") as csvfile:
    content = list(csv.reader(csvfile))
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
        for j in range(NB_METHOD):
            if len(data[content[index][i]]) == 0:
                data[content[index][i]] = [[] for _ in range(NB_METHOD)]
            if content[index][i] == "time":
                data[content[index][i]][j].append(float(content[index + j + 1][i]))
            elif content[index][i] == "method":
                data[content[index][i]][j].append(content[index + j + 1][i])
            else:
                data[content[index][i]][j].append(int(content[index + j + 1][i]))
    index += NB_METHOD + 1

pprint(data)
data["stopped"] = [False] * len(data["time"][1])
was_stopped = False
for i in range(len(data["time"][1])):
    if data["time"][1][i] >= TIME_STOP:
        data["stopped"][i] = True
        was_stopped = True


def reorder(list1, list2):
    zipped_lists = zip(list1, list2)
    sorted_pairs = sorted(zipped_lists)
    tuples = zip(*sorted_pairs)
    f_list1, f_list2 = [list(tuple) for tuple in tuples]
    return f_list1, f_list2


def pretty_K_plot(v):
    if v > 1000:
        return f"{round(v/1000,1)}K"
    else:
        return v


X = [data[variation][i] for i in range(len(data[variation]))]
labels = [
    str(data[variation][i]) + "\n" + str(pretty_K_plot(data["|S|"][i]))
    for i in range(len(data[variation]))
]
plot = [data[metric][i] for i in range(NB_METHOD)]
width = 10  # the width of the bars
x = np.arange(
    start=0, stop=(NB_METHOD) * width * len(labels), step=(NB_METHOD + 0.5) * width
)  # the label locations

_, sorted_labels = reorder(X, labels)
for i in range(NB_METHOD):
    _, plot[i] = reorder(X, plot[i])

fig, ax = plt.subplots()
if metric == "distortion":
    fig, ax = plt.subplots(figsize=(4, 5))
    plt.xticks(fontsize=14)
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))

if data["method"][0][0] == "minimize_sum_unfrequent_distance_to_tau":
    names = ["G-HEU", "G-ILP", "FIXED", "RANDOM"]
    rects = []
    for i in range(NB_METHOD):
        rects.append(
            ax.bar(
                x - width * (NB_METHOD / 2) + width / 2 + i * width,
                plot[i],
                width,
                label=names[i],
            )
        )
else:
    names = ["TPM", "ILP", "HEU"]
    rects = []
    for i in range(NB_METHOD):
        rects.append(
            ax.bar(x - width / 2 + i * (width / 2), plot[i], width / 2, label=names[i])
        )

offset = max(plot[1]) * 0.08
if was_stopped:
    _, stopped = reorder(X, data["stopped"])
    for i in range(len(stopped)):
        if stopped[i]:
            plt.plot(x[i] - width / 2, plot[1][i] + offset, "kX")

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel(metric.capitalize())


sym_nb_hashes = "|P|"
if data["method"][0][0] == "minimize_sum_unfrequent_distance_to_tau":
    sym_nb_hashes = "δ"
symbols = {"tau": "τ", "nb_sensitive": "# sensitive patterns"}

if variation in symbols:
    ax.set_xlabel(symbols[variation] + "\n" + sym_nb_hashes)
else:
    ax.set_xlabel(variation + "\n" + sym_nb_hashes)

# ax.set_title(metric.capitalize() + " vs " + variation + " - " + dataset)
ax.set_xticks(x)
ax.set_xticklabels(sorted_labels)
ax.legend()


def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate(
            "{}".format(pretty_K_plot(height)),
            xy=(rect.get_x() + rect.get_width() / 2, height),
            xytext=(0, -1),  # 3 points vertical offset
            fontsize=4,
            # backgroundcolor="w",
            textcoords="offset points",
            ha="center",
            va="bottom",
        )


if metric == "ghosts":
    for rect in rects:
        autolabel(rect)

# plt.legend(loc="upper center", bbox_to_anchor=(0.5, 1.5), ncol=4)
# plt.legend().remove()
plt.legend(prop={"size": 14})
if metric == "distortion":
    plt.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, 1.25),
        ncol=NB_METHOD,
        prop={"size": 12},
    )
if data["method"][0][0] == "minimize_sum_unfrequent_distance_to_tau":
    plt.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, 1.25),
        ncol=NB_METHOD,
        prop={"size": 12},
    )
plt.tight_layout()
plt.gray()
# plt.show()
fig.savefig(
    "data/figures/"
    + file_name.split("/")[-1].split(".")[0]
    + "_"
    + variation
    + "_"
    + metric
    + ".svg",
    dpi=fig.dpi,
    bbox_inches="tight",
)
