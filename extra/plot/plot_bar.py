import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import sys
from pprint import pprint
import csv
from matplotlib import rc

plt.rcParams.update({"font.size": 16})


file_name = sys.argv[1]
variation = sys.argv[2]
metric = sys.argv[3]
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

pprint(data)


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
else:
    X = [data["nb_sensitive"][i] for i in range(len(data["k"]))]
    labels = [
        str(data["nb_sensitive"][i]) + "\n" + str(data["|S|"][i])
        for i in range(len(data["nb_sensitive"]))
    ]

if metric == "ghosts":
    plot = [data["ghosts"][i] for i in range(3)]
else:
    plot = [data["distortion"][i] for i in range(3)]

x = np.arange(start=0, stop=3 * len(labels), step=3)  # the label locations
width = 1.6  # the width of the bars

_, sorted_labels = reorder(X, labels)
_, plot0 = reorder(X, plot[0])
_, plot1 = reorder(X, plot[1])
_, plot2 = reorder(X, plot[2])

fig, ax = plt.subplots()
if metric == "distortion":
    fig, ax = plt.subplots(figsize=(4, 5))
    plt.xticks(fontsize=14)
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))

rects1 = ax.bar(x - width / 2, plot0, width / 2, label="TPM")
rects2 = ax.bar(x, plot1, width / 2, label="ILP")
rects3 = ax.bar(x + width / 2, plot2, width / 2, label="HEU")  # h2
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
ax.legend()


def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate(
            "{}".format(height),
            xy=(rect.get_x() + rect.get_width() / 2, height),
            xytext=(0, -1),  # 3 points vertical offset
            fontsize=12,
            # backgroundcolor="w",
            textcoords="offset points",
            ha="center",
            va="bottom",
        )


if metric == "ghosts":
    autolabel(rects1)
    autolabel(rects2)
    autolabel(rects3)

# plt.legend(loc="upper center", bbox_to_anchor=(0.5, 1.5), ncol=4)
# plt.legend().remove()
plt.legend(prop={"size": 14})
if metric == "distortion":
    plt.legend(
        loc="upper center", bbox_to_anchor=(0.5, 1.25), ncol=3, prop={"size": 12},
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
