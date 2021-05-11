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
metric = "time"
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

        if content[index][i] == "n":
            data[content[index][i]].append(content[index + 1][i])
        else:
            data[content[index][i]].append(int(content[index + 1][i]))
    index += 2
    for i in range(len(content[index])):
        if content[index][i] == "":
            break
        for j in range(3):
            if len(data[content[index][i]]) == 0:
                data[content[index][i]] = [[] for _ in range(4)]
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
_, sorted_labels = reorder(x, labels)
_, plot0 = reorder(x, plot[0])
_, plot1 = reorder(x, plot[1])
_, plot2 = reorder(x, plot[2])
print(sorted_labels)
plt.plot(sorted_labels, plot0, "o-", label="TPM", linewidth=2.5)
plt.plot(sorted_labels, plot1, "s-", label="ILP", linewidth=2.5)
plt.plot(sorted_labels, plot2, "^-", label="HEU", linewidth=2.5)
ax.set_ylabel("Runime (s)")
if variation == "k":
    ax.set_xlabel("k\n|P|")
elif variation == "tau":
    ax.set_xlabel("τ" + "\n" + "|P|", fontsize=19)
elif variation == "n":
    ax.set_xlabel("n" + "\n" + "|P|")
else:
    ax.set_xlabel("# sensitive patterns\n|P|")

plt.xticks(fontsize=14)
# ax.set_title(metric.capitalize() + " vs " + variation + " - " + dataset)
ax.legend()


# plt.legend(loc="upper center", bbox_to_anchor=(0.5, 1.5), ncol=4)
# plt.legend().remove()
plt.legend(prop={"size": 14})
plt.tight_layout()
plt.gray()
# plt.show()
fig.savefig(
    file_name.split(".")[0] + "_" + variation + "_" + metric + ".svg",
    dpi=fig.dpi,
    bbox_inches="tight",
)
