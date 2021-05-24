import sys

name = sys.argv[1]
outname = sys.argv[2]
in_file = open(name, "r")
out_file = open(outname, "w")
parameter_file = open("data/parameters/" + outname.split("/")[-1], "w")

nb_hash = 0
if name.split(".")[-1] == "txt":
    txt = in_file.read()
    nb_hash = txt.count("#")
    out_file.write(txt)
else:
    assert name.split(".")[-1] in ["fasta", "fa"]
    while True:
        name = in_file.readline().rstrip()
        seq = in_file.readline().rstrip()
        if not seq:
            break
        if "N" in seq:
            nb_hash += 1
            hash_seq = seq.replace("N", "#")
            out_file.write(hash_seq)

param = [0, 0, 0, nb_hash]
for p in param:
    parameter_file.write(str(p) + "\n")
parameter_file.close()

out_file.close()
in_file.close()
