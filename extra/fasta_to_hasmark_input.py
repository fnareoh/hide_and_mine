import sys

name = sys.argv[1]
outname = sys.argv[2]
assert name.split(".")[-1] in ["fasta", "fa"]
in_file = open(name, "r")
out_file = open(outname, "w")

while True:
    name = in_file.readline().rstrip()
    seq = in_file.readline().rstrip()
    if not seq:
        break
    if "N" in seq:
        hash_seq = seq.replace("N", "#")
        out_file.write(hash_seq)


out_file.close()
in_file.close()
