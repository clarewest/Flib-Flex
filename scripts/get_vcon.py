import sys
import csv

def writeOut(listoflists, outfile):
    with open(outfile, "w") as f:
        wr = csv.writer(f, delimiter=" ")
        wr.writerows(listoflists)

target=sys.argv[1]
terminus=sys.argv[3]
#length=int(sys.argv[3])
missing=int(sys.argv[2])

confile=target+".metapsicov_stage1"
vconfile=target+".vcon"
fasta=target+".fasta"
with open(fasta, "r") as fin:
    fin.readline()
    length=len(fin.readline().strip())

vcons = []

if terminus == "C":
    begin = 1
    end = length - missing + 1
elif terminus == "N":
    begin = missing + 1
    end = length 
else:
    print("Not a valid terminus")

with open(confile, "r") as incon:
    for predcon in incon:
        con = predcon.strip().split()
        if ((int(con[0]) or int(con[1])) not in range(begin, end+1)) or (int(con[3])==1):
            vcons.append(con[:3])
writeOut(vcons, vconfile)


