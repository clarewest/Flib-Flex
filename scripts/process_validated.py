import sys
import re
import os

outfile=sys.argv[1]+"_proc"
try:
    os.remove(outfile)
except OSError:
    pass
newlib=[]
with open(sys.argv[1], "r") as f:
    for line in f:
        cols=line.strip().split()
        if float(cols[-1]) > 0:
            newline=("\t").join(cols[:-2]+["-"+cols[-1]])+"\n"
        else:
            newline=("\t").join(cols[:-1])+"\n"
        newlib.append(newline)
with open(outfile, "a") as f:
    for line in newlib:
        f.write(line)


