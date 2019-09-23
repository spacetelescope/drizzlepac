#!/usr/bin/env python
import fileinput
import pdb
import sys
"""Simple little script to reset source flags in the last column of HLA classic and HAP sorucelists to zero"""



for line in fileinput.input(sys.argv[1], inplace=True):
    if not line.endswith(",0\n") and not line[0].isdigit():
        line=line.replace(line.split(",")[-1],"0\n")
        print(line.strip())
    else:
        print(line.strip())


# pdb.set_trace()
