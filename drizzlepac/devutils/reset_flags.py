#!/usr/bin/env python
import fileinput
import pdb
import sys
"""Simple little script to reset source flags in the last column of HLA classic and HAP sorucelists to zero"""



for line in fileinput.input(sys.argv[1], inplace=True):
    if sys.argv[1].endswith("daophot.txt"):
        if not line.startswith("X-Center") and not line.endswith(",0\n"):
            line=line.replace(line.split(",")[-1],"0\n")
            print(line.strip())
        else:
            print(line.strip())
    if sys.argv[1].endswith("sexphot.txt"):
        if not line.startswith("X_IMAGE") and not line.endswith(",0\n"):
            line=line.replace(line.split(",")[-1],"0\n")
            print(line.strip())
        else:
            print(line.strip())


# pdb.set_trace()
