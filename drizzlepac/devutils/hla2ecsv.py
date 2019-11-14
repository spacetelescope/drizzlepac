#!/usr/bin/env python


import sys
from astropy.table import Table
import numpy as np

infilename = sys.argv[1]
outfilename = infilename.replace(".txt",".ecsv")
source_cat = Table.read(infilename, format='ascii')

source_cat.write(outfilename, format="ascii.ecsv")
