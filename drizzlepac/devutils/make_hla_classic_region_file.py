#!/usr/bin/env python


import sys
from astropy.table import Table
import numpy as np


infilename = sys.argv[1]
outfilename = infilename.replace(".txt",".reg")
source_cat = Table.read(infilename, format='ascii')


out_table = source_cat.copy()

# Remove all other columns besides 'X-Center and Y-Center
out_table.keep_columns(['X-Center', 'Y-Center'])
# Add offset of 1.0 in X and Y to line up sources in region file with image displayed in ds9.
# out_table['X-Center'].data[:] += np.float64(1.0)
# out_table['Y-Center'].data[:] += np.float64(1.0)


out_table.write(outfilename, format="ascii")
print("Wrote region file '{}' containing {} sources".format(outfilename, len(out_table)))