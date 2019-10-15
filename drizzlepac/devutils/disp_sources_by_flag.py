#!/usr/bin/env python

import pdb
import sys

from astropy.io import fits as fits
from astropy.table import Table
import numpy as np
import pyds9



# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
def deconstruct_flag(flagval):
    """Breaks down an integer flag value into individual component bit values.

    Parameters
    ----------
    flagval : int
        Flag value to deconstruct

    Returns
    -------
    out_idx_list : list
        a 9-element numpy array of 0s and 1s. Each element of the array represents the presence of a particular
        bit value (element 0 = bit 0, element 1 = bit 1, ..., element 3 = bit 4 and so on...)
    """
    bit_list = [1, 2, 4, 8, 16, 32, 64, 128]
    flagval = int(flagval)
    # out_bit_list = []
    out_idx_list = np.zeros(9, dtype=int)
    if flagval == 0:
        # out_bit_list = [0]
        out_idx_list[0] = 1
    if flagval > 0:
        idx = 1
        for bit in bit_list:
            if flagval & bit > 0:
                # out_bit_list.append(bit)
                out_idx_list[idx] = 1
            if bit > flagval: break
            idx += 1
    return (out_idx_list)

# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-

def make_region_files(sl_name):
    table_data = Table.read(sl_name, format='ascii.ecsv')
    for table_line in table_data:
        x = table_line[0]
        y = table_line[1]
        flagval  = table_line[-1]
        flag_bits = deconstruct_flag((flagval))
        print(table_line,x,y,flagval,flag_bits)
# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
def run(imgname,slanme):
    make_region_files(slname)





# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
if __name__ == "__main__":
    imgname = "hst_10265_01_acs_wfc_f606w_j92c01_drc.fits"#sys.argv[1]
    slname = "hst_10265_01_acs_wfc_f606w_j92c01_point-cat.ecsv"#sys.argv[2]

    bit_list = [0, 1, 2, 4, 8, 16, 32, 64, 128]
    run(imgname,slname)