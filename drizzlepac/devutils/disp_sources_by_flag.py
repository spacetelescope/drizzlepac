#!/usr/bin/env python
"""
This script breaks down the flag values in the specified sourcelist by flag bit value and plots source locations by
bit value

NOTE: This script requires the "pyds9" library.
"""
import sys
# import pdb
from astropy.io import fits as fits
from astropy.table import Table
import numpy as np
import pyds9

bit_list = [0, 1, 2, 4, 8, 16, 32, 64, 128]
flag_meanings = ['Point Source', 'Extended Source', 'Single-Pixel Saturation', 'Multi-Pixel Saturation',
               'Faint Magnitude Limit', 'Hot Pixel', 'Swarm Detection', 'Edge and Chip Gap',
               'Bleeding and Cosmic Rays']

# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-


def display_regions(imgname, reg_dict, flag_counts, n_sources):
    """
    Display the input input image overplot the positions of flagged sources

    parameters
    ----------
    imgname : string
        input image name

    reg_dict : dictionary
        dictionary containing region info keyed by flag bit

    flag_counts : list
        overall flag population stats broken down by individual flag bit

    n_sources : int
        total number of sources in sourcelist

    returns
    -------
    nothing.
    """
    imghdu = fits.open(imgname)
    d = pyds9.DS9()

    for ctr in range(0, len(bit_list)):
        bit_val = bit_list[ctr]
        # padding0 = 2 - len(str(ctr + 1))
        padding1 = 6 - len(str(bit_val))
        padding2 = 27 - len(flag_meanings[ctr])
        print("Frame {}: Bit value {}{}{}{}{} sources flagged ({}% of all sources)".format(ctr + 1,
                            bit_val, "." * padding1,
                            flag_meanings[ctr],
                            padding2 * ".",
                            flag_counts[ctr],
                            100.0 * (float(flag_counts[ctr]) / float(n_sources))))
        if ctr != 0:
            d.set("frame new")
        d.set_fits(imghdu)
        d.set("scale zscale")
        if reg_dict[bit_val]:
            d.set('regions', 'physical; {}'.format(reg_dict[bit_val]))
    # d.set("frame first")
    # d.set("frame delete")


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
    bitlist = [1, 2, 4, 8, 16, 32, 64, 128]
    flagval = int(flagval)
    # out_bit_list = []
    out_idx_list = np.zeros(9, dtype=int)
    if flagval == 0:
        # out_bit_list = [0]
        out_idx_list[0] = 1
    if flagval > 0:
        idx = 1
        for bit in bitlist:
            if flagval & bit > 0:
                # out_bit_list.append(bit)
                out_idx_list[idx] = 1
            if bit > flagval: break
            idx += 1
    return out_idx_list


# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~


def make_regions(sl_name):
    """
    generate a dictionary of dictionary region info to plot in ds9.

    parameters
    ----------
    sl_name : string
        sourcelist name

    returns
    -------
    reg_dict : dictionary
        dictionary containing region info keyed by flag bit

    flag_counts : list
        overall flag population stats broken down by individual flag bit

    n_sources : int
        total number of sources in sourcelist
    """
    table_data = Table.read(sl_name, format='ascii.ecsv')
    flag_counts = np.zeros(9, dtype=int)
    reg_dict = {}
    for bit_val in bit_list:
        reg_dict[bit_val] = ""

    for table_line in table_data:
        if sl_name.endswith("point-cat.ecsv"):
            x = table_line["X-Center"]
            y = table_line["Y-Center"]
        elif sl_name.endswith("segment-cat.ecsv"):
            x = table_line["X-Centroid"]
            y = table_line["Y-Centroid"]
        else:
            sys.exit("ERROR! Unrecognized catalog filetype!")
        flagval = table_line["Flags"]
        flag_bits = deconstruct_flag((flagval))
        flag_counts += flag_bits
        for bit_val, flag_element in zip(bit_list, flag_bits):
            if flag_element == 1:
                reg_dict[bit_val] += "point({}, {}) #point=circle color=red; \n".format(x + 1.0, y + 1.0)

    for bit_val in bit_list:
        if len(reg_dict[bit_val]) > 0:
            reg_dict[bit_val] = reg_dict[bit_val][:-2]
    n_sources = len(table_data)
    return reg_dict, flag_counts, n_sources


# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
def write_regions(reg_dict, sl_name):
    """Write out separate regions files for each bit flagged

    The region files written out can be read in directly by DS9 by the user, as an interactive
    alternative to using the Python interface pyds9.

    parameters
    -----------
    reg_dict : dict
        dictionary containing region info keyed by flag bit

    sl_name : string
        sourcelist name

    """
    for bit in reg_dict:
        fname = sl_name.replace("-cat.ecsv", "-bit{:03d}.reg".format(int(bit)))
        regions = reg_dict[bit]
        with open(fname, 'w') as bitfile:
            bitfile.writelines(regions)

# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
if __name__ == "__main__":
    imgname = sys.argv[1]
    slname = sys.argv[2]


    reg_dict, flag_counts, n_sources = make_regions(slname)

    display_regions(imgname, reg_dict, flag_counts, n_sources)
