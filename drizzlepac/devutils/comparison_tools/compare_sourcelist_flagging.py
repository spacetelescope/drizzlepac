#!/usr/bin/env python
"""
This script allows the user to compare the flagging in two filter sourcelists on a bit-by-bit basis. It breaks down
the flag values in the specified sourcelists by flag bit value and plots source locations by bit value

NOTE: This script requires the "pyds9" library.
"""
import sys
from astropy.io import fits as fits
from astropy.table import Table
import numpy as np
import pyds9
import pdb
bit_list = [0, 1, 2, 4, 8, 16, 32, 64, 128]
flag_meanings = ['Point Source', 'Extended Source', 'Single-Pixel Saturation', 'Multi-Pixel Saturation',
                 'Faint Magnitude Limit', 'Hot Pixel', 'Swarm Detection', 'Edge and Chip Gap',
                 'Bleeding and Cosmic Rays']

# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-


def display_regions(imgname, reg_dict_list, flag_counts_list, n_sources_list):
    """
    Display the input input image overplot the positions of flagged sources

    parameters
    ----------
    imgname : string
        input image name

    reg_dict_list : list
        list of dictionaries containing region info keyed by flag bit for both sourcelists

    flag_counts_list : list
        list of overall flag population stats broken down by individual flag bit for both sourcelists

    n_sources_list : list
        list total number of sources in each sourcelist

    returns
    -------
    nothing.
    """
    imghdu = fits.open(imgname)
    d = pyds9.DS9()
    print("{}Cat. #1{}Cat. #2".format(" " * 52, " " * 3))
    for ctr in range(0, len(bit_list)):
        bit_val = bit_list[ctr]
        padding1 = 6 - len(str(bit_val))
        padding2 = 27 - len(flag_meanings[ctr])
        print("Frame {}: Bit value {}{}{}{}{}{}{}".format(ctr+1,
                                                          bit_val,
                                                          "."*padding1,
                                                          flag_meanings[ctr],
                                                          padding2*".",
                                                          flag_counts_list[0][ctr],
                                                          (10-len(str(flag_counts_list[0][ctr])))*".",
                                                          flag_counts_list[1][ctr]))
        if ctr != 0:
            d.set("frame new")
        d.set_fits(imghdu)
        d.set("scale zscale")
        for reg_dict in reg_dict_list:
            if reg_dict[bit_val]:
                d.set('regions', 'physical; {}'.format(reg_dict[bit_val]))

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
            if bit > flagval:
                break
            idx += 1
    return out_idx_list


# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~


def make_regions(sl_name, shape, color):
    """
    generate a dictionary of dictionary region info to plot in ds9. Assumes that X and Y coords are stored the first
    two (leftmost) columns, and flag values are stored in the last (rightmost) column.

    parameters
    ----------
    sl_name : string
        sourcelist name

    shape : string
        desired region shape

    color : string
        desired region color

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
        flagval  = table_line["Flags"]
        flag_bits = deconstruct_flag((flagval))
        flag_counts += flag_bits
        for bit_val, flag_element in zip(bit_list, flag_bits):
            if flag_element == 1:
                reg_dict[bit_val] += "point({}, {}) #point={} color={}; ".format(x+1.0, y+1.0, shape, color)

    for bit_val in bit_list:
        if len(reg_dict[bit_val]) > 0:
            reg_dict[bit_val] = reg_dict[bit_val][:-2]
    n_sources = len(table_data)
    return reg_dict, flag_counts, n_sources

# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~


if __name__ == "__main__":
    imgname = sys.argv[1]
    sl_name_a = sys.argv[2]
    sl_name_b = sys.argv[3]
    if sl_name_a.split("_")[5] == "total" or sl_name_b.split("_")[5] == "total":
        sys.exit("Invalid catalog input. This script only compares filter catalogs, not total catalogs.")
    shapes = ['box', 'circle']
    colors = ['red', 'green']

    reg_dict_list = ["", ""]
    flag_counts_list = ["", ""]
    n_sources_list = ["", ""]
    print("Creating regions for catalog " + sl_name_a)
    reg_dict_list[0], flag_counts_list[0], n_sources_list[0] = make_regions(sl_name_a, shapes[0], colors[0])
    print("Creating regions for catalog " + sl_name_b+"\n")
    reg_dict_list[1], flag_counts_list[1], n_sources_list[1] = make_regions(sl_name_b, shapes[1], colors[1])
    len_list = [len("{} {}".format(colors[0], shapes[0])), len("{} {}".format(colors[1], shapes[1]))]
    max_width = max(len_list)
    print("Catalog #1 ({} {}):{} {}".format(colors[0], shapes[0], (max_width - len_list[0])*" ", sl_name_a))
    print("Catalog #2 ({} {}):{} {}".format(colors[1], shapes[1], (max_width - len_list[1])*" ", sl_name_b))

    display_regions(imgname, reg_dict_list, flag_counts_list, n_sources_list)
