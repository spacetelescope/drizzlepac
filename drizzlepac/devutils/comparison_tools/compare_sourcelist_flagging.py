#!/usr/bin/env python
"""
This script allows the user to compare the flagging in two filter sourcelists on a bit-by-bit basis. It breaks down
the flag values in the specified sourcelists by flag bit value and plots source locations by bit value

NOTE: This script requires the "pyds9" library.
"""
import argparse
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


def display_regions(imgname, coordsys, reg_dict_list, flag_counts_list, n_sources_list):
    """
    Display the input input image overplot the positions of flagged sources

    parameters
    ----------
    imgname : string
        input image name

    coordsys : string
        coordinate system to use

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
    print("{}Catalog #1{}Catalog #2".format(" " * 56, " " * 10))
    for ctr in range(0, len(bit_list)):
        bit_val = bit_list[ctr]
        padding1 = 6 - len(str(bit_val))
        padding2 = 27 - len(flag_meanings[ctr])


        text_stuff = "Frame {}: Bit value {}{}{}{}".format(ctr+1,bit_val,"."*padding1,flag_meanings[ctr],padding2*".")
        pct_1 = (float(flag_counts_list[0][ctr])/float(n_sources_list[0]))*100.0
        pct_2 = (float(flag_counts_list[1][ctr]) / float(n_sources_list[1])) * 100.0
        out_string = "{:9d} {:8.3f}% {:9d} {:8.3f}%".format(flag_counts_list[0][ctr],pct_1,flag_counts_list[1][ctr],pct_2)

        print("{}{}".format(text_stuff,out_string.replace(" ",".")))
        if ctr != 0:
            d.set("frame new")
        d.set_fits(imghdu)
        d.set("scale zscale")
        if coordsys == "radec":
            d.set("wcs fk5")
            d.set("wcs degrees")
        for reg_dict in reg_dict_list:
            if reg_dict[bit_val]:
                if coordsys == "xy":
                    d.set('regions', 'physical; {}'.format(reg_dict[bit_val]))
                if coordsys == "radec":
                    d.set('regions', 'fk5; {}'.format(reg_dict[bit_val]))
    out_string = ".......{:9d}{}{:9d}".format(n_sources_list[0], "." * 11 , n_sources_list[1])
    print("{}TOTAL CATALOG LENGTH{}".format(" " * 25,out_string.replace(" ",".")))
    d.set("lock frame image")
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


def make_regions(sl_name, coordsys, shape, color):
    """
    generate a dictionary of dictionary region info to plot in ds9. Assumes that X and Y coords are stored the first
    two (leftmost) columns, and flag values are stored in the last (rightmost) column.

    parameters
    ----------
    sl_name : string
        sourcelist name

    coordsys : string
        coordinate system to use (either "xy" or radec)

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
    if sl_name.endswith(".ecsv"):
        table_data = Table.read(sl_name, format='ascii.ecsv')
    elif (sl_name.endswith("daophot.txt") or sl_name.endswith("daophot_corrected.txt")):
        table_data = Table.read(sl_name, format='ascii')
    elif (sl_name.endswith("sexphot.txt") or sl_name.endswith("sexphot_corrected.txt")):
        table_data = Table.read(sl_name, format='ascii')
    else:
        sys.exit("ERROR! Unrecognized catalog filetype!")
    flag_counts = np.zeros(9, dtype=int)
    reg_dict = {}
    for bit_val in bit_list:
        reg_dict[bit_val] = ""

    table_col_titles = {
        "point-cat.ecsv": {"xy": ("X-Center", "Y-Center", "Flags"),
                           "radec": ("RA", "DEC", "Flags")},
        "segment-cat.ecsv": {"xy": ("X-Centroid", "Y-Centroid", "Flags"),
                             "radec": ("RA", "DEC", "Flags")},
        "daophot.txt": {"xy": ("X-Center", "Y-Center", "Flags"),
                        "radec": ("RA", "DEC", "Flags")},
        "daophot_corrected.txt": {"xy": ("X-Center", "Y-Center", "Flags"),
                                  "radec": ("RA", "DEC", "Flags")},
        "sexphot.txt": {"xy": ("X_IMAGE", "Y_IMAGE", "FLAGS"),
                        "radec": ("RA", "DEC", "FLAGS")},
        "sexphot_corrected.txt": {"xy": ("X_IMAGE", "Y_IMAGE", "FLAGS"),
                                  "radec": ("RA", "DEC", "FLAGS")}}
    found_ending = False
    fileendinglist = ""
    for filenameending in table_col_titles.keys():
        fileendinglist += filenameending+"\n"
        if sl_name.endswith(filenameending):
            x_col_name = table_col_titles[filenameending][coordsys][0]
            y_col_name = table_col_titles[filenameending][coordsys][1]
            flag_col_name = table_col_titles[filenameending][coordsys][2]
            found_ending = True
            break
    if not found_ending:
        errmsg = "'{}' is not a valid filetype. Valid filetypes are as follows:\n{}".format(sl_name,fileendinglist)
        log.error(errmsg)
        raise Exception(errmsg)
    for table_line in table_data:
        x = table_line[x_col_name]
        y = table_line[y_col_name]
        flagval = table_line[flag_col_name]
        flag_bits = deconstruct_flag((flagval))
        flag_counts += flag_bits
        for bit_val, flag_element in zip(bit_list, flag_bits):
            if flag_element == 1:
                if coordsys == "xy":
                    reg_dict[bit_val] += "point({}, {}) #point={} color={}; ".format(x+1.0, y+1.0, shape, color)
                if coordsys == "radec":
                    reg_dict[bit_val] += "point({}, {}) #point={} color={}; ".format(x, y, shape, color)
    for bit_val in bit_list:
        if len(reg_dict[bit_val]) > 0:
            reg_dict[bit_val] = reg_dict[bit_val][:-2]
    n_sources = len(table_data)
    return reg_dict, flag_counts, n_sources

# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description='Compare Sourcelists')
    # required positional input arguments
    PARSER.add_argument('sourcelistNames', nargs=2,help='A space-separated pair of sourcelists to compare. The first sourcelist is assumed to be the reference sourcelist that the second is being compared to.')
    PARSER.add_argument('-c', '--coordsys',required=False,default='xy',choices=['xy','radec'],help='Coordinate system to use for plotting sources. Choices are either "xy" for simple x-y coords, or "radec" for fk5 right ascention and declination values in degrees. Default is "xy".')
    PARSER.add_argument('-i', '--imageName', required=True, help='Image to overplot flags on in ds9')
    ARGS = PARSER.parse_args()

    imgname = ARGS.imageName
    sl_name_a = ARGS.sourcelistNames[0]
    sl_name_b = ARGS.sourcelistNames[1]
    coordsys = ARGS.coordsys

    if sl_name_a.split("_")[5] == "total" or sl_name_b.split("_")[5] == "total":
        sys.exit("Invalid catalog input. This script only compares filter catalogs, not total catalogs.")
    shapes = ['box', 'circle']
    colors = ['red', 'green']

    reg_dict_list = ["", ""]
    flag_counts_list = ["", ""]
    n_sources_list = ["", ""]
    print("Creating regions for catalog " + sl_name_a)
    reg_dict_list[0], flag_counts_list[0], n_sources_list[0] = make_regions(sl_name_a, coordsys, shapes[0], colors[0])
    print("Creating regions for catalog " + sl_name_b+"\n")
    reg_dict_list[1], flag_counts_list[1], n_sources_list[1] = make_regions(sl_name_b, coordsys, shapes[1], colors[1])
    len_list = [len("{} {}".format(colors[0], shapes[0])), len("{} {}".format(colors[1], shapes[1]))]
    max_width = max(len_list)
    print("Catalog #1 ({} {}):{} {}".format(colors[0], shapes[0], (max_width - len_list[0])*" ", sl_name_a))
    print("Catalog #2 ({} {}):{} {}".format(colors[1], shapes[1], (max_width - len_list[1])*" ", sl_name_b))

    display_regions(imgname, coordsys, reg_dict_list, flag_counts_list, n_sources_list)
