#!/usr/bin/env python
"""Reports the name(s) of the projection cell(s)/skycell(s) that the observations in the user-specified input
 file occupy."""

import argparse
import glob
import os
import sys

from drizzlepac.haputils import cell_utils


__version__ = 0.1
__version_date__ = '08-Sept-2021'
# ------------------------------------------------------------------------------------------------------------


def report_skycells(img_list):
    """reports the names of the projection cell(s)/skycell(s) that the observations in the user-specified
    input list occupy.

    Parameters
    ----------
    img_list : list
        list of the images to process

    Returns
    -------
    Nothing!
    """
    skycell_dict = cell_utils.get_sky_cells(img_list)
    print("\n")
    for skycell_name in skycell_dict.keys():
        print("Skycell {} contains {} image(s)".format(skycell_name, len(skycell_dict[skycell_name].members)))
        for imgname in skycell_dict[skycell_name].members:
            print("     {}".format(imgname))
        print("\n")

# ------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='identifies the projection cell(s) and skycell(s) that the '
                                                 'user-specified observations occupy')
    parser.add_argument('-i', '--input_list', required=False, default="NONE",
                        help='Name of a file containing a list of calibrated fits files (ending with '
                             '"_flt.fits" or "_flc.fits") to process. If not explicitly specified, all '
                             'flc/flt fits files in the current working directory will be processed.')
    in_args = parser.parse_args()
    if in_args.input_list is "NONE":
        img_list = glob.glob("*_fl?.fits")
        if len(img_list) > 0:
            report_skycells(img_list)
        else:
            sys.exit("ERROR: No flc/flt fits files found in current path {}/".format(os.getcwd()))
    elif os.path.exists(in_args.input_list):
        with open(in_args.input_list) as fin:
            img_list = fin.read().splitlines()
        report_skycells(img_list)
    else:
        sys.exit("ERROR: Input file {} does not exist!".format(in_args.input_list))
