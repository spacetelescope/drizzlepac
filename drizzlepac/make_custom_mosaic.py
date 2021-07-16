#!/usr/bin/env python

""" make_custom_mosaic.py - Module to control processing of user-defined custom mosaics

USAGE:

- python drizzlepac/make_custom_mosaic.py <search pattern enclosed in quotes>
- python drizzlepac/make_custom_mosaic.py <text file with list of input files>

Python USAGE:
    python
    from drizzlepac import make_custom_mosaic
    make_custom_mosaic.perform(<list file or search pattern)
"""

import argparse
import datetime
import glob
import logging
import math
import pdb
import os
import sys

from astropy.io import fits
from astropy.table import Table
import numpy as np
import drizzlepac

from drizzlepac.haputils import cell_utils

from stsci.tools import logutil
from stwcs import wcsutil

__taskname__ = 'make_custom_mosaic'
MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)

__version__ = 0.1
__version_date__ = '14-July-2021'
# ------------------------------------------------------------------------------------------------------------

def create_input_image_list(user_input):
    """Create list of input images based in user input from command-line
    
    Parameters
    ----------
    user_input : str
        Search pattern to be used to identify images to process or the name of a text file containing a list 
        of images to process
    
    Returns
    -------
    img_list : list
        list of images to process
    """
    # Get files from list in user-specified text file
    if os.path.isfile(user_input) and not user_input.endswith(".fits"):
        with open(user_input, "r") as fin:
            img_list = fin.readlines()
        # clean up any stray carrage returns
        for ctr in range(0, len(img_list)):
            img_list[ctr] = img_list[ctr].strip()
        search_method_string = "in list file {}".format(user_input)
    # Assume user specified a search pattern
    else:
        img_list = glob.glob(user_input)
        search_method_string = "using search pattern '{}'".format(user_input)
    if len(img_list) == 1:
        plural_string = ""
    else:
        plural_string = "s"
    log.info("Found {} input image{} {}:".format(len(img_list), plural_string, search_method_string))
    for item in img_list:
        log.info("{}".format(item))
    return img_list

# ------------------------------------------------------------------------------------------------------------

def determine_projection_cell(img_list):
    """Determine which projection cell should be used as the basis for the WCS of the output mosaic
    product(s) based on which projection cell center is closest to the center of the observations.

    Parameters
    ----------
    img_list : list
        A list of images to process

    Returns
    -------
    Not sure yet: Unknown type
        I really just don't know at this point.
    """
    # Create dictionary containing information about the skycells that contain input images
    skycell_dict = cell_utils.get_sky_cells(img_list)

    # recast skycell_dict into nested dictionary with the top-most layer keyed by projection cell name
    proj_cell_dict = {}
    for key in skycell_dict.keys():
        proj_cell = key[9:13]
        if proj_cell not in proj_cell_dict.keys():
            proj_cell_dict[proj_cell] = {}
        proj_cell_dict[proj_cell][key] = skycell_dict[key]

    # Determine which skycell's WCS information should be used as the basis for WCS of the output product(s)
    if len(proj_cell_dict.keys()) == 1:
        log.info("Observations are present in only a single projection cell.")
        best_pc = list(proj_cell_dict)[0]
    else:
        log.info("Observations are present in multiple projection cells.")
        log.info("Output WCS will be based on WCS from the projection cell whose center is closest to the center of the input observations.")
        # Determine which projection cell's center is closest to the center of the observations
        pcell_filename = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'pars', 'allsky_cells.fits')
        sc_nxy = fits.getval(pcell_filename, "SC_NXY", ext=0) # Projection cell side length (in skycells) in X or Y
        min_dist = math.sqrt(sc_nxy**2 + sc_nxy**2) + 1.0
        best_pc = ""
        for pc in proj_cell_dict.keys():
            dist_ra = np.empty(len(proj_cell_dict[pc].keys()))
            for i, skycell_name in zip(range(0, len(proj_cell_dict[pc].keys())), proj_cell_dict[pc].keys()):
                pc_center_length = math.ceil(21/2.0)
                x_index = proj_cell_dict[pc][skycell_name].x_index
                y_index = proj_cell_dict[pc][skycell_name].y_index
                dist_ra[i] = math.sqrt((pc_center_length - x_index)**2 + (pc_center_length - y_index)**2)
            if min_dist > dist_ra.mean():
                min_dist = dist_ra.mean()
                best_pc = pc
    log.info("Output WCS will be based on WCS from projection cell {}".format(best_pc))


    return proj_cell_dict[best_pc]
# ------------------------------------------------------------------------------------------------------------


def perform(input_image_source):
    """Main calling subroutine

    Parameters
    ----------
    input_image_source : str
        Search pattern to be used to identify images to process or the name of a text file containing a list
        of images to process.

    Returns
    -------
    return_value : int
        A simple status value. '0' for a successful run or '1' for a failed run.
    """
    # optimistically pre-set return value to 0.
    return_value = 0
    log.setLevel(logging.DEBUG)
    # Get list input fits files from input args, and raise an exception if no input images can be found.
    img_list = create_input_image_list(input_image_source)
    if not img_list:
        err_msg = "ERROR: No input images were found. Please double-check the search pattern or contents of the input list text file."
        log.critical(err_msg)
        raise Exception(err_msg)

    # get list of skycells/projection cells that observations are in
    # figure out which projection cell center is closest to the center of the observations, use that projection cell as basis for WCS
    proj_cell_dict = determine_projection_cell(img_list)

    # Create MVM poller file

    # Generate custom MVM config .json file and insert relevant WCS info from proj_cell_dict

    # Execute hapmultisequencer.run_mvm_processing() with poller file, custom config file
    # TODO: PROBLEM: This will still cut off image at skycell boundry.


    # use cell_utils.bounded_wcs() to crop down image size to just around the mosaic footprint.

    pdb.set_trace()
    return return_value

# ------------------------------------------------------------------------------------------------------------


def main():
    """Command-line interface

    Parameters
    ----------
    None

    Returns
    -------
    Nothing.
    """
    # Parse command-line inputs
    parser = argparse.ArgumentParser(description='Create custom mosaic based on user-specified images and '
                                                 'world coordinate system (WCS) information')
    parser.add_argument('input_image_source',
                        help='Search pattern to be used to identify images to process (NOTE: Pattern must be '
                             'enclosed in single or double quotes) or alternately, the '
                             'name of a text file containing a list of images to process')
    user_args = parser.parse_args()

    rv = perform(user_args.input_image_source)
# ------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    main()
