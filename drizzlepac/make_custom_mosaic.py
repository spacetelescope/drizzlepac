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
import glob
import math
import pdb
import os
import sys
import tempfile
import traceback

from astropy.io import fits
import numpy as np

from drizzlepac import hapmultisequencer
from drizzlepac.haputils import cell_utils
from drizzlepac.haputils import make_poller_files

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

def calc_skycell_dist(x, y, x_ref, y_ref):
    """Calculate distance from one skyframe to another

    Parameters
    ----------
    x : int
        skyframe x index

    y : int
        skyframe y index

    x_ref : int
        reference x index

    y_ref : int
        reference y index

    Returns
    -------
    dist : float
        distance from (x, y) to (x_ref, y_ref)
    """
    dist = math.sqrt(float((x - x_ref)**2) + float((y - y_ref)**2))
    return(dist)
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


def create_poller_file(img_list, proj_cell_dict):
    """Subroutine that executes make_poller_files.generate_poller_file() to generate custom MVM poller file
    for execution of hapmultisequencer.run_mvm_processing().

    Parameters
    ----------
    img_list : list
        list of images to process

    proj_cell_dict : dictionary
        Dictionary containing projection cell information

    Returns
    -------
    poller_filename : str
        Name of the newly created poller_filename
    """
    # Locate bottom-left most skycell
    pcell_filename = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'pars', 'allsky_cells.fits')
    sc_nxy = fits.getval(pcell_filename, "SC_NXY",
                         ext=0)  # Projection cell side length (in skycells) in X or Y
    closest_skycell = ""
    min_dist = calc_skycell_dist(sc_nxy, sc_nxy, 1, 1) + 2.0
    for sc in proj_cell_dict.keys():
        x_index = proj_cell_dict[sc].x_index
        y_index = proj_cell_dict[sc].y_index

        dist = calc_skycell_dist(x_index, y_index, 1, 1)
        if dist < min_dist:
            min_dist = dist
            closest_skycell = sc

    # Create poller file
    skycell_name = closest_skycell.replace("skycell-", "")
    poller_filename = "temp_{}_mvm.out".format(skycell_name)
    rootname_list = []
    for imgname in img_list:
        rootname_list.append(imgname.split("_")[0]+"\n")
    tf = tempfile.NamedTemporaryFile(mode='w+t', dir=os.getcwd())
    with open(tf.name, 'w') as f:
        f.writelines(rootname_list)
        f.close()

        make_poller_files.generate_poller_file(tf.name, input_file_path=os.getcwd(),
                                               poller_file_type="mvm",
                                               output_poller_filename=poller_filename,
                                               skycell_name=skycell_name)
    return poller_filename

# ------------------------------------------------------------------------------------------------------------

def compute_mosaic_limits(proj_cell_dict):
    """Compute min and max limits of the rectangle that encloses the mosaic observations

    Parameters
    ----------
    proj_cell_dict : dictionary
        Dictionary containing projection cell information

    Returns
    -------
    mosaic_limits : list
        4-element list containing the mosaic bounding rectangle X min and max and Y min and max values
    """
    # set up storage arrays
    array_size = len(proj_cell_dict.keys()) * 4
    ra_values = np.empty(array_size)
    dec_values = np.empty(array_size)
    x_values = np.empty(array_size)
    y_values = np.empty(array_size)

    # Use wcs.calc_footprint() to get skycell corners, convert RA, Dec to X, Y in projection cell frame of reference
    i = 0
    for sc_name in proj_cell_dict.keys():
        for ra_dec in proj_cell_dict[sc_name].wcs.calc_footprint().tolist():
            x_y = proj_cell_dict[sc_name].projection_cell.wcs.all_world2pix([ra_dec], 0).tolist()[0]
            ra_values[i] = ra_dec[0]
            dec_values[i] = ra_dec[1]
            x_values[i] = x_y[0]
            y_values[i] = x_y[1]
            i += 1

    # Return mosaic bounding rectangle X min and max and Y min and max values
    mosaic_limits = [x_values.min(), x_values.max(), y_values.min(), y_values.max()]
    limit_labels = ["X_min", "X_max", "Y_min", "Y_max"]
    for i in range(0, 4):
        mosaic_limits[i] = int(np.rint(mosaic_limits[i]))
        print("{}: {}".format(limit_labels[i], mosaic_limits[i]))
    return mosaic_limits
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
        sc_nxy = fits.getval(pcell_filename, "SC_NXY", ext=0)  # Projection cell side length (in skycells) in X or Y
        min_dist = calc_skycell_dist(sc_nxy, sc_nxy, 1, 1) + 2.0
        best_pc = ""
        for pc in proj_cell_dict.keys():
            dist_ra = np.empty(len(proj_cell_dict[pc].keys()))
            for i, skycell_name in zip(range(0, len(proj_cell_dict[pc].keys())), proj_cell_dict[pc].keys()):
                pc_center_length = math.ceil(21/2.0)
                x_index = proj_cell_dict[pc][skycell_name].x_index
                y_index = proj_cell_dict[pc][skycell_name].y_index
                dist_ra[i] = calc_skycell_dist(x_index, y_index, pc_center_length, pc_center_length)
            if min_dist > dist_ra.mean():
                min_dist = dist_ra.mean()
                best_pc = pc
    log.info("Output WCS will be based on WCS from projection cell {}".format(best_pc))

    return proj_cell_dict[best_pc]

# ------------------------------------------------------------------------------------------------------------


def perform(input_image_source, log_level='info'):
    """Main calling subroutine

    Parameters
    ----------
    input_image_source : str
        Search pattern to be used to identify images to process or the name of a text file containing a list
        of images to process.

    log_level : str, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the
        .log file. The level of verboseness from left to right, and includes all log statements with a
        log_level left of the specified level. Specifying "critical" will only record/display "critical" log
        statements, and specifying "error" will record/display both "error" and "critical" log statements,
        and so on. Unless explicitly set, the default value is 'info'.

    Returns
    -------
    return_value : int
        A simple status value. '0' for a successful run or '1' for a failed run.
    """
    try:
        # optimistically pre-set return value to 0.
        return_value = 0
        log_level_dict = {"critical": logutil.logging.CRITICAL,
                          "error": logutil.logging.ERROR,
                          "warning": logutil.logging.WARNING,
                          "info": logutil.logging.INFO,
                          "debug": logutil.logging.DEBUG}
        logging_level = log_level_dict[log_level]
        log.setLevel(logging_level)
        temp_files_to_delete = []
        # Get list input fits files from input args, and raise an exception if no input images can be found.
        img_list = create_input_image_list(input_image_source)
        if not img_list:
            err_msg = "ERROR: No input images were found. Please double-check the search pattern or contents of the input list text file."
            log.critical(err_msg)
            raise Exception(err_msg)

        # get list of skycells/projection cells that observations are in
        # figure out which projection cell center is closest to the center of the observations, use that projection cell as basis for WCS
        proj_cell_dict = determine_projection_cell(img_list)

        # compute bounding rectangle limits for mosaic observations
        custom_limits = compute_mosaic_limits(proj_cell_dict)

        # Create MVM poller file
        poller_filename = create_poller_file(img_list, proj_cell_dict)
        temp_files_to_delete.append(poller_filename)

        # Execute hapmultisequencer.run_mvm_processing() with poller file
        log.debug("Creating custom mosaic from the following {} input images".format(len(img_list)))
        for item in img_list:
            log.debug(" {}".format(item))
        log.info("Mosaic bounding box limits")
        for limit_name, limit_value in zip(["X_min", "X_max", "Y_min", "Y_max"], custom_limits):
            log.debug("{}: {}".format(limit_name, int(np.rint(limit_value))))
        return_value = hapmultisequencer.run_mvm_processing(poller_filename,
                                                            custom_limits=custom_limits,
                                                            log_level=logging_level)

    except Exception:
        if return_value == 0:
            return_value = 1
        print("\a\a\a")
        exc_type, exc_value, exc_tb = sys.exc_info()
        traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)
    finally:
        if temp_files_to_delete and log_level != "debug":
            log.info("Time delete some temporary files...")
            for filename in temp_files_to_delete:
                if os.path.exists(filename):
                    os.remove(filename)
                    log.info("Removed temporary file {}.".format(filename))


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
    parser.add_argument('-l', '--log_level', required=False, default='info',
                        choices=['critical', 'error', 'warning', 'info', 'debug'],
                        help='The desired level of verboseness in the log statements displayed on the screen '
                             'and written to the .log file. The level of verboseness from left to right, and '
                             'includes all log statements with a log_level left of the specified level. '
                             'Specifying "critical" will only record/display "critical" log statements, and '
                             'specifying "error" will record/display both "error" and "critical" log '
                             'statements, and so on.')
    user_args = parser.parse_args()

    rv = perform(user_args.input_image_source, user_args.log_level)
# ------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    main()
