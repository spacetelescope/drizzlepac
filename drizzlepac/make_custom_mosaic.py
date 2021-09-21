#!/usr/bin/env python

""" make_custom_mosaic.py - Module to control the generation of user-defined multi-skycell mosaics. This
script extends the capabilities of the HAP multi-visit mosaics (MVM) processing code (hapmultisequencer.py).
Based on user input, this script will produce drizzle-combined mosaics that will be as large as necessary to
incorporate all input images. Mosaic images can span multiple skycells, but NOT multiple projection cells. In
the exceedingly rare case that a custom mosaic contains observations that fall on multiple adjacent
projection cells, the bounds of the mosaic will be automatically clipped at the projection cell edge.

As this script is simply a wrapper around hapmultisequencer.py, the final output products are the exact same
as what one would expect from a pipeline MVM run. For a complete details about all hapmultisequencer.py
inputs and output products, please refer to the hapmultisequencer.py documentation. # TODO: Make hapmultisequencer.py documentation

hapmultisequencer.py  was explicitly written to process data from only a single skycell. However, this script
gives users the ability to create mosaics spanning multiple skycells. As such, all output products produced
by hapmultisequencer.py normally include both the projection cell and skycell names in the all output
products. To avoid any confusion that might arise during the creation of mosaics that span multiple skycells,
all custom mosaic output product filenames will by default include projection cell name and the projection
cell reference right ascension and declination instead of the skycell name. Users can also optionally
specify an output product filename prefix of their choosing.

The output world cooordinate system (WCS) information is based on that of the projection cell in which the
observations reside. If the input observations happen to fall in a region of the sky where more than one
projection cells overlap, the WCS information of the output products will be based on the projection cell
whose center is closest to the geometric center of the input observations.

For all detectors, output products will be generated with a plate scale of 0.04 arcseconds per pixel. The
script produces an additional set up output products for WFC3/IR observations. This additional set of
products is generated using the native detector platescale of 0.12 arcseconds per pixel. These files end in
"_coarse_all_flt.fits".

The output world cooordinate system (WCS) information is based on that of the projection cell in which the
observations reside. If the input observations happen to fall in a region of the sky where more than one
projection cells overlap, the WCS information of the output products will be based on the projection cell
whose center is closest to the geometric center of the input observations.

USAGE:

- python drizzlepac/make_custom_mosaic.py <search pattern enclosed in quotes>
- python drizzlepac/make_custom_mosaic.py <text file with list of input files>

Python USAGE:
    python
    from drizzlepac import make_custom_mosaic
    make_custom_mosaic.perform(<list file or search pattern>)
"""

import argparse
import glob
import math
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
        rootname_list.append(imgname+"\n")
    tf = tempfile.NamedTemporaryFile(mode='w+t', dir=os.getcwd())
    with open(tf.name, 'w') as f:
        f.writelines(rootname_list)
        f.close()

        make_poller_files.generate_poller_file(tf.name, poller_file_type="mvm",
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

    # determine min and max X and Y values, convert them from floating points to integers
    x_min = int(np.rint(x_values.min()))
    x_max = int(np.rint(x_values.max()))
    y_min = int(np.rint(y_values.min()))
    y_max = int(np.rint(y_values.max()))

    # Make sure all min and max X and Y values are within the projection cell. If not, reset to the projection cell limit.
    ps_x_max = proj_cell_dict[sc_name].projection_cell.wcs.pixel_shape[0]
    ps_y_max = proj_cell_dict[sc_name].projection_cell.wcs.pixel_shape[1]
    if x_min < 0:
        log.warning("Custom mosaic minimum X value '{}' is beyond the projection cell minimum X value of '{}' Resetting custom mosaic minimum X value to projection cell minimum X value.".format(x_min, 0))
        x_min = 0
    if x_max > ps_x_max:
        log.warning("Custom mosaic maximum X value '{}' is beyond the projection cell maximum X value of '{}' Resetting custom mosaic maximum X value to projection cell maximum X value.".format(x_max, ps_x_max))
        x_max = ps_x_max
    if y_min < 0:
        log.warning("Custom mosaic minimum Y value '{}' is beyond the projection cell minimum Y value of '{}' Resetting custom mosaic minimum Y value to projection cell minimum Y value.".format(y_min, 0))
        y_min = 0
    if y_max > ps_y_max:
        log.warning("Custom mosaic maximum Y value '{}' is beyond the projection cell maximum Y value of '{}' Resetting custom mosaic maximum Y value to projection cell maximum Y value.".format(y_max, ps_y_max))
        y_max = ps_y_max

    mosaic_limits = [x_min, x_max, y_min, y_max]

    if log.level <= 10:  # only display x and y limits if logging level is set to "debug".
        for mosaic_limit, limit_label in zip(mosaic_limits, ["minimum X", "maximum X", "minimum Y", "maximum Y"]):
            log.debug("Custom mosaic {} value (in projection cell pixels) : {}".format(limit_label, mosaic_limit))

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
    best_pc_dict : dict
        Dictionary of SkyCell objects for the input observations. These SkyCell objects contain projection
        cell information from the projection cell whose center is closest to the geometric center of the
        input observations.
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
        log.info("Output WCS will be based on WCS from the projection cell whose center is closest to the "
                 "center of the input observations.")
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

    best_pc_dict = proj_cell_dict[best_pc]
    return best_pc_dict

# ------------------------------------------------------------------------------------------------------------


def perform(input_image_source, log_level='info', output_file_prefix=None, skip_gaia_alignment=False):
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

    output_file_prefix : str, optional
        'Text string that will be used as the filename prefix all files created by hapmultisequencer.py
        during the MVM custom mosaic generation process. If not explicitly specified, all output files will
        start with the following formatted text string:
        "hst-skycell-p<pppp>-ra<##>d<####>-dec<n|s><##>d<####>", where p<pppp> is the projection cell ID,
        ra<##>d<####> are the whole-number and decimal portions of the right ascention, respectively, and
        dec<n|s><##>d<####> are the whole-number and decimal portions of the declination, respectively. Note
        that the "<n|s>" denotes if the declination is north (positive) or south (negative). Example: For
        skycell = 1974, ra = 201.9512, and dec = +26.0012, The filename prefix would be
        "skycell-p1974-ra201d9512-decn26d0012".

    skip_gaia_alignment : bool, optional
        Skip alignment of all input images to known Gaia/HSC sources in the input image footprint? If set to
        'True', The existing input image alignment solution will be used instead. The default is False.

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
                                                            log_level=logging_level,
                                                            output_file_prefix=output_file_prefix,
                                                            skip_gaia_alignment=skip_gaia_alignment)

    except Exception:
        if return_value == 0:
            return_value = 1
        print("\a\a\a")
        exc_type, exc_value, exc_tb = sys.exc_info()
        traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)
    finally:
        # remove any temp files like the poller file and
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
    return_value : int
        return value from the run. 0 for successful run, something else otherwise.
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
    parser.add_argument('-o', '--output_file_prefix', required=False,
                        help='Text string that will be used as the filename prefix all files created by '
                             'hapmultisequencer.py during the MVM custom mosaic generation process. If not '
                             'explicitly specified, all output files will start with the following formatted '
                             'text string: "hst-skycell-p<pppp>-ra<##>d<####>-dec<n|s><##>d<####>", where '
                             'p<pppp> is the projection cell ID, ra<##>d<####> are the whole-number and '
                             'decimal portions of the right ascention, respectively, and dec<n|s><##>d<####> '
                             'are the whole-number and decimal portions of the declination, respectively. '
                             'Note that the "<n|s>" denotes if the declination is north (positive) or south '
                             '(negative). Example: For skycell = 1974, ra = 201.9512, and dec = +26.0012, '
                             'The filename prefix would be "skycell-p1974-ra201d9512-decn26d0012".')
    parser.add_argument('-s', '--skip_gaia_alignment', required=False, action='store_true',
                        help='Skip alignment of all input images to known Gaia/HSC sources in the input '
                             'image footprint? If this option is turned on, the existing input image '
                             'alignment solution will be used instead. The default is False.')
    user_args = parser.parse_args()
    return_value = perform(user_args.input_image_source,
                           log_level=user_args.log_level,
                           output_file_prefix=user_args.output_file_prefix,
                           skip_gaia_alignment=user_args.skip_gaia_alignment)
# ------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    main()
