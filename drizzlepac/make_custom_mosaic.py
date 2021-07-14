#!/usr/bin/env python

"""
MAIN DOCSTRING GOES HERE!
"""
# TODO: Add main docstring

import argparse
import datetime
import glob
import logging
import pdb
import os
import sys


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
    # Assume user specifed a search pattern
    else:
        img_list = glob.glob(user_input)

    pdb.set_trace()
    return img_list

# ------------------------------------------------------------------------------------------------------------


def perform(input_image_source, output_wcs_source=None):
    """Main calling subroutine

    Parameters
    ----------
    input_image_source : str
        Search pattern to be used to identify images to process or the name of a text file containing a list
        of images to process.

    output_wcs_source : str, optional.
        Name of a file that contains the desired world coordinate system (WCS) information for the final
        output mosaic product(s). Users should use one of the following valid filetypes: 1) a calibrated fits
        image file (_flt.fits or _flc.fits fits) 2) a headerlet fits file (_hlet.fits) or 3) a text file
        containing all values of a complete WCS solution, one value per line, in the following specific
        order: CRVAL terms, CRPIX terms, Platescale value, and Orientation. If a source file for the output
        WCS information is not explicitly specified, the WCS information of the skycell containing the
        largest fraction of the input observations will be used.'

    Returns
    -------
    return_value : int
        A simple status value. '0' for a successful run or '1' for a failed run.
    """
    # optimistically pre-set return value to 0.
    return_value = 0

    # Get list input fits files from input args
    img_list = create_input_image_list(input_image_source)

    # Generate WCS object based on user-specified WCS info source file (or lack there of)


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
    parser.add_argument('-i', '--output_wcs_source', required=False, default=None,
                        help='Name of a file that contains the desired world coordinate system (WCS) '
                             'information for the final output mosaic product(s). Users should use one of '
                             'the following valid filetypes: 1) a calibrated fits image file (_flt.fits or '
                             '_flc.fits fits) 2) a headerlet fits file (_hlet.fits) or 3) a text file '
                             'containing all values of a complete WCS solution, one value per line, in the '
                             'following specific order: CRVAL terms, CRPIX terms, Platescale value, and '
                             'Orientation. If a source file for the output WCS information is not '
                             'explicitly specified, the WCS information of the skycell containing the '
                             'largest fraction of the input observations will be used.')
    user_args = parser.parse_args()

    rv = perform(user_args.input_images, output_wcs_soruce = user_args.output_wcs_source)
# ------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    main()
