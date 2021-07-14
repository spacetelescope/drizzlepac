#!/usr/bin/env python

"""
MAIN DOCSTRING GOES HERE!
"""
# TODO: Add main docstring

import argpase
import datetime
import glob
import logging
import pdb
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


def perform():
    """Main calling subroutine

    Parameters
    ----------

    Returns
    -------
    return_value : int
        Return value.
    """

    return_value = 0  # TODO: Remove hard-wired return value
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
    parser = argparse.ArgumentParser(description='Create custom mosaic based on user-specified images and '
                                                 'world coordinate system (WCS) information')
    parser.add_argument('input_images',
                        help='Search pattern to be used to identify images to process or alternately, the '
                             'name of a text file containing a list of images to process')
    parser.add_argument('-i', '--output_wcs_info_source_file', required=False, default=None,
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

    rv = perform()
# ------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    main()
