#!/usr/bin/env python

""" runsinglehap.py - Module to control processing of single-visit mosaics

:License: :doc:`LICENSE`

USAGE: runsinglehap [-d] inputFilename

The '-d' option will run this task in DEBUG mode producing additional outputs.

Python USAGE:
    python
    from drizzlepac import runsinglehap
    runsinglehap.perform(inputFilename,debug=True)

"""
# Import standard Python modules
import argparse

import sys
import logging

# THIRD-PARTY
from stsci.tools import logutil

from drizzlepac import hapsequencer

__taskname__ = "runsinglehap"

# Local variables
__version__ = "0.1.1"
__version_date__ = "(16-Oct-2019)"
#
# These lines (or something similar) will be needed in the HAP processing code
#
MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout, 
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)
# Any module which uses 'util.with_logging' should be added separately here...
# logging.getLogger('astrodrizzle').addHandler(log)
# logging.getLogger('alignimages').addHandler(log)

# ----------------------------------------------------------------------------------------------------------------------

def perform(input_filename, **kwargs):
    """
    Main calling subroutine.

    Parameters
    ----------
    input_filename : string
        Name of the input csv file containing information about the files to
        be processed

    debug : Boolean
        display all tracebacks, and debug information?

    Updates
    -------
    return_value : list
        a simple status value. '0' for a successful run and '1' for a failed
        run
    """
    log.info("Starting single-visit processing of {}".format(input_filename))
    return_value = hapsequencer.run_hap_processing(input_filename, **kwargs)
    return return_value

# ----------------------------------------------------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description='Process images, produce drizzled images and sourcelists')
    parser.add_argument('input_filename', help='Name of the input csv file containing information about the files to '
                        'be processed')
    parser.add_argument('-d', '--diagnostic_mode', required=False, action='store_true', help='If this option is turned '
                        'on, additional log messages will be displayed and additional files will be created during the '
                        'course of the run.')
    user_args = parser.parse_args()

    print("Single-visit processing started for: {}".format(user_args.input_filename))
    rv = perform(user_args.input_filename, diagnostic_mode=user_args.diagnostic_mode)
    print("Return Value: ", rv)
    return rv

if __name__ == '__main__':
    main()
