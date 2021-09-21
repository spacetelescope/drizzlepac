#!/usr/bin/env python

""" runmultihap.py - Module to control processing of multi-visit mosaics

:License: :doc:`LICENSE`

USAGE: runmultihap [-d] inputFilename

The '-d' option will run this task in DEBUG mode producing additional outputs.

Python USAGE:
    python
    from drizzlepac import runmultihap
    runmultihap.perform(inputFilename,debug=True)

"""
# Import standard Python modules
import argparse
import sys

# THIRD-PARTY
from stsci.tools import logutil

from drizzlepac import hapmultisequencer

__taskname__ = "runmultihap"

# Local variables
__version__ = "0.1.1"
__version_date__ = "(16-Oct-2019)"
#
# These lines (or something similar) will be needed in the HAP processing code
#
MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)

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

    log_level : string
        The desired level of verboseness in the log statements displayed on the screen and written to the
        .log file.

    Updates
    -------
    return_value : list
        a simple status value. '0' for a successful run and '1' for a failed
        run
    """
    # set up log_level as an input to hapmultisequencer.run_mvm_processing().
    log_level_dict = {"critical": logutil.logging.CRITICAL,
                      "error": logutil.logging.ERROR,
                      "warning": logutil.logging.WARNING,
                      "info": logutil.logging.INFO,
                      "debug": logutil.logging.DEBUG}

    if 'log_level' in kwargs:
        kwargs['log_level'] = kwargs['log_level'].lower()
        if kwargs['log_level'] in log_level_dict.keys():
            kwargs['log_level'] = log_level_dict[kwargs['log_level']]
        else:
            print("Log level set to default level 'log.info'.")
            kwargs['log_level'] = logutil.logging.INFO
    else:
        print("Log level set to default level 'log.info'.")
        kwargs['log_level'] = logutil.logging.INFO

    # execute hapmultisequencer.run_mvm_processing()
    return_value = hapmultisequencer.run_mvm_processing(input_filename, **kwargs)

    return return_value

# ------------------------------------------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description='Process images, produce mosaics and sourcelists')
    parser.add_argument('input_filename',
                        help='Name of the input csv file containing information about the files to be '
                             'processed')
    parser.add_argument('-d', '--diagnostic_mode', required=False, action='store_true',
                        help='If this option is turned on, additional log messages will be displayed and '
                             'additional files will be created during the course of the run.')
    parser.add_argument('-l', '--log_level', required=False, default='info',
                        choices=['critical', 'error', 'warning', 'info', 'debug'],
                        help='The desired level of verboseness in the log statements displayed on the screen '
                             'and written to the .log file. The level of verboseness from left to right, and '
                             'includes all log statements with a log_level left of the specified level. '
                             'Specifying "critical" will only record/display "critical" log statements, and '
                             'specifying "error" will record/display both "error" and "critical" log '
                             'statements, and so on.')
    parser.add_argument('-s', '--skip_gaia_alignment', required=False, action='store_true',
                        help='Skip alignment of all input images to known Gaia/HSC sources in the input '
                             'image footprint? If this option is turned on, the existing input image '
                             'alignment solution will be used instead. The default is False.')
    user_args = parser.parse_args()

    print("Multi-visit processing started for: {}".format(user_args.input_filename))
    rv = perform(user_args.input_filename,
                 diagnostic_mode=user_args.diagnostic_mode,
                 log_level=user_args.log_level,
                 skip_gaia_alignment=user_args.skip_gaia_alignment)
    print("Return Value: ", rv)
    return rv


# ------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    main()
