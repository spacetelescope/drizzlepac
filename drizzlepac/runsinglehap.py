#!/usr/bin/env python

""" runsinglehap.py - Module to control processing of single-visit mosaics

:License: :doc:`LICENSE`

USAGE:
    >>> runsinglehap [-cdl] inputFilename

    - The '-c' option allows the user to specify a customized configuration JSON file which has been tuned for
      specialized processing.  This file should contain ALL the input parameters necessary for processing. If
      not specified, default configuration values will be used.

    - The '-d' option will run this task in a more verbose diagnostic mode producing additional log messages
      will be displayed and additional files will be created.

    - The '-l' option allows the user to set the desired level of verboseness in the log statements displayed on
      the screen and written to the .log file. Specifying "critical" will only record/display "critical" log
      statements, and specifying "error" will record/display both "error" and "critical" log statements, and so
      on. Valid inputs: 'critical', 'error', 'warning', 'info', or 'debug'.

Python USAGE:
    >>> python
    >>> from drizzlepac import runsinglehap
    >>> runsinglehap.perform(inputFilename, debug=False, input_custom_pars_file=None, log_level='info')
"""
# Import standard Python modules
import argparse
import sys

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
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)

# ----------------------------------------------------------------------------------------------------------------------

def perform(input_filename, **kwargs):
    """Main calling subroutine for the ``runsinglehap`` task.

    Parameters
    ----------
    input_filename : str
        Name of the input csv file containing information about the files to
        be processed

    debug : bool, optional
        display all tracebacks, and debug information? If not specified, the default value is Boolean 'False'.

    input_custom_pars_file : str, optional
        Represents a fully specified input filename of a configuration JSON file which has been
        customized for specialized processing. This file should contain ALL the input parameters necessary
        for processing. If not specified, default configuration parameter values will be used.

    log_level : str, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the
        .log file. Valid inputs: 'critical', 'error', 'warning', 'info', or 'debug'. If not specified, the
        default value is 'info'.

    Updates
    -------
    return_value : list
        a simple status value. '0' for a successful run and '1' for a failed run
    """
    # set up log_level as an input to hapsequencer.run_hap_processing().
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

    # execute hapsequencer.run_hap_processing()
    return_value = hapsequencer.run_hap_processing(input_filename, **kwargs)

    return return_value

# ----------------------------------------------------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description='Process images, produce drizzled images and sourcelists')
    parser.add_argument('input_filename', help='Name of the input csv file containing information about the '
                        'files to be processed')
    parser.add_argument('-c', '--input_custom_pars_file', required=False, default=None, help='filename of a '
                        'configuration JSON file which has been customized for specialized processing. This '
                        'file should contain ALL the input parameters necessary for processing. If not '
                        'specified, default configuration parameter values will be used.')
    parser.add_argument('-d', '--diagnostic_mode', required=False, action='store_true', help='If this option '
                        'is turned on, additional log messages will be displayed and additional files will '
                        'be created during the course of the run.')
    parser.add_argument('-l', '--log_level', required=False, default='info',
                        choices=['critical', 'error', 'warning', 'info', 'debug'], help='The desired level '
                        'of verboseness in the log statements displayed on the screen and written to the '
                        '.log file. The level of verboseness from left to right, and includes all log '
                        'statements with a log_level left of the specified level. Specifying "critical" will '
                        'only record/display "critical" log statements, and specifying "error" will '
                        'record/display both "error" and "critical" log statements, and so on.')
    user_args = parser.parse_args()

    print("Single-visit processing started for: {}".format(user_args.input_filename))
    rv = perform(user_args.input_filename, input_custom_pars_file=user_args.input_custom_pars_file,
                 diagnostic_mode=user_args.diagnostic_mode, log_level=user_args.log_level)
    print("Return Value: ", rv)
    return rv

if __name__ == '__main__':
    main()
