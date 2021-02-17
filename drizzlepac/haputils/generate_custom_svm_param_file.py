#!/usr/bin/env python

""" generate_custom_svm_param_file - a module to create a template SVM processing pipeline parameter file
based on the observations present in the current working directory for the user to customize

Command-line USAGE:
    >>> python generate_custom_svm_param_file.py [-clo] poller_filename

    - poller_filename: (required) the name of the input SVM poller file

    - -c: (optional) If turned on, existing files with the same name as the output custom SVM parameter file
      created by this script will be overwritten.

    - -l: (optional) The desired level of verboseness in the log statements displayed on the screen and
      written to the .log file. Valid input options are 'critical', 'error', 'warning', 'info', or 'debug'.

    - -o: (optional) Name of the output configuration JSON file which will be created for specialized
      processing. This file will contain ALL the input parameters necessary for processing. If not explicitly
      specified, the default name for the output file is "custom_svm_params.json".

Python USAGE:
    >>> python
    >>> from drizzlepac.haputils import generate_custom_svm_param_file
    >>> generate_custom_svm_param_file.make_svm_input_file(input_filename,
                                                           output_custom_pars_file='custom_svm_params.json',
                                                           clobber=False,
                                                           log_level=logutil.logging.INFO)

.. note:: This script only generates a template input parameter file populated with default values based on
          the observations present in the current working directory. It is entirely incumbent on the user to
          modify this file with non-default values.
"""

import argparse
import datetime
import json
import logging
import os
import sys
import traceback

from drizzlepac.haputils import config_utils
from drizzlepac.haputils import ci_table
from drizzlepac.haputils import poller_utils

from stsci.tools import logutil

__taskname__ = 'generate_custom_svm_param_file'
MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)
__version__ = 0.1
__version_date__ = '04-Dec-2020'
# ----------------------------------------------------------------------------------------------------------------------


def make_svm_input_file(input_filename, output_custom_pars_file='custom_svm_params.json', clobber=False,
                        log_level=logutil.logging.INFO):
    """
    create a custom SVM processing pipeline parameter file based on the observations present in the current
    working directory using config_utils.HapConfig() and optionally update_ci_values() to adjust CI upper and
    lower limits for filter products

    Parameters
    ----------
    input_filename: str
        The 'poller file' where each line contains information regarding an exposures taken
        during a single visit.

    output_custom_pars_file: str, optional
        Fully specified output filename which contains all the configuration parameters
        available during the processing session. Default is 'custom_svm_params.json'.

    clobber : Bool, optional
        If set to Boolean 'True', existing files with the same name as *output_custom_pars_file*, the output
        custom SVM parameter file created by this script will be overwritten. Default value is Boolean
        'False'.

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the
        .log file. Default value is 20, or 'info'.

    RETURNS
    -------
    Nothing.
    """
    log.setLevel(log_level)
    if not clobber:
        if os.path.exists(output_custom_pars_file):
            msg = "A file named '{}' already exists. Please choose a unique name for the custom SVM parameter file.".format(output_custom_pars_file)
            log.critical(msg)
            sys.exit()
    # Define trailer file (log file) that will contain the log entries for all processing
    if isinstance(input_filename, str):  # input file is a poller file -- easy case
        logname = input_filename.replace('.out', '_svm_partam_gen.log')

    else:
        logname = 'svm_param_gen.log'

    # Initialize total trailer filename as temp logname
    logging.basicConfig(filename=logname, format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)
    # start processing
    starting_dt = datetime.datetime.now()
    log.info("Run start time: {}".format(str(starting_dt)))

    try:
        # Parse the poller file and generate the the obs_info_dict, as well as the total detection
        # product lists which contain the ExposureProduct, FilterProduct, and TotalProduct objects
        # A poller file contains visit data for a single instrument.  The TotalProduct discriminant
        # is the detector.  A TotalProduct object is comprised of FilterProducts and ExposureProducts
        # where its FilterProduct is distinguished by the filter in use, and the ExposureProduct
        # is the atomic exposure data.
        log.info("Parse the poller and determine what exposures need to be combined into separate products.\n")
        obs_info_dict, total_obj_list = poller_utils.interpret_obset_input(input_filename, log_level)

        # Update all of the product objects with their associated configuration information.
        for total_item in total_obj_list:
            log.info("Preparing configuration parameter values for total product {}".format(total_item.drizzle_filename))
            total_item.configobj_pars = config_utils.HapConfig(total_item,
                                                               log_level=log_level,
                                                               output_custom_pars_file=output_custom_pars_file)
            for filter_item in total_item.fdp_list:
                log.info("Preparing configuration parameter values for filter product {}".format(filter_item.drizzle_filename))
                filter_item.configobj_pars = config_utils.HapConfig(filter_item,
                                                                    log_level=log_level,
                                                                    output_custom_pars_file=output_custom_pars_file)
                update_ci_values(filter_item, output_custom_pars_file, log_level)

            for expo_item in total_item.edp_list:
                log.info("Preparing configuration parameter values for exposure product {}".format(expo_item.drizzle_filename))
                expo_item.configobj_pars = config_utils.HapConfig(expo_item,
                                                                  log_level=log_level,
                                                                  output_custom_pars_file=output_custom_pars_file)
                # Housekeeping: remove those pesky renamed copies of the input flc.fits/flt.fits files
                # generated by drizzlepac.haputils.product()
                if expo_item.drizzle_filename.endswith("_drc.fits"):
                    file_to_remove = expo_item.drizzle_filename.replace("_drc.fits", "_flc.fits")
                if expo_item.drizzle_filename.endswith("_drz.fits"):
                    file_to_remove = expo_item.drizzle_filename.replace("_drz.fits", "_flt.fits")
                if os.path.exists(file_to_remove):
                    os.remove(file_to_remove)
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)
        err_msg = "Something went wrong!"
        log.error(err_msg)
        raise Exception(err_msg)


# ----------------------------------------------------------------------------------------------------------------------


def update_ci_values(filter_item, config_filename, log_level=logutil.logging.INFO):
    """
    if the 'lookup_ci_limits_from_table' parameter is set to Boolean 'True', update the custom parameter
    file with concentration index (CI) upper and lower limit values from ci_table.get_ci_from_file() for the
    specified filter object.

    Parameters
    ----------
    filter_item : drizzlepac.haputils.product.FilterProduct object
        object containing all the relevant info for the drizzled filter product

    config_filename : str
        Name of the output custom SVM parameter file

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the
        .log file. Default value is 20, or 'info'.

    Returns
    -------
    Nothing!
    """
    log.setLevel(log_level)
    # Check to see if the CI values need to be updated. If not, simply skip over all of this.
    update_ci_vals = False
    for phot_mode in ['aperture', 'segment']:
        if filter_item.configobj_pars.pars['quality control'].outpars['ci filter'][phot_mode]['lookup_ci_limits_from_table']:
            update_ci_vals = True
    if update_ci_vals:
        # read in the custom SVM param file data
        with open(config_filename) as f:
            json_data = json.load(f)
        for phot_mode in ['aperture', 'segment']:
            if filter_item.configobj_pars.pars['quality control'].outpars['ci filter'][phot_mode]['lookup_ci_limits_from_table']:
                log.info("NOTE: The 'lookup_ci_limits_from_table' setting in the 'quality control'>'{}' section of the parameters for filter image {} is set to 'True'. This means that any custom user-tuned values for 'ci_upper_limit' and 'ci_lower_limit' will be overwritten. To prevent this, please set 'lookup_ci_limits_from_table' to 'False' in the custom parameter file {}".format(phot_mode, filter_item.drizzle_filename, config_filename))
                # set up inputs to ci_table.get_ci_from_file() and execute to get new CI values
                drizzled_image = filter_item.drizzle_filename
                ci_lookup_file_path = "default_parameters/any"
                diagnostic_mode = False
                ci_lower_limit = filter_item.configobj_pars.pars['quality control'].outpars['ci filter'][phot_mode]['ci_lower_limit']
                ci_upper_limit = filter_item.configobj_pars.pars['quality control'].outpars['ci filter'][phot_mode]['ci_upper_limit']
                ci_dict = ci_table.get_ci_from_file(drizzled_image, ci_lookup_file_path, log_level,
                                                    diagnostic_mode=diagnostic_mode, ci_lower=ci_lower_limit,
                                                    ci_upper=ci_upper_limit)
                log.debug("{} {} CI upper limit updated from {} to {}".format(filter_item.drizzle_filename,
                                                                              phot_mode,
                                                                              ci_upper_limit,
                                                                              ci_dict["ci_upper_limit"]))
                log.debug("{} {} CI lower limit updated from {} to {}\n".format(filter_item.drizzle_filename,
                                                                                phot_mode,
                                                                                ci_lower_limit,
                                                                                ci_dict["ci_lower_limit"]))
                # update CI values
                json_data[drizzled_image[:-9]]["default_values"]["quality control"]["ci filter"][phot_mode]["ci_lower_limit"] = ci_dict["ci_lower_limit"]
                json_data[drizzled_image[:-9]]["default_values"]["quality control"]["ci filter"][phot_mode]["ci_upper_limit"] = ci_dict["ci_upper_limit"]

        # Write out the updated custom SVM param file data back to file.
        with open(config_filename, 'w') as f:
            json.dump(json_data, f, indent=4)
    else:
        log.debug("Using existing concentration index limits from parameter file")


# ----------------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process images, produce drizzled images and sourcelists')
    parser.add_argument('poller_filename', help='Name of the input csv file containing information about the '
                                                'files to be processed')
    parser.add_argument('-c', '--clobber', required=False, action='store_true',
                        help='If turned on, existing files with the same name as the output custom SVM '
                             'parameter file created by this script will be overwritten.')
    parser.add_argument('-l', '--log_level', required=False, default='info',
                        choices=['critical', 'error', 'warning', 'info', 'debug'], help='The desired level '
                        'of verboseness in the log statements displayed on the screen and written to the '
                        '.log file. The level of verboseness from left to right, and includes all log '
                        'statements with a log_level left of the specified level. Specifying "critical" will '
                        'only record/display "critical" log statements, and specifying "error" will '
                        'record/display both "error" and "critical" log statements, and so on.')
    parser.add_argument('-o', '--output_custom_pars_file', required=False, default="custom_svm_params.json",
                        help='Name of the output configuration JSON file which will be created for '
                             'specialized processing. This file will contain ALL the input parameters '
                             'necessary for processing. If not explicitly specified, the default name for '
                             'the output file is "custom_svm_params.json".')
    user_args = parser.parse_args()

    log_level_dict = {"critical": logutil.logging.CRITICAL,
                      "error": logutil.logging.ERROR,
                      "warning": logutil.logging.WARNING,
                      "info": logutil.logging.INFO,
                      "debug": logutil.logging.DEBUG}

    make_svm_input_file(user_args.poller_filename, output_custom_pars_file=user_args.output_custom_pars_file,
                        clobber=user_args.clobber, log_level=log_level_dict[user_args.log_level])
