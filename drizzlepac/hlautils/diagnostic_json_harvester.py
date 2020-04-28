#!/usr/bin/env python

"""This script 'harvests' information stored in the .json files produced by
drizzlepac/hlautils/svm_quality_analysis.py and stores it as a Pandas DataFrame"""

# Standard library imports
import glob
import json
import os
import pdb
import sys

# Related third party imports
import pandas as pd

# Local application imports
import drizzlepac.hlautils.diagnostic_utils as du
from stsci.tools import logutil

__taskname__ = 'diagnostic_json_harvester'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)


# ------------------------------------------------------------------------------------------------------------


def get_json_files(search_path="", log_level=logutil.logging.INFO):
    """use grep to create a list of json files to harvest

    Parameters
    ----------
    search_path : str, optional
        directory path to search for .json files. Default value is the current working directory.

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the
        .log file. Default value is 'INFO'.

    Returns
    -------
    sorted_json_list : list
        list of json files to harvest, sorted first alphabetically by filename
    """

    # set up search string and use glob to get list of files
    log.setLevel(log_level)
    if search_path != "" and not search_path.endswith("/"):
        search_path += "/"
    search_string = search_path + "*_svm_*.json"
    json_list = glob.glob(search_string)

    # Fail gracefully if no .json files were found
    if not json_list:
        err_msg = "No .json files were found!"
        log.error(err_msg)
        raise Exception(err_msg)

    # sort list first alphabetically by filename
    sorted_json_list = sorted(json_list)

    return sorted_json_list


# ------------------------------------------------------------------------------------------------------------


def json_harvester(log_level=logutil.logging.INFO):
    """Main calling function

    Parameters
    ----------
    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the
        .log file. Default value is 'INFO'.

    Returns
    -------
    TBD
    """
    log.setLevel(log_level)
    json_list = get_json_files(log_level=log_level)
    # json_list = [json_list[1]]
    for json_filename in json_list:
        json_data = du.read_json_file(json_filename)
        print(json_filename)
        print(json_data.keys())
        for item in json_data['data'].keys():
            print(">>", item)
            if hasattr(json_data['data'][item],"keys"):
                for item2 in json_data['data'][item].keys():
                    print(">>>>>{}: {}".format(item2,json_data['data'][item][item2]))
            else:
                print(">> {}: {}".format(item,json_data['data'][item]))
        input("\n")



# ======================================================================================================================


if __name__ == "__main__":
    # Testing
    json_harvester(log_level=logutil.logging.DEBUG)
