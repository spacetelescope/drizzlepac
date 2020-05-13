#!/usr/bin/env python

"""This script 'harvests' information stored in the .json files produced by
drizzlepac/hlautils/svm_quality_analysis.py and stores it as a Pandas DataFrame"""

# Standard library imports
import collections
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
    out_json_dict : ordered dictionary
        dictionary containing lists of all identified json files, grouped by and  keyed by Pandas Dataframe
        index value
    """
    log.setLevel(log_level)

    # set up search string and use glob to get list of files
    search_string = os.path.join(search_path, "*_svm_*.json")
    json_list = glob.glob(search_string)

    # store json filenames in a dictionary keyed by Pandas DataFrame index value
    if json_list:
        out_json_dict = collections.OrderedDict()
        for json_filename in sorted(json_list):
            json_data = du.read_json_file(json_filename)
            dataframe_idx = json_data['general information']['dataframe_index']
            if dataframe_idx in out_json_dict.keys():
                out_json_dict[dataframe_idx].append(json_filename)
            else:
                out_json_dict[dataframe_idx]=[json_filename]
            del(json_data)  # Housekeeping!

    # Fail gracefully if no .json files were found
    else:
        err_msg = "No .json files were found!"
        log.error(err_msg)
        raise Exception(err_msg)

    return out_json_dict


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

    # Get sorted list of json files
    json_dict = get_json_files(log_level=log_level)
    master_dataframe = None
    for idx in json_dict.keys():
        ingest_dict = make_dataframe_line(json_dict[idx], idx, log_level=log_level)
        if ingest_dict:
            if master_dataframe is not None:
                print("APPENDED DATAFRAME")
                master_dataframe = master_dataframe.append(pd.DataFrame(ingest_dict, index=[idx]))

            else:
                print("CREATED DATAFRAME")
                master_dataframe = pd.DataFrame(ingest_dict, index=[idx])
    if master_dataframe is not None:
        out_csv_filename = "master_dataframe.csv"
        if os.path.exists(out_csv_filename):
            os.remove(out_csv_filename)

        master_dataframe.to_csv(out_csv_filename)
        print("Wrote "+out_csv_filename)


# ------------------------------------------------------------------------------------------------------------


def make_dataframe_line(json_filename_list, idx, log_level=logutil.logging.INFO):
    allowed_json_types = ['_svm_gaia_distribution_characterization.json',
                          '_svm_num_sources.json',
                          '_photometry.json',
                          '_svm_wcs.json']


    dataframe_line = None
    header_ingested = True # TODO: RESET TO FALSE
    gen_info_ingested = False
    ingest_dict = collections.OrderedDict()

    print(idx)  # TODO: REMOVE
    for json_filename in json_filename_list:
        file_status = False
        for ending in allowed_json_types:
            if json_filename.endswith(ending):
                file_status = True
                break
        if file_status: # TODO: REMOVE ONCE THIGNS are working and de-indent everything by one level.
            json_data = du.read_json_file(json_filename)
            if not header_ingested:
                for header_item in json_data['header'].keys():
                    ingest_dict["header."+header_item] = json_data['header'][header_item]
                header_ingested = True
            if not gen_info_ingested:
                for gi_item in json_data['general information'].keys():
                    ingest_dict["gen_info."+gi_item] = json_data['general information'][gi_item]
                gen_info_ingested = True
            print(json_filename)
            flattened_data = flatten_dict(json_data['data'])
            for fd_key in flattened_data.keys():
                ingest_dict[fd_key.replace(" ","_")] = flattened_data[fd_key]


        print("   {} {}".format(file_status,json_filename))  # TODO: REMOVE

    return ingest_dict
# ------------------------------------------------------------------------------------------------------------

def flatten_dict(dd, separator ='.', prefix =''):
    return { prefix + separator + k if prefix else k : v
             for kk, vv in dd.items()
             for k, v in flatten_dict(vv, separator, kk).items()
             } if isinstance(dd, dict) else { prefix : dd }


# ======================================================================================================================


if __name__ == "__main__":
    #  Testing
    json_harvester(log_level=logutil.logging.DEBUG)

# TODO: automagically ingest data columns into single dataframe cell
# TODO: join individual
