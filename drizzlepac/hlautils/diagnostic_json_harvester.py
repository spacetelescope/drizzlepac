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


def flatten_dict(dd, separator='.', prefix=''):
    """Recursive subroutine to flatten nested dictionaries down into a single-layer dictionary.
    Borrowed from https://www.geeksforgeeks.org/python-convert-nested-dictionary-into-flattened-dictionary/

    Parameters
    ----------
    dd : dict
        dictionary to flatten

    separator : str, optional
        separator character used in constructing flattened dictionary key names from multiple recursive
        elements. Default value is '.'

    prefix : str, optional
        flattened dictionary key prefix. Default value is an empty string ('').

    Returns
    -------
    a version of input dictionary *dd* that has been flattened by one layer
    """
    return {prefix + separator + k if prefix else k: v
            for kk, vv in dd.items()
            for k, v in flatten_dict(vv, separator, kk).items()
            } if isinstance(dd, dict) else {prefix: dd}

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
    master_dataframe : Pandas DataFrame
        pandas DataFrame containing all information harvested from the json files.
    """
    log.setLevel(log_level)

    # Get sorted list of json files
    json_dict = get_json_files(log_level=log_level)
    master_dataframe = None
    for idx in json_dict.keys():
        ingest_dict = make_dataframe_line(json_dict[idx], log_level=log_level)
        if ingest_dict:
            if master_dataframe is not None:
                log.debug("APPENDED DATAFRAME")
                master_dataframe = master_dataframe.append(pd.DataFrame(ingest_dict, index=[idx]))
            else:
                log.debug("CREATED DATAFRAME")
                master_dataframe = pd.DataFrame(ingest_dict, index=[idx])
    if master_dataframe is not None:
        out_csv_filename = "master_dataframe.csv"
        if os.path.exists(out_csv_filename):
            os.remove(out_csv_filename)

        master_dataframe.to_csv(out_csv_filename)
        print("Wrote "+out_csv_filename)

    return master_dataframe


# ------------------------------------------------------------------------------------------------------------


def make_dataframe_line(json_filename_list, log_level=logutil.logging.INFO):
    """extracts information from the json files specified by the input list *json_filename_list*.

    Parameters
    ----------
    json_filename_list : list
        list of json files to process

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the
        .log file. Default value is 'INFO'.

    Returns
    -------
    ingest_dict : collections.OrderedDict
        ordered dictionary containing all information extracted from json files specified by the input list
        *json_filename_list*.
    """
    header_ingested = False
    gen_info_ingested = False
    ingest_dict = collections.OrderedDict()
    for json_filename in json_filename_list:
        if json_filename.endswith("_point-cat_svm_compare_sourcelists.json"):
            title_suffex = "hap_vs_hla_point_"
        elif json_filename.endswith("_segment-cat_svm_compare_sourcelists.json"):
            title_suffex = "hap_vs_hla_segment_"
        else:
            title_suffex = ""
        json_data = du.read_json_file(json_filename)
        if not header_ingested:
            for header_item in json_data['header'].keys():
                ingest_dict["header."+header_item] = json_data['header'][header_item]
            header_ingested = True
        if not gen_info_ingested:
            for gi_item in json_data['general information'].keys():
                ingest_dict["gen_info."+gi_item] = json_data['general information'][gi_item]
            gen_info_ingested = True
        flattened_data = flatten_dict(json_data['data'])
        for fd_key in flattened_data.keys():
            json_data_item = flattened_data[fd_key]
            ingest_key = fd_key.replace(" ","_")
            if str(type(json_data_item)) == "<class 'astropy.table.table.Table'>":
                for coltitle in json_data_item.colnames:
                    ingest_value = json_data_item[coltitle].tolist()
                    ingest_dict[title_suffex + ingest_key + "." + coltitle] = [ingest_value]
            else:
                ingest_value = json_data_item
                ingest_dict[title_suffex + ingest_key] = ingest_value
    return ingest_dict



# ======================================================================================================================


if __name__ == "__main__":
    #  Testing
    json_harvester(log_level=logutil.logging.DEBUG)
