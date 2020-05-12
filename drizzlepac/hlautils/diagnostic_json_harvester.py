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
    df_line = make_dataframe_line(json_dict, log_level=log_level)

    # master_dataframe = None
    #
    # for json_filename in json_list:
    #     master_dataframe = json_ingest(master_dataframe, json_filename, log_level=log_level)
    # if len(master_dataframe) > 0:
    #     master_dataframe = pd.concat(master_dataframe)
    # if master_dataframe is not None:
    #     out_csv_filename = "master_dataframe.csv"
    #     if os.path.exists(out_csv_filename):
    #         os.remove(out_csv_filename)
    #
    #     master_dataframe.to_csv(out_csv_filename)
    #     print("Wrote "+out_csv_filename)


# ------------------------------------------------------------------------------------------------------------


def make_dataframe_line(json_dict, log_level=logutil.logging.INFO):
    allowed_json_types = ['_svm_gaia_distribution_characterization.json',
                          '_svm_num_sources.json',
                          '_photometry.json',
                          '_svm_wcs.json']
    allowed_json_types=['_photometry.json']
    for idx in json_dict.keys():
        dataframe_line = None
        header_ingested = True # TODO: RESET TO FALSE
        gen_info_ingested = False
        ingest_dict = collections.OrderedDict()

        print(idx)  # TODO: REMOVE
        for json_filename in json_dict[idx]:

            file_status = False
            for ending in allowed_json_types:
                if json_filename.endswith(ending):
                    file_status = True
                    break
            if file_status: # TODO: REMOVE ONCE THIGNS are working and de-indent everything by one level.
                json_data = du.read_json_file(json_filename)
                if not header_ingested:
                    for header_item in json_data['header'].keys():
                        ingest_dict["header-"+header_item] = json_data['header'][header_item]
                    header_ingested = True
                if not gen_info_ingested:
                    dataframe_index = json_data['general information']["dataframe_index"]
                    for gi_item in json_data['general information'].keys():
                        ingest_dict["gen_info-"+gi_item] = json_data['general information'][gi_item]
                    gen_info_ingested = True
                print(json_filename)

                foo = flatten_dict(json_data['data'], '', {})
                for item in foo.keys():
                    print(item,foo[item])






            print("   {} {}".format(file_status,json_filename))  # TODO: REMOVE

    return 0
# ------------------------------------------------------------------------------------------------------------


def json_ingest(master_dataframe, json_filename, log_level=logutil.logging.INFO):
    """ingests data from specified json file into a pandas dataframe

    Parameters
    ----------
    master_dataframe : pandas DataFrame
        The pandas DataFrame that information from the specified json file will be appended to

    json_filename : str
        The json file to ingest into a master_datagrame

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the
        .log file. Default value is 'INFO'.

    Returns
    -------
    master_dataframe : pandas DataFrame
        an updated version the input master_dataframe that now includes information harvested from the json
        file specified in json_filename
    """
    # Generate pandas dataframe index string
    pdindex = json_filename.split("_svm_")[0]
    for pdindex_ending in ["point-cat", "segment-cat"]:
        if pdindex.endswith(pdindex_ending):
            pdindex = pdindex.replace("_" + pdindex_ending, "")
            break
    # NEXT TWO LINES ARE FOR TESTING. REMOVE!
    if not json_filename.endswith("crossmatch.json"):  # TODO: REMOVE!
        return master_dataframe  # TODO: REMOVE!

    print("-----------------------", json_filename, pdindex, "-----------------------")
    # ingest data from json file into the dataframe
    json_data = du.read_json_file(json_filename)
    json_header = json_data['header']
    json_data = json_data['data']
    ingest_dict = collections.OrderedDict()
    for data_item in json_data.keys():
        for diag_key in json_data[data_item].keys():
            new_key = "{}-{}".format(data_item.replace(" ","_"), diag_key.replace(" ","_"))
            print(new_key, json_data[data_item][diag_key])
            ingest_dict[new_key] = json_data[data_item][diag_key]

    if master_dataframe is not None:
        print("APPENDED DATAFRAME")
        master_dataframe.append(pd.DataFrame(ingest_dict, index=[pdindex]))

    else:
        print("CREATED DATAFRAME")
        master_dataframe = [pd.DataFrame(ingest_dict, index=[pdindex])]

    # print(json_data.keys())
    # for item in json_data['data'].keys():
    #     print(">>", item)
    #     if hasattr(json_data['data'][item],"keys"):
    #         for item2 in json_data['data'][item].keys():
    #             print(">>>>>{}: {}".format(item2,json_data['data'][item][item2]))
    #     else:
    #         print(">> {}: {}".format(item,json_data['data'][item]))
    # input("\n")
    return master_dataframe


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Mike's new stuff
def flatten_dict(current, key, result):
    """Flatten nested dictionaries into a non-nested dictionary. Assumes that there are no non-unique keys.
    Code credit: https://stackoverflow.com/questions/24448543/how-would-i-flatten-a-nested-dictionary-in-python-3
    Solution submitted by user 'Matthew Franglen'.
    """
    print(current,isinstance(current, dict))
    print(key)
    print(result)
    print("============================================\n")
    if isinstance(current, dict):
        for k in current:
            new_key = "{1}".format(key, k) if len(key) > 0 else k
            flatten_dict(current[k], new_key, result)
    else:
        result[key] = current
    input()
    return result

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# ======================================================================================================================


if __name__ == "__main__":
    #  Testing
    json_harvester(log_level=logutil.logging.DEBUG)

# TODO: automagically ingest data columns into single dataframe cell
# TODO: join individual
