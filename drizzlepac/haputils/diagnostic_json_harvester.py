#!/usr/bin/env python

"""This script 'harvests' information stored in the .json files produced by
drizzlepac/haputils/svm_quality_analysis.py and stores it as a Pandas DataFrame"""

# Standard library imports
import argparse
import collections
import glob
import os
import pdb
import sys

# Related third party imports
import pandas as pd

# Local application imports
import drizzlepac.haputils.diagnostic_utils as du
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


def get_json_files(search_path=os.getcwd(), 
                   search_patterns = ["*_svm_*.json", "*_mvm_*.json", "*_cal_qa_*.json"],
                   log_level=logutil.logging.INFO):
    """use glob to create a list of json files to harvest
    
    This function looks for all the json files containing qa test results generated
    by `runastrodriz` and `runsinglehap`.  The search starts in the directory 
    specified in the `search_path` parameter, but will look in immediate
    sub-directories as well if no json files are located in the directory 
    specified by `search_path`.

    Parameters
    ----------
    search_path : str, optional
        directory path to search for .json files. Default value is the current working directory.
        This serves as the starting directory for finding the .json files, with the
        search expanding to immediate sub-directories if no .json files are found
        in this directory.

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the
        .log file. Default value is 'INFO'.

    Returns
    -------
    out_json_dict : ordered dictionary
        dictionary containing lists of all identified json files, grouped by and  keyed by Pandas DataFrame
        index value
    """
    log.setLevel(log_level)

    # set up search string and use glob to get list of files
    json_list = []
    for search_pattern in search_patterns:
        search_string = os.path.join(search_path, search_pattern)
        search_results = glob.glob(search_string)
        if len(search_results) == 0:
            # Try another directory lower
            search_string = os.path.join(search_path, '*', search_pattern)
            search_results = glob.glob(search_string)
            
        log.info("{} files found: {}".format(search_pattern, len(search_results)))
        if len(search_results) > 0:
            json_list += search_results

    # store json filenames in a dictionary keyed by Pandas DataFrame index value
    if json_list:
        out_json_dict = collections.OrderedDict()
        for json_filename in sorted(json_list):
            json_data = du.read_json_file(json_filename)
            dataframe_idx = json_data['general information']['dataframe_index']
            if dataframe_idx in out_json_dict.keys():
                out_json_dict[dataframe_idx].append(json_filename)
            else:
                out_json_dict[dataframe_idx] = [json_filename]
            del(json_data)  # Housekeeping!

    # Fail gracefully if no .json files were found
    else:
        err_msg = "No .json files were found!"
        log.error(err_msg)
        raise Exception(err_msg)
    return out_json_dict

# ------------------------------------------------------------------------------------------------------------


def h5load(filename):
    """Load pandas DataFrame and metadata stored in specified HDF5 file. If filename does not end with a '.h5'
    extension,  an exception will be raised.
    Code borrowed from https://stackoverflow.com/questions/29129095/save-additional-attributes-in-pandas-dataframe/29130146#29130146

    Parameters
    ----------
    filename : str
        Name of the HDF5 file to open

    Returns
    -------
    data : Pandas DataFrame
        Pandas Dataframe read in from specified HDF5 file

    metadata : dict
        dictionary containing dataFrame metadata
    """
    if os.path.exists(filename):
        if not filename.endswith(".h5"):
            errmsg = "Incorrect input file extension. The expected input filetype is HDF5, and should have "
            errmsg += "the extension '.h5', not '{}'.".format(os.path.splitext(filename)[1])
            log.error(errmsg)
            raise Exception(errmsg)
        with pd.HDFStore(filename) as store:
            data = store['mydata']
            metadata = store.get_storer('mydata').attrs.metadata
    else:
        errmsg = "HDF5 file {} not found!".format(filename)
        log.error(errmsg)
        raise Exception(errmsg)
    return data, metadata


# ------------------------------------------------------------------------------------------------------------

def h5store(filename, df, **kwargs):
    """Write a pandas Dataframe and metadata to a HDF5 file.
    Code borrowed from https://stackoverflow.com/questions/29129095/save-additional-attributes-in-pandas-dataframe/29130146#29130146

    Parameters
    ----------
    filename : str
        Name of the output HDF5 file

    df : Pandas DataFrame
        Pandas DataFrame to write to output file

    Returns
    -------
    Nothing.
    """
    store = pd.HDFStore(filename)
    store.put('mydata', df)
    store.get_storer('mydata').attrs.metadata = kwargs
    store.close()


# ------------------------------------------------------------------------------------------------------------

def json_harvester(json_search_path=os.getcwd(), 
                   json_patterns = ["*_svm_*.json", "*_mvm_*.json", "*_cal_qa_*.json"],
                   log_level=logutil.logging.INFO,
                   output_filename_basename="svm_qa_dataframe"):
    """Main calling function

    Parameters
    ----------
    json_search_path : str, optional
        The full path of the directory that will be searched for json files to process. If not explicitly
        specified, the current working directory will be used.

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the
        .log file. Default value is 'INFO'.

    output_filename_basename : str, optional
        Name of the output file basename (filename without the extension) for the Hierarchical Data Format
        version 5 (HDF5) .h5 file that the Pandas DataFrame will be written to. If not explicitly specified,
        the default filename basename that will be used is "svm_qa_dataframe". The default location that the
        output file will be written to is the current working directory
    """
    log.setLevel(log_level)

    # Get sorted list of json files
    json_dict = get_json_files(search_path=json_search_path,
                               search_patterns=json_patterns,
                               log_level=log_level)
    master_dataframe = None
    # extract all information from all json files related to a specific Pandas DataFrame index value into a
    # single line in the master dataframe
    num_json = len(json_dict)
    for n,idx in enumerate(json_dict.keys()):
        if ((n/num_json) % 0.1) == 0:
            log.info("Harvested {}% of the JSON files".format((n/num_json)*100))

        ingest_dict = make_dataframe_line(json_dict[idx], log_level=log_level)
        if ingest_dict:
            if master_dataframe is not None:
                log.debug("APPENDED DATAFRAME")
                master_dataframe = master_dataframe.append(pd.DataFrame(ingest_dict["data"], index=[idx]))
                for key_item in ingest_dict['descriptions'].keys():
                    if key_item not in metadata['descriptions'].keys(): # only have to check descriptions dict because units dict has same keys
                        metadata['descriptions'][key_item] = ingest_dict['descriptions'][key_item]
                        metadata['units'][key_item] = ingest_dict['units'][key_item]
            else:
                log.debug("CREATED DATAFRAME")
                master_dataframe = pd.DataFrame(ingest_dict["data"], index=[idx])
                metadata = {}
                metadata['descriptions'] = ingest_dict['descriptions']
                metadata['units'] = ingest_dict['units']

    # Write master_dataframe out to a HDF5 .hdfile
    if master_dataframe is not None:
        output_h5_filename = output_filename_basename + ".h5"
        if os.path.exists(output_h5_filename):
            os.remove(output_h5_filename)
        h5store(output_h5_filename, master_dataframe, **metadata)
        log.info("Wrote dataframe and metadata to HDF5 file {}".format(output_h5_filename))
        #optionally also write dataframe out to .csv file for human inspection
        output_csv_filename = output_filename_basename + ".csv"
        if log_level == logutil.logging.DEBUG:
            if os.path.exists(output_csv_filename):
                os.remove(output_csv_filename)
            master_dataframe.to_csv(output_csv_filename)
            log.debug("Wrote dataframe to csv file {}".format(output_csv_filename))

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
    log.setLevel(log_level)
    header_ingested = False
    gen_info_ingested = False
    ingest_dict = collections.OrderedDict()
    ingest_dict['data'] = collections.OrderedDict()
    ingest_dict['descriptions'] = collections.OrderedDict()
    ingest_dict['units'] = collections.OrderedDict()
    for json_filename in json_filename_list:
        short_json_filename = os.path.basename(json_filename)
        # This is to differentiate point catalog compare_sourcelists columns from segment catalog
        # compare_sourcelists columns in the dataframe
        if json_filename.endswith("_point-cat_svm_compare_sourcelists.json"):
            title_suffix = "hap_vs_hla_point_"
        elif json_filename.endswith("_segment-cat_svm_compare_sourcelists.json"):
            title_suffix = "hap_vs_hla_segment_"
        else:
            title_suffix = ""
        json_data = du.read_json_file(json_filename)
        # add information from "header" section to ingest_dict just once
        if not header_ingested:
            # filter out ALL header keywords not included in 'header_keywords_to_keep'
            header_keywords_to_keep = ['APERTURE',
                                       'CHINJECT',
                                       'DATE-OBS',
                                       'DEC_TARG',
                                       'EXPTIME',
                                       'FGSLOCK',
                                       'GYROMODE',
                                       'MTFLAG',
                                       'OBSKEY',
                                       'OBSTYPE',
                                       'RA_TARG',
                                       'SCAN_TYP',
                                       'SUBARRAY',
                                       'TIME-OBS']
            for header_item in json_data['header'].keys():
                if header_item in header_keywords_to_keep:
                    ingest_dict["data"]["header."+header_item] = json_data['header'][header_item]
            header_ingested = True

        # add information from "general information" section to ingest_dict just once
        if not gen_info_ingested:
            for gi_item in json_data['general information'].keys():
                ingest_dict["data"]["gen_info."+gi_item] = json_data['general information'][gi_item]
            gen_info_ingested = True

        # recursively flatten nested "data" section dictionaries and build ingest_dict
        flattened_data = flatten_dict(json_data['data'])
        flattened_descriptions = flatten_dict(json_data['descriptions'])
        flattened_units = flatten_dict(json_data['units'])
        for fd_key in flattened_data.keys():
            json_data_item = flattened_data[fd_key]
            ingest_key = fd_key.replace(" ", "_")
            if str(type(json_data_item)) == "<class 'astropy.table.table.Table'>":
                for coltitle in json_data_item.colnames:
                    ingest_value = json_data_item[coltitle].tolist()
                    id_key = title_suffix + ingest_key + "." + coltitle
                    ingest_dict["data"][id_key] = [ingest_value]
                    # print(">>>>",short_json_filename,id_key, title_suffix + fd_key + "." + coltitle)
                    try:
                        ingest_dict["descriptions"][id_key] = flattened_descriptions[fd_key + "." + coltitle]
                        log.debug("Added Description {} = {}".format(id_key, flattened_descriptions[fd_key + "." + coltitle]))
                    except:  # insert placeholders if the code runs into trouble getting descriptions
                        log.warning("Descriptions not found for {} in file{}. Using placeholder value '>>>UNDEFINED<<<' instead.".format(id_key, short_json_filename))
                        ingest_dict["descriptions"][id_key] = ">>>UNDEFINED<<<"

                    try:
                        ingest_dict["units"][id_key] = flattened_units[fd_key + "." + coltitle]
                        log.debug("Added Unit {} = {}".format(id_key, flattened_units[fd_key + "." + coltitle]))
                    except:  # insert placeholders if the code runs into trouble getting units
                        log.warning("Units not found for {} in file {}. Using placeholder value '>>>UNDEFINED<<<' instead.".format(id_key, short_json_filename))
                        ingest_dict["units"][id_key] = ">>>UNDEFINED<<<"

            else:
                ingest_value = json_data_item
                id_key = title_suffix + ingest_key
                if str(type(ingest_value)) == "<class 'list'>":
                    ingest_dict["data"][id_key] = [ingest_value]
                else:
                    ingest_dict["data"][id_key] = ingest_value
                try:
                    ingest_dict["descriptions"][id_key] = flattened_descriptions[fd_key]
                    log.debug("Added Description {} = {}".format(id_key, flattened_descriptions[fd_key]))
                except:  # insert placeholders if the code runs into trouble getting the descriptions
                    log.warning("Descriptions not found for {} in file {}. Using placeholder value '>>>UNDEFINED<<<' instead.".format(id_key, short_json_filename))
                    ingest_dict["descriptions"][id_key] = ">>>UNDEFINED<<<"

                try:
                    ingest_dict["units"][id_key] = flattened_units[fd_key]
                    log.debug("Added unit {} = {}".format(id_key, flattened_units[fd_key]))
                except:  # insert placeholders if the code runs into trouble getting units
                    log.warning("Units not found for {} in file {}. Using placeholder value '>>>UNDEFINED<<<' instead.".format(id_key, short_json_filename))
                    ingest_dict["units"][id_key] = ">>>UNDEFINED<<<"

    return ingest_dict

# ======================================================================================================================


if __name__ == "__main__":
    # process command-line inputs with argparse
    parser = argparse.ArgumentParser(description='ingest all SVM QA json files into a Pandas DataFrame and'
                                                 'and write DataFrame to an .csv file.')
    parser.add_argument('-j', '--json_search_path', required=False, default=os.getcwd(),
                        help='The full path of the directory that will be searched for json files to '
                             'process. If not explicitly specified, the current working directory will be '
                             'used.')
    parser.add_argument('-l', '--log_level', required=False, default='info',
                        choices=['critical', 'error', 'warning', 'info', 'debug'],
                        help='The desired level of verboseness in the log statements displayed on the screen '
                             'and written to the .log file. The level of verboseness from left to right, and '
                             'includes all log statements with a log_level left of the specified level. '
                             'Specifying "critical" will only record/display "critical" log statements, and '
                             'specifying "error" will record/display both "error" and "critical" log '
                             'statements, and so on. If the log_level is set to "debug", a human-readable'
                             '.csv will be generated along with the .h5 file.')
    parser.add_argument('-o', '--output_filename_basename', required=False, default="svm_qa_dataframe",
                        help='Name of the output file basename (filename without the extension) for the  '
                             'Hierarchical Data Format version 5 (HDF5) .h5 file that the Pandas DataFrame '
                             'will be written to. If not explicitly specified, the default filename basename '
                             'that will be used is "svm_qa_dataframe". The default location that the output '
                             'file will be written to is the current working directory')
    user_args = parser.parse_args()

    # set up logging
    log_dict = {"critical": logutil.logging.CRITICAL,
                "error": logutil.logging.ERROR,
                "warning": logutil.logging.WARNING,
                "info": logutil.logging.INFO,
                "debug": logutil.logging.DEBUG}
    log_level = log_dict[user_args.log_level]
    log.setLevel(log_level)

    json_harvester(json_search_path=user_args.json_search_path,
                   log_level=log_level,
                   output_filename_basename=user_args.output_filename_basename)
