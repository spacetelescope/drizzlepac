#!/usr/bin/env python

import argparse
import pdb

import pandas as pd

def make_search_string(exposure, skycell, config, spec):
    """Create ssearch string that will be used to query observation tables

    Parameters
    ----------
    exposure : str
        Image name

    skycell : str
        Skycell name

    config : str
        Instrument/Detector name. All caps, seperted by a "/" (e.g. WFC3/IR).

    spec : str
        Filter name. All caps.py

    Returns
    -------
    search_string : str
        Properly formatted seach string
    """
    search_string = ""

def query_dataframe(csv_data_file, search_string):
    dataframe = pd.DataFrame.from_csv(csv_data_file, header=0, index_col=0)
    results = dataframe.query(search_string)
    print(results)
    pdb.set_trace()

# ----------------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    # Argparse input
    parser = argparse.ArgumentParser(description='Search all observations')
    parser.add_argument('-c', '--csv_data_file', required=False, default="allexposures.csv",
                        help='Name of the input .csv file containing comma-separated columns index #, '
                             'exposure, skycell, config, and spec. Default value is "allexposures.csv".')
    parser.add_argument('-i', '--exposure', required=False, default="None", help='Image name.')
    parser.add_argument('-s', '--skycell', required=False, default="None", help='Skycell name.')
    parser.add_argument('-d', '--config', required=False, default="None", choices=["ACS/HRC", "ACS/SBC",
                                                                                   "ACS/WFC", "WFC3/IR",
                                                                                   "WFC3/UVIS", "None"],
                        help='Instrument/Detector name. All caps, seperted by a "/" (e.g. WFC3/IR).')

    parser.add_argument('-f', '--spec', required=False, default="None", help='Filter name.')
    in_args = parser.parse_args()
    pdb.set_trace()
    arg_dict = {}
    for item in in_args.__dict__.keys():
        if item not "csv_data_file":
            arg_dict[item] = in_args[item]

    # Reformat input args
    for item in arg_dict.keys():
        if arg_dict[item] is "None":
            arg_dict[item] = None

    if arg_dict['spec']:
        args.spec = args.spec.upper()



    pdb.set_trace()
    search_string = "skycell == 'p1911x04y12' and spec == 'f625w'"
    query_dataframe(args.csv_data_file, search_string)