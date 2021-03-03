#!/usr/bin/env python

""" This script allows the user to query a master observations .csv file by any combination of exposure name,
    skycell name, instrument/detector, and/or filter name. The results are printed the screen and optionally
    written to a user-specified output .csv file. The results will be sorted first by skycell name, then by
    instrument/detector, filter name and finally by image name."""

import argparse
import sys
import pandas as pd


def make_search_string(arg_dict):
    """Create search string that will be used to query observation tables

    Parameters
    ----------
    arg_dict: dict
        dictionary with the search arguments. They are as follows: 'exposure': Image name,
        'skycell': Skycell name, 'config': Instrument/Detector name (Uppercase, separated by a "/"
        (e.g. WFC3/IR)), 'spec': Filter name. lowercase.

    Returns
    -------
    search_string : str
        Properly formatted search string
    """
    search_string = ""
    for item in arg_dict:
        if not arg_dict[item]:
            continue
        substring = "{} == '{}'".format(item, arg_dict[item])
        if len(search_string) == 0:
            query_delimiter = ""
        else:
            query_delimiter = " and "
        search_string = "{}{}{}".format(search_string, query_delimiter, substring)
    return search_string


# ------------------------------------------------------------------------------------------------------------

def query_dataframe(master_observations_file, search_string, out_file=None):
    """Read in master observations file to pandas dataframe, run query and display (and optionally) write out
    the query results. The results will be sorted first by skycell name, then by instrument/detector, filter
    name and finally by image name.

    Parameters
    ----------
    master_observations_file : str
        Name of the master observations .csv file containing comma-separated columns index #, exposure,
        skycell, config, and spec.

    search_string : str
        properly formatted search string to use in the query of the master observations file.

    out_file : str, optional
        Optional name of an output .csv file to write the query results to.

    Returns
    -------
    Nothing.
    """
    dataframe = pd.DataFrame.from_csv(master_observations_file, header=0, index_col=0)
    results = dataframe.query(search_string)
    results = results.sort_values(by=['skycell', 'config', 'spec', 'exposure'])
    results = results[['exposure', 'skycell']]
    n_lines = results.index.size
    print("\n")
    print(results.to_string())
    print("\n")
    if n_lines == 1:
        print('Query "{}" found 1 result.'.format(search_string))
    else:
        print('Query "{}" found {} results.'.format(search_string, n_lines))
    if out_file:
        results.to_csv(out_file)
        print("Wrote query results to {}".format(out_file))


# ------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    # Argparse input
    parser = argparse.ArgumentParser(description='Search all observations')
    parser.add_argument('-m', '--master_observations_file', required=False, default="allexposures.csv",
                        help='Name of the master observations .csv file containing comma-separated columns '
                             'index #, exposure, skycell, config, and spec. '
                             'Default value is "allexposures.csv".')
    parser.add_argument('-i', '--exposure', required=False, default="None", help='Image name.')
    parser.add_argument('-s', '--skycell', required=False, default="None", help='Skycell name.')
    parser.add_argument('-d', '--config', required=False, default="None", choices=["ACS/HRC", "ACS/SBC",
                                                                                   "ACS/WFC", "WFC3/IR",
                                                                                   "WFC3/UVIS", "None"],
                        help='Instrument/Detector name. All caps, separated by a "/" (e.g. WFC3/IR).')

    parser.add_argument('-f', '--spec', required=False, default="None", help='Filter name.')
    
    parser.add_argument('-o', '--output_results_file', required=False, default="None",
                        help='Optional name of an output .csv file to write the query results to. If not '
                             'explicitly specified, no output file will be written.')
    in_args = parser.parse_args()
    arg_dict = {}
    for item in in_args.__dict__.keys():
        if item not in ["master_observations_file", "output_results_file"]:
            arg_dict[item] = in_args.__dict__[item]

    # Reformat input args
    blank_entry_ctr = 0
    for item in arg_dict.keys():
        if arg_dict[item] == "None":
            arg_dict[item] = None
            blank_entry_ctr += 1
    if arg_dict['spec']:
        arg_dict['spec'] = arg_dict['spec'].lower()
    if in_args.output_results_file == "None":
        in_args.output_results_file = None

    # Bail out if user didn't enter any search criteria
    if blank_entry_ctr == 4:
        sys.exit("ERROR: Search results too broad. No search criteria were entered.")

    # Bail out if in_args.master_observations_file == in_args.output_results_file so the master observations table doesn't get overwritten
    if in_args.master_observations_file == in_args.output_results_file:
        sys.exit("ERROR: The output results file cannot have the same name is the input master observations table")
    search_string = make_search_string(arg_dict)

    query_dataframe(in_args.master_observations_file, search_string, out_file=in_args.output_results_file)
