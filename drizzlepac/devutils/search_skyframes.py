#!/usr/bin/env python

""" This script allows the user to query a master observations .csv file by any combination of exposure name,
    skycell name, instrument/detector, and/or filter name. The results are printed the screen and optionally
    written to a user-specified output .csv file. The results will be sorted first by skycell name, then by
    instrument/detector, filter name and finally by image name."""

import argparse
import os
import pdb
import sys

from astropy.io import fits
import pandas as pd

from drizzlepac.haputils import cell_utils
from drizzlepac.haputils import make_poller_files


# ------------------------------------------------------------------------------------------------------------

def augment_results(results):
    """Add additional searchable columns to query results. Columns added are date-obs and  full flt/flc fits
    file path

    Parameters
    ----------
    results : pandas.DataFrame
        query results to augment

    Returns
    -------
    results : pandas.DataFrame
        Augmented query results
    """
    dateobs_list = []
    path_list = []
    basepath = os.getcwd()
    for idx in results.index:
        rootname = results.exposure[idx]
        filestub="{}/{}".format(basepath, rootname)
        if os.path.exists("{}_flc.fits".format(filestub)):
            imgname = "{}_flc.fits".format(filestub)
        else:
            imgname = "{}_flt.fits".format(filestub)
        dateobs_list.append(fits.getval(imgname, "DATE-OBS"))
        path_list.append(imgname)
    results['dateobs'] = dateobs_list
    results['dateobs'] = pd.to_datetime(results.dateobs)
    results['filename'] = path_list
    pdb.set_trace()
    return results


# ------------------------------------------------------------------------------------------------------------

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
        if item is "skycell":
            substring = "{} == '{}'".format(item, arg_dict[item])
        else:
            substring = "{}.str.contains('{}')".format(item, arg_dict[item])
        if len(search_string) == 0:
            query_delimiter = ""
        else:
            query_delimiter = " and "
        search_string = "{}{}{}".format(search_string, query_delimiter, substring)
    return search_string


# ------------------------------------------------------------------------------------------------------------

def query_dataframe(master_observations_file, search_string, output_columns=None, output_sorting=None,
                    output_filename=None):
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

    output_columns : list, optional
        Optional list describing which query result columns to display, and in what order. If not explicitly
        specified, all columns will be displayed.

    output_sorting : list, optional
        Optional list describing how, and in what order query results should be sorted. If not explicitly
        specified, no sorting will be performed on query results.

    output_filename : str, optional
        Optional name of an output .csv file to write the query results to.

    Returns
    -------
    results : Pandas DataFrame class
        A Pandas DataFrame continaing the query results.
    """
    dataframe = pd.read_csv(master_observations_file, header=0, index_col=0)
    results = dataframe.query(search_string, engine='python')
    results = augment_results(results)
    ret_results = results.copy()
    if output_sorting:
        results = results.sort_values(by=output_sorting)
    if output_columns:
        results = results[output_columns]
    n_lines = results.index.size
    print("\n")
    print(results.to_string())
    print("\n")
    if n_lines == 1:
        print('Query "{}" found 1 result.'.format(search_string))
    else:
        print('Query "{}" found {} results.'.format(search_string, n_lines))
    if output_sorting:
        col_sorting_str = ", ".join(output_sorting)
    else:
        col_sorting_str = "No sorting"
    print("Column sorting order: {}".format(col_sorting_str))
    if output_filename:
        results.to_csv(output_filename)
        print("Wrote query results to {}".format(output_filename))

    return ret_results

# ------------------------------------------------------------------------------------------------------------


def visualize_footprints(results):
    """Visualize footprints of skycells and exposures in query result

    Parameters
    ----------
    results : Pandas DataFrame class
        Pandas DataFrame containing query results.

    Returns
    -------
    Nothing.
    """
    # Get list of all the unique skycells in the query results
    unique_skycell_list = []
    for idx in results.index:
        unique_skycell_list.append(results.skycell[idx])
    unique_skycell_list = list(set(unique_skycell_list))  # remove duplicate items from skycell list

    for skycell_name in unique_skycell_list:
        skycell = cell_utils.SkyCell.from_name("skycell-{}".format(skycell_name))
        footprint = cell_utils.SkyFootprint(meta_wcs=skycell.wcs)
        footprint.build(img_path_dict[skycell_name])
        footprint_imgname = "skycell-{}-footprint.fits".format(skycell_name)
        foo = footprint.get_footprint_hdu(filename=footprint_imgname)
        print("Skycell footprint image {} contains {} individual exposures.".format(footprint_imgname,
                                                                                    len(img_path_dict[skycell_name])))
    print("\a\a\a")

# ------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    # Argparse input
    parser = argparse.ArgumentParser(description='Search all observations')
    parser.add_argument('-m', '--master_observations_file', required=False,
                        default="~/Documents/HLAtransition/mvm_testing/allexp/allexposures.csv",
                        help='Name of the master observations .csv file containing comma-separated columns '
                             'index #, exposure, skycell, config, and spec. '
                             'Default value is "allexposures.csv".')
    parser.add_argument('-i', '--exposure', required=False, default="None",
                        help='Image name. Does not have to be a full 9-character iamge name. Partials are '
                             'acceptable.')
    parser.add_argument('-s', '--skycell', required=False, default="None",
                        help='Skycell name. Only full skycell names are accepted.')
    parser.add_argument('-d', '--config', required=False, default="None", choices=["ACS/HRC", "ACS/SBC",
                                                                                   "ACS/WFC", "WFC3/IR",
                                                                                   "WFC3/UVIS", "None"],
                        help='Instrument/Detector name. All caps, separated by a "/" (e.g. WFC3/IR).')

    parser.add_argument('-f', '--spec', required=False, default="None", help='Filter name.')
    parser.add_argument('-o', '--output_results_file', required=False, default="None",
                        help='Optional name of an output .csv file to write the query results to. If not '
                             'explicitly specified, no output file will be written.')
    parser.add_argument('-v', '--visualize_footprints', required=False, action='store_true',
                        help='If turned on, the footprints of the skycell and the exposures returned by the '
                             'query will be displayed.')
    parser.add_argument('--output_columns', required=False, default="None", nargs='?',
                        help="Columns to display in the query results. Needs to be some combination of "
                             "'dateobs', 'exposure', 'skycell', 'config' and/or 'spec'. If not explicitly "
                             "specified, columns will be displayed in the query results.")
    parser.add_argument('--output_sorting', required=False, default="None", nargs='?',
                        help="Order (if any) in which to sort columns of the query results. Needs to be some "
                             "combination of 'dateobs', 'exposure', 'skycell', 'config' and/or 'spec'. If "
                             "not explicitly specified, query results will not be sorted. Recommended "
                             "sorting order is 'skycell,config,spec,exposure'.")
    in_args = parser.parse_args()
    arg_dict = {}
    for item in in_args.__dict__.keys():
        if item not in ["master_observations_file", "output_results_file",
                        "output_columns", "output_sorting", "visualize_footprints"]:
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

    if in_args.output_columns == "None":
        in_args.output_columns = None
    else:
        in_args.output_columns = in_args.output_columns.split(",")

    if in_args.output_sorting == "None":
        in_args.output_sorting = None
    else:
        in_args.output_sorting = in_args.output_sorting.split(",")

    # FAULT TOLERANCE. Exit if something in the input arguments isn't quite right.
    # Exit if user didn't enter any search criteria
    if blank_entry_ctr == 4:
        sys.exit("ERROR: Search results too broad. No search criteria were entered.")

    # Exit if in_args.master_observations_file == in_args.output_results_file so the master observations table doesn't get overwritten
    if in_args.master_observations_file == in_args.output_results_file:
        sys.exit("ERROR: The output results file cannot have the same name is the input master observations table")

    # Exit if named columns don't match the names of the available columns.
    valid_cols = ['dateobs', 'exposure', 'skycell', 'config', 'spec']
    if in_args.output_columns is not None:
        for item in in_args.output_columns:
            if item not in valid_cols:
                sys.exit("ERROR: {} is not a valid column name. Valid column names are 'dateobs', 'exposure', 'skycell', 'config', 'spec'.".format(item))

    if in_args.output_sorting is not None:
        for item in in_args.output_sorting:
            if item not in valid_cols:
                sys.exit("ERROR: {} is not a valid column name. Valid column names are 'dateobs', 'exposure', 'skycell', 'config', 'spec'.".format(item))

    # Generate search string and run query
    search_string = make_search_string(arg_dict)

    results = query_dataframe(in_args.master_observations_file,
                              search_string,
                              output_columns=in_args.output_columns,
                              output_sorting=in_args.output_sorting,
                              output_filename=in_args.output_results_file)

    # Optionally visualize footprints of skycells and exposures returned by the query.
    if in_args.visualize_footprints:
        if results.index.size > 0:
            visualize_footprints(results)
        else:
            print("WARNING: Null query result. No footprints to visualize.")