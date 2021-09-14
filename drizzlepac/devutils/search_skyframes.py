#!/usr/bin/env python

""" This script allows users to query a master observations .csv file by any combination of exposure name,
    skycell name, instrument/detector, and/or filter name. The results are printed the screen and optionally
    written to a user-specified output .csv file. The results will be sorted first by skycell name, then by
    instrument/detector, filter name and finally by image name. Optionally, footprint fits files will be
    generated that show the footprints of the images in the query result in the skycell in which they
    reside. These multi-detector, multi-filter composite footprint images use the naming format "
    "skycell-<SKYCELL NAME>_all_footprint.fits". Additional footprint images are created for each
    skycell-instrument/detector/filter combination present in the query results. These
    skycell/detector/filter-specific footprint images use the naming format
    "skycell-<SKYCELL NAME>_<INSTRUMENT NAME>-<DETECTOR NAME>-<FILTER NAME_footprint.fits".


    USAGE:

        >>> search_skyframes [-cdefmosv --output_columns --output_sorting]

        - The '-c' option allows users to specify an instrument/detector configuration to search for. The
        expected format is instrument name, then detector name in all caps, separated by a "/"; example:
        "-c WFC3/IR"

        - The '-d' option allows users to specify the Min and max date (inclusive) of dateobs search
        window. The exected date format is "YYYY-MM-DD". Min amd max dates should be separated by a space;
        example: "-d 2011-06-15 2019-04-21"

        - The '-e' option allows users to specify an exposure name to search for. The expected format is a
        lower-case full or partial image name; examples: "-e j92", "-e j92c01" or "-e j92c01b4q"

        - The '-f' option allows users to specify filter name(s) to search for. To search for ACS
        observations that use two spectral elements, enter the names of both spectral elements in any order
        seperated by a dash. The expected format is filter name(s), lower-case; Example ACS two-spectral
        element input: "-f f606w-pol60v", Example single-filter input: "-f f814w"

        - The '-m' option allows users to specify the name of the master observations .csv file containing
        comma-separated columns "index #", "exposure", "skycell", "config", and "spec" that are searched to
        produce query results.

        - The '-o' option allows users to specify an output file to write the query results to. The output
        file format is '.csv'.

        - The '-s' option allows users to specify the skycell name to search for. Only full skycell names are
        accepted. Expected format is lowercase; Example: "-s p0797x13y05"

        - The '-v' opition allows users to visualize the footprints of the images returned by the query
        results. If turned on, footprint fits files will be generated that show the footprints of the images
        returned by the query result in the skycell in which they reside. Additional footprint images are
        created for each skycell-instrument/detector/filter combination present in the query results. This
        input arguemnt does not expect any additinal information. Simply specifying "-v" turns on the
        footprint generation.

        - The '--output_columns' option allows users to specify which columns to display (and optionally
        write to the output file, and in what order.  Valid column names are 'dateobs', 'exposure',
        'filename', 'skycell', 'config', and 'spec'; Example: "--output_columns skycell,config,spec,exposure"

        - The '--output_sorting' option allows users Order in which to sort columns of the query results.
        The expected input is to be some combination of 'dateobs', 'exposure', 'skycell', 'config', 'spec'
        and 'filename'. If not explicitly specified, query results will not be sorted; Example:
        "--output_sorting skycell,config,spec"

    EXAMPLE RUN:

    Search for all ACS/WFC F814W observations located in skycell p0123x45y67 taken during 2015; return
    date and exposure name; sort by date; write results to file 'query_results.csv'; generate footprint
    visualization footprint fits imaages

        >>> search_skyframes -c ACS/WFC -d 2015-01-01 2015-12-31 -f f814w -s p0123x45y67 -o query_results.csv -v --output_columns dateobs, exposure --output_sorting dateobs
    """

import argparse
from datetime import datetime
import os
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

    # populate dateobs_list, path_list
    for idx in results.index:
        rootname = results.exposure[idx]
        imgname = make_poller_files.locate_fitsfile(rootname)
        dateobs_list.append(fits.getval(imgname, "DATE-OBS"))
        path_list.append(imgname)

    # Add 'dateobs' (as datetime objects) and 'filename' columns to results dataframe.
    results['dateobs'] = dateobs_list
    results['dateobs'] = pd.to_datetime(results.dateobs)
    results['filename'] = path_list
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
        if not arg_dict[item] or item is "date_range":
            continue
        if item is "skycell":
            substring = "{} == '{}'".format(item, arg_dict[item])
        elif item is "spec":
            parse_spec = arg_dict[item].split("-")
            substring = "{}.str.contains('{}')".format(item, parse_spec[0])
            if len(parse_spec) == 2:  # Search for two-filter observations (ACS)
                substring = "{} and {}.str.contains('{}')".format(substring, item, parse_spec[1])
        else:
            substring = "{}.str.contains('{}')".format(item, arg_dict[item])
        if len(search_string) == 0:
            query_delimiter = ""
        else:
            query_delimiter = " and "
        search_string = "{}{}{}".format(search_string, query_delimiter, substring)
    return search_string


# ------------------------------------------------------------------------------------------------------------

def query_dataframe(master_observations_file, search_string, date_range=None, output_columns=None,
                    output_sorting=None, output_filename=None):
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

    date_range : list
        Optional two-element list containing the start and end dates (inclusive) of the date range to
        query. Date format: YYYY-MM-DD

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
        A Pandas DataFrame containing the query results.
    """
    dataframe = pd.read_csv(master_observations_file, header=0, index_col=0)
    results = dataframe.query(search_string, engine='python')
    results = augment_results(results)
    if date_range:
        search_string = "{} and '{}' <= dateobs <= '{}'".format(search_string, date_range[0], date_range[1])
        results = results[(results['dateobs'] >= date_range[0]) & (results['dateobs'] <= date_range[1])]
    ret_results = results.copy()
    if output_sorting:
        results = results.sort_values(by=output_sorting)
    if output_columns:
        results = results[output_columns]
    n_lines = results.index.size

    print("\n")
    print(results.to_string(index=False))
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
    """Visualize footprints of skycells and exposures in query result. footprint fits files will be
    generated that show the footprints of the images in the query result in the skycell in which they
    reside. These multi-detector, multi-filter composite footprint images use the naming format "
    "skycell-<SKYCELL NAME>_all_footprint.fits". Additional footprint images are created for each
    skycell-instrument/detector/filter combination present in the query results. These
    skycell/detector/filter-specific footprint images use the naming format
    "skycell-<SKYCELL NAME>_<INSTRUMENT NAME>-<DETECTOR NAME>-<FILTER NAME_footprint.fits".

    Parameters
    ----------
    results : Pandas DataFrame class
        Pandas DataFrame containing query results.

    Returns
    -------
    Nothing.
    """
    print("Building footprint images. Please be patient. This may take some time...")
    for skycell_name in list(set(results.skycell.values.tolist())):
        footprint_imgname = "skycell-{}_all_footprint.fits".format(skycell_name)
        img_list = results.query('skycell == "{}"'.format(skycell_name)).filename.values.tolist()
        make_footprint_fits_file(skycell_name, img_list, footprint_imgname)
        for config_name in list(set(results.query('skycell == "{}"'.format(skycell_name)).config.values.tolist())):
            for spec_name in list(set(results.query('skycell == "{}" and config == "{}"'.format(skycell_name, config_name)).spec.values.tolist())):
                img_list = results.query('skycell == "{}" and config == "{}" and spec == "{}"'.format(skycell_name, config_name, spec_name)).filename.values.tolist()
                footprint_imgname = "skycell-{}_{}_{}_footprint.fits".format(skycell_name, config_name.replace("/", "-"), spec_name)
                make_footprint_fits_file(skycell_name, img_list, footprint_imgname)


# ------------------------------------------------------------------------------------------------------------


def make_footprint_fits_file(skycell_name, img_list, footprint_imgname):
    """Generates footprint fits file

    Parameters
    ----------
    skycell_name : str
        Name of the skycell

    img_list : list
        list containing the images to process into a footprint fits files

    footprint_imgname : str
        name of the footprint fits file to generate

    Returns
    -------
    Nothing.
    """
    skycell = cell_utils.SkyCell.from_name("skycell-{}".format(skycell_name))
    footprint = cell_utils.SkyFootprint(meta_wcs=skycell.wcs)
    footprint.build(img_list)
    foo = footprint.get_footprint_hdu(filename=footprint_imgname)
    n_images = len(img_list)
    if n_images > 1:
        add_an_s = "s"
    else:
        add_an_s = ""
    print("Skycell footprint image {} contains {} individual exposure{}.".format(footprint_imgname, n_images,
                                                                                 add_an_s))


# ------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    # Argparse input
    parser = argparse.ArgumentParser(description='Search all observations')
    parser.add_argument('-c', '--config', required=False, default="None", choices=["ACS/HRC", "ACS/SBC",
                                                                                   "ACS/WFC", "WFC3/IR",
                                                                                   "WFC3/UVIS", "None"],
                        help='Instrument/Detector configuration to search for. All caps, separated by a "/" '
                             '(e.g. WFC3/IR).')
    parser.add_argument('-d', '--date_range', required=False, default=None, nargs=2,
                        help='Min and max date (inclusive) of dateobs search window. FORMAT: Date format is '
                             '"YYYY-MM-DD". min amd max dates should be separated by a space.')
    parser.add_argument('-e', '--exposure', required=False, default="None",
                        help='Exposure name to search for. Does not have to be a full 9-character exposure '
                             'name. Partials are acceptable.')
    parser.add_argument('-f', '--spec', required=False, default="None",
                        help='Filter name(s) to search for. To search for ACS observations that use two '
                             'spectral elements, enter the names of both spectral elements in any order '
                             'separated by a dash. Example two-spectral element input: f606w-pol60v')
    parser.add_argument('-m', '--master_observations_file', required=False,
                        default=os.getenv("ALL_EXP_FILE"),
                        help='Name of the master observations .csv file containing comma-separated columns '
                             'index #, exposure, skycell, config, and spec. '
                             'Default value is "allexposures.csv".')
    parser.add_argument('-o', '--output_results_file', required=False, default="None",
                        help='Optional name of an output .csv file to write the query results to. If not '
                             'explicitly specified, no output file will be written.')
    parser.add_argument('-s', '--skycell', required=False, default="None",
                        help='Skycell name to search for. Only full skycell names are accepted.')
    parser.add_argument('-v', '--visualize_footprints', required=False, action='store_true',
                        help='If turned on, footprint fits files will be generated that show the footprints '
                             'of the images returned by the query result in the skycell in which they '
                             'reside.')
    parser.add_argument('--output_columns', required=False, default="None", nargs='?',
                        help="Columns to display in the query results. Needs to be some combination of "
                             "'dateobs', 'exposure', 'skycell', 'config' and/or 'spec'. If not explicitly "
                             "specified, all columns will be displayed in the query results.")
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
    if blank_entry_ctr == 5:
        sys.exit("ERROR: Search results too broad. No search criteria were entered.")
    # Exit if specified dates don't match the format "YYYY-MM-DD".
    if arg_dict["date_range"]:
        date_range_dt = []
        for date_string, date_label in zip(arg_dict["date_range"], ["Start", "End"]):
            try:
                date_range_dt.append(datetime.strptime(date_string, "%Y-%m-%d"))
            except ValueError:
                sys.exit("{} date value '{}' is not properly formatted. Please use the format 'YYYY-MM-DD'.".format(date_label, date_string))
        # flip order of dates if first is later than the second.
        if date_range_dt[1] < date_range_dt[0]:
            arg_dict["date_range"].reverse()

    # Exit if in_args.master_observations_file is not found.
    if not in_args.master_observations_file:
        sys.exit("ERROR: Undefined master observations table path. Please set environment variable 'ALL_EXP_FILE'")

    # Exit if the specified master observations table doesn't exist
    if not os.path.exists(in_args.master_observations_file):
        sys.exit("ERROR: Specified master observations table '{}' does not exist.".format(in_args.master_observations_file.replace(os.getenv("ALL_EXP_FILE"), "$ALL_EXP_FILE/")))

    # Exit if in_args.master_observations_file == in_args.output_results_file so the master observations table doesn't get overwritten
    if in_args.master_observations_file == in_args.output_results_file:
        sys.exit("ERROR: The output results file cannot have the same name is the input master observations table")

    # Exit if named columns don't match the names of the available columns.
    valid_cols = ['dateobs', 'exposure', 'filename', 'skycell', 'config', 'spec']
    if in_args.output_columns is not None:
        for item in in_args.output_columns:
            if item not in valid_cols:
                sys.exit("ERROR: {} is not a valid column name. Valid column names are 'dateobs', 'exposure', 'filename', 'skycell', 'config', 'spec'.".format(item))
    if in_args.output_sorting is not None:
        for item in in_args.output_sorting:
            if item not in valid_cols:
                sys.exit("ERROR: {} is not a valid column name. Valid column names are 'dateobs', 'exposure', 'filename', 'skycell', 'config', 'spec'.".format(item))

    # Generate search string and run query
    search_string = make_search_string(arg_dict)

    results = query_dataframe(in_args.master_observations_file,
                              search_string,
                              date_range=arg_dict["date_range"],
                              output_columns=in_args.output_columns,
                              output_sorting=in_args.output_sorting,
                              output_filename=in_args.output_results_file)

    # Optionally visualize footprints of skycells and exposures returned by the query.
    if in_args.visualize_footprints:
        if results.index.size > 0:
            visualize_footprints(results)
        else:
            print("WARNING: Null query result. No footprints to visualize.")
