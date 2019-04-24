#!/usr/bin/env python

"""Simple script to combine the ecsv files from individual test runs into a single monolithic file for easier
review.

"""
import argparse
from astropy.io import ascii
from astropy.table import Table, vstack
import glob
import os
import pdb
import sys
import traceback

def find_files(input_file_basepath):
    """Find the ecsv files.

    Parameters
    ----------
    input_file_basepath : string
        Path to start recursive search for the .ecsv files.

    Returns
    -------
    file_list: list
        List of ecsv full filenames
    """
    search_string = "{}popen-gw?/*[!current]/*.ecsv".format(input_file_basepath) #ignore '_current' sym links to existing directories
    file_list = glob.glob(search_string)

    if len(file_list) == 0: sys.exit("No .ecsv files found. Exiting...")

    return(file_list)

# ------------------------------------------------------------------------------------------------------------

def generate_output_file(ecsv_file_list,output_filename,clobber):
    """Generate combined output ecsv file.

    Parameters
    ----------
    ecsv_file_list: list
        List of ecsv full filenames

    output_filename: string
         Name of the output combined .ecsv file.

    Returns
    -------
    Nothing.
    """
    try:
        file_start = True
        for ecsv_filename in ecsv_file_list:
            table_data = ascii.read(ecsv_filename, format='ecsv') # Read ecsv file

            dataset = os.path.basename(ecsv_filename)[:-5].lower() # scrape dataset name out of ecsv filename
            dataset_column=Table.Column(name='datasetName', data=[dataset]*len(table_data)) # make new column
            table_data.add_column(dataset_column,index=0) # add dataset column to table data to append.

            if not file_start: #append out_data with ecsv file data for all files after the first in the list.
                out_data = vstack([out_data, table_data])
            else: # use the data from the first ecsv file to initialize out_data
                file_start = False
                out_data = table_data.copy()

        ascii.write(out_data,output_filename,format='ecsv',overwrite=clobber) #write output file.
        print("Wrote {}".format(output_filename))
        out_data.pprint(max_width=-1)

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout) #display traceback


# ------------------------------------------------------------------------------------------------------------

def run_combine(clobber=False,input_file_basepath="",output_filename=""):
    """Main calling subroutine.

    Parameters
    ----------
    clobber : Boolean
        Overwrite existing files with the same name as output_filename?

    input_file_basepath : string
        Path to start recursive search for the .ecsv files.

    output_filename: string
         Name of the output combined .ecsv file.

    Returns
    -------
    Nothing.
    """
    # 0: set up input arg defaults, make sure input_file_basepath always ends with a "/".
    if input_file_basepath == "":
        input_file_basepath = os.getcwd()

    if not input_file_basepath.endswith("/"):
        input_file_basepath += "/"

    if output_filename == "":
        output_filename = "{}/{}.ecsv".format(os.getcwd(),os.getcwd().split("/")[-1])

    # 1: create list of ecsv files to be combined.
    ecsv_file_list = find_files(input_file_basepath)


    # 2: combine individual ecsv files into a single monolithic file.
    generate_output_file(ecsv_file_list,output_filename,clobber)


# ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Find and combine all individual ecsv files into a single '
                                                 'file.')
    parser.add_argument('-c', '--clobber', required=False, action='store_true',help='If this option is turned '
                        'on, any existing file with same name as the output_fileanme will be overwritten.')

    parser.add_argument('-i', '--input_file_basepath', required=False, default="",help='path to '
                        'start recursive search for the .ecsv files. If not specified, the current working '
                        'directory will be used.')

    parser.add_argument('-o', '--output_filename', required=False, default="",help='Name of the '
                        'output combined .ecsv file. This may include a full file path. If not specified, the '
                        'file be named "<CURRENT WORKING DIRECTORY>/<CURRENT WORKING DIRECTORY>.ecsv".')
    args = parser.parse_args()


    run_combine(clobber=args.clobber,input_file_basepath=args.input_file_basepath,output_filename=args.output_filename)
