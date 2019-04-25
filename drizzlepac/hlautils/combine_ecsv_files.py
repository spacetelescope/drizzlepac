#!/usr/bin/env python

"""Simple script to combine the ecsv files from individual test runs into a single monolithic file for easier
review.

"""
import argparse
from astropy.io import ascii
from astropy.table import Table, vstack
import glob
import os
from stsci.tools import logutil
import sys
import traceback



__taskname__ = 'combine_ecsv_files'

log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout)

# ------------------------------------------------------------------------------------------------------------

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
    #Search for ecsv files ignoring '_current' sym links to existing directories
    file_list = glob.glob("{}**/*[!current]/*.ecsv".format(input_file_basepath), recursive=True)
    n_found = len(file_list)
    if n_found == 0:
        log.info("No .ecsv files found. Exiting...")
        sys.exit("No .ecsv files found. Exiting...")
    elif n_found == 1:
        log.info("{} ecsv file found.".format(n_found))
    else:
        log.info("{} ecsv files found.".format(n_found))
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

    clobber : Boolean
        Overwrite existing files with the same name as output_filename?

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
        n_found = len(ecsv_file_list)
        if n_found == 1:
            file_plural_string = ""
        else:
            file_plural_string = "s"
        total_rows = len(out_data) #display total number of rows in output file.
        if total_rows == 1:
            row_plural_string = ""
        else:
            row_plural_string = "s"
        log.info("Wrote {} row{} from {} input file{} to output file {}".format(total_rows,
                                                                                row_plural_string,
                                                                                n_found,
                                                                                file_plural_string,
                                                                                output_filename))
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout) #display traceback


# ------------------------------------------------------------------------------------------------------------

def run_ecsv_combine(clobber=False,input_file_basepath=None,output_filename=None):
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
    # 0a: set up input arg defaults,

    if not input_file_basepath:
        input_file_basepath = os.getcwd()

    if not output_filename:
        output_filename = "{}/{}.ecsv".format(os.getcwd(),os.getcwd().split("/")[-1])

    if clobber == False and os.path.exists(output_filename) == True:
        sys.exit("Output file {} already exists. Rename the file or rerun with the clobber option on (-c) to "
              "overwrite.".format(output_filename))

    # 0c: make sure input_file_basepath always ends with a "/".
    if not input_file_basepath.endswith("/"):
        input_file_basepath += "/"

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

    parser.add_argument('-i', '--input_file_basepath', required=False, default=None,help='path to '
                        'start recursive search for the .ecsv files. If not specified, the current working '
                        'directory will be used.')

    parser.add_argument('-o', '--output_filename', required=False, default=None,help='Name of the '
                        'output combined .ecsv file. This may include a full file path. If not specified, the '
                        'file be named "<CURRENT WORKING DIRECTORY>/<CURRENT WORKING DIRECTORY>.ecsv".')
    args = parser.parse_args()


    run_ecsv_combine(clobber=args.clobber,
                     input_file_basepath=args.input_file_basepath,
                     output_filename=args.output_filename)
