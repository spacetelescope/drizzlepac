""" This module tests full pipeline SVM processing as a demonstration template.

"""
import datetime
import glob
import os
import pytest
import numpy as np

from drizzlepac.haputils import astroquery_utils as aqutils
from drizzlepac import runsinglehap
from astropy.io import fits, ascii
from pathlib import Path

"""
    test_svm_demo.py

    This test file can be executed in the following manner:
        $ pytest -s --basetemp=/internal/hladata/yourUniqueDirectoryHere test_svm.py >& test_svm.log &
        $ tail -f test_svm.log
      * Note: When running this test, the `--basetemp` directory should be set to a unique
        existing directory to avoid deleting previous test output.
      * The POLLER_FILE exists in the tests/hap directory.

"""

WCS_SUB_NAME = "FIT_SVM_GAIA"
POLLER_FILE = "acs_hrc_sbc_input.out"

def read_csv_for_filenames():
    # Read the CSV poller file residing in the tests directory to extract the individual visit FLT/FLC filenames
    path = os.path.join(os.path.dirname(__file__), POLLER_FILE)
    table = ascii.read(path, format="no_header")
    filename_column = table.colnames[0]
    filenames = list(table[filename_column])
    print("\nread_csv_for_filenames. Filesnames from poller: {}".format(filenames))

    return filenames


def gather_data_for_processing(tmp_path_factory):
    # create working directory specified for the test
    curdir = tmp_path_factory.mktemp(os.path.basename(__file__))
    os.chdir(curdir)

    # Establish FLC/FLT lists and obtain the requested data
    flc_flag = ""
    flt_flag = ""
    # In order to obtain individual FLC or FLT images from MAST (if the files are not reside on disk) which
    # may be part of an ASN, use only IPPPSS with a wildcard.  The unwanted images have to be removed
    # after-the-fact.
    filenames = read_csv_for_filenames()

    for fn in filenames:
        if fn.lower().endswith("flc.fits") and flc_flag == "":
            flc_flag = fn[0:6] + "*"
        elif fn.lower().endswith("flt.fits") and flt_flag == "":
            flt_flag = fn[0:6] + "*"

        # If both flags have been set, then break out the loop early.  It may be
        # that all files have to be checked which means the for loop continues
        # until its natural completion.
        if flc_flag and flt_flag:
            break

    # Get test data through astroquery - only retrieve the pipeline processed FLC and/or FLT files
    # (e.g., j*_flc.fits) as necessary. The logic here and the above for loop is an attempt to
    # avoid downloading too many images which are not needed for processing.
    flcfiles = []
    fltfiles = []
    if flc_flag:
        flcfiles = aqutils.retrieve_observation(flc_flag, suffix=["FLC"], product_type="pipeline")
    if flt_flag:
        fltfiles = aqutils.retrieve_observation(flt_flag, suffix=["FLT"], product_type="pipeline")

    flcfiles.extend(fltfiles)

    # Keep only the files which exist in BOTH lists for processing
    files_to_process = set(filenames).intersection(set(flcfiles))

    # Identify unwanted files from the download list and remove from disk
    files_to_remove = set(filenames).symmetric_difference(set(flcfiles))
    try:
        for ftr in files_to_remove:
            os.remove(ftr)
    except Exception as x_cept:
        print("")
        print("Exception encountered: {}.".format(x_cept))
        print("The file {} could not be deleted from disk. ".format(ftr))
        print("Remove files which are not used for processing from disk manually.")

    print("\ngather_data_for_processing. Gathered data: {}".format(files_to_process))

    return list(files_to_process)


def gather_output_data(manifest_filename):
    # Determine the filenames of all the output files from the manifest
    print(f"\nManifest Filename: {manifest_filename}")
    files = []
    with open(manifest_filename, 'r') as fout:
        for line in fout.readlines():
            files.append(line.rstrip('\n'))
    print("\ngather_output_data. Output data files: {}".format(files))

    return files


def construct_manifest_filename(filenames):
    # Construct the output manifest filename from input file keywords
    inst = fits.getval(filenames[0], "INSTRUME", ext=0).lower()
    root = fits.getval(filenames[0], "ROOTNAME", ext=0).lower()
    tokens_tuple = (inst, root[1:4], root[4:6], "manifest.txt")
    manifest_filename = "_".join(tokens_tuple)
    print("\nconstruct_manifest_filename. Manifest filename: {}".format(manifest_filename))

    return manifest_filename


def test_driver(tmp_path_factory):
    # Act: Process the input data by executing runsinglehap - time consuming activity

    current_dt = datetime.datetime.now()
    print(str(current_dt))

    # Read the "poller file" and download the input files, as necessary
    input_names = gather_data_for_processing(tmp_path_factory)

    # Construct the manifest filename for later
    manifest_filename = construct_manifest_filename(input_names)

    # Run the SVM processing
    path = os.path.join(os.path.dirname(__file__), POLLER_FILE)
    try:
        status = runsinglehap.perform(path, log_level="debug")

        output_files = gather_output_data(manifest_filename)

        # Check the output primary WCSNAME includes FIT_SVM_GAIA as part of the string value
        tdp_files = [files for files in output_files if
                     files.lower().find("total") > -1 and files.lower().endswith(".fits")]

        for tdp in tdp_files:
            wcsname = fits.getval(tdp, "WCSNAME", ext=1).upper()
            print("\ntest_svm_wcs.  WCSNAME: {} Output file: {}".format(wcsname, tdp))
            assert WCS_SUB_NAME in wcsname, f"WCSNAME is not as expected for file {tdp}."

    # Catch anything that happens and report it.  This is meant to catch unexpected errors and
    # generate sufficient output exception information so algorithmic problems can be addressed.
    except Exception as except_details:
        print(except_details)
        pytest.fail("\nsvm_setup. Exception Visit: {}\n", path)

    current_dt = datetime.datetime.now()
    print(str(current_dt))
