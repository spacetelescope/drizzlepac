""" This module tests full pipeline SVM processing - visit has one detector and uses one filter.

"""
import datetime
import glob
import os
import pytest

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

    *** The --basetemp does not seem to be working at this time!!!

"""

WCS_SUB_NAME = "FIT_SVM_GAIA"
POLLER_FILE = "acs_8f6_61_input.out"

@pytest.fixture(scope="module")
def read_csv_for_filenames():
    # Read the CSV poller file residing in the tests directory to extract the individual visit FLT/FLC filenames
    path = os.path.join(os.path.dirname(__file__), POLLER_FILE)
    table = ascii.read(path, format="no_header")
    filename_column = table.colnames[0]
    filenames = list(table[filename_column])
    print("\nread_csv_for_filenames. Filesnames from poller: {}".format(filenames))

    return filenames


@pytest.fixture(scope="module")
def gather_data_for_processing(read_csv_for_filenames):
    # Establish FLC/FLT lists and obtain the requested data 
    flc_list = []
    flt_list = []
    flcfiles = []
    fltfiles = []
    for fn in read_csv_for_filenames:
        if fn.lower().endswith("flc.fits"):
            flc_list.append(fn.split("_")[0])
        elif fn.lower().endswith("flt.fits"):
            flt_list.append(fn.split("_")[0])

    # Get test data through astroquery - only retrieve the pipeline processed FLC/FLT files
    # (e.g., j*_flc.fits).  Soon the SVM and MVM products will be pipeline produced!
    if flc_list:
        flcfiles = aqutils.retrieve_observation(flc_list, suffix=["FLC"], product_type="pipeline")
    if flt_list:
        fltfiles = aqutils.retrieve_observation(flt_list, suffix=["FLT"], product_type="pipeline")

    flcfiles.extend(fltfiles)
    print("\ngather_data_for_processing. Gathered data: {}".format(flcfiles))

    return flcfiles


@pytest.fixture(scope="module")
def gather_output_data(construct_manifest_filename):
    # Determine the filenames of all the output files from the manifest
    table = ascii.read(construct_manifest_filename, format="no_header")
    file_col = table.colnames[0]
    files = list(table[file_col])
    print("\ngather_output_data. Output data files: {}".format(files))

    return files


@pytest.fixture(scope="module")
def construct_manifest_filename(read_csv_for_filenames):
    # Construct the output manifest filename from input file keywords
    inst = fits.getval(read_csv_for_filenames[0], "INSTRUME", ext=0).lower()
    root = fits.getval(read_csv_for_filenames[0], "ROOTNAME", ext=0).lower()
    tokens_tuple = (inst, root[1:4], root[4:6], "manifest.txt")
    manifest_filename = "_".join(tokens_tuple)
    print("\nconstruct_manifest_filename. Manifest filename: {}".format(manifest_filename))

    return manifest_filename


@pytest.fixture(scope="module", autouse=True)
def svm_setup(gather_data_for_processing):
    # Act: Process the input data by executing runsinglehap - time consuming activity

    current_dt = datetime.datetime.now()
    print(str(current_dt))
    print("\nsvm_setup fixture")

    # Read the "poller file" and download the input files, as necessary
    input_names = gather_data_for_processing

    # Run the SVM processing
    try:
        status = runsinglehap.perform(POLLER_FILE)

    # Catch anything that happens and report it.  This is meant to catch unexpected errors and
    # generate sufficient output exception information so algorithmic problems can be addressed.
    except Exception as except_details:
        print(except_details)
        pytest.fail("\nsvm_setup. Exception Visit: {}\n", POLLER_FILE)

    current_dt = datetime.datetime.now()
    print(str(current_dt))


# TESTS

#def test_svm_manifest_name(gather_data_for_processing):
def test_svm_manifest_name(construct_manifest_filename):
    # Construct the manifest filename from the header of an input file in the list and check it exists.
    path = Path(construct_manifest_filename)
    print("\ntest_svm_manifest. Filename: {}".format(path))

    # Ensure the manifest file uses the proper naming convention
    assert(path.is_file())


def test_svm_wcs(gather_output_data):
    # Check the output primary WCSNAME includes FIT_SVM_GAIA as part of the string value
    tdp_files = [files for files in gather_output_data if files.lower().find("total") and files.lower().endswith(".fits")]

    for tdp in tdp_files:
        wcsname = fits.getval(tdp, "WCSNAME", ext=1).upper()
        print("\ntest_svm_wcs.  WCSNAME: {} Output file: {}".format(wcsname, tdp))
        assert WCS_SUB_NAME in wcsname, f"WCSNAME is not as expected for file {tdp}."


def test_svm_cat_sources(gather_output_data):
    # Check the output catalogs should contain > 0 measured sources
    cat_files = [files for files in gather_output_data if files.lower().endswith("-cat.ecsv")]

    for cat in cat_files:
        table_length = len(ascii.read(cat, format="ecsv"))
        print("\ntest_svm_cat_sources. Number of sources in catalog {} is {}.".format(cat, table_length))
        assert table_length > 0, f"Catalog file {cat} is unexpectedly empty"
