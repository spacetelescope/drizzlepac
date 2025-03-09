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
    test_svm_u2as01.py

    This test file can be executed in the following manner:
        $ pytest -s --basetemp=/internal/hladata/yourUniqueDirectoryHere test_svm_u2as01.py >& test_svm_u2as01.log &
        $ tail -f test_svm_u2as01.log
      * Note: When running this test, the `--basetemp` directory should be set to a unique
        existing directory to avoid deleting previous test output.
      * The POLLER_FILE exists in the tests/hap directory.

"""

POLLER_FILE = "wfpc2_2as_01_input.out"


@pytest.fixture(scope="module")
def read_csv_for_filenames():
    # Read the CSV poller file residing in the tests directory to extract the individual visit FLT filenames
    path = os.path.join(os.path.dirname(__file__), POLLER_FILE)
    table = ascii.read(path, format="no_header")
    filename_column = table.colnames[0]
    filenames = list(table[filename_column])
    print("\nread_csv_for_filenames. Filesnames from poller: {}".format(filenames))

    return filenames


@pytest.fixture(scope="module")
def gather_data_for_processing(read_csv_for_filenames, tmp_path_factory):
    # create working directory specified for the test
    curdir = tmp_path_factory.mktemp(os.path.basename(__file__))
    os.chdir(curdir)

    # Establish FLT lists and obtain the requested data
    flt_flag = ""
    # In order to obtain individual FLT images from MAST (if the files do not reside on disk) which
    # may be part of an ASN, use only IPPPSS with a wildcard.  The unwanted images have to be removed
    # after-the-fact.
    for fn in read_csv_for_filenames:
        if fn.lower().endswith("flt.fits") and flt_flag == "":
            flt_flag = fn[0:6] + "*"

    # Get test data through astroquery - only retrieve the pipeline processed FLT files
    # (e.g., j*_flt.fits) as necessary. The logic here and the above for loop is an attempt to
    # avoid downloading too many images which are not needed for processing.
    fltfiles = []
    if flt_flag:
        fltfiles = aqutils.retrieve_observation(flt_flag, suffix=["FLT"], product_type="pipeline")

    print("\ngather_data_for_processing. Gathered data: {}".format(fltfiles))

    return fltfiles


@pytest.fixture(scope="module")
def gather_output_data(construct_manifest_filename):
    # Determine the filenames of all the output files from the manifest
    print(f"\nManifest Filename: {construct_manifest_filename}")
    files = []
    with open(construct_manifest_filename, 'r') as fout:
        for line in fout.readlines():
            files.append(line.rstrip('\n'))
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
    path = os.path.join(os.path.dirname(__file__), POLLER_FILE)
    try:
        status = runsinglehap.perform(path)

    # Catch anything that happens and report it.  This is meant to catch unexpected errors and
    # generate sufficient output exception information so algorithmic problems can be addressed.
    except Exception as except_details:
        print(except_details)
        pytest.fail("\nsvm_setup. Exception Visit: {}\n", path)

    current_dt = datetime.datetime.now()
    print(str(current_dt))


# TESTS

def test_svm_manifest_name(construct_manifest_filename):
    # Construct the manifest filename from the header of an input file in the list and check it exists.
    path = Path(construct_manifest_filename)
    print("\ntest_svm_manifest. Filename: {}".format(path))

    # Ensure the manifest file uses the proper naming convention
    assert (path.is_file())


def test_svm_empty_cats(gather_output_data):
    # Check the output catalogs should contain > 0 measured sources
    cat_files = [files for files in gather_output_data if files.lower().endswith("-cat.ecsv")]

    valid_tables = {}
    for cat in cat_files:
        table_length = len(ascii.read(cat, format="ecsv"))
        print("\ntest_svm_cat_sources. Number of sources in catalog {} is {}.".format(cat, table_length))
        valid_tables[cat] = table_length > 0
    bad_tables = [cat for cat in cat_files if not valid_tables[cat]]
    assert len(bad_tables) == 0, f"Catalog file(s) {bad_tables} is/are unexpectedly empty"
