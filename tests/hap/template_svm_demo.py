""" This module tests full pipeline SVM processing as a demonstration template.

"""
import datetime
import glob
import os
import pytest

from drizzlepac import runsinglehap
from astropy.io import fits, ascii
from pathlib import Path
from ci_watson.artifactory_helpers import get_bigdata

"""
    template_svm_demo.py

    This test file can be executed in the following manner:
        $ pytest -s --basetemp=/internal/hladata/yourUniqueDirectoryHere template_svm_demo.py >& template_svm_demo.log &
        $ tail -f template_svm_demo.log
      * Note: When running this test, the `--basetemp` directory should be set to a unique
        existing directory to avoid deleting previous test output.
      * The POLLER_FILE exists in the tests/hap directory.
      * If running manually with `--basetemp`, the template_svm_demo.log file will still be written to the 
        originating directory.

"""

WCS_SUB_NAME = "FIT_SVM_GAIA"
POLLER_FILE = "acs_e28_1u_input.out"

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
def gather_data_for_processing(read_csv_for_filenames, tmp_path_factory):
    # Create working directory specified for the test
    curdir = tmp_path_factory.mktemp(os.path.basename(__file__)) 
    os.chdir(curdir)

    # Get the data from Artifactory
    inputs = [os.path.basename(get_bigdata('drizzlepac', 'dev', 'acs', 'input', i))
              for i in read_csv_for_filenames]

    files_to_process = read_csv_for_filenames
    print("\ngather_data_for_processing. Gathered data: {}".format(files_to_process))

    return files_to_process


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


# TESTS - Avoid doing an assert in a loop as the test could exit before all the data has processed

def test_svm_manifest_name(construct_manifest_filename):
    # Construct the manifest filename from the header of an input file in the list and check it exists.
    path = Path(construct_manifest_filename)
    print("\ntest_svm_manifest. Filename: {}".format(path))

    # Ensure the manifest file uses the proper naming convention
    assert(path.is_file())


def test_svm_wcs(gather_output_data):
    # Check the output primary WCSNAME includes FIT_SVM_GAIA as part of the string value
    tdp_files = [files for files in gather_output_data if files.lower().find("total") > -1 and files.lower().endswith(".fits")]

    # This check is for all total data products which have the same "type" of WCSNAME -
    # in this case a name akin to *-FIT_SVM_GAIA*.
    wcsnames = [fits.getval(tdp, "WCSNAME", ext=1).upper().split('-')[1] for tdp in tdp_files]
    assert len(set(wcsnames)) == 1, f"WCSNAMES are not all the same: {wcsnames}"


# Due to the way the catalogs are filtered, check the size of the total catalog and one of the filter
# catalogs separately.  The total catalog has the row removed for each source where the constituent 
# filter catalogs *ALL* have flag>5.
def test_svm_point_total_cat(gather_output_data):
    # Check the output catalogs should contain the correct number of sources
    tdp_files = [files for files in gather_output_data if files.lower().find("total") > -1 and files.lower().endswith("point-cat.ecsv")]

    valid_tables = {}
    for cat in tdp_files:
        table_length = len(ascii.read(cat, format="ecsv"))
        print("\ntest_svm_point_total_cat. Number of sources in catalog {} is {}.".format(cat, table_length))
        valid_tables[cat] = table_length > 0
    bad_tables = [cat for cat in cat_files if not valid_tables[cat]]
    assert len(bad_tables) == 0, f"Catalog file(s) {bad_tables} is/are unexpectedly empty"


