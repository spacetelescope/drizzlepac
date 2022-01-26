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

# Gather expected values for pass/fail criteria here
expected_total_point_sources = {'hrc': 268, 'sbc': 65}
expected_total_segment_sources = {'hrc': 642, 'sbc': 283}
expected_filter_point_sources = {'hrc': 269, 'sbc': 208}
expected_filter_segment_sources = {'hrc': 642, 'sbc': 399}
tolerance = 0.25 


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
    # create working directory specified for the test
    curdir = tmp_path_factory.mktemp(os.path.basename(__file__))
    os.chdir(curdir)

    # Establish FLC/FLT lists and obtain the requested data
    flc_flag = ""
    flt_flag = ""
    # In order to obtain individual FLC or FLT images from MAST (if the files are not reside on disk) which
    # may be part of an ASN, use only IPPPSS with a wildcard.  The unwanted images have to be removed
    # after-the-fact.
    for fn in read_csv_for_filenames:
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
    files_to_process = set(read_csv_for_filenames).intersection(set(flcfiles))

    # Identify unwanted files from the download list and remove from disk
    files_to_remove = set(read_csv_for_filenames).symmetric_difference(set(flcfiles))
    try:
        for ftr in files_to_remove:
            os.remove(ftr)
    except Exception as x_cept:
        print("")
        print("Exception encountered: {}.".format(x_cept))
        print("The file {} could not be deleted from disk. ".format(ftr))
        print("Remove files which are not used for processing from disk manually.")

    print("\ngather_data_for_processing. Gathered data: {}".format(files_to_process))

    return files_to_process


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


def test_svm_wcs(gather_output_data):
    # Check the output primary WCSNAME includes FIT_SVM_GAIA as part of the string value
    tdp_files = [files for files in gather_output_data if
                 files.lower().find("total") > -1 and files.lower().endswith(".fits")]

    for tdp in tdp_files:
        wcsname = fits.getval(tdp, "WCSNAME", ext=1).upper()
        print("\ntest_svm_wcs.  WCSNAME: {} Output file: {}".format(wcsname, tdp))
        assert WCS_SUB_NAME in wcsname, f"WCSNAME is not as expected for file {tdp}."


def test_svm_samewcs(gather_output_data):
    # Check that products for both detectors are aligned to the same catalog
    # The assumption is that if they are all aligned to the same catalog, they are
    # correctly aligned to each other.
    tdp_files = [files for files in gather_output_data if
                 files.lower().find("total") > -1 and files.lower().endswith(".fits")]

    print(f'TDP_FILES: \n{tdp_files}')

    wcsnames = [fits.getval(tdp, "WCSNAME", ext=1).upper().split('-')[1] for tdp in tdp_files]
    assert len(set(wcsnames)) == 1, f"WCSNAMES are not all the same: {wcsnames}"


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


# Due to the way the catalogs are filtered, check the size of the total catalog and one of the filter
# catalogs separately.  The total catalog has the row removed for each source where the constituent 
# filter catalogs *ALL* have flag>5 for the source.  Rows are NOT removed from the filter table based on
# flag values.
def test_svm_point_total_cat(gather_output_data):
    # Check the output catalogs should contain the correct number of sources -- allows for a broad tolerance
    print("\ntest_svm_point_total_cat.")
    tdp_files = [files for files in gather_output_data if files.lower().find("total") > -1 and files.lower().endswith("point-cat.ecsv")]

    num_sources = {tdp:len(ascii.read(tdp, format="ecsv")) for tdp in tdp_files}
    valid_cats = {}
    for tdp in expected_total_point_sources.keys():
        for file in tdp_files:
            if tdp in file:
                tol_limit = tolerance * expected_total_point_sources[tdp]
                valid_cats[tdp] = (file, np.isclose(expected_total_point_sources[tdp], num_sources[file], atol=tol_limit))
                break
    bad_cats = [cat for cat in valid_cats if not valid_cats[cat][1]]
    assert len(bad_cats) == 0,  f"Total Point Catalog(s) {bad_cats} had {valid_cats} sources, expected {expected_total_point_sources}"


def test_svm_segment_total_cat(gather_output_data):
    # Check the output catalogs should contain the correct number of sources -- allows for a broad tolerance
    print("\ntest_svm_segment_total_cat.")
    tdp_files = [files for files in gather_output_data if files.lower().find("total") > -1 and files.lower().endswith("segment-cat.ecsv")]

    num_sources = {tdp:len(ascii.read(tdp, format="ecsv")) for tdp in tdp_files}
    valid_cats = {}
    for tdp in expected_total_segment_sources.keys():
        for file in tdp_files:
            if tdp in file:
                tol_limit = tolerance * expected_total_segment_sources[tdp]
                valid_cats[tdp] = (file, np.isclose(expected_total_segment_sources[tdp], num_sources[file], atol=tol_limit))
                break
    bad_cats = [cat for cat in valid_cats if not valid_cats[cat][1]]
    assert len(bad_cats) == 0,  f"Total Segment Catalog(s) {bad_cats} had {valid_cats} sources, expected {expected_total_segment_sources}"
