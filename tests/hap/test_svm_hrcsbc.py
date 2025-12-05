""" This module tests full pipeline SVM processing as a demonstration template.

"""
import pytest
import numpy as np

from ci_watson.plugin import _jail
from drizzlepac import runsinglehap
from astropy.io import fits, ascii
from pathlib import Path

from . import svm_test_utils as svm_utils

pytestmark = pytest.mark.bigdata

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
expected_total_segment_sources = {'hrc': 642, 'sbc': 250}
tolerance = 0.25 


@pytest.fixture(scope="module")
def read_csv_for_filenames():
    return svm_utils.load_poller_filenames(POLLER_FILE)

@pytest.fixture(scope="function")
def gather_data_for_processing(_jail, read_csv_for_filenames, pytestconfig):
    """Retrieve required files inside the temporary jail directory."""
    return svm_utils.retrieve_data_for_processing(
        read_csv_for_filenames,
        pytestconfig=pytestconfig,
    )


@pytest.fixture(scope="function")
def svm_run(gather_data_for_processing):
    """Execute the SVM pipeline after inputs are staged in the jail."""
    svm_utils.run_svm_pipeline(POLLER_FILE, runsinglehap.perform)
    return gather_data_for_processing


@pytest.fixture(scope="function")
def construct_manifest_filename(read_csv_for_filenames, svm_run):
    return svm_utils.build_manifest_name(read_csv_for_filenames[0])


@pytest.fixture(scope="function")
def gather_output_data(construct_manifest_filename, svm_run):
    return svm_utils.read_manifest(construct_manifest_filename)


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
# flag values. NOTE: Filtered catalogs are actually not checked by these tests.
@pytest.mark.skip(reason="Modifying tests and cannot reproduce failed result at this time - need for RC.")
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


@pytest.mark.skip(reason="Modifying tests and cannot reproduce failed result at this time. - need for RC")
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
