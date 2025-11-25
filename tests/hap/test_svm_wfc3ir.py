""" This module tests full pipeline SVM processing as a demonstration template.

"""
import pytest
import numpy as np

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

WCS_SUB_NAME = "FIT_SVM"
POLLER_FILE = "wfc3_ir_ib6807_input.out"

# Gather all expected values used for determining pass/fail criteria here
expected_point_sources = {'ir': 316}
expected_seg_sources = {'ir': 299}


@pytest.fixture(scope="module")
def read_csv_for_filenames():
    return svm_utils.load_poller_filenames(POLLER_FILE)


@pytest.fixture(scope="module")
def gather_data_for_processing(read_csv_for_filenames, tmp_path_factory):
    svm_utils.change_to_temp_working_dir(tmp_path_factory, __file__)
    return svm_utils.retrieve_data_for_processing(read_csv_for_filenames)


@pytest.fixture(scope="module")
def gather_output_data(construct_manifest_filename):
    return svm_utils.read_manifest(construct_manifest_filename)


@pytest.fixture(scope="module")
def construct_manifest_filename(read_csv_for_filenames):
    return svm_utils.build_manifest_name(read_csv_for_filenames[0])


@pytest.fixture(scope="module", autouse=True)
def svm_setup(gather_data_for_processing):
    svm_utils.run_svm_pipeline(POLLER_FILE, runsinglehap.perform)


# TESTS

def test_svm_manifest_name(construct_manifest_filename):
    # Construct the manifest filename from the header of an input file in the list and check it exists.
    path = Path(construct_manifest_filename)
    print("\ntest_svm_manifest. Filename: {}".format(path))

    # Ensure the manifest file uses the proper naming convention
    assert (path.is_file())


def test_svm_wcs(gather_output_data):
    # Check the output primary WCSNAME includes FIT_SVM as part of the string value
    tdp_files = [files for files in gather_output_data if
                 files.lower().find("total") > -1 and files.lower().endswith(".fits")]

    for tdp in tdp_files:
        wcsname = fits.getval(tdp, "WCSNAME", ext=1).upper()
        print("\ntest_svm_wcs.  WCSNAME: {} Output file: {}".format(wcsname, tdp))
        assert WCS_SUB_NAME in wcsname, f"WCSNAME is not as expected for file {tdp}."


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


@pytest.mark.skip(reason="Skipping test as missing the 'science commits' for this release - need for RC.")
def test_svm_point_cats(gather_output_data):
    # Check that the point catalogs have the expected number of sources
    cat_files = [files for files in gather_output_data if files.lower().endswith("point-cat.ecsv")]

    num_sources = {cat:len(ascii.read(cat, format="ecsv")) for cat in cat_files}
    valid_cats = {}
    for cat in expected_point_sources.keys():
        for file in cat_files:
            if cat in file and "total" in file:
                valid_cats[cat] = (np.isclose(num_sources[file], expected_point_sources[cat], rtol=0.25), num_sources[file])
                break
    bad_cats = [cat for cat in valid_cats if not valid_cats[cat][0]]
    assert len(bad_cats) == 0,  f"Point Catalog(s) {bad_cats} had {valid_cats} sources, expected {expected_point_sources}"


@pytest.mark.skip(reason="Skipping test as missing the 'science commits' for this release - need for RC.")
def test_svm_segment_cats(gather_output_data):
    # Check that the point catalogs have the expected number of sources
    cat_files = [files for files in gather_output_data if files.lower().endswith("segment-cat.ecsv")]

    num_sources = {cat: len(ascii.read(cat, format="ecsv")) for cat in cat_files}
    valid_cats = {}
    for cat in expected_seg_sources.keys():
        for file in cat_files:
            if cat in file and "total" in file:
                valid_cats[cat] = (np.isclose(num_sources[file], expected_seg_sources[cat], rtol=0.25), num_sources[file])
                break
    bad_cats = [cat for cat in valid_cats if not valid_cats[cat][0]]
    assert len(bad_cats) == 0, f"Segment Catalog(s) {bad_cats} had {valid_cats} sources, expected {expected_seg_sources}"
