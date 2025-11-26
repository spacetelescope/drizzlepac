""" This module tests full pipeline SVM processing for an ACS WFC, full-frame, one filter dataset.

"""
import pytest

from drizzlepac import runsinglehap
from astropy.io import fits, ascii
from pathlib import Path

from . import svm_test_utils as svm_utils

pytestmark = pytest.mark.bigdata

"""
    test_svm_je281u.py

    This test file can be executed in the following manner:
        $ pytest -s --basetemp=/internal/hladata/yourUniqueDirectoryHere test_svm_je281u.py >& je281u.log &
        $ tail -f je281u.log
      * Note: When running this test, the `--basetemp` directory should be set to a unique
        existing directory to avoid deleting previous test output.
      * The POLLER_FILE exists in the tests/hap directory.
      * If running manually with `--basetemp`, the je281u.log file will still be written to the 
        originating directory.

"""

WCS_SUB_NAME = "FIT_SVM_GAIA"
POLLER_FILE = "acs_e28_1u_input.out"

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
    assert(path.is_file())


def test_svm_wcs(gather_output_data):
    # Check the output primary WCSNAME includes FIT_SVM_GAIA as part of the string value
    tdp_files = [files for files in gather_output_data if files.lower().find("total") > -1 and files.lower().endswith(".fits")]

    for tdp in tdp_files:
        wcsname = fits.getval(tdp, "WCSNAME", ext=1).upper()
        print("\ntest_svm_wcs.  WCSNAME: {} Output file: {}".format(wcsname, tdp))
        assert WCS_SUB_NAME in wcsname, f"WCSNAME is not as expected for file {tdp}."


def test_svm_point_cat_cols(gather_output_data):
    # Check the total catalog product does not contain any unexpected, non-filter dependent columns
    tdp_files = [files for files in gather_output_data if files.lower().find("total") > -1 and files.lower().endswith("point-cat.ecsv")]

    ref_strings = ["ID", "Center", "RA", "DEC"]
    for tdp in tdp_files:
        table = ascii.read(tdp, format="ecsv")
        sub_columns = []
        for c in table.colnames:
            if "_f" not in c:
                strip_c = c.lstrip("XY-")
                sub_columns.append(strip_c)

        for c in sub_columns:
            if c not in ref_strings:
                assert 0, f"Unexpected column, {c}, found in Total Point Catalog file"


def test_svm_segment_cat_cols(gather_output_data):
    # Check the total catalog product does not contain any unexpected, non-filter dependent columns
    tdp_files = [files for files in gather_output_data if files.lower().find("total") > -1 and files.lower().endswith("segment-cat.ecsv")]

    ref_strings = ["ID", "Centroid", "RA", "DEC"]
    for tdp in tdp_files:
        table = ascii.read(tdp, format="ecsv")
        sub_columns = []
        for c in table.colnames:
            if "_f" not in c:
                strip_c = c.lstrip("XY-")
                sub_columns.append(strip_c)

        for c in sub_columns:
            if c not in ref_strings:
                assert 0, f"Unexpected column, {c}, found in Total Segment Catalog file"


def test_svm_cat_sources(gather_output_data):
    # Check the output catalogs should contain > 0 measured sources
    cat_files = [files for files in gather_output_data if files.lower().endswith("-cat.ecsv")]

    for cat in cat_files:
        table_length = len(ascii.read(cat, format="ecsv"))
        print("\ntest_svm_cat_sources. Number of sources in catalog {} is {}.".format(cat, table_length))
        assert table_length > 0, f"Catalog file {cat} is unexpectedly empty"
