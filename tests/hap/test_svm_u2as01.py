""" This module tests full pipeline SVM processing as a demonstration template.

"""
import pytest

from drizzlepac import runsinglehap
from astropy.io import fits, ascii
from pathlib import Path

from . import svm_test_utils as svm_utils

pytestmark = pytest.mark.bigdata

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
    return svm_utils.load_poller_filenames(POLLER_FILE)


@pytest.fixture(scope="module")
def gather_data_for_processing(read_csv_for_filenames, tmp_path_factory):
    svm_utils.change_to_temp_working_dir(tmp_path_factory, __file__)
    return svm_utils.retrieve_data_for_processing(read_csv_for_filenames, suffixes=("FLT",))


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
