""" This module tests full pipeline SVM processing of ACS WFC Grism and Direct image data.

"""
import numpy as np
import pytest

from ci_watson.plugin import _jail
from drizzlepac import runsinglehap
from astropy.io import fits, ascii
from pathlib import Path

from . import svm_test_utils as svm_utils

pytestmark = pytest.mark.bigdata

"""
    test_svm_j97e06.py

    This test file can be executed in the following manner:
        $ pytest -s --basetemp=/internal/hladata/yourUniqueDirectoryHere test_svm_j97e06.py >& j97e06.log &
        $ tail -f j97e06.log
      * Note: When running this test, the `--basetemp` directory should be set to a unique
        existing directory to avoid deleting previous test output.
      * The POLLER_FILE exists in the tests/hap directory with the PyTests.
      * If running manually with `--basetemp`, the j97e06.log file will still be written to the 
        originating directory.

"""

# Expectation values used directly or indirectly for the test assert statements
WCS_SUB_NAME = "IDC_4BB1536OJ"
POLLER_FILE = "acs_97e_06_input.out"
EXPECTED_POINT_SOURCES = {"wfc": 2}
EXPECTED_SEG_SOURCES = {"wfc": 6}
MEAN_CAT_MAGAP1_POINT = {
"hst_10374_06_acs_wfc_f814w_j97e06_point-cat.ecsv": 17.55,
"hst_10374_06_acs_wfc_f625w_j97e06_point-cat.ecsv": 17.86,
"hst_10374_06_acs_wfc_f775w_j97e06_point-cat.ecsv": 17.53,
"hst_10374_06_acs_wfc_f555w_j97e06_point-cat.ecsv": 18.14,
"hst_10374_06_acs_wfc_f850lp_j97e06_point-cat.ecsv": 17.75,
"hst_10374_06_acs_wfc_f606w_j97e06_point-cat.ecsv": 17.94}
MEAN_CAT_MAGAP1_SEGMENT = {
"hst_10374_06_acs_wfc_f850lp_j97e06_segment-cat.ecsv": 20.36,
"hst_10374_06_acs_wfc_f555w_j97e06_segment-cat.ecsv": 21.80,
"hst_10374_06_acs_wfc_f775w_j97e06_segment-cat.ecsv": 20.48,
"hst_10374_06_acs_wfc_f814w_j97e06_segment-cat.ecsv": 20.38,
"hst_10374_06_acs_wfc_f625w_j97e06_segment-cat.ecsv": 21.26,
"hst_10374_06_acs_wfc_f606w_j97e06_segment-cat.ecsv": 21.28}
POINT_DIFF = 0.5
SEGMENT_DIFF = 0.5

@pytest.fixture(scope="module")
def read_csv_for_filenames():
    return svm_utils.load_poller_filenames(POLLER_FILE)


@pytest.fixture(scope="function")
def gather_data_for_processing(_jail, read_csv_for_filenames, pytestconfig):
    """Retrieve observation files inside the temporary jail directory."""
    return svm_utils.retrieve_data_for_processing(
            read_csv_for_filenames,
            pytestconfig=pytestconfig,
        )


@pytest.fixture(scope="function")
def svm_run(gather_data_for_processing):
    """Execute SVM pipeline once inputs are staged."""
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
    assert(path.is_file())


def test_svm_wcs(gather_output_data):
    # Check the output primary WCSNAME includes IDC_* of the string value
    tdp_files = [files for files in gather_output_data if files.lower().find("total") > -1 and files.lower().endswith(".fits")]

    for tdp in tdp_files:
        wcsname = fits.getval(tdp, "WCSNAME", ext=1).upper()
        print("\ntest_svm_wcs.  WCSNAME: {} Output file: {}".format(wcsname, tdp))
        assert WCS_SUB_NAME in wcsname, f"WCSNAME is not as expected for file {tdp}."

@pytest.mark.skip(reason="Need to update logic as changes in the small number of catalog sources is too sensitive.")
def test_svm_point_cat_numsources(gather_output_data):
   # Check that the point catalogs have the expected number of sources
    cat_files = [files for files in gather_output_data if files.lower().endswith("point-cat.ecsv")]

    num_sources = {cat:len(ascii.read(cat, format="ecsv")) for cat in cat_files}
    valid_cats = {}
    for cat in EXPECTED_POINT_SOURCES.keys():
        for file in cat_files:
            if cat in file and "total" in file:
                valid_cats[cat] = (np.isclose(num_sources[file], EXPECTED_POINT_SOURCES[cat], atol=1), num_sources[file])
                break
    bad_cats = [cat for cat in valid_cats if not valid_cats[cat][0]]
    assert len(bad_cats) == 0,  f"Point Catalog(s) {bad_cats} had {valid_cats} sources, expected {EXPECTED_POINT_SOURCES}"

@pytest.mark.skip(reason="Need to update logic as changes in the small number of catalog sources is too sensitive.")
def test_svm_segment_cat_numsources(gather_output_data):
   # Check that the point catalogs have the expected number of sources
    cat_files = [files for files in gather_output_data if files.lower().endswith("segment-cat.ecsv")]

    num_sources = {cat: len(ascii.read(cat, format="ecsv")) for cat in cat_files}
    valid_cats = {}
    for cat in EXPECTED_SEG_SOURCES.keys():
        for file in cat_files:
            if cat in file and "total" in file:
                valid_cats[cat] = (np.isclose(num_sources[file], EXPECTED_SEG_SOURCES[cat], atol=1), num_sources[file])
                break
    bad_cats = [cat for cat in valid_cats if not valid_cats[cat][0]]
    assert len(bad_cats) == 0, f"Segment Catalog(s) {bad_cats} had {valid_cats} sources, expected {EXPECTED_SEG_SOURCES}"


def test_svm_point_cat_meanmag(gather_output_data):
    cat_files = [files.lower() for files in gather_output_data if files.lower().endswith("point-cat.ecsv") and files.lower().find("total") < 0]

    # Compute the mean of the MagAp1 in the filtered catalogs and do not include flagged bad data
    Mag1_mean = {}
    for cat in cat_files:
        table = ascii.read(cat)
        Mag1_array = table['MagAp1'].data
        Mag1_mean[cat] = -9999.0
        if len(Mag1_array[Mag1_array > -9999.0]) > 0:
            Mag1_mean[cat] = Mag1_array[Mag1_array > -9999.0].mean()

    good_cats = {}
    for cat in MEAN_CAT_MAGAP1_POINT.keys():
        for file in cat_files:
            if cat == file:
                good_cats[cat] = (np.isclose(MEAN_CAT_MAGAP1_POINT[cat], Mag1_mean[cat], rtol=POINT_DIFF), Mag1_mean[cat])
                break

    bad_cats = [cat for cat in good_cats if not good_cats[cat][0]]
    assert len(bad_cats) == 0, f"Point Catalog(s) {bad_cats} had {good_cats} sources, expected {MEAN_CAT_MAGAP1_POINT}"


def test_svm_segment_cat_meanmag(gather_output_data):
    cat_files = [files.lower() for files in gather_output_data if files.lower().endswith("segment-cat.ecsv") and files.lower().find("total") < 0]

    # Compute the mean of the MagAp1 in the filtered catalogs and do not include flagged bad data
    Mag1_mean = {}
    for cat in cat_files:
        table = ascii.read(cat)
        Mag1_array = table['MagAp1'].data
        Mag1_mean[cat] = -9999.0
        if len(Mag1_array[Mag1_array > -9999.0]) > 0:
            Mag1_mean[cat] = Mag1_array[Mag1_array > -9999.0].mean()

    good_cats = {}
    for cat in MEAN_CAT_MAGAP1_SEGMENT.keys():
        for file in cat_files:
            if cat == file:
                good_cats[cat] = (np.isclose(MEAN_CAT_MAGAP1_SEGMENT[cat], Mag1_mean[cat], rtol=SEGMENT_DIFF), Mag1_mean[cat])
                break

    bad_cats = [cat for cat in good_cats if not good_cats[cat][0]]
    assert len(bad_cats) == 0, f"Segment Catalog(s) {bad_cats} had {good_cats} sources, expected {MEAN_CAT_MAGAP1_SEGMENT}"
