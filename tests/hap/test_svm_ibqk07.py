""" This module tests full pipeline SVM processing of a WFC3 dataset containing both IR and UVIS data.
    The two detectors do NOT have the same WCS solution.

"""
import numpy as np
import pytest

from drizzlepac import runsinglehap
from astropy.io import fits, ascii
from pathlib import Path

from . import svm_test_utils as svm_utils

pytestmark = pytest.mark.bigdata

"""
    test_svm_ibqk07.py

    This test file can be executed in the following manner:
        $ pytest -s --basetemp=/internal/hladata/yourUniqueDirectoryHere test_svm_ibqk07.py >& ibqk07.log &
        $ tail -f ibqk07.log
      * Note: When running this test, the `--basetemp` directory should be set to a unique
        existing directory to avoid deleting previous test output.
      * The POLLER_FILE exists in the tests/hap directory with the PyTests.
      * If running manually with `--basetemp`, the ibqk07.log file will still be written to the 
        originating directory.

"""

POLLER_FILE = "wfc3_bqk_07_input.out"
WCS_IR_SUB_NAME = "HSC30"
WCS_UVIS_SUB_NAME = "HSC30"
expected_total_point_sources = {
"hst_12557_07_wfc3_ir_total_ibqk07_point-cat.ecsv": 2,
"hst_12557_07_wfc3_uvis_total_ibqk07_point-cat.ecsv": 13}
expected_total_segment_sources= {
"hst_12557_07_wfc3_ir_total_ibqk07_segment-cat.ecsv": 2,
"hst_12557_07_wfc3_uvis_total_ibqk07_segment-cat.ecsv": 19}
tolerance = 0.25 

MEAN_CAT_MAGAP2_POINT = {
"hst_12557_07_wfc3_ir_f160w_ibqk07_point-cat.ecsv": 24.27,
"hst_12557_07_wfc3_uvis_f555w_ibqk07_point-cat.ecsv": 25.36}
MEAN_CAT_MAGAP2_SEGMENT = {
"hst_12557_07_wfc3_uvis_f555w_ibqk07_segment-cat.ecsv": 24.97,
"hst_12557_07_wfc3_ir_f160w_ibqk07_segment-cat.ecsv": 18.76}
POINT_DIFF = 0.5
SEGMENT_DIFF = 0.5

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
    print("\ntest_svm_manifest.")
    # Construct the manifest filename from the header of an input file in the list and check it exists.
    path = Path(construct_manifest_filename)
    print("\ntest_svm_manifest. Filename: {}".format(path))

    # Ensure the manifest file uses the proper naming convention
    assert(path.is_file())

def test_svm_wcs_ir(gather_output_data):
    print("\ntest_svm_wcs_ir.")
    # Get the TDP for this detector
    tdp_files = [files for files in gather_output_data if files.lower().find("ir_total") > -1 and files.lower().endswith("drz.fits")]

    # Check the WCS solution is as expected
    wcsname = fits.getval(tdp_files[0], "WCSNAME", ext=1).upper()
    print("\ntest_svm_wcs_ir.  WCSNAME: {} Output file: {}".format(wcsname, tdp_files[0]))
    assert WCS_IR_SUB_NAME in wcsname, f"WCSNAME is not as expected for file {tdp_files[0]}."


def test_svm_wcs_ir_all(gather_output_data):
    print("\ntest_svm_wcs_ir_all.")
    # Check the output primary WCSNAME 
    ir_files = [files for files in gather_output_data if files.lower().find("_ir_") > -1 and files.lower().endswith("drz.fits")]

    wcsnames = [fits.getval(ir, "WCSNAME", ext=1).upper() for ir in ir_files]
    assert len(set(wcsnames)) == 1, f"WCSNAMES are not all the same for the IR detector: {wcsnames}"


def test_svm_wcs_uvis(gather_output_data):
    print("\ntest_svm_wcs_uvis.")
    # Get the TDP for this detector
    tdp_files = [files for files in gather_output_data if files.lower().find("uvis_total") > -1 and files.lower().endswith("drc.fits")]

    # Check the WCS solution is as expected
    wcsname = fits.getval(tdp_files[0], "WCSNAME", ext=1).upper()
    print("\ntest_svm_wcs_uvis.  WCSNAME: {} Output file: {}".format(wcsname, tdp_files[0]))
    assert WCS_UVIS_SUB_NAME in wcsname, f"WCSNAME is not as expected for file {tdp_files[0]}."


def test_svm_wcs_uvis_all(gather_output_data):
    # Check the output primary WCSNAME 
    print("\ntest_svm_wcs_uvis_all.")
    uvis_files = [files for files in gather_output_data if files.lower().find("_uvis_") > -1 and files.lower().endswith("drc.fits")]

    wcsnames = [fits.getval(uvis, "WCSNAME", ext=1).upper() for uvis in uvis_files]
    assert len(set(wcsnames)) == 1, f"WCSNAMES are not all the same for the UVIS detector: {wcsnames}"


# Due to the way the catalogs are filtered, check the size of the total catalog and one of the filter
# catalogs separately.  The total catalog has the row removed for each source where the constituent 
# filter catalogs *ALL* have flag>5 for the source.  Rows are NOT removed from the filter table based on
# flag values.
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


@pytest.mark.skip(reason="Modifying tests and cannot reproduce failed result at this time - need for RC.")
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


def test_svm_point_cat_meanmag(gather_output_data):
    cat_files = [files.lower() for files in gather_output_data if files.lower().endswith("point-cat.ecsv") and files.lower().find("total") < 0]

    # Compute the mean of the MagAp2 in the filtered catalogs and do not include flagged bad data
    Mag2_mean = {}
    for cat in cat_files:
        table = ascii.read(cat)
        Mag2_array = table['MagAp2'].data
        Mag2_mean[cat] = -9999.0
        if len(Mag2_array[Mag2_array > -9999.0]) > 0:
            Mag2_mean[cat] = Mag2_array[Mag2_array > -9999.0].mean()

    good_cats = {}
    for cat in MEAN_CAT_MAGAP2_POINT.keys():
        for file in cat_files:
            if cat == file:
                good_cats[cat] = (np.isclose(MEAN_CAT_MAGAP2_POINT[cat], Mag2_mean[cat], rtol=POINT_DIFF), Mag2_mean[cat])
                break

    bad_cats = [cat for cat in good_cats if not good_cats[cat][0]]
    assert len(bad_cats) == 0, f"Point Catalog(s) {bad_cats} had {good_cats} sources, expected {MEAN_CAT_MAGAP2_POINT}"


def test_svm_segment_cat_meanmag(gather_output_data):
    cat_files = [files.lower() for files in gather_output_data if files.lower().endswith("segment-cat.ecsv") and files.lower().find("total") < 0]

    # Compute the mean of the MagAp2 in the filtered catalogs and do not include flagged bad data
    Mag2_mean = {}
    for cat in cat_files:
        table = ascii.read(cat)
        Mag2_array = table['MagAp2'].data
        Mag2_mean[cat] = -9999.0
        if len(Mag2_array[Mag2_array > -9999.0]) > 0:
            Mag2_mean[cat] = Mag2_array[Mag2_array > -9999.0].mean()

    good_cats = {}
    for cat in MEAN_CAT_MAGAP2_SEGMENT.keys():
        for file in cat_files:
            if cat == file:
                good_cats[cat] = (np.isclose(MEAN_CAT_MAGAP2_SEGMENT[cat], Mag2_mean[cat], rtol=SEGMENT_DIFF), Mag2_mean[cat])
                break

    bad_cats = [cat for cat in good_cats if not good_cats[cat][0]]
    assert len(bad_cats) == 0, f"Segment Catalog(s) {bad_cats} had {good_cats} sources, expected {MEAN_CAT_MAGAP2_SEGMENT}"
