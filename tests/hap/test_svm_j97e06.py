""" This module tests full pipeline SVM processing of ACS WFC Grism and Direct image data.

"""
import datetime
import glob
import os
import math
import numpy as np
import pytest

from drizzlepac.haputils import astroquery_utils as aqutils
from drizzlepac import runsinglehap
from astropy.io import fits, ascii
from astropy.table import Table
from pathlib import Path

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
EXPECTED_SEG_SOURCES = {"wfc": 2}
MEAN_CAT_MAGAP1_POINT = {
"hst_10374_06_acs_wfc_f814w_j97e06_point-cat.ecsv": 18.79,
"hst_10374_06_acs_wfc_f625w_j97e06_point-cat.ecsv": 18.86,
"hst_10374_06_acs_wfc_f775w_j97e06_point-cat.ecsv": 18.78,
"hst_10374_06_acs_wfc_f555w_j97e06_point-cat.ecsv": 19.05,
"hst_10374_06_acs_wfc_f850lp_j97e06_point-cat.ecsv": 18.84,
"hst_10374_06_acs_wfc_f606w_j97e06_point-cat.ecsv": 18.97}
MEAN_CAT_MAGAP1_SEGMENT = {
"hst_10374_06_acs_wfc_f850lp_j97e06_segment-cat.ecsv": 20.56,
"hst_10374_06_acs_wfc_f555w_j97e06_segment-cat.ecsv": 21.73,
"hst_10374_06_acs_wfc_f775w_j97e06_segment-cat.ecsv": 20.71,
"hst_10374_06_acs_wfc_f814w_j97e06_segment-cat.ecsv": 20.73,
"hst_10374_06_acs_wfc_f625w_j97e06_segment-cat.ecsv": 21.08,
"hst_10374_06_acs_wfc_f606w_j97e06_segment-cat.ecsv": 21.17}
POINT_DIFF = 0.5
SEGMENT_DIFF = 0.5

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
    files_to_process= set(read_csv_for_filenames).intersection(set(flcfiles))

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
    files = []
    with open(construct_manifest_filename, 'r') as fout:
        for line in fout.readlines():
            files.append(line.rstrip("\n"))
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
    assert(path.is_file())


def test_svm_wcs(gather_output_data):
    # Check the output primary WCSNAME includes IDC_* of the string value
    tdp_files = [files for files in gather_output_data if files.lower().find("total") > -1 and files.lower().endswith(".fits")]

    for tdp in tdp_files:
        wcsname = fits.getval(tdp, "WCSNAME", ext=1).upper()
        print("\ntest_svm_wcs.  WCSNAME: {} Output file: {}".format(wcsname, tdp))
        assert WCS_SUB_NAME in wcsname, f"WCSNAME is not as expected for file {tdp}."


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
