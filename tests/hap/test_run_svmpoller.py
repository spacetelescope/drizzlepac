""" This module is the high-level wrapper for running a test on a list of input
    poller files (e.g., ib4606.out) in order to process the specified data through
    the Single Visit Mosaic pipeline.
    """
import datetime
import glob
import os
import shutil
import sys
import traceback
import numpy as np
import pytest

from astropy.io import ascii
from drizzlepac import runsinglehap
from astropy.table import Table
from drizzlepac.haputils import astroquery_utils as aqutils


def pytest_generate_tests(metafunc):
    """Get the command line options."""

    # The input svm_list option should be given a fully qualified filename
    svm_list = metafunc.config.option.svm_list

    # The list containing the poller filenames must exist in the current working directory
    if not os.path.exists(svm_list):
        # If not, copy the default file from the installation directory
        # Find where this module has been installed.
        install_dir = os.path.dirname(__file__)
        default_file = os.path.join(install_dir, svm_list)
        if os.path.exists(default_file):
            # Copy the file
            shutil.copy2(default_file, os.getcwd())

    # Read the svm_list to get the poller filenames
    table = Table.read(svm_list, format="ascii.fast_no_header")
    poller_file_list = table["col1"].tolist()

    print("Input file: {}".format(svm_list))
    print("List of poller files: {}".format(poller_file_list))
    metafunc.parametrize('dataset', poller_file_list)


@pytest.mark.bigdata
@pytest.mark.slow
@pytest.mark.unit
def test_run_svmpoller(tmpdir, dataset):
    """ Tests to read a series of poller files and process the contents of each as Single Visit Mosaic

        Characteristics of these tests:

        Success Criteria:
            The SVM processing returns a value of 0: Success or 1: Failure

        The input svm_list file is a list of poller filenames, one filename per line.
        Each poller file must be obtained from a specified directory and read to obtain the
        names of the data files which need to be processed.

        This test file can be executed in the following manner:
            $ pytest -n # -s --basetemp=/internal/hladata/yourUniqueDirectoryHere --bigdata --slow
              --svm_list svm_input.lst test_run_svmpoller.py >& test_svmpoller_output.txt &
            $ tail -f test_svmpoller_output.txt
          * The `-n #` option can be used to run tests in parallel if `pytest-xdist` has
            been installed where `#` is the number of cpus to use. THIS IS NOT ADVISED FOR USE.
          * Note: When running this test, the `--basetemp` directory should be set to a unique
            existing directory to avoid deleting previous test output.
          * A default master list, svm_input.lst, exists in the tests/hla directory and contains 3 datasets.
            This specific list may NOT the list you want to use, but it allows you to see what this file
            should contain.  Please note the PyTests should be kept to runtimes which are not
            excessive.

    """
    print("TEST_RUN_SVMPOLLER. Dataset: ", dataset)

    current_dt = datetime.datetime.now()
    print(str(current_dt))

    subdir = ""
    prevdir = os.getcwd()

    # create working directory specified for the test
    if not tmpdir.ensure(subdir, dir=True):
        curdir = tmpdir.mkdir(subdir).strpath
    else:
        curdir = tmpdir.join(subdir).strpath
    os.chdir(curdir)

    return_value = 1

    try:

        # Read the CSV poller file residing in the tests directory to extract the individual visit FLT/FLC filenames
        path = os.path.join(os.path.dirname(__file__), dataset)
        table = ascii.read(path, format="no_header")
        filename_column = table.colnames[0]
        filenames = list(table[filename_column])
        print("\nread_csv_for_filenames. Filesnames from poller: {}".format(filenames))

        # Establish FLC/FLT lists and obtain the requested data
        flc_flag = ""
        flt_flag = ""
        # In order to obtain individual FLC or FLT images from MAST (if the files are not reside on disk) which
        # may be part of an ASN, use only IPPPSS with a wildcard.  The unwanted images have to be removed
        # after-the-fact.
        for fn in filenames:
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
        files_to_process= set(filenames).intersection(set(flcfiles))

        # Identify unwanted files from the download list and remove from disk
        files_to_remove = set(filenames).symmetric_difference(set(flcfiles))

        try:
            for ftr in files_to_remove:
               os.remove(ftr)
        except Exception as x_cept:
            print("")
            print("Exception encountered: {}.".format(x_cept))
            print("The file {} could not be deleted from disk. ".format(ftr))
            print("Remove files which are not used for processing from disk manually.")

        # Run the SVM processing
        path = os.path.join(os.path.dirname(__file__), dataset)

        return_value = runsinglehap.perform(path)

    # Catch anything that happens and report it.  This is meant to catch unexpected errors and
    # generate sufficient output exception information so algorithmic problems can be addressed.
    except Exception as except_details:
        traceback.print_exc()
        pytest.fail("TEST_RUN_SVMPOLLER. Exception Dataset: {}\n", dataset)
        return_value = 1

    assert return_value == 0

    current_dt = datetime.datetime.now()
    print(str(current_dt))

    # Return to original directory
    os.chdir(prevdir)
