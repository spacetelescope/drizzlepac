""" This module is the high-level wrapper for running a test on a list of input
    poller files (e.g., ib4606.out) in order to process the specified data through
    the Single Visit Mosaic pipeline.
    """
import datetime
import traceback
import os
import shutil
import sys
import numpy as np
from astropy.table import Table
import pytest
import logging
from importlib import reload
import glob

from stsci.tools import logutil
from astropy.table import Table

from drizzlepac import runsinglehap


log = logutil.create_logger('test_run_svmpoller', level=logutil.logging.INFO, stream=sys.stdout)


def pytest_generate_tests(metafunc):
    """Get the command line options."""

    # The input master_list option should be given a fully qualified filename
    master_list = metafunc.config.option.master_list

    # The list containing the poller filenames must exist in the current working directory
    if not os.path.exists(master_list):
        # If not, copy the default file from the installation directory
        # Find where this module has been installed.
        install_dir = os.path.dirname(__file__)
        default_file = os.path.join(install_dir, master_list)
        if os.path.exists(default_file):
            # Copy the file
            shutil.copy2(default_file, os.getcwd())

    # Read the master_list to get the poller fully qualified filenames
    table = Table.read(master_list, format="ascii.fast_no_header")
    poller_fqfile_list = table["col1"].tolist()

    print("Input file: {}".format(master_list))
    print("List of poller files: {}".format(poller_fqfile_list))
    metafunc.parametrize('dataset', poller_fqfile_list)


@pytest.mark.bigdata
@pytest.mark.slow
@pytest.mark.unit
@pytest.mark.skip
def test_run_svmpoller(tmpdir, dataset):
    """ Tests to read a series of poller files and process the contents of each as Single Visit Mosaic

        Characteristics of these tests:

        Success Criteria:
            The SVM processing returns a value of 0: Success or 1: Failure

        The input master_list file is a list of poller filenames, one filename per line.
        Each poller file must be obtained from a specified directory and read to obtain the
        names of the data files which need to be processed.

        This test file can be executed in the following manner:
            $ pytest -n # -s --basetemp=/internal/hladata/yourUniqueDirectoryHere --bigdata --slow
              --master_list /internal/hladata/input/master_poller_list.txt test_run_svmpoller.py >&
              test_svmpoller_output.txt &
            $ tail -f test_svmpoller_output.txt
          * The `-n #` option can be used to run tests in parallel if `pytest-xdist` has
            been installed where `#` is the number of cpus to use. THIS IS NOT ADVISED FOR USE.
          * Note: When running this test, the `--basetemp` directory should be set to a unique
            existing directory to avoid deleting previous test output.
          * A default master list exists in the tests/hla directory and contains 121 datasets.  This
            is probably NOT the list you want to use, but it allows you to see what this file should
            contain.

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

        # The dataset variable is a fully qualified poller file (CSV) with no header
        # Copy the poller file to the current working directory
        shutil.copy2(dataset, ".")

        # Get the poller file path, as well as its simple name
        poller_path = os.path.dirname(dataset)
        poller_file = os.path.basename(dataset)
        table = Table.read(poller_file, format="ascii")

        # Column "col8" contains fully qualified constituent filenames for the single visit, and
        # the fully qualified path designated in the poller file is the shared cache.
        file_list = table["col8"].tolist()

        # Check if the files to be processed are in the same directory as poller
        # directory, otherwise they need to be copied from the on-line cache
        for full_filename in file_list:
            filename = os.path.basename(full_filename)
            log.info("Looking for file {}".format(filename))
            local_path = os.path.join(poller_path, filename)
            if os.path.exists(local_path):
                shutil.copy2(local_path, ".")
            else:
                shutil.copy2(full_filename, ".")

        log.info("Obtained all input files for dataset {}.".format(dataset))

        # Run SVM pipeline processing
        return_value = runsinglehap.perform(poller_file)

    # Catch anything that happens as this dataset will be considered a failure, but
    # the processing of datasets should continue.  This is meant to catch
    # unexpected errors and generate sufficient output exception
    # information so algorithmic problems can be addressed.
    except Exception as except_details:
        traceback.print_exc()
        pytest.fail("TEST_RUN_SVMPOLLER. Exception Dataset: {}\n", dataset)
        return_value = 1

    assert return_value == 0

    # Return to original directory
    os.chdir(prevdir)
