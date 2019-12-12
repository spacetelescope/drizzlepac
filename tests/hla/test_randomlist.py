""" This module is the high-level wrapper for running a test on a list of input
    datasets in order to collect statistics on alignment of each dataset to an
    astrometric catalog."""
import datetime
import os
import shutil
import numpy as np
from astropy.table import Table
import pytest

from drizzlepac import align as alignimages

def pytest_generate_tests(metafunc):
    """Get the command line option."""
    start_row = metafunc.config.option.start_row
    num_rows = metafunc.config.option.num_rows
    master_list = metafunc.config.option.master_list

    # Check to see if the specified file exists in the current working directory
    if not os.path.exists(master_list):
        # If not, copy the default file from the installation directory
        # Find where this module has been installed.
        install_dir = os.path.dirname(__file__)
        default_file = os.path.join(install_dir, master_list)
        if os.path.exists(default_file):
            # Copy the file
            shutil.copy2(default_file, os.getcwd())

    # Read a randomized table
    data_table = Table.read(master_list, format='ascii.csv')
    data_list = get_dataset_list(data_table)

    # Extract the subset rows
    start_row = int(start_row)
    end_row = start_row + int(num_rows)
    print("\nTEST_RANDOM. Start row: {}   Number of rows to process: {}.".format(start_row, num_rows))
    print("MASTER_TABLE: {}".format(master_list))
    random_candidate_table = data_list[start_row:end_row]
    print(random_candidate_table)
    metafunc.parametrize('dataset', random_candidate_table)


@pytest.mark.bigdata
@pytest.mark.slow
@pytest.mark.unit
def test_randomlist(tmpdir, dataset):
    """ Tests which validate whether mosaics can be aligned to an astrometric standard.

        Characteristics of these tests:
          * A reference WCS is generated based upon all the input images for
            the field.
          * A source astrometric catalog is created using the Photutils
            package to detect explicitly sources in the images.
          * An astrometric catalog is created to extract astrometric positions
            for the found sources in the input images' field-of-view using
            GAIADR2 (preferred) or GAIADR1.
          * Cross matching/fitting is done between found sources and catalog
            coordinates with the Tweakwcs package.
          * The quality of the fit is evaluated against a minimum threshold and
            potentially another fit algorithm is invoked or an alternative
            catalog is used in an effort to obtain a better quality fit.
          * If the option is set, the WCS information is updated for the
            input exposures. The default is False.
          * No mosaic is generated.
          * An output table containing characterizations of the process and
            associated fit is generated.

        Success Criteria:
          * Success criterion hard-coded for this test represents whether a
            statistical sample (70%) of ACS and WFC3 datasets were able to be
            aligned to within 10mas RMS.
              * RMS values are extracted from the table output from `perform_align`
              * This criterion will need to be determined by the user after the test has been
                run based on how many datasets were run and skipped.

        The input master_list CSV file is
        output from a database and lists associations and singletons for ACS
        and WFC3 instruments randomly sorted.  The actual data files are
        downloaded from MAST via astroquery.

        This test file can be executed in the following manner:
            $ pytest -n # -s --basetemp=/internal/hladata/yourUniqueDirectoryHere --bigdata --slow
              --master_list ACSWFC3ListDefault50.csv --start_row 0 --num_rows 50 test_randomlist.py >&
              test_random_output.txt &
            $ tail -f test_random_output.txt
          * The `-n #` option can be used to run tests in parallel if `pytest-xdist` has
            been installed where `#` is the number of cpus to use.
          * Note: When running this test, the `--basetemp` directory should be set to a unique
            existing directory to avoid deleting previous test output.
          * The default master list exists in the tests/hla directory and contains 50 datasets.  The
            full master list of thousands of datasets resides in Artifactory as ACSWFC3List.csv
            (https://bytesalad.stsci.edu/artifactory/hst-hla-pipeline/dev/master_lists).

    """
    print("TEST_RANDOM. Dataset: ", dataset)
    output_name = dataset + '.ecsv'

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

    try:

        dataset_table = alignimages.perform_align([dataset], archive=False,
                                                  clobber=True, debug=False,
                                                  update_hdr_wcs=False,
                                                  print_fit_parameters=True,
                                                  print_git_info=False,
                                                  output=False)

        # Filtered datasets
        if dataset_table['doProcess'].sum() == 0:
            pytest.skip("TEST_RANDOM. Filtered Dataset: {}.".format(dataset))
        # Datasets to process
        elif dataset_table['doProcess'].sum() > 0:
            # Determine images in dataset to be processed and the number of images
            # This is in case an image was filtered out (e.g., expotime = 0)
            index = np.where(dataset_table['doProcess'] == 1)[0]
            fit_qual = dataset_table['fit_qual'][index[0]]

            # Update the table with the dataset_key which is really just a counter
            dataset_table['completed'][:] = True
            dataset_table.write(output_name, format='ascii.ecsv')

            if fit_qual > 4:
                pytest.fail("TEST_RANDOM. Unsuccessful Dataset (fit_qual = 5): {}.".format(dataset))
            else:
                assert 0 < fit_qual <= 4

    # Catch anything that happens as this dataset will be considered a failure, but
    # the processing of datasets should continue.  This is meant to catch
    # unexpected errors and generate sufficient output exception
    # information so algorithmic problems can be addressed.
    except Exception as except_details:
        print(except_details)
        pytest.fail("TEST_RANDOM. Exception Dataset: {}\n", dataset)

    finally:
        # Perform some clean up
        if os.path.isfile('ref_cat.ecsv'):
            os.remove('ref_cat.ecsv')
        if os.path.isfile('refcatalog.cat'):
            os.remove('refcatalog.cat')
        for filename in os.listdir():
            if filename.endswith('flt.fits') or filename.endswith('flc.fits'):
                os.remove(filename)

    # Return to original directory
    os.chdir(prevdir)


def get_dataset_list(table_name):
    """ Standalone function to read the Astropy table and get the dataset names

    Parameters
    ==========
    table_name : str
        Filename of the input master CSV file containing individual
        images or association names, as well as observational
        information regarding the images

    Returns
    =======
    dataset_names: list
        List of individual image or association base (IPPSSOOT) names
    """

    dataset_names = []

    # Determine if the data is part of an association or is an individual image
    for imgid, asnid in zip(table_name['observationID'], table_name['asnID']):

        # Protect against incomplete lines, missing asnID values, ...
        # Happens when lines like "(236319 rows affected)" are included
        if isinstance(asnid, np.ma.core.MaskedConstant):
            continue

        # If the asnID is the string NONE, this is an individual image,
        # and it is necessary to get the individual image dataset name.
        # Otherwise, this is an association dataset, so just add the asnID.
        if asnid.upper() == "NONE":
            dataset_names.append(imgid)
        else:
            dataset_names.append(asnid)

    # Turn into a set to remove duplicate ASNID entries, only want 1 per ASN
    dataset = set()
    # Return results as a list of unique dataset names while retaining the original order
    return [x for x in dataset_names if x not in dataset and not dataset.add(x)]
