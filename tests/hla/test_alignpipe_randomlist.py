""" This module is the high-level wrapper for running a test on a list of input
    datasets in order to process the dataset in manner similar to the pipeline
    processing (updatewcs, alignment, astrodrizzle).
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
from astropy.io import fits

from drizzlepac.hlautils import astroquery_utils as aqutils
from drizzlepac import runastrodriz


log = logutil.create_logger('test_alignpipe_randomlist', level=logutil.logging.INFO, stream=sys.stdout)

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
def test_alignpipe_randomlist(tmpdir, dataset):
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
              --master_list ACSWFC3ListDefault50.csv --start_row 0 --num_rows 50 test_pipe_randomlist.py >&
              test_pipe_random_output.txt &
            $ tail -f test_pipe_random_output.txt
          * The `-n #` option can be used to run tests in parallel if `pytest-xdist` has
            been installed where `#` is the number of cpus to use.
          * Note: When running this test, the `--basetemp` directory should be set to a unique
            existing directory to avoid deleting previous test output.
          * The default master list exists in the tests/hla directory and contains 50 datasets.  The
            full master list of thousands of datasets resides in Artifactory as ACSWFC3List.csv
            (https://bytesalad.stsci.edu/artifactory/hst-hla-pipeline/dev/master_lists).

    """
    # Start by resetting the logging for all the modules used
    rl = logging.getLogger('stwcs.wcsutil.headerlet')
    if len(rl.handlers) > 1: del rl.handlers[-1]

    print("TEST_RANDOM_ALIGN. Dataset: ", dataset)
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

    return_value = 1

    try:

        # The dataset name is an ipppssoot only.
        dataset = dataset.lower()

        # Download the data necessary for processing - The dataset list contains only the ipppssoot of a
        # singleton or an ASN.  Need to get the FLT/FLC files as runastrodriz expects them to be present.
        input_parameter = "RAW"
        filename_for_runastrodriz = dataset.lower() + "_raw.fits"
        if dataset.endswith("0"):
            input_parameter = "ASN"
            filename_for_runastrodriz = dataset.lower() + "_asn.fits"
        log.info("Input parameter: {}".format(input_parameter))
        files_on_disk = check_disk_get_data([dataset], suffix=input_parameter, archive=False, clobber=False)

        log.info("Dataset: {}".format(dataset))
        if os.path.exists('mastDownload'):
            shutil.rmtree('mastDownload')
        log.info("\nFiles: {}".format(files_on_disk))


        # Insure environment variables are set for full processing
        os.environ['ASTROMETRY_STEP_CONTROL'] = 'on'
        os.environ['ASTROMETRY_COMPUTE_APOSTERIORI'] = 'on'
        os.environ['ASTROMETRY_APPLY_APRIORI'] = 'on'

        flts = sorted(glob.glob('*fl?.fits'))
        jref_dir = 'jref.old/' if '16r12191j_mdz' in fits.getval(flts[0], 'mdriztab') else 'jref/'
        os.environ['jref'] = os.path.join(os.environ['crrefer'], jref_dir)
        log.info("JREF: {}".format(os.environ['jref']))

        # Run pipeline processing using
        # runastrodriz accepts *_raw.fits or *_asn.fits, but it assumes the *_fl[t|c].fits files
        # are also present
        runastrodriz.process(filename_for_runastrodriz, force=True)

        return_value = 0

    # Catch anything that happens as this dataset will be considered a failure, but
    # the processing of datasets should continue.  This is meant to catch
    # unexpected errors and generate sufficient output exception
    # information so algorithmic problems can be addressed.
    except Exception as except_details:
        traceback.print_exc()
        pytest.fail("TEST_RANDOM. Exception Dataset: {}\n", dataset)
        return_value = 1

    """
    finally:
        # Perform some clean up
        if os.path.isfile('ref_cat.ecsv'):
            os.remove('ref_cat.ecsv')
        if os.path.isfile('refcatalog.cat'):
            os.remove('refcatalog.cat')
        for filename in os.listdir():
            if filename.endswith('flt.fits') or filename.endswith('flc.fits'):
                os.remove(filename)
        """

    assert return_value == 0
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


def check_disk_get_data(input_list, **pars):
    """Verify that all specified files are present. If not, retrieve them from MAST.

    Parameters
    ----------
    input_list : list
        List of one or more calibrated fits images that will be used for catalog generation.

    Returns
    =======
    total_input_list: list
        list of full filenames

    """
    reload(aqutils)

    empty_list = []
    retrieve_list = []    # Actual files retrieved via astroquery and resident on disk
    candidate_list = []   # File names gathered from *_asn.fits file
    ipppssoot_list = []   # ipppssoot names used to avoid duplicate downloads
    total_input_list = []  # Output full filename list of data on disk
    member_suffix = '_flc.fits'

    # Get the suffix values
    suffix_to_check = pars.get("suffix")
    # List set up with FLT before FLC to ensure both are retrieved if they both exist
    suffix_to_retrieve = ["ASN", "FLT", "FLC"]
    if suffix_to_check == "RAW":
        suffix_to_retrieve = ["RAW", "FLT", "FLC"]

    # Loop over the input_list to determine if the item in the input_list is a full association file
    # (*_asn.fits), a full individual image file (aka singleton, *_flt.fits), or a root name specification
    # (association or singleton, ipppssoot).
    for input_item in input_list:
        log.info('Input item: {}'.format(input_item))
        indx = input_item.find('_')

        # Input with a suffix (_xxx.fits)
        if indx != -1:
            lc_input_item = input_item.lower()
            suffix = lc_input_item[indx + 1:indx + 4]
            log.info('file: {}'.format(lc_input_item))
            # For an association, need to open the table and read the image names as this could
            # be a custom association.  The assumption is this file is on local disk when specified
            # in this manner (vs just the ipppssoot of the association).
            # This "if" block just collects the wanted full file names.
            if suffix == 'asn':
                try:
                    asntab = Table.read(input_item, format='fits')
                except FileNotFoundError:
                    log.error('File {} not found.'.format(input_item))
                    return(empty_list)
                for row in asntab:
                    if row['MEMTYPE'].startswith('PROD'):
                        continue
                    memname = row['MEMNAME'].lower().strip()
                    # Need to check if the MEMNAME is a full filename or an ipppssoot
                    if memname.find('_') != -1:
                        candidate_list.append(memname)
                    else:
                        # Define suffix for all members based on what files are present
                        if not os.path.exists(memname + member_suffix):
                            member_suffix = '_flt.fits'

                        candidate_list.append(memname + member_suffix)
            elif suffix in ['flc', 'flt']:
                if lc_input_item not in candidate_list:
                    candidate_list.append(lc_input_item)
            else:
                log.error(
                    'Inappropriate file suffix: {}.  Looking for "asn.fits", '
                    '"flc.fits", or "flt.fits".'.format(
                        suffix))
                return (empty_list)

        # Input is an ipppssoot (association or singleton), nine characters by definition.
        # This "else" block actually downloads the data specified as ipppssoot.
        elif len(input_item) == 9:
            try:
                if input_item not in ipppssoot_list:
                    # An ipppssoot of an individual file which is part of an association cannot be
                    # retrieved from MAST
                    log.info("Collect data: {} Suffix: {}".format(input_item, suffix_to_retrieve))
                    for filetype in suffix_to_retrieve:
                        retrieve_list += aqutils.retrieve_observation(input_item, suffix=filetype)
                    log.info("Collected data: {}".format(retrieve_list))

                    # If the retrieved list is not empty, add filename(s) to the total_input_list.
                    # Also, update the ipppssoot_list so we do not try to download the data again.  Need
                    # to do this since retrieve_list can be empty because (1) data cannot be acquired (error)
                    # or (2) data is already on disk (ok).
                    if retrieve_list:
                        total_input_list += retrieve_list
                        ipppssoot_list.append(input_item)
                    else:
                        log.error('File {} cannot be retrieved from MAST.'.format(input_item))
                        return(empty_list)
            except Exception:
                log.info("Exception in check_disk_get_data")
                exc_type, exc_value, exc_tb = sys.exc_info()
                traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)

    # Only the retrieve_list files via astroquery have been put into the total_input_list thus far.
    # Now check candidate_list to detect or acquire the requested files from MAST via astroquery.
    for file in candidate_list:
        # If the file is found on disk, add it to the total_input_list and continue
        if glob.glob(file):
            total_input_list.append(file)
            continue
        else:
            log.error('File {} cannot be found on the local disk.'.format(file))
            return(empty_list)

    log.info("TOTAL INPUT LIST: {}".format(total_input_list))
    return(total_input_list)
