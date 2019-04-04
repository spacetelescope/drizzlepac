import sys
import traceback
import os
import datetime
import random
import pytest

import numpy as np
from astropy.table import Table, vstack
from astropy.io import ascii

from .base_test import BaseHLATest
from drizzlepac import alignimages
import drizzlepac.hlautils.catalog_utils as catutils


def _random_select_from_csv(table_name, num_entries, seed_value):
    """
    Function to extract random entries (lines) from a CSV file

    The function, _random_select_from_csv, allows the user to specify the
    desired number of entries to obtain from a comma-separated values (CSV)
    file which contains data extracted from the STScI archive for ACS, WFC3,
    and WFPC2 instruments. The data are comprised of the following information
    as identified by the database name: observationID, trgName, trgposRA,
    trgPosDec, detector, aperture, filters, enrBandpassName, timMin, timMax,
    timExposure,and asnID.  The entries from the input table are chosen
    at random with no duplication, and returned as an Astropy table.

    Parameters
    ==========
    table_name: str
        Filename of the input master CSV file containing individual
        images or association names, as well as observational
        information regarding the images

    num_entries : int
        Number of entries/rows to extract from the master input CSV file

    seed_value : int
        Value used to initialize the random number generator for the
        selection of random entries

    Returns
    =======
    output_table : object
        Astropy Table object

    """
    # Initialize the random number generator
    random.seed(seed_value)

    # Get the contents of the table
    data_table = Table.read(table_name, format='ascii.csv')

    # Generate a sequence of integers the size of the table, and then
    # obtain a random subset of the sequence with no duplicate selections
    subset = random.sample(range(len(data_table)), num_entries)

    # Extract the subset rows...
    output_table = data_table[subset]

    # Returns the outputTable which is an Astropy Table object
    return output_table


@pytest.mark.bigdata
class TestAlignMosaic(BaseHLATest):
    """ Tests which validate whether mosaics can be aligned to an astrometric standard.

        Characeteristics of these tests:
          * A single astrometric catalog was obtained with both GAIA and non-GAIA
            (PanSTARRS?) sources for the entire combined field-of-view using the GSSS
            server.
              * These tests assume combined input exposures for each test have enough
                astrometric sources from the external catalog to perform a high quality
                fit.
          * This test only determines the fit between sources
            extracted from the images by Tweakreg and the source positions included in
            the astrometric catalog.
          * The WCS information for the input exposures do not get updated in this test.
          * No mosaic gets generated.

        Success Criteria:
          * Success criteria hard-coded for this test represents 10mas RMS for the
            WFC3 images based on the fit of source positions to the astrometric catalog
            source positions.
              * RMS values are extracted from optional shiftfile output from `tweakreg`
              * Number of stars used for the fit and other information is not available
                with the current version of `tweakreg`.

        The environment variable needs to be set in the following manner:
            export TEST_BIGDATA=https://bytesalad.stsci.edu/artifactory/
            OR
            export TEST_BIGDATA=/Users/YourNameHere/TestDataDirectory/

        For this test, the TEST_BIGDATA defines the root of the location  where the CSV
        file is stored.  The full path is TEST_BIGDATA/self.input_repo/self.tree/self.input_loc.
        The CSV is output from a database and lists associations and singletons for ACS
        and WFC3 instruments randomly sorted.  The actual data files are downloaded from MAST
        via astroquery.

        This test file can be executed in the following manner:
            $ pytest -s --bigdata test_randomlist.py >& test_random_output.txt &
            $ tail -f test_random_output.txt
    """

    @pytest.mark.xfail
    @pytest.mark.slow
    def test_align_randomFields(self):
        """ Wrapper to set up the test for aligning a large number of randomly
            selected fields (aka datasets) from a input ascii file (CSV).

            The wrapper provides the parameter settings for the underlying test,
            as well as implements the criterion for the overall success or failure
            of the test.
        """

        inputListFile = 'ACSWFC3List.csv'

        # Desired number of random entries for testing
        inputNumEntries = 50

        # Seed for random number generator
        inputSeedValue = 1

        # Obtain the full path to the file containing the dataset field names
        self.input_repo = 'hst-hla-pipeline'
        self.tree = 'dev'
        self.input_loc = 'master_lists'
        input_file_path = self.get_data(inputListFile)

        # Randomly select a subset of field names (each field represented by a row) from
        # the master CSV file and return as an Astropy table
        randomCandidateTable = catutils._random_select_from_csv(
            input_file_path[0], inputNumEntries, inputSeedValue
        )

        # Invoke the methods which will handle acquiring/downloading the data from
        # MAST and perform the alignment
        percentSuccess = 0.0
        try:
            percentSuccess = self.align_randomFields (randomCandidateTable)
        except Exception:
            pass

        assert(percentSuccess >= 0.70)

    def align_randomFields(self, randomTable):
        """ Process randomly selected fields (aka datasets) stored in an Astropy table.

            Each field is used as input to determine if it can be aligned to an
            astrometric standard.  The success or fail status for each test is retained
            as the overall success or fail statistic is the necessary output from
            this test.
        """

        numSuccess = 0
        numUnsuccessful = 0
        numException = 0
        numProcessedDatasets = 0

        # Read the table and extract a list of each dataset name in IPPSSOOT format
        # which is either an association ID or an individual filename
        dataset_list = get_dataset_list(randomTable)

        numProcessedDatasets = len(dataset_list)
        numStartTests  = numProcessedDatasets

        # Process the dataset names in the list
        #
        # If the dataset name represents an association ID, the multiplicity
        # of images within the association need to be processed.  Otherwise,
        # the dataset is a single image.
        #
        # If the "alignment" of a field/dataset fails for any reason, trap
        # the exception and keep going.
        allDatasetTable = Table()
        datasetKey = -1
        print("TEST_RANDOM. Dataset List: ", dataset_list)
        for dataset in dataset_list:
            datasetKey += 1
            outputName = dataset + '.ecsv'

            print("TEST_RANDOM. Dataset: ", dataset)
            currentDT = datetime.datetime.now()
            print(str(currentDT))

            try:

                datasetTable = alignimages.perform_align([dataset],archive=False,clobber=True,debug=False,
                    update_hdr_wcs=False,print_fit_parameters=True,print_git_info=False,output=False)

                # Filtered datasets
                if datasetTable['doProcess'].sum() == 0:
                    print("TEST_RANDOM. Filtered Dataset: ", dataset, "\n")
                    numProcessedDatasets -= 1;
                # Datasets to process
                elif datasetTable['doProcess'].sum() > 0:
                    # Determine images in dataset to be processed and the number of images
                    # This is in case an image was filtered out (e.g., expotime = 0)
                    index = np.where(datasetTable['doProcess']==1)[0]
                    sumOfStatus = datasetTable['status'][index].sum()

                    # Update the table with the datasetKey which is really just a counter
                    datasetTable['datasetKey'][:] = datasetKey
                    datasetTable['completed'][:] = True
                    datasetTable.write(outputName, format='ascii.ecsv')
                    #datasetTable.pprint(max_width=-1)

                    # Successful datasets
                    if (sumOfStatus == 0):
                        print("TEST_RANDOM. Successful Dataset: ", dataset, "\n")
                        numSuccess += 1
                    # Unsuccessful datasets
                    else:
                        print("TEST_RANDOM. Unsuccessful Dataset: ", dataset, "\n")
                        numUnsuccessful += 1

                # Append the latest dataset table to the summary table
                allDatasetTable = vstack([allDatasetTable, datasetTable])

            # Catch anything that happens as this dataset will be considered a failure, but
            # the processing of datasets should continue.  Generate sufficient output exception
            # information so problems can be addressed.
            except Exception:

                exc_type, exc_value, exc_tb = sys.exc_info()
                traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)
                print("TEST_RANDOM. Exception Dataset: ", dataset, "\n")
                numException += 1
                continue

        # Perform some clean up
        if os.path.isfile('ref_cat.ecsv'):
            os.remove('ref_cat.ecsv')
        if os.path.isfile('refcatalog.cat'):
            os.remove('refcatalog.cat')
        for filename in os.listdir():
            if filename.endswith('flt.fits') or filename.endswith('flc.fits'):
                os.unlink(filename)

        # Write out the table
        allDatasetTable.write('resultsBigTest.ecsv', format='ascii.ecsv')

        # Determine the percent success over all datasets processed
        percentSuccess = numSuccess/numProcessedDatasets
        print('TEST_RANDOM. Number of tests started: ', numStartTests)
        print('TEST_RANDOM. Number of tests (excluding filtered): ', numProcessedDatasets)
        print('TEST_RANDOM. Number of successful tests: ', numSuccess)
        print('TEST_RANDOM. Number of unsuccessful tests: ', numUnsuccessful)
        print('TEST_RANDOM. Number of exception tests: ', numException)
        print('TEST_RANDOM. Percentage success/numberOfTests: ', numSuccess/numProcessedDatasets*100.0)

        return percentSuccess

def get_dataset_list(tableName):
    """ Standalone function to read the Astropy table and get the dataset names

    Parameters
    ==========
    tableName : str
        Filename of the input master CSV file containing individual
        images or association names, as well as observational
        information regarding the images

    Returns
    =======
    datasetNames: list
        List of individual image or association base (IPPSSOOT) names
    """

    #dataFromTable = Table.read(filename, format='ascii')
    datasetIDs = tableName['observationID']
    asnIDs     = tableName['asnID']

    datasetNames = []

    # Determine if the data is part of an association or is an individual image
    for imgid,asnid in zip(datasetIDs,asnIDs):

        # If the asnID is the string NONE, this is an individual image,
        # and it is necessary to get the individual image dataset name.
        # Otherwise, this is an association dataset, so just add the asnID.
        if (asnid.upper() == "NONE"):
            datasetNames.append(imgid)
        else:
            datasetNames.append(asnid)

    return datasetNames
