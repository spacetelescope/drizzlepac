import sys
import traceback
import os
import datetime
import pytest
import numpy as np
from astropy.table import Table, vstack
from astropy.io import ascii

from .base_test import BaseHLATest
from drizzlepac import alignimages
import drizzlepac.hlautils.catalog_utils as catutils

class TestAlignMosaic(BaseHLATest):
    """ Process a large sample of ACS and WFC3 datasets to determine if they can be
        aligned to an astrometric standard.
    """

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
        randomCandidateTable = catutils.randomSelectFromCSV(input_file_path[0],
            inputNumEntries, inputSeedValue)

        # Invoke the methods which will handle acquiring/downloading the data from
        # MAST and perform the alignment
        percentSuccess = 0.0
        try:
            percentSuccess = self.align_randomFields (randomCandidateTable)
        except Exception:
            pass

        return(percentSuccess)

    def align_randomFields(self, randomTable):
        """ Process randomly selected fields (aka datasets) stored in an Astropy table.

            Each field is used as input to determine if it can be aligned to an
            astrometric standard.  The success or fail status for each test is retained
            as the overall success or fail statistic is the necessary output from
            this test.
        """

        numSuccess = 0
        numQualSuccess = 0
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
                    fitQual = datasetTable['fit_qual'][index[0]]

                    # Update the table with the datasetKey which is really just a counter
                    datasetTable['datasetKey'][:] = datasetKey
                    datasetTable['completed'][:] = True
                    datasetTable.write(outputName, format='ascii.ecsv')
                    #datasetTable.pprint(max_width=-1)
   
                    # Successful datasets
                    if (fitQual <= 2):
                        print("TEST_RANDOM. Successful Dataset (fit_qual <= 2): ", dataset, "\n")
                        numSuccess += 1
                    elif 2 < fitQual <= 4:
                        print("TEST_RANDOM. Qualified Successful Dataset (2 < fit_qual <= 4): ", dataset, "\n")
                        numQualSuccess += 1
                    # Unsuccessful datasets
                    else:
                        print("TEST_RANDOM. Unsuccessful Dataset (fit_qual = 5): ", dataset, "\n")
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
        print('TEST_RANDOM. Number of qualified successful tests: ', numQualSuccess)
        print('TEST_RANDOM. Number of unsuccessful tests: ', numUnsuccessful)
        print('TEST_RANDOM. Number of exception tests: ', numException)
        print('TEST_RANDOM. Percentage success/numberOfTests: ', numSuccess/numProcessedDatasets*100.0)
        print('TEST_RANDOM. Percentage success+qualsuccess/numberOfTests: ', (numSuccess+numQualSuccess)/numProcessedDatasets*100.0)
 
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
