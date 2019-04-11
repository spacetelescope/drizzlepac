""" This module processes a test on a list of input datasets to collect
    statistics on the quality of the alignment of each dataset to an
    astrometric catalog."""
import datetime
import time
import os
import random
import numpy as np
from astropy.table import Table, vstack

from drizzlepac import alignimages
from .base_test import BaseHLATest

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

    Notes
    =====
    This routine is for future use.

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

class TestAlignMosaic(BaseHLATest):
    """ Process a large sample of ACS and WFC3 datasets to determine if they can be
        aligned to an astrometric standard.
    """

    def test_align_randomfields(self, start_row, num_rows):
        """ Wrapper to set up the test for aligning a large number of randomly
            selected fields (aka datasets) from a input ascii file (CSV).

            The wrapper provides the parameter settings for the underlying test,
            as well as implements the criterion for the overall success or failure
            of the test.
        """

        input_list_file = 'ACSWFC3List.csv'

        # Obtain the full path to the file containing the dataset field names
        self.input_repo = 'hst-hla-pipeline'
        self.tree = 'dev'
        self.input_loc = 'master_lists'
        input_file_path = self.get_data(input_list_file)

        # FOR FUTURE
        # Desired number of random entries for testing
        # num_rows = 50

        # Seed for random number generator
        # input_seed_value = 1

        # Randomly select a subset of field names (each field represented by a row) from
        # the master CSV file and return as an Astropy table
        # NOTE: Preserved for use in the future.
        # random_candidate_table = _random_select_from_csv(input_file_path[0],
        #                                                 num_rows,
        #                                                 input_seed_value)

        # Read a randomized table
        data_table = Table.read(input_file_path[0], format='ascii.csv')

        # Extract the subset rows
        start_row = int(start_row)
        end_row = start_row + int(num_rows)
        print("\nTEST_RANDOM. Start row: {}   Number of rows to process: {}.".format(start_row, num_rows))
        random_candidate_table = data_table[start_row:end_row]

        # Invoke the methods which will handle acquiring/downloading the data from
        # MAST and perform the alignment.  If an exception happens, just abort
        # out - no further analysis needed.
        percent_success = 0.0
        try:
            percent_success = self.align_randomfields(random_candidate_table)
        except Exception:
            pass

        return percent_success

    def align_randomfields(self, random_table):
        """ Process randomly selected fields (aka datasets) stored in an Astropy table.

            Each field is used as input to determine if it can be aligned to an
            astrometric standard.  The success or fail status for each test is retained
            as the overall success or fail statistic is the necessary output from
            this test.
        """

        num_success = 0
        num_qual_success = 0
        num_unsuccessful = 0
        num_exception = 0

        # Read the table and extract a list of each dataset name in IPPSSOOT format
        # which is either an association ID or an individual filename
        dataset_list = get_dataset_list(random_table)

        num_processed_datasets = len(dataset_list)
        print('TEST_RANDOM. Number of tests started: ', num_processed_datasets)

        # Process the dataset names in the list
        #
        # If the dataset name represents an association ID, the multiplicity
        # of images within the association need to be processed.  Otherwise,
        # the dataset is a single image.
        #
        # If the "alignment" of a field/dataset fails for any reason, trap
        # the exception and keep going.
        all_dataset_table = Table()
        dataset_key = -1
        print("TEST_RANDOM. Dataset List: ", dataset_list)
        for dataset in dataset_list:
            dataset_key += 1
            output_name = dataset + '.ecsv'

            print("TEST_RANDOM. Dataset: ", dataset)
            current_dt = datetime.datetime.now()
            print(str(current_dt))

            try:

                dataset_table = alignimages.perform_align([dataset], archive=False,
                                                          clobber=True, debug=False,
                                                          update_hdr_wcs=False,
                                                          print_fit_parameters=True,
                                                          print_git_info=False,
                                                          output=False)

                # Filtered datasets
                if dataset_table['doProcess'].sum() == 0:
                    print("TEST_RANDOM. Filtered Dataset: ", dataset, "\n")
                    num_processed_datasets -= 1
                # Datasets to process
                elif dataset_table['doProcess'].sum() > 0:
                    # Determine images in dataset to be processed and the number of images
                    # This is in case an image was filtered out (e.g., expotime = 0)
                    index = np.where(dataset_table['doProcess'] == 1)[0]
                    fit_qual = dataset_table['fit_qual'][index[0]]

                    # Update the table with the dataset_key which is really just a counter
                    dataset_table['datasetKey'][:] = dataset_key
                    dataset_table['completed'][:] = True
                    dataset_table.write(output_name, format='ascii.ecsv')

                    # Successful datasets
                    if 0 < fit_qual <= 2:
                        print("TEST_RANDOM. Successful Dataset (fit_qual <= 2): ", dataset, "\n")
                        num_success += 1
                    elif 2 < fit_qual <= 4:
                        print("TEST_RANDOM. Qualified Successful Dataset (2 < fit_qual <= 4): ",
                              dataset, "\n")
                        num_qual_success += 1
                    # Unsuccessful datasets
                    else:
                        print("TEST_RANDOM. Unsuccessful Dataset (fit_qual = 5): ", dataset, "\n")
                        num_unsuccessful += 1

                # Append the latest dataset table to the summary table
                all_dataset_table = vstack([all_dataset_table, dataset_table])

            # Catch anything that happens as this dataset will be considered a failure, but
            # the processing of datasets should continue.  This is meant to catch
            # unexpected errors and generate sufficient output exception
            # information so algorithmic problems can be addressed.
            except Exception as except_details:
                print(except_details)
                print("TEST_RANDOM. Exception Dataset: ", dataset, "\n")
                num_exception += 1
                continue

        # Perform some clean up
        if os.path.isfile('ref_cat.ecsv'):
            os.remove('ref_cat.ecsv')
        if os.path.isfile('refcatalog.cat'):
            os.remove('refcatalog.cat')
        for filename in os.listdir():
            if filename.endswith('flt.fits') or filename.endswith('flc.fits'):
                os.remove(filename)

        # Write out the summary table for all processed datasets - generate a unique output
        # name based on seconds since the epoch in units of seconds
        all_dataset_table.write('randomResults{}.ecsv'.format(int(time.time())),
                                format='ascii.ecsv')

        # Determine the percent success over all datasets processed
        percent_success = num_success / num_processed_datasets
        print('TEST_RANDOM. Number of tests (excluding filtered): ', num_processed_datasets)
        print('TEST_RANDOM. Number of successful tests: ', num_success)
        print('TEST_RANDOM. Number of qualified successful tests: ', num_qual_success)
        print('TEST_RANDOM. Number of unsuccessful tests: ', num_unsuccessful)
        print('TEST_RANDOM. Number of exception tests: ', num_exception)
        print('TEST_RANDOM. Percentage success/numberOfTests: ',
              num_success / num_processed_datasets * 100.0)
        print('TEST_RANDOM. Percentage success+qualsuccess/numberOfTests: ',
              (num_success + num_qual_success) / num_processed_datasets * 100.0)

        return percent_success

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

        # If the asnID is the string NONE, this is an individual image,
        # and it is necessary to get the individual image dataset name.
        # Otherwise, this is an association dataset, so just add the asnID.
        if asnid.upper() == "NONE":
            dataset_names.append(imgid)
        else:
            dataset_names.append(asnid)

    return dataset_names
