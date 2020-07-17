""" Definition of PandasDFReader class for reading the SVM harvester file.

    This class is responsible for reading the SVM "harvester file" which
    contains a Pandas dataframe stored in an HDF5 file.  The harvester file
    is the collection of JSON output generated by the SVM quality analysis
    tests.

"""
import array
from bokeh.layouts import row
from bokeh.plotting import figure, output_file, show
from bokeh.models import ColumnDataSource, Label
import csv
import glob
import json
import logging
import os
import sys
import pandas as pd
import pathlib
import shutil
import traceback

from stsci.tools import logutil
from astropy.io import fits
import numpy as np

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)
__version__ = 0.1
__version_date__ = '08-Jun-2020'

DETECTOR_LEGEND = {'UVIS': 'magenta', 'IR': 'red', 'WFC': 'blue',
                   'SBC': 'yellow', 'HRC': 'black'}

def get_pandas_data(pandas_filename, data_columns, log_level=logutil.logging.NOTSET):
    """Load the harvested data, stored in a CSV file, into local arrays.

    Parameters
    ==========
    pandas_filename : str
        Name of the CSV file created by the harvester.
        
    data_columns : list
        List of column names to be extracted from the input dataframe.

    Returns
    =======
    data_colsDF : Pandas dataframe
        Dataframe which is a subset of the input Pandas dataframe written out as
        a CSV file.  The subset dataframe consists of only the requested columns
        and rows where all of the requested columns did not contain NaNs.

    """
    
    # Instantiate a Pandas Dataframe Reader (lazy instantiation)
    # df_handle = PandasDFReader_CSV("svm_qa_dataframe.csv")
    df_handle = PandasDFReader(pandas_filename, log_level=log_level)

    # In this particular case, the names of the desired columns do not
    # have to be further manipulated, for example, to add dataset specific
    # names.
    # 
    # Get the relevant column data, eliminating all rows which have NaNs
    # in any of the relevant columns.
    if pandas_filename.endswith('.h5'):
        data_colsDF = df_handle.get_columns_HDF5(data_columns)
    else:
        data_colsDF = df_handle.get_columns_CSV(data_columns)

    return data_colsDF




class PandasDFReader:
    """ A base class to read a Pandas dataframe generated by the HAP harvester.
    """
    def __init__(self, harvester_filename, log_level):
        # set logging level to user-specified level
        log.setLevel(log_level)
        self.log_level = log_level

        # Check for the existence of the input file.  There is no race condition
        # here, so this style of checking is fine.
        path = pathlib.Path(harvester_filename)
        if not path.is_file():

            raise ValueError('Pandas dataframe file {} does not exist.')

        self.harvester_filename = harvester_filename

        # Lazy attribute creation
        self.dataframe = pd.DataFrame()

        # Lists to contain all of the 'header' and 'general information' data, respectively
        self.header_cols = []
        self.gen_info_cols = []

    def get_columns_CSV(self, column_names):
        """ Method to do the actual reading of dataframe and get the data in the
            specified columns.

            Parameters
            ----------
            column_names : list of str
            A list of the column names which specify the desired data

            Returns
            -------
            column_data : Pandas dataframe
            A Pandas dataframe containing only the specified named columns

            *** OUT OF DATE
        """
        # Only read the input dataframe once
        if self.dataframe.empty:
            self.dataframe = pd.read_csv(self.harvester_filename)

        # Get the requested columns and eliminate all rows which have
        # Generate a new column in the HAP dataframe, 'inst_det'.  Also append
        # this column to the user requested columns.
        self.dataframe['inst_det'] = self.dataframe['gen_info.instrument'] + '/' + self.dataframe['gen_info.detector']

        column_data = self.extract_columns(column_names)

        return column_data

    def get_columns_HDF5(self, column_names, do_drop=False):
        """ Method to do the actual reading of dataframe and get and return the
            data in the specified columns.

            Parameters
            ----------
            column_names : list
            A list of the column names which specify the desired data

            do_drop : bool, optional
            Indicates whether or not rows with NaNs should be dropped.

            Returns
            -------
            column_dataDF : Pandas dataframe
            A Pandas dataframe containing only the specified named columns where
            any rows containing NaNs have been eliminated.
        """
        # Only read the input dataframe once
        if self.dataframe.empty:
            hdf5 = pd.HDFStore(self.harvester_filename, mode="r")

            # Get the zeroth key from the file as there is really only one dataframe
            # stored in the file - just do not assume its key.
            key0 = hdf5.keys()[0]
            self.dataframe = hdf5.get(key0)

            hdf5.close()

            # Generate a new column in the HAP dataframe so the instrument and detector
            # can be reported in the same entry as a HoverTool tooltip (e.g., ACS/SBC)
            self.dataframe['inst_det'] = self.dataframe['gen_info.instrument'] + '/' + self.dataframe['gen_info.detector']

            # Generate a new column in the HAP dataframe which contains a color associated
            # with each instrument/detector combination.  These colors can then used for the
            # graphics so the data is consistently represented by the same set of colors.
            self.dataframe['colormap'] = self.dataframe['gen_info.detector']  # Just create a new column
            for key, value in DETECTOR_LEGEND.items():
                self.dataframe.loc[self.dataframe['gen_info.detector'] == key, 'colormap'] = value

            # TO DO: Is getting all the header info wise as the values can have NaNs
            # Always get all the "header" and "general information" columns
            self.header_cols = [hd_cols for hd_cols in self.dataframe if 'header' in hd_cols]
            self.gen_info_cols = [hd_cols for hd_cols in self.dataframe if 'gen_info' in hd_cols]

        # Set up for the final columns to extract from the dataframe and remove any duplicates
        tmp_columns = list(column_names) + self.header_cols + self.gen_info_cols + ['inst_det', 'colormap']
        final_columns_to_extract = list(set(tmp_columns))

        column_dataDF = self.extract_columns(final_columns_to_extract)

        return column_dataDF

    def extract_columns(self, column_names):
        """ Helper method to get the requested columns

            Parameters
            ----------
            column_names : list of str
            List of column names to extract from the Pandas dataframe

            Returns
            -------
            column_dataDF : Pandas dataframe
            A Pandas dataframe containing the specified named columns
        """
        column_dataDF = pd.DataFrame()
        try:
            column_dataDF = self.dataframe.loc[:, column_names]
        # Columns may not be present in the dataframe for legitimate reasons
        except KeyError as err:
            log.warning("Column data missing from the Pandas dataframe. {}".format(err))
        # Some other possibly more serious problem
        except Exception:
            log.critical("Unable to extract column data from the Pandas dataframe, {}.\n".format(self.harvester_filename))
            sys.exit(1)

        return column_dataDF
