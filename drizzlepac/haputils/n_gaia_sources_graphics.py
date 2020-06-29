#!/usr/bin/env python

"""This script harvests """


# Standard library imports
import argparse
import logging
import os
import pdb
import sys

from bokeh.layouts import gridplot, row
from bokeh.plotting import figure, output_file, show, save
from bokeh.models import ColumnDataSource, Label
from bokeh.models.tools import HoverTool

# Local application imports
from drizzlepac.haputils.pandas_utils import PandasDFReader
from drizzlepac.haputils.diagnostic_json_harvester import h5load
from stsci.tools import logutil

__taskname__ = 'n_gaia_sources_graphics'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)

# ------------------------------------------------------------------------------------------------------------
# Global variable declaration
new_df_colnames = {'gen_info.instrument': 'instrument',
                   'gen_info.detector': 'detector',
                   'gen_info.filter': 'filter',
                   'gen_info.dataset': 'dataset',
                   'gen_info.proposal_id': 'proposal_id',
                   'gen_info.visit': 'visit',
                   'Number_of_GAIA_sources.Number_of_GAIA_sources': 'Number_of_GAIA_sources'}

# ------------------------------------------------------------------------------------------------------------


def does_file_exist(dataframe_filename, log_level=logutil.logging.INFO):
    """Checks that the specified file exists. If it does not, an exception is raised.
    
    Parameters
    ----------
    dataframe_filename : str
        Name of the storage file for the Pandas dataframe created by the harvester.
    
    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the 
        .log file. Default value is 20, or 'info'.

    Returns
    -------
    Nothing.
    """
    log.setLevel(log_level)

    # Verify the input file exists
    if not os.path.exists(dataframe_filename):
        err_msg = "Harvester HDF5 File {} does not exist.".format(dataframe_filename)
        log.critical(err_msg)
        raise Exception(err_msg)

# ------------------------------------------------------------------------------------------------------------

def get_data(dataframe_filename, log_level=logutil.logging.INFO):
    """Extract releavant information from dataframe stored in user-specified .h5 file and return it as a
    Pandas DataFrame
    
    Parameters
    ----------
    dataframe_filename : str
        Name of the storage file for the Pandas dataframe created by the harvester.
    
    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the 
        .log file. Default value is 20, or 'info'.

    Returns
    -------
    return_dataframe : Pandas DataFrame
        Pandas DataFrame containing just the data relevant for plotting.
    """
    log.setLevel(log_level)
    
    # check to make sure that the input file does indeed exist
    does_file_exist(dataframe_filename, log_level=log_level)
    
    # Instantiate a Pandas Dataframe Reader (lazy instantiation)
    df_handle = PandasDFReader(dataframe_filename, log_level=logutil.logging.NOTSET)

    # Get the relevant column data, eliminating all rows which have NaNs
    # in any of the relevant columns.
    return_dataframe = df_handle.get_columns_HDF5(new_df_colnames.keys())
    # log.info("Photometry_graphics. Photometric data has been retrieved from the storage Pandas dataframe: {}.\n".format(storage_filename))

    # Rename the columns to abbreviated text as the graph titles further
    # document the information.
    for old_col_name, new_col_name in new_df_colnames.items():
        return_dataframe.rename(columns={old_col_name: new_col_name}, inplace=True)

    # Keep only lines for total-level products. Discard filter-level and exposure-level lines from the
    # dataframe.
    initial_nrows = len(return_dataframe.index)
    if initial_nrows == 1:
        initial_row_descriptor = "row"
    else:
        initial_row_descriptor = "rows"
    nrows_removed = 0
    for item in return_dataframe.index:
        if "total" not in item.split("_"):
            log.debug("Removed dataframe line with index value {}".format(item))
            return_dataframe = return_dataframe.drop([item])
            nrows_removed += 1
    final_nrows = len(return_dataframe.index)
    if final_nrows == 1:
        final_row_descriptor = "row"
    else:
        final_row_descriptor = "rows"
    log.info("Initial dataframe length: {} {}".format(initial_nrows, initial_row_descriptor))
    log.info("Number of row(s) removed: {}".format(nrows_removed))
    log.info("Final dataframe length:   {} {}".format(final_nrows, final_row_descriptor))

    return return_dataframe

# ------------------------------------------------------------------------------------------------------------

def n_gaia_sources_graphics_driver(dataframe_filename, output_base_filename='n_gaia_sources_graphics', log_level=logutil.logging.INFO):
    """This is the primary driver subroutine for this script.

    Parameters
    ----------
    dataframe_filename : str
        Name of the storage file for the Pandas dataframe created by the harvester.

    output_base_filename : str
        Base name for the HMTL file generated by Bokeh.

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the .log 
        file. Default value is 20, or 'info'.

    Returns
    -------
    Nothing
    """
    log.setLevel(log_level)

    # 1: get relevant data values from user-specified input file
    dataframe = get_data(dataframe_filename, log_level=log_level)
    
    # 2: generate visualization with bokeh


# ======================================================================================================================


if __name__ == "__main__":
    # Process command-line inputs with argparse
    parser = argparse.ArgumentParser(description='Read the harvested Pandas dataframe stored as and HDF5 file.')

    parser.add_argument('dataframe_filename', help='File which holds the Pandas dataframe in an HDF5 file.')
    parser.add_argument('-o', '--output_base_filename', required=False, default="n_gaia_sources_graphics",
                        help='Name of the output base filename (filename without the extension) for the  '
                             'HTML file generated by Bokeh. Unless explicitly specified, the default value '
                             'is "n_gaia_sources_graphics".'
                             ' ')
    parser.add_argument('-l', '--log_level', required=False, default='info',
                        choices=['critical', 'error', 'warning', 'info', 'debug'],
                        help='The desired level of verboseness in the log statements displayed on the screen '
                             'and written to the .log file. The level of verboseness from left to right, and '
                             'includes all log statements with a log_level left of the specified level. '
                             'Specifying "critical" will only record/display "critical" log statements, and '
                             'specifying "error" will record/display both "error" and "critical" log '
                             'statements, and so on.')
    user_args = parser.parse_args()

    # Set up logging
    log_dict = {"critical": logutil.logging.CRITICAL,
                "error": logutil.logging.ERROR,
                "warning": logutil.logging.WARNING,
                "info": logutil.logging.INFO,
                "debug": logutil.logging.DEBUG}
    log_level = log_dict[user_args.log_level]
    log.setLevel(log_level)

    # Verify the input file exists
    does_file_exist(user_args.dataframe_filename, log_level=log_level)

    print(user_args.output_base_filename)

    n_gaia_sources_graphics_driver(user_args.dataframe_filename,
                                   output_base_filename=user_args.output_base_filename,
                                   log_level=log_level)