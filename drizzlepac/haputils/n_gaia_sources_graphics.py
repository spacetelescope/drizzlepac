#!/usr/bin/env python

"""This script harvests the hdf5 .h5 file generated by diagnostic_json_harvester.py and generates a
histogram of the number of GAIA sources as a function of the number of total-level products present in the
dataframe. (x axis: Number of GAIA sources; y axis: total-level products)"""

# Standard library imports
import argparse
import os
import sys

# Related third party imports
from bokeh.plotting import figure, output_file, show
import numpy as np

# Local application imports
from drizzlepac.haputils.pandas_utils import PandasDFReader
from stsci.tools import logutil

__taskname__ = 'n_gaia_sources_graphics'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)

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


def generate_histogram(dataframe, output_base_filename, log_level=logutil.logging.INFO):
    """Generate a histogram of the 'number of GAIA sources' data and write it to an .html file.

    Parameters
    ----------
    dataframe : Pandas DataFrame
        Pandas DataFrame containing just the data relevant for plotting.

    output_base_filename : str
        Name of the output base filename (filename without the extension) for the HTML file generated by Bokeh

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the
        .log file. Default value is 20, or 'info'.

    Returns
    -------
    Nothing.
    """
    log.setLevel(log_level)
    out_filename = output_base_filename+".html"
    # clobber existing file if it exists
    if os.path.exists(out_filename):
        os.remove(out_filename)

    # Generate histogram data
    hist, edges = np.histogram(dataframe['Number_of_GAIA_sources'].values)

    # initialize plot
    p = figure(title="Histogram of the number of GAIA sources",
               x_axis_label="Number of GAIA sources",
               y_axis_label="Number of total-level exposures")

    # generate histogram plot
    p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:], line_color="white")

    # write to out_filename
    output_file(out_filename)
    log.info("Histogram written to {}".format(out_filename))
    if log_level == logutil.logging.DEBUG:
        show(p)


# ------------------------------------------------------------------------------------------------------------

def get_data(dataframe_filename, log_level=logutil.logging.INFO):
    """Extract relevant information from dataframe stored in user-specified .h5 file and return it as a
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
    new_df_colnames = {'gen_info.instrument': 'instrument',
                       'gen_info.detector': 'detector',
                       'gen_info.filter': 'filter',
                       'gen_info.dataset': 'dataset',
                       'gen_info.proposal_id': 'proposal_id',
                       'gen_info.visit': 'visit',
                       'Number_of_GAIA_sources.Number_of_GAIA_sources': 'Number_of_GAIA_sources'}

    return_dataframe = df_handle.get_columns_HDF5(new_df_colnames.keys())

    # Rename the columns to abbreviated text as the graph titles further
    # document the information.
    for old_col_name, new_col_name in new_df_colnames.items():
        return_dataframe.rename(columns={old_col_name: new_col_name}, inplace=True)

    # Discard rows that contain entries for filter-level and exposure-level products from dataframe. Only
    # rows containing entries for total-level products will remain.
    initial_nrows = len(return_dataframe.index)
    nrows_removed = 0
    for item in return_dataframe.index:
        if "total" not in item.split("_"):
            log.debug("Removed dataframe line with index value {}".format(item))
            return_dataframe = return_dataframe.drop([item])
            nrows_removed += 1
    final_nrows = len(return_dataframe.index)

    # report initial number of rows in dataframe, number of rows removed, and the number of rows remaining.
    if initial_nrows == 1:
        initial_row_descriptor = "row"
    else:
        initial_row_descriptor = "rows"
    if final_nrows == 1:
        final_row_descriptor = "row"
    else:
        final_row_descriptor = "rows"
    log.info("Initial dataframe length: {} {}".format(initial_nrows, initial_row_descriptor))
    log.info("Number of row(s) removed: {}".format(nrows_removed))
    log.info("Final dataframe length:   {} {}".format(final_nrows, final_row_descriptor))

    # bail if no rows remain that contain entries for total-level products
    if final_nrows == 0:
        err_msg = "There are no rows left remaining to plot! Exiting...".format(dataframe_filename)
        log.critical(err_msg)
        raise Exception(err_msg)
    return return_dataframe

# ------------------------------------------------------------------------------------------------------------


def n_gaia_sources_graphics_driver(dataframe_filename, output_base_filename='n_gaia_sources_graphics',
                                   log_level=logutil.logging.INFO):
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
    generate_histogram(dataframe, output_base_filename, log_level=log_level)

# ======================================================================================================================


if __name__ == "__main__":
    # Process command-line inputs with argparse
    parser = argparse.ArgumentParser(description='Read the harvested Pandas dataframe stored as and HDF5 '
                                                 'file.')

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
                             'statements, and so on. If log_level is set to "debug", the histogram plot .html'
                             'file will be displayed in the default web browser.')
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

    # execute main driver subroutine
    n_gaia_sources_graphics_driver(user_args.dataframe_filename,
                                   output_base_filename=user_args.output_base_filename,
                                   log_level=log_level)
