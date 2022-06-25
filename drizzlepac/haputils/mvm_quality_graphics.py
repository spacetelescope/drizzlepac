"""Code that evaluates the quality of the MVM products generated by the drizzlepac package.

The JSON files generated here can be converted directly into a Pandas DataFrame
using the syntax:

>>> import json
>>> import pandas as pd
>>> with open("<rootname>_astrometry_resids.json") as jfile:
>>>     resids = json.load(jfile)
>>> pdtab = pd.DataFrame(resids)

These DataFrames can then be concatenated using:

>>> allpd = pdtab.concat([pdtab2, pdtab3])

where 'pdtab2' and 'pdtab3' are DataFrames generated from other datasets.  For
more information on how to merge DataFrames, see

https://pandas.pydata.org/pandas-docs/stable/user_guide/merging.html

Visualization of these Pandas DataFrames with Bokeh can follow the example
from:

https://programminghistorian.org/en/lessons/visualizing-with-bokeh

PUBLIC INTERFACE FOR THIS MODULE
 - build_svm_plots(HDF5_FILE, output_basename='', display_plot=False):

"""

# Standard library imports
import argparse
import collections
from datetime import datetime
import glob
import json
import logging
import os
import pdb
import pickle
import re
import sys
import time

# Related third party imports
from astropy.coordinates import SkyCoord
from astropy.io import ascii, fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from itertools import chain
import numpy as np
import pandas as pd
from scipy.spatial import KDTree

from bokeh.layouts import row, column, gridplot
from bokeh.plotting import figure, output_file, save, show
from bokeh.models import ColumnDataSource, Label, CDSView, Div
from bokeh.models.tools import HoverTool

# Local application imports
from drizzlepac import util, wcs_functions
from drizzlepac.haputils.graph_utils import HAPFigure, build_tooltips
from drizzlepac.haputils.pandas_utils import PandasDFReader, get_pandas_data
from drizzlepac.haputils.pipeline_graphics import build_vector_plot
from stsci.tools import logutil
from stwcs import wcsutil
from stwcs.wcsutil import HSTWCS

__taskname__ = 'mvm_quality_graphics'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)


# ----------------------------------------------------------------------------------------------------------------------
# Module level variables

# ====================================================================================
# GAIA plots: number of GAIA sources, mean distance to neighbors, centroid/offset/std
# of sources in field
# ====================================================================================


def build_mvm_plots(data_source, output_basename='', display_plot=False, log_level=logutil.logging.INFO):
    """Create all the plots for the results generated by these comparisons

    Parameters
    ----------
    data_source : str
        Filename for master data file which contains all the results.  This will
        typically be an HSF5 file generated by the JSON harvester.

    output_base_filename : str
        Base name for the HMTL file generated by Bokeh.

    display_plot : bool, optional
        Option to display the plot to the screen
        Default: False

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the
        .log file. Default value is 20, or 'info'.

    Returns
    -------
    Nothing.
    """
    log.setLevel(log_level)
    if output_basename == '':
        output_basename = "mvm_qa"
    else:
        output_basename = "{}_mvm_qa".format(output_basename)
    build_overlap_crossmatch_plots(data_source)

# ------------------------------------------------------------------------------------------------------------


def build_overlap_crossmatch_plots(data_source):
    details_column_basenames = ["overlap_region_size",
                                "reference_catalog_name",
                                "comparison_catalog_name",
                                "total_reference_catalog_size",
                                "total_comparison_catalog_size",
                                "number_of_reference_sources_available_for_crossmatch",
                                "number_of_comparison_sources_available_for_crossmatch",
                                "number_of_crossmatched_sources"]

    difference_column_basenames = ["X-axis_differences",
                                   "Y-axis_differences",
                                   "On-sky_separation_(X-Y)",
                                   "On-sky_separation_(RA-Dec)"]

    stats_column_basenames = ["Non-clipped_minimum",
                              "Non-clipped_maximum",
                              "Non-clipped_mean",
                              "Non-clipped_median",
                              "Non-clipped_standard_deviation",
                              "3x3-sigma_clipped_mean",
                              "3x3-sigma_clipped_median",
                              "3x3-sigma_clipped_standard_deviation",
                              "Percent_of_all_diff_values_within_1-sigma_of_the_clipped_mean",
                              "Percent_of_all_diff_values_within_2-sigma_of_the_clipped_mean",
                              "Percent_of_all_diff_values_within_3-sigma_of_the_clipped_mean"]
    data_table_column_basename = "crossmatched_reference_X,_Y,_RA,_Dec,_and_crossmatched_comparison_-_reference_difference_values"
    data_table_colnames = ["X-Skycell",
                           "Y-Skycell",
                           "RA",
                           "DEC",
                           "X-axis differences",
                           "Y-axis differences",
                           "On-sky separation (X-Y)",
                           "On-sky separation (RA-Dec)"]

    n_layers_colname = 'gen_info.number of overlap regions present'
    num_layers = get_pandas_data(data_source, [n_layers_colname])[n_layers_colname]
    # columns_to_retrieve = []
    # for layer_ctr in range(1, max(num_layers.values)+1):
    #     column_basename = "overlap_region_#{}".format(layer_ctr)
    #     # add "overlap details" columns
    #     for details_colname in details_column_basenames:
    #         columns_to_retrieve.append("{}_details.{}".format(column_basename, details_colname))
    #     # add stats columns for each difference type
    #     for diff_type in difference_column_basenames:
    #         for stats_colname in stats_column_basenames:
    #             columns_to_retrieve.append("{}_{}.{}".format(column_basename, diff_type, stats_colname))
    #     # add all the data table columns
    #     for data_table_colname in data_table_colnames:
    #         columns_to_retrieve.append("{}_{}.{}".format(column_basename, data_table_column_basename, data_table_colname))
    # overlap_dataframe = get_pandas_data(data_source, columns_to_retrieve)

    # overlap_dataframe =overlap_dataframe[overlap_dataframe.
    # create blank dataframe restacked_overlap_dataframe"
    restacked_overlap_dataframe = pd.DataFrame()

    for df_indexname, layer_val in zip(num_layers.index.values, num_layers.values):
        for layer_ctr in range(1, layer_val + 1):
            columns_to_retrieve = []
            column_basename = "overlap_region_#{}".format(layer_ctr)
            # add "overlap details" columns
            for details_colname in details_column_basenames:
                columns_to_retrieve.append("{}_details.{}".format(column_basename, details_colname))
            # add stats columns for each difference type
            for diff_type in difference_column_basenames:
                for stats_colname in stats_column_basenames:
                    columns_to_retrieve.append("{}_{}.{}".format(column_basename, diff_type, stats_colname))
            # add all the data table columns
            for data_table_colname in data_table_colnames:
                columns_to_retrieve.append(
                    "{}_{}.{}".format(column_basename, data_table_column_basename, data_table_colname))
            overlap_dataframe = get_pandas_data(data_source, columns_to_retrieve)
            overlap_dataframe = overlap_dataframe[overlap_dataframe['gen_info.dataframe_index'] == df_indexname]
            overlap_dataframe['gen_info.dataframe_index'] = "{}_{}".format(overlap_dataframe['gen_info.dataframe_index'],layer_ctr)
            col_rename_dict = {}
            for colname in columns_to_retrieve:
                col_rename_dict[colname] = colname.replace(column_basename, "overlap_region")
            overlap_dataframe = overlap_dataframe.rename(columns=col_rename_dict)
            restacked_overlap_dataframe.append(overlap_dataframe)

    # use following line of code to filter out all other DF rows besides the one we want:
    #
    # add overlap number to then end of dataframe index
    # rename portions of column titles with "overlap_region_#N" to simply "overlap_region"
    # append this updated one-row dataframe to restacked_overlap_dataframe

    pdb.set_trace()


# ------------------------------------------------------------------------------------------------------------


if __name__ == "__main__":
    """Simple command-line interface. That is all.
    """
    # process command-line inputs with argparse
    parser = argparse.ArgumentParser(description='Generate MVM quality assessment plots based on information'
                                                 ' stored in the user-specified .h5 file.')
    parser.add_argument('input_filename',
                        help='Input .h5 file produced by diagnostic_json_harvester.py that holds the '
                             'information to plot.')
    parser.add_argument('-d', '--display_plot', required=False, action='store_true',
                        help='If specified, plots will be automatically opened in the default web browser as '
                             'they are generated. Otherwise, .html plot files will be generated but not '
                             'opened.')
    parser.add_argument('-l', '--log_level', required=False, default='info',
                        choices=["critical", "error", "warning", "info", "debug", "notset"],
                        help='The desired level of verboseness in the log statements displayed on the screen '
                             'and written to the .log file.')
    parser.add_argument('-o', '--output_basename', required=False, default='',
                        help='Base name for the HMTL file generated by Bokeh.')
    user_args = parser.parse_args()

    # verify that input file exists
    if not os.path.exists(user_args.input_filename):
        err_msg = "File {} doesn't exist.".format(user_args.input_filename)
        log.critical(err_msg)
        sys.exit(err_msg)

    # set up logging
    log_dict = {"critical": logutil.logging.CRITICAL,
                "error": logutil.logging.ERROR,
                "warning": logutil.logging.WARNING,
                "info": logutil.logging.INFO,
                "debug": logutil.logging.DEBUG,
                "notset": logutil.logging.NOTSET}
    log_level = log_dict[user_args.log_level]
    log.setLevel(log_level)

    # execute plot generation!
    build_mvm_plots(user_args.input_filename, output_basename=user_args.output_basename,
                    display_plot=user_args.display_plot, log_level=log_level)
