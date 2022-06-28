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



def build_mvm_plots(data_source, output_basename='mvm_qa', display_plot=False, log_level=logutil.logging.INFO):
    """Create all the plots for the results generated by these comparisons

    Parameters
    ----------
    data_source : str
        name of the .h5 file produced by diagnostic_json_harvester.py that holds the information to plot.

    output_basename : str
        text string that will be used as the basis for all .html files generated by this script. Unless
        explicitly specified, the default value is 'mvm_qa'. If a value is explicitly specified, the text
        string '_mvm_qa' will be automatically appended to the end of the user-specified value.

    display_plot : bool, optional
        If set to Boolean 'True', plots .html files will be automatically opened in the default web browser
        as they are generated. Unless explicitly specified, the default value is Boolean 'False', meaning
        that plots will only be written to output .html files but not displayed on-screen.

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

    # Generate overlap crossmatch plots
    try:
        overlap_xmatch_filename = build_overlap_crossmatch_plots(data_source,
                                                                 display_plot,
                                                                 output_basename=output_basename,
                                                                 log_level=log_level)
    except Exception:
        log.warning("Overlap crossmatch plot generation encountered a problem.")
        log.exception("message")
        log.warning("Continuing to next plot...")
#     -      -     -      -     -      -     -      -     -      -     -      -     -      -     -      -


# ------------------------------------------------------------------------------------------------------------


def build_overlap_crossmatch_plots(data_source, display_plot=False, output_basename='mvm_qa', log_level=logutil.logging.INFO):
    """retrieve required data and reformat the dataframe in preperation for the generation of overlap
    crossmatch plots.

    Parameters
    ----------
    data_source : str
        name of the .h5 file produced by diagnostic_json_harvester.py that holds the information to plot.

    display_plot : bool, optional.
        If set to Boolean 'True', plots .html files will be automatically opened in the default web browser
        as they are generated. Unless explicitly specified, the default value is Boolean 'False', meaning
        that plots will only be written to output .html files but not displayed on-screen.

    output_basename : str, optional.
        text string that will be used as the basis for all .html files generated by this script. Unless
        explicitly specified, the default value is 'mvm_qa'.

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the
        .log file. Default value is 20, or 'info'.

    Returns
    -------
    output : str
        Name of HTML file where the plot was saved.
    """
    log.setLevel(log_level)
    # lists of column titles that will be used.
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

    # retrieve relevant data and "restack" and rename dataframe so that all the information for each overlap
    # region is stored in discrete rows
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
            overlap_dataframe['gen_info.dataframe_index'] = "{}_overlap_region_{}".format(overlap_dataframe['gen_info.dataframe_index'][0], layer_ctr)
            col_rename_dict = {}
            for colname in columns_to_retrieve:
                col_rename_dict[colname] = colname.replace(column_basename, "overlap_region")
            overlap_dataframe = overlap_dataframe.rename(columns=col_rename_dict)
            restacked_overlap_dataframe = restacked_overlap_dataframe.append(overlap_dataframe)
    # Sort columns alphabetically to make it more human-friendly
    restacked_overlap_dataframe = restacked_overlap_dataframe[overlap_dataframe.columns.sort_values()]
    # optionally write dataframe to .csv file.
    if log_level == logutil.logging.DEBUG:
        output_csv_filename = "testout.csv"
        if os.path.exists(output_csv_filename):
            os.remove(output_csv_filename)
        restacked_overlap_dataframe.to_csv(output_csv_filename)
        log.debug("Wrote restacked dataframe to csv file {}".format(output_csv_filename))

    # generate plots!
    output_file = generate_overlap_crossmatch_graphics(restacked_overlap_dataframe,
                                                       display_plot=display_plot,
                                                       output_basename=output_basename,
                                                       log_level=log_level)

    return(output_file)


# ------------------------------------------------------------------------------------------------------------

def generate_overlap_crossmatch_graphics(dataframe, display_plot=False, output_basename='mvm_qa', log_level=logutil.logging.INFO):
    """Generate plots to statistically quantify the quality of the alignment of crossmatched sources in
    regions where observations from different proposal/vists overlap.

    Parameters
    ----------
    dataframe : pandas dataframe
        dataframe containing results from the overlap crossmatch(s) to plot.

    display_plot : bool, optional.
        If set to Boolean 'True', plots .html files will be automatically opened in the default web browser
        as they are generated. Unless explicitly specified, the default value is Boolean 'False', meaning
        that plots will only be written to output .html files but not displayed on-screen.

    output_basename : str, optional.
        text string that will be used as the basis for all .html files generated by this script. Unless
        explicitly specified, the default value is 'mvm_qa'.

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the
        .log file. Default value is 20, or 'info'.

    Returns
    -------
    output : str
        Name of HTML file where the plot was saved.
    """
    log.setLevel(log_level)
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
    xmatch_cds = ColumnDataSource(dataframe)
    # generate plots of x vs. y components of various stat. measures for each difference
    output_basename = "{}_overlap_crossmatch_stats_plots".format(output_basename)
    if not output_basename.endswith('.html'):
        output = output_basename + '.html'
    else:
        output = output_basename
    # Set the output file immediately as advised by Bokeh.
    output_file(output)

    # Generate the graphic-specific tooltips - be mindful of the default tooltips defined in graph_utils.py
    # TODO: Add more detail to hovertips

    # Define the graphics
    # Create title text at the top of the html file
    html_title_text = Div(text="""<h1>Distribution characteristics of crossmatched sources identified in regions of overlapping observations in the MVM product</h1>""")

    # Scatter plots! #TODO: add more detail to hover_tips
    p0 = HAPFigure(title='Minimum difference value',
                   x_label='X minimum difference (pixels)',
                   y_label='Y minimum difference (pixels)')#,hover_tips=gaia_tips)
    p0.build_glyph('circle',
                   x='overlap_region_X-axis_differences.Non-clipped_minimum',
                   y='overlap_region_Y-axis_differences.Non-clipped_minimum',
                   sourceCDS=xmatch_cds,
                   glyph_color='colormap',
                   legend_group='inst_det')

    p1 = HAPFigure(title='Maximum difference value',
                   x_label='X maximum difference (pixels)',
                   y_label='Y maximum difference (pixels)')#,hover_tips=gaia_tips)
    p1.build_glyph('circle',
                   x='overlap_region_X-axis_differences.Non-clipped_maximum',
                   y='overlap_region_Y-axis_differences.Non-clipped_maximum',
                   sourceCDS=xmatch_cds,
                   glyph_color='colormap',
                   legend_group='inst_det')
    row1 = row(p0.fig, p1.fig)

    p2 = HAPFigure(title='Median difference value',
                   x_label='X median difference (pixels)',
                   y_label='Y median difference (pixels)')#,hover_tips=gaia_tips)
    p2.build_glyph('circle',
                   x='overlap_region_X-axis_differences.Non-clipped_median',
                   y='overlap_region_Y-axis_differences.Non-clipped_median',
                   sourceCDS=xmatch_cds,
                   glyph_color='colormap',
                   legend_group='inst_det')

    p3 = HAPFigure(title='Mean difference value',
                   x_label='X mean difference (pixels)',
                   y_label='Y mean difference (pixels)')#,hover_tips=gaia_tips)
    p3.build_glyph('circle',
                   x='overlap_region_X-axis_differences.Non-clipped_mean',
                   y='overlap_region_Y-axis_differences.Non-clipped_mean',
                   sourceCDS=xmatch_cds,
                   glyph_color='colormap',
                   legend_group='inst_det')
    row2 = row(p2.fig, p3.fig)

    p4 = HAPFigure(title='Difference standard deviation value',
                   x_label='X standard deviation difference (pixels)',
                   y_label='Y standard deviation difference (pixels)')#,hover_tips=gaia_tips)
    p4.build_glyph('circle',
                   x='overlap_region_X-axis_differences.Non-clipped_standard_deviation',
                   y='overlap_region_Y-axis_differences.Non-clipped_standard_deviation',
                   sourceCDS=xmatch_cds,
                   glyph_color='colormap',
                   legend_group='inst_det')

    p5 = HAPFigure(title='3x3 sigma-clippped median difference value',
                   x_label='X sigma-clipped median difference (pixels)',
                   y_label='Y sigma-clipped median difference (pixels)')#,hover_tips=gaia_tips)
    p5.build_glyph('circle',
                   x='overlap_region_X-axis_differences.3x3-sigma_clipped_median',
                   y='overlap_region_Y-axis_differences.3x3-sigma_clipped_median',
                   sourceCDS=xmatch_cds,
                   glyph_color='colormap',
                   legend_group='inst_det')
    row3 = row(p4.fig, p5.fig)

    p6 = HAPFigure(title='3x3 sigma-clipped mean difference value',
                   x_label='X sigma-clipped mean difference (pixels)',
                   y_label='Y sigma-clipped mean difference (pixels)')#,hover_tips=gaia_tips)
    p6.build_glyph('circle',
                   x='overlap_region_X-axis_differences.3x3-sigma_clipped_mean',
                   y='overlap_region_Y-axis_differences.3x3-sigma_clipped_mean',
                   sourceCDS=xmatch_cds,
                   glyph_color='colormap',
                   legend_group='inst_det')

    p7 = HAPFigure(title='3x3 sigma-clippped Difference standard deviation value',
                   x_label='X sigma-clipped difference standard_deviation (pixels)',
                   y_label='Y sigma-clipped difference standard_deviation (pixels)')#,hover_tips=gaia_tips)
    p7.build_glyph('circle',
                   x='overlap_region_X-axis_differences.3x3-sigma_clipped_standard_deviation',
                   y='overlap_region_Y-axis_differences.3x3-sigma_clipped_standard_deviation',
                   sourceCDS=xmatch_cds,
                   glyph_color='colormap',
                   legend_group='inst_det')
    row4 = row(p6.fig, p7.fig)
    # Display and save
    row_list = [row1, row2, row3, row4]
    if display_plot:
        show(column(row_list))
    # Just save
    else:
        save(column(row_list))
    log.info("Output HTML graphic file {} has been written.\n".format(output))

    # generate quad resid (x vs. dx, y vs. dx, x vs. dy, y vs. dy) plots for each DF row
    for index, row in dataframe.iterrows():
        pdb.set_trace()
    # Set the output file immediately as advised by Bokeh.


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
