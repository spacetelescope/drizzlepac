# This routine is designed specifically to support plotting of the
# SVM data quality output, report_wcs().  As such, it expects
# the input Pandas Dataframe to have a particular structure and the
# data to contain specific measurements.
#
# To use:
# 1) Start a Python session
# 2) import photometry_graphics as pg
# 3) You will need to provide a fully qualified pathname to the directory
#    above all of the output directories which contain the photometry
#    JSON output files.
# 4) Invoke the reader/graphics routine
#    >>> pg.photometry_graphics_driver(storage_filename)
#        This will generate a Bokeh plot to the browser, as
#        well as generate an HTML file.

# Standard library imports
import argparse
import logging
import os
import re
import sys
import pandas as pd

from bokeh.layouts import gridplot, row
from bokeh.plotting import figure, output_file, show, save
from bokeh.models import ColumnDataSource, Label
#from bokeh.models.tools import HoverTool, WheelZoomTool, PanTool, ZoomInTool
from bokeh.models.tools import HoverTool

# Local application imports
from drizzlepac.haputils.pandas_utils import PandasDFReader
from drizzlepac.haputils.graph_utils import HAPFigure
from stsci.tools import logutil

# Observations taken after Oct 2017 will most likely NOT have any apriori
# solutions since they were observed using GAIA-based coordinates for the
# guide stars by default.  In a sense, the default IDC_<rootname> WCS is
# the apriori solution.
# WCS columns

# TO DO:  Need to keep gen_info and header until fix made in pandas_utils.py
""" 
WCS_COLUMNS = {'gen_info.instrument': 'gen_info.instrument',
               'gen_info.detector': 'gen_info.detector',
               'gen_info.inst_det': 'gen_info.inst_det',
               'gen_info.filter': 'gen_info.filter',
               'gen_info.dataset': 'gen_info.dataset',
               'gen_info.proposal_id': 'gen_info.proposal_id',
               'gen_info.imgname': 'gen_info.imgname',
               'header.ASN_ID': 'header.ASN_ID',
""" 
WCS_COLUMNS = {'PrimaryWCS.primary_wcsname': 'prim_wcsname',
               'PrimaryWCS.crpix1': 'prim_crpix1',
               'PrimaryWCS.crpix2': 'prim_crpix2',
               'PrimaryWCS.crval1': 'prim_crval1',
               'PrimaryWCS.crval2': 'prim_crval2',
               'PrimaryWCS.scale': 'prim_scale',
               'PrimaryWCS.orientation': 'prim_orient',
               'PrimaryWCS.exposure': 'prim_exponame',
               'AlternateWCS_default.alt_wcsname': 'alt_default_wcsname',
               'AlternateWCS_default.crpix1': 'alt_default_crpix1',
               'AlternateWCS_default.crpix2': 'alt_default_crpix2',
               'AlternateWCS_default.crval1': 'alt_default_crval1',
               'AlternateWCS_default.crval2': 'alt_default_crval2',
               'AlternateWCS_default.scale': 'alt_default_scale',
               'AlternateWCS_default.orientation': 'alt_default_orient',
               'AlternateWCS_default.exposure': 'alt_default_exponame',
               'DeltaWCS_default.delta_wcsname': 'del_default_wcsname',
               'DeltaWCS_default.d_crpix1': 'del_default_crpix1',
               'DeltaWCS_default.d_crpix2': 'del_default_crpix2',
               'DeltaWCS_default.d_crval1': 'del_default_crval1',
               'DeltaWCS_default.d_crval2': 'del_default_crval2',
               'DeltaWCS_default.d_scale': 'del_default_scale',
               'DeltaWCS_default.d_orientation': 'del_default_orient',
               'DeltaWCS_default.exposure': 'del_default_exponame',
               'AlternateWCS_apriori.alt_wcsname': 'alt_apriori_wcsname',
               'AlternateWCS_apriori.crpix1': 'alt_apriori_crpix1',
               'AlternateWCS_apriori.crpix2': 'alt_apriori_crpix2',
               'AlternateWCS_apriori.crval1': 'alt_apriori_crval1',
               'AlternateWCS_apriori.crval2': 'alt_apriori_crval2',
               'AlternateWCS_apriori.scale': 'alt_apriori_scale',
               'AlternateWCS_apriori.orientation': 'alt_apriori_orient',
               'AlternateWCS_apriori.exposure': 'alt_apriori_exponame',
               'DeltaWCS_apriori.delta_wcsname': 'del_apriori_wcsname',
               'DeltaWCS_apriori.d_crpix1': 'del_apriori_crpix1',
               'DeltaWCS_apriori.d_crpix2': 'del_apriori_crpix2',
               'DeltaWCS_apriori.d_crval1': 'del_apriori_crval1',
               'DeltaWCS_apriori.d_crval2': 'del_apriori_crval2',
               'DeltaWCS_apriori.d_scale': 'del_apriori_scale',
               'DeltaWCS_apriori.d_orientation': 'del_apriori_orient',
               'DeltaWCS_apriori.exposure': 'del_apriori_exponame',
               'AlternateWCS_aposteriori.alt_wcsname': 'alt_aposteriori_wcsname',
               'AlternateWCS_aposteriori.crpix1': 'alt_aposteriori_crpix1',
               'AlternateWCS_aposteriori.crpix2': 'alt_aposteriori_crpix2',
               'AlternateWCS_aposteriori.crval1': 'alt_aposteriori_crval1',
               'AlternateWCS_aposteriori.crval2': 'alt_aposteriori_crval2',
               'AlternateWCS_aposteriori.scale': 'alt_aposteriori_scale',
               'AlternateWCS_aposteriori.orientation': 'alt_aposteriori_orient',
               'AlternateWCS_aposteriori.exposure': 'alt_aposteriori_exponame',
               'DeltaWCS_aposteriori.delta_wcsname': 'del_aposteriori_wcsname',
               'DeltaWCS_aposteriori.d_crpix1': 'del_aposteriori_crpix1',
               'DeltaWCS_aposteriori.d_crpix2': 'del_aposteriori_crpix2',
               'DeltaWCS_aposteriori.d_crval1': 'del_aposteriori_crval1',
               'DeltaWCS_aposteriori.d_crval2': 'del_aposteriori_crval2',
               'DeltaWCS_aposteriori.d_scale': 'del_aposteriori_scale',
               'DeltaWCS_aposteriori.d_orientation': 'del_aposteriori_orient',
               'DeltaWCS_aposteriori.exposure': 'del_aposteriori_exponame'}

__taskname__ = 'wcs_graphics'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)
# ----------------------------------------------------------------------------------------------------------------------


def compute_global_stats(value_array):
    """Compute the mean of the input array values.

    Parameters
    ==========
    value_array : array
        Values for which a mean is to be computed

    Returns
    =======
    mean_value : float
        Computed mean of the input array
    """

    mean_value = value_array.mean()
    return mean_value


def get_data(storage_filename):
    """Load the harvested data, stored in a storage file, into local arrays.

    Parameters
    ==========
    storage_filename : str
        Name of the storage file for the Pandas dataframe created by the harvester.

    Returns
    =======
    wcs_dataDF : Pandas dataframe
        Dataframe which is a subset of the input Pandas dataframe which
        consists of only the requested columns and rows where all of the requested
        columns do not contain NaNs.
    """

    # Instantiate a Pandas Dataframe Reader (lazy instantiation)
    df_handle = PandasDFReader(storage_filename, log_level=logutil.logging.NOTSET)

    # Get the relevant column data, eliminating all rows which have NaNs
    # in any of the relevant columns.
    wcs_column_type = ['apriori', 'aposteriori']
    windex = -1
    wcs_dataDF = pd.DataFrame()
    while wcs_dataDF.empty:
        wcs_dataDF = df_handle.get_columns_HDF5(WCS_COLUMNS.keys())

        # If no dataframe were returned, there was a KeyError because columns were
        # not present in the original dataframe versus the columns contained NaNs.
        entries_to_remove = []
        if wcs_dataDF.empty and windex < len(wcs_column_type) - 1:
            log.info("Columns are missing from the data file {}. Remove some requested column names\n"
                     "from the list and try again.\n".format(storage_filename))

            # Identify missing columns and remove them from a copy of the original dictionary
            windex += 1
            for key in WCS_COLUMNS:
                if re.match(r'.+{}.+'.format(wcs_column_type[windex]), key):
                    entries_to_remove.append(key)

            for key in entries_to_remove:
                WCS_COLUMNS.pop(key, None)

        elif wcs_dataDF.empty:
            log.critical("Critical columns not found in storage Pandas dataframe: {}.\n".format(storage_filename))
            sys.exit(1)

    log.info("wcs_graphics. WCS data has been retrieved from the storage Pandas dataframe: {}.\n".format(storage_filename))

    # Rename the columns to abbreviated text as the graph titles further
    # document the information.
    for old_col_name, new_col_name in WCS_COLUMNS.items():
        wcs_dataDF.rename(columns={old_col_name: new_col_name}, inplace=True)

    return wcs_dataDF


# Generate the actual plot for the "svm_graphic_type" data.
def generate_graphic(wcs_dataDF, output_base_filename, display_plot, log_level):
    """Generate the graphics associated with this particular type of data.

    Parameters
    ==========
    wcs_dataDF : Pandas dataframe
        A subset of the input Dataframe consisting only of the WCS_COLUMNS
        and Aperture 2

    output_base_filename : str
        Base name for the HMTL file generated by Bokeh.

    display_plot : bool
        Option to display the plot to the screen
        Default: False

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the .log file.
        Default: 20 or 'info'.

    Returns
    =======
    Nothing

    HTML file is generated and written to the current directory.
    """

    # Set the output file immediately as advised by Bokeh.
    output_file(output_base_filename + '.html')

    # Setup the source of the data to be plotted so the axis variables can be
    # referenced by column name in the Pandas dataframe
    sourceCDS = ColumnDataSource(wcs_dataDF)
    num_of_datasets = len(wcs_dataDF.index)
    print('Number of datasets: {}'.format(num_of_datasets))

    # TO DO: Alternate no longer works.  $name is not interpreted. See original code.
    wcs_tips = [("Primary WCS", "@prim_wcsname"),
                ("Alternate WCS", "@$name")]

    # Instantiate the figure objects 
    text = 'WCS Component Differences (Active - Alternate)'
    figure1 = HAPFigure(title = text,
                        x_label = 'Delta CRPIX1 (pixels)',
                        y_label = 'Delta CRPIX2 (pixels)',
                        grid_line_color = 'gainsboro',
                        hover_tips = wcs_tips)
    figure2 = HAPFigure(title = text,
                        x_label = 'Delta CRVAL1 (pixels)',
                        y_label = 'Delta CRVAL2 (pixels)',
                        grid_line_color = 'gainsboro',
                        hover_tips = wcs_tips)
    figure3 = HAPFigure(title = text,
                        x_label = 'Delta Scale (pixels/arcseconds)',
                        y_label = 'Delta Orientation (degrees)',
                        grid_line_color = 'gainsboro',
                        hover_tips = wcs_tips)

    # Figure 1 delta_crpix1 vs delta_crpix2
    # Figure 2 delta_crval1 vs delta_crval2
    # Figure 3 delta_scale vs delta_orientation
    wcs_type_names = ['default', 'apriori', 'aposteriori']
    wcs_components = []
    alt_wcs_names = []
    for wcs_type in wcs_type_names:
        wcs_components.append('del_{}_wcsname'.format(wcs_type))
        alt_wcs_names.append('alt_{}_wcsname'.format(wcs_type))

    # TO DO: Need way to specify any "shape" glyph
    # There are three distinct figures in this graphic layout, each figure can have up to
    # three datasets plotted with the circular glyph
    wcs_type_colors = ['blue', 'green', 'purple']
    wcs_type_colors = ['orange', 'orange', 'orange']
    for i, wcs_component in enumerate(wcs_components):
        if wcs_component in wcs_dataDF.columns:
            slist = wcs_component.rsplit('_')
            figure1.build_glyph('triangle', x='del_{}_crpix1'.format(slist[1]),
                                       y='del_{}_crpix2'.format(slist[1]),
                                       sourceCDS=sourceCDS,
                                       marker_color=wcs_type_colors[i],
                                       legend_label='Delta(Active-'+wcs_type_names[i].capitalize()+')',
                                       name=alt_wcs_names[i])

            """

            figure2.build_circle_glyph(x='del_{}_crval1'.format(slist[1]),
                                       y='del_{}_crval2'.format(slist[1]),
                                       sourceCDS=sourceCDS,
                                       marker_color=wcs_type_colors[i],
                                       legend_label='Delta(Active-'+wcs_type_names[i].capitalize()+')',
                                       name=alt_wcs_names[i])

            figure3.build_circle_glyph(x='del_{}_scale'.format(slist[1]),
                                       y='del_{}_orient'.format(slist[1]),
                                       sourceCDS=sourceCDS,
                                       marker_color=wcs_type_colors[i],
                                       legend_label='Delta(Active-'+wcs_type_names[i].capitalize()+')',
                                       name=alt_wcs_names[i])
            """

    """
    # Create the the HTML output and optionally display the plot
    grid = gridplot([[figure1.fig, figure2.fig], [figure3.fig, None]], plot_width=500, plot_height=500)
    if display_plot:
        show(grid)
    log.info("WCS_graphics. Output HTML graphic file {} has been written.\n".format(output_base_filename + ".html"))
    """
    show(figure1.fig)


def wcs_graphics_driver(storage_filename, output_base_filename='wcs_graphics', display_plot=False, log_level=logutil.logging.INFO):
    """Driver to load the data from the storage file and generate the graphics.

    Parameters
    ==========
    storage_filename : str
        Name of the storage file for the Pandas dataframe created by the harvester.

    output_base_filename : str
        Base name for the HMTL file generated by Bokeh.

    display_plot : bool
        Option to display the plot in addition to writing out the file.

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the .log file.
        Default value is 20, or 'info'.

    Returns
    =======
    Nothing
    """

    # Retrieve the relevant dataframe columns
    log.info('Retrieve Pandas dataframe from file {}.\n'.format(storage_filename))
    wcs_dataDF = get_data(storage_filename)

    # Generate the WCS graphic
    generate_graphic(wcs_dataDF, output_base_filename, display_plot, log_level)


# ======================================================================================================================


if __name__ == "__main__":
    # Process command-line inputs with argparse
    parser = argparse.ArgumentParser(description='Read the harvested Pandas dataframe stored as and HDF5 file.')

    parser.add_argument('harvester_filename', help='File which holds the Pandas dataframe in an HDF5 file.')
    parser.add_argument('-o', '--output_base_filename', required=False, default="wcs_graphics",
                        help='Name of the output base filename (filename without the extension) for the  '
                             'HTML file generated by Bokeh.')
    parser.add_argument('-d', '--display_plot', required=False, default=False,
                        help='Option to display the plot to the screen in addition to writing the output file.')
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
    if not os.path.exists(user_args.harvester_filename):
        err_msg = "Harvester HDF5 File {} does not exist.".format(user_args.harvester_filename)
        log.critical(err_msg)
        raise Exception(err_msg)

    log.info("Harvester file {} found.\n".format(user_args.harvester_filename))

    wcs_graphics_driver(user_args.harvester_filename,
                        output_base_filename=user_args.output_base_filename,
                        display_plot=user_args.display_plot,
                        log_level=log_level)
