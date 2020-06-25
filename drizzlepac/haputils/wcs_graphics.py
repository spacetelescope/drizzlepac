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
from bokeh.models.tools import HoverTool

# Local application imports
from drizzlepac.haputils.pandas_utils import PandasDFReader
from stsci.tools import logutil


# Observations taken after Oct 2017 will most likely NOT have any apriori
# solutions since they were observed using GAIA-based coordinates for the
# guide stars by default.  In a sense, the default IDC_<rootname> WCS is
# the apriori solution.
# WCS columns
WCS_COLUMNS = {'gen_info.instrument': 'Instrument',
               'gen_info.detector': 'Detector',
               'gen_info.filter': 'Filter',
               'gen_info.dataset': 'Dataset',
               'gen_info.proposal_id': 'Proposal ID',
               'PrimaryWCS.primary_wcsname': 'prim_wcsname',
               'PrimaryWCS.crpix1': 'prim_crpix1',
               'PrimaryWCS.crpix2': 'prim_crpix2',
               'PrimaryWCS.crval1': 'prim_crval1',
               'PrimaryWCS.crval2': 'prim_crval2',
               'PrimaryWCS.scale': 'prim_scale',
               'PrimaryWCS.orientation': 'prim_orient',
               'PrimaryWCS.exposure': 'prim_exponame',
               'AlternateWCS_default.alt_wcsname': 'alt_def_wcsname',
               'AlternateWCS_default.crpix1': 'alt_def_crpix1',
               'AlternateWCS_default.crpix2': 'alt_def_crpix2',
               'AlternateWCS_default.crval1': 'alt_def_crval1',
               'AlternateWCS_default.crval2': 'alt_def_crval2',
               'AlternateWCS_default.scale': 'alt_def_scale',
               'AlternateWCS_default.orientation': 'alt_def_orient',
               'AlternateWCS_default.exposure': 'alt_def_exponame',
               'DeltaWCS_default.delta_wcsname': 'del_def_wcsname',
               'DeltaWCS_default.d_crpix1': 'del_def_crpix1',
               'DeltaWCS_default.d_crpix2': 'del_def_crpix2',
               'DeltaWCS_default.d_crval1': 'del_def_crval1',
               'DeltaWCS_default.d_crval2': 'del_def_crval2',
               'DeltaWCS_default.d_scale': 'del_def_scale',
               'DeltaWCS_default.d_orientation': 'del_def_orient',
               'DeltaWCS_default.exposure': 'del_def_exponame',
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
               'AlternateWCS_aposteriori.alt_wcsname': 'alt_apost_wcsname',
               'AlternateWCS_aposteriori.crpix1': 'alt_apost_crpix1',
               'AlternateWCS_aposteriori.crpix2': 'alt_apost_crpix2',
               'AlternateWCS_aposteriori.crval1': 'alt_apost_crval1',
               'AlternateWCS_aposteriori.crval2': 'alt_apost_crval2',
               'AlternateWCS_aposteriori.scale': 'alt_apost_scale',
               'AlternateWCS_aposteriori.orientation': 'alt_apost_orient',
               'AlternateWCS_aposteriori.exposure': 'alt_apost_exponame',
               'DeltaWCS_aposteriori.delta_wcsname': 'del_apost_wcsname',
               'DeltaWCS_aposteriori.d_crpix1': 'del_apost_crpix1',
               'DeltaWCS_aposteriori.d_crpix2': 'del_apost_crpix2',
               'DeltaWCS_aposteriori.d_crval1': 'del_apost_crval1',
               'DeltaWCS_aposteriori.d_crval2': 'del_apost_crval2',
               'DeltaWCS_aposteriori.d_scale': 'del_apost_scale',
               'DeltaWCS_aposteriori.d_orientation': 'del_apost_orient',
               'DeltaWCS_aposteriori.exposure': 'del_apost_exponame'}

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

    if wcs_dataDF.empty:
        log.critical("Critical columns not found in storage Pandas dataframe: {}.\n".format(storage_filename))
        sys.exit(1)

    log.info("wcs_graphics. WCS data has been retrieved from the storage Pandas dataframe: {}.\n".format(storage_filename))

    # Rename the columns to abbreviated text as the graph titles further
    # document the information.
    for old_col_name, new_col_name in WCS_COLUMNS.items():
        wcs_dataDF.rename(columns={old_col_name: new_col_name}, inplace=True)

    return wcs_dataDF


# Generate the actual plot for the "svm_graphic_type" data.
def generate_graphic(wcs_dataDF, output_base_filename, log_level):
    """Generate the graphics associated with this particular type of data.

    Parameters
    ==========
    wcs_dataDF : Pandas dataframe
    A subset of the input Dataframe consisting only of the WCS_COLUMNS
    and Aperture 2

    output_base_filename : str
    Base name for the HMTL file generated by Bokeh.

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the .log file.
        Default value is 20, or 'info'.

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

    TOOLS = "box_zoom,wheel_zoom,box_select,lasso_select,reset,help"

    # Define a figure object
    # Aperture 1
    p1 = figure(tools=TOOLS, toolbar_location="right")
    p2 = figure(tools=TOOLS, toolbar_location="right")
    p3 = figure(tools=TOOLS, toolbar_location="right")

    # Get the existing/remaining columns of the WCS_COLUMNS to know which
    # glyphs to generate

    # Figure 1 delta_crpix1 vs delta_crpix2
    # Figure 2 delta_crval1 vs delta_crval2
    wcs_components = ['del_def_wcsname', 'del_apriori_wcsname', 'del_apost_wcsname']
    wcs_type_colors = ['blue', 'green', 'purple']
    wcs_type_names = ['default', 'apriori', 'aposteriori']
    alt_wcs_names = ['alt_def_wcsname', 'alt_apriori_wcsname', 'alt_apost_wcsname']
    for i, wcs_component in enumerate(wcs_components):
        if wcs_component in wcs_dataDF.columns:
            slist = wcs_component.rsplit('_')
            xcol_pix = 'del_{}_crpix1'.format(slist[1])
            ycol_pix = 'del_{}_crpix2'.format(slist[1])
            p1.circle(x=xcol_pix, y=ycol_pix, source=sourceCDS,
                      fill_alpha=0.5, line_alpha=0.5, size=10, color=wcs_type_colors[i],
                      legend_label='Delta(Active-'+wcs_type_names[i].capitalize()+')', name=alt_wcs_names[i])

            xcol_val = 'del_{}_crval1'.format(slist[1])
            ycol_val = 'del_{}_crval2'.format(slist[1])
            p2.circle(x=xcol_val, y=ycol_val, source=sourceCDS,
                      fill_alpha=0.5, line_alpha=0.5, size=10, color=wcs_type_colors[i],
                      legend_label='Delta(Active-'+wcs_type_names[i].capitalize()+')', name=alt_wcs_names[i])

            xcol_val = 'del_{}_scale'.format(slist[1])
            ycol_val = 'del_{}_orient'.format(slist[1])
            p3.circle(x=xcol_val, y=ycol_val, source=sourceCDS,
                      fill_alpha=0.5, line_alpha=0.5, size=10, color=wcs_type_colors[i],
                      legend_label='Delta(Active-'+wcs_type_names[i].capitalize()+')', name=alt_wcs_names[i])

    p1.legend.click_policy = 'hide'
    p1.title.text = 'WCS Component Differences (Active - Alternate)'
    p1.xaxis.axis_label = 'Delta CRPIX1 (pixels)'
    p1.yaxis.axis_label = 'Delta CRPIX2 (pixels)'

    hover_p1 = HoverTool()
    # Make the name of the graphic the actual WCS
    # This hover tool shows how to use a variable for one of the tooltips.  The glyph for the
    # alternate WCS data could represent one of three different options.  The glyph has the "name"
    # parameter set to the key of the appropriate option, and the key is interpreted as
    hover_p1.tooltips = [("Instrument", "@Instrument"),
                         ("Detector", "@Detector"),
                         ("Filter", "@Filter"),
                         ("Dataset", "@Dataset"),
                         ("Primary WCS", "@prim_wcsname"),
                         ("Alternate WCS", "@$name")]
    p1.add_tools(hover_p1)

    p2.legend.click_policy = 'hide'
    p2.title.text = 'WCS Component Differences (Active - Alternate)'
    p2.xaxis.axis_label = 'Delta CRVAL1 (pixels)'
    p2.yaxis.axis_label = 'Delta CRVAL2 (pixels)'
    p2.add_tools(hover_p1)

    p3.legend.click_policy = 'hide'
    p3.title.text = 'WCS Component Differences (Active - Alternate)'
    p3.xaxis.axis_label = 'Delta Scale (pixels/arcseconds)'
    p3.yaxis.axis_label = 'Delta Orientation (degrees)'
    p3.add_tools(hover_p1)

    # Create the the HTML output, but do not display at this time

    grid = gridplot([[p1, p2], [p3, None]], plot_width=500, plot_height=500)
    show(grid)
    log.info("WCS_graphics. Output HTML graphic file {} has been written.\n".format(output_base_filename + ".html"))


def wcs_graphics_driver(storage_filename, output_base_filename='wcs_graphics', log_level=logutil.logging.INFO):
    """Driver to load the data from the storage file and generate the graphics.

    Parameters
    ==========
    storage_filename : str
    Name of the storage file for the Pandas dataframe created by the harvester.

    output_base_filename : str
    Base name for the HMTL file generated by Bokeh.

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
    generate_graphic(wcs_dataDF, output_base_filename, log_level)


# ======================================================================================================================


if __name__ == "__main__":
    # Process command-line inputs with argparse
    parser = argparse.ArgumentParser(description='Read the harvested Pandas dataframe stored as and HDF5 file.')

    parser.add_argument('harvester_filename', help='File which holds the Pandas dataframe in an HDF5 file.')
    parser.add_argument('-o', '--output_base_filename', required=False, default="wcs_graphics",
                        help='Name of the output base filename (filename without the extension) for the  '
                             'HTML file generated by Bokeh.')
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
                        log_level=log_level)
