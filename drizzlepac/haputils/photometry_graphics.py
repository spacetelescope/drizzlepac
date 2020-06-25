# This routine is designed specifically to support plotting of the
# SVM data quality output, compare_photometry().  As such, it expects
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
#    By default, the routine will save the graphic to an HTML file
#    in the working directory.
#
# Alternative execution:
# bash> python photometry_graphics.py harvested_dataframe.h5 -o name_of_my_graphic_output
#    where HTML will be appended to the name_of_my_graphic_output by Bokeh.
#

# Standard library imports
import argparse
import logging
import os
import sys

from bokeh.layouts import gridplot, row
from bokeh.plotting import figure, output_file, show, save
from bokeh.models import ColumnDataSource, Label
from bokeh.models.tools import HoverTool

# Local application imports
from drizzlepac.haputils.pandas_utils import PandasDFReader
from stsci.tools import logutil

# Abbreviated dataframe column names
NEW_PHOT_COLUMNS = ['Ap1 Mean Differences',
                    'Ap1 Standard Deviation',
                    'Ap1 Median Differences',
                    'Ap2 Mean Differences',
                    'Ap2 Standard Deviation',
                    'Ap2 Median Differences',
                    'instrument',
                    'detector',
                    'filter',
                    'dataset']

# Dictionary mapping the original dataframe column name to their new names
PHOT_COLUMNS = {'Statistics_MagAp1.Delta_MagAp1.Mean': NEW_PHOT_COLUMNS[0],
                'Statistics_MagAp1.Delta_MagAp1.StdDev':  NEW_PHOT_COLUMNS[1],
                'Statistics_MagAp1.Delta_MagAp1.Median':  NEW_PHOT_COLUMNS[2],
                'Statistics_MagAp2.Delta_MagAp2.Mean':  NEW_PHOT_COLUMNS[3],
                'Statistics_MagAp2.Delta_MagAp2.StdDev': NEW_PHOT_COLUMNS[4],
                'Statistics_MagAp2.Delta_MagAp2.Median': NEW_PHOT_COLUMNS[5],
                'gen_info.instrument': NEW_PHOT_COLUMNS[6],
                'gen_info.detector': NEW_PHOT_COLUMNS[7],
                'gen_info.filter': NEW_PHOT_COLUMNS[8],
                'gen_info.dataset': NEW_PHOT_COLUMNS[9]}

__taskname__ = 'photometry_graphics'

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


def get_data(storage_filename, log_level):
    """Load the harvested data, stored in a storage file, into local arrays.

    Parameters
    ==========
    storage_filename : str
    Name of the storage file for the Pandas dataframe created by the harvester.

    log_level : int
        The desired level of verboseness in the log statements displayed on the screen and written to the .log file.

    Returns
    =======
    phot_dataDF : Pandas dataframe
    Dataframe which is a subset of the input Pandas dataframe which
    consists of only the requested columns and rows where all of the requested
    columns do not contain NaNs.

    (mean_dMagAp1_mean, mean_dMagAp1_median) : tuple
    Aperture 1 mean of means and mean of medians

    (mean_dMagAp2_mean, mean_dMagAp2_median) : tuple
    Aperture 2 mean of means and mean of medians
    """

    # Instantiate a Pandas Dataframe Reader (lazy instantiation)
    df_handle = PandasDFReader(storage_filename, log_level=logutil.logging.NOTSET)

    # Get the relevant column data, eliminating all rows which have NaNs
    # in any of the relevant columns.
    phot_dataDF = df_handle.get_columns_HDF5(PHOT_COLUMNS.keys())
    log.info("Photometry_graphics. Photometric data has been retrieved from the storage Pandas dataframe: {}.\n".format(storage_filename))

    # Rename the columns to abbreviated text as the graph titles further
    # document the information.
    for old_col_name, new_col_name in PHOT_COLUMNS.items():
        phot_dataDF.rename(columns={old_col_name: new_col_name}, inplace=True)

    # Generate a general index array and add it to the dataframe
    x_index = list(range(0, len(phot_dataDF.index)))
    phot_dataDF['x_index'] = x_index
    x_index.clear()

    mean_dMagAp1mean = compute_global_stats(phot_dataDF['Ap1 Mean Differences'])
    mean_dMagAp2mean = compute_global_stats(phot_dataDF['Ap1 Median Differences'])

    mean_dMagAp1median = compute_global_stats(phot_dataDF['Ap2 Mean Differences'])
    mean_dMagAp2median = compute_global_stats(phot_dataDF['Ap2 Median Differences'])

    return phot_dataDF, (mean_dMagAp1mean, mean_dMagAp1median), \
        (mean_dMagAp2mean, mean_dMagAp2median)


# Generate the actual plot for the "svm_graphic_type" data.
def generate_graphic(phot_dataDF, stat_Ap1, stat_Ap2, output_base_filename, log_level):
    """Generate the graphics associated with this particular type of data.

    Parameters
    ==========
    phot_dataDF : Pandas dataframe
    Dataframe consisting of the Magnitude statistics of mean, std, and median for Aperture 1
    and Aperture 2

    stat_Ap1 : tuple
    Tuple containing the average of the means and the average of the medians for Aperture 1

    stat_Ap2 : tuple
    Tuple containing the average of the means and the average of the medians for Aperture 2

    output_base_filename : str
    Base name for the HMTL file generated by Bokeh.

    log_level : int
        The desired level of verboseness in the log statements displayed on the screen and written to the .log file.

    Returns
    =======
    Nothing

    HTML file is generated and written to the current directory.
    """

    # Set the output file immediately as advised by Bokeh.
    output_file(output_base_filename + '.html')

    # Setup the source of the data to be plotted so the axis variables can be
    # referenced by column name in the Pandas dataframe
    sourceCDS = ColumnDataSource(phot_dataDF)
    num_of_datasets = len(phot_dataDF.index)
    log.info('Photometry_graphics. Number of datasets: {}.\n'.format(num_of_datasets))

    TOOLS = "box_zoom, wheel_zoom, box_select, lasso_select, pan, reset, save, help"

    # Define a figure object
    # Aperture 1
    p1 = figure(tools=TOOLS, toolbar_location="right")

    # Add the glyphs
    p1.circle(x='x_index', y=NEW_PHOT_COLUMNS[0], source=sourceCDS, fill_alpha=0.5, line_alpha=0.5, size=10, color='green',
              legend_label='Mean Aperture 1 Magnitude Differences', name="bee")
    p1.cross(x='x_index', y=NEW_PHOT_COLUMNS[2], source=sourceCDS, size=10, color='blue', legend_label='Median Aperture 1 Magnitude Differences')
    info_text = 'Number of datasets: ' + str(num_of_datasets)
    p1.legend.click_policy = 'hide'

    p1.title.text = 'Differences (Point - Segment) Aperture 1 Magnitude'
    p1.xaxis.axis_label = 'Index     ' + info_text
    p1.yaxis.axis_label = 'Difference (ABMag)'

    hover_p1 = HoverTool()
    hover_p1.tooltips = [("Instrument", "@instrument"),
                         ("Detector", "@detector"),
                         ("Filter", "@filter"),
                         ("Dataset", "@dataset"),
                         ("Mean", "@{Ap1 Mean Differences}"),
                         ("StdDev", "@{Ap1 Standard Deviation}"),
                         ("Median", "@{Ap1 Median Differences}")]
    p1.add_tools(hover_p1)

    stat_text = ('<DeltaMagAp1_Mean>: {:6.2f}     <DeltaMagAp1_Median>: {:6.2f}'.format(stat_Ap1[0], stat_Ap1[1]))
    stat_label = Label(x=20, y=20, x_units='screen', y_units='screen', text=stat_text)
    p1.add_layout(stat_label)

    # Aperture 2
    p2 = figure(tools=TOOLS, toolbar_location="right", x_range=p1.x_range, y_range=p1.y_range,)
    p2.circle(x='x_index', y=NEW_PHOT_COLUMNS[3], source=sourceCDS, fill_alpha=0.5, line_alpha=0.5, size=10, color='green',
              legend_label='Mean Aperture 2 Magnitude Differences')
    p2.cross(x='x_index', y=NEW_PHOT_COLUMNS[5], source=sourceCDS, size=10, color='blue', legend_label='Median Aperture 2 Magnitude Differences')
    p2.legend.click_policy = 'hide'

    p2.title.text = 'Differences (Point - Segment) Aperture 2 Magnitude'
    p2.xaxis.axis_label = 'Index     ' + info_text
    p2.yaxis.axis_label = 'Difference (ABMag)'

    hover_p2 = HoverTool()
    hover_p2.tooltips = [("Instrument", "@instrument"),
                         ("Detector", "@detector"),
                         ("Filter", "@filter"),
                         ("Dataset", "@dataset"),
                         ("Mean", "@{Ap2 Mean Differences}"),
                         ("StdDev", "@{Ap2 Standard Deviation}"),
                         ("Median", "@{Ap2 Median Differences}")]
    p2.add_tools(hover_p2)

    stat_text = ('<DeltaMagAp2_Mean>: {:6.2f}     <DeltaMagAp2_Median>: {:6.2f}'.format(stat_Ap2[0], stat_Ap2[1]))
    stat_label = Label(x=20, y=20, x_units='screen', y_units='screen', text=stat_text)
    p2.add_layout(stat_label)
    p = gridplot([[p1, p2]])

    # Create the the HTML output, but do not display at this time
    save(p)
    log.info("Photometry_graphics. Output HTML graphic file {} has been written.\n".format(output_base_filename + ".html"))


def photometry_graphics_driver(storage_filename, output_base_filename='photometry_graphics', log_level=logutil.logging.NOTSET):
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
    log.setLevel(log_level)

    # Retrieve the relevant dataframe and statistics (mean of
    # means and mean of medians) for Aperture 1 and Aperture 2
    log.info('Photometry_graphics. Retrieve Pandas dataframe from file {}.\n'.format(storage_filename))
    phot_dataDF, stat_Ap1, stat_Ap2 = get_data(storage_filename, log_level)

    # Generate the photometric graphic
    generate_graphic(phot_dataDF, stat_Ap1, stat_Ap2, output_base_filename, log_level)


# ======================================================================================================================


if __name__ == "__main__":
    # Process command-line inputs with argparse
    parser = argparse.ArgumentParser(description='Read the harvested Pandas dataframe stored as and HDF5 file.')

    parser.add_argument('harvester_filename', help='File which holds the Pandas dataframe in an HDF5 file.')
    parser.add_argument('-o', '--output_base_filename', required=False, default="photometry_graphics",
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
        err_msg = "Photometry_graphics. Harvester HDF5 File {} does not exist.".format(user_args.harvester_filename)
        log.critical(err_msg)
        raise Exception(err_msg)

    print(user_args.output_base_filename)

    photometry_graphics_driver(user_args.harvester_filename,
                               output_base_filename=user_args.output_base_filename,
                               log_level=log_level)
