# This routine is designed specifically to support plotting the
# SVM data quality output.  As such, it expects the input Pandas
# Dataframe to have a particular structure and the data to contain
# specific measurements.
#
# This is NOT a general plotting utility and is a prototype.
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

from drizzlepac.haputils.pandas_utils import PandasDFReader
from stsci.tools import logutil

from bokeh.layouts import row
from bokeh.plotting import figure, output_file, show
from bokeh.models import ColumnDataSource, Label
from bokeh.models.tools import HoverTool
import csv
import glob
import json
import logging
import os
import sys

PHOT_COLUMNS = ['Statistics_MagAp1.Delta_MagAp1.Mean', 'Statistics_MagAp1.Delta_MagAp1.StdDev', 
                'Statistics_MagAp1.Delta_MagAp1.Median', 'Statistics_MagAp2.Delta_MagAp2.Mean', 
                'Statistics_MagAp2.Delta_MagAp2.StdDev', 'Statistics_MagAp2.Delta_MagAp2.Median']
DPHOT_COLUMNS = {'Ap1MeanDiff': 'Statistics_MagAp1.Delta_MagAp1.Mean', 
                 'Ap1StdDev': 'Statistics_MagAp1.Delta_MagAp1.StdDev',
                 'Ap1MedianDiff': 'Statistics_MagAp1.Delta_MagAp1.Median',
                 'Ap2MeanDiff': 'Statistics_MagAp2.Delta_MagAp2.Mean',
                 'Ap2StdDev': 'Statistics_MagAp2.Delta_MagAp2.StdDev',
                 'Ap2MedianDiff': 'Statistics_MagAp2.Delta_MagAp2.Median'}
"""
# FIX MDD: Fix hovertool 
# FIX MDD: Where should this code live?  Join with svm_quality_analysis.py?
"""


MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)

def compute_global_stats(value_array):
    """Compute the mean of the input array values.

    Parameters
    ==========
    value_array: array
    Values for which a mean is to be computed

    Returns
    =======
    mean_value: float
    Computed mean of the input array
    """

    mean_value = value_array.mean()
    return mean_value


def get_data(storage_filename):
    """Load the harvested data, stored in a storage file, into local arrays.

    Parameters
    ==========
    storage_filename: str
    Name of the storage file for the Pandas dataframe created by the harvester.

    Returns
    =======
    phot_data: Pandas dataframe
    Dataframe which is a subset of the input Pandas dataframe which
    consists of only the requested columns and rows where all of the requested 
    columns do not contain NaNs.

    (mean_dMagAp1_mean, mean_dMagAp1_median): tuple
    Aperture 1 mean of means and mean of medians

    (mean_dMagAp2_mean, mean_dMagAp2_median): tuple
    Aperture 2 mean of means and mean of medians
    """
    
    # Instantiate a Pandas Dataframe Reader (lazy instantiation)
    df_handle = PandasDFReader(storage_filename, log_level=logutil.logging.NOTSET)

    # In this particular case, the names of the desired columns do not
    # have to be further manipulated, for example, to add dataset specific
    # names.
    # 
    # Get the relevant column data, eliminating all rows which have NaNs
    # in any of the relevant columns.
    phot_data = df_handle.get_columns_HDF5(PHOT_COLUMNS)

    # Generate a general index array and add it to the dataframe
    x_index = list(range(0, len(phot_data.index)))
    phot_data['x_index'] = x_index
    x_index.clear()

    mean_dMagAp1mean = compute_global_stats(phot_data[DPHOT_COLUMNS['Ap1MeanDiff']])
    mean_dMagAp2mean = compute_global_stats(phot_data[DPHOT_COLUMNS['Ap2MeanDiff']])

    mean_dMagAp1median = compute_global_stats(phot_data[DPHOT_COLUMNS['Ap1MedianDiff']])
    mean_dMagAp2median = compute_global_stats(phot_data[DPHOT_COLUMNS['Ap2MedianDiff']])

    return phot_data, (mean_dMagAp1mean, mean_dMagAp1median), \
        (mean_dMagAp2mean, mean_dMagAp2median)


# Generate the actual plot for the "svm_graphic_type" data.
def generate_graphic(phot_data, stat_Ap1, stat_Ap2):
    """Generate the graphics associated with this particular type of data.

    Parameters
    ==========
    phot_data: Pandas dataframe
    Dataframe consisting of the Magnitude statistics of mean, std, and median for Aperture 1
    and Aperture 2

    stat_Ap1: tuple
    Tuple containing the average of the means and the average of the medians for Aperture 1

    stat_Ap2: tuple
    Tuple containing the average of the means and the average of the medians for Aperture 2
    """

    # Set the output file immediately as advised by Bokeh.
    output_file('photometry_grahics.html')

    # Setup the source of the data to be plotted so the axis variables can be
    # referenced by column name in the Pandas dataframe
    sourceDF = ColumnDataSource(phot_data)
    num_of_datasets = len(phot_data.index)
    print('Number of datasets: {}'.format(num_of_datasets))

    # Define a figure object
    p1 = figure()

    # Add the glyphs
    p1.circle(x='x_index', y=PHOT_COLUMNS[0], source=sourceDF, fill_alpha=0.5, line_alpha=0.5, size=10, color='green',
              legend_label='Mean Aperture 1 Magnitude Differences')
    p1.cross(x='x_index', y=PHOT_COLUMNS[2], source=sourceDF, size=10, color='blue', legend_label='Median Aperture 1 Magnitude Differences')
    info_text = 'Number of datasets: ' + str(num_of_datasets)
    p1.legend.click_policy = 'hide'
    
    p1.title.text = 'Differences Point - Segment Ap1 Magnitude'
    p1.xaxis.axis_label = 'Index     ' + info_text
    p1.yaxis.axis_label = 'Difference (Magnitudes)'

    # hover_p1 = HoverTool()
    # col0 = DPHOT_COLUMNS['C1']
    # print(col0)
    # hover_p1.tooltips=[('Mean', @'DPHOT_COLUMNS["C1"]')]
    # p1.add_tools(hover_p1)

    stat_text = ('Mean deltaMAp1_mean: {:6.2f}     Mean deltaMAp1_median: {:6.2f}'.format(stat_Ap1[0], stat_Ap1[1]))
    stat_label = Label(x=20, y=20, x_units='screen', y_units='screen', text=stat_text)
    p1.add_layout(stat_label)

    p2 = figure()
    p2.circle(x='x_index', y=PHOT_COLUMNS[3], source=sourceDF, fill_alpha=0.5, line_alpha=0.5, size=10, color='green',
              legend_label='Mean Aperture 2 Magnitude Differences')
    p2.cross(x='x_index', y=PHOT_COLUMNS[5], source=sourceDF, size=10, color='blue', legend_label='Median Aperture 2 Magnitude Differences')
    p2.legend.click_policy = 'hide'

    p2.title.text = 'Differences Point - Segment Ap2 Magnitude'
    p2.xaxis.axis_label = 'Index     ' + info_text
    p2.yaxis.axis_label = 'Difference (Magnitudes)'

    stat_text = ('Mean deltaMAp2_mean: {:6.2f}     Mean deltaMAp2_median: {:6.2f}'.format(stat_Ap2[0], stat_Ap2[1]))
    stat_label = Label(x=20, y=20, x_units='screen', y_units='screen', text=stat_text)
    p2.add_layout(stat_label)

    # Display!
    show(row(p1, p2))


def photometry_graphics_driver(storage_filename, log_level=logutil.logging.INFO):
    """Driver to load the data from the storage file and generate the graphics.

    Parameters
    ==========
    storage_filename: str
    Name of the storage file for the Pandas dataframe created by the harvester.

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the .log file.
        Default value is 20, or 'info'.
    """

    # Retrieve the relevant dataframe and statistics (mean of
    # means and mean of medians) for Aperture 1 and Aperture 2
    log.info('Retrieve Pandas dataframe from file {}.\n'.format(storage_filename))
    phot_data, stat_Ap1, stat_Ap2 = get_data(storage_filename)

    # Generate the photometric graphic
    generate_graphic(phot_data, stat_Ap1, stat_Ap2)


