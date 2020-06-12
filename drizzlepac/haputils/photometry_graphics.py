# This routine is designed specifically to support plotting the
# SVM data quality output.  As such, it expects the input Pandas
# Dataframe to have a particular structure and the data to contain
# specific measurements.
#
# This is NOT a general plotting utility and is a prototype.
# To use:
# 1) Start a Python session
# 2) import photometry_graphics as pg
# 3) You need to provide a fully qualified pathname to the directory
#    above all of the output directories which contain the photometry
#    JSON output files.
# 4) In the Python session, invoke the harvester
#    >>> generate_photometry_df(input_directory, output_csv_filename='svm_photometry_stats.csv')
#        This will generate an output CSV file where a Pandas
#        dataframe is stored.
# 5) Once the CSV file is generated, invoke the graphics routine
#    >>> pg.photometry_graphics_driver(CSV_filename)
#        This will generate a Bokeh plot to the browsers, as
#        well as generate an HTML file.
#
# SVM Photometry file has the following data values where
# differences are (POINT - SEGMENT).
# Delta_MagAp1-Mean Difference: -0.04963999999999993
# Delta_MagAp1-Standard Deviation: 0.0869423778526139
# Delta_MagAp1-Median Difference: -0.0259999999999998
# Delta_MagAp2-Mean Difference: 0.009571428571428555
# Delta_MagAp2-Standard Deviation: 0.02693804431083513
# Delta_MagAp2-Median Difference: 0.006000000000000227
# ISSUES:
# Too many columns from header in JSON cluttering up dataframe.
# Hierachical information in JSON not tightly coupled to values once in dataframe.
# To use hover, MUST have all data columns coming from ColumnDataSource.
# ColumnDataSource needs to use columns names, so use succinct names.
# Think about information useful for hover.
####  Fix the docs above.

from drizzlepac.haputils.pandas_utils import PandasDFReader_CSV
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
PHOT_COLUMNS = ["AP1_Mean","AP1_StdDev","AP1_Median","AP2_Mean","AP2_StdDev","AP2_Median"]
DPHOT_COLUMNS = {"C1": "AP1_Mean","C2": "AP1_StdDev","C3": "AP1_Median","C4": "AP2_Mean","C5": "AP2_StdDev","C6": "AP2_Median"}

"""
PHOT_COLUMNS = ["Photometry_Statistics_MagAp1.Delta_MagAp1=Point_MagAp1_-_Segment_MagAp1.Mean Difference",
                "Photometry_Statistics_MagAp1.Delta_MagAp1=Point_MagAp1_-_Segment_MagAp1.Standard Deviation",
                "Photometry_Statistics_MagAp1.Delta_MagAp1=Point_MagAp1_-_Segment_MagAp1.Median Difference",
                "Photometry_Statistics_MagAp2.Delta_MagAp2=Point_MagAp2_-_Segment_MagAp2.Mean Difference",
                "Photometry_Statistics_MagAp2.Delta_MagAp2=Point_MagAp2_-_Segment_MagAp2.Standard Deviation",
                "Photometry_Statistics_MagAp2.Delta_MagAp2=Point_MagAp2_-_Segment_MagAp2.Median Difference"]
DPHOT_COLUMNS = {"Ap1MeanDiff": "Photometry_Statistics_MagAp1.Delta_MagAp1=Point_MagAp1_-_Segment_MagAp1.Mean Difference",
                 "Ap1StdDev": "Photometry_Statistics_MagAp1.Delta_MagAp1=Point_MagAp1_-_Segment_MagAp1.Standard Deviation",
                 "Ap1MedianDiff": "Photometry_Statistics_MagAp1.Delta_MagAp1=Point_MagAp1_-_Segment_MagAp1.Median Difference",
                 "Ap2MeanDiff": "Photometry_Statistics_MagAp2.Delta_MagAp2=Point_MagAp2_-_Segment_MagAp2.Mean Difference",
                 "Ap2StdDev": "Photometry_Statistics_MagAp2.Delta_MagAp2=Point_MagAp2_-_Segment_MagAp2.Standard Deviation",
                 "Ap2MedianDiff": "Photometry_Statistics_MagAp2.Delta_MagAp2=Point_MagAp2_-_Segment_MagAp2.Median Difference"}
"""

# FIX MDD: Create an enumeration and use a map function
# FIX MDD: Fix names again svm_quality_analysis.py in JSON generation
# FIX MDD: Need to use the real harvester for data and metadata
# FIX MDD: Fix hovertool 
# FIX MDD: Where should this code live?  Join with svm_quality_analysis.py?


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


def get_data(csv_filename):
    """Load the harvested data, stored in a CSV file, into local arrays.

    Parameters
    ==========
    csv_filename: str
    Name of the CSV file created by the harvester.

    Returns
    =======
    phot_data: Pandas dataframe
    Dataframe which is a subset of the input Pandas dataframe written out as
    a CSV file.  The subset dataframe consists of only the requested columns
    and rows where all of the requested columns did not contain NaNs.

    (mean_dMagAp1_mean, mean_dMagAp1_median): tuple
    Aperture 1 mean of means and mean of medians

    (mean_dMagAp2_mean, mean_dMagAp2_median): tuple
    Aperture 2 mean of means and mean of medians
    """
    
    # Instantiate a Pandas Dataframe Reader (lazy instantiation)
    # df_handle = PandasDFReader_CSV("svm_qa_dataframe.csv")
    df_handle = PandasDFReader_CSV(csv_filename, log_level=logutil.logging.NOTSET)

    # In this particular case, the names of the desired columns do not
    # have to be further manipulated, for example, to add dataset specific
    # names.
    # 
    # Get the relevant column data, eliminating all rows which have NaNs
    # in any of the relevant columns.
    phot_data = df_handle.get_columns(PHOT_COLUMNS)

    # Generate a general index array and add it to the dataframe
    x_index = list(range(0, len(phot_data.index)))
    phot_data['x_index'] = x_index
    x_index.clear()

    # mean_dMagAp1mean = compute_global_stats(phot_data[PHOTC[0]])
    # mean_dMagAp2mean = compute_global_stats(phot_data[PHOTC[3]])

    # mean_dMagAp1median = compute_global_stats(phot_data[PHOT[2]])
    # mean_dMagAp2median = compute_global_stats(phot_data[PHOT[5]])
    mean_dMagAp1mean = 10.0
    mean_dMagAp2mean = 11.0
    mean_dMagAp1median = 12.0
    mean_dMagAp2median = 13.0

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


def photometry_graphics_driver(csv_filename, log_level=logutil.logging.INFO):
    """Driver to load the data from the CSV file and generate the graphics.

    Parameters
    ==========
    csv_filename: str
    Name of the CSV file created by the harvester.

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the .log file.
        Default value is 20, or 'info'.
    """

    # Retrieve the relevant dataframe and statistics (mean of
    # means and mean of medians) for Aperture 1 and Aperture 2
    log.info('Retrieve the dfjslf')
    phot_data, stat_Ap1, stat_Ap2 = get_data(csv_filename)

    # Generate the photometric graphic
    generate_graphic(phot_data, stat_Ap1, stat_Ap2)


