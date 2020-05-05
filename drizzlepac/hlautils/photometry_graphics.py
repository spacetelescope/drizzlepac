# This routine is designed specifically to support plotting the
# SVM data quality output.  As such, it expects the input Pandas
# Dataframe to have a particular structure and the data to contain
# specific measurements.
# This is NOT a general plotting utility.
#
# SVM Photometry file has the following data values where
# differences are POINT - SEGMENT.
# Delta_MagAp1-Mean Difference: -0.04963999999999993
# Delta_MagAp1-Standard Deviation: 0.0869423778526139
# Delta_MagAp1-Median Difference: -0.0259999999999998
# Delta_MagAp2-Mean Difference: 0.009571428571428555
# Delta_MagAp2-Standard Deviation: 0.02693804431083513
# Delta_MagAp2-Median Difference: 0.006000000000000227
# ISSUES: 
# Too many columns from header cluttering up the file.
# Hierachical information not tightly coupled to values

import qa_harvest as qa

import pandas as pd

import array
import csv
from bokeh.layouts import row
from bokeh.plotting import figure, output_file, show
from bokeh.models import ColumnDataSource, CustomJS
from bokeh.models.tools import HoverTool
from bokeh import events

def load_data_to_arrays(csv_filename):
    # col_list = df.columns.tolist()

    df = pd.read_csv(csv_filename)
    key_MagAp1_mean = 'Delta_MagAp1-Mean Difference'
    key_MagAp2_mean = 'Delta_MagAp2-Mean Difference'
    key_MagAp1_std = 'Delta_MagAp1-Standard Deviation'
    key_MagAp2_std = 'Delta_MagAp2-Standard Deviation'
    key_MagAp1_median = 'Delta_MagAp1-Median Difference'
    key_MagAp2_median = 'Delta_MagAp2-Median Difference'

    d_MagAp1_mean = df[key_MagAp1_mean]
    d_MagAp2_mean = df[key_MagAp2_mean]
    d_MagAp1_std = df[key_MagAp1_std]
    d_MagAp2_std = df[key_MagAp2_std]
    d_MagAp1_median = df[key_MagAp1_median]
    d_MagAp2_median = df[key_MagAp2_median]

    return (d_MagAp1_mean, d_MagAp1_std, d_MagAp1_median), (d_MagAp2_mean, d_MagAp2_std, d_MagAp2_median)


# Generate the actual plot for the "svm_graphic_type" data.
def generate_graphic(ap1, ap2):

    # Set the output file immediately as advised by Bokeh.
    output_file('proto.html')

    # Generate a general index array
    x_index = array.array('i', (i for i in range(0,len(ap1[0])))) 
    num_of_datasets = len(x_index)
    print('Number of datasets: {}'.format(num_of_datasets))

    # Define a figure object
    p1 = figure()

    # Add the glyphs
    p1.circle(x_index, ap1[0], size=10, color='green', legend_label='Mean Aperture 1 Magnitude Differences')
    p1.triangle(x_index, ap1[2], size=10, color='blue', legend_label='Median Aperture 1 Magnitude Differences')
    info_text = 'Number of datasets: ' + str(num_of_datasets)
    p1.legend.click_policy='hide'

    p1.title.text = 'Differences Point - Segment Ap1 Magnitude'
    p1.xaxis.axis_label = 'Index     ' + info_text
    p1.yaxis.axis_label = 'Difference (Magnitudes)'

    p2 = figure()
    p2.circle(x_index, ap2[0], size=10, color='green', legend_label='Mean Aperture 2 Magnitude Differences')
    p2.triangle(x_index, ap2[2], size=10, color='blue', legend_label='Median Aperture 2 Magnitude Differences')
    p2.legend.click_policy='hide'

    p2.title.text = 'Differences Point - Segment Ap2 Magnitude'
    p2.xaxis.axis_label = 'Index     ' + info_text
    p2.yaxis.axis_label = 'Difference (Magnitudes)'

    # Display!
    show(row(p1, p2))


def photometry_graphics_driver(csv_filename):
   
    # Get the data back as two tuples which contain mean, std, and median
    ap1, ap2 = load_data_to_arrays(csv_filename)

    generate_graphic(ap1, ap2)


def generate_photometry_df(input_directory, output_csv_filename='svm_photometry_stats.csv'):
    """Generate a CSV file containing the Pandas dataframe for statistics.

    This routine searches the named directory for a specified JSON file 
    and extracts the statistics to generate a Pandas dataframe.  The dataframe
    is written out to a CSV file for safe keeping.  In this way, if there is
    an error generating the graphics, the summary file still exists.

    Parameters
    ==========
    input_directory: str
    Directory to search for *.json files
    """

    df = qa.make_master_df(input_directory, pattern='*photometry.json')
    df.to_csv(output_csv_filename)
