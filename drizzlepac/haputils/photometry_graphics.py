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

import array
from bokeh.layouts import row
from bokeh.plotting import figure, output_file, show
from bokeh.models import ColumnDataSource, Label
import csv
import glob
import json
import os
import pandas as pd


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


def load_data_to_arrays(csv_filename):
    """Load the harvested data, stored in a CSV file, into local arrays.

    Parameters
    ==========
    csv_filename: str
    Name of the CSV file created by generate_photometry_df.

    Returns
    =======
    (d_MagAp1_mean, d_MagAp1_std, d_MagAp1_median): tuple
    Aperture 1 statistics

    (d_MagAp2_mean, d_MagAp2_std, d_MagAp2_median): tuple
    Aperture 2 statistics

    (mean_dMagAp1_mean, mean_dMagAp1_median): tuple
    Aperture 1 mean of means and mean of medians

    (mean_dMagAp2_mean, mean_dMagAp2_median): tuple
    Aperture 2 mean of means and mean of medians
    """

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

    mean_dMagAp1_mean = compute_global_stats(d_MagAp1_mean)
    mean_dMagAp2_mean = compute_global_stats(d_MagAp2_mean)

    mean_dMagAp1_median = compute_global_stats(d_MagAp1_median)
    mean_dMagAp2_median = compute_global_stats(d_MagAp2_median)

    return (d_MagAp1_mean, d_MagAp1_std, d_MagAp1_median), \
        (d_MagAp2_mean, d_MagAp2_std, d_MagAp2_median), \
        (mean_dMagAp1_mean, mean_dMagAp1_median), \
        (mean_dMagAp2_mean, mean_dMagAp2_median)


# Generate the actual plot for the "svm_graphic_type" data.
def generate_graphic(ap1, ap2, means_ap1, means_ap2):
    """Generate the graphics associated with this particular type of data.

    Parameters
    ==========
    ap1: tuple
    Tuple containing the mean, std, and median statistics for Aperture 1

    ap2: tuple
    Tuple containing the mean, std, and median statistics for Aperture 2

    means_ap1: tuple
    Tuple containing the average of the means and the average of the medians for Aperture 1

    means_ap2: tuple
    Tuple containing the average of the means and the average of the medians for Aperture 2
    """

    # Set the output file immediately as advised by Bokeh.
    output_file('photometry_grahics.html')

    # Generate a general index array
    x_index = array.array('i', (i for i in range(0, len(ap1[0]))))
    num_of_datasets = len(x_index)
    print('Number of datasets: {}'.format(num_of_datasets))

    # Define a figure object
    p1 = figure()

    # Add the glyphs
    p1.circle(x_index, ap1[0], fill_alpha=0.5, line_alpha=0.5, size=10, color='green',
              legend_label='Mean Aperture 1 Magnitude Differences')
    p1.cross(x_index, ap1[2], size=10, color='blue', legend_label='Median Aperture 1 Magnitude Differences')
    info_text = 'Number of datasets: ' + str(num_of_datasets)
    p1.legend.click_policy = 'hide'

    p1.title.text = 'Differences Point - Segment Ap1 Magnitude'
    p1.xaxis.axis_label = 'Index     ' + info_text
    p1.yaxis.axis_label = 'Difference (Magnitudes)'

    stat_text = ('Mean deltaMAp1_mean: {:6.2f}     Mean deltaMAp1_median: {:6.2f}'.format(means_ap1[0], means_ap1[1]))
    stat_label = Label(x=20, y=20, x_units='screen', y_units='screen', text=stat_text)
    p1.add_layout(stat_label)

    p2 = figure()
    p2.circle(x_index, ap2[0], fill_alpha=0.5, line_alpha=0.5, size=10, color='green',
              legend_label='Mean Aperture 2 Magnitude Differences')
    p2.cross(x_index, ap2[2], size=10, color='blue', legend_label='Median Aperture 2 Magnitude Differences')
    p2.legend.click_policy = 'hide'

    p2.title.text = 'Differences Point - Segment Ap2 Magnitude'
    p2.xaxis.axis_label = 'Index     ' + info_text
    p2.yaxis.axis_label = 'Difference (Magnitudes)'

    stat_text = ('Mean deltaMAp2_mean: {:6.2f}     Mean deltaMAp2_median: {:6.2f}'.format(means_ap2[0], means_ap2[1]))
    stat_label = Label(x=20, y=20, x_units='screen', y_units='screen', text=stat_text)
    p2.add_layout(stat_label)

    # Display!
    show(row(p1, p2))


def photometry_graphics_driver(csv_filename):
    """Driver to load the data from the CSV file and generate the graphics.

    Parameters
    ==========
    csv_filename: str
    Name of the CSV file created by generate_photometry_df.
    """

    # Get the data back as two tuples which contain mean, std, and median
    ap1, ap2, means_ap1, means_ap2 = load_data_to_arrays(csv_filename)

    generate_graphic(ap1, ap2, means_ap1, means_ap2)


def generate_photometry_df(input_directory, output_csv_filename='svm_photometry_stats.csv'):
    """Generate a CSV file containing the Pandas dataframe for statistics.

    This routine searches the named directory for a specified JSON file
    and extracts the statistics to generate a Pandas dataframe.  The dataframe
    is written out to a CSV file for safe keeping.  In this way, if there were to
    be an error generating the graphics, the summary file still exists and the
    harvesting does not have to be redone.

    Parameters
    ==========
    input_directory: str
    Directory to search for *.json files

    output_csv_filename: str
    Filename for the output CSV file
    """

    df = make_master_df(input_directory, pattern='*photometry.json')
    df.to_csv(output_csv_filename)


def make_dataset_df(dirname, pattern='*.json'):
    """Convert dir full of JSON files into a DataFrame"""

    jpatt = os.path.join(dirname, pattern)
    hdr = None

    pdtabs = []
    for jfilename in sorted(glob.glob(jpatt)):
        with open(jfilename) as jfile:
            resids = json.load(jfile)
        pdindx = None
        if hdr is None:
            hdr = resids['header']
        rootname = hdr['FILENAME'].replace('.fits', '')
        k = resids['data'].keys()
        for key in k:
            dat = resids['data'][key]['data']

            det = dat['detector']
            filtname = dat['filter_name']
            del dat['detector']
            del dat['filter_name']
            if pdindx is None:
                pdindx = '-'.join([rootname, det, filtname])
            for dk in dat.keys():
                for di in dat[dk].keys():
                    hdr.update(dict([('-'.join([dk.split(" ")[2], di]), dat[dk][di])]))

        pdtabs.append(pd.DataFrame(hdr, index=[pdindx]))
    if len(pdtabs) == 0:
        allpd = None
    else:
        allpd = pd.concat(pdtabs)
    return allpd


def make_master_df(dirname, pattern='*.json', num=None):
    dirs = sorted(glob.glob(os.path.join(dirname, '*')))
    allpd = None
    for d in dirs[:num]:
        pdtab = make_dataset_df(d, pattern=pattern)
        if pdtab is not None:
            if allpd is None:
                allpd = pdtab
            else:
                allpd = allpd.append(pdtab)

    return allpd
