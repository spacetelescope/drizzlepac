#!/usr/bin/env python
import argparse
import pdb

from astropy.stats import sigma_clipped_stats
from astropy.table import Table, vstack
import matplotlib.pyplot as plt
import numpy as np

from drizzlepac.devutils.comparison_tools import compare_sourcelists
from drizzlepac.haputils.diagnostic_utils import read_json_file
"""
1: read in ci_fwhm file(s)
2: smash all input datasets into a single dataset
3: determine min, max CI values
4: set up CI value bins based on dataset min/max values and specified bin size
5: for each bin, compute resistant mean FWHM value, sigma value
6: plot bin CI value vs. mean FWHM values
"""
# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
def get_data(input_files):
    """read in data from one or more files. In the case that there are multiple files, the data from all files
    will be combined into a single CI vs FWHM dataset.

    Parameters
    ----------
    input_files : list
        list of CI vs. FWHM .csv files to process

    Returns
    -------
    data_table: : numpy.ndarray
        A 2 x n sized numpy array. Column 1: CI values; Column 2: FWHM values
    """
    data_table = Table()
    for filename in input_files:
        if filename.endswith(".csv"):
            data_table = vstack([data_table,Table.read(filename,format='ascii.csv')])
        if filename.endswith(".json"):
            json_data = read_json_file(filename)
            data_table = vstack([data_table, json_data['data']['CI_FWHM']])
    return data_table

# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-

def run(input_files,bin_size,ci_limits,plot_title):
    """Main calling subroutine

    Parameters
    ----------
    input_files : list
        list of CI vs. FWHM .csv files to process

    bin_size : float
        Size of the bin to use for CI values.

    Returns
    -------
    Nothing.
    """
    data_table = get_data(input_files)
    processed_data_table = process_data(data_table,bin_size)
    plot_data(processed_data_table,bin_size,ci_limits=ci_limits,plot_title=plot_title)
# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-

def plot_data(data_table,bin_size,ci_limits=[-1.0, -1.0],plot_title=None):
    fig = plt.figure(figsize=(11, 8.5))
    ax1 = fig.add_subplot(111)
    plt.scatter(data_table['CI'], data_table['FWHM'], marker=".", s=10, color="blue")
    if ci_limits[0] != -1.0 and ci_limits[1] != -1.0:
        ax1.axvline(x=ci_limits[0], color='k', linestyle='--')
        ax1.axvline(x=ci_limits[1], color='k', linestyle='--')
    full_plot_title = "Binned CI vs Mean FWHM value; CI binsize = {}".format(bin_size)
    if plot_title:
        full_plot_title = "{}\n{}".format(plot_title,full_plot_title)
    ax1.set_title(full_plot_title)
    ax1.set_xlabel("CI")
    ax1.set_ylabel("Mean FWHM Value (Pixels)")
    ax1.grid(True)


    plt.show()
# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-

def process_data(data_table,bin_size):
    """Bin CI values according to specified bin_size, Compute mean FWHM for each bin from only individual
    fWHM values in the bin

    Parameters
    ----------
    data_table: : numpy.ndarray
        A 2 x n sized numpy array. Column 1: CI values; Column 2: FWHM values

    bin_size : float
        Size of the bin to use for CI values.

    Returns
    -------
    processed_data_table : numpy.ndarray
        A 2 x n sized numpy array. Column 1: Binned CI values; Column 2: Mean FWHM values
    """
    processed_data_table = Table(names=("CI","FWHM"))
    start_value = compare_sourcelists.round2ArbatraryBase(min(data_table['CI']), "down",bin_size)
    stop_value = compare_sourcelists.round2ArbatraryBase(max(data_table['CI']), "down",bin_size)
    print(start_value, stop_value)
    while start_value <= stop_value:
        bin_upper_limit = start_value + bin_size
        ix = np.where((data_table['CI'] >= start_value) & (data_table['CI'] < bin_upper_limit))
        ix = ix[0]
        if len(ix >0):
            clippedStats = sigma_clipped_stats([data_table['FWHM'][ix]])
            processed_data_table.add_row([start_value,clippedStats[0]])
        start_value += bin_size
    return processed_data_table

# =======================================================================================================================
if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description='Plot binned CI values vs. mean FWHM values')
    PARSER.add_argument('input_files', nargs='+',help='one or more space-separated ci vs fwhm csv files')
    PARSER.add_argument('-b', '--bin_size', required=False, default=0.01, type=float,
                        help = "Size of the bin to use for CI values. Default value is 0.01")
    PARSER.add_argument('-c', '--ci_limits',nargs=2, required=False, type=float, default=[-1.0, -1.0],
                        help = "Optional values for ci_lower_limit and ci_upper_limit to be plotted as vertical lines in scatter plot")
    PARSER.add_argument('-t', '--plot_title', required=False, default=None, help="Optional plot title")
    ARGS = PARSER.parse_args()
    run(ARGS.input_files,ARGS.bin_size,ARGS.ci_limits,ARGS.plot_title)