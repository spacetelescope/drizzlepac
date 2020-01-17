#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4 ai :
"""This script compares sources common to two user-specified sourcelists and displays various measures of their
differences. The differences can be calculated using one of three user-selectable methods:

* Absolute: Just a simple difference e.g. Comparision - Reference.
* Overall mean percent difference:  Percent difference based on the **overall mean** reference value e.g. ((Comparision-Reference)/mean_overall(Reference)) x 100.
* Dynamic percent difference: Percent difference values are calculated discretely for each pair of comparision/reference values e.g. ((Comparision[n]-Reference[n])/Reference[n]) x 100.

3x3-sigma clipped mean, median and standard deviation, and non-sigma clipped min and max values are
computed for the following:

* Image X position (in pixels)
* Image Y position (in pixels)
* Right Ascension (in arcseconds)
* Declination (in arcseconds)
* Flux (Outer Aperture)
* Magnitude (Inner Aperture)
* Magnitude (Outer Aperture)
* Magnitude Error (Inner Aperture)
* Magnitude Error (Outer Aperture)
* MSKY
* STDDEV
* CI (Concentration Index)

Absolute bit-wise comparisons are also performed for the following item:

* Flag Value

.. note::
    Not all sourcelist types compatible with this comparision script contain all of the columns listed above.
    Statistics (and optionally plots) can only be generated for columns common to both user-specified sourcelists.
    Thus, it is to be expected that not all runs will yield comparisons for all columns listed above.



Regression Testing
------------------
**All** of the following criteria must be met for the test to be declared "successful":

* All non-flag (linear) column comparisons: Less than 5% of all comparison - reference difference values are greater than 3 sigma from the sigma-clipped mean
* X position:
* Y position:
* Right Ascension:
* Declination:
* Flux (Inner Aperture):
* Flux (Outer Aperture):
* Magnitude (Inner Aperture):
* Magnitude (Outer Aperture):
* Magnitude error (Inner Aperture):
* Magnitude error (Outer Aperture):
* Flag Value: Less than 5% of all matched sources have differing flag values.

.. note::
    Sigma-clipped values for mean, sigma, and median are computed using the astropy.stats.sigma_clipped_stats() routine with three rounds of three-sigma clipping.


Plots
-----
Plots will be generated if the optional plot input is set to 'screen' or 'file'. If set to 'screen', the plots will be
simply displayed on-screen, and not written to disk. If set to 'file', a plot and statistical summary (on a second
page) is generated for each valid comparison. All plots and summaries saved to a single multiple-page pdf file. The
default name of the combined plot pdf file is 'combined_plot.pdf'. If a user-specified filename prefix string was
specified, it will be prepended to the default name.

The code generates up to three different types of plots:

* Difference histogram plots: These are the most commonly generated plot products are are produced for all column
comparisons except for the bit-wise comparison.
* Bit value barcharts
* X-Y absolute difference vector plot


Path
----
drizzlepac/drizzlepac/devutils/comparison_tools/compare_sourcelists.py

Dependencies
------------
- drizzlepac/drizzlepac/devutils/comparison_tools/starmatch_hist.py
- The `PyPDF2 <https://pypi.org/project/PyPDF2/>`_ python library

Inputs
------
* Required input
    1. *sourcelistNames*
        * A space-separated pair of sourcelists to compare. The first sourcelist is assumed to be the reference sourcelist that the second is being compared to.

* Optional inputs:
    #. -d *debugMode*
        * Perform additional match quality diagnostics
        * Input choices: "True" or "False"
        * Default value: False

    #. -i *imageNames*
        * A space-separated list of the fits images that were used to generate the input sourcelists. The first image corresponds to the first listed sourcelist, and so in. These will be used to improve the sourcelist alignment and matching.

    #. -m *diffMode*
        * How should the comp-ref difference be calculated? "absolute" is simply the straight comp-ref difference. "peman" is the mean percent difference ((C-R)/avg(R)) x 100. "pdynamic" is the dynamic percent difference ((C-R)/R) x 100
        * Input choices: "absolute", "pmean" or "pdynamic"
        * Default value: "pmean"

    #. -p *plotGen*
        * Generate plots?
        * Input choices: "True" or "False"
        * Default value: False

    #. -s *plotfile_prefix_string*
        * Text string that will prepend the plot files generated if plots are written to files ***REQUIRES the -p option set to 'file'***
        * Default value: blank text string ('')

    #. -v *verbose*
        * Display verbose output?
        * Input choices: "True" or "False"
        * Default value: True

Classes and Functions
---------------------
"""
import argparse
from datetime import datetime
import os
import pdb
import random
import sys

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from PyPDF2 import PdfFileMerger

from drizzlepac.devutils.comparison_tools import starmatch_hist
from stsci.tools import logutil
from stwcs import wcsutil

__taskname__ = 'compare_sourcelists'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout, format=SPLUNK_MSG_FORMAT,
                            datefmt=MSG_DATEFMT)


# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-


def check_match_quality(matched_x_list, matched_y_list):
    """Creates region file to check quality of source matching.

    Parameters
    ----------
    matched_x_list : list
        list of ref and comp x coords for matched sources

    matched_y_list : list
        list of ref and comp y coords for matched sources

    Returns
    -------
    Nothing.
    """
    out_filename = "match_check.reg"
    num_display = 5000  # Number of pairs to plot
    list_length = len(matched_x_list[0])
    if num_display > list_length:  # if the list of matched sources is smaller than num_display, just use all matched pairs, rather than a randomly selected subset.
        index_list = np.arange(list_length)
    else:
        index_list = random.sample(range(1, list_length), num_display)
    with open(out_filename, "w") as fout:
        for index_no in index_list:
            fout.write("circle({},{},10)  # color=green\n".format(matched_x_list[0][index_no], matched_y_list[0][
                index_no]))  # write ref source circle
            fout.write("circle({},{},10) # color=red\n".format(matched_x_list[1][index_no],
                                                               matched_y_list[1][index_no]))  # write comp source circle
            fout.write(
                "line({},{},{},{}) # color=blue\n".format(matched_x_list[0][index_no], matched_y_list[0][index_no],
                                                          matched_x_list[1][index_no],
                                                          matched_y_list[1][index_no]))  # write line connecting the two
    log.info("Wrote region file {}".format(out_filename))


# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
def computeFlagStats(matchedRA, max_diff, plotGen, plot_title, plotfile_prefix, verbose):
    """Compute and report statistics on the differences in flagging.

    Parameters
    ----------
    matchedRA : numpy.ndarray
        A 2 x len(refLines) sized numpy array. Column 1: matched reference values.
        Column 2: The corresponding matched comparision values

    max_diff : float
        Maximum allowable percentage of all matched sources with differences in their flag values for comparison to be
        declared a success

    plotGen : bool
        Generate plots and display them to the screen (True/False)?

    plot_title : str
        text string that will be used in plot title.

    plotfile_prefix : str
        text string that will prepend the plot files generated if plots are written to files

    verbose : bool
        display verbose output?

    Returns
    -------
    regTestStatus : str
        overall test result and statistics
    """
    pdf_file_list = []
    log.info(">>>>>> Comparision - reference sourcelist {} differences <<<<<<".format(plot_title))
    # compute overall percentage of matched sources with flagging differences
    flag_diff_list = list(matchedRA[0] - matchedRA[1])
    n_total = len(matchedRA[0])
    n_unchanged = flag_diff_list.count(0)
    n_changed = n_total - n_unchanged
    pct_changed = (float(n_changed) / float(n_total)) * 100.0

    # set up arrays to count stuff up
    bit_list = [0, 1, 2, 4, 8, 16, 32, 64, 128]
    refFlagBreakdown = np.zeros(9, dtype=int)
    compFlagBreakdown = np.zeros(9, dtype=int)
    unchangedFlagBreakdown = np.zeros(9, dtype=int)
    on_off_FlagFlips = np.zeros(9, dtype=int)
    off_on_FlagFlips = np.zeros(9, dtype=int)
    for refFlagVal, compFlagVal in zip(matchedRA[0], matchedRA[1]):
        # break down each flag value into component bit values, add values to totals
        refFlagRA = deconstruct_flag(refFlagVal)
        refFlagBreakdown += refFlagRA
        compFlagRA = deconstruct_flag(compFlagVal)
        compFlagBreakdown += compFlagRA
        # find differences in flagging, total up which bits were turned on, which were turned off.
        diffFlagRA = compFlagRA - refFlagRA
        if not np.array_equal(refFlagRA, compFlagRA):
            off_on_FlagFlips[np.where(diffFlagRA == 1)] += 1  # bits that are off in ref but on in comp
            on_off_FlagFlips[np.where(diffFlagRA == -1)] += 1  # bits that are on in ref but off in comp
            unchangedFlagBreakdown[np.where((refFlagRA == 1) & (
                        compFlagRA == 1))] += 1  # takes care of the case were comp and ref have differing bits, but also have additional bits that are unchanged.
        if np.array_equal(refFlagRA, compFlagRA):  # if there are absolutely no differences
            unchangedFlagBreakdown += refFlagRA
    regTestStatus = "OK     "
    if pct_changed > max_diff:
        regTestStatus = "FAILURE"
    log_output_string_list = []

    tot_str_len = len(str(n_total))
    uch_padding = tot_str_len - len(str(n_unchanged))
    ch_padding = tot_str_len - len(str(n_changed))
    if ((verbose is True) or (regTestStatus == "FAILURE")):
        # Generate result tables
        n = np.sum(refFlagBreakdown, dtype=float)
        log_output_string_list.append("        Statistical Breakdown of Reference List Flagging Differences")
        log_output_string_list.append("  Flagging differences by number         Flagging differences by percentage")
        log_output_string_list.append(
            "--------------------------------------------------------------------------------")
        log_output_string_list.append("%5s%9s%12s%10s%5s%5s %9s %12s %10s" % ("FLAG", "# TOTAL", "# UNCHANGED", "# ON->OFF", "  |  ", "FLAG", "% TOTAL", "% UNCHANGED", "% ON->OFF"))
        for ctr in range(0, len(bit_list)):
            log_output_string_list.append("%5d%9d%12d%10d%5s%5d  %8.4f %12.4f %10.4f" % (bit_list[ctr],
                                                                                         refFlagBreakdown[ctr],
                                                                                         unchangedFlagBreakdown[ctr],
                                                                                         on_off_FlagFlips[ctr],
                                                                                         "  |  ",
                                                                                         bit_list[ctr],
                                                                                         (float(refFlagBreakdown[ctr]) / n) * 100.0,
                                                                                         (float(unchangedFlagBreakdown[ctr]) / n) * 100.0,
                                                                                         (float(on_off_FlagFlips[ctr]) / n) * 100.0))
        log_output_string_list.append("%5s%9d%12d%10d%5s%5s  %8.4f %12.4f %10.4f" % ("TOTAL",
                                                                                     np.sum(refFlagBreakdown),
                                                                                     np.sum(unchangedFlagBreakdown),
                                                                                     np.sum(on_off_FlagFlips),
                                                                                     "  |  ",
                                                                                     "TOTAL",
                                                                                     (float(np.sum(refFlagBreakdown)) / n) * 100.0,
                                                                                     (float(np.sum(unchangedFlagBreakdown)) / n) * 100.0,
                                                                                     (float(np.sum(on_off_FlagFlips)) / n) * 100.0))
        log_output_string_list.append("\n")
        n = np.sum(compFlagBreakdown, dtype=float)
        log_output_string_list.append("        Statistical Breakdown of Comparison List Flagging Differences")
        log_output_string_list.append("  Flagging differences by number         Flagging differences by percentage")
        log_output_string_list.append("--------------------------------------------------------------------------------")
        log_output_string_list.append("%5s%9s%12s%10s%5s%5s %9s %12s %10s" % ("FLAG", "# TOTAL", "# UNCHANGED", "# OFF->ON", "  |  ", "FLAG", "% TOTAL", "% UNCHANGED", "% OFF->ON"))
        for ctr in range(0, len(bit_list)):
            log_output_string_list.append("%5d%9d%12d%10d%5s%5d  %8.4f %12.4f %10.4f" % (bit_list[ctr], compFlagBreakdown[ctr], unchangedFlagBreakdown[ctr], off_on_FlagFlips[ctr], "  |  ", bit_list[ctr], (float(compFlagBreakdown[ctr]) / n) * 100.0, (float(unchangedFlagBreakdown[ctr]) / n) * 100.0, (float(off_on_FlagFlips[ctr]) / n) * 100.0))
        log_output_string_list.append("%5s%9d%12d%10d%5s%5s  %8.4f %12.4f %10.4f" % ("TOTAL", np.sum(compFlagBreakdown), np.sum(unchangedFlagBreakdown), np.sum(off_on_FlagFlips), "  |  ", "TOTAL", (float(np.sum(compFlagBreakdown)) / n) * 100.0, (float(np.sum(unchangedFlagBreakdown)) / n) * 100.0, (float(np.sum(off_on_FlagFlips)) / n) * 100.0))
        log_output_string_list.append("\n")
        log_output_string_list.append(
            "Total flag bit differences......... {}\n".format(np.sum([off_on_FlagFlips, on_off_FlagFlips])))
    log_output_string_list.append(
        "{}{}".format(" " * 8, "Overall Percentage of Matched Sources with Flagging Differences"))
    log_output_string_list.append(
        "Total number of matched sources with unchanged flag values...... {}{}".format(" " * uch_padding, n_unchanged))
    log_output_string_list.append(
        "Total number of matched sources with flag value differences..... {}{}".format(" " * ch_padding, n_changed))
    log_output_string_list.append("Total number of matched sources................................. {}".format(n_total))
    log_output_string_list.append(
        "Percentage of all matched sources with flag value differences... {:7.4f}%".format(pct_changed))
    log_output_string_list.append("Regression test status............. {}".format(regTestStatus))

    for log_line in log_output_string_list:
        if log_line == "\n":
            log.info("")
        else:
            log.info(log_line)

    if plotGen != "none":
        idx = np.arange(9)
        x_title_list = []
        for bit in bit_list:
            x_title_list.append(str(bit))

        # plot flag breakdown by bit for all matched sources in the reference and comparison sourcelists
        width = 0.35
        fig = plt.figure(figsize=(11, 8.5))
        ax1 = fig.add_subplot(111)
        plt.bar(idx - width / 2, refFlagBreakdown, width, label='Reference')
        plt.bar(idx + width / 2, compFlagBreakdown, width, label='Comparison')
        plt.legend()
        plt.xticks(idx, x_title_list)
        plt.xlabel("Flag Bit Value")
        plt.ylabel("Number of matched sources")
        fullPlotTitle = "Flag Breakdown by Bit"
        plt.title(fullPlotTitle)

        if plotGen == "screen":
            plt.show()
            plotFileName = ""
            stat_file_name = ""
        if plotGen == "file":
            # Put timestamp and plotfile_prefix text string in lower left corner below plot
            timestamp = "Generated {}".format(datetime.now().strftime("%m/%d/%Y %H:%M:%S"))
            plt.text(0.0, -0.081, timestamp, horizontalalignment='left', verticalalignment='center', fontsize=5,
                     transform=ax1.transAxes)
            plt.text(0.0, -0.105, plotfile_prefix, horizontalalignment='left', verticalalignment='center', fontsize=5,
                     transform=ax1.transAxes)
            plotFileName = "{}_{}.pdf".format(plotfile_prefix, fullPlotTitle.replace(" ", "_"))
            if plotFileName.startswith("_"):
                plotFileName = plotFileName[1:]
            plt.savefig(plotFileName)
            plt.close()
            # generate second pdf page with statistics
            stat_file_name = plotFileName.replace(".pdf", "_stats.pdf")
            fig = plt.figure(figsize=(11, 8.5))
            fig.text(0.5, 0.92, fullPlotTitle, transform=fig.transFigure, size=12, ha="center")
            stat_text_blob = ""
            for log_line in log_output_string_list:
                if log_line != "\n":
                    stat_text_blob = "{}{}\n".format(stat_text_blob, log_line)
                else:
                    stat_text_blob += "\n"
            stat_text_blob += "\n" + timestamp + "\n"
            stat_text_blob += plotfile_prefix
            fig.text(0.5, 0.5, stat_text_blob, transform=fig.transFigure, size=10, ha="center", va="center",
                     multialignment="left", family="monospace")
            fig.savefig(stat_file_name)
            plt.close()

        pdf_file_list += [plotFileName, stat_file_name]

        # plot flag changes broken down by bit
        fig = plt.figure(figsize=(11, 8.5))
        ax2 = fig.add_subplot(111)
        p_unchanged = plt.bar(idx, unchangedFlagBreakdown)
        p_offOn = plt.bar(idx, off_on_FlagFlips, bottom=unchangedFlagBreakdown)
        p_onOff = plt.bar(idx, on_off_FlagFlips, bottom=off_on_FlagFlips + unchangedFlagBreakdown)
        plt.xticks(idx, x_title_list)
        plt.xlabel("Flag Bit Value")
        plt.ylabel("Number of matched sources")
        fullPlotTitle = "Flagging Differences by Bit"
        plt.title(fullPlotTitle)
        plt.legend((p_onOff[0], p_offOn, p_unchanged[0]), ("On -> Off", "Off -> On", "Unchanged"))
        if plotGen == "screen":
            plt.show()
            plotFileName = ""
        if plotGen == "file":
            # Put timestamp and plotfile_prefix text string in lower left corner below plot
            timestamp = "Generated {}".format(datetime.now().strftime("%m/%d/%Y %H:%M:%S"))
            plt.text(0.0, -0.081, timestamp, horizontalalignment='left', verticalalignment='center', fontsize=5,
                     transform=ax2.transAxes)
            plt.text(0.0, -0.105, plotfile_prefix, horizontalalignment='left', verticalalignment='center', fontsize=5,
                     transform=ax2.transAxes)
            plotFileName = "{}_{}.pdf".format(plotfile_prefix, fullPlotTitle.replace(" ", "_"))
            if plotFileName.startswith("_"):
                plotFileName = plotFileName[1:]
            plt.savefig(plotFileName)
            plt.close()
        pdf_file_list.append(plotFileName)

    regTestStatus = "%s %11.7f" % (regTestStatus, pct_changed)
    return (regTestStatus, pdf_file_list)


# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
def computeLinearStats(matchedRA, max_diff, plotGen, plot_title, plotfile_prefix, verbose, plate_scale=None):
    """Compute stats on the quantities with differences that can be computed with simple subtraction 
    (X, Y, RA, Dec, Flux, and Magnitude).

    Parameters
    ----------  
    matchedRA : numpy.ndarray
        A 2 x len(refLines) sized numpy array. Column 1: matched reference values. Column 2: The corresponding matched
    comparision values

    max_diff : float
        Maximum allowable value for comparison test to be declared a success

    plotGen : str
        Generate plots?

    plot_title : str
        text string that will be used in plot title.

    plotfile_prefix : str
        text string that will prepend the plot files generated if plots are written to files

    verbose : bool
        display verbose output?

    plate_scale : float
        plate scale, in arcseconds/pixel


    Returns
    -------
    regTestStatus : str
        overall test result and statistics
    """
    log.info(">>>>>> Comparision - reference sourcelist {} absolute differences <<<<<<".format(plot_title))
    # remove any "inf" or "nan" values in matchedRA.
    nanIdx = np.where(np.isnan(matchedRA) == True)[1]
    if len(nanIdx) > 0:
        if verbose:
            log.info("{} np.nan values will be removed from input list. New input list contains {} values.".format(len(nanIdx), len(matchedRA[0, :]) - len(nanIdx)))
        matchedRA = np.delete(matchedRA, nanIdx, axis=1)

    infIdx = np.where(np.isinf(matchedRA) == True)[1]
    if len(infIdx) > 0:
        if verbose:
            log.info("{} np.inf values will be removed from input list. New input list contains {} values.".format(len(infIdx), len(matchedRA[0, :]) - len(infIdx)))
        matchedRA = np.delete(matchedRA, infIdx, axis=1)

    # 'sigma' and 'iters' input values used for various np.sigma_clipped_stats() runs
    sigVal = 3
    intersVal = 3

    diffRA = matchedRA[1, :] - matchedRA[0, :]
    if plate_scale:  # Convert X and Y differences from pixels to arcseconds
        diffRA = diffRA * plate_scale
    clippedStats = sigma_clipped_stats(diffRA, sigma=sigVal, maxiters=intersVal)
    sigma_percentages = []
    for sig_val in [1.0, 2.0, 3.0]:
        sigma_percentages.append((float(np.shape(np.where((diffRA >= (clippedStats[0] - sig_val * clippedStats[2])) & (diffRA <= (clippedStats[0] + sig_val * clippedStats[2]))))[1])/float(np.shape(diffRA)[0])) * 100.0)

    out_stats = "%11.7f %11.7f  %11.7f  %11.7f  %11.7f " % (clippedStats[0], clippedStats[1], clippedStats[2], sigma_percentages[2], 100.0-sigma_percentages[2])

    if ((sigma_percentages[2] >=95.0) and (abs(clippedStats[0])) <= max_diff):  # success condition: Greater than or equal to 95% of all difference values within 3-sigma of sigma-clipped mean
        regTestStatus = "OK      "
    else:
        regTestStatus = "FAILURE "
    if ((verbose == True) or (regTestStatus == "FAILURE ")):
        log_output_string_list = []
        log_output_string_list.append("            Non-Clipped Statistics")
        log_output_string_list.append("Non-clipped minimum....................... {}".format(np.min(diffRA)))
        log_output_string_list.append("Non-clipped maximum....................... {}".format(np.max(diffRA)))
        log_output_string_list.append("Non-clipped mean.......................... {}".format(np.mean(diffRA)))
        log_output_string_list.append("Non-clipped median........................ {}".format(np.median(diffRA)))
        log_output_string_list.append("Non-clipped standard deviation............ {}".format(np.std(diffRA)))
        log_output_string_list.append(
            "Non-clipped mean in units of SD........... {}\n".format(np.divide(np.mean(diffRA), np.std(diffRA))))
        log_output_string_list.append("       Sigma-clipped Statistics; Sigma = {}, # steps = {}".format(sigVal, intersVal))
        log_output_string_list.append("Sigma-clipped mean........................ {}".format(clippedStats[0]))
        log_output_string_list.append("Sigma-clipped median...................... {}".format(clippedStats[1]))
        log_output_string_list.append("Sigma-clipped standard deviation.......... {}".format(clippedStats[2]))
        for sig_val, pct_val in zip([1.0, 2.0, 3.0], sigma_percentages):
            log_output_string_list.append("% all diff values within {} sigma of mean...{}".format(sig_val,pct_val))

    log_output_string_list.append("Regression test status.................... {}".format(regTestStatus))
    log_output_string_list.append("\n")


    if plotGen == "none":
        pdf_files = []
        log.info("\n")
        for log_line in log_output_string_list:
            if log_line == "\n":
                log.info("")
            else:
                log.info(log_line)

    else:

        xAxisString = "{}".format(plot_title.split(" ")[0])
        plotCutoff = (10.0 * np.abs(clippedStats[2])) + np.abs(clippedStats[0])
        if plotCutoff != 0.0:
            origSize = len(diffRA)
            log_output_string_list.append("Plot cutoff: {}".format(plotCutoff))
            goodIdx = np.where(np.abs(diffRA) <= plotCutoff)
            diffRA = diffRA[goodIdx]
            log_output_string_list.append("%d values (%7.4f percent) clipped from plot." % (origSize - len(diffRA), (float(origSize - len(diffRA)) / float(origSize)) * 100.0))
        fig = plt.figure(figsize=(11, 8.5))
        ax1 = fig.add_subplot(111)
        fullPlotTitle = "Comparision - reference sourcelist %s absolute differences" % (plot_title)
        plt.title(fullPlotTitle)
        bins = "auto"
        ax1.hist(diffRA, bins=bins)
        ax1.axvline(x=clippedStats[0], color='k', linestyle='--')
        ax1.axvline(x=clippedStats[1], color='r', linestyle='--')
        ax1.axvline(x=clippedStats[0] + 3.0*clippedStats[2], color='k', linestyle=':')
        ax1.axvline(x=clippedStats[0] - 3.0*clippedStats[2], color='k', linestyle=':')

        ax1.set_xlabel("$\Delta %s$" % (xAxisString))
        ax1.set_ylabel("Number of matched sources")

        ax2 = ax1.twinx()
        ax2.hist(diffRA, bins=bins, cumulative=-1, density=True, histtype='step', color='r')
        ax2.set_ylabel("Fraction of all matched sources", color='r')
        for tl in ax2.get_yticklabels():
            tl.set_color('r')

        log.info("\n")
        for log_line in log_output_string_list:
            if log_line == "\n":
                log.info("")
            else:
                log.info(log_line)

        if plotGen == "screen":
            plt.show()
            plotFileName = ""
            stat_file_name = ""
        if plotGen == "file":
            # Put timestamp and plotfile_prefix text string in lower left corner below plot
            timestamp = "Generated {}".format(datetime.now().strftime("%m/%d/%Y %H:%M:%S"))
            plt.text(0.0, -0.081, timestamp, horizontalalignment='left', verticalalignment='center', fontsize=5,
                     transform=ax1.transAxes)
            plt.text(0.0, -0.105, plotfile_prefix, horizontalalignment='left', verticalalignment='center', fontsize=5,
                     transform=ax1.transAxes)
            plotFileName = "{}_{}.pdf".format(plotfile_prefix, plot_title.replace(" ", "_"))
            if plotFileName.startswith("_"):
                plotFileName = plotFileName[1:]
            plt.savefig(plotFileName)
            plt.close()
            # generate second pdf page with statistics
            stat_file_name = plotFileName.replace(".pdf", "_stats.pdf")
            fig = plt.figure(figsize=(11, 8.5))
            fig.text(0.5, 0.87, fullPlotTitle[:-1] + " statistics", transform=fig.transFigure, size=12, ha="center")
            stat_text_blob = ""
            for log_line in log_output_string_list:
                if log_line != "\n":
                    stat_text_blob += log_line + "\n"
                else:
                    stat_text_blob += "\n"
            stat_text_blob += "\n" + timestamp + "\n"
            stat_text_blob += plotfile_prefix
            fig.text(0.5, 0.5, stat_text_blob, transform=fig.transFigure, size=10, ha="center", va="center",
                     multialignment="left", family="monospace")
            fig.savefig(plotFileName.replace(".pdf", "_stats.pdf"))
            plt.close()
        pdf_files = [plotFileName, stat_file_name]

    log.info("\n")

    return (regTestStatus + out_stats, pdf_files)


# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
def deconstruct_flag(flagval):
    """Breaks down an integer flag value into individual component bit values.

    Parameters
    ----------
    flagval : int
        Flag value to deconstruct
    
    Returns
    -------
    out_idx_list : list
        a 9-element numpy array of 0s and 1s. Each element of the array represents the presence of a particular 
        bit value (element 0 = bit 0, element 1 = bit 1, ..., element 3 = bit 4 and so on...)
    """
    bit_list = [1, 2, 4, 8, 16, 32, 64, 128]
    flagval = int(flagval)
    # out_bit_list = []
    out_idx_list = np.zeros(9, dtype=int)
    if flagval == 0:
        # out_bit_list = [0]
        out_idx_list[0] = 1
    if flagval > 0:
        idx = 1
        for bit in bit_list:
            if flagval & bit > 0:
                # out_bit_list.append(bit)
                out_idx_list[idx] = 1
            if bit > flagval:
                break
            idx += 1
    return (out_idx_list)


# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
def extractMatchedLines(col2get, refData, compData, refLines, compLines):
    """Extracts only matching lines of data for a specific column of refData and compData. Returns empty list if the
    specified column is not found in both tables.

    Parameters
    ----------
    col2get : str
        Title of the column to return

    refData : astropy Table object
        reference data table

    compData : astropy Table object
        comparison data table

    refLines : numpy.ndarray
        List of matching refData line numbers

    compLines : numpy.ndarray
        List of matching compData line numbers

    Returns
    -------
    return_ra : numpy ndarray
        A 2 x len(refLines) sized numpy array. Column 1: matched reference values. Column 2: The corresponding matched
        comparision values
    """
    return_ra = []
    if col2get in list(refData.keys()) and col2get in list(compData.keys()):
        matching_refData = refData[col2get][refLines].data
        matching_compData = compData[col2get][compLines].data
        return_ra = np.stack((matching_refData, matching_compData))
    return (return_ra)


# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
def getMatchedLists(slNames, imgNames, slLengths, log_level):
    """run starmatch_hist to get the indices of matching sources that are common to both input source catalogs

    Parameters
    ----------
    slNames : list
        list of input source lists

    imgNames : list
        list of input images

    slLengths : list
        list of integer sourcelist lengths

    log_level : int
        The desired level of verboseness in the log statements displayed on the screen and written to the .log file.

    Returns
    -------
    matching_lines_ref : list
        A list of the indices of reference sourcelist sources that match comparison sourcelist sources

    matching_lines_img : list
        A corresponding list of the indices of comparison sourcelist sources that match reference sourcelist sources
    """
    source_list_dict = {}
    equal_flag = False
    if slLengths[0] == slLengths[1]:
        slLengths[0] += 1
        equal_flag = True

    for ctr, slName in enumerate(slNames, 0):
        source_list_dict[slName] = slLengths[ctr]

    try:
        fh = fits.open(imgNames[0])
        xref = fh[0].header['crpix1']
        yref = fh[0].header['crpix2']
        fh.close()
    except Exception:
        log.info(
            "WARNING: Unable to fetch values for xref and yref from fits file headers. Using xref = 0.0 and yref = 0.0.")
        xref = 0.0
        yref = 0.0
    log.info("Summary of input catalog lengths")
    for source_list in source_list_dict.keys():
        log.info("{}: {}".format(os.path.basename(source_list), source_list_dict[source_list]))
    out_dict = starmatch_hist.run(source_list_dict, log_level, xref=xref, yref=yref)
    matching_lines_ref = out_dict[slNames[0]]
    matching_lines_img = out_dict[slNames[1]]

    if equal_flag:
        slLengths[0] -= 1
    # Report number and percentage of the total number of detected ref and comp sources that were matched
    log.info("Sourcelist Matching Results")
    log.info(
        "Reference sourcelist:  {} of {} total sources matched ({} %)".format(len(matching_lines_ref), slLengths[0],
                                                                              100.0 * (float(
                                                                                  len(matching_lines_ref)) / float(
                                                                                  slLengths[0]))))
    log.info(
        "Comparison sourcelist: {} of {} total sources matched ({} %)".format(len(matching_lines_img), slLengths[1],
                                                                              100.0 * (float(
                                                                                  len(matching_lines_img)) / float(
                                                                                  slLengths[1]))))

    return (matching_lines_ref, matching_lines_img)


# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
def makeVectorPlot(x, y, plotDest, plotfile_prefix, binThresh=10000, binSize=250):
    """Generate vector plot of dx and dy values vs. reference (x,y) positions

    Parameters
    ----------
    x : numpy.ndarray
        A 2 x n sized numpy array. Column 1: matched reference X values. Column 2: The corresponding matched
        comparision X values

    y : numpy.ndarray
        A 2 x n sized numpy array. Column 1: matched reference Y values. Column 2: The corresponding matched
        comparision Y values

    plotDest : str
        plot destination; screen or file

    plotfile_prefix : str
        text string that will prepend the plot files generated if plots are written to files

    binThresh : int, optional
        Minimum size of list *x* and *y* that will trigger generation of a binned vector plot. Default value = 10000.

    binSize : int, optional
        Size of binning box in pixels. When generating a binned vector plot, mean dx and dy values are computed by
        taking the mean of all points located within the box. Default value = 250.

    Returns
    -------
    nothing
    """
    dx = x[1, :] - x[0, :]
    dy = y[1, :] - y[0, :]
    if len(dx) > binThresh:  # if the input list is larger than binThresh, a binned vector plot will be generated.
        binStatus = "%d x %d Binning" % (binSize, binSize)
        log.info("Input list length greater than threshold length value. Generating binned vector plot using %d pixel x %d pixel bins" % (binSize, binSize))
        if min(x[0, :]) < 0.0:
            xmin = min(x[0, :])
        else:
            xmin = 0.0
        if min(y[0, :]) < 0.0:
            ymin = min(y[0, :])
        else:
            ymin = 0.0

        p_x = np.empty(shape=[0])
        p_y = np.empty(shape=[0])
        p_dx = np.empty(shape=[0])
        p_dy = np.empty(shape=[0])
        color_ra = []
        for xBinCtr in range(int(round2ArbatraryBase(xmin, "down", binSize)),
                             int(round2ArbatraryBase(max(x[0, :]), "up", binSize)), binSize):
            for yBinCtr in range(int(round2ArbatraryBase(ymin, "down", binSize)),
                                 int(round2ArbatraryBase(max(y[0, :]), "up", binSize)), binSize):
                # define bin box x,y upper and lower bounds
                xBinMin = xBinCtr
                xBinMax = xBinMin + binSize
                yBinMin = yBinCtr
                yBinMax = yBinMin + binSize
                # get indicies of x and y within bounding box
                ix0 = np.where((x[0, :] >= xBinMin) & (x[0, :] < xBinMax) & (y[0, :] >= yBinMin) & (y[0, :] < yBinMax))
                if len(dx[ix0]) > 0 and len(dy[ix0]) > 0:  # ignore empty bins
                    p_x = np.append(p_x, xBinCtr + 0.5 * binSize)  # X and Y position at center of bin.
                    p_y = np.append(p_y, yBinCtr + 0.5 * binSize)
                    mean_dx = np.mean(dx[ix0])
                    p_dx = np.append(p_dx, mean_dx)  # compute mean dx, dy values
                    mean_dy = np.mean(dy[ix0])
                    p_dy = np.append(p_dy, mean_dy)
                    avg_npts = (float(len(dx[ix0])) + float(
                        len(dy[ix0]))) / 2.0  # keep an eye out for mean values computed from less than 10 samples.
                    if (
                            avg_npts < 10.0):  # if less than 10 samples were used in mean calculation, color the vector red.
                        color_ra.append('r')
                    else:
                        color_ra.append('k')
        lowSampleWarning = ""
        if "r" in color_ra:
            lowSampleWarning = "; Red Vectors were computed with less than 10 values"
    else:
        log.info("Generating unbinned vector plot")
        binStatus = "Unbinned"
        lowSampleWarning = ""
        color_ra = ["k"]
        p_x = x[0, :]
        p_y = y[0, :]
        p_dx = dx
        p_dy = dy
    plt_mean = np.mean(np.hypot(p_dx, p_dy))
    e = np.log10(5.0 * plt_mean).round()
    plt_scaleValue = 10 ** e
    fig = plt.figure(figsize=(11, 8.5))
    ax1 = fig.add_subplot(111)
    if len(dx) > binThresh:
        Q = plt.quiver(p_x, p_y, p_dx, p_dy, color=color_ra, units="xy")
    else:
        Q = plt.quiver(p_x, p_y, p_dx, p_dy)
    plt.quiverkey(Q, 0.75, 0.05, plt_scaleValue, r'%5.3f' % (plt_scaleValue), labelpos='S', coordinates='figure',
                  color="k")
    plot_title = "Comparision - reference $\Delta X$, $\Delta Y$ values vs. $(X_{ref}, Y_{ref})$ positions\n%s%s" % (binStatus, lowSampleWarning)
    plt.title(plot_title)
    plt.xlabel(r"$X_{ref}$")
    plt.ylabel(r"$Y_{ref}$")
    if plotDest == "screen":
        plt.show()
        plotFileName = ""
    if plotDest == "file":
        # Put timestamp and plotfile_prefix text string in lower left corner below plot
        timestamp = "Generated {}".format(datetime.now().strftime("%m/%d/%Y %H:%M:%S"))
        plt.text(0.0, -0.081, timestamp, horizontalalignment='left', verticalalignment='center', fontsize=5,
                 transform=ax1.transAxes)
        plt.text(0.0, -0.105, plotfile_prefix, horizontalalignment='left', verticalalignment='center', fontsize=5,
                 transform=ax1.transAxes)
        plotFileName = "{}_xy_vector_plot.pdf".format(plotfile_prefix, plot_title.replace(" ", "_"))
        if plotFileName.startswith("_"):
            plotFileName = plotFileName[1:]
        plt.savefig(plotFileName)
        plt.close()
        log.info("Vector plot saved to file {}".format(plotFileName))
        log.info("\n")
    return plotFileName


# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
def round2ArbatraryBase(value, direction, roundingBase):
    """Round value up or down to arbitrary base

    Parameters
    ----------
    value : float
        Value to be rounded.

    direction : str
        Round up, down to the nearest base (choices: "up","down","nearest")

    roundingBase : int
        rounding base. (Example: if base = 5, values will be rounded to the nearest multiple of 5 or 10.)

    Returns
    -------
    rv : int
        rounded value.
    """
    if direction.lower().startswith("u"):
        rv = value + (roundingBase - value % roundingBase)  # round up to nearest base
    elif direction.lower().startswith("d"):
        rv = value - value % roundingBase  # round down to nearest base
    else:
        rv = int(roundingBase * round(float(value) / roundingBase))  # round up or down to nearest base
    return rv


# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-

def comparesourcelists(slNames, imgNames, plotGen=None, plotfile_prefix=None, verbose=False,
                       log_level=logutil.logging.NOTSET, debugMode=False):
    """Main calling subroutine to compare sourcelists.

    Parameters
    ----------
    slNames : list
        list of input source lists

    imgNames : list
        optional list of input images that starmatch_hist will use to improve sourcelist matching

    plotfile_prefix : str, optional
        text string that will prepend the plot files generated if plots are written to files

    plotGen : bool, optional
        Generate plots and display them to the screen (True/False)?

    diffMode : str, optional
        method used to compute comp-ref difference in computeLinearStats().

    verbose : bool, optional
        display verbose output?

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the .log file.
        Default value is 'NOTSET'.

    debugMode : bool, optional
        executes subroutine check_match_quality(), which the writes the matched sources (x, y) coordinates of the
        comparision and reference source lists to ds9 region files for follow-up human visual inspection.

    Returns
    -------
    overallStatus : str
        "OK" if all tests were passed, or "FAILURE" if inconsistencies were found.
    """
    log.setLevel(log_level)
    if not plotfile_prefix:
        plotfile_prefix = ""
    regressionTestResults = {}
    colTitles = []
    pdf_file_list = []

    # 0: define dictionary of max allowable mean sigma-clipped difference values
    max_diff_dict = {"X Position": 0.1,
                     "Y Position": 0.1,
                     "RA Position" : 0.1,
                     "DEC Position": 0.1,
                     "Flux (Inner Aperture)": 999.0,
                     "Flux (Outer Aperture)": 999.0,
                     "Magnitude (Inner Aperture)": 999.0,
                     "Magnitude (Inner Aperture) Error": 999.0,
                     "Magnitude (Outer Aperture)": 999.0,
                     "Magnitude (Outer Aperture) Error": 999.0,
                     "MSKY value": 999.0,
                     "STDEV value": 999.0,
                     "CI": 999.0,
                     "Flag Value": 5.0}
    # 1: Read in sourcelists files into astropy table or 2-d array so that individual columns from each sourcelist can be easily accessed later in the code.
    refData, compData = slFiles2dataTables(slNames)
    log.info("Valid reference data columns:   {}".format(list(refData.keys())))
    log.info("Valid comparision data columns: {}".format(list(compData.keys())))
    log.info("\n")
    log.info("Data columns to be compared:")
    for listItem in sorted(list(set(refData.keys()).intersection(set(compData.keys())))):
        log.info(listItem)
    log.info("\n")
    # 2: Run starmatch_hist to get list of matched sources common to both input sourcelists
    slLengths = [len(refData['X']), len(compData['X'])]
    matching_lines_ref, matching_lines_img = getMatchedLists(slNames, imgNames, slLengths, log_level)
    if len(matching_lines_ref) == 0 or len(matching_lines_img) == 0:
        log.critical("*** Comparisons cannot be computed. No matching sources were found. ***")
        return ("ERROR")

    # 3: Compute and display statistics on X position differences for matched sources
    # Get platescale
    plate_scale = wcsutil.HSTWCS(imgNames[0], ext=('sci', 1)).pscale

    matched_values = extractMatchedLines("X", refData, compData, matching_lines_ref, matching_lines_img)
    if len(matched_values) > 0:
        formalTitle = "X Position"
        rt_status, pdf_files = computeLinearStats(matched_values, max_diff_dict[formalTitle], plotGen, formalTitle, plotfile_prefix,
                                                  verbose, plate_scale=plate_scale)
        if plotGen == "file":
            pdf_file_list = pdf_files
        regressionTestResults[formalTitle] = rt_status
        colTitles.append(formalTitle)
        matchedXValues = matched_values.copy()
    # 4: Compute and display statistics on Y position differences for matched sources
    matched_values = extractMatchedLines("Y", refData, compData, matching_lines_ref, matching_lines_img)
    if len(matched_values) > 0:
        formalTitle = "Y Position"
        rt_status, pdf_files = computeLinearStats(matched_values, max_diff_dict[formalTitle], plotGen, formalTitle, plotfile_prefix,
                                                  verbose, plate_scale=plate_scale)
        if plotGen == "file":
            pdf_file_list += pdf_files
        regressionTestResults[formalTitle] = rt_status
        colTitles.append(formalTitle)
        matchedYValues = matched_values.copy()
        if plotGen != "none":
            pdf_file_name = makeVectorPlot(matchedXValues, matchedYValues, plotGen, plotfile_prefix)
            if plotGen == "file":
                pdf_file_list.append(pdf_file_name)
    if debugMode:
        check_match_quality(matchedXValues, matchedYValues)
    # 5: Compute and display statistics on RA position differences for matched sources
    matched_values = extractMatchedLines("RA", refData, compData, matching_lines_ref, matching_lines_img)
    if len(matched_values) > 0:
        formalTitle = "RA Position"
        rt_status, pdf_files = computeLinearStats(matched_values, max_diff_dict[formalTitle], plotGen, formalTitle, plotfile_prefix,
                                                  verbose)
        if plotGen == "file":
            pdf_file_list += pdf_files
        regressionTestResults[formalTitle] = rt_status
        colTitles.append(formalTitle)

    # 6: Compute and display statistics on DEC position differences for matched sources
    matched_values = extractMatchedLines("DEC", refData, compData, matching_lines_ref, matching_lines_img)
    if len(matched_values) > 0:
        formalTitle = "DEC Position"
        rt_status, pdf_files = computeLinearStats(matched_values, max_diff_dict[formalTitle], plotGen, formalTitle, plotfile_prefix,
                                                  verbose)
        if plotGen == "file":
            pdf_file_list += pdf_files
        regressionTestResults[formalTitle] = rt_status
        colTitles.append(formalTitle)

    # 7: Compute and display statistics on flux differences for matched sources
    matched_values = extractMatchedLines("FLUX1", refData, compData, matching_lines_ref, matching_lines_img)
    if len(matched_values) > 0:
        formalTitle = "Flux (Inner Aperture)"
        rt_status, pdf_files = computeLinearStats(matched_values, max_diff_dict[formalTitle], plotGen, formalTitle,
                                                  plotfile_prefix, verbose)
        if plotGen == "file":
            pdf_file_list += pdf_files
        regressionTestResults[formalTitle] = rt_status
        colTitles.append(formalTitle)

    matched_values = extractMatchedLines("FLUX2", refData, compData, matching_lines_ref, matching_lines_img)
    if len(matched_values) > 0:
        formalTitle = "Flux (Outer Aperture)"
        rt_status, pdf_files = computeLinearStats(matched_values, max_diff_dict[formalTitle], plotGen, formalTitle,
                                                  plotfile_prefix, verbose)
        if plotGen == "file":
            pdf_file_list += pdf_files
        regressionTestResults[formalTitle] = rt_status
        colTitles.append(formalTitle)

    # 8: Compute and display statistics on magnitude differences for matched sources
    matched_values = extractMatchedLines("MAGNITUDE1", refData, compData, matching_lines_ref, matching_lines_img)
    if len(matched_values) > 0:
        formalTitle = "Magnitude (Inner Aperture)"
        rt_status, pdf_files = computeLinearStats(matched_values, max_diff_dict[formalTitle], plotGen, formalTitle,
                                                  plotfile_prefix, verbose)
        if plotGen == "file":
            pdf_file_list += pdf_files
        regressionTestResults[formalTitle] = rt_status
        colTitles.append(formalTitle)

    matched_values = extractMatchedLines("MERR1", refData, compData, matching_lines_ref, matching_lines_img)
    if len(matched_values) > 0:
        formalTitle = "Magnitude (Inner Aperture) Error"
        rt_status, pdf_files = computeLinearStats(matched_values, max_diff_dict[formalTitle], plotGen, formalTitle, plotfile_prefix,
                                                  verbose)
        if plotGen == "file":
            pdf_file_list += pdf_files
        regressionTestResults[formalTitle] = rt_status
        colTitles.append(formalTitle)

    matched_values = extractMatchedLines("MAGNITUDE2", refData, compData, matching_lines_ref, matching_lines_img)
    if len(matched_values) > 0:
        formalTitle = "Magnitude (Outer Aperture)"
        rt_status, pdf_files = computeLinearStats(matched_values, max_diff_dict[formalTitle], plotGen, formalTitle,
                                                  plotfile_prefix, verbose)
        if plotGen == "file":
            pdf_file_list += pdf_files
        regressionTestResults[formalTitle] = rt_status
        colTitles.append(formalTitle)

    matched_values = extractMatchedLines("MERR2", refData, compData, matching_lines_ref, matching_lines_img)
    if len(matched_values) > 0:
        formalTitle = "Magnitude (Outer Aperture) Error"
        rt_status, pdf_files = computeLinearStats(matched_values, max_diff_dict[formalTitle], plotGen, formalTitle, plotfile_prefix,
                                                  verbose)
        if plotGen == "file":
            pdf_file_list += pdf_files
        regressionTestResults[formalTitle] = rt_status
        colTitles.append(formalTitle)

    # 9: Compute and display statistics on differences in background sky values
    matched_values = extractMatchedLines("MSKY", refData, compData, matching_lines_ref, matching_lines_img)
    if len(matched_values) > 0:
        formalTitle = "MSKY value"
        rt_status, pdf_files = computeLinearStats(matched_values, max_diff_dict[formalTitle], plotGen, formalTitle, plotfile_prefix,
                                                  verbose)
        if plotGen == "file":
            pdf_file_list += pdf_files
        regressionTestResults[formalTitle] = rt_status
        colTitles.append(formalTitle)

    matched_values = extractMatchedLines("STDEV", refData, compData, matching_lines_ref, matching_lines_img)
    if len(matched_values) > 0:
        formalTitle = "STDEV value"
        rt_status, pdf_files = computeLinearStats(matched_values, max_diff_dict[formalTitle], plotGen, formalTitle, plotfile_prefix,
                                                  verbose)
        if plotGen == "file":
            pdf_file_list += pdf_files
        regressionTestResults[formalTitle] = rt_status
        colTitles.append(formalTitle)

    # 10: Compute and display statistics on differences in concentration index  for matched sources
    matched_values = extractMatchedLines("CI", refData, compData, matching_lines_ref, matching_lines_img)
    if len(matched_values) > 0:
        formalTitle = "CI"
        rt_status, pdf_files = computeLinearStats(matched_values, max_diff_dict[formalTitle], plotGen, formalTitle, plotfile_prefix,
                                                  verbose)
        if plotGen == "file":
            pdf_file_list += pdf_files
        regressionTestResults[formalTitle] = rt_status
        colTitles.append(formalTitle)

    # 11: Compute and display statistics on differences in flag populations for matched sources
    matched_values = extractMatchedLines("FLAGS", refData, compData, matching_lines_ref, matching_lines_img)
    if len(matched_values) > 0:
        formalTitle = "Source Flagging"
        rt_status, flag_pdf_list = computeFlagStats(matched_values, max_diff_dict[formalTitle], plotGen, formalTitle, plotfile_prefix,
                                                    verbose)
        if plotGen == "file":
            pdf_file_list += flag_pdf_list
        regressionTestResults[formalTitle] = rt_status
        colTitles.append(formalTitle)

    log.info("\n")
    log_output_string_list = []

    lenList = []
    for item in colTitles:
        lenList.append(len(item))
    totalPaddedSize = max(lenList) + 3
    log_output_string_list.append("{}{}".format(" " * 35, "REGRESSION TESTING SUMMARY"))
    log_output_string_list.append("-" * (70 + totalPaddedSize))
    log_output_string_list.append("{}{}".format(" " * (totalPaddedSize + 46), "% within     % beyond"))
    log_output_string_list.append(
        "COLUMN{}STATUS   MEAN        MEDIAN       STD DEV     3\u03C3 of mean   3\u03C3 of mean".format(
            " " * (totalPaddedSize - 6)))
    overallStatus = "OK"
    for colTitle in colTitles:
        log_output_string_list.append(
            "%s%s%s" % (colTitle, "." * (totalPaddedSize - len(colTitle)), regressionTestResults[colTitle]))
        if not regressionTestResults[colTitle].startswith("OK"):
            overallStatus = "FAILURE"
    log_output_string_list.append("-" * (70 + totalPaddedSize))
    log_output_string_list.append("OVERALL TEST STATUS{}{}".format("." * (totalPaddedSize - 19), overallStatus))
    for log_line in log_output_string_list:
        log.info(log_line)

    if plotGen == "file":

        # generate final overall summary pdf page
        stat_summary_file_name = "stats_summary.pdf"
        if plotfile_prefix is not None:
            stat_summary_file_name = "{}_{}".format(plotfile_prefix, stat_summary_file_name)
        fig = plt.figure(figsize=(11, 8.5))
        # fig.text(0.5, 0.87, fullPlotTitle[:-1] + " statistics", transform=fig.transFigure, size=12, ha="center")
        stat_text_blob = ""
        for log_line in log_output_string_list:
            if log_line != "\n":
                stat_text_blob += log_line + "\n"
            else:
                stat_text_blob += "\n"

        stat_text_blob += "\n\nGenerated {}\n".format(datetime.now().strftime("%m/%d/%Y %H:%M:%S"))
        stat_text_blob += plotfile_prefix
        fig.text(0.5, 0.5, stat_text_blob, transform=fig.transFigure, size=10, ha="center", va="center",
                 multialignment="left", family="monospace")
        fig.savefig(stat_summary_file_name)
        plt.close()
        pdf_file_list.append(stat_summary_file_name)
        # combine all individual plot files into a single multiple-page pdf file
        final_plot_filename = "comparision_plots.pdf"
        if plotfile_prefix:
            final_plot_filename = "{}_{}".format(plotfile_prefix, final_plot_filename)

        pdf_merger(final_plot_filename, pdf_file_list)
        log.info("Sourcelist comparison plots saved to file {}.".format(final_plot_filename))
    return (overallStatus)


# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
def slFiles2dataTables(slNames):
    """Reads in data from sourcelists, returns some or all of the following data columns in a pair of properly
    formatted astropy.table objects:

    * X Position
    * Y Position
    * Right Ascension
    * Declination
    * Flux (Inner Aperture)
    * Flux (Outer Aperture)
    * Magnitude (Inner Aperture)
    * Magnitude (Outer Aperture)
    * Flag Value

    .. note::
        'X position' and 'Y position' data columns will always be returned. However, not every sourcelist or catalog
        file is guaranteed to have any of the seven remaining data columns in the above list.

    Parameters
    ----------
    slNames : list
        A list containing the reference sourcelist filename and the comparison sourcelist filename, in that order.

    Returns
    -------
    refData : astropy.table object
        data from the reference sourcelist

    compData : astropy.table object
        data from the comparision sourcelist
    """
    if slNames[0].endswith(".ecsv"):
        refData_in = Table.read(slNames[0], format='ascii.ecsv')
    else:
        try:
            refData_in = Table.read(slNames[0], format='ascii.daophot')
        except Exception:
            refData_in = Table.read(slNames[0], format='ascii')
    if slNames[1].endswith(".ecsv"):
        compData_in = Table.read(slNames[1], format='ascii.ecsv')
    else:
        try:
            compData_in = Table.read(slNames[1], format='ascii.daophot')
        except Exception:
            compData_in = Table.read(slNames[1], format='ascii')
    titleSwapDict_dao1 = {"X": "X-Center", "Y": "Y-Center", "RA": "RA", "DEC": "DEC", "FLUX1": "n/a",
                          "FLUX2": "Flux(0.15)", "MAGNITUDE1": "MagAp(0.05)", "MAGNITUDE2": "MagAp(0.15)",
                          "MERR1": "MagErr(0.05)", "MERR2": "MagErr(0.15)", "MSKY": "MSky(0.15)",
                          "STDEV": "Stdev(0.15)", "FLAGS": "Flags", "ID": "ID", "CI": "CI"}
    titleSwapDict_dao2 = {"X": "X-Center", "Y": "Y-Center", "RA": "RA", "DEC": "DEC", "FLUX1": "n/a",
                          "FLUX2": "Flux(0.45)", "MAGNITUDE1": "MagAp(0.15)", "MAGNITUDE2": "MagAp(0.45)",
                          "MERR1": "MagErr(0.15)", "MERR2": "MagErr(0.45)", "MSKY": "MSky(0.45)",
                          "STDEV": "Stdev(0.45)", "FLAGS": "Flags", "ID": "ID", "CI": "CI"}
    titleSwapDict_dao3 = {"X": "X-Center", "Y": "Y-Center", "RA": "RA", "DEC": "DEC", "FLUX1": "n/a",
                          "FLUX2": "Flux(0.125)", "MAGNITUDE1": "MagAp(0.03)", "MAGNITUDE2": "MagAp(0.125)",
                          "MERR1": "MagErr(0.03)", "MERR2": "MagErr(0.125)", "MSKY": "MSky(0.125)",
                          "STDEV": "Stdev(0.125)", "FLAGS": "Flags", "ID": "ID", "CI": "CI"}
    titleSwapDict_point = {"X": "X-Center", "Y": "Y-Center", "RA": "RA", "DEC": "DEC", "FLUX1": "n/a",
                           "FLUX2": "FluxAp2", "MAGNITUDE1": "MagAp1", "MAGNITUDE2": "MagAp2", "MERR1": "MagErrAp1",
                           "MERR2": "MagErrAp2", "MSKY": "MSkyAp2", "STDEV": "StdevAp2", "FLAGS": "Flags", "ID": "ID",
                           "CI": "CI"}
    titleSwapDict_segment = {"X": "X-Centroid", "Y": "Y-Centroid", "RA": "RA", "DEC": "DEC", "FLUX1": "n/a",
                             "FLUX2": "FluxAp2", "MAGNITUDE1": "MagAp1", "MAGNITUDE2": "MagAp2", "MERR1": "MagErrAp1",
                             "MERR2": "MagErrAp2", "MSKY": "n/a", "STDEV": "n/a", "FLAGS": "Flags", "ID": "ID",
                             "CI": "CI"}
    titleSwapDict_daoTemp = {"X": "XCENTER", "Y": "YCENTER", "RA": "n/a", "DEC": "n/a", "FLUX1": "FLUX1",
                             "FLUX2": "FLUX2", "MAGNITUDE1": "MAG1", "MAGNITUDE2": "MAG2", "FLAGS": "n/a", "ID": "ID",
                             "MERR1": "MERR1", "MERR2": "MERR2", "MSKY": "MSKY", "STDEV": "STDEV"}
    titleSwapDict_sourceX = {"X": "X_IMAGE", "Y": "Y_IMAGE", "RA": "RA", "DEC": "DEC", "FLUX1": "FLUX_APER1",
                             "FLUX2": "FLUX_APER2", "MAGNITUDE1": "MAG_APER1", "MAGNITUDE2": "MAG_APER2",
                             "FLAGS": "FLAGS", "ID": "NUMBER"}
    titleSwapDict_cooNew = {"X": "col1", "Y": "col2", "RA": "n/a", "DEC": "n/a", "FLUX1": "n/a", "FLUX2": "n/a",
                            "MAGNITUDE1": "n/a", "MAGNITUDE2": "n/a", "FLAGS": "n/a", "ID": "col7"}
    # titleSwapDict_cooOld = {"X": "XCENTER", "Y": "YCENTER", "RA": "n/a", "DEC": "n/a", "FLUX1": "n/a", "FLUX2": "n/a",
    #                         "MAGNITUDE1": "n/a", "MAGNITUDE2": "n/a", "FLAGS": "n/a", "ID": "ID"}

    titleSwapDict_cooOld2 = {"X": "XCENTER", "Y": "YCENTER", "RA": "RA", "DEC": "DEC", "FLUX1": "FLUX_0.05",
                             "FLUX2": "FLUX_0.15", "MAGNITUDE1": "MAG_0.05", "MAGNITUDE2": "MAG_0.15", "FLAGS": "n/a",
                             "ID": "ID", "MERR1": "MERR_0.05", "MERR2": "MERR_0.15", "MSKY": "MSKY", "STDEV": "STDEV"}

    titleSwapDict_daorep = {"X": "X", "Y": "Y", "RA": "n/a", "DEC": "n/a", "FLUX1": "flux_0", "FLUX2": "flux_1",
                            "MAGNITUDE1": "mag_0", "MAGNITUDE2": "mag_1", "FLAGS": "n/a", "ID": "n/a"}
    ctr = 1
    for dataTable in [refData_in, compData_in]:
        if "X-Center" in list(dataTable.keys()):
            if (("MagAp(0.05)" in list(dataTable.keys())) and (
                    "MagAp(0.15)" in list(dataTable.keys()))):  # ACS/WFC, WFC3/UVIS
                log.info("titleSwapDict_dao1")
                titleSwapDict = titleSwapDict_dao1
            elif (("MagAp(0.15)" in list(dataTable.keys())) and ("MagAp(0.45)" in list(dataTable.keys()))):  # WFC3/IR
                log.info("titleSwapDict_dao2")
                titleSwapDict = titleSwapDict_dao2
            elif (("MagAp(0.03)" in list(dataTable.keys())) and ("MagAp(0.125)" in list(dataTable.keys()))):  # ACS/HRC
                log.info("titleSwapDict_dao3")
                titleSwapDict = titleSwapDict_dao3
            elif "MagAp1" in list(dataTable.keys()):
                log.info("titleSwapDict_point")
                titleSwapDict = titleSwapDict_point
            else:
                sys.exit("ERROR: Unrecognized format. Exiting...")
        elif ("XCENTER" in list(dataTable.keys()) and "FLUX1" in list(dataTable.keys())):
            log.info("titleSwapDict_daoTemp")
            titleSwapDict = titleSwapDict_daoTemp
        elif "X_IMAGE" in list(dataTable.keys()):
            log.info("titleSwapDict_sourceX")
            titleSwapDict = titleSwapDict_sourceX
        elif "col1" in list(dataTable.keys()):
            log.info("titleSwapDict_cooNew")
            titleSwapDict = titleSwapDict_cooNew
        elif "XCENTER" in list(dataTable.keys()):
            log.info("titleSwapDict_cooOld2")
            titleSwapDict = titleSwapDict_cooOld2
        elif "X" in list(dataTable.keys()):
            log.info("titleSwapDict_daorep")
            titleSwapDict = titleSwapDict_daorep
        elif "X-Centroid" in list(dataTable.keys()) and "MagAp1" in list(dataTable.keys()):
            log.info("titleSwapDict_segment")
            titleSwapDict = titleSwapDict_segment
        else:
            sys.exit("ERROR: Unrecognized format. Exiting...")
        outTable = Table()
        for swapKey in list(titleSwapDict.keys()):
            if titleSwapDict[swapKey] != "n/a":
                try:
                    col2add = Table.Column(name=swapKey, data=dataTable[titleSwapDict[swapKey]])
                except TypeError:
                    col2add = Table.MaskedColumn(name=swapKey, data=dataTable[titleSwapDict[swapKey]])
                outTable.add_column(col2add)
        if ctr == 1:
            refData = outTable
        if ctr == 2:
            compData = outTable
        ctr += 1

    return (refData, compData)


# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-

def pdf_merger(output_path, input_paths):
    """Merges multiple pdf files into a single multi-page pdf file
    
    Parameters
    ----------
    output_path : str
        name of output multipage pdf file
        
    input_paths : list
        list of pdf files to combine
    
    Returns
    -------
    nothing.
    """
    pdf_merger = PdfFileMerger()

    for path in input_paths:
        pdf_merger.append(path)

    with open(output_path, 'wb') as fileobj:
        pdf_merger.write(fileobj)

    for path in input_paths:
        os.remove(path)


# =======================================================================================================================
if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description='Compare Sourcelists')
    # required positional input arguments
    PARSER.add_argument('sourcelistNames', nargs=2,
                        help='A space-separated pair of sourcelists to compare. The first sourcelist is assumed to be the reference sourcelist that the second is being compared to.')
    # optional input arguments
    PARSER.add_argument('-d', '--debugMode', required=False, choices=["True", "False"], default="False",
                        help="Turn on debug mode? Default value is False.")
    PARSER.add_argument('-i', '--imageNames', required=False, nargs=2,
                        help='A space-separated list of the fits images that were used to generate the input sourcelists. The first image corresponds to the first listed sourcelist, and so in. These will be used to improve the sourcelist alignment and matching.')
    PARSER.add_argument('-p', '--plotGen', required=False, choices=["screen", "file", "none"], default="none",
                        help='Generate Plots? "screen" displays plots on-screen. "file" saves them to a .pdf file, and "none" skips all plot generation.')
    PARSER.add_argument('-s', '--plotfile_prefix_string', required=False, default="",
                        help="text string that will prepend the plot files generated if plots are written to files ***REQUIRES the -p option set to 'file'***")
    PARSER.add_argument('-v', '--verbose', required=False, choices=["True", "False"], default="True",
                        help='Display verbose output? Default value is "True".')
    ARGS = PARSER.parse_args()
    print(ARGS)

    if ARGS.verbose == "True":
        ARGS.verbose = True
    else:
        ARGS.verbose = False

    if ARGS.debugMode == "True":
        ARGS.debugMode = True
    else:
        ARGS.debugMode = False

    runStatus = comparesourcelists(ARGS.sourcelistNames, ARGS.imageNames, plotGen=ARGS.plotGen,
                                   verbose=ARGS.verbose, log_level=logutil.logging.INFO, debugMode=ARGS.debugMode,
                                   plotfile_prefix=ARGS.plotfile_prefix_string)


# TODO: fix PEP 8 violations
# TODO: compute RA and Dec differences properly using astropy.skycoordinate (https://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html)
# TODO: Compute combined angular seperation based on both RA and DEC values for a given matched pair, then break it down into seperate componants
# TODO: Compute Magnitude differences properly
# TODO: Add flag-value based rejection (good flag bits = 0, 1, and/or 4, ignore lines containing flag bits other than these). This could also be an optional user-specified input.