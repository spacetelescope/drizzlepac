#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4 ai :

"""
This script compares two sourcelists and displays various measures of their differences. 3x3-sigma clipped mean, median and standard deviatnion, and  non-sigma clipped min and max values are computed for the following:

* X position
* Y position
* Right Ascension
* Declination
* Flux (Inner Aperture)
* Flux (Outer Aperture)
* Magnitude (Inner Aperture)
* Magnitude (Outer Aperture)

Bit-wise comparisons are also performed for the following item:

* Flag Value

.. note::
    Statistics (and optionally plots) for 'X position' and 'Y position' differences will always be displayed. However, not every sourcelist or catalog file is guaranteed to have any of the seven remaining data columns in the above list. Results for these seven columns will only be displayed if data columns are found in both input files.

Regression Testing
------------------
The following criteria must be met for the test to be declared "successful":

* X position: The sigma-clipped mean of all comparision - reference difference values is less than 1 sigma from zero.
* Y position: The sigma-clipped mean of all comparision - reference difference values is less than 1 sigma from zero.
* Right Ascension: The sigma-clipped mean of all comparision - reference difference values is less than 1 sigma from zero.
* Declination: The sigma-clipped mean of all comparision - reference difference values is less than 1 sigma from zero.
* Flux (Inner Aperture): The sigma-clipped mean of all comparision - reference difference values is less than 1 sigma from zero.
* Flux (Outer Aperture): The sigma-clipped mean of all comparision - reference difference values is less than 1 sigma from zero.
* Magnitude (Inner Aperture): The sigma-clipped mean of all comparision - reference difference values is less than 1 sigma from zero.
* Magnitude (Outer Aperture): The sigma-clipped mean of all comparision - reference difference values is less than 1 sigma from zero.
* Flag Value: The total number of differing flag bits is less than 5% of the total number of reference flag bits.

.. note::
    Sigma-clipped values for mean, sigma, and median are computed using the astropy.stats.sigma_clipped_stats() routine with three rounds of three-sigma clipping.

Path
----
HLApipeline/regression_testing/compare_sourcelists.py

Dependencies
------------
HLApipeline/regression_testing/starmatch_hist.py

Inputs
------
* Required input
    1: *sourcelistNames*
        * A space-separated pair of sourcelists to compare. The first sorucelist is assumed to be the reference sourcelist that the second is being compared to.

* Optional inputs:
    1: -i *imageNames*
        * A space-separated list of the fits images that were used to generate the input sourcelists. The first image corresponds to the first listed sourcelist, and so in. These will be used to imporove the sourcelist alignment and matching.

    2: -m *diffMode*
        * How should the comp-ref difference be calculated? "absolute" is simply the straight comp-ref difference. "peman" is the mean percent difference ((C-R)/avg(R)) x 100. "pdynamic" is the dynamic percent difference ((C-R)/R) x 100
        * Input choices: "absolute", "pmean" or "pdynamic"
        * Default value: "pmean"

    3: -p *plotGen*
        * Generate plots?
        * Input choices: "True" or "False"
        * Default value: False

    4: -p *verbose*
        * Display verbose output?
        * Input choices: "True" or "False"
        * Default value: True

Examples
--------

1: Compare j92c01b7q_single_sci.fits.coo to reference list j92c01b5q_single_sci.fits.coo using the corresponding single_sci.fits images to imporve alignment; display verbose output and generate plots::

    $HLApipeline/regression_testing/compare_sourcelists.py j92c01b5q_single_sci.fits.coo j92c01b7q_single_sci.fits.coo -i j92c01b5q_single_sci.fits j92c01b7q_single_sci.fits -p True

2: Compare V99.8/hst_10265_01_acs_wfc_f606w_daophot.txt to reference sourcelist V10.0/hst_10265_01_acs_wfc_f606w_daophot.txt. Return statistics on the differencs, but no plots or verbose output::

    $HLApipeline/regression_testing/compare_sourcelists.py V10.0/hst_10265_01_acs_wfc_f606w_daophot.txt V99.8/hst_10265_01_acs_wfc_f606w_daophot.txt -v False

Classes and Functions
---------------------
"""
import argparse
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import pdb
# from . import starmatch_hist
import starmatch_hist
import sys,os

from drizzlepac import util
from stsci.tools import fileutil, logutil

log =logutil.create_logger('compare_sourcelists', level=logutil.logging.INFO, stream=sys.stdout)
#-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
def computeFlagStats(matchedRA,plotGen,plot_title,verbose):
    """
    Compute and report statistics on the differences in flagging.

    :param matchedRA: A 2 x len(refLines) sized numpy array. Column 1: matched reference values. Column 2: The corresponding matched comparision values
    :param plotGen: Generate plots and display them to the screen (True/False)?
    :param plot_title: text string that will be used in plot title.
    :param verbose: display verbose output?
    :type matchedRA: numpy.ndarray
    :type plotGen: Boolean
    :type plot_title: string
    :type verbose: Boolean
    :return: Regression test status: Either "OK" or "FAILURE".
    """
    log.info(">>>>>> Comparision - reference sourcelist {} differences <<<<<<".format(plot_title))
    #set up arrays to count stuff up
    bit_list = [0, 1, 2, 4, 8, 16, 32, 64, 128]
    refFlagBreakdown=np.zeros(9,dtype=int)
    compFlagBreakdown=np.zeros(9,dtype=int)
    unchangedFlagBreakdown=np.zeros(9,dtype=int)
    on_off_FlagFlips=np.zeros(9,dtype=int)
    off_on_FlagFlips=np.zeros(9,dtype=int)
    for refFlagVal, compFlagVal in zip(matchedRA[0], matchedRA[1]):
        #break down each flag value into componant bit values, add values to totals
        refFlagRA=deconstruct_flag(refFlagVal)
        refFlagBreakdown += refFlagRA
        compFlagRA = deconstruct_flag(compFlagVal)
        compFlagBreakdown += compFlagRA
        #find differences in flagging, total up which bits were turned on, which were turned off.
        diffFlagRA=compFlagRA-refFlagRA
        if not np.array_equal(refFlagRA,compFlagRA):
            off_on_FlagFlips[np.where(diffFlagRA == 1)] += 1 #bits that are off in ref but on in comp
            on_off_FlagFlips[np.where(diffFlagRA == -1)] += 1 #bits that are on in ref but off in comp
            unchangedFlagBreakdown[np.where((refFlagRA ==1) & (compFlagRA ==1))] += 1 #takes care of the case were comp and ref have differing bits, but also have additional bits that are unchanged.
        if np.array_equal(refFlagRA, compFlagRA): #if there are absolutly no differences
            unchangedFlagBreakdown+=refFlagRA



    regTestStatus = "OK     "
    pct_diff_refbits=(np.sum([off_on_FlagFlips,on_off_FlagFlips],dtype=float)/np.sum(refFlagBreakdown, dtype=float))*100.0
    if pct_diff_refbits >=5.0:
        regTestStatus = "FAILURE"
    if ((verbose == True) or (regTestStatus == "FAILURE")):
        #Generate result tables
        n = np.sum(refFlagBreakdown, dtype=float)
        print("        Statistical Breakdown of Reference List Flagging Differences")
        print("  Flagging differences by number         Flagging differences by percentage")
        print("--------------------------------------------------------------------------------")
        print("%5s%9s%12s%10s%5s%5s %9s %12s %10s" % ("FLAG", "# TOTAL", "# UNCHANGED", "# ON->OFF","  |  ","FLAG", "% TOTAL", "% UNCHANGED", "% ON->OFF"))
        for ctr in range(0, len(bit_list)):
            print("%5d%9d%12d%10d%5s%5d  %8.4f %12.4f %10.4f" % (bit_list[ctr], refFlagBreakdown[ctr], unchangedFlagBreakdown[ctr],on_off_FlagFlips[ctr],"  |  ",bit_list[ctr],(float(refFlagBreakdown[ctr]) / n) * 100.0,(float(unchangedFlagBreakdown[ctr]) / n) * 100.0,(float(on_off_FlagFlips[ctr]) / n) * 100.0))
        print("%5s%9d%12d%10d%5s%5s  %8.4f %12.4f %10.4f"%("TOTAL",np.sum(refFlagBreakdown), np.sum(unchangedFlagBreakdown),np.sum(on_off_FlagFlips),"  |  ","TOTAL",(float(np.sum(refFlagBreakdown)) / n) * 100.0,(float(np.sum(unchangedFlagBreakdown)) / n) * 100.0,(float(np.sum(on_off_FlagFlips)) / n) * 100.0))
        print()
        n = np.sum(compFlagBreakdown,dtype=float)
        print("        Statistical Breakdown of Comparison List Flagging Differences")
        print("  Flagging differences by number         Flagging differences by percentage")
        print("--------------------------------------------------------------------------------")
        print("%5s%9s%12s%10s%5s%5s %9s %12s %10s" % ("FLAG", "# TOTAL", "# UNCHANGED", "# OFF->ON", "  |  ", "FLAG", "% TOTAL", "% UNCHANGED", "% OFF->ON"))
        for ctr in range(0, len(bit_list)):
            print("%5d%9d%12d%10d%5s%5d  %8.4f %12.4f %10.4f" % (
            bit_list[ctr], compFlagBreakdown[ctr], unchangedFlagBreakdown[ctr], off_on_FlagFlips[ctr], "  |  ",bit_list[ctr], (float(compFlagBreakdown[ctr]) / n) * 100.0, (float(unchangedFlagBreakdown[ctr]) / n) * 100.0,(float(off_on_FlagFlips[ctr]) / n) * 100.0))
        print("%5s%9d%12d%10d%5s%5s  %8.4f %12.4f %10.4f" % ("TOTAL", np.sum(compFlagBreakdown), np.sum(unchangedFlagBreakdown), np.sum(off_on_FlagFlips), "  |  ","TOTAL", (float(np.sum(compFlagBreakdown)) / n) * 100.0, (float(np.sum(unchangedFlagBreakdown)) / n) * 100.0,(float(np.sum(off_on_FlagFlips)) / n) * 100.0))
        print()
        print("Total flag bit differences......... ",np.sum([off_on_FlagFlips,on_off_FlagFlips]))
        print("Percentage change of all ref bits.. ",pct_diff_refbits,"%")

    print("Regression test status............. ",regTestStatus)



    if plotGen != "none":
        idx=np.arange(9)
        x_title_list=[]
        for bit in bit_list:x_title_list.append(str(bit))

        #plot flag breakdown by bit for all matched sources in the reference and comparison sourcelists
        width = 0.35
        p1=plt.bar(idx - width/2,refFlagBreakdown,width,label='Reference')
        p2 = plt.bar(idx + width/2, compFlagBreakdown,width, label='Comparison')
        plt.legend()
        plt.xticks(idx, x_title_list)
        plt.xlabel("Flag Bit Value")
        plt.ylabel("Number of matched sources")
        fullPlotTitle = "Flag Breakdown by Bit"
        plt.title(fullPlotTitle)

        if plotGen == "screen":
            plt.show()
        if plotGen == "file":
            plotFileName = fullPlotTitle.replace(" ","_")+".pdf"
            plt.savefig(plotFileName)
            plt.close()
            print("{} plot saved to file {}.".format(fullPlotTitle, plotFileName))

        #plot flag changes brokeken down by bit
        p_unchanged=plt.bar(idx,unchangedFlagBreakdown)
        p_offOn=plt.bar(idx,off_on_FlagFlips,bottom=unchangedFlagBreakdown)
        p_onOff=plt.bar(idx,on_off_FlagFlips,bottom=off_on_FlagFlips+unchangedFlagBreakdown)
        plt.xticks(idx,x_title_list)
        plt.xlabel("Flag Bit Value")
        plt.ylabel("Number of matched sources")
        fullPlotTitle = "Flagging Differences by Bit"
        plt.title(fullPlotTitle)
        plt.legend((p_onOff[0],p_offOn,p_unchanged[0]),("On -> Off","Off -> On","Unchanged"))

        if plotGen == "screen":
            plt.show()
        if plotGen == "file":
            plotFileName = fullPlotTitle.replace(" ","_")+".pdf"
            plt.savefig(plotFileName)
            plt.close()
            print("{} plot saved to file {}.".format(fullPlotTitle, plotFileName))

    regTestStatus = "%s %11.7f"%(regTestStatus,pct_diff_refbits)
    return (regTestStatus)
# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
def computeLinearStats(matchedRA,plotGen,diffMode,plot_title,verbose):
    """
    Compute stats on the quantities with differences that can be computed with simple subtraction (X, Y, RA, Dec, Flux, and Magnitude).

    :param matchedRA: A 2 x len(refLines) sized numpy array. Column 1: matched reference values. Column 2: The corresponding matched comparision values
    :param plotGen: Generate plots?
    :param diffMode: Method used to compute diffRA
    :param plot_title: text string that will be used in plot title.
    :param verbose: display verbose output?
    :type matchedRA: numpy.ndarray
    :type plotGen: string
    :type diffMode: string
    :type plot_title: string
    :type verbose: Boolean
    :return: Regression test status: Either "OK" or "FAILURE".
    """
    print(">>>>>> Comparision - reference sourcelist %s differences <<<<<<" % (plot_title))
    #remove any "inf" or "nan" values in matchedRA.
    nanIdx = np.where(np.isnan(matchedRA) == True)[1]
    if len(nanIdx) >0:
        if verbose: print("{} np.nan values will be removed from input list. New input list contains {} values.".format(len(nanIdx),len(matchedRA[0,:]) - len(nanIdx)))
        matchedRA = np.delete(matchedRA, nanIdx, axis=1)

    infIdx = np.where(np.isinf(matchedRA) == True)[1]
    if len(infIdx) >0:
        if verbose: print("{} np.inf values will be removed from input list. New input list contains {} values.".format(len(infIdx),len(matchedRA[0,:]) - len(infIdx)))
        matchedRA = np.delete(matchedRA, infIdx, axis=1)

    # 'sigma' and 'iters' input values used for various np.sigma_clipped_stats() runs
    sigVal= 3
    intersVal= 3

    if diffMode == "absolute":
        diffRA=matchedRA[1, :] - matchedRA[0, :] # simple difference

    else:
        if diffMode == "pmean":
            normValue = np.abs(sigma_clipped_stats(matchedRA[0, :], sigma=sigVal, maxiters=intersVal)[0]) #divide all comp-ref values by a single sigma-clipped mean ref value
            if verbose: print("normValue: ", normValue)
        if diffMode == "pdynamic":
            normValue = matchedRA[0, :] #divide each comp-ref value by the corresponding ref value

        diffRA = ((matchedRA[1, :] - matchedRA[0, :]) / normValue) * 100.0

    clippedStats = sigma_clipped_stats(diffRA, sigma=sigVal, maxiters=intersVal)
    pct_1sig = (float(np.shape(np.where(abs(diffRA) <= clippedStats[2]))[1]) / float(np.shape(diffRA)[0])) * 100.0
    pct_five = (float(np.shape(np.where(abs(diffRA) <= 5.0))[1]) / float(np.shape(diffRA)[0])) * 100.0

    out_stats="%11.7f %11.7f  %11.7f  %11.7f  %11.7f "%(clippedStats[0],clippedStats[1],clippedStats[2],pct_five,pct_1sig)
    if diffMode == "absolute":
        if abs(clippedStats[0]) <= abs(clippedStats[2]): #success conditon: sigma clippped mean less then 1 sigma from Zero.
            regTestStatus = "OK      "
        else: regTestStatus = "FAILURE "
    else:
        if np.abs(clippedStats[0]) <= 5.0: #success condition: sigma-clipped mean less then 5%
            regTestStatus =   "OK      "
        else:
            regTestStatus = "FAILURE "
    if ((verbose == True) or (regTestStatus == "FAILURE ")):
        print("       Sigma-clipped Statistics; Sigma = {}, # steps = {}".format(sigVal,intersVal))
        print("Sigma-clipped mean........................ ",clippedStats[0])
        print("Sigma-clipped median...................... ",clippedStats[1])
        print("Sigma-clipped standard deviation.......... ",clippedStats[2])
        print("Sigma-clipped mean in units of SD......... ",np.divide(clippedStats[0],clippedStats[2]))
        print()
        print()
        print("            Non-Clipped Statistics")
        print("Non-clipped mean.......................... ",np.mean(diffRA))
        print("Non-clipped median........................ ",np.median(diffRA))
        print("Non-clipped standard deviation............ ",np.std(diffRA))
        print("Non-clipped mean in units of SD........... ",np.divide(np.mean(diffRA), np.std(diffRA)))
        print("Non-clipped minimum....................... ",np.min(diffRA))
        print("Non-clipped maximum....................... ",np.max(diffRA))
        print("% all diff values within 1 sigma of 0.0... ", pct_1sig)
        print("% all diff values within 5% of 0.0........ ",pct_five)
    print("Regression test status.................... ",regTestStatus)

    if plotGen != "none":
        # plt.hist(diffRA,bins='auto')
        # plt.axvline(x=clippedStats[0], color='k', linestyle='--')
        # plt.axvline(x=clippedStats[0] + clippedStats[2], color='r', linestyle=':')
        # plt.axvline(x=clippedStats[0] - clippedStats[2], color='r', linestyle=':')
        # plt.xlabel("$\Delta %s$"%(plot_title.split(" ")[0]))
        # plt.ylabel("Number of matched sources")
        # plt.title("Comparision - reference sourcelist %s differences"%(plot_title))
        # plt.show()

        if diffMode.startswith("p"): xAxisString="{} (percent)".format(plot_title.split(" ")[0])
        else: xAxisString="{}".format(plot_title.split(" ")[0])
        plotCutoff=(10.0*np.abs(clippedStats[2]))+np.abs(clippedStats[0])
        if plotCutoff != 0.0:
            origSize=len(diffRA)
            print("Plot cutoff: ",plotCutoff)
            goodIdx=np.where(np.abs(diffRA)<=plotCutoff)
            diffRA=diffRA[goodIdx]
            print("%d values (%7.4f percent) clipped from plot."%(origSize-len(diffRA),(float(origSize-len(diffRA))/float(origSize))*100.0))
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        fullPlotTitle = "Comparision - reference sourcelist %s differences" % (plot_title)
        plt.title(fullPlotTitle)
        ax1.hist(diffRA,bins='auto')
        ax1.axvline(x=clippedStats[0], color='k', linestyle='--')
        ax1.axvline(x=clippedStats[0] + clippedStats[2], color='k', linestyle=':')
        ax1.axvline(x=clippedStats[0] - clippedStats[2], color='k', linestyle=':')

        ax1.axvline(x=np.mean(diffRA)+np.std(diffRA), color='g', linestyle=':')
        ax1.axvline(x=np.mean(diffRA)-np.std(diffRA), color='g', linestyle=':')
        ax1.set_xlabel("$\Delta %s$" % (xAxisString))
        ax1.set_ylabel("Number of matched sources")

        ax2 = ax1.twinx()
        ax2.hist(diffRA, bins='auto', cumulative=-1, density=True, histtype='step', color='r')
        ax2.set_ylabel("Fraction of all matched sources",color='r')
        for tl in ax2.get_yticklabels():
            tl.set_color('r')
        if plotGen == "screen":
            plt.show()
        if plotGen == "file":
            plotFileName = plot_title.replace(" ","_")+".pdf"
            plt.savefig(plotFileName)
            plt.close()
            print("{} plot saved to file {}.".format(fullPlotTitle, plotFileName))
    print()
    return(regTestStatus+out_stats)
# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
def deconstruct_flag(flagval):
    """
    Breaks down an integer flag value into individual component bit values.

    :param flagval: Flag value to deconstruct
    :type flagval: integer
    :return: a 9-element numpy array of 0s and 1s. Each element of the array represents the presence of a particular bit value (element 0 = bit 0, element 1 = bit 1, ..., element 3 = bit 4 and so on...)
    """
    bit_list = [1, 2, 4, 8, 16, 32, 64, 128]
    flagval = int(flagval)
    #out_bit_list = []
    out_idx_list = np.zeros(9,dtype=int)
    if flagval == 0:
        #out_bit_list = [0]
        out_idx_list[0] =1
    if flagval > 0:
        idx =1
        for bit in bit_list:
            if flagval & bit > 0:
                #out_bit_list.append(bit)
                out_idx_list[idx] = 1
            if bit > flagval: break
            idx += 1
    return (out_idx_list)
#-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
def extractMatchedLines(col2get,refData,compData,refLines,compLines):
    """
    Extracts only matching lines of data for a specific column of refData and compData. Returns empty list if the specified column is not found in both tables.

    :param col2get: Title of the column to return
    :param refData: reference data table
    :param compData: comparison data table
    :param refLines: List of matching refData line numbers
    :param compLines: List of matching compData line numbers
    :type col2get: string
    :type refData: astropy Table object
    :type compData: astropy Table object
    :type refLines: numpy.ndarray
    :type compLines: numpy.ndarray
    :return: A 2 x len(refLines) sized numpy array. Column 1: matched reference values. Column 2: The corresponding matched comparision values
    """
    return_ra=[]
    if col2get in list(refData.keys()) and col2get in list(compData.keys()):
        matching_refData = refData[col2get][refLines].data
        matching_compData = compData[col2get][compLines].data
        return_ra=np.stack((matching_refData,matching_compData))
    return(return_ra)
#-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
def getMatchedLists(slNames,imgNames,slLengths):
    """
    run starmatch_hist to get the indices of matching sources that are common to both input source catalogs

    :param slNames: list of input source lists
    :param imgNames: list of input images
    :param slLengths: list of integer sourcelist lengths
    :type slNames: list
    :type imgNames: list
    :type slLengths: list
    :return: A list of the indices of reference sourcelist sources that match comparison sourcelist sources , and a corresponding list of the indices of comparison sourcelist sources that match reference sourcelist sources
    """
    source_list_dict = {}
    equal_flag=False
    if slLengths[0] == slLengths[1]:
        slLengths[0] +=1
        equal_flag=True

    for ctr, slName in enumerate(slNames, 0):
        source_list_dict[slName] = slLengths[ctr]

    try:
        fh = fits.open(imgNames[0])
        xref = fh[0].header['crpix1']
        yref = fh[0].header['crpix2']
        fh.close()
    except:
        print("WARNING: Unable to fetch values for xref and yref from fits file headers. Using xref = 0.0 and yref = 0.0.")
        xref=0.0
        yref=0.0
    print("source_list_dict: ",source_list_dict)
    out_dict = starmatch_hist.run(source_list_dict, xref=xref, yref=yref)
    matching_lines_ref = out_dict[slNames[0]]
    matching_lines_img = out_dict[slNames[1]]

    if equal_flag: slLengths[0] -=1
    #Report number and percentage of the total number of detected ref and comp sources that were matched
    print("Sourcelist Matching Results")
    print("Reference sourcelist:  {} of {} total sources matched ({} %)".format(len(matching_lines_ref),slLengths[0],100.0*(float(len(matching_lines_ref))/float(slLengths[0]))))
    print("Comparison sourcelist: {} of {} total sources matched ({} %)".format(len(matching_lines_img),slLengths[1],100.0*(float(len(matching_lines_img))/float(slLengths[1]))))

    return(matching_lines_ref,matching_lines_img)
#-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
def makeVectorPlot(x,y,plotDest,binThresh = 10000,binSize=250):
    """
    Generate vector plot of dx and dy values vs. reference (x,y) positions

    :param x: A 2 x n sized numpy array. Column 1: matched reference X values. Column 2: The corresponding matched comparision X values
    :param y: A 2 x n sized numpy array. Column 1: matched reference Y values. Column 2: The corresponding matched comparision Y values
    :param plotDest: plot destination; screen or file
    :param binThresh: Minimum size of list *x* and *y* that will trigger generation of a binned vector plot. Default value = 10000.
    :param binSize: Size of binning box in pixels. When generating a binned vector plot, mean dx and dy values are computed by taking the mean of all points located within the box. Default value = 250.
    :type x: numpy.ndarray
    :type y: numpy.ndarray
    :type plotDest: string
    :type binThresh: integer
    :type binSize: integer
    :return: nothing
    """
    dx = x[1, :] - x[0, :]
    dy = y[1, :] - y[0, :]
    if len(dx)>binThresh:# if the input list is larger than binThresh, a binned vector plot will be generated.
        binStatus = "%d x %d Binning"%(binSize,binSize)
        print("Input list length greater than threshold length value. Generating binned vector plot using %d pixel x %d pixel bins"%(binSize,binSize))
        if min(x[0,:])<0.0: xmin=min(x[0,:])
        else: xmin = 0.0
        if min(y[0,:])<0.0: ymin=min(y[0,:])
        else: ymin = 0.0

        p_x=np.empty(shape=[0])
        p_y=np.empty(shape=[0])
        p_dx=np.empty(shape=[0])
        p_dy=np.empty(shape=[0])
        color_ra=[]
        for xBinCtr in range(int(round2ArbatraryBase(xmin,"down",binSize)),int(round2ArbatraryBase(max(x[0,:]),"up",binSize)),binSize):
            for yBinCtr in range(int(round2ArbatraryBase(ymin, "down", binSize)),
                                 int(round2ArbatraryBase(max(y[0, :]), "up", binSize)), binSize):
                #define bin box x,y upper and lower bounds
                xBinMin=xBinCtr
                xBinMax=xBinMin+binSize
                yBinMin=yBinCtr
                yBinMax=yBinMin+binSize
                #get indicies of x and y withen bounding box
                ix0 = np.where((x[0,:] >= xBinMin) & (x[0,:] < xBinMax) & (y[0,:] >= yBinMin) & (y[0,:] < yBinMax))
                if len(dx[ix0]) > 0 and len(dy[ix0]) > 0: #ignore empty bins
                    p_x=np.append(p_x, xBinCtr + 0.5 * binSize) #X and Y posotion at center of bin.
                    p_y=np.append(p_y, yBinCtr + 0.5 * binSize)
                    mean_dx=np.mean(dx[ix0])
                    p_dx=np.append(p_dx, mean_dx) #compute mean dx, dy values
                    mean_dy = np.mean(dy[ix0])
                    p_dy=np.append(p_dy,mean_dy)
                    avg_npts=(float(len(dx[ix0]))+float(len(dy[ix0])))/2.0 #keep an eye out for mean values computed from less than 10 samples.
                    if (avg_npts<10.0): #if less than 10 samples were used in mean calculation, color the vector red.
                        color_ra.append('r')
                    else:
                        color_ra.append('k')
        lowSampleWarning = ""
        if "r" in color_ra: lowSampleWarning = "; Red Vectors were computed with less than 10 values"
    else:
        print("Generating unbinned vector plot")
        binStatus = "Unbinned"
        lowSampleWarning = ""
        color_ra=["k"]
        p_x=x[0,:]
        p_y = y[0, :]
        p_dx=dx
        p_dy=dy
    plt_mean = np.mean(np.hypot(p_dx, p_dy))
    e = np.log10(5.0*plt_mean).round()
    plt_scaleValue=10**e
    if len(dx) > binThresh: Q = plt.quiver(p_x, p_y, p_dx, p_dy,color=color_ra,units="xy")
    else: Q = plt.quiver(p_x, p_y, p_dx, p_dy)
    plt.quiverkey(Q, 0.75, 0.05, plt_scaleValue, r'%5.3f'%(plt_scaleValue), labelpos='S', coordinates='figure', color="k")
    plot_title="Comparision - reference $\Delta X$, $\Delta Y$ values vs. $(X_{ref}, Y_{ref})$ positions\n%s%s" % (binStatus,lowSampleWarning)
    plt.title(plot_title)
    plt.xlabel(r"$X_{ref}$")
    plt.ylabel(r"$Y_{ref}$")
    if plotDest == "screen":
        plt.show()
    if plotDest == "file":
        plt.savefig("xy_vector_plot.pdf")
        plt.close()
        print("Vector plot saved to file xy_vector_plot.pdf")
# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
def round2ArbatraryBase(value,direction,roundingBase):
    """
    Round value up or down to arbitrary base

    :param value: Value to be rounded.
    :param direction: Round up, down to the nearest base (choices: "up","down","nearest")
    :param roundingBase: rounding base. (Example: if base = 5, values will be rounded to the nearest multiple of 5 or 10.)
    :type value: float
    :type direction: string
    :type roundingBase: integer
    :return: rounded value.
    """
    if direction.lower().startswith("u"):
        rv=value+(roundingBase-value%roundingBase) #round up to nearest base
    elif direction.lower().startswith("d"):
        rv=value-value%roundingBase #round down to nearest base
    else:
        rv=int(roundingBase * round(float(value)/roundingBase)) #round up or down to nearest base
    return rv
#-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
@util.with_logging
def comparesourcelists(slNames,imgNames,plotGen,diffMode,verbose):
    """
    Main calling subroutine to compare sourcelists.

    :param slNames: list of input source lists
    :param imgNames: optional list of input images that starmatch_hist will use to improve sourcelist matching
    :param plotGen: Generate plots and display them to the screen (True/False)?
    :param diffMode: method used to compute comp-ref difference in computeLinearStats().
    :param verbose: display verbose output?
    :type slNames: list
    :type imgNames: list
    :type plotGen: Boolean
    :type diffMode: string
    :type verbose: Boolean
    :return: "OK" if all tests were passed, or "FAILURE" if inconsistencies were found.
    """
    regressionTestResults={}
    colTitles=[]
    # 1: Read in sourcelists fiels into astropy table or 2-d array so that individual columns from each sourcelist can be easily accessed later in the code.
    refData,compData=slFiles2dataTables(slNames)
    print("Valid reference data columns:   ",list(refData.keys()))
    print("Valid comparision data columns: ",list(compData.keys()))
    print()
    print("Data columns to be compared:")
    for listItem in sorted(list(set(refData.keys()).intersection(set(compData.keys())))): print(listItem)
    print()
    # 2: Run starmatch_hist to get list of matched sources common to both input sourcelists
    slLengths=[len(refData['X']),len(compData['X'])]
    matching_lines_ref, matching_lines_img=getMatchedLists(slNames,imgNames,slLengths)
    # 3: Compute and display statistics on X position differences for matched sources
    matched_values=extractMatchedLines("X",refData,compData,matching_lines_ref, matching_lines_img)
    if len(matched_values) >0:
        rt_status=computeLinearStats(matched_values,plotGen,diffMode,"X position",verbose)
        regressionTestResults["X Position"]=rt_status
        colTitles.append("X Position")
        if plotGen != "none":matchedXValues=matched_values
    # 4: Compute and display statistics on Y position differences for matched sources
    matched_values=extractMatchedLines("Y",refData,compData,matching_lines_ref, matching_lines_img)
    if len(matched_values) >0:
        rt_status=computeLinearStats(matched_values,plotGen,diffMode,"Y position",verbose)
        regressionTestResults["Y Position"]=rt_status
        colTitles.append("Y Position")
        if plotGen != "none" and diffMode == "absolute":
            matchedYValues=matched_values
            makeVectorPlot(matchedXValues,matchedYValues,plotGen)

    # 5: Compute and display statistics on RA position differences for matched sources
    matched_values=extractMatchedLines("RA",refData,compData,matching_lines_ref, matching_lines_img)
    if len(matched_values) >0:
        rt_status=computeLinearStats(matched_values,plotGen,diffMode,"RA position",verbose)
        regressionTestResults["RA Position"]=rt_status
        colTitles.append("RA Position")

    # 6: Compute and display statistics on DEC position differences for matched sources
    matched_values=extractMatchedLines("DEC",refData,compData,matching_lines_ref, matching_lines_img)
    if len(matched_values) >0:
        rt_status=computeLinearStats(matched_values,plotGen,diffMode,"DEC position",verbose)
        regressionTestResults["DEC Position"]=rt_status
        colTitles.append("DEC Position")

    # 7: Compute and display statistics on flux differences for matched sources
    matched_values=extractMatchedLines("FLUX1",refData,compData,matching_lines_ref, matching_lines_img)
    if len(matched_values) >0:
        rt_status=computeLinearStats(matched_values,plotGen,diffMode,"Flux (Inner Aperture)",verbose)
        regressionTestResults["Flux (Inner Aperture)"]=rt_status
        colTitles.append("Flux (Inner Aperture)")

    matched_values=extractMatchedLines("FLUX2",refData,compData,matching_lines_ref, matching_lines_img)
    if len(matched_values) >0:
        rt_status=computeLinearStats(matched_values,plotGen,diffMode,"Flux (Outer Aperture)",verbose)
        regressionTestResults["Flux (Outer Aperture)"]=rt_status
        colTitles.append("Flux (Outer Aperture)")

    # 8: Compute and display statistics on magnitude differences for matched sources
    matched_values=extractMatchedLines("MAGNITUDE1",refData,compData,matching_lines_ref, matching_lines_img)
    if len(matched_values) >0:
        rt_status=computeLinearStats(matched_values,plotGen,diffMode,"Magnitude (Inner Aperture)",verbose)
        regressionTestResults["Magnitude (Inner Aperture)"]=rt_status
        colTitles.append("Magnitude (Inner Aperture)")

    matched_values=extractMatchedLines("MERR1",refData,compData,matching_lines_ref, matching_lines_img)
    if len(matched_values) >0:
        formalTitle = "Magnitude (Inner Aperture) Error"
        rt_status=computeLinearStats(matched_values,plotGen,diffMode,formalTitle,verbose)
        regressionTestResults[formalTitle]=rt_status
        colTitles.append(formalTitle)

    matched_values=extractMatchedLines("MAGNITUDE2",refData,compData,matching_lines_ref, matching_lines_img)
    if len(matched_values) >0:
        rt_status=computeLinearStats(matched_values,plotGen,diffMode,"Magnitude (Outer Aperture)",verbose)
        regressionTestResults["Magnitude (Outer Aperture)"]=rt_status
        colTitles.append("Magnitude (Outer Aperture)")

    matched_values=extractMatchedLines("MERR2",refData,compData,matching_lines_ref, matching_lines_img)
    if len(matched_values) >0:
        formalTitle = "Magnitude (Outer Aperture) Error"
        rt_status=computeLinearStats(matched_values,plotGen,diffMode,formalTitle,verbose)
        regressionTestResults[formalTitle]=rt_status
        colTitles.append(formalTitle)

    # 9: Compute and display statistics on differences in background sky values
    matched_values=extractMatchedLines("MSKY",refData,compData,matching_lines_ref, matching_lines_img)
    if len(matched_values) >0:
        formalTitle = "MSKY value"
        rt_status=computeLinearStats(matched_values,plotGen,diffMode,formalTitle,verbose)
        regressionTestResults[formalTitle]=rt_status
        colTitles.append(formalTitle)

    matched_values=extractMatchedLines("STDEV",refData,compData,matching_lines_ref, matching_lines_img)
    if len(matched_values) >0:
        formalTitle = "STDEV value"
        rt_status=computeLinearStats(matched_values,plotGen,diffMode,formalTitle,verbose)
        regressionTestResults[formalTitle]=rt_status
        colTitles.append(formalTitle)

    # 10: Compute and display statistics on differences in concentration index  for matched sources
    matched_values=extractMatchedLines("CI",refData,compData,matching_lines_ref, matching_lines_img)
    if len(matched_values) >0:
        formalTitle = "CI"
        rt_status=computeLinearStats(matched_values,plotGen,diffMode,formalTitle,verbose)
        regressionTestResults[formalTitle]=rt_status
        colTitles.append(formalTitle)

    # 11: Compute and display statistics on differences in flag populations for matched sources
    matched_values=extractMatchedLines("FLAGS",refData,compData,matching_lines_ref, matching_lines_img)
    if len(matched_values) >0:
        rt_status=computeFlagStats(matched_values,plotGen,"Source Flagging",verbose)
        regressionTestResults["Source Flagging"]=rt_status
        colTitles.append("Source Flagging")
    overallStatus="OK"
    print("\n"*2)
    print(">---> REGRESSION TESTING SUMMARY <---<")
    print("                                                                                 % within     % within")
    print("COLUMN                             STATUS   MEAN        MEDIAN       STD DEV     5% of 0.     1 SD of 0.")
    lenList=[]
    for item in colTitles:
        lenList.append(len(item))
    totalPaddedSize=max(lenList)+3

    for colTitle in colTitles:
        print("%s%s%s"%(colTitle,"."*(totalPaddedSize-len(colTitle)),regressionTestResults[colTitle]))
        if not regressionTestResults[colTitle].startswith("OK"):overallStatus="FAILURE"
    return(overallStatus)
#-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
def slFiles2dataTables(slNames):
    """
    Reads in data from sourcelists, returns some or all of the following data columns in a pair of properly formatted astropy.table objects:

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
        'X position' and 'Y position' data columns will always be returned. However, not every sourcelist or catalog file is guaranteed to have any of the seven remaining data columns in the above list.

    :param slNames: A list containing the reference sourcelist filename and the comparison sourcelist filename, in that order.
    :type slNames: list
    :return: a pair of astropy.table objects the first is the data from the reference sourcelist, and the second is from the comparison sourcelist.
    """
    if slNames[0].endswith(".ecsv"):
        refData_in = Table.read(slNames[0], format='ascii.ecsv')
    else:
        try:
            refData_in=Table.read(slNames[0], format='ascii.daophot')
        except:
            refData_in = Table.read(slNames[0], format='ascii')
    if slNames[1].endswith(".ecsv"):
        compData_in = Table.read(slNames[1], format='ascii.ecsv')
    else:
        try:
            compData_in = Table.read(slNames[1], format='ascii.daophot')
        except:
            compData_in=Table.read(slNames[1], format='ascii')
    titleSwapDict_dao1={"X":"X-Center","Y":"Y-Center","RA":"RA","DEC":"DEC","FLUX1":"n/a","FLUX2":"Flux(0.15)","MAGNITUDE1":"MagAp(0.05)","MAGNITUDE2":"MagAp(0.15)","MERR1":"MagErr(0.05)","MERR2":"MagErr(0.15)","MSKY":"MSky(0.15)","STDEV":"Stdev(0.15)","FLAGS":"Flags","ID":"ID","CI":"CI"}
    titleSwapDict_dao2 = {"X": "X-Center", "Y": "Y-Center", "RA": "RA", "DEC": "DEC", "FLUX1": "n/a",
                          "FLUX2": "Flux(0.45)", "MAGNITUDE1": "MagAp(0.15)", "MAGNITUDE2": "MagAp(0.45)",
                          "MERR1": "MagErr(0.15)", "MERR2": "MagErr(0.45)", "MSKY": "MSky(0.45)",
                          "STDEV": "Stdev(0.45)", "FLAGS": "Flags", "ID": "ID", "CI": "CI"}
    titleSwapDict_dao3 = {"X": "X-Center", "Y": "Y-Center", "RA": "RA", "DEC": "DEC", "FLUX1": "n/a",
                          "FLUX2": "Flux(0.125)", "MAGNITUDE1": "MagAp(0.03)", "MAGNITUDE2": "MagAp(0.125)",
                          "MERR1": "MagErr(0.03)", "MERR2": "MagErr(0.125)", "MSKY": "MSky(0.125)",
                          "STDEV": "Stdev(0.125)", "FLAGS": "Flags", "ID": "ID", "CI": "CI"}
    titleSwapDict_daoTemp={"X":"XCENTER","Y":"YCENTER","RA":"n/a","DEC":"n/a","FLUX1":"FLUX1","FLUX2":"FLUX2","MAGNITUDE1":"MAG1","MAGNITUDE2":"MAG2","FLAGS":"n/a","ID":"ID","MERR1":"MERR1","MERR2":"MERR2","MSKY":"MSKY","STDEV":"STDEV"}
    titleSwapDict_sourceX={"X":"X_IMAGE","Y":"Y_IMAGE","RA":"RA","DEC":"DEC","FLUX1":"FLUX_APER1","FLUX2":"FLUX_APER2","MAGNITUDE1":"MAG_APER1","MAGNITUDE2":"MAG_APER2","FLAGS":"FLAGS","ID":"NUMBER"}
    titleSwapDict_cooNew={"X":"col1","Y":"col2","RA":"n/a","DEC":"n/a","FLUX1":"n/a","FLUX2":"n/a","MAGNITUDE1":"n/a","MAGNITUDE2":"n/a","FLAGS":"n/a","ID":"col7"}
    titleSwapDict_cooOld = {"X":"XCENTER","Y":"YCENTER","RA":"n/a","DEC":"n/a","FLUX1":"n/a","FLUX2":"n/a","MAGNITUDE1":"n/a","MAGNITUDE2":"n/a","FLAGS":"n/a","ID":"ID"}

    titleSwapDict_cooOld2 = {"X": "XCENTER", "Y": "YCENTER", "RA": "RA", "DEC": "DEC", "FLUX1": "FLUX_0.05", "FLUX2": "FLUX_0.15",
                            "MAGNITUDE1": "MAG_0.05", "MAGNITUDE2": "MAG_0.15", "FLAGS": "n/a", "ID": "ID","MERR1":"MERR_0.05","MERR2":"MERR_0.15","MSKY":"MSKY","STDEV":"STDEV"}


    titleSwapDict_daorep = {"X": "X", "Y": "Y", "RA": "n/a", "DEC": "n/a", "FLUX1": "flux_0", "FLUX2": "flux_1","MAGNITUDE1": "mag_0", "MAGNITUDE2": "mag_1", "FLAGS": "n/a", "ID": "n/a"}
    ctr=1
    for dataTable in [refData_in,compData_in]:
        if "X-Center" in list(dataTable.keys()):
            if (("MagAp(0.05)" in list(dataTable.keys())) and ("MagAp(0.15)" in list(dataTable.keys()))): #ACS/WFC, WFC3/UVIS
                print("titleSwapDict_dao1")
                titleSwapDict = titleSwapDict_dao1
            elif (("MagAp(0.15)" in list(dataTable.keys())) and ("MagAp(0.45)" in list(dataTable.keys()))): #WFC3/IR
                print("titleSwapDict_dao2")
                titleSwapDict = titleSwapDict_dao2
            elif (("MagAp(0.03)" in list(dataTable.keys())) and ("MagAp(0.125)" in list(dataTable.keys()))): #ACS/HRC
                print("titleSwapDict_dao3")
                titleSwapDict = titleSwapDict_dao3
            else:
                sys.exit("ERROR: Unrecognized format. Exiting...")
        elif ("XCENTER" in list(dataTable.keys()) and "FLUX1" in list(dataTable.keys())):
            print("titleSwapDict_daoTemp")
            titleSwapDict = titleSwapDict_daoTemp
        elif "X_IMAGE" in list(dataTable.keys()):
            print("titleSwapDict_sourceX")
            titleSwapDict = titleSwapDict_sourceX
        elif "col1" in list(dataTable.keys()):
            print("titleSwapDict_cooNew")
            titleSwapDict = titleSwapDict_cooNew
        elif "XCENTER" in list(dataTable.keys()):
            print("titleSwapDict_cooOld2")
            titleSwapDict = titleSwapDict_cooOld2
        elif "X" in list(dataTable.keys()):
            print("titleSwapDict_daorep")
            titleSwapDict = titleSwapDict_daorep
        else: sys.exit("ERROR: Unrecognized format. Exiting...")
        outTable=Table()
        for swapKey in list(titleSwapDict.keys()):
            if titleSwapDict[swapKey] != "n/a":
                col2add=Table.Column(name=swapKey,data=dataTable[titleSwapDict[swapKey]])
                outTable.add_column(col2add)
        if ctr == 1: refData = outTable
        if ctr == 2: compData = outTable
        ctr+=1

    return(refData, compData)
#=======================================================================================================================
if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description='Compare Sourcelists')
    # required positional input arguments
    PARSER.add_argument('sourcelistNames', nargs=2,help='A space-separated pair of sourcelists to compare. The first sorucelist is assumed to be the reference sourcelist that the second is being compared to.')
    # optional input arguments
    PARSER.add_argument('-i', '--imageNames', required = False, nargs=2,help='A space-seperated list of the fits images that were used to generate the input sourcelists. The first image corresponds to the first listed sourcelist, and so in. These will be used to imporove the sourcelist alignment and matching.')
    PARSER.add_argument('-m', '--diffMode', required=False, choices=["absolute", "pmean","pdynamic"], default="pmean",
                        help='How should the comp-ref difference be calculated? "absolute" is simply the stright comp-ref difference. "peman" is the mean percent difference ((C-R)/avg(R)) x 100. "pdynamic" is the dynamic percent difference ((C-R)/R) x 100. Default value is "pmean".')
    PARSER.add_argument('-p', '--plotGen', required=False, choices=["screen","file","none"], default="none",help='Generate Plots? "screen" displays plots on-screen. "file" saves them to a .pdf file, and "none" skips all plot generation.')
    PARSER.add_argument('-v', '--verbose', required=False, choices=["True", "False"], default="True",
                        help='Display verbose output? Default value is "True".')
    ARGS = PARSER.parse_args()
    print(ARGS)

    if ARGS.verbose == "True":
        ARGS.verbose = True
    else: ARGS.verbose = False

    runStatus=comparesourcelists(ARGS.sourcelistNames,ARGS.imageNames,ARGS.plotGen,ARGS.diffMode,ARGS.verbose)
