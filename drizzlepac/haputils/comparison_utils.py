#!/usr/bin/env python
"""A collection of functions that assist with sourcelist comparison"""

# Standard library imports
import os
import sys

# Related third party imports
from astropy.table import Table
import numpy as np
from PyPDF2 import PdfFileMerger

# Local application imports
from drizzlepac.haputils import starmatch_hist
from stsci.tools import logutil

__taskname__ = 'comparison_utils'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)


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
        a 9-element numpy array of 0s and 1s. Each element of the array represents the presence of a
        particular bit value (element 0 = bit 0, element 1 = bit 1, ..., element 3 = bit 4 and so on...)
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

def extractMatchedLines(col2get, refData, compData, refLines, compLines, bitmask=[]):
    """Extracts only matching lines of data for a specific column of refData and compData. Returns empty list
    if the specified column is not found in both tables.

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

    bitmask : numpy.ndarray, optional
        list of True/False values where False corresponds to values to keep, and True corresponds to values
        to remove

    Returns
    -------
    return_ra : numpy ndarray
        A 2 x len(refLines) sized numpy array. Column 1: matched reference values. Column 2: The
        corresponding matched comparison values
    """
    if col2get in list(refData.keys()) and col2get in list(compData.keys()):
        matching_refData = refData[col2get][refLines].data
        matching_compData = compData[col2get][compLines].data
        if bitmask != []:
            bitmask = bitmask.astype(int)
            matching_refData = np.ma.array(matching_refData, mask=bitmask)
            matching_compData = np.ma.array(matching_compData, mask=bitmask)
            matching_refData = matching_refData.compressed()
            matching_compData = matching_compData.compressed()
        return_ra = np.stack((matching_refData, matching_compData))
    else:
        return_ra = np.empty((2, 0))
    return return_ra


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
        The desired level of verboseness in the log statements displayed on the screen and written to the .log
        file.

    Returns
    -------
    matching_lines_ref : list
        A list of the indices of reference sourcelist sources that match comparison sourcelist sources

    matching_lines_img : list
        A corresponding list of the indices of comparison sourcelist sources that match reference sourcelist
        sources
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
        log.info("WARNING: Unable to fetch values for xref and yref from fits file headers. Using xref = 0.0 "
                 "and yref = 0.0.")
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
        "Reference sourcelist:  {} of {} total sources matched ({} %)".format(len(matching_lines_ref),
                                                                              slLengths[0],
                                                                              100.0 * (float(
                                                                                  len(matching_lines_ref)) /
                                                                                       float(slLengths[0]))))
    log.info(
        "Comparison sourcelist: {} of {} total sources matched ({} %)".format(len(matching_lines_img),
                                                                              slLengths[1],
                                                                              100.0 * (float(
                                                                                  len(matching_lines_img)) /
                                                                                       float(slLengths[1]))))

    return (matching_lines_ref, matching_lines_img)


# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-

def make_flag_mask(matched_flag_values, good_flag_sum, missing_mask):
    """Returns a list of the array index values to mask based on user-specified good flag value, and missing
    mask

    Parameters
    ----------
    matched_flag_values : numpy.ndarray
        A 2 x len(refLines) sized numpy array. Column 1: matched reference values.
        Column 2: The corresponding matched comparison values

    good_flag_sum : int
        sum of flag bit values that should be considered "good" for masking purposes

    missing_mask : numpy.ndarray
        Updated list of True/False values where False corresponds to values to keep, and True corresponds to
        values to remove

    Returns
    -------
    masked_index_list : numpy list
        list of the array index values to mask
    """
    full_refFlag_list = []
    full_compFlag_list = []
    bitmask = missing_mask  # np.full(len(matched_flag_values[0]),0,dtype=bool)
    if good_flag_sum != 255:
        good_bit_list = deconstruct_flag(good_flag_sum)  # break good bit sum into list of component bits
        good_bit_list[0] = 1
        bad_bit_list = np.invert(good_bit_list.astype(bool))  # invert good bit list to make bad bit list
    ctr = 0
    for refFlagVal, compFlagVal in zip(matched_flag_values[0], matched_flag_values[1]):
        refFlag_list = deconstruct_flag(refFlagVal)  # break ref flag bit sum into list of component bits
        full_refFlag_list.append(refFlag_list)
        compFlag_list = deconstruct_flag(compFlagVal)  # break comp flag bit sum into list of component bits
        full_compFlag_list.append(compFlag_list)
        if good_flag_sum != 255:
            merged_flag_val = np.logical_or(refFlag_list, compFlag_list)  # merge comp and ref flag lists
            bitmask[ctr] = np.any(np.logical_and(merged_flag_val, bad_bit_list))  # generate mask value by checking to see if any of the bad bits are found in the merged comp+ref bit list

        ctr += 1

    masked_index_list = np.where(bitmask == True)
    log.info("{} of {} ({} %) values masked.".format(np.shape(masked_index_list)[1], ctr,
                                                     100.0*(float(np.shape(masked_index_list)[1])/float(ctr))))
    log.info("{} remain.".format(ctr-np.shape(masked_index_list)[1]))
    return bitmask


# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-

def mask_missing_values(refData, compData, refLines, compLines, columns_to_compare):
    """Update the bitmask to include lines where any values are missing, nan, or inf from any column in
    either the comp or ref matched catalogs

    Parameters
    ----------
    refData : astropy Table object
        reference data table

    compData : astropy Table object
        comparison data table

    refLines : numpy.ndarray
        List of matching refData line numbers

    compLines : numpy.ndarray
        List of matching compData line numbers

    columns_to_compare : list
        list of columns that are common to both comparison and reference catalogs and will be used in the
        comparisons

    Returns
    -------
    out_mask : numpy.ndarray, optional
        Updated list of True/False values where False corresponds to values to keep, and True corresponds to
        values to remove
    """
    out_mask = np.full(np.shape(refLines), 0, dtype=bool)

    for col_title in columns_to_compare:
        matching_refData = refData[col_title][refLines].data
        matching_compData = compData[col_title][compLines].data
        for data_set in [matching_refData, matching_compData]:  # merge together all input mask arrays
            if hasattr(data_set, "mask"):
                out_mask = np.logical_or(out_mask, data_set.mask)
            inf_nan_idx = np.where((np.isnan(data_set) == True) | (np.isinf(data_set) == True))  # identify any 'nan' or 'inf' values and flag them out as well
            for mask_idx in inf_nan_idx:
                out_mask[mask_idx] = True
    return out_mask


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
        'X position' and 'Y position' data columns will always be returned. However, not every sourcelist or
        catalog file is guaranteed to have any of the seven remaining data columns in the above list.

    Parameters
    ----------
    slNames : list
        A list containing the reference sourcelist filename and the comparison sourcelist filename, in that
        order.

    Returns
    -------
    refData : astropy.table object
        data from the reference sourcelist

    compData : astropy.table object
        data from the comparison sourcelist
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
                           "FLUX2": "FluxAp2", "MAGNITUDE1": "MagAp1", "MAGNITUDE2": "MagAp2",
                           "MERR1": "MagErrAp1", "MERR2": "MagErrAp2", "MSKY": "MSkyAp2",
                           "STDEV": "StdevAp2", "FLAGS": "Flags", "ID": "ID",
                           "CI": "CI"}
    titleSwapDict_segment = {"X": "X-Centroid", "Y": "Y-Centroid", "RA": "RA", "DEC": "DEC", "FLUX1": "n/a",
                             "FLUX2": "FluxAp2", "MAGNITUDE1": "MagAp1", "MAGNITUDE2": "MagAp2",
                             "MERR1": "MagErrAp1", "MERR2": "MagErrAp2", "MSKY": "n/a", "STDEV": "n/a",
                             "FLAGS": "Flags", "ID": "ID", "CI": "CI"}
    titleSwapDict_daoTemp = {"X": "XCENTER", "Y": "YCENTER", "RA": "n/a", "DEC": "n/a", "FLUX1": "FLUX1",
                             "FLUX2": "FLUX2", "MAGNITUDE1": "MAG1", "MAGNITUDE2": "MAG2", "FLAGS": "n/a",
                             "ID": "ID", "MERR1": "MERR1", "MERR2": "MERR2", "MSKY": "MSKY", "STDEV": "STDEV"}
    titleSwapDict_sourceX = {"X": "X_IMAGE", "Y": "Y_IMAGE", "RA": "RA", "DEC": "DEC", "FLUX1": "FLUX_APER1",
                             "FLUX2": "FLUX_APER2", "MAGNITUDE1": "MAG_APER1", "MAGNITUDE2": "MAG_APER2",
                             "FLAGS": "FLAGS", "ID": "NUMBER"}
    titleSwapDict_cooNew = {"X": "col1", "Y": "col2", "RA": "n/a", "DEC": "n/a", "FLUX1": "n/a",
                            "FLUX2": "n/a", "MAGNITUDE1": "n/a", "MAGNITUDE2": "n/a", "FLAGS": "n/a",
                            "ID": "col7"}
    # titleSwapDict_cooOld = {"X": "XCENTER", "Y": "YCENTER", "RA": "n/a", "DEC": "n/a", "FLUX1": "n/a",
    #                         "FLUX2": "n/a", "MAGNITUDE1": "n/a", "MAGNITUDE2": "n/a", "FLAGS": "n/a",
    #                         "ID": "ID"}

    titleSwapDict_cooOld2 = {"X": "XCENTER", "Y": "YCENTER", "RA": "RA", "DEC": "DEC", "FLUX1": "FLUX_0.05",
                             "FLUX2": "FLUX_0.15", "MAGNITUDE1": "MAG_0.05", "MAGNITUDE2": "MAG_0.15",
                             "FLAGS": "n/a", "ID": "ID", "MERR1": "MERR_0.05", "MERR2": "MERR_0.15",
                             "MSKY": "MSKY", "STDEV": "STDEV"}

    titleSwapDict_daorep = {"X": "X", "Y": "Y", "RA": "n/a", "DEC": "n/a", "FLUX1": "flux_0",
                            "FLUX2": "flux_1", "MAGNITUDE1": "mag_0", "MAGNITUDE2": "mag_1", "FLAGS": "n/a",
                            "ID": "n/a"}
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
