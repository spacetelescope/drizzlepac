#!/usr/bin/env python

"""This script is a modernized implementation of tweakreg.

"""

import datetime
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u
from collections import OrderedDict
from drizzlepac import updatehdr
import glob
import math
import numpy as np
import os
import pdb
import pickle
from stsci.tools import fileutil
from stwcs.wcsutil import headerlet
from stwcs.wcsutil import HSTWCS
import sys
import tweakwcs
try:
    from drizzlepac.hlautils import astrometric_utils as amutils
    from drizzlepac.hlautils import astroquery_utils as aqutils
    from drizzlepac.hlautils import analyze as filter
    from drizzlepac.hlautils import get_git_rev_info
except:
    from drizzlepac.hlautils import astrometric_utils as amutils
    from drizzlepac.hlautils import astroquery_utils as aqutils
    from drizzlepac.hlautils import analyze as filter
    from drizzlepac.hlautils import get_git_rev_info

MIN_CATALOG_THRESHOLD = 3
MIN_OBSERVABLE_THRESHOLD = 10
MIN_CROSS_MATCHES = 3
MIN_FIT_MATCHES = 6
MAX_FIT_RMS = 10 # RMS now in mas, 1.0
MAX_FIT_LIMIT = 1000 # Maximum RMS that a result is useful
MAX_SOURCES_PER_CHIP = 250  # Maximum number of sources per chip to include in source catalog

# Module-level dictionary contains instrument/detector-specific parameters used later on in the script.
detector_specific_params = {"acs":
                                {"hrc":
                                     {"fwhmpsf": 0.073,
                                      "classify": True,
                                      "threshold": None},
                                 "sbc":
                                     {"fwhmpsf": 0.065,
                                      "classify": False,
                                      "threshold": 2.0},
                                 "wfc":
                                     {"fwhmpsf": 0.13, #0.076,
                                      "classify": True,
                                      "threshold": -1.1}},
                            "wfc3":
                                {"ir":
                                     {"fwhmpsf": 0.14,
                                      "classify": False,
                                      "threshold": None},
                                 "uvis":
                                     {"fwhmpsf": 0.076,
                                      "classify": True,
                                      "threshold": None}}} # fwhmpsf in units of arcsec


# ----------------------------------------------------------------------------------------------------------------------


def check_and_get_data(input_list,**pars):
    """Verify that all specified files are present. If not, retrieve them from MAST.

    Parameters
    ----------
    imglist : list
        List of one or more calibrated fits images that will be used for catalog generation.

    Returns
    =======
    input_file_list : list
        list of full filenames

    """
    totalInputList = []
    retrieve_list = [] # list of already retrieved ASNs, only request each once
    for input_item in input_list:
        filelist = aqutils.retrieve_observation(input_item,**pars)
        if len(filelist) == 0:
            # look for local copy of the file
            input_rootname = input_item[:input_item.find('_')]
            fitsfilenames = sorted(glob.glob("{}_fl?.fits".format(input_rootname)))
            if len(fitsfilenames) > 0:
                imghdu = fits.open(fitsfilenames[0])
                imgprimaryheader = imghdu[0].header
                try:
                    asnid = imgprimaryheader['ASN_ID'].strip().lower()
                    if asnid == 'NONE':
                        asnid = None
                except KeyError:
                    asnid = None
                if asnid and asnid not in retrieve_list:
                    filelist = aqutils.retrieve_observation(asnid,**pars)
                    retrieve_list.append(asnid) # record asn's already retrieved

        if len(filelist) > 0:
            totalInputList += filelist

    if len(filelist) > 0: totalInputList = sorted(
        list(set(totalInputList)))  # remove duplicate list elements, sort resulting list of unique elements
    print("TOTAL INPUT LIST: ",totalInputList)
    # TODO: add trap to deal with non-existent (incorrect) rootnames
    # TODO: Address issue about how the code will retrieve association information if there isn't a local file to get 'ASN_ID' header info
    return(totalInputList)


# ----------------------------------------------------------------------------------------------------------------------


def convert_string_tf_to_boolean(invalue):
    """Converts string 'True' or 'False' value to Boolean True or Boolean False.

    :param invalue: string
        input true/false value

    :return: Boolean
        converted True/False value
    """
    outvalue = False
    if invalue == 'True':
        outvalue = True
    return(outvalue)


# ----------------------------------------------------------------------------------------------------------------------


def perform_align(input_list, archive=False, clobber=False, debug=True, update_hdr_wcs=False,
                  print_fit_parameters=True, print_git_info=False, output=True): # TODO: set 'debug' and 'output' back to 'False' before release.
    """Main calling function.

    Parameters
    ----------
    input_list : list
        List of one or more IPPSSOOTs (rootnames) to align.

    archive : Boolean
        Retain copies of the downloaded files in the astroquery created sub-directories?

    clobber : Boolean
        Download and overwrite existing local copies of input files?

    debug : Boolean
        Attempt to use saved sourcelists stored in pickle files if they exist, or if they do not exist, save
        sourcelists in pickle files for reuse so that step 4 can be skipped for faster subsequent debug/development
        runs??

    update_hdr_wcs : Boolean
        Write newly computed WCS information to image image headers?

    print_fit_parameters : Boolean
        Specify whether or not to print out FIT results for each chip.

    print_git_info : Boolean
        Display git repository information?

    output : Boolean
        Should utils.astrometric_utils.create_astrometric_catalog() generate file 'ref_cat.ecsv' and should
        generate_source_catalogs() generate the .reg region files for every chip of every input image and should
        generate_astrometric_catalog() generate file 'refcatalog.cat'?

    Returns
    -------
    int value 0 if successful, int value 1 if unsuccessful

    """

    # Define astrometric catalog list in priority order
    catalogList = ['GAIADR2', 'GSC241']

    # 0: print git info
    if print_git_info:
        print("-------------------- STEP 0: Display Git revision info  --------------------")
        full_path = os.path.dirname(__file__)+"/utils"
        repo_path=None
        if "hlapipeline/hlapipeline" in full_path:
            repo_path = full_path.split("hlapipeline/hlapipeline")[0]+"hlapipeline"
        elif "hlapipeline" in full_path:
            repo_path = full_path.split("hlapipeline")[0]+"hlapipeline"
        else:
            pass
        if not os.path.exists(repo_path): repo_path = None # protect against non-existent paths
        if repo_path:
            get_git_rev_info.print_rev_id(repo_path) # Display git repository information
        else:
            print("WARNING: Unable to display Git repository revision information.")

    # 1: Interpret input data and optional parameters
    print("-------------------- STEP 1: Get data --------------------")
    zeroDT = startingDT = datetime.datetime.now()
    print(str(startingDT))
    imglist = check_and_get_data(input_list, archive=archive, clobber=clobber)
    print("\nSUCCESS")

    currentDT = datetime.datetime.now()
    deltaDT = (currentDT - startingDT).total_seconds()
    print('Processing time of [STEP 1]: {} sec'.format(deltaDT))
    startingDT = currentDT
    # 2: Apply filter to input observations to insure that they meet minimum criteria for being able to be aligned
    print("-------------------- STEP 2: Filter data --------------------")
    filteredTable = filter.analyze_data(imglist)

    # Check the table to determine if there is any viable data to be aligned.  The
    # 'doProcess' column (bool) indicates the image/file should or should not be used
    # for alignment purposes.  For filtered data, 'doProcess=0' and 'status=9999' in the table
    # (the status value by default), so there is no need to update the filteredTable here.
    if filteredTable['doProcess'].sum() == 0:
        print("No viable images in filtered table - no processing done.\n")
        return(filteredTable)

    # Get the list of all "good" files to use for the alignment
    processList = filteredTable['imageName'][np.where(filteredTable['doProcess'])]
    processList = list(processList) #Convert processList from numpy list to regular python list
    print("\nSUCCESS")

    currentDT = datetime.datetime.now()
    deltaDT = (currentDT - startingDT).total_seconds()
    print('Processing time of [STEP 2]: {} sec'.format(deltaDT))
    startingDT = currentDT
    # 3: Build WCS for full set of input observations
    print("-------------------- STEP 3: Build WCS --------------------")
    refwcs = amutils.build_reference_wcs(processList)
    print("\nSUCCESS")


    currentDT = datetime.datetime.now()
    deltaDT = (currentDT - startingDT).total_seconds()
    print('Processing time of [STEP 3]: {} sec'.format(deltaDT))
    startingDT = currentDT
    # 4: Extract catalog of observable sources from each input image
    print("-------------------- STEP 4: Source finding --------------------")
    if debug:
        pickle_filename = "{}.source_catalog.pickle".format(
            processList[0])  # TODO: All this pickle stuff is only here for debugging. <START>
        if os.path.exists(pickle_filename):
            pickle_in = open(pickle_filename, "rb")
            extracted_sources = pickle.load(pickle_in)
            print("Using sourcelist extracted from {} generated during the last run to save time.".format(
                pickle_filename))
        else:
            extracted_sources = generate_source_catalogs(processList,
                                                         centering_mode='starfind',nlargest=MAX_SOURCES_PER_CHIP,output=output)
            pickle_out = open(pickle_filename, "wb")
            pickle.dump(extracted_sources, pickle_out)
            pickle_out.close()
            print("Wrote ", pickle_filename)
    else:
        extracted_sources = generate_source_catalogs(processList,
                                                     centering_mode='starfind',nlargest=MAX_SOURCES_PER_CHIP,output=output)

    for imgname in extracted_sources.keys():
        table=extracted_sources[imgname]["catalog_table"]

        # Get the location of the current image in the filtered table
        index = np.where(filteredTable['imageName']==imgname)[0][0]

        # First ensure sources were found
        if table[1] == None:
            print("No sources found in image {}".format(imgname))
            filteredTable[:]['status'] = 1
            filteredTable[:]['processMsg'] = "No sources found"
            return(filteredTable)

        # The catalog of observable sources must have at least MIN_OBSERVABLE_THRESHOLD entries to be useful
        total_num_sources = 0
        for chipnum in table.keys():
            total_num_sources += len(table[chipnum])

        # Update filtered table with number of found sources
        filteredTable[index]['foundSources'] = total_num_sources

        if total_num_sources < MIN_OBSERVABLE_THRESHOLD:
            print("Not enough sources ({}) found in image {}".format(total_num_sources,imgname))
            filteredTable[:]['status'] = 1
            filteredTable[:]['processMsg'] = "Not enough sources found"
            return(filteredTable)
    print("\nSUCCESS")
    currentDT = datetime.datetime.now()
    deltaDT = (currentDT - startingDT).total_seconds()
    print('Processing time of [STEP 4]: {} sec'.format(deltaDT))
    startingDT = currentDT
    # 5: Retrieve list of astrometric sources from database

    # Convert input images to tweakwcs-compatible NDData objects and
    # attach source catalogs to them.
    imglist = []
    for group_id, image in enumerate(processList):
        img = amutils.build_nddata(image, group_id,
                                   extracted_sources[image]['catalog_table'])
        # add the name of the image to the imglist object
        for im in img:
            im.meta['name'] = image
        imglist.extend(img)

    best_fit_rms = -99999.0
    best_fitStatusDict={}
    best_fitQual = 5
    fit_algorithm_list= [match_2dhist_fit,match_default_fit]
    for catalogIndex in range(0, len(catalogList)): #loop over astrometric catalog
        print("-------------------- STEP 5: Detect astrometric sources --------------------")
        print("Astrometric Catalog: ",catalogList[catalogIndex])
        reference_catalog = generate_astrometric_catalog(processList, catalog=catalogList[catalogIndex], output=output)

        currentDT = datetime.datetime.now()
        deltaDT = (currentDT - startingDT).total_seconds()
        print('Processing time of [STEP 5]: {} sec'.format(deltaDT))
        startingDT = currentDT

        if len(reference_catalog) < MIN_CATALOG_THRESHOLD:
            print("Not enough sources found in catalog " + catalogList[catalogIndex])
            fitQual = 5
            if catalogIndex < len(catalogList) -1:
                print("Try again with other catalog")
            else:
                print("ERROR! No astrometric sources found in any catalog. Exiting...") #bail out if not enough sources can be found any of the astrometric catalogs
                filteredTable['status'][:] = 1
                filteredTable['processMsg'][:] = "No astrometric sources found"
                return (filteredTable)
        else:
            print("-------------------- STEP 5b: Cross matching and fitting --------------------")
            for algorithm_name in fit_algorithm_list: #loop over fit algorithm type
                print("------------------ Catalog {} matched using {} ------------------ ".format(catalogList[catalogIndex],algorithm_name.__name__))
                try:
                    #execute the correct fitting/matching algorithm
                    imglist = algorithm_name(imglist, reference_catalog)

                    # determine the quality of the fit
                    fit_rms, fit_num, fitQual, filteredTable, fitStatusDict = determine_fit_quality(imglist,filteredTable, print_fit_parameters=print_fit_parameters)

                    # Figure out which fit solution to go with based on fitQual value and maybe also total_rms
                    if fitQual < 5:
                        if fitQual == 1: #valid, non-comprimised solution with total rms < 10 mas...go with this solution.
                            best_fit_rms = fit_rms
                            best_fit_num = fit_num
                            for item in imglist:
                                item.best_meta = item.meta.copy()
                            best_fitStatusDict = fitStatusDict.copy()
                            break #break out of while loop
                        elif fitQual < best_fitQual: # better solution found. keep looping but with the better solution as "best" for now.
                            print("Better solution found!")
                            best_fit_rms = fit_rms
                            best_fit_num = fit_num
                            for item in imglist:
                                item.best_meta = item.meta.copy()
                            best_fitStatusDict = fitStatusDict.copy()
                            best_fitQual = fitQual
                        elif fitQual == best_fitQual: # new solution same level of fitQual. Choose whichever one has the lowest total rms as "best" and keep looping.
                            if best_fit_rms >= 0.:
                                if fit_rms < best_fit_rms:
                                    best_fit_rms = fit_rms
                                    best_fit_num = fit_num
                                    for item in imglist:
                                        item.best_meta = item.meta.copy()
                                    best_fitStatusDict = fitStatusDict.copy()
                        else: # new solution has worse fitQual. discard and continue looping.
                            continue
                except Exception:
                    print("WARNING: Catastrophic fitting failure with catalog {} and matching algorithm {}.".format(catalogList[catalogIndex],algorithm_name.__name__))
                    filteredTable['status'][:] = 1
                    filteredTable['processMsg'][:] = "Fitting failure"
                    # It may be there are additional catalogs and algorithms to try, so keep going
                    continue

                if fitQual == 1:  # break out of inner fit algorithm loop
                    break
        if fitQual == 1: #break out of outer astrometric catalog loop
            break
    currentDT = datetime.datetime.now()
    deltaDT = (currentDT - startingDT).total_seconds()
    print('Processing time of [STEP 5b]: {} sec'.format(deltaDT))
    startingDT = currentDT
    # 6: Populate the filteredTable
    print("-------------------- STEP 6: Collect up information and populate the filtered table --------------------")
    if best_fit_rms > 0 and best_fit_rms < MAX_FIT_RMS:
        print("The fitting process was successful with a best fit total rms of {} mas".format(best_fit_rms))
    else:
        print("The fitting process was unsuccessful with a best fit total rms of {} mas".format(best_fit_rms))

    if best_fit_rms > 0 and best_fit_rms < MAX_FIT_LIMIT:
        # update to the meta information with the lowest rms if it is reasonable
        for item in imglist:
            item.meta = item.best_meta.copy()
        filteredTable['status'][:] = 0
        fitStatusDict = best_fitStatusDict.copy()

        # Protect the writing of the table within the best_fit_rms
        info_keys = OrderedDict(imglist[0].meta['tweakwcs_info']).keys()
        # Update filtered table with number of matched sources and other information
        for item in imglist:
            imgname = item.meta['name']
            index = np.where(filteredTable['imageName'] == imgname)[0][0]

            if item.meta['tweakwcs_info']['status'].startswith("FAILED") != True:
                for tweakwcs_info_key in info_keys:
                    if not tweakwcs_info_key.startswith("matched"):
                        if tweakwcs_info_key.lower() == 'rms':
                            filteredTable[index]['rms_x'] = item.meta['tweakwcs_info'][tweakwcs_info_key][0]
                            filteredTable[index]['rms_y'] = item.meta['tweakwcs_info'][tweakwcs_info_key][1]

                filteredTable[index]['catalog'] = item.meta['tweakwcs_info']['catalog']
                filteredTable[index]['catalogSources'] = len(reference_catalog)
                filteredTable[index]['matchSources'] = item.meta['tweakwcs_info']['nmatches']
                filteredTable[index]['rms_ra'] = item.meta['tweakwcs_info']['RMS_RA'].value
                filteredTable[index]['rms_dec'] = item.meta['tweakwcs_info']['RMS_DEC'].value
                filteredTable[index]['fit_rms'] = item.meta['tweakwcs_info']['FIT_RMS']
                filteredTable[index]['total_rms'] = item.meta['tweakwcs_info']['TOTAL_RMS']
                # populate filteredTable fields "status", "compromised" and
                # "processMsg" with fitStatusDict fields "valid", "compromised"
                # and "reason".
                explicitDictKey ="{},{}".format(item.meta['name'], item.meta['chip'])
                if fitStatusDict[explicitDictKey]['valid'] == True:
                    filteredTable[index]['status'] = 0
                else:
                    filteredTable[index]['status'] = 1
                if fitStatusDict[explicitDictKey]['compromised'] == False:
                    filteredTable['compromised'] = 0
                else:
                    filteredTable['compromised'] = 1
                if fitStatusDict[explicitDictKey]['reason'] != "":
                    filteredTable[index]['processMsg'] = fitStatusDict[explicitDictKey]['reason']

    currentDT = datetime.datetime.now()
    deltaDT = (currentDT - startingDT).total_seconds()
    print('Processing time of [STEP 6]: {} sec'.format(deltaDT))
    startingDT = currentDT
    # 7: Write new fit solution to input image headers
    print("-------------------- STEP 7: Update image headers with new WCS information --------------------")
    if best_fit_rms > 0 and update_hdr_wcs:
        headerlet_dict = update_image_wcs_info(imglist)
        for tableIndex in range(0,len(filteredTable)):
            filteredTable[tableIndex]['headerletFile'] = headerlet_dict[filteredTable[tableIndex]['imageName']]
        print("\nSUCCESS")
    else:
        print("\n STEP SKIPPED")

    currentDT = datetime.datetime.now()
    deltaDT = (currentDT - startingDT).total_seconds()
    print('Processing time of [STEP 7]: {} sec'.format(deltaDT))
    print('TOTAL Processing time of {} sec'.format((currentDT- zeroDT).total_seconds()))
    print(best_fitStatusDict)
    return (filteredTable)


# ----------------------------------------------------------------------------------------------------------------------


def match_default_fit(imglist, reference_catalog):
    """Perform cross-matching and final fit using 2dHistogram matching

    Parameters
    ----------
    imglist : list
        List of input image NDData objects with metadata and source catalogs

    reference_catalog : Table
        Astropy Table of reference sources for this field

    Returns
    --------
    imglist : list
        List of input image NDData objects with metadata and source catalogs

    """
    # Specify matching algorithm to use
    match = tweakwcs.TPMatch(searchrad=250, separation=0.1,
                             tolerance=100, use2dhist=False)
    # Align images and correct WCS
    tweakwcs.tweak_image_wcs(imglist, reference_catalog, match=match, expand_refcat=False) #TODO: turn on 'expand_refcat' option in future development

    # Interpret RMS values from tweakwcs
    interpret_fit_rms(imglist, reference_catalog)

    return imglist


# ----------------------------------------------------------------------------------------------------------------------


def match_2dhist_fit(imglist, reference_catalog):
    """Perform cross-matching and final fit using 2dHistogram matching

    Parameters
    ----------
    imglist : list
        List of input image NDData objects with metadata and source catalogs

    reference_catalog : Table
        Astropy Table of reference sources for this field

    Returns
    --------
    imglist : list
        List of input image NDData objects with metadata and source catalogs

    """
    print("-------------------- STEP 5b: (match_2dhist_fit) Cross matching and fitting --------------------")
    # Specify matching algorithm to use
    match = tweakwcs.TPMatch(searchrad=75, separation=0.1,
                             tolerance=2.0, use2dhist=True)
    # Align images and correct WCS
    tweakwcs.tweak_image_wcs(imglist, reference_catalog, match=match, expand_refcat=False) #TODO: turn on 'expand_refcat' option in future development

    # Interpret RMS values from tweakwcs
    interpret_fit_rms(imglist, reference_catalog)

    return imglist


# ----------------------------------------------------------------------------------------------------------------------


def determine_fit_quality(imglist,filteredTable, print_fit_parameters=True):
    """Determine the quality of the fit to the data

    Parameters
    ----------
    imglist : list
        output of interpret_fits. Contains sourcelist tables, newly computed WCS info, etc. for every chip of every valid
        input image.  This list should have been  updated, in-place, with the new RMS values;
        specifically,

            * 'FIT_RMS': RMS of the separations between fitted image positions and reference positions
            * 'TOTAL_RMS': mean of the FIT_RMS values for all observations
            * 'NUM_FITS': number of images/group_id's with successful fits included in the TOTAL_RMS

        These entries are added to the 'tweakwcs_info' dictionary.

    filteredTable : object
        Astropy Table object containing data pertaining to the associated dataset, including
        the doProcess bool.  It is intended this table is updated by subsequent functions for
        bookkeeping purposes.

    print_fit_parameters : bool
        Specify whether or not to print out FIT results for each chip

    Returns
    -------
    max_rms_val : float
        The best Total rms determined from all of the images

    num_xmatches: int
        The number of stars used in matching the data


    fitQual : int
        fit quality catagory:
        1 = valid solution with rms < 10 mas;
        2 = Valid but compromised solution with rms < 10 mas;
        3 = Valid solution with RMS >= 10 mas;
        4 = Valid but compromised solution with RMS >= 10 mas;
        5 = Not valid solution

    filteredTable : object
        modified filteredTable objecgt

    fitStatusDict : dictionary
        Dictionary containing the following:
            overall fit validity (Boolean)
            total (visit-level) RMS value in mas (float)
            number of matched sources (int)
            fit compromised status (Boolean)
            reason fit is considered 'compromised' (only populated if 'compromised' field is "True")
    """
    tweakwcs_info_keys = OrderedDict(imglist[0].meta['tweakwcs_info']).keys()
    max_rms_val = 1e9
    num_xmatches = 0
    fitStatusDict={}
    xshifts=[]
    yshifts=[]
    overall_valid = True
    overall_comp = False
    for item in imglist:
        xshifts.append(item.meta['tweakwcs_info']['shift'][0])
        yshifts.append(item.meta['tweakwcs_info']['shift'][1])

    for item in imglist:
        image_name = item.meta['name']
        chip_num = item.meta['chip']
        #Handle fitting failures (no matches found)
        if item.meta['tweakwcs_info']['status'].startswith("FAILED") == True:
                print("No cross matches found in any catalog for {} - no processing done.".format(image_name))
                continue
        fit_rms_val = item.meta['tweakwcs_info']['FIT_RMS']
        max_rms_val = item.meta['tweakwcs_info']['TOTAL_RMS']
        num_xmatches = item.meta['tweakwcs_info']['nmatches']
        if num_xmatches < MIN_CROSS_MATCHES:
            if catalogIndex < numCatalogs-1:
                print("Not enough cross matches found between astrometric catalog and sources found in {}".format(image_name))
                continue

        # Build fitStatusDict entry
        dictKey = "{},{}".format(image_name, chip_num)
        fitStatusDict[dictKey] = {'valid': False,
                                  'max_rms': max_rms_val,
                                  'num_matches': num_xmatches,
                                  'compromised': False,
                                  'reason': ""} # Initialize dictionary entry for current image/chip
        # Execute checks
        nmatchesCheck = False
        if num_xmatches > 4:
            nmatchesCheck = True

        radialOffsetCheck = False
        radialOffset = math.sqrt(float(item.meta['tweakwcs_info']['shift'][0])**2 +
                                 float(item.meta['tweakwcs_info']['shift'][0])**2)*item.wcs.pscale #radial offset in arssec
        if float(num_xmatches) * 0.36 > 0.8 + (radialOffset/10.0)**8:
            radialOffsetCheck = True

        largeRmsCheck = True
        if fit_rms_val > 150. or max_rms_val > 150.:
            largeRmsCheck = False

        # fitRmsCheck = False
        # if fit_rms_val < max_rms_val:
        #     fitRmsCheck = True

        consistencyCheck = True
        rms_limit = max(item.meta['tweakwcs_info']['TOTAL_RMS'], 10.)
        if not math.sqrt(np.std(np.asarray(xshifts)) ** 2 + np.std(np.asarray(yshifts)) ** 2) <= (
                    rms_limit / 1000.0) / (item.wcs.pscale):
            consistencyCheck = False

        # Decide if fit solutions are valid based on checks
        if consistencyCheck == False: # Failed consistency check
            fitStatusDict[dictKey]['valid'] = False
            fitStatusDict[dictKey]['compromised'] = False
            fitStatusDict[dictKey]['reason'] = "Consistency violation!"
        elif largeRmsCheck == False: # RMS value(s) too large
            fitStatusDict[dictKey]['valid'] = False
            fitStatusDict[dictKey]['compromised'] = False
            fitStatusDict[dictKey]['reason'] = "RMS too large (>150 mas)!"
        elif radialOffsetCheck == False: # Failed radial offset check
            fitStatusDict[dictKey]['valid'] = False
            fitStatusDict[dictKey]['compromised'] = False
            fitStatusDict[dictKey]['reason'] = "Radial offset value too large!"
        elif nmatchesCheck == False: # Too few matches
            fitStatusDict[dictKey]['valid'] = True
            fitStatusDict[dictKey]['compromised'] = True
            fitStatusDict[dictKey]['reason'] = "Too few matches!"
        else: # all checks passed. Valid solution.
            fitStatusDict[dictKey]['valid'] = True
            fitStatusDict[dictKey]['compromised'] = False
            fitStatusDict[dictKey]['reason'] = ""

        # for now, generate overall valid and compromised values. Basically, if any of the entries for "valid" is False,
        # treat the whole dataset as not valid. Same goes for compromised.
        if fitStatusDict[dictKey]['valid'] == False:
            overall_valid = False
        if fitStatusDict[dictKey]['compromised'] == True:
            overall_comp = True

        print('RESULTS FOR {} Chip {}: FIT_RMS = {} mas, TOTAL_RMS = {} mas, NUM =  {}'.format(image_name, item.meta['chip'], fit_rms_val, max_rms_val, num_xmatches))
        # print fit params to screen
        if print_fit_parameters:
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIT PARAMETERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print("image: {}".format(image_name))
            print("chip: {}".format(item.meta['chip']))
            print("group_id: {}".format(item.meta['group_id']))
            for tweakwcs_info_key in tweakwcs_info_keys:
                if not tweakwcs_info_key.startswith("matched"):
                    print("{} : {}".format(tweakwcs_info_key,item.meta['tweakwcs_info'][tweakwcs_info_key]))
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print("nmatchesCheck: {} radialOffsetCheck: {} largeRmsCheck: {}, consistencyCheck: {}".format(nmatchesCheck,radialOffsetCheck,largeRmsCheck,consistencyCheck))



    # determine which fit quality category this latest fit falls into
    if overall_valid == False:
        fitQual = 5
        print("FIT SOLUTION REJECTED")
        filteredTable['status'][:] = 1
        for ctr in range(0, len(filteredTable)):
            filteredTable[ctr]['processMsg'] = fitStatusDict[filteredTable[ctr]['imageName'] + ",1"]["reason"]
    else:
        if overall_comp == False and max_rms_val < 10.:
            print("Valid solution with RMS < 10 mas found!")
            fitQual = 1
        elif overall_comp == True and max_rms_val < 10.:
            print("Valid but compromised solution with RMS < 10 mas found!")
            fitQual = 2
        elif overall_comp == False and max_rms_val >= 10.:
            print("Valid solution with RMS >= 10 mas found!")
            fitQual = 3
        else:
            print("Valid but compromised solution with RMS >= 10 mas found!")
            fitQual = 4

    if print_fit_parameters:
        for item in imglist: print(fitStatusDict["{},{}".format(item.meta['name'], item.meta['chip'])])

    if max_rms_val > MAX_FIT_RMS:
        print("Total fit RMS value = {} mas greater than the maximum threshold value {}.".format(max_rms_val, MAX_FIT_RMS))
    if not overall_valid:
        print("The fit solution for some or all of the images is not valid.")
    if max_rms_val > MAX_FIT_RMS or overall_valid == False:
        print("Try again with the next catalog")
    else:
        print("Fit calculations successful.")
    return max_rms_val, num_xmatches, fitQual, filteredTable, fitStatusDict


# ----------------------------------------------------------------------------------------------------------------------


def generate_astrometric_catalog(imglist, **pars):
    """Generates a catalog of all sources from an existing astrometric catalog are in or near the FOVs of the images in
        the input list.

    Parameters
    ----------
    imglist : list
        List of one or more calibrated fits images that will be used for catalog generation.

    Returns
    =======
    ref_table : object
        Astropy Table object of the catalog

    """
    # generate catalog
    temp_pars = pars.copy()
    if pars['output'] == True:
        pars['output'] = 'ref_cat.ecsv'
    else:
        pars['output'] = None
    out_catalog = amutils.create_astrometric_catalog(imglist,**pars)
    pars = temp_pars.copy()
    #if the catalog has contents, write the catalog to ascii text file
    if len(out_catalog) > 0 and pars['output']:
        catalog_filename = "refcatalog.cat"
        out_catalog.write(catalog_filename, format="ascii.fast_commented_header")
        print("Wrote reference catalog {}.".format(catalog_filename))

    return(out_catalog)


# ----------------------------------------------------------------------------------------------------------------------


def generate_source_catalogs(imglist, **pars):
    """Generates a dictionary of source catalogs keyed by image name.

    Parameters
    ----------
    imglist : list
        List of one or more calibrated fits images that will be used for source detection.

    Returns
    -------
    sourcecatalogdict : dictionary
        a dictionary (keyed by image name) of two element dictionaries which in tern contain 1) a dictionary of the
        detector-specific processing parameters and 2) an astropy table of position and photometry information of all
        detected sources
    """
    output = pars.get('output', False)
    sourcecatalogdict = {}
    for imgname in imglist:
        print("Image name: ", imgname)

        sourcecatalogdict[imgname] = {}

        # open image
        imghdu = fits.open(imgname)
        imgprimaryheader = imghdu[0].header
        instrument = imgprimaryheader['INSTRUME'].lower()
        detector = imgprimaryheader['DETECTOR'].lower()

        # get instrument/detector-specific image alignment parameters
        if instrument in detector_specific_params.keys():
            if detector in detector_specific_params[instrument].keys():
                detector_pars = detector_specific_params[instrument][detector]
                # to allow generate_source_catalog to get detector specific parameters
                detector_pars.update(pars)
                sourcecatalogdict[imgname]["params"] = detector_pars
            else:
                sys.exit("ERROR! Unrecognized detector '{}'. Exiting...".format(detector))
        else:
            sys.exit("ERROR! Unrecognized instrument '{}'. Exiting...".format(instrument))

        # Identify sources in image, convert coords from chip x, y form to reference WCS sky RA, Dec form.
        imgwcs = HSTWCS(imghdu, 1)
        fwhmpsf_pix = sourcecatalogdict[imgname]["params"]['fwhmpsf']/imgwcs.pscale #Convert fwhmpsf from arsec to pixels

        sourcecatalogdict[imgname]["catalog_table"] = amutils.generate_source_catalog(imghdu, fwhm=fwhmpsf_pix, **detector_pars)

        # write out coord lists to files for diagnostic purposes. Protip: To display the sources in these files in DS9,
        # set the "Coordinate System" option to "Physical" when loading the region file.
        imgroot = os.path.basename(imgname).split('_')[0]
        numSci = amutils.countExtn(imghdu)
        # Allow user to decide when and how to write out catalogs to files
        if output:
            for chip in range(1,numSci+1):
                regfilename = "{}_sci{}_src.reg".format(imgroot, chip)
                out_table = Table(sourcecatalogdict[imgname]["catalog_table"][chip])
                out_table.write(regfilename, include_names=["xcentroid", "ycentroid"], format="ascii.fast_commented_header")
                print("Wrote region file {}\n".format(regfilename))
        imghdu.close()
    return(sourcecatalogdict)


# ----------------------------------------------------------------------------------------------------------------------


def update_image_wcs_info(tweakwcs_output):
    """Write newly computed WCS information to image headers and write headerlet files

        Parameters
        ----------
        tweakwcs_output : list
            output of tweakwcs. Contains sourcelist tables, newly computed WCS info, etc. for every chip of every valid
            input image.

        Returns
        -------
        out_headerlet_list : dictionary
            a dictionary of the headerlet files created by this subroutine, keyed by flt/flc fits filename.
        """
    out_headerlet_dict = {}
    for item in tweakwcs_output:
        imageName = item.meta['filename']
        chipnum = item.meta['chip']
        if chipnum == 1:
            chipctr = 1
            hdulist = fits.open(imageName, mode='update')
            num_sci_ext = amutils.countExtn(hdulist)

            # generate wcs name for updated image header, headerlet
            if not hdulist['SCI',1].header['WCSNAME'] or hdulist['SCI',1].header['WCSNAME'] =="": #Just in case header value 'wcsname' is empty.
                wcsName = "FIT_{}".format(item.meta['catalog_name'])
            else:
                wname = hdulist['sci', 1].header['wcsname']
                if "-" in wname:
                    wcsName = '{}-FIT_{}'.format(wname[:wname.index('-')], item.meta['tweakwcs_info']['catalog'])
                else:
                    wcsName = '{}-FIT_{}'.format(wname, item.meta['tweakwcs_info']['catalog'])

            # establish correct mapping to the science extensions
            sciExtDict = {}
            for sciExtCtr in range(1, num_sci_ext + 1):
                sciExtDict["{}".format(sciExtCtr)] = fileutil.findExtname(hdulist,'sci',extver=sciExtCtr)

        # update header with new WCS info
        updatehdr.update_wcs(hdulist, sciExtDict["{}".format(item.meta['chip'])], item.wcs, wcsname=wcsName,
                                 reusename=True, verbose=True)
        if chipctr == num_sci_ext:
            # Close updated flc.fits or flt.fits file
            print("CLOSE {}\n".format(imageName))  # TODO: Remove before deployment
            hdulist.flush()
            hdulist.close()

            # Create headerlet
            out_headerlet = headerlet.create_headerlet(imageName, hdrname=wcsName, wcsname=wcsName)

            # Update headerlet
            update_headerlet_phdu(item, out_headerlet)

            # Write headerlet
            if imageName.endswith("flc.fits"):
                headerlet_filename = imageName.replace("flc", "flt_hlet")
            if imageName.endswith("flt.fits"):
                headerlet_filename = imageName.replace("flt", "flt_hlet")
            out_headerlet.writeto(headerlet_filename, clobber=True)
            print("Wrote headerlet file {}.\n\n".format(headerlet_filename))
            out_headerlet_dict[imageName] = headerlet_filename

        chipctr +=1
    return (out_headerlet_dict)


# ----------------------------------------------------------------------------------------------------------------------
def update_headerlet_phdu(tweakwcs_item, headerlet):
    """Update the primary header data unit keywords of a headerlet object in-place

    Parameters
    ==========
    tweakwc_item :
        Basically the output from tweakwcs which contains the cross match and fit
        information for every chip of every valid input image.

    headerlet :
        object containing WCS information
    """

    # Get the data to be used as values for FITS keywords
    rms_ra = tweakwcs_item.meta['tweakwcs_info']['RMS_RA'].value
    rms_dec = tweakwcs_item.meta['tweakwcs_info']['RMS_DEC'].value
    fit_rms = tweakwcs_item.meta['tweakwcs_info']['FIT_RMS']
    nmatch = tweakwcs_item.meta['tweakwcs_info']['nmatches']
    catalog = tweakwcs_item.meta['tweakwcs_info']['catalog']

    x_shift = (tweakwcs_item.meta['tweakwcs_info']['shift'])[0]
    y_shift = (tweakwcs_item.meta['tweakwcs_info']['shift'])[1]
    rot = tweakwcs_item.meta['tweakwcs_info']['rot']
    scale = tweakwcs_item.meta['tweakwcs_info']['scale'][0]
    skew = tweakwcs_item.meta['tweakwcs_info']['skew']

    # Update the existing FITS keywords
    primary_header = headerlet[0].header
    primary_header['RMS_RA'] = rms_ra
    primary_header['RMS_DEC'] = rms_dec
    primary_header['NMATCH'] = nmatch
    primary_header['CATALOG'] = catalog

    # Create a new FITS keyword
    primary_header['FIT_RMS'] = (fit_rms, 'RMS (mas) of the 2D fit of the headerlet solution')

    # Create the set of HISTORY keywords
    primary_header['HISTORY'] = '~~~~~ FIT PARAMETERS ~~~~~'
    primary_header['HISTORY'] = '{:>15} : {:9.4f} "/pixels'.format('platescale', tweakwcs_item.wcs.pscale)
    primary_header['HISTORY'] = '{:>15} : {:9.4f} pixels'.format('x_shift', x_shift)
    primary_header['HISTORY'] = '{:>15} : {:9.4f} pixels'.format('y_shift', y_shift)
    primary_header['HISTORY'] = '{:>15} : {:9.4f} degrees'.format('rotation', rot)
    primary_header['HISTORY'] = '{:>15} : {:9.4f}'.format('scale', scale)
    primary_header['HISTORY'] = '{:>15} : {:9.4f}'.format('skew', skew)


# ----------------------------------------------------------------------------------------------------------------------


def interpret_fit_rms(tweakwcs_output, reference_catalog):
    """Interpret the FIT information to convert RMS to physical units

    Parameters
    ----------
    tweakwcs_output : list
        output of tweakwcs. Contains sourcelist tables, newly computed WCS info, etc. for every chip of every valid
        input image.  This list gets updated, in-place, with the new RMS values;
        specifically,

            * 'FIT_RMS': RMS of the separations between fitted image positions and reference positions
            * 'TOTAL_RMS': mean of the FIT_RMS values for all observations
            * 'NUM_FITS': number of images/group_id's with successful fits included in the TOTAL_RMS

        These entries are added to the 'tweakwcs_info' dictionary.

    reference_catalog : astropy.Table
        Table of reference source positions used for the fit

    Returns
    -------
    Nothing
    """
    # Start by collecting information by group_id
    group_ids = [info.meta['group_id'] for info in tweakwcs_output]
    group_dict = {'avg_RMS':None}
    obs_rms = []
    for group_id in group_ids:
        for item in tweakwcs_output:
            if item.meta['tweakwcs_info']['status'].startswith('FAILED'):
                continue
            if item.meta['group_id'] == group_id and \
               group_id not in group_dict:
                    group_dict[group_id] = {'ref_idx':None, 'FIT_RMS':None}
                    tinfo = item.meta['tweakwcs_info']
                    ref_idx = tinfo['fit_ref_idx']
                    group_dict[group_id]['ref_idx'] = ref_idx
                    ref_RA = reference_catalog[ref_idx]['RA']
                    ref_DEC = reference_catalog[ref_idx]['DEC']
                    img_coords = SkyCoord(tinfo['fit_RA'], tinfo['fit_DEC'],
                                          unit='deg',frame='icrs')
                    ref_coords = SkyCoord(ref_RA, ref_DEC, unit='deg',frame='icrs')
                    dra, ddec = img_coords.spherical_offsets_to(ref_coords)
                    ra_rms = np.std(dra.to(u.mas))
                    dec_rms = np.std(ddec.to(u.mas))
                    fit_rms = np.std(Angle(img_coords.separation(ref_coords), unit=u.mas)).value
                    group_dict[group_id]['FIT_RMS'] = fit_rms
                    group_dict[group_id]['RMS_RA'] = ra_rms
                    group_dict[group_id]['RMS_DEC'] = dec_rms
                    obs_rms.append(fit_rms)
    # Compute RMS for entire ASN/observation set
    total_rms = np.mean(obs_rms)
    #total_rms = np.sqrt(np.sum(np.array(obs_rms)**2))

    # Now, append computed results to tweakwcs_output
    for item in tweakwcs_output:
        group_id = item.meta['group_id']
        if group_id in group_dict:
            fit_rms = group_dict[group_id]['FIT_RMS']
            ra_rms = group_dict[group_id]['RMS_RA']
            dec_rms = group_dict[group_id]['RMS_DEC']
        else:
            fit_rms = None
            ra_rms = None
            dec_rms = None

        item.meta['tweakwcs_info']['FIT_RMS'] = fit_rms
        item.meta['tweakwcs_info']['TOTAL_RMS'] = total_rms
        item.meta['tweakwcs_info']['NUM_FITS'] = len(group_ids)
        item.meta['tweakwcs_info']['RMS_RA'] = ra_rms
        item.meta['tweakwcs_info']['RMS_DEC'] = dec_rms
        item.meta['tweakwcs_info']['catalog'] = reference_catalog.meta['catalog']


# ----------------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    import argparse
    PARSER = argparse.ArgumentParser(description='Align images')
    PARSER.add_argument('raw_input_list', nargs='+', help='A space-separated list of fits files to align, or a simple '
                    'text file containing a list of fits files to align, one per line')

    PARSER.add_argument( '-a', '--archive', required=False,choices=['True','False'],default='False',help='Retain '
                    'copies of the downloaded files in the astroquery created sub-directories? Unless explicitly set, '
                    'the default is "False".')

    PARSER.add_argument( '-c', '--clobber', required=False,choices=['True','False'],default='False',help='Download and '
                    'overwrite existing local copies of input files? Unless explicitly set, the default is "False".')

    PARSER.add_argument( '-d', '--debug', required=False,choices=['True','False'],default='False',help='Attempt to use '
                    'saved sourcelists stored in pickle files if they exist, or if they do not exist, save sourcelists'
                    ' in pickle files for reuse so that step 4 can be skipped for faster subsequent debug/development '
                    'runs?? Unless explicitly set, the default is "False".')

    PARSER.add_argument( '-u', '--update_hdr_wcs', required=False,choices=['True','False'],default='False',help='Write '
                    'newly computed WCS information to image image headers and create headerlet files? Unless explicitly '
                    'set, the default is "False".')

    PARSER.add_argument( '-p', '--print_fit_parameters', required=False,choices=['True','False'],default='True',help=''
                    'Specify whether or not to print out FIT results for each chip. Unless explicitly set, the default '
                    'is "True".')

    PARSER.add_argument( '-g', '--print_git_info', required=False,choices=['True','False'],default='False',help='Display '
                    'git repository information? Unless explicitly set, the default is "False".')

    PARSER.add_argument( '-o', '--output', required=False,choices=['True','False'],default='False',help='Should '
                    'utils.astrometric_utils.create_astrometric_catalog() generate file "ref_cat.ecsv" and should '
                    'generate_source_catalogs() generate the .reg region files for every chip of every input image and '
                    'should generate_astrometric_catalog() generate file "refcatalog.cat"? Unless explicitly set, the '
                    'default is "False".')
    'Should the code generate '

    ARGS = PARSER.parse_args()

    # Build list of input images
    input_list = []
    for item in ARGS.raw_input_list:
        if os.path.exists(item):
            with open(item, 'r') as infile:
                fileLines = infile.readlines()
            for fileLine in fileLines:
                input_list.append(fileLine.strip())
        else:
            input_list.append(item)

    # Convert input args from text strings to Boolean True/False values
    archive = convert_string_tf_to_boolean(ARGS.archive)

    clobber = convert_string_tf_to_boolean(ARGS.clobber)

    debug = convert_string_tf_to_boolean(ARGS.debug)

    update_hdr_wcs = convert_string_tf_to_boolean(ARGS.update_hdr_wcs)

    print_fit_parameters = convert_string_tf_to_boolean(ARGS.print_fit_parameters)

    print_git_info = convert_string_tf_to_boolean(ARGS.print_git_info)

    output = convert_string_tf_to_boolean(ARGS.output)

    # Get to it!
    return_value = perform_align(input_list,archive,clobber,debug,update_hdr_wcs,print_fit_parameters,print_git_info)

    print(return_value)
