#!/usr/bin/env python

"""This script is a modernized implementation of tweakreg.

"""
import copy
import datetime
import sys
import glob
import math
import os
import pickle
from collections import OrderedDict
import traceback

import numpy as np
from astropy.table import Table
from astropy.io import fits

from stsci.tools import logutil

from . import util
from .hlautils import astrometric_utils as amutils
from .hlautils import astroquery_utils as aqutils
from .hlautils import get_git_rev_info
from .hlautils import align_utils
from .hlautils import config_utils

__taskname__ = 'align'

MIN_CATALOG_THRESHOLD = 3
MIN_OBSERVABLE_THRESHOLD = 10
MIN_CROSS_MATCHES = 3
MIN_FIT_MATCHES = 6
MAX_FIT_RMS = 10  # RMS now in mas, 1.0
MAX_FIT_LIMIT = 150  # Maximum RMS that a result is useful
MAX_SOURCES_PER_CHIP = 250  # Maximum number of sources per chip to include in source catalog
# MAX_RMS_RATIO = 1.0  # Maximum ratio between RMS in RA and DEC which still represents a valid fit
MAS_TO_ARCSEC = 1000.  # Conversion factor from milli-arcseconds to arcseconds


MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)

__version__ = 0.0
__version_date__ = '21-Aug-2019'

# ----------------------------------------------------------------------------------------------------------
def check_and_get_data(input_list, **pars):
    """Verify that all specified files are present. If not, retrieve them from MAST.

    Parameters
    ----------
    input_list : list
        List of one or more calibrated fits images that will be used for catalog generation.

    Returns
    =======
    total_input_list: list
        list of full filenames

    """
    empty_list = []
    retrieve_list = []    # Actual files retrieved via astroquery and resident on disk
    candidate_list = []   # File names gathered from *_asn.fits file
    ipppssoot_list = []   # ipppssoot names used to avoid duplicate downloads
    total_input_list = []  # Output full filename list of data on disk
    member_suffix = '_flc.fits'

    # Loop over the input_list to determine if the item in the input_list is a full association file
    # (*_asn.fits), a full individual image file (aka singleton, *_flt.fits), or a root name specification
    # (association or singleton, ipppssoot).
    for input_item in input_list:
        log.info('Input item: {}'.format(input_item))
        indx = input_item.find('_')

        # Input with a suffix (_xxx.fits)
        if indx != -1:
            lc_input_item = input_item.lower()
            suffix = lc_input_item[indx + 1:indx + 4]
            log.info('file: {}'.format(lc_input_item))
            # For an association, need to open the table and read the image names as this could
            # be a custom association.  The assumption is this file is on local disk when specified
            # in this manner (vs just the ipppssoot of the association).
            # This "if" block just collects the wanted full file names.
            if suffix == 'asn':
                try:
                    asntab = Table.read(input_item, format='fits')
                except FileNotFoundError:
                    log.error('File {} not found.'.format(input_item))
                    return(empty_list)
                for row in asntab:
                    if row['MEMTYPE'].startswith('PROD'):
                        continue
                    memname = row['MEMNAME'].lower().strip()
                    # Need to check if the MEMNAME is a full filename or an ipppssoot
                    if memname.find('_') != -1:
                        candidate_list.append(memname)
                    else:
                        # Define suffix for all members based on what files are present
                        if not os.path.exists(memname + member_suffix):
                            member_suffix = '_flt.fits'

                        candidate_list.append(memname + member_suffix)
            elif suffix in ['flc', 'flt']:
                if lc_input_item not in candidate_list:
                    candidate_list.append(lc_input_item)
            else:
                log.error(
                    'Inappropriate file suffix: {}.  Looking for "asn.fits", '
                    '"flc.fits", or "flt.fits".'.format(
                        suffix))
                return (empty_list)

        # Input is an ipppssoot (association or singleton), nine characters by definition.
        # This "else" block actually downloads the data specified as ipppssoot.
        elif len(input_item) == 9:
            try:
                if input_item not in ipppssoot_list:
                    # An ipppssoot of an individual file which is part of an association cannot be
                    # retrieved from MAST
                    retrieve_list = aqutils.retrieve_observation(input_item, **pars)

                    # If the retrieved list is not empty, add filename(s) to the total_input_list.
                    # Also, update the ipppssoot_list so we do not try to download the data again.  Need
                    # to do this since retrieve_list can be empty because (1) data cannot be acquired (error)
                    # or (2) data is already on disk (ok).
                    if retrieve_list:
                        total_input_list += retrieve_list
                        ipppssoot_list.append(input_item)
                    else:
                        log.error('File {} cannot be retrieved from MAST.'.format(input_item))
                        return(empty_list)
            except Exception:
                exc_type, exc_value, exc_tb = sys.exc_info()
                traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)

    # Only the retrieve_list files via astroquery have been put into the total_input_list thus far.
    # Now check candidate_list to detect or acquire the requested files from MAST via astroquery.
    for file in candidate_list:
        # If the file is found on disk, add it to the total_input_list and continue
        if glob.glob(file):
            total_input_list.append(file)
            continue
        else:
            log.error('File {} cannot be found on the local disk.'.format(file))
            return(empty_list)

    log.info("TOTAL INPUT LIST: {}".format(total_input_list))
    return(total_input_list)

# ------------------------------------------------------------------------------------------------------------
def perform_align(input_list, archive=False, clobber=False, debug=False, update_hdr_wcs=False, result=None,
              runfile=None, print_fit_parameters=True, print_git_info=False, output=False, num_sources=500,
              headerlet_filenames=None, catalog_list=['GAIADR2', 'GAIADR1'], fit_label=None,
              **alignment_pars):
    """Actual Main calling function.

    Parameters
    ----------
    input_list : list
        List of one or more IPPSSOOTs (rootnames) to align.

    archive : Boolean
        Retain copies of the downloaded files in the astroquery created
        sub-directories?

    clobber : Boolean
        Download and overwrite existing local copies of input files?

    debug : Boolean
        Attempt to use saved sourcelists stored in pickle files if they exist, or if they do not exist, save
        sourcelists in pickle files for reuse so that step 4 can be skipped for faster subsequent
        debug/development runs??

    update_hdr_wcs : Boolean
        Write newly computed WCS information to image image headers?

    result: Table
        name of variable to be updated by subroutine run.

    runfile : string
        log file name

    print_fit_parameters : Boolean
        Specify whether or not to print out FIT results for each chip.

    print_git_info : Boolean
        Display git repository information?

    output : Boolean
        Should utils.astrometric_utils.create_astrometric_catalog() generate file 'ref_cat.ecsv' and should
        the alignment source catalogs get written out to files?

    num_sources : int, optional
        Maximum number of **brightest sources per chip** which will be used for cross-matching and fitting.
        If set to None, all sources will be used.

    headerlet_filenames : dictionary, optional
        dictionary that maps the flt/flc.fits file name to the corresponding custom headerlet filename.

    Updates
    -------
    filtered_table: Astropy Table
        Table which contains processing information and alignment results for every raw image evaluated

    """
    log.info("*** HLAPIPELINE Processing Version {!s} ({!s}) started at: {!s} ***\n".format(__version__, __version_date__, util._ptime()[0]))

    # 0: print git info
    if print_git_info:
        log.info("{} STEP 0: Display Git revision info  {}".format("-" * 20, "-" * 49))
        full_path = os.path.dirname(__file__)
        repo_path = None
        if "drizzlepac/drizzlepac" in full_path:
            repo_path = full_path.split("drizzlepac/drizzlepac")[0] + "drizzlepac"
        elif "hlapipeline" in full_path:
            repo_path = full_path.split("drizzlepac")[0] + "drizzlepac"
        else:
            pass
        if not os.path.exists(repo_path): repo_path = None  # protect against non-existent paths
        if repo_path:
            get_git_rev_info.print_rev_id(repo_path)  # Display git repository information
        else:
            log.warning("WARNING: Unable to display Git repository revision information.")

    try:
        # Initialize key variables
        filtered_table = None

        # 1: Interpret input data and optional parameters
        log.info("{} STEP 1: Get data {}".format("-" * 20, "-" * 66))
        zero_dt = starting_dt = datetime.datetime.now()
        log.info(str(starting_dt))
        imglist = check_and_get_data(input_list, archive=archive, clobber=clobber)
        log.info("SUCCESS")

        log.info(make_label('Processing time of [STEP 1]', starting_dt))
        starting_dt = datetime.datetime.now()

        # Get default alignment parameters if not provided by the user...
        inst = fits.getval(imglist[0], 'instrume')
        det = fits.getval(imglist[0], 'detector')
        apars = get_default_pars(inst, det)
        alignment_pars.update(apars)

        # Instantiate AlignmentTable class with these input files
        alignment_table = align_utils.AlignmentTable(imglist, **alignment_pars)
        process_list = alignment_table.process_list

        # Define fitting algorithm list in priority order
        # The match_relative_fit algorithm must have more than one image as the first image is
        # the reference for the remaining images.
        fit_algorithm_list = alignment_table.get_fit_methods()

        if len(process_list) == 1:
            fit_algorithm_list.remove("relative")

        log.info(make_label('Processing time of [STEP 2]', starting_dt))
        starting_dt = datetime.datetime.now()
        # 3: Build WCS for full set of input observations
        log.info("{} STEP 3: Build WCS {}".format("-" * 20, "-" * 65))
        # refwcs = amutils.build_reference_wcs(process_list)
        log.info("SUCCESS")

        log.info(make_label('Processing time of [STEP 3]', starting_dt))
        starting_dt = datetime.datetime.now()
        # 4: Extract catalog of observable sources from each input image
        log.info(
            "{} STEP 4: Source finding {}".format("-" * 20, "-" * 60))
        if debug:
            pickle_filename = "{}.source_catalog.pickle".format(process_list[0])
            if os.path.exists(pickle_filename):
                pickle_in = open(pickle_filename, "rb")
                alignment_table.extracted_sources = pickle.load(pickle_in)
                log.info("Using sourcelist extracted from {} generated during the last run to save time.".format(
                    pickle_filename))
            else:
                alignment_table.find_alignment_sources(output=True)

                pickle_out = open(pickle_filename, "wb")
                pickle.dump(alignment_table.extracted_sources, pickle_out)
                pickle_out.close()
                log.info("Wrote {}".format(pickle_filename))
        else:
            alignment_table.find_alignment_sources(output=output)


        alignment_table.configure_fit()

        for imgname in alignment_table.extracted_sources.keys():
            table = alignment_table.extracted_sources[imgname]

            # Get the location of the current image in the filtered table
            index = np.where(alignment_table.filtered_table['imageName'] == imgname)[0][0]

            # First ensure sources were found

            if table is None or not table[1]:
                log.warning("No sources found in image {}".format(imgname))
                alignment_table.filtered_table[:]['status'] = 1
                alignment_table.filtered_table[:]['processMsg'] = "No sources found"
                log.info(make_label('Processing time of [STEP 4]', starting_dt))
                alignment_table.close()
                return

            # The catalog of observable sources must have at least MIN_OBSERVABLE_THRESHOLD entries to be useful
            total_num_sources = 0
            for chipnum in table.keys():
                total_num_sources += len(table[chipnum])

            # Update filtered table with number of found sources
            alignment_table.filtered_table[index]['foundSources'] = total_num_sources

            if total_num_sources < MIN_OBSERVABLE_THRESHOLD:
                log.warning("Not enough sources ({}) found in image {}".format(total_num_sources, imgname))
                alignment_table.filtered_table[:]['status'] = 1
                alignment_table.filtered_table[:]['processMsg'] = "Not enough sources found"
                log.info(make_label('Processing time of [STEP 4]', starting_dt))
                alignment_table.close()
                return
        log.info("SUCCESS")
        log.info(make_label('Processing time of [STEP 4]', starting_dt))
        starting_dt = datetime.datetime.now()
        # 5: Retrieve list of astrometric sources from database

        # Convert input images to tweakwcs-compatible FITSWCS objects and
        # attach source catalogs to them.
        imglist = []
        for group_id, image in enumerate(process_list):
            img = amutils.build_wcscat(image, group_id,
                                       alignment_table.extracted_sources[image])
            imglist.extend(img)
        # store mapping of group_id to filename/chip
        group_id_dict = {}
        for image in imglist:
            group_id_dict["{}_{}".format(image.meta["filename"], image.meta["chip"])] = image.meta["group_id"]

        best_fit_rms = -99999.0
        best_fit_status_dict = {}
        best_fit_qual = 5
        best_fit_label = [None, None]
        # create pristine copy of imglist that will be used to restore imglist back so it always starts exactly the same
        # for each run.
        orig_imglist = copy.deepcopy(imglist)
        # create dummy list that will be used to preserve imglist best_meta information through the imglist reset process
        # best_imglist = []
        fit_info_dict = OrderedDict()
        for algorithm_name in fit_algorithm_list:  # loop over fit algorithm type
            log.info("Applying {} fit method".format(algorithm_name))
            for catalog_index, catalog_name in enumerate(catalog_list):  # loop over astrometric catalog
                log.info("{} STEP 5: Detect astrometric sources {}".format("-" * 20, "-" * 48))
                log.info("Astrometric Catalog: {}".format(catalog_name))
                # store reference catalogs in a dictionary so that generate_astrometric_catalog() doesn't
                #  execute unnecessarily after it's been run once for a given astrometric catalog.
                if catalog_name in alignment_table.reference_catalogs:
                    log.info("Using {} reference catalog from earlier this run.".format(catalog_name))
                    reference_catalog = alignment_table.reference_catalogs[catalog_name]
                else:
                    log.info("Generating new reference catalog for {};"
                             " Storing it for potential re-use later this run.".format(catalog_name))
                    reference_catalog = generate_astrometric_catalog(process_list,
                                                                     catalog=catalog_name,
                                                                     output=output)
                    alignment_table.reference_catalogs[catalog_name] = reference_catalog
                log.info(make_label('Processing time of [STEP 5]', starting_dt))
                starting_dt = datetime.datetime.now()

                if len(reference_catalog) < MIN_CATALOG_THRESHOLD:
                    log.warning("Not enough sources found in catalog {}".format(catalog_name))
                    fit_quality = 5
                    if catalog_index < len(catalog_list) - 1:
                        log.info("Try again with other catalog")
                    else:
                        # bail out if not enough sources can be found any of the astrometric catalogs
                        log.warning("ERROR! No astrometric sources found in any catalog. Exiting...")
                        alignment_table.filtered_table['status'][:] = 1
                        alignment_table.filtered_table['processMsg'][:] = "No astrometric sources found"
                        alignment_table.filtered_table['fit_qual'][:] = fit_quality
                        current_dt = datetime.datetime.now()
                        delta_dt = (current_dt - starting_dt).total_seconds()
                        log.info('Processing time of [STEP 5]: {} sec'.format(delta_dt))
                        alignment_table.close()
                        return (alignment_table.filtered_table)
                else:
                    log.info("{} Cross matching and "
                         "fitting {}".format("-" * 20, "-" * 47))
                    imglist = copy.deepcopy(orig_imglist)  # reset imglist to pristine state

                    log.info(
                        "{} Catalog {} matched using {} {}".format("-" * 18,
                                                                   catalog_name,
                                                                   algorithm_name, "-" * 18))
                    try:
                        # restore group IDs to their pristine state prior to each run.
                        alignment_table.reset_group_id(len(reference_catalog))

                        # execute the correct fitting/matching algorithm
                        imglist = alignment_table.perform_fit(algorithm_name, catalog_name, reference_catalog)

                        # determine the quality of the fit
                        fit_rms, fit_num, fit_quality, filtered_table, fit_status_dict = \
                            determine_fit_quality(
                                alignment_table.imglist,
                                alignment_table.filtered_table,
                                (catalog_index < (len(catalog_list) - 1)),
                                print_fit_parameters=print_fit_parameters)
                        alignment_table.filtered_table = filtered_table

                        # save fit algorithm name to dictionary key "fit method" in imglist.
                        for imglist_ctr in range(0, len(imglist)):
                            table_fit = alignment_table.fit_dict[(catalog_name, algorithm_name)]
                            table_fit[imglist_ctr].meta['fit method'] = algorithm_name
                            table_fit[imglist_ctr].meta['fit quality'] = fit_quality
                        #print(alignment_table.fit_dict[(catalog_name, algorithm_name)][0].meta)

                        # populate fit_info_dict
                        fit_info_dict["{} {}".format(catalog_name, algorithm_name)] = \
                            fit_status_dict[next(iter(fit_status_dict))]
                        fit_info_dict["{} {}".format(catalog_name,
                            algorithm_name)]['fit_qual'] = fit_quality

                        # Figure out which fit solution to go with based on fit_quality value and maybe also total_rms
                        if fit_quality < 5:
                            if fit_quality == 1:  # valid, non-comprimised solution with total rms < 10 mas...go with this solution.
                                best_fit_rms = fit_rms
                                best_fit_label = (catalog_name, algorithm_name)

                                best_fit_status_dict = fit_status_dict.copy()
                                best_fit_qual = fit_quality
                                break  # break out of while loop
                            elif fit_quality < best_fit_qual:  # better solution found. keep looping but with the better solution as "best" for now.
                                log.info("Better solution found!")
                                best_fit_rms = fit_rms
                                best_fit_label = (catalog_name, algorithm_name)

                                best_fit_status_dict = fit_status_dict.copy()
                                best_fit_qual = fit_quality
                            elif fit_quality == best_fit_qual:  # new solution same level of fit_quality. Choose whichever one has the lowest total rms as "best" and keep looping.
                                if best_fit_rms >= 0.:
                                    if fit_rms < best_fit_rms:
                                        best_fit_rms = fit_rms
                                        best_fit_label = (catalog_name, algorithm_name)

                                        best_fit_status_dict = fit_status_dict.copy()
                                        best_fit_qual = fit_quality
                            else:  # new solution has worse fit_quality. discard and continue looping.
                                continue

                    except Exception:
                        exc_type, exc_value, exc_tb = sys.exc_info()
                        traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)
                        log.warning(
                            "WARNING: Catastrophic fitting failure with catalog {} and matching "
                            "algorithm {}.".format(catalog_name,
                                                   algorithm_name))
                        alignment_table.filtered_table['status'][:] = 1
                        alignment_table.filtered_table['processMsg'][:] = "Fitting failure"
                        # It may be there are additional catalogs and algorithms to try, so keep going
                        fit_quality = 5  # Flag this fit with the 'bad' quality value
                        alignment_table.filtered_table['fit_qual'][:] = fit_quality
                        continue
                    if fit_quality == 1:  # break out of inner  astrometric catalog loop
                        break
            # break out of outer fit algorithm loop
            # either with a fit_rms < 10 or a 'valid' relative fit
            if fit_quality == 1 or (best_fit_qual in [2, 3, 4] and
                "relative" in algorithm_name):
                break
        log.info("best_fit found to be: {}".format(best_fit_label))
        log.info("FIT_DICT: {}".format(alignment_table.fit_dict.keys()))
        # Reset imglist to point to best solution...
        alignment_table.select_fit(best_fit_label[0], best_fit_label[1])
        imglist = alignment_table.selected_fit
        filtered_table = alignment_table.filtered_table

        # Report processing time for this step
        log.info(make_label('Processing time of [STEP 5b]', starting_dt))
        starting_dt = datetime.datetime.now()

        # 6: Populate the filtered_table
        log.info(
            "{} STEP 6: Collect up information and populate the filtered table "
            "{}".format("-" * 20, "-" * 20))
        if 0 < best_fit_rms < MAX_FIT_RMS:
            log.info("The fitting process was successful with a best fit total "
                     "rms of {} mas".format(best_fit_rms))
        else:
            log.info(
                "The fitting process was unsuccessful with a best fit total rms "
                "of {} mas".format(best_fit_rms))

        if 0 < best_fit_rms < MAX_FIT_LIMIT:
            # Update filtered table with best fit results
            filtered_table['status'][:] = 0
            fit_status_dict = best_fit_status_dict.copy()

            for item in imglist:
                imgname = item.meta['name']
                index = np.where(filtered_table['imageName'] == imgname)[0][0]

                # populate self.filtered_table fields "status", "compromised" and
                # "processMsg" with fit_status_dict fields "valid", "compromised"
                # and "reason".
                explicit_dict_key = "{},{}".format(item.meta['name'], item.meta['chip'])
                if fit_status_dict[explicit_dict_key]['valid'] is True:
                    filtered_table[index]['status'] = 0
                else:
                    filtered_table[index]['status'] = 1
                if fit_status_dict[explicit_dict_key]['compromised'] is False:
                    filtered_table['compromised'] = 0
                else:
                    filtered_table['compromised'] = 1

                filtered_table[index]['processMsg'] = fit_status_dict[explicit_dict_key]['reason']
                filtered_table['fit_qual'][index] = item.meta['fit quality']

        log.info(make_label('Processing time of [STEP 6]', starting_dt))
        starting_dt = datetime.datetime.now()
        # 7: Write new fit solution to input image headers
        log.info("{} STEP 7: Update image headers with new WCS information "
                 "{}".format("-" * 20, "-" * 29))
        if (0 < best_fit_rms < 9999.) and update_hdr_wcs:
            # determine the quality of the fit
            alignment_table.apply_fit(headerlet_filenames)
            log.info("SUCCESS")
        else:
            log.info(" STEP SKIPPED")

        log.info(make_label('Processing time of [STEP 7]', starting_dt))
        log.info('TOTAL Processing time of {} sec'.format(
                 (datetime.datetime.now() - zero_dt).total_seconds()))
        log.info(best_fit_status_dict)
        log.info("-" * 104)

        log.info("-" * 104)
        log.info("                             SUMMARY OF ALL FIT ATTEMPTS")
        for item in fit_info_dict.keys():
            log.info("{} {}".format(item, fit_info_dict[item]))
        log.info("-" * 104)

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)

    finally:
        # Always make sure that all file handles are closed
        alignment_table.close()

        # Now update the result with the filtered_table contents
        if result:
            result.meta = filtered_table.meta
            for col in filtered_table.colnames:
                result.add_column(filtered_table[col], name=col)
        if filtered_table is not None:
            filtered_table.pprint(max_width=-1)
    return filtered_table



# ----------------------------------------------------------------------------------------------------------

def make_label(label, starting_dt):
    """Create a time-stamped label for use in log messages"""
    current_dt = datetime.datetime.now()
    delta_dt = (current_dt - starting_dt).total_seconds()
    return '{}: {} sec'.format(label, delta_dt)



# ----------------------------------------------------------------------------------------------------------


def determine_fit_quality(imglist, filtered_table, catalogs_remaining, print_fit_parameters=True):
    """Determine the quality of the fit to the data

    Parameters
    ----------
    imglist : list
        output of interpret_fits. Contains sourcelist tables, newly computed WCS info, etc. for every chip of
        every valid input image.  This list should have been  updated, in-place, with the new RMS values;
        specifically,

            * 'FIT_RMS': RMS of the separations between fitted image positions and reference positions
            * 'TOTAL_RMS': mean of the FIT_RMS values for all observations
            * 'NUM_FITS': number of images/group_id's with successful fits included in the TOTAL_RMS

        These entries are added to the 'fit_info' dictionary.

    filtered_table : object
        Astropy Table object containing data pertaining to the associated dataset, including
        the doProcess bool.  It is intended this table is updated by subsequent functions for
        bookkeeping purposes.

    catalogs_remaining : bool
        Specify whether additional catalogs remain to be fit against.

    print_fit_parameters : bool
        Specify whether or not to print out FIT results for each chip

    Returns
    -------
    max_rms_val : float
        The best Total rms determined from all of the images

    num_xmatches: int
        The number of stars used in matching the data


    fit_quality : int
        fit quality category:
            * 1 = valid solution with rms < 10 mas
            * 2 = Valid but compromised solution with rms < 10 mas
            * 3 = Valid solution with RMS >= 10 mas
            * 4 = Valid but compromised solution with RMS >= 10 mas
            * 5 = Not valid solution

    filtered_table : object
        modified filtered_table object

    fit_status_dict : dictionary
        Dictionary containing the following:
            * overall fit validity (Boolean)
            * total (visit-level) RMS value in mas (float)
            * number of matched sources (int)
            * fit compromised status (Boolean)
            * reason fit is considered 'compromised' (only populated if "compromised" field is "True")
    """
    tweakwcs_info_keys = OrderedDict(imglist[0].meta['fit_info']).keys()
    max_rms_val = 1e9
    num_xmatches = 0
    fit_status_dict = {}
    xshifts = []
    yshifts = []
    overall_valid = True
    overall_comp = False
    for item in imglist:
        if item.meta['fit_info']['status'].startswith('FAILED') is False:
            xshifts.append(item.meta['fit_info']['shift'][0])
            yshifts.append(item.meta['fit_info']['shift'][1])

    for item in imglist:
        image_name = item.meta['name']
        chip_num = item.meta['chip']

        # Build fit_status_dict entry
        dict_key = "{},{}".format(image_name, chip_num)
        fit_status_dict[dict_key] = {'valid': False,
                                  'max_rms': max_rms_val,
                                  'num_matches': num_xmatches,
                                  'compromised': False,
                                  'reason': ""}
        # Handle fitting failures (no matches found)
        if item.meta['fit_info']['status'].startswith("FAILED") is True:
            log.warning("No cross matches found in any catalog for {} "
                        "- no processing done.".format(image_name))
            overall_valid = False
            continue
        fit_rms_val = item.meta['fit_info']['FIT_RMS']
        max_rms_val = item.meta['fit_info']['TOTAL_RMS']
        # fit_rms_ra = item.meta['fit_info']['RMS_RA']
        # fit_rms_dec = item.meta['fit_info']['RMS_DEC']
        # rms_ratio = abs(fit_rms_ra - fit_rms_dec) / min(fit_rms_ra, fit_rms_dec)
        num_xmatches = item.meta['fit_info']['nmatches']
        fit_status_dict[dict_key]['max_rms'] = max_rms_val
        fit_status_dict[dict_key]['num_matches'] = num_xmatches

        if num_xmatches < MIN_CROSS_MATCHES:
            if catalogs_remaining:
                log.warning(
                    "Not enough cross matches found between astrometric"
                    "catalog and sources found in {}".format(image_name))
                overall_valid = False
                continue

        # Compute correlation between input and GAIA magnitudes
        if num_xmatches < max(0.1 * item.meta['num_ref_catalog'], 10):
            cross_match_check = amutils.check_mag_corr([item])[0]
            log.info("Cross-match check: {} on {} ref sources".format(cross_match_check,
                                                                      item.meta['num_ref_catalog']))
        else:
            cross_match_check = True

        # Execute checks
        nmatches_check = False
        if num_xmatches > 4 or (num_xmatches > 2 and fit_rms_val > 0.5):
            nmatches_check = True

        radial_offset_check = False
        radial_offset = math.sqrt(
            float(item.meta['fit_info']['shift'][0])**2 +
            float(item.meta['fit_info']['shift'][0])**2) * item.wcs.pscale  # radial offset in arssec
        if float(num_xmatches) * 0.36 > 0.8 + (radial_offset / 10.0)**8:
            radial_offset_check = True

        large_rms_check = True
        if fit_rms_val > 150. or max_rms_val > 150.:
            large_rms_check = False

        # fitRmsCheck = False
        # if fit_rms_val < max_rms_val:
        #     fitRmsCheck = True

        consistency_check = True
        rms_limit = max(item.meta['fit_info']['TOTAL_RMS'], 10.)
        if not math.sqrt(np.std(np.asarray(xshifts)) ** 2 + np.std(
                         np.asarray(yshifts)) ** 2) <= (rms_limit / MAS_TO_ARCSEC) / (item.wcs.pscale):  # \
                         # or rms_ratio > MAX_RMS_RATIO:
            consistency_check = False

        # Decide if fit solutions are valid based on checks
        if not consistency_check:  # Failed consistency check
            fit_status_dict[dict_key]['valid'] = False
            fit_status_dict[dict_key]['compromised'] = False
            fit_status_dict[dict_key]['reason'] = "Consistency violation!"
        elif not large_rms_check:  # RMS value(s) too large
            fit_status_dict[dict_key]['valid'] = False
            fit_status_dict[dict_key]['compromised'] = False
            fit_status_dict[dict_key]['reason'] = "RMS too large (>150 mas)!"
        elif not radial_offset_check:  # Failed radial offset check
            fit_status_dict[dict_key]['valid'] = False
            fit_status_dict[dict_key]['compromised'] = False
            fit_status_dict[dict_key]['reason'] = "Radial offset value too large!"
        elif not nmatches_check:  # Too few matches
            fit_status_dict[dict_key]['valid'] = False
            fit_status_dict[dict_key]['compromised'] = True
            fit_status_dict[dict_key]['reason'] = "Too few matches!"
        elif not cross_match_check:
            fit_status_dict[dict_key]['valid'] = True
            fit_status_dict[dict_key]['compromised'] = True
            fit_status_dict[dict_key]['reason'] = "Cross-match magnitudes not correlated!"
        else:  # all checks passed. Valid solution.
            fit_status_dict[dict_key]['valid'] = True
            fit_status_dict[dict_key]['compromised'] = False
            fit_status_dict[dict_key]['reason'] = ""
        # for now, generate overall valid and compromised values. Basically, if any of the entries for "valid" is False,
        # "valid" is False, treat the whole dataset as not valid. Same goes for compromised.
        if not fit_status_dict[dict_key]['valid']:
            overall_valid = False
        if fit_status_dict[dict_key]['compromised']:
            overall_comp = True

        log.info('RESULTS FOR {} Chip {}: FIT_RMS = {} mas, TOTAL_RMS = {}'
                 ' mas, NUM =  {}'.format(image_name,
                                          item.meta['chip'],
                                          fit_rms_val,
                                          max_rms_val,
                                          num_xmatches))
        # print fit params to screen
        if print_fit_parameters:
            log_info_keys = ['status', 'fitgeom', 'eff_minobj', 'matrix', 'shift', 'center', 'rot', 'proper',
                'rotxy', 'scale', 'skew', 'rmse', 'mae', 'nmatches', 'FIT_RMS', 'TOTAL_RMS', 'NUM_FITS',
                'RMS_RA', 'RMS_DEC', 'catalog']
            log.info("{} FIT PARAMETERS {}".format("~" * 35, "~" * 34))
            log.info("image: {}".format(image_name))
            log.info("chip: {}".format(item.meta['chip']))
            log.info("group_id: {}".format(item.meta['group_id']))
            for tweakwcs_info_key in log_info_keys:
                log.info("{} : {}".format(tweakwcs_info_key, item.meta['fit_info'][tweakwcs_info_key]))
            log.info("~" * 84)
            log.info("nmatches_check: {} radial_offset_check: {}"
                     " large_rms_check: {},"
                     " consistency_check: {}".format(nmatches_check,
                                                    radial_offset_check,
                                                    large_rms_check,
                                                    consistency_check))


    # determine which fit quality category this latest fit falls into
    if overall_valid is False:
        fit_quality = 5
        log.info("FIT SOLUTION REJECTED")
        filtered_table['status'][:] = 1
        for ctr in range(0, len(filtered_table)):
            filtered_table[ctr]['processMsg'] = fit_status_dict[filtered_table[ctr]['imageName'] + ",1"]["reason"]
    else:
        for ctr in range(0, len(filtered_table)):
            filtered_table[ctr]['processMsg'] = ""
        if overall_comp is False and max_rms_val < 10.:
            log.info("Valid solution with RMS < 10 mas found!")
            fit_quality = 1
        elif overall_comp is True and max_rms_val < 10.:
            log.info("Valid but compromised solution with RMS < 10 mas found!")
            fit_quality = 2
        elif overall_comp is False and 1000. >= max_rms_val >= 10.:
            log.info("Valid solution with RMS >= 10 mas found!")
            fit_quality = 3
        else:
            log.info("Valid but compromised solution with RMS >= 10 mas found!")
            fit_quality = 4

    if print_fit_parameters:
        for item in imglist: log.info(fit_status_dict["{},{}".format(item.meta['name'], item.meta['chip'])])

    if max_rms_val > MAX_FIT_RMS:
        log.info("Total fit RMS value = {} mas greater than the maximum threshold value {}.".format(max_rms_val, MAX_FIT_RMS))
    if not overall_valid:
        log.info("The fit solution for some or all of the images is not valid.")
    if max_rms_val > MAX_FIT_RMS or not overall_valid:
        log.info("Try again with the next catalog")
    else:
        log.info("Fit calculations successful.")
    return max_rms_val, num_xmatches, fit_quality, filtered_table, fit_status_dict


# ----------------------------------------------------------------------------------------------------------------------


def generate_astrometric_catalog(imglist, **pars):
    """Generates a catalog of all sources from an existing astrometric catalog are
       in or near the FOVs of the images in the input list.

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
    if pars['output'] is True:
        pars['output'] = 'ref_cat.ecsv'
    else:
        pars['output'] = None
    out_catalog = amutils.create_astrometric_catalog(imglist, **pars)
    pars = temp_pars.copy()
    # if the catalog has contents, write the catalog to ascii text file
    if len(out_catalog) > 0 and pars['output']:
        catalog_filename = "refcatalog.cat"
        out_catalog.write(catalog_filename, format="ascii.fast_commented_header")
        log.info("Wrote reference catalog {}.".format(catalog_filename))

    return(out_catalog)

# ----------------------------------------------------------------------------------------------------------------------

def update_headerlet_phdu(tweakwcs_item, headerlet):
    """Update the primary header data unit keywords of a headerlet object in-place

    Parameters
    ==========
    tweakwcs_item :
        Basically the output from tweakwcs which contains the cross match and fit information for every chip
        of every valid input image.

    headerlet : headerlet object
        object containing WCS information
    """

    # Get the data to be used as values for FITS keywords
    rms_ra = tweakwcs_item.meta['fit_info']['RMS_RA'].value
    rms_dec = tweakwcs_item.meta['fit_info']['RMS_DEC'].value
    fit_rms = tweakwcs_item.meta['fit_info']['FIT_RMS']
    nmatch = tweakwcs_item.meta['fit_info']['nmatches']
    catalog = tweakwcs_item.meta['fit_info']['catalog']
    fit_method = tweakwcs_item.meta['fit method']

    x_shift = (tweakwcs_item.meta['fit_info']['shift'])[0]
    y_shift = (tweakwcs_item.meta['fit_info']['shift'])[1]
    rot = tweakwcs_item.meta['fit_info']['rot']
    scale = tweakwcs_item.meta['fit_info']['scale'][0]
    skew = tweakwcs_item.meta['fit_info']['skew']

    log.info("Headerlet being updated with RMS_RA={},  RMS_DEC={}".format(rms_ra, rms_dec))
    # Update the existing FITS keywords
    primary_header = headerlet[0].header
    primary_header['RMS_RA'] = rms_ra
    primary_header['RMS_DEC'] = rms_dec
    primary_header['NMATCH'] = nmatch
    primary_header['CATALOG'] = catalog
    primary_header['FITMETH'] = fit_method

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



def get_default_pars(instrument, detector, step='alignment',
                     condition=['filter_basic']):

    step_list = config_utils.step_title_list
    if step not in step_list:
        log.error("{} not valid!  Needs to be one of: {}".format(step,
                  step_list))
        raise ValueError

    full_cfg_index, pars_dir = config_utils.read_index(instrument, detector)

    par_class = config_utils.step_name_list[step_list.index(step)]
    apars = par_class(full_cfg_index[step], condition,
                      pars_dir, step, True, None)

    return apars.outpars
