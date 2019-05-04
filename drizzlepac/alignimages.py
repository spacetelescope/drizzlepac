#!/usr/bin/env python

"""This script is a modernized implementation of tweakreg.

"""
import argparse
import copy
import datetime
import sys
import glob
import os
import pickle
from itertools import product
from collections import OrderedDict
import traceback

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u
from stsci.tools import fileutil, logutil
from stwcs.wcsutil import headerlet, HSTWCS
import tweakwcs

from drizzlepac import updatehdr
from drizzlepac import util
from drizzlepac.hlautils.astrometric_utils import (
    build_wcscat, create_astrometric_catalog, generate_source_catalog,
    countExtn
)
from drizzlepac.hlautils import astroquery_utils as aqutils
from drizzlepac.hlautils import analyze
from drizzlepac.hlautils import get_git_rev_info


__taskname__ = 'alignimages'

__version__ = 0.1
__version_date__ = '15-Feb-2019'


_MIN_CATALOG_THRESHOLD = 3
_MIN_OBSERVABLE_THRESHOLD = 10
_MIN_CROSS_MATCHES = 3
_MIN_FIT_MATCHES = 6
_MAX_FIT_RMS = 10  # RMS now in mas, 1.0
_MAX_FIT_LIMIT = 1000  # Maximum RMS that a result is useful
# Maximum number of sources per chip to include in source catalog:
_MAX_SOURCES_PER_CHIP = 250

# Module-level dictionary contains instrument/detector-specific parameters used
# later on in the script.
_DETECTOR_SPECIFIC_PARAMS = {
    "acs": {"hrc": {"fwhmpsf": 0.152,  # 0.073
                    "classify": True,
                    "threshold": None},
            "sbc": {"fwhmpsf": 0.13,  # 0.065
                    "classify": False,
                    "threshold": 2.0},
            "wfc": {"fwhmpsf": 0.13,  # 0.076,
                    "classify": True,
                    "threshold": -1.1}
            },

    "wfc3": {"ir": {"fwhmpsf": 0.25,  # 0.14
                    "classify": False,
                    "threshold": None},
             "uvis": {"fwhmpsf": 0.152,  # 0.076
                      "classify": True,
                      "threshold": None}
             }
}


log = logutil.create_logger('alignimages', level=logutil.logging.INFO,
                            stream=sys.stdout)


def check_and_get_data(input_list, **pars):
    """ Verify that all specified files are present. If not, retrieve them
    from MAST.

    Parameters
    ----------
    input_list : list
        List of one or more calibrated fits images that will be used for
        catalog generation.

    Returns
    =======
    total_input_list: list
        list of full filenames

    """
    retrieve_list = []  # Files retrieved via astroquery and resident on disk
    candidate_list = []  # File names gathered from *_asn.fits file
    ipppssoot_list = []  # ipppssoot names used to avoid duplicate downloads
    total_input_list = []  # Output full filename list of data on disk
    member_suffix = '_flc.fits'

    # Loop over the input_list to determine if the item in the input_list is a
    # full association file (*_asn.fits), a full individual image file
    # (aka singleton, *_flt.fits), or a root name specification
    # (association or singleton, ipppssoot).
    for input_item in input_list:
        log.info('Input item: {}'.format(input_item))
        indx = input_item.find('_')

        if indx != -1:  # Input with a suffix (_xxx.fits)
            lc_input_item = input_item.lower()
            suffix = lc_input_item[indx + 1:indx + 4]
            log.info('file: {}'.format(lc_input_item))
            # For an association, need to open the table and read the image
            # names as this could be a custom association.  The assumption'
            # is this file is on local disk when specified in this manner
            # (vs just the ipppssoot of the association). This "if" block
            # just collects the wanted full file names.
            if suffix == 'asn':
                try:
                    asntab = Table.read(input_item, format='fits')
                except FileNotFoundError:
                    log.error('File {} not found.'.format(input_item))
                    return []

                for row in asntab:
                    if row['MEMTYPE'].startswith('PROD'):
                        continue
                    memname = row['MEMNAME'].lower().strip()
                    # Need to check if the MEMNAME is a full filename
                    # or an ipppssoot
                    if memname.find('_') == -1:
                        memname += + '_flc.fits'
                    candidate_list.append(memname)

            elif suffix in ['flc', 'flt']:
                if lc_input_item not in candidate_list:
                    candidate_list.append(lc_input_item)

            else:
                log.error(
                    'Inappropriate file suffix: {}.  Looking for "asn.fits", '
                    '"flc.fits", or "flt.fits".'.format(
                        suffix))
                return []

        # Input is an ipppssoot (association or singleton), nine characters by
        # definition. This "else" block actually downloads the data specified
        # as ipppssoot.
        elif len(input_item) == 9:
            try:
                if input_item not in ipppssoot_list:
                    # An ipppssoot of an individual file which is part of an
                    # association cannot be retrieved from MAST
                    retrieve_list = aqutils.retrieve_observation(input_item,
                                                                 **pars)

                    # If the retrieved list is not empty, add filename(s) to
                    # the total_input_list. Also, update the ipppssoot_list
                    # so we do not try to download the data again. Need
                    # to do this since retrieve_list can be empty because
                    # (1) data cannot be acquired (error) or (2) data
                    # is already on disk (ok).
                    if retrieve_list:
                        total_input_list += retrieve_list
                        ipppssoot_list.append(input_item)
                    else:
                        log.error('File {} cannot be retrieved from MAST.'
                                  .format(input_item))
                        return []

            except Exception:
                exc_type, exc_value, exc_tb = sys.exc_info()
                traceback.print_exception(exc_type, exc_value, exc_tb,
                                          file=sys.stdout)

    # Only the retrieve_list files via astroquery have been put into the
    # total_input_list thus far. Now check candidate_list to detect or
    # acquire the requested files from MAST via astroquery.
    for file in candidate_list:
        # If the file is found on disk, add it to the total_input_list
        # and continue
        if glob.glob(file):
            total_input_list.append(file)
        else:
            log.error('File {} cannot be found on the local disk.'
                      .format(file))
            return []

    log.info('TOTAL INPUT LIST: {}'.format(total_input_list))
    return total_input_list


def _log_step(description, pad_char='-'):
    log.info("{:s} {:s} {:s}".format(
        20 * pad_char, description, max(0, 82 - len(description)) * pad_char))


@util.with_logging
def perform_align(input_list, archive=False, clobber=False, debug=False,
                  update_hdr_wcs=False, result=None, runfile=None,
                  print_fit_parameters=True, print_git_info=False,
                  output=False, num_sources=250):
    """Main calling function.


    Parameters
    ----------
    input_list : list, optional
        List of one or more IPPSSOOTs (rootnames) to align.

    archive : bool, optional
        Retain copies of the downloaded files in the astroquery created
        sub-directories?

    clobber : bool, optional
        Download and overwrite existing local copies of input files?

    debug : bool, optional
        Attempt to use saved sourcelists stored in pickle files if they exist,
        or if they do not exist, save sourcelists in pickle files for re-use
        so that step 4 can be skipped for faster subsequent debug/development
        runs??

    update_hdr_wcs : bool, optional
        Write newly computed WCS information to image image headers?

    runfile : string, optional
        log file name

    print_fit_parameters : bool, optional
        Specify whether or not to print out FIT results for each chip.

    print_git_info : bool, optional
        Display git repository information?

    output : bool, optional
        Should ``utils.astrometric_utils.create_astrometric_catalog()``
        generate file ``'ref_cat.ecsv'`` and should
        ``generate_source_catalogs()`` generate the ``.reg`` region files for
        every chip of every input image and should
        ``generate_astrometric_catalog()`` generate file ``'refcatalog.cat'``?

    num_sources : int, optional
        Maximum number of **brightest sources per chip** which will be used
        for cross-matching and fitting. If set to `None`, all sources will
        be used.


    Updates
    -------
    fitab: astropy.table.Table
        Table which contains processing information and alignment results for
        every raw image evaluated.

    """
    log.info(
        "*** HLAPIPELINE Processing Version {!s} ({!s}) started at: {!s} ***\n"
        .format(__version__, __version_date__, util._ptime()[0])
    )
    fitab = Table()

    # Define astrometric catalog list in priority order
    catalogs = ['GAIADR2', 'GAIADR1']
    assert len(set(catalogs)) == len(catalogs)  # catalogues are unique

    # 0: print git info
    if print_git_info:
        _log_step('STEP 0: Display Git revision info')
        full_path = os.path.dirname(__file__)

        if "drizzlepac/drizzlepac" in full_path:
            repo_path = (
                full_path.split("drizzlepac/drizzlepac")[0] + "drizzlepac"
            )

        elif "hlapipeline" in full_path:
            repo_path = full_path.split("drizzlepac")[0] + "drizzlepac"

        else:
            repo_path = None

        if repo_path and os.path.exists(repo_path):
            # Display git repository information
            get_git_rev_info.print_rev_id(repo_path)
        else:
            log.warning("WARNING: Unable to display Git repository revision "
                        "information.")

    log.info(input_list)

    def _log_run_time(start_time, step):
        current_time = datetime.datetime.now()
        runtime = (current_time - start_time).total_seconds()
        log.info('Processing time of [STEP {}]: {} sec'.format(step, runtime))
        return current_time

    try:
        # 1: Interpret input data and optional parameters
        _log_step('STEP 1: Get data')

        begin_time = start_time = datetime.datetime.now()
        log.info(str(start_time))
        imglist = check_and_get_data(input_list, archive=archive,
                                     clobber=clobber)
        log.info("SUCCESS")

        start_time = _log_run_time(start_time, '1')

        # 2: Apply filter to input observations to insure that they meet
        # minimum criteria for being able to be aligned
        _log_step('STEP 2: Filter data')

        fitab = analyze.analyze_data(imglist)

        # Check the table to determine if there is any viable data to be
        # aligned. The 'doProcess' column (bool) indicates the image/file
        # should or should not be used for alignment purposes. For filtered
        # data, 'doProcess=0' and 'status=9999' in the table (the status value
        # by default), so there is no need to update the fitab here.
        if fitab['doProcess'].sum() == 0:
            log.warning("No viable images in filtered table - no processing "
                        "done.\n")
            _log_run_time(start_time, '2')
            return fitab

        # Get the list of all "good" files to use for the alignment
        process_list = list(
            fitab['imageName'][np.where(fitab['doProcess'])]
        )  # Convert process_list from numpy list to regular python list
        log.info("SUCCESS")

        start_time = _log_run_time(start_time, '2')

        # 3: Build WCS for full set of input observations
        _log_step('STEP 3: Build WCS')
        log.info("SUCCESS")
        start_time = _log_run_time(start_time, '2')

        # 4: Extract catalog of observable sources from each input image
        _log_step('STEP 4: Source finding')

        if debug:
            pickle_filename = ("{}.source_catalog.pickle"
                               .format(process_list[0]))
            if os.path.exists(pickle_filename):
                pickle_in = open(pickle_filename, "rb")
                extracted_sources = pickle.load(pickle_in)
                log.info("Using sourcelist extracted from {} generated during "
                         "the last run to save time.".format(pickle_filename))

            else:
                extracted_sources = generate_source_catalogs(
                    process_list, centering_mode='starfind',
                    nlargest=num_sources, output=output
                )
                pickle_out = open(pickle_filename, "wb")
                pickle.dump(extracted_sources, pickle_out)
                pickle_out.close()
                log.info("Wrote {}".format(pickle_filename))

        else:
            extracted_sources = generate_source_catalogs(
                process_list, centering_mode='starfind',
                nlargest=num_sources, output=output
            )

        for imgname in extracted_sources.keys():
            table = extracted_sources[imgname]["catalog_table"]

            # Get the location of the current image in the filtered table
            index = np.where(fitab['imageName'] == imgname)[0][0]

            # First ensure sources were found
            if table is None or not table[1]:
                log.warning("No sources found in image {}".format(imgname))
                fitab[:]['status'] = 1
                fitab[:]['processMsg'] = "No sources found"
                _log_run_time(start_time, '4')
                return fitab

            # The catalog of observable sources must have at least
            # _MIN_OBSERVABLE_THRESHOLD entries to be useful:
            total_num_sources = sum(map(len, table.values()))

            # Update filtered table with number of found sources
            fitab[index]['foundSources'] = total_num_sources

            if total_num_sources < _MIN_OBSERVABLE_THRESHOLD:
                log.warning("Not enough sources ({}) found in image {}"
                            .format(total_num_sources, imgname))
                fitab[:]['status'] = 1
                fitab[:]['processMsg'] = "Not enough sources found"
                _log_run_time(start_time, '4')
                return fitab

        log.info("SUCCESS")
        start_time = _log_run_time(start_time, '4')

        # 5: Retrieve list of astrometric sources from database
        # Convert input images to tweakwcs-compatible FITSWCS objects and
        # attach source catalogs to them.
        imglist = []
        for group_id, image in enumerate(process_list):
            img = build_wcscat(
                image, group_id, extracted_sources[image]['catalog_table']
            )
            # add the name of the image to the imglist object
            for im in img:
                log.info('im.meta[name] = {}'.format(im.meta['name']))
            imglist.extend(img)

        # store mapping of group_id to filename/chip
        group_id_dict = {}
        for image in imglist:
            k = "{}_{}".format(image.meta["filename"], image.meta["chip"])
            group_id_dict[k] = image.meta["group_id"]

        best_fit_rms = -99999.0
        best_fit_status_dict = {}
        best_fit_qual = 5
        # create pristine copy of imglist that will be used to restore imglist
        # back so it always starts exactly the same for each run.
        orig_imglist = copy.deepcopy(imglist)

        # Define fitting algorithm list in priority order.
        # The match_relative_fit algorithm must have more than one image
        # as the first image is the reference for the remaining images.
        fit_algorithms = [match_2dhist_fit, match_default_fit]
        if len(process_list) > 1:
            fit_algorithms.insert(0, match_relative_fit)

        # create dummy list that will be used to preserve imglist best_meta
        # information through the imglist reset process
        temp_imglist = []
        fit_info_dict = OrderedDict()
        reference_catalog_dict = {}

        for algorithm, catalog in product(fit_algorithms, catalogs):
            _log_step('STEP 5: Detect astrometric sources')
            log.info("Astrometric Catalog: {}".format(catalog))

            # store reference catalogs in a dictionary so that
            # generate_astrometric_catalog() doesn't execute unnecessarily
            # after it's been run once for a given astrometric catalog.
            if catalog in reference_catalog_dict:
                log.info(
                    "Using {} reference catalog from earlier this run."
                    .format(catalog)
                )
                reference_catalog = reference_catalog_dict[catalog]

            else:
                log.info(
                    "Generating new reference catalog for {}; Storing it "
                    "for potential re-use later this run."
                    .format(catalog)
                )
                reference_catalog = generate_astrometric_catalog(
                    process_list, catalog=catalog, output=output
                )
                reference_catalog_dict[catalog] = reference_catalog

            start_time = _log_run_time(start_time, '5')

            if len(reference_catalog) < _MIN_CATALOG_THRESHOLD:
                log.warning(
                    "Not enough sources found in catalog {}"
                    .format(catalog)
                )
                fit_quality = 5

                if catalog == catalogs[-1]:
                    log.info("Try again with other catalog")
                    continue

                else:
                    # bail out if not enough sources can be found any of
                    # the astrometric catalogs
                    log.warning("ERROR! No astrometric sources found in "
                                "any catalog. Exiting...")
                    fitab['status'][:] = 1
                    fitab['processMsg'][:] = "No astrometric sources found"
                    fitab['fit_qual'][:] = fit_quality
                    _log_run_time(start_time, '5')
                    return fitab

            else:
                _log_step('STEP 5b: Cross matching and fitting')
                # reset imglist to pristine state
                imglist = copy.deepcopy(orig_imglist)

                if temp_imglist:
                    # migrate best_meta to new imglist
                    for temp_item, item in zip(temp_imglist, imglist):
                        item.best_meta = temp_item.best_meta.copy()

                log.info(
                    "{} Catalog {} matched using {} {}"
                    .format('-' * 18, catalog,
                            algorithm.__name__, '-' * 18)
                )

                # restore group IDs to their pristine state prior to
                # each run.
                for image in imglist:
                    image.meta["group_id"] = group_id_dict[
                        "{}_{}".format(
                            image.meta["filename"], image.meta["chip"]
                        )
                    ]

                try:
                    # execute the correct fitting/matching algorithm
                    imglist = algorithm(imglist, reference_catalog)
                except Exception:
                    exc_type, exc_value, exc_tb = sys.exc_info()
                    traceback.print_exception(
                        exc_type, exc_value, exc_tb, file=sys.stdout
                    )
                    log.warning(
                        "WARNING: Catastrophic fitting failure with "
                        "catalog {} and matching algorithm {}."
                        .format(catalog, algorithm.__name__)
                    )
                    fitab['status'][:] = 1
                    fitab['processMsg'][:] = "Fitting failure"
                    # It may be there are additional catalogs and
                    # algorithms to try, so keep going but flag this fit
                    # with the 'bad' quality value
                    fit_quality = 5
                    fitab['fit_qual'][:] = fit_quality
                    continue

                # determine the quality of the fit
                fit_rms, fit_num, fit_quality, fitab, fit_status_dict = \
                    determine_fit_quality(
                        imglist, fitab, catalog == catalogs[-1],
                        print_fit_parameters=print_fit_parameters
                    )

                # save fit algorithm name to dictionary key
                # "fit method" in imglist.
                for img in imglist:
                    img.meta['fit method'] = algorithm.__name__
                    img.meta['fit quality'] = fit_quality

                # populate fit_info_dict
                kwd = "{} {}".format(
                    catalog, algorithm.__name__
                )
                fit_info_dict[kwd] = fit_status_dict[
                    next(iter(fit_status_dict))
                ]
                fit_info_dict[kwd]['fit_qual'] = fit_quality

                # Figure out which fit solution to go with based on
                # fit_quality value and maybe also total_rms
                if fit_quality >= 5:
                    pass

                elif fit_quality == 1:
                    # valid, non-comprimised solution with total
                    # rms < 10 mas...go with this solution.
                    best_fit_rms = fit_rms
                    for item in imglist:
                        item.best_meta = item.meta.copy()
                    best_fit_status_dict = fit_status_dict.copy()
                    break

                elif fit_quality < best_fit_qual:
                    # better solution found. keep looping but with
                    # the better solution as "best" for now.
                    log.info("Better solution found!")
                    best_fit_rms = fit_rms
                    for item in imglist:
                        item.best_meta = item.meta.copy()
                    best_fit_status_dict = fit_status_dict.copy()
                    best_fit_qual = fit_quality
                    # preserve best fit solution to be inserted into
                    # a reinitialized imglist next time through.
                    temp_imglist = copy.deepcopy(imglist)

                elif fit_quality == best_fit_qual:
                    # new solution same level of fit_quality.
                    # Choose whichever one has the lowest total
                    # rms as "best" and keep looping.
                    if fit_rms < best_fit_rms:
                        best_fit_rms = fit_rms
                        for item in imglist:
                            item.best_meta = item.meta.copy()
                        best_fit_status_dict = fit_status_dict.copy()
                    # preserve best fit solution to be inserted into
                    # a reinitialized imglist next time through.
                    temp_imglist = copy.deepcopy(imglist)

                else:
                    # new solution has worse fit_quality. discard
                    continue

            if fit_quality == 1:
                break

        start_time = _log_run_time(start_time, '5b')

        # 6: Populate the fitab
        _log_step('STEP 6: Collect up information and populate the filtered '
                  'table')

        if 0 < best_fit_rms < _MAX_FIT_RMS:
            log.info("The fitting process was successful with a best fit "
                     "total rms of {} mas".format(best_fit_rms))
        else:
            log.info("The fitting process was unsuccessful with a best fit "
                     "total rms of {} mas".format(best_fit_rms))

        if 0 < best_fit_rms < _MAX_FIT_LIMIT:
            # update to the meta information with the lowest rms if it is
            # reasonable
            for item in imglist:
                item.meta.update(item.best_meta)
            fitab['status'][:] = 0
            fit_status_dict = best_fit_status_dict .copy()

            # Protect the writing of the table within the best_fit_rms
            info_keys = OrderedDict(imglist[0].meta['fit_info']).keys()

            # Update filtered table with number of matched sources and other
            # information
            for item in imglist:
                imgname = item.meta['name']
                index = np.where(fitab['imageName'] == imgname)[0][0]

                fi = item.meta['fit_info']

                if fi['status'].startswith("FAILED"):
                    continue

                for tweakwcs_info_key in info_keys:
                    if tweakwcs_info_key.lower() == 'rms':
                        rx, ry = fi[tweakwcs_info_key]
                        fitab[index]['rms_x'] = rx
                        fitab[index]['rms_y'] = ry

                fitab[index]['fit_method'] = item.meta['fit method']
                fitab[index]['catalog'] = fi['catalog']
                fitab[index]['catalogSources'] = len(reference_catalog)
                fitab[index]['matchSources'] = fi['nmatches']
                fitab[index]['rms_ra'] = fi['RMS_RA'].value
                fitab[index]['rms_dec'] = fi['RMS_DEC'].value
                fitab[index]['fit_rms'] = fi['FIT_RMS']
                fitab[index]['total_rms'] = fi['TOTAL_RMS']
                shx, shy = fi['shift']
                fitab[index]['offset_x'] = shx
                fitab[index]['offset_y'] = shy
                fitab[index]['scale'] = fi['scale'][0]
                fitab[index]['rotation'] = fi['rot']

                # populate fitab fields "status", "compromised" and
                # "processMsg" with fit_status_dict fields "valid",
                # "compromised" and "reason".
                explicit_dict_key = "{},{}".format(
                    item.meta['name'], item.meta['chip']
                )
                fitab[index]['status'] = int(
                    not fit_status_dict[explicit_dict_key]['valid']
                )

                fitab['compromised'] = int(
                    fit_status_dict[explicit_dict_key]['compromised']
                )

                fitab[index]['processMsg'] = \
                    fit_status_dict[explicit_dict_key]['reason']
                fitab['fit_qual'][index] = item.meta['fit quality']

        start_time = _log_run_time(start_time, '6')

        # 7: Write new fit solution to input image headers
        _log_step('STEP 7: Update image headers with new WCS information')

        if 0 < best_fit_rms < 9999. and update_hdr_wcs:
            headerlet_dict = update_image_wcs_info(imglist)
            for tab in fitab:
                tab['headerletFile'] = headerlet_dict[tab['imageName']]
            log.info("SUCCESS")
        else:
            log.info(" STEP SKIPPED")

        current_time = _log_run_time(start_time, '7')
        log.info(
            'TOTAL Processing time of {} sec'
            .format((current_time - begin_time).total_seconds())
        )
        log.info(best_fit_status_dict)
        log.info(104 * '-')
        log.info(104 * '-')
        log.info("SUMMARY OF ALL FIT ATTEMPTS".rjust(56))
        for item in fit_info_dict.keys():
            log.info("{} {}".format(item, fit_info_dict[item]))
        log.info(104 * '-')

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)

    finally:
        fitab.pprint(max_width=-1)

    return fitab


def match_relative_fit(imglist, reference_catalog):
    """Perform cross-matching and final fit using relative matching algorithm

    Parameters
    ----------
    imglist : list
        List of input image `~tweakwcs.tpwcs.FITSWCS` objects with metadata
        and source catalogs

    reference_catalog : astropy.table.Table
        `~astropy.table.Table` of reference sources for this field

    Returns
    --------
    imglist : list
        List of input image `~tweakwcs.tpwcs.FITSWCS` objects with metadata
        and source catalogs.

    """
    _log_step('STEP 5b: (match_relative_fit) Cross matching and fitting')
    # 0: Specify matching algorithm to use
    match = tweakwcs.TPMatch(searchrad=75, separation=0.1, tolerance=2,
                             use2dhist=True)

    # Align images and correct WCS
    # NOTE: this invocation does not use an astrometric catalog. This call
    # allows all the input images to be aligned in a relative way using the
    # first input image as the reference.
    # 1: Perform relative alignment
    tweakwcs.align_wcs(imglist, None, match=match, expand_refcat=True)

    # Set all the group_id values to be the same so the various images/chips
    # will be aligned to the astrometric reference catalog as an ensemble.
    # astrometric reference catalog as an ensemble.
    #
    # BEWARE: If additional iterations of solutions are to be done,
    # the group_id values need to be restored.
    for image in imglist:
        image.meta["group_id"] = 1234567

    # 2: Perform absolute alignment
    tweakwcs.align_wcs(imglist, reference_catalog, match=match)

    # 3: Interpret RMS values from tweakwcs
    interpret_fit_rms(imglist, reference_catalog)

    return imglist


def match_default_fit(imglist, reference_catalog):
    """Perform cross-matching and final fit using default tolerance matching

    Parameters
    ----------
    imglist : list
        List of input image `~tweakwcs.tpwcs.FITSWCS` objects with metadata
        and source catalogs

    reference_catalog : astropy.table.Table
        `~astropy.table.Table` of reference sources for this field.

    Returns
    --------
    imglist : list
        List of input image `~tweakwcs.tpwcs.FITSWCS` objects with metadata
        and source catalogs

    """
    _log_step('STEP 5b: (match_default_fit) Cross matching and fitting')
    # Specify matching algorithm to use
    match = tweakwcs.TPMatch(searchrad=250, separation=0.1, tolerance=100,
                             use2dhist=False)

    # Align images and correct WCS
    tweakwcs.align_wcs(imglist, reference_catalog, match=match,
                       expand_refcat=False)

    # Interpret RMS values from tweakwcs
    interpret_fit_rms(imglist, reference_catalog)

    return imglist


def match_2dhist_fit(imglist, reference_catalog):
    """Perform cross-matching and final fit using 2D histogram matching.

    Parameters
    ----------
    imglist : list
        List of input image `~tweakwcs.tpwcs.FITSWCS` objects with metadata
        and source catalogs

    reference_catalog : astropy.table.Table
        `astropy.table.Table` of reference sources for this field.

    Returns
    --------
    imglist : list
        List of input image `~tweakwcs.tpwcs.FITSWCS` objects with metadata and
        source catalogs.

    """
    _log_step('STEP 5b: (match_2dhist_fit) Cross matching and fitting')

    # Specify matching algorithm to use
    match = tweakwcs.TPMatch(searchrad=75, separation=0.1, tolerance=2.0,
                             use2dhist=True)

    # Align images and correct WCS
    tweakwcs.align_wcs(imglist, reference_catalog, match=match,
                       expand_refcat=False)

    # Interpret RMS values from tweakwcs
    interpret_fit_rms(imglist, reference_catalog)

    return imglist


def determine_fit_quality(imglist, fitab, catalogs_remaining,
                          print_fit_parameters=True):
    """Determine the quality of the fit to the data

    Parameters
    ----------
    imglist: list
        output of ``interpret_fits``. Contains sourcelist tables,
        newly computed WCS info, etc. for every chip of every valid input
        image. This list should have been  updated, in-place, with the new
        RMS values; specifically,

            * ``'FIT_RMS'``: RMS of the separations between fitted image
              positions and reference positions;
            * ``'TOTAL_RMS'``: mean of the FIT_RMS values for all observations;
            * ``'NUM_FITS'``: number of images/group_id's with successful fits
              included in the ``'TOTAL_RMS'``;

        These entries are added to the ``'fit_info'`` dictionary.

    fitab: astropy.table.Table
        `~astropy.table.Table` object containing data pertaining to the
        associated dataset, including the ``doProcess`` bool. It is intended
        this table is updated by subsequent functions for bookkeeping purposes.

    catalogs_remaining: bool
        Specify whether additional catalogs remain to be fit against.

    print_fit_parameters: bool, optional
        Specify whether or not to print out FIT results for each chip.

    Returns
    -------
    max_rms_val: float
        The best Total rms determined from all of the images

    num_xmatches: int
        The number of stars used in matching the data

    fit_quality: int
        fit quality category:
            * 1 = valid solution with rms < 10 mas
            * 2 = Valid but compromised solution with rms < 10 mas
            * 3 = Valid solution with RMS >= 10 mas
            * 4 = Valid but compromised solution with RMS >= 10 mas
            * 5 = Not valid solution

    fitab: astropy.table.Table
        modified fitab object

    fit_status_dict: dictionary
        Dictionary containing the following:
            * overall fit validity (Boolean)
            * total (visit-level) RMS value in mas (float)
            * number of matched sources (int)
            * fit compromised status (Boolean)
            * reason fit is considered 'compromised' (only populated if
              "compromised" field is "True")
    """
    tweakwcs_info_keys = OrderedDict(imglist[0].meta['fit_info']).keys()
    max_rms_val = 1e9
    num_xmatches = 0
    fit_status_dict = {}
    overall_valid = True
    overall_comp = False

    xsh = []
    ysh = []
    for item in imglist:
        if not item.meta['fit_info']['status'].startswith('FAILED'):
            xsh.append(item.meta['fit_info']['shift'][0])
            ysh.append(item.meta['fit_info']['shift'][1])

    for item in imglist:
        image_name = item.meta['name']
        chip_num = item.meta['chip']

        # Build fit_status_dict entry
        dict_key = "{},{}".format(image_name, chip_num)
        fit_status_dict[dict_key] = {
            'valid': False,
            'max_rms': max_rms_val,
            'num_matches': num_xmatches,
            'compromised': False,
            'reason': ''
        }

        fit_info = item.meta['fit_info']

        # Handle fitting failures (no matches found)
        if fit_info['status'].startswith("FAILED"):
            log.warning("No cross matches found in any catalog for {} "
                        "- no processing done.".format(image_name))
            continue

        fit_rms_val = fit_info['FIT_RMS']
        max_rms_val = fit_info['TOTAL_RMS']
        num_xmatches = fit_info['nmatches']
        fit_status_dict[dict_key]['max_rms'] = max_rms_val
        fit_status_dict[dict_key]['num_matches'] = num_xmatches

        if num_xmatches < _MIN_CROSS_MATCHES and catalogs_remaining:
            log.warning(
                "Not enough cross matches found between astrometric"
                "catalog and sources found in {}".format(image_name)
            )
            continue

        # Execute checks
        nmatches_check = num_xmatches > 4
        rshift = np.linalg.norm(fit_info['shift']) * item.wcs.pscale
        radial_offset_check = 0.36 * num_xmatches > 0.8 + (rshift / 10.0)**8
        large_rms_check = fit_rms_val <= 150. and max_rms_val <= 150.

        rms_limit = max(fit_info['TOTAL_RMS'], 10.) / 1000.0 / item.wcs.pscale
        consistency_check = np.sqrt(np.var(xsh) + np.var(ysh)) <= rms_limit

        # Decide if fit solutions are valid based on checks
        if not consistency_check:
            # Failed consistency check
            fit_status_dict[dict_key]['valid'] = False
            fit_status_dict[dict_key]['compromised'] = False
            fit_status_dict[dict_key]['reason'] = "Consistency violation!"

        elif not large_rms_check:
            # RMS value(s) too large
            fit_status_dict[dict_key]['valid'] = False
            fit_status_dict[dict_key]['compromised'] = False
            fit_status_dict[dict_key]['reason'] = "RMS too large (>150 mas)!"

        elif not radial_offset_check:
            # Failed radial offset check
            fit_status_dict[dict_key]['valid'] = False
            fit_status_dict[dict_key]['compromised'] = False
            fit_status_dict[dict_key]['reason'] = "Radial offset too large!"

        elif not nmatches_check:
            # Too few matches
            fit_status_dict[dict_key]['valid'] = True
            fit_status_dict[dict_key]['compromised'] = True
            fit_status_dict[dict_key]['reason'] = "Too few matches!"

        else:
            # all checks passed. Valid solution.
            fit_status_dict[dict_key]['valid'] = True
            fit_status_dict[dict_key]['compromised'] = False
            fit_status_dict[dict_key]['reason'] = ""

        # for now, generate overall valid and compromised values. Basically,
        # if any of the entries for "valid" is False, "valid" is False,
        # treat the whole dataset as not valid. Same goes for compromised.
        if not fit_status_dict[dict_key]['valid']:
            overall_valid = False

        if fit_status_dict[dict_key]['compromised']:
            overall_comp = True

        log.info(
            'RESULTS FOR {} Chip {}: FIT_RMS = {} mas, TOTAL_RMS = {}'
            ' mas, NUM =  {}'.format(
                image_name, item.meta['chip'], fit_rms_val, max_rms_val,
                num_xmatches
            )
        )

        if print_fit_parameters:
            log.info("{} FIT PARAMETERS {}".format(35 * '~', 34 * '~'))
            log.info("image: {}".format(image_name))
            log.info("chip: {}".format(item.meta['chip']))
            log.info("group_id: {}".format(item.meta['group_id']))
            for tweakwcs_info_key in tweakwcs_info_keys:
                if not tweakwcs_info_key.startswith("matched"):
                    log.info("{} : {}".format(tweakwcs_info_key,
                                              fit_info[tweakwcs_info_key]))
            log.info("~" * 84)
            log.info(
                "nmatches_check: {} radial_offset_check: {}"
                " large_rms_check: {}, consistency_check: {}"
                .format(nmatches_check, radial_offset_check, large_rms_check,
                        consistency_check)
            )

    # determine which fit quality category this latest fit falls into
    if overall_valid is False:
        fit_quality = 5
        log.info("FIT SOLUTION REJECTED")
        fitab['status'][:] = 1
        for ctr in range(0, len(fitab)):
            fitab[ctr]['processMsg'] = fit_status_dict[
                fitab[ctr]['imageName'] + ',1']['reason']
    else:
        for ctr in range(0, len(fitab)):
            fitab[ctr]['processMsg'] = ''

        if overall_comp is False and max_rms_val < 10.:
            log.info("Valid solution with RMS < 10 mas found!")
            fit_quality = 1

        elif overall_comp is True and max_rms_val < 10.:
            log.info("Valid but compromised solution with RMS < 10 mas found!")
            fit_quality = 2

        elif overall_comp is False and max_rms_val >= 10.:
            log.info("Valid solution with RMS >= 10 mas found!")
            fit_quality = 3

        else:
            log.info("Valid but compromised solution with RMS >= 10 "
                     "mas found!")
            fit_quality = 4

    if print_fit_parameters:
        for item in imglist:
            log.info(fit_status_dict["{},{}".format(item.meta['name'],
                                                    item.meta['chip'])])

    if max_rms_val > _MAX_FIT_RMS:
        log.info("Total fit RMS value = {} mas greater than the maximum "
                 "threshold value {}.".format(max_rms_val, _MAX_FIT_RMS))

    if not overall_valid:
        log.info("The fit solution for some or all of the images is not "
                 "valid.")

    if max_rms_val > _MAX_FIT_RMS or not overall_valid:
        log.info("Try again with the next catalog")
    else:
        log.info("Fit calculations successful.")

    return max_rms_val, num_xmatches, fit_quality, fitab, fit_status_dict


def generate_astrometric_catalog(imglist, **pars):
    """Generates a catalog of all sources from an existing astrometric
    catalog are in or near the FOVs of the images in the input list.

    Parameters
    ----------
    imglist : list
        List of one or more calibrated fits images that will be used for
        catalog generation.

    Returns
    =======
    ref_table : astropy.table.Table
        Catalog as a `~astropy.table.Table`

    """
    temp_pars = pars.copy()
    pars['output'] = 'ref_cat.ecsv' if pars['output'] else None
    out_catalog = create_astrometric_catalog(imglist, **pars)
    pars = temp_pars.copy()

    if out_catalog and pars['output']:
        catalog_filename = "refcatalog.cat"
        out_catalog.write(catalog_filename,
                          format="ascii.fast_commented_header")
        log.info("Wrote reference catalog {}.".format(catalog_filename))

    return(out_catalog)


def generate_source_catalogs(imglist, **pars):
    """Generates a dictionary of source catalogs keyed by image name.

    Parameters
    ----------
    imglist : list
        List of one or more calibrated fits images that will be used for
        source detection.

    Returns
    -------
    src_cat : dictionary
        a dictionary (keyed by image name) of two-element dictionaries which
        contain the following:
            * a dictionary of the detector-specific processing parameters;
            * an astropy table of position and photometry information of
              all detected sources.
    """
    output = pars.get('output', False)
    src_cat = {}
    for imgname in imglist:
        log.info("Image name: {}".format(imgname))

        src_cat[imgname] = {}

        # open image
        with fits.open(imgname) as f:
            phdr = f[0].header
            instrument = phdr['INSTRUME'].lower()
            detector = phdr['DETECTOR'].lower()

            # get instrument/detector-specific image alignment parameters
            try:
                detector_info = _DETECTOR_SPECIFIC_PARAMS[instrument]
                try:
                    detector_pars = detector_info[detector]
                    # to allow generate_source_catalog to get detector
                    # specific parameters
                    detector_pars.update(pars)
                    src_cat[imgname]["params"] = detector_pars

                except KeyError:
                    log.exit(
                        "ERROR! Unrecognized detector '{}'. Exiting..."
                        .format(detector)
                    )
                    os.error(
                        "ERROR! Unrecognized detector '{}'. Exiting..."
                        .format(detector)
                    )

            except KeyError:
                log.exit("ERROR! Unrecognized instrument '{}'. Exiting..."
                         .format(instrument))
                os.error("ERROR! Unrecognized instrument '{}'. Exiting..."
                         .format(instrument))

            # Identify sources in image, convert coords from chip x, y form
            # to reference WCS sky RA, Dec form.
            imgwcs = HSTWCS(f, 1)

            # Convert fwhmpsf from arsec to pixels
            fwhmpsf_pix = src_cat[imgname]["params"]['fwhmpsf'] / imgwcs.pscale
            src_cat[imgname]["catalog_table"] = generate_source_catalog(
                f, fwhm=fwhmpsf_pix, **detector_pars
            )

            # write out coord lists to files for diagnostic purposes. Protip:
            # To display the sources in these files in DS9, set the
            # "Coordinate System" option to "Physical" when loading the region
            # file.
            imgroot = os.path.basename(imgname).split('_')[0]
            num_sci = countExtn(f)
            if output:
                for chip in range(1, num_sci + 1):
                    chip_cat = src_cat[imgname]["catalog_table"][chip]
                    if chip_cat:
                        regfilename = "{}_sci{}_src.reg".format(imgroot, chip)
                        out_table = Table(chip_cat)

                        # To align with positions of sources in DS9/IRAF
                        out_table['xcentroid'] += 1
                        out_table['ycentroid'] += 1
                        out_table.write(
                            regfilename,
                            include_names=["xcentroid", "ycentroid"],
                            format="ascii.fast_commented_header"
                        )
                        log.info("Wrote region file {}\n".format(regfilename))

    return src_cat


def update_image_wcs_info(tweakwcs_output):
    """Write newly computed WCS information to image headers and write
    headerlet files.

    Parameters
    ----------
    tweakwcs_output : list
        output of tweakwcs. Contains sourcelist tables, newly computed WCS
        info, etc. for every chip of every valid every valid input image.

    Returns
    -------
    out_headerlet_list : dictionary
        a dictionary of the headerlet files created by this subroutine,
        keyed by ``flt/flc`` FITS filename.

    """
    out_headerlet_dict = {}

    # group "tweaked_output" by file name. Then we can safely open/close
    # that file and process all its chips at once.
    g_tweak_output = OrderedDict()
    for item in tweakwcs_output:
        fname = item.meta['filename']
        if fname not in g_tweak_output:
            g_tweak_output[fname] = []
        g_tweak_output[fname].append(item)

    # sort chips:
    for fname, chip_data in g_tweak_output.items():
        g_tweak_output[fname] = sorted(chip_data, lambda x: x.meta['chip'])

    for fname, chip_data in g_tweak_output.items():
        assert fname == chip_data[0].meta['filename']
        with fits.open(fname, mode='update') as hdulist:
            meta0 = chip_data[0].meta

            # generate wcs name for updated image header, headerlet
            # Just in case header value 'wcs_name' is empty.
            if item.meta['fit method'] == 'match_relative_fit':
                fit_method = 'REL'
            else:
                fit_method = 'IMG'

            if hdulist['SCI', 1].header['WCSNAME']:
                wcs_name = "FIT_{}_{}".format(fit_method, meta0['catalog'])

            else:
                wname = hdulist['sci', 1].header['wcsname']
                if "-" in wname:
                    wcs_name = '{}-FIT_{}_{}'.format(
                        wname[:wname.index('-')], fit_method,
                        meta0['fit_info']['catalog']
                    )
                else:
                    wcs_name = '{}-FIT_{}_{}'.format(
                        wname, fit_method, meta0['fit_info']['catalog']
                    )

            for item in chip_data:
                sci_extn = fileutil.findExtname(
                    hdulist, 'sci', extver=item.meta['chip']
                )

                # update header with new WCS info
                updatehdr.update_wcs(
                    hdulist, sci_extn, item.wcs, wcsname=wcs_name,
                    reusename=True, verbose=True
                )

        # Create headerlet
        out_headerlet = headerlet.create_headerlet(
            fname, hdrname=wcs_name, wcsname=wcs_name
        )

        # Update headerlet
        update_headerlet_phdu(item, out_headerlet)

        # Write headerlet
        if fname.endswith('flc.fits'):
            headerlet_filename = fname.replace('flc', 'flt_hlet')

        elif fname.endswith('flt.fits'):
            headerlet_filename = fname.replace('flt', 'flt_hlet')

        out_headerlet.writeto(headerlet_filename, clobber=True)
        log.info("Wrote headerlet file {}.\n\n".format(headerlet_filename))
        out_headerlet_dict[fname] = headerlet_filename

        # Attach headerlet as HDRLET extension
        headerlet.attach_headerlet(fname, headerlet_filename)

    return (out_headerlet_dict)


def update_headerlet_phdu(tpwcs, headerlet):
    """Update the primary header data unit keywords of a headerlet
    object in-place

    Parameters
    ----------
    tpwcs: tweakwcs.tpwcs.TPWCS
        Basically the output from ``tweakwcs`` which contains the cross
        match and fit information as well as an aligned WCS.

    headerlet: headerlet object
        Object containing WCS information

    """
    # Get the data to be used as values for FITS keywords
    rms_ra = tpwcs.meta['fit_info']['RMS_RA'].value
    rms_dec = tpwcs.meta['fit_info']['RMS_DEC'].value
    fit_rms = tpwcs.meta['fit_info']['FIT_RMS']
    nmatch = tpwcs.meta['fit_info']['nmatches']
    catalog = tpwcs.meta['fit_info']['catalog']
    fit_method = tpwcs.meta['fit method']

    x_shift, y_shift = tpwcs.meta['fit_info']['shift']
    rot = tpwcs.meta['fit_info']['rot']
    scale = tpwcs.meta['fit_info']['scale'][0]
    skew = tpwcs.meta['fit_info']['skew']

    # Update the existing FITS keywords
    primary_header = headerlet[0].header
    primary_header['RMS_RA'] = rms_ra
    primary_header['RMS_DEC'] = rms_dec
    primary_header['NMATCH'] = nmatch
    primary_header['CATALOG'] = catalog
    primary_header['FITMETH'] = fit_method

    # Create a new FITS keyword
    primary_header['FIT_RMS'] = (
        fit_rms, 'RMS (mas) of the 2D fit of the headerlet solution'
    )

    # Create the set of HISTORY keywords
    primary_header['HISTORY'] = '~~~~~ FIT PARAMETERS ~~~~~'
    primary_header['HISTORY'] = '{:>15} : {:9.4f} "/pixels'.format(
        'platescale', tpwcs.wcs.pscale
    )
    primary_header['HISTORY'] = '{:>15} : {:9.4f} pixels'.format(
        'x_shift', x_shift
    )
    primary_header['HISTORY'] = '{:>15} : {:9.4f} pixels'.format(
        'y_shift', y_shift
    )
    primary_header['HISTORY'] = '{:>15} : {:9.4f} degrees'.format('rotation',
                                                                  rot)
    primary_header['HISTORY'] = '{:>15} : {:9.4f}'.format('scale', scale)
    primary_header['HISTORY'] = '{:>15} : {:9.4f}'.format('skew', skew)


def _compute_fit_sky_rms(fit_info, refcat):
    ref_idx = fit_info['matched_ref_idx']
    fitmask = fit_info['fitmask']
    ref_RA = refcat[ref_idx]['RA'][fitmask]
    ref_DEC = refcat[ref_idx]['DEC'][fitmask]
    img_coords = SkyCoord(fit_info['fit_RA'], fit_info['fit_DEC'],
                          unit='deg', frame='icrs')
    ref_coords = SkyCoord(ref_RA, ref_DEC, unit='deg', frame='icrs')
    dra, ddec = img_coords.spherical_offsets_to(ref_coords)
    fit_rms = np.std(
        Angle(img_coords.separation(ref_coords), unit=u.mas)
    ).value
    rms_ra = np.std(dra.to(u.mas))
    rms_dec = np.std(ddec.to(u.mas))
    return fit_rms, rms_ra, rms_dec


def interpret_fit_rms(tpwcs_lst, reference_catalog):
    """Interpret the FIT information to convert RMS to physical units

    Parameters
    ----------
    tpwcs_lst : list of tweakwcs.tpwcs.TPWCS
        Each item in the list contains sourcelist tables, newly computed WCS
        info, etc. for every chip of every valid input image. This list gets
        updated, in-place, with the new RMS values; specifically,

            * ``'FIT_RMS'``: RMS of the separations between fitted image
              positions and reference positions
            * ``'TOTAL_RMS'``: mean of the FIT_RMS values for all observations
            * ``'NUM_FITS'``: number of images/group_id's with successful fits
              included in the ``'TOTAL_RMS'``

        These entries are added to the 'fit_info' dictionary.

    reference_catalog : astropy.table.Table
        Table of reference source positions used for the fit.

    """
    # Start by collecting information by group_id
    grouped_tpwcs = OrderedDict([(None, [])])
    for tpwcs in tpwcs_lst:
        group_id = tpwcs.meta['group_id']
        if group_id in grouped_tpwcs:
            grouped_tpwcs[group_id].append(tpwcs)
        else:
            grouped_tpwcs[group_id] = [tpwcs]
    non_grouped = grouped_tpwcs.pop(None)

    num_fits = 0
    total_rms = 0.0

    # iterate over "groups":
    for tp_items in grouped_tpwcs.values():
        tpwcs = tp_items[0]
        tinfo = tpwcs.meta['fit_info']
        if tinfo['status'] == 'SUCCESS':
            log.info("fit_info: {}".format(tinfo))
            fit_rms, ra_rms, dec_rms = _compute_fit_sky_rms(
                tinfo, reference_catalog
            )
            num_fits += 1
            total_rms += fit_rms

        else:
            # fit info is not complete when fit status is other than 'SUCCESS'
            fit_rms, ra_rms, dec_rms = None

        for tpwcs in tp_items:
            tpwcs.meta['fit_info']['FIT_RMS'] = fit_rms
            tpwcs.meta['fit_info']['RMS_RA'] = ra_rms
            tpwcs.meta['fit_info']['RMS_DEC'] = dec_rms

    # iterate over non-grouped/stand-alone WCS:
    for tpwcs in non_grouped:
        tinfo = tpwcs.meta['fit_info']
        if tinfo['status'] == 'SUCCESS':
            log.info("fit_info: {}".format(tinfo))
            fit_rms, ra_rms, dec_rms = _compute_fit_sky_rms(
                tinfo, reference_catalog
            )
            num_fits += 1
            total_rms += fit_rms

        else:
            # fit info is not complete when fit status is other than 'SUCCESS'
            fit_rms, ra_rms, dec_rms = None

        tpwcs.meta['fit_info']['FIT_RMS'] = fit_rms
        tpwcs.meta['fit_info']['RMS_RA'] = ra_rms
        tpwcs.meta['fit_info']['RMS_DEC'] = dec_rms

    # Compute RMS for entire ASN/observation set
    total_rms /= num_fits

    # Now, append computed results to tpwcs_lst
    for tpwcs in tpwcs_lst:
        tpwcs.meta['fit_info']['TOTAL_RMS'] = total_rms
        tpwcs.meta['fit_info']['NUM_FITS'] = num_fits
        tpwcs.meta['fit_info']['catalog'] = reference_catalog.meta['catalog']


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Align images')
    parser.add_argument('raw_input_list', nargs='+', help='The Images one '
                        'wishes to align. Valid input formats: '
                        '1. An association name; Example; j92c12345. '
                        '2. A space-separated list of flc.fits (or flt.fits) '
                        'files to align; Example: aaa_flc.fits bbb_flc.fits  '
                        'ccc_flc.fits 3. a simple text file containing a list '
                        'of fits files to align, one per line; '
                        'Example: input_list.txt')

    parser.add_argument('-a', '--archive', required=False, action='store_true',
                        help='Turning on this option will retain '
                             'copies of the downloaded files in the '
                             'astroquery created sub-directories.')

    parser.add_argument('-c', '--clobber', required=False, action='store_true',
                        help='If this option is turned on, the '
                             'program will download new copies of the input '
                             'files, overwriting any existing local copies in '
                             'the working directory')

    parser.add_argument('-d', '--debug', required=False, action='store_true',
                        help='If this option is turned on, the program will '
                             'attempt to use saved sourcelists stored in a '
                             'pickle file generated during a previous run. '
                             'Using a saved sorucelist instead of generating '
                             'new sourcelists greatly reduces overall run '
                             'time. If the pickle file does not exist, the '
                             'program will generate new sourcelists and save '
                             'them in a pickle file named after the first '
                             'input file.')

    parser.add_argument('-g', '--print_git_info', required=False,
                        action='store_true',
                        help='Turning on this option will '
                             'display git repository information at the start '
                             'of the run.')

    parser.add_argument('-o', '--output', required=False, action='store_true',
                        help='If turned on, '
                             'utils.astrometric_utils.create_astrometric_'
                             'catalog() generate file "ref_cat.ecsv", '
                             'generate_source_catalogs() generate the .reg '
                             'region files for every chip of every input '
                             'image and generate_astrometric_catalog() '
                             'generate file "refcatalog.cat".')

    parser.add_argument('-p', '--print_fit_parameters', required=False,
                        action='store_true', help='Turning on this option '
                               'will print out fit results for each chip.')

    parser.add_argument('-u', '--update_hdr_wcs', required=False,
                        action='store_true',
                        help='Turning on this option will write newly '
                             'computed WCS information to image image headers '
                             'and create headerlet files.')
    args = parser.parse_args()

    # Build list of input images
    input_list = []
    for item in args.raw_input_list:
        if os.path.exists(item):
            if item.endswith(".fits"):
                input_list.append(item)
            else:
                with open(item, 'r') as infile:
                    file_lines = infile.readlines()
                for file_line in file_lines:
                    input_list.append(file_line.strip())
        else:
            log.info("{} not found in working directory!".format(item))
            input_list.append(item)

    # Get to it!
    return_value = perform_align(
        input_list,
        archive=args.archive,
        clobber=args.clobber,
        debug=args.debug,
        update_hdr_wcs=args.update_hdr_wcs,
        print_fit_parameters=args.print_fit_parameters,
        print_git_info=args.print_git_info,
        output=args.output
    )
    log.info(return_value)
