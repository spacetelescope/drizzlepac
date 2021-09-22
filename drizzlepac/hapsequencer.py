#!/usr/bin/env python

""" This script defines the HST Advanced Products (HAP) generation portion of the
    calibration pipeline.  This portion of the pipeline produces mosaic and catalog
    products. This script provides similar functionality as compared to the Hubble
    Legacy Archive (HLA) pipeline in that it acts as the controller for the overall
    sequence of processing.

    Note regarding logging...
    During instantiation of the log, the logging level is set to NOTSET which essentially
    means all message levels (debug, info, etc.) will be passed from the logger to
    the underlying handlers, and the handlers will dispatch all messages to the associated
    streams.  When the command line option of setting the logging level is invoked, the
    logger basically filters which messages are passed on to the handlers according the
    level chosen. The logger is acting as a gate on the messages which are allowed to be
    passed to the handlers.

    Creation of source catalogs can be controlled through the use of environment variables:

      - SVM_CATALOG_HRC
      - SVM_CATALOG_SBC
      - SVM_CATALOG_WFC
      - SVM_CATALOG_UVIS
      - SVM_CATALOG_IR

    These variables can be defined using values of:

      - 'on', 'true', 'yes' : Create catalogs
      - 'off', 'false', 'no' : Turn off generation of catalogs

    The output products can be evaluated to determine the quality of the alignment and
    output data through the use of the environment variable:

    - SVM_QUALITY_TESTING : Turn on quality assessment processing.  This environment
      variable, if found with an affirmative value, will turn on processing to generate a JSON
      file which contains the results of evaluating the quality of the generated products.

    NOTE: Step 9 compares the output HAP products to the Hubble Legacy Archive (HLA)
    products. In order for step 9 (run_sourcelist_comparison()) to run, the following
    environment variables need to be set:
    - HLA_CLASSIC_BASEPATH
    - HLA_BUILD_VER

    Alternatively, if the HLA classic path is unavailable, The comparison can be run using
    locally stored HLA classic files. The relevant HLA classic imagery and sourcelist files
    must be placed in a subdirectory of the current working directory called 'hla_classic'.
"""
import datetime
import fnmatch
import logging
import os
import pickle
import shutil
import sys
import traceback

import numpy as np
from astropy.io import fits
from astropy.table import Table

import drizzlepac
from drizzlepac import util
from drizzlepac.haputils import config_utils
from drizzlepac.haputils import diagnostic_utils
from drizzlepac.haputils import hla_flag_filter
from drizzlepac.haputils import poller_utils
from drizzlepac.haputils import product
from drizzlepac.haputils import processing_utils as proc_utils
from drizzlepac.haputils import svm_quality_analysis as svm_qa
from drizzlepac.haputils.catalog_utils import HAPCatalogs

from stsci.tools.fileutil import countExtn
from stsci.tools import logutil
from stwcs import updatewcs
from stwcs import wcsutil

__taskname__ = 'hapsequencer'
MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)
__version__ = 0.1
__version_date__ = '07-Nov-2019'

# Environment variable which controls the quality assurance testing
# for the Single Visit Mosaic processing.
envvar_bool_dict = {'off': False, 'on': True, 'no': False, 'yes': True, 'false': False, 'true': True}
envvar_qa_svm = "SVM_QUALITY_TESTING"

envvar_cat_svm = {"SVM_CATALOG_SBC": 'on',
                  "SVM_CATALOG_HRC": 'on',
                  "SVM_CATALOG_WFC": 'on',
                  "SVM_CATALOG_UVIS": 'on',
                  "SVM_CATALOG_IR": 'on'}
envvar_cat_str = "SVM_CATALOG_{}"

# --------------------------------------------------------------------------------------------------------------


def create_catalog_products(total_obj_list, log_level, diagnostic_mode=False, phot_mode='both',
                            catalog_switches=None):
    """This subroutine utilizes haputils/catalog_utils module to produce photometric sourcelists for the specified
    total drizzle product and it's associated child filter products.

    PARAMETERS
    -----------
    total_obj_list : `drizzlepac.haputils.product.TotalProduct`
                    total drizzle product that will be processed by catalog_utils.
                    catalog_utils will also create photometric
                    sourcelists for the child filter products of this total product.
    log_level : int, optional
                The desired level of verboseness in the log statements displayed on the screen and written to the .log file.
    diagnostic_mode : bool, optional
                      generate ds9 region file counterparts to the photometric sourcelists? Default value is False.
    phot_mode : str, optional
                Which algorithm should be used to generate the sourcelists? 'aperture' for aperture photometry;
                'segment' for segment map photometry; 'both' for both 'segment' and 'aperture'. Default value is 'both'.
    catalog_switches : dict, optional
                       Specify which, if any, catalogs should be generated at all, based on detector.  This dictionary
                       needs to contain values for all instruments; namely:

                       SVM_CATALOG_HRC, SVM_CATALOG_SBC, SVM_CATALOG_WFC, SVM_CATALOG_UVIS, SVM_CATALOG_IR

                       These variables can be defined with values of 'on'/'off'/'yes'/'no'/'true'/'false'.

    Returns
    -------
    product_list : list
        list of all catalogs generated.
    """
    product_list = []
    log.info("Generating total product source catalogs")
    phot_mode = phot_mode.lower()
    input_phot_mode = phot_mode

    for total_product_obj in total_obj_list:
        # Need to have direct exposures for catalog generation
        if total_product_obj.edp_list:
            cat_sw_name = envvar_cat_str.format(total_product_obj.detector.upper())
            if catalog_switches[cat_sw_name] is False:
                log.info("Catalog generation turned OFF for {}".format(total_product_obj.detector.upper()))
                continue

            # Make sure this is re-initialized for the new total product
            phot_mode = input_phot_mode

            # Generate an "n" exposure mask which has the image footprint set to the number
            # of exposures which constitute each pixel.
            total_product_obj.generate_footprint_mask()

            # Instantiate filter catalog product object
            total_product_catalogs = HAPCatalogs(total_product_obj.drizzle_filename,
                                                 total_product_obj.configobj_pars.get_pars('catalog generation'),
                                                 total_product_obj.configobj_pars.get_pars('quality control'),
                                                 total_product_obj.mask,
                                                 log_level,
                                                 types=input_phot_mode,
                                                 diagnostic_mode=diagnostic_mode)

            # Identify sources in the input image and delay writing the total detection
            # catalog until the photometric measurements have been done on the filter
            # images and some of the measurements can be appended to the total catalog
            total_product_catalogs.identify(mask=total_product_obj.mask)

            # Build dictionary of total_product_catalogs.catalogs[*].sources to use for
            # filter photometric catalog generation
            sources_dict = {}
            filter_catalogs = {}
            source_mask = {}
            for cat_type in total_product_catalogs.catalogs.keys():
                sources_dict[cat_type] = {}
                sources_dict[cat_type]['sources'] = total_product_catalogs.catalogs[cat_type].sources
                # For the segmentation source finding, both the segmentation image AND the segmentation catalog
                # computed for the total object need to be provided to the filter objects
                if cat_type == "segment":
                    sources_dict['segment']['kernel'] = total_product_catalogs.catalogs['segment'].kernel
                    sources_dict['segment']['source_cat'] = total_product_catalogs.catalogs['segment'].source_cat

            # Get parameter from config files for CR rejection of catalogs
            cr_residual = total_product_obj.configobj_pars.get_pars('catalog generation')['cr_residual']
            n1_exposure_time = 0

            log.info("Generating filter product source catalogs")
            for filter_product_obj in total_product_obj.fdp_list:

                # Load a dictionary with the filter subset table for each catalog...
                subset_columns_dict = {}

                # Instantiate filter catalog product object
                filter_product_catalogs = HAPCatalogs(filter_product_obj.drizzle_filename,
                                                      total_product_obj.configobj_pars.get_pars('catalog generation'),
                                                      total_product_obj.configobj_pars.get_pars('quality control'),
                                                      total_product_obj.mask,
                                                      log_level,
                                                      types=phot_mode,
                                                      diagnostic_mode=diagnostic_mode,
                                                      tp_sources=sources_dict)

                flag_trim_value = filter_product_catalogs.param_dict['flag_trim_value']

                # Perform photometry
                # The measure method also copies a specified portion of the filter table into
                # a filter "subset" table which will be combined with the total detection table.
                filter_name = filter_product_obj.filters
                filter_product_catalogs.measure(filter_name)
                log.info("Flagging sources in filter product catalog")
                filter_product_catalogs = run_sourcelist_flagging(filter_product_obj,
                                                                  filter_product_catalogs,
                                                                  log_level,
                                                                  diagnostic_mode)

                if total_product_obj.detector.upper() not in ['IR', 'SBC']:
                    # Apply cosmic-ray threshold criteria used by HLA to determine whether or not to reject
                    # the catalogs.
                    tot_exposure_time = 0
                    n1_factor = 0.0
                    n1_dict = {}
                    for edp in total_product_obj.edp_list:
                        if edp.filters not in n1_dict:
                            n1_dict[edp.filters] = {'n': 1, 'texptime': edp.exptime}
                        else:
                            n1_dict[edp.filters]['n'] += 1
                            n1_dict[edp.filters]['texptime'] += edp.exptime

                    for edp in total_product_obj.edp_list:
                        tot_exposure_time += edp.exptime
                        if n1_dict[edp.filters]['n'] == 1:
                            n1_exposure_time += edp.exptime
                            n1_factor += cr_residual

                    # Insure that n1_factor only improves the threshold, not make it worse.
                    n1_factor = min(n1_factor, 1.0)

                    # Account for the influence of the single-image cosmic-ray identification
                    # This fraction represents the residual number of cosmic-rays after single-image identification
                    if n1_exposure_time < tot_exposure_time:
                        n1_exposure_time *= n1_factor

                # write out CI and FWHM values to file (if IRAFStarFinder was used instead of DAOStarFinder) for hla_flag_filter parameter optimization.
                if diagnostic_mode and phot_mode in ['aperture', 'both']:
                    if "fwhm" in total_product_catalogs.catalogs['aperture'].sources.colnames:
                        diag_obj = diagnostic_utils.HapDiagnostic(log_level=log_level)
                        diag_obj.instantiate_from_hap_obj(filter_product_obj,
                                                          data_source=__taskname__,
                                                          description="CI vs. FWHM values")
                        output_table = Table([filter_product_catalogs.catalogs['aperture'].source_cat['CI'], total_product_catalogs.catalogs['aperture'].sources['fwhm']], names=("CI", "FWHM"))

                        diag_obj.add_data_item(output_table, "CI_FWHM")
                        diag_obj.write_json_file(filter_product_obj.point_cat_filename.replace(".ecsv", "_ci_fwhm.json"))
                        del output_table
                        del diag_obj

                # Replace zero-value total-product catalog 'Flags' column values with meaningful filter-product catalog
                # 'Flags' column values
                for cat_type in filter_product_catalogs.catalogs.keys():
                    filter_product_catalogs.catalogs[cat_type].subset_filter_source_cat[
                       'Flags_{}'.format(filter_product_obj.filters)] = \
                       filter_product_catalogs.catalogs[cat_type].source_cat['Flags']
                    source_mask[cat_type] = None

                filter_catalogs[filter_product_obj.drizzle_filename] = filter_product_catalogs

                # Determine which rows should be removed from each type of catalog based on Flag values
                # Any source with Flag > 5 in any filter product catalog will be marked for removal from
                # all catalogs.
                # This requires collating results for each type of catalog from all filter products.
                for cat_type in filter_product_catalogs.catalogs.keys():
                    catalog_mask = filter_product_catalogs.catalogs[cat_type].source_cat['Flags'] > flag_trim_value
                    if source_mask[cat_type] is None:
                        source_mask[cat_type] = catalog_mask
                    else:
                        # Combine masks for all filters for this catalog type
                        source_mask[cat_type] = np.bitwise_or(source_mask[cat_type], catalog_mask)

                    # Trim based on user-specified/default flag limit 'flag_trim_value' specified in parameter file
                    trimmed_rows = np.where(source_mask[cat_type])[0].tolist()
                    filter_product_catalogs.catalogs[cat_type].source_cat.remove_rows(trimmed_rows)
                    filter_product_catalogs.catalogs[cat_type].subset_filter_source_cat.remove_rows(trimmed_rows)

                    subset_columns_dict[cat_type] = {}
                    subset_columns_dict[cat_type]['subset'] = \
                        filter_product_catalogs.catalogs[cat_type].subset_filter_source_cat

                # ...and append the filter columns to the total detection product catalog.
                log.info("create catalogs. combine() columns: {}".format(subset_columns_dict))
                total_product_catalogs.combine(subset_columns_dict)

            # At this point the total product catalog contains all columns contributed
            # by each filter catalog. However, some of the columns originating in one or more of
            # the filter catalogs contain no measurements for a particular source.  Remove all
            # rows which contain empty strings (masked values) for all measurements for all
            # of the filter catalogs.
            for cat_type in total_product_catalogs.catalogs.keys():
                good_rows_index = []
                if cat_type == 'aperture':
                    all_columns = total_product_catalogs.catalogs[cat_type].sources.colnames
                    table_filled = total_product_catalogs.catalogs[cat_type].sources.filled(-9999.9)
                else:
                    all_columns = total_product_catalogs.catalogs[cat_type].source_cat.colnames
                    table_filled = total_product_catalogs.catalogs[cat_type].source_cat.filled(-9999.9)
                flag_columns = [colname for colname in all_columns if "Flags_" in colname]
                filled_flag_columns = table_filled[flag_columns]
                for i, trow in enumerate(filled_flag_columns):
                    for tcol in trow:
                        if tcol != -9999:
                            good_rows_index.append(i)
                            break
                if cat_type == 'aperture':
                    total_product_catalogs.catalogs[cat_type].sources = total_product_catalogs.catalogs[cat_type].sources[good_rows_index]
                else:
                    total_product_catalogs.catalogs[cat_type].source_cat = total_product_catalogs.catalogs[cat_type].source_cat[good_rows_index]

            # Determine whether any catalogs should be written out at all based on comparison to expected
            # rate of cosmic-ray contamination for the total detection product
            reject_catalogs = total_product_catalogs.verify_crthresh(n1_exposure_time)

            if diagnostic_mode:
                # If diagnostic mode, we want to inspect the original full source catalogs
                reject_catalogs = False

            for filter_product_obj in total_product_obj.fdp_list:
                filter_product_catalogs = filter_catalogs[filter_product_obj.drizzle_filename]

                # Now write the catalogs out for this filter product
                log.info("Writing out filter product catalog")
                # Write out photometric (filter) catalog(s)
                filter_product_catalogs.write(reject_catalogs)

                # append filter product catalogs to list
                if phot_mode in ['aperture', 'both']:
                    product_list.append(filter_product_obj.point_cat_filename)
                if phot_mode in ['segment', 'both']:
                    product_list.append(filter_product_obj.segment_cat_filename)

            log.info("Writing out total product catalog")
            # write out list(s) of identified sources
            total_product_catalogs.write(reject_catalogs)

            # append total product catalogs to manifest list
            if phot_mode in ['aperture', 'both']:
                product_list.append(total_product_obj.point_cat_filename)
            if phot_mode in ['segment', 'both']:
                product_list.append(total_product_obj.segment_cat_filename)
    return product_list


# ----------------------------------------------------------------------------------------------------------------------

def create_drizzle_products(total_obj_list):
    """
    Run astrodrizzle to produce products specified in the total_obj_list.

    Parameters
    ----------
    total_obj_list : list
        List of TotalProduct objects, one object per instrument/detector combination is
        a visit.  The TotalProduct objects are comprised of FilterProduct and ExposureProduct
        objects.

    RETURNS
    -------
    product_list : list
        A list of output products
    """
    log.info("Processing with astrodrizzle version {}".format(drizzlepac.astrodrizzle.__version__))
    # Get rules files
    rules_files = {}

    # Generate list of all input exposure filenames that are to be processed
    edp_names = []
    for t in total_obj_list:
        edp_names += [e.full_filename for e in t.edp_list]

    # Define dataset-specific rules filenames for each input exposure
    for imgname in edp_names:
        rules_files[imgname] = proc_utils.get_rules_file(imgname)

    print('Generated RULES_FILE names of: \n{}\n'.format(rules_files))

    # Keep track of all the products created for the output manifest
    product_list = []

    # For each detector (as the total detection product are instrument- and detector-specific),
    # create the drizzle-combined filtered image, the drizzled exposure (aka single) images,
    # and finally the drizzle-combined total detection image.
    for total_obj in total_obj_list:
        # Need to have direct exposures to drizzle
        if total_obj.edp_list:
            log.info("~" * 118)
            # Get the common WCS for all images which are part of a total detection product,
            # where the total detection product is detector-dependent.
            meta_wcs = total_obj.generate_metawcs()

            # Create drizzle-combined filter image as well as the single exposure drizzled image
            for filt_obj in total_obj.fdp_list:
                log.info("~" * 118)
                filt_obj.rules_file = rules_files[filt_obj.edp_list[0].full_filename]

                log.info("CREATE DRIZZLE-COMBINED FILTER IMAGE: {}\n".format(filt_obj.drizzle_filename))
                filt_obj.wcs_drizzle_product(meta_wcs)
                product_list.append(filt_obj.drizzle_filename)
                product_list.append(filt_obj.trl_filename)

                # Create individual single drizzled images
                for exposure_obj in filt_obj.edp_list:
                    log.info("~" * 118)
                    exposure_obj.rules_file = rules_files[exposure_obj.full_filename]

                    log.info("CREATE SINGLE DRIZZLED IMAGE: {}".format(exposure_obj.drizzle_filename))
                    exposure_obj.wcs_drizzle_product(meta_wcs)
                    product_list.append(exposure_obj.drizzle_filename)
                    product_list.append(exposure_obj.full_filename)
                    # product_list.append(exposure_obj.headerlet_filename)
                    product_list.append(exposure_obj.trl_filename)

            # Create drizzle-combined total detection image after the drizzle-combined filter image and
            # drizzled exposure images in order to take advantage of the cosmic ray flagging.
            log.info("CREATE DRIZZLE-COMBINED TOTAL IMAGE: {}\n".format(total_obj.drizzle_filename))
            total_obj.rules_file = total_obj.fdp_list[0].rules_file
            total_obj.wcs_drizzle_product(meta_wcs)
            product_list.append(total_obj.drizzle_filename)
            product_list.append(total_obj.trl_filename)

    # Ensure that all drizzled products have headers that are to specification
    try:
        log.info("Updating these drizzle products for CAOM compatibility:")
        fits_files = fnmatch.filter(product_list, "*dr?.fits")
        for filename in fits_files:
            log.info("    {}".format(filename))
            proc_utils.refine_product_headers(filename, total_obj_list)
    except Exception:
        log.critical("Trouble updating drizzle products for CAOM.")
        exc_type, exc_value, exc_tb = sys.exc_info()
        traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)
        logging.exception("message")
    # Remove rules files copied to the current working directory
    for rules_filename in list(rules_files.values()):
        log.info("Removed rules file {}".format(rules_filename))
        os.remove(rules_filename)

    # Add primary header information to all objects
    for total_obj in total_obj_list:
        if total_obj.edp_list:
            total_obj = poller_utils.add_primary_fits_header_as_attr(total_obj)
            for filt_obj in total_obj.fdp_list:
                filt_obj = poller_utils.add_primary_fits_header_as_attr(filt_obj)
                for exposure_obj in filt_obj.edp_list:
                    exposure_obj = poller_utils.add_primary_fits_header_as_attr(exposure_obj)

    # Return product list for creation of pipeline manifest file
    return product_list

# ----------------------------------------------------------------------------------------------------------------------


def run_hap_processing(input_filename, diagnostic_mode=False, input_custom_pars_file=None,
                       output_custom_pars_file=None, phot_mode="both", log_level=logutil.logging.INFO):
    """
    Run the HST Advanced Products (HAP) generation code.  This routine is the sequencer or
    controller which invokes the high-level functionality to process the single visit data.

    Parameters
    ----------
    input_filename : str
        The 'poller file' where each line contains information regarding an exposures taken
        during a single visit.

    diagnostic_mode : bool, optional
        Allows printing of additional diagnostic information to the log.  Also, can turn on
        creation and use of pickled information.

    input_custom_pars_file : str, optional
        Represents a fully specified input filename of a configuration JSON file which has been
        customized for specialized processing. This file should contain ALL the input parameters necessary
        for processing. If not specified, default configuration parameter values will be used. The default is
        None.

    output_custom_pars_file : str, optional
        Fully specified output filename which contains all the configuration parameters
        available during the processing session.  The default is None.

    phot_mode : str, optional
        Which algorithm should be used to generate the sourcelists? 'aperture' for aperture photometry;
        'segment' for segment map photometry; 'both' for both 'segment' and 'aperture'. Default value is
        'both'.

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the
        .log file. Default value is 20, or 'info'.


    RETURNS
    -------
    return_value: int
        A return exit code used by the calling Condor/OWL workflow code: 0 (zero) for success, 1 for error
    """
    # This routine needs to return an exit code, return_value, for use by the calling
    # Condor/OWL workflow code: 0 (zero) for success, 1 for error condition
    return_value = 0
    log.setLevel(log_level)
    # Define trailer file (log file) that will contain the log entries for all processing
    logname = proc_utils.build_logname(input_filename)

    # Initialize total trailer filename as temp logname
    logging.basicConfig(filename=logname, format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)
    # start processing
    starting_dt = datetime.datetime.now()
    log.info("Run start time: {}".format(str(starting_dt)))

    # Start by reading in any environment variable related to catalog generation that has been set
    cat_switches = {sw: _get_envvar_switch(sw, default=envvar_cat_svm[sw]) for sw in envvar_cat_svm}

    total_obj_list = []
    manifest_name = ""
    try:
        # Parse the poller file and generate the the obs_info_dict, as well as the total detection
        # product lists which contain the ExposureProduct, FilterProduct, and TotalProduct objects
        # A poller file contains visit data for a single instrument.  The TotalProduct discriminant
        # is the detector.  A TotalProduct object is comprised of FilterProducts and ExposureProducts
        # where its FilterProduct is distinguished by the filter in use, and the ExposureProduct
        # is the atomic exposure data. Note: the TotalProduct was enhanced to also be comprised
        # of an GrismExposureProduct list which is exclusive to the TotalProduct.
        log.info("Parse the poller and determine what exposures need to be combined into separate products.\n")
        obs_info_dict, total_obj_list = poller_utils.interpret_obset_input(input_filename, log_level)

        # Generate the name for the manifest file which is for the entire visit.  It is fine
        # to use only one of the Total Products to generate the manifest name as the name is not
        # dependent on the detector.
        # Example: instrument_programID_obsetID_manifest.txt (e.g.,wfc3_b46_06_manifest.txt)
        manifest_name = total_obj_list[0].manifest_name
        log.info("Generate the manifest name for this visit.")
        log.info("The manifest will contain the names of all the output products.")

        # The product_list is a list of all the output products which will be put into the manifest file
        product_list = []

        # It is possible the total_obj_list output from the poller_utils contains only Grism/Prism
        # data and no direct images, so no further processing should be done.  If this is the case,
        # there is actually nothing to be done for the visit, except write out a manifest and
        # a log file.  Check every item in the total data product list.
        found_data = False
        for total_item in total_obj_list:
            if total_item.edp_list and not found_data:
                found_data = True
            elif total_item.grism_edp_list:
                no_data_trl = total_item.trl_filename
        if not found_data:
            log.warning("")
            log.warning("There are no viable direct images in any Total Data Product for this visit. No processing can be done.")
            log.warning("No SVM processing is done for the Grism/Prism data - no SVM output products are generated.")
            product_list += [no_data_trl]
            sys.exit(0)

        # Update all of the product objects with their associated configuration information.
        for total_item in total_obj_list:
            log.info("Preparing configuration parameter values for total product {}".format(total_item.drizzle_filename))

            total_item.configobj_pars = config_utils.HapConfig(total_item,
                                                               log_level=log_level,
                                                               input_custom_pars_file=input_custom_pars_file,
                                                               output_custom_pars_file=output_custom_pars_file)
            for filter_item in total_item.fdp_list:
                log.info("Preparing configuration parameter values for filter product {}".format(filter_item.drizzle_filename))
                filter_item.configobj_pars = config_utils.HapConfig(filter_item,
                                                                    log_level=log_level,
                                                                    input_custom_pars_file=input_custom_pars_file,
                                                                    output_custom_pars_file=output_custom_pars_file)
            for expo_item in total_item.edp_list:
                log.info("Preparing configuration parameter values for exposure product {}".format(expo_item.drizzle_filename))
                expo_item.configobj_pars = config_utils.HapConfig(expo_item,
                                                                  log_level=log_level,
                                                                  input_custom_pars_file=input_custom_pars_file,
                                                                  output_custom_pars_file=output_custom_pars_file)
                expo_item = poller_utils.add_primary_fits_header_as_attr(expo_item, log_level)

            log.info("The configuration parameters have been read and applied to the drizzle objects.")

            # Check to ensure the object contains a list of direct images (protects against processing Grism/Prism data)
            if total_item.edp_list:
                reference_catalog = run_align_to_gaia(total_item, log_level=log_level, diagnostic_mode=diagnostic_mode)
                if reference_catalog:
                    product_list += reference_catalog

            # Need to delete the Ramp filter Exposure objects from the *Product lists as
            # these images should not be processed beyond the alignment to Gaia (run_align_to_gaia).
            ramp_product_list = []
            _ = delete_ramp_exposures(total_item.fdp_list, 'FILTER')
            ramp_product_list = delete_ramp_exposures(total_item.edp_list, 'EXPOSURE')
            product_list += ramp_product_list

            # If there are Grism/Prism images present in this visit, as well as corresponding direct images
            # for the same detector, update the primary WCS in the direct and/or Grism/Prism images as
            # appropriate to be an 'a priori' or the pipeline default (fallback) solution.  Note: there
            # is no need to delete the Grism/Prism lists after the WCS update as the Grism/Prism
            # exposures are stored in an list ignored by a majority of the processing.
            grism_product_list = []
            if total_item.grism_edp_list and total_item.edp_list:
                grism_product_list = update_wcs_in_visit(total_item)
                product_list += grism_product_list

            # Need a trailer file to log the situation in this special case of a total data product in the visit with
            # only Grism/Prism and no direct images.  No processing is actually done in this case.
            if total_item.grism_edp_list and not total_item.edp_list:
                log.warning("This Total Data Product only has Grism/Prism data and no direct images: {}".format(total_item.drizzle_filename))
                log.warning("No SVM processing is done for the Grism/Prism data - no SVM output products are generated.")
                product_list += [total_item.trl_filename]

        # Run AstroDrizzle to produce drizzle-combined products
        log.info("\n{}: Create drizzled imagery products.".format(str(datetime.datetime.now())))
        driz_list = create_drizzle_products(total_obj_list)
        product_list += driz_list

        # Create source catalogs from newly defined products (HLA-204)
        log.info("{}: Create source catalog from newly defined product.\n".format(str(datetime.datetime.now())))
        if "total detection product 00" in obs_info_dict.keys():
            catalog_list = create_catalog_products(total_obj_list, log_level,
                                                   diagnostic_mode=diagnostic_mode,
                                                   phot_mode=phot_mode,
                                                   catalog_switches=cat_switches)
            product_list += catalog_list
        else:
            log.warning("No total detection product has been produced. The sourcelist generation step has been skipped")

        # Store total_obj_list to a pickle file to speed up development
        if log_level <= logutil.logging.DEBUG:
            pickle_filename = "total_obj_list_full.pickle"
            if os.path.exists(pickle_filename):
                os.remove(pickle_filename)
            pickle_out = open(pickle_filename, "wb")
            pickle.dump(total_obj_list, pickle_out)
            pickle_out.close()
            log.info("Successfully wrote total_obj_list to pickle file {}!".format(pickle_filename))

        # Quality assurance portion of the processing - done only if the environment
        # variable, SVM_QUALITY_TESTING, is set to 'on', 'yes', or 'true'.
        qa_switch = _get_envvar_switch(envvar_qa_svm)

        # If requested, generate quality assessment statistics for the SVM products
        if qa_switch:
            log.info("SVM Quality Assurance statistics have been requested for this dataset, {}.".format(input_filename))
            svm_qa.run_quality_analysis(total_obj_list, log_level=log_level)

        # 10: Return exit code for use by calling Condor/OWL workflow code: 0 (zero) for success, 1 for error condition
        return_value = 0
    except Exception:
        return_value = 1
        print("\a\a\a")
        exc_type, exc_value, exc_tb = sys.exc_info()
        traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)
        logging.exception("message")

    finally:
        # Try to ensure there is a manifest filename if the code exits in the poller_utils.py
        # so the shutdown is tidy
        ntokens = len(input_filename.split("_"))
        if manifest_name == "" and input_filename.lower().endswith("input.out") and ntokens == 4:
            manifest_name = input_filename.lower().replace("input.out", "manifest.txt")
        else:
            manifest_name = "manifest.txt"

        # Write out manifest file listing all products generated during processing
        log.info("Creating manifest file {}.".format(manifest_name))
        log.info("  The manifest contains the names of products generated during processing.")
        with open(manifest_name, mode='w') as catfile:
            if total_obj_list:
                [catfile.write("{}\n".format(name)) for name in product_list]

        end_dt = datetime.datetime.now()
        log.info('Processing completed at {}'.format(str(end_dt)))
        log.info('Total processing time: {} sec'.format((end_dt - starting_dt).total_seconds()))
        log.info("Return exit code for use by calling Condor/OWL workflow code: 0 (zero) for success, 1 for error ")
        log.info("Return condition {}".format(return_value))
        logging.shutdown()

        # The Grism/Prism SVM FLT/FLC images which have had their WCS reconciled with the
        # corresponding direct images need trailer files.  This is also true of the Ramp images
        # which have only been processed through the "align to Gaia" stage.  At this time, just
        # copy the total trailer file, and rename it appropriately.
        if total_obj_list:
            for tot_obj in total_obj_list:
                for gitem in grism_product_list:
                    if gitem.endswith('trl.txt'):
                        shutil.copy(logname, gitem)
                for ritem in ramp_product_list:
                    if ritem.endswith('trl.txt'):
                        shutil.copy(logname, ritem)

        # Append total trailer file (from astrodrizzle) to all total log files
        if total_obj_list:
            for tot_obj in total_obj_list:
                proc_utils.append_trl_file(tot_obj.trl_filename, logname, clean=False)
                if tot_obj.edp_list:
                    # Update DRIZPARS keyword value with new logfile name in ALL drizzle products
                    tot_obj.update_drizpars()

        # Now remove single temp log file
        if os.path.exists(logname):
            os.remove(logname)
        else:
            print("Master log file not found.  Please check logs to locate processing messages.")

        return return_value


# ------------------------------------------------------------------------------------------------------------

def run_align_to_gaia(tot_obj, log_level=logutil.logging.INFO, diagnostic_mode=False):
    # Run align.py on all input images sorted by overlap with GAIA bandpass
    log.info("\n{}: Align the all filters to GAIA with the same fit".format(str(datetime.datetime.now())))
    gaia_obj = None
    headerlet_filenames = []

    # Start by creating a FilterProduct instance which includes ALL input exposures
    for exp_obj in tot_obj.edp_list:
        if gaia_obj is None:
            prod_list = exp_obj.info.split("_")
            prod_list[4] = "metawcs"
            gaia_obj = product.FilterProduct(prod_list[0], prod_list[1], prod_list[2],
                                             prod_list[3], prod_list[4], "all",
                                             prod_list[5][0:3], log_level)
            gaia_obj.configobj_pars = tot_obj.configobj_pars
        gaia_obj.add_member(exp_obj)

    log.info("\n{}: Combined all filter objects in gaia_obj".format(str(datetime.datetime.now())))

    # Now, perform alignment to GAIA with 'match_relative_fit' across all inputs
    # Need to start with one filt_obj.align_table instance as gaia_obj.align_table
    #  - append imglist from each filt_obj.align_table to the gaia_obj.align_table.imglist
    #  - reset group_id for all members of gaia_obj.align_table.imglist to the unique incremental values
    #  - run gaia_obj.align_table.perform_fit() with 'match_relative_fit' only
    #  - migrate updated WCS solutions to exp_obj instances, if necessary (probably not?)
    #  - re-run tot_obj.generate_metawcs() method to recompute total object meta_wcs based on updated
    #    input exposure's WCSs
    catalog_list = gaia_obj.configobj_pars.pars['alignment'].pars_multidict['all']['run_align']['mosaic_catalog_list']
    align_table, filt_exposures = gaia_obj.align_to_gaia(catalog_list=catalog_list,
                                                         output=diagnostic_mode,
                                                         fit_label='SVM')

    tot_obj.generate_metawcs()

    log.info("\n{}: Finished aligning gaia_obj to GAIA".format(str(datetime.datetime.now())))
    log.info("ALIGNED WCS: \n{}".format(tot_obj.meta_wcs))

    # Clean up if the align_table does not exist - the metawcs file
    # may already be deleted, but make sure here.
    if align_table is None:
        headerlet_filenames = []
        try:
            os.remove(gaia_obj.refname)
            gaia_obj.refname = None
        except (OSError, TypeError):
            pass
    else:
        # Get names of all headerlet files written out to file
        headerlet_filenames = [f for f in align_table.filtered_table['headerletFile'] if f != "None"]

    return headerlet_filenames
# ----------------------------------------------------------------------------------------------------------------------


def run_sourcelist_flagging(filter_product_obj, filter_product_catalogs, log_level, diagnostic_mode=False):
    """
    Super-basic and profoundly inelegant interface to hla_flag_filter.py.

    Execute haputils.hla_flag_filter.run_source_list_flaging() to populate the "Flags" column in the catalog tables
    generated by HAPcatalogs.measure().

    Parameters
    ----------
    filter_product_obj : drizzlepac.haputils.product.FilterProduct object
        object containing all the relevant info for the drizzled filter product

    filter_product_catalogs : drizzlepac.haputils.catalog_utils.HAPCatalogs object
        drizzled filter product catalog object

    log_level : int
        The desired level of verboseness in the log statements displayed on the screen and written to the .log file.

    diagnostic_mode : Boolean, optional.
        create intermediate diagnostic files? Default value is False.

    Returns
    -------
    filter_product_catalogs : drizzlepac.haputils.catalog_utils.HAPCatalogs object
        updated version of filter_product_catalogs object with fully populated source flags

    """
    drizzled_image = filter_product_obj.drizzle_filename
    flt_list = []
    for edp_obj in filter_product_obj.edp_list:
        flt_list.append(edp_obj.full_filename)
    param_dict = filter_product_obj.configobj_pars.as_single_giant_dict()
    plate_scale = wcsutil.HSTWCS(drizzled_image, ext=('sci', 1)).pscale
    median_sky = filter_product_catalogs.image.bkg_median

    # Create mask array that will be used by hla_flag_filter.hla_nexp_flags() for both point and segment catalogs.
    if not hasattr(filter_product_obj, 'hla_flag_msk'):
        filter_product_obj.hla_flag_msk = hla_flag_filter.make_mask_array(drizzled_image)

    if filter_product_obj.configobj_pars.use_defaults:
        ci_lookup_file_path = "default_parameters/any"
    else:
        ci_lookup_file_path = "user_parameters/any"
    output_custom_pars_file = filter_product_obj.configobj_pars.output_custom_pars_file
    for cat_type in filter_product_catalogs.catalogs.keys():
        exptime = filter_product_catalogs.catalogs[cat_type].image.imghdu[0].header['exptime']  # TODO: This works for ACS. Make sure that it also works for WFC3. Look at "TEXPTIME"
        catalog_name = filter_product_catalogs.catalogs[cat_type].sourcelist_filename
        catalog_data = filter_product_catalogs.catalogs[cat_type].source_cat
        drz_root_dir = os.getcwd()
        log.info("Run source list flagging on catalog file {}.".format(catalog_name))

        # TODO: REMOVE BELOW CODE ONCE FLAGGING PARAMS ARE OPTIMIZED
        write_flag_filter_pickle_file = False
        if write_flag_filter_pickle_file:
            pickle_dict = {"drizzled_image": drizzled_image,
                           "flt_list": flt_list,
                           "param_dict": param_dict,
                           "exptime": exptime,
                           "plate_scale": plate_scale,
                           "median_sky": median_sky,
                           "catalog_name": catalog_name,
                           "catalog_data": catalog_data,
                           "cat_type": cat_type,
                           "drz_root_dir": drz_root_dir,
                           "hla_flag_msk": filter_product_obj.hla_flag_msk,
                           "ci_lookup_file_path": ci_lookup_file_path,
                           "output_custom_pars_file": output_custom_pars_file,
                           "log_level": log_level,
                           "diagnostic_mode": diagnostic_mode}
            out_pickle_filename = catalog_name.replace("-cat.ecsv", "_flag_filter_inputs.pickle")
            pickle_out = open(out_pickle_filename, "wb")
            pickle.dump(pickle_dict, pickle_out)
            pickle_out.close()
            log.info("Wrote hla_flag_filter param pickle file {} ".format(out_pickle_filename))
        # TODO: REMOVE ABOVE CODE ONCE FLAGGING PARAMS ARE OPTIMIZED
        if catalog_data is not None and len(catalog_data) > 0:
             source_cat = hla_flag_filter.run_source_list_flagging(drizzled_image,
                                                                   flt_list,
                                                                   param_dict,
                                                                   exptime,
                                                                   plate_scale,
                                                                   median_sky,
                                                                   catalog_name,
                                                                   catalog_data,
                                                                   cat_type,
                                                                   drz_root_dir,
                                                                   filter_product_obj.hla_flag_msk,
                                                                   ci_lookup_file_path,
                                                                   output_custom_pars_file,
                                                                   log_level,
                                                                   diagnostic_mode)
        else:
            source_cat = catalog_data

        filter_product_catalogs.catalogs[cat_type].source_cat = source_cat

    return filter_product_catalogs


def _get_envvar_switch(envvar_name, default=None):
    """
    This private routine interprets any environment variable, such as SVM_QUALITY_TESTING.

    PARAMETERS
    -----------
    envvar_name : str
        name of environment variable to be interpreted

    default : str or None
        Value to be used in case environment variable was not defined or set.

    .. note :
    This is a copy of the routine in runastrodriz.py.  This code should be put in a common place.

    """
    if envvar_name in os.environ:
        val = os.environ[envvar_name].lower()
        if val not in envvar_bool_dict:
            msg = "ERROR: invalid value for {}.".format(envvar_name)
            msg += "  \n    Valid Values: on, off, yes, no, true, false"
            raise ValueError(msg)
        switch_val = envvar_bool_dict[val]
    else:
        switch_val = envvar_bool_dict[default] if default else None

    return switch_val

# ------------------------------------------------------------------------------


def update_wcs_in_visit(tdp):
    """Examine entries in the total data product for Grism/Prism data

    This routine is only invoked if the total data product for a particular detector
    for the visit is comprised of both an ExposureProduct list AND a
    GrismExposureProduct list which contain entries.

    *** If the visit contains Grism/Prism data, then either the 'a priori' WCS for
    the Grism/Prism data or the pipeline default WCS for the direct images will be
    set as the primary/active WCS for all the Grism/Prism and direct image data for
    that detector - whichever WCS is common to all the images.

    Note: The direct exposure image active WCS solutions will be modified in-place.

    Parameters
    ----------
    tdp : total data product (TotalProduct) object
        Object comprised of FilterProduct, ExposureProduct, and GrismExposureProduct

    Returns
    -------
    grism_product_list : list
        List of all the SVM Grism/Prism FLT/FLC files updated with the common WCS
        and their corresponding trailer filenames

    """
    log.info("\n***** Grism/Prism Image Processing *****")
    # The TotalProduct (tdp) for this instrument/dectector has both a Grism/Prism
    # exposure list and a direct exposure list - both with contents.
    # Every image should (!) have an IDC_?????????-GSC240 solution.
    wcs_preference = ['IDC_?????????-GSC240', 'IDC_?????????']

    # Grism output product list for the manifest
    grism_product_list = []

    grism_wcs_set, skip_grism_list, g_keyword_wcs_names_dict, grism_dict = collect_wcs_names(tdp.grism_edp_list, 'GRISM')
    log.info("WCS solutions common to all viable Grism/Prism images: {}".format(grism_wcs_set))

    # There is a preference for the active WCS for the viable images in the visit
    # Check the grism_wcs_set for the preferential solutions
    # If a common name is found here, wait to report it as the active WCS may still need to
    # be rolled back depending upon the WCS solutions available in the direct images.
    match_list = []
    grism_wcsname = ''
    for wcs_item in wcs_preference:
        match_list = fnmatch.filter(grism_wcs_set, wcs_item)
        if match_list:
            grism_wcsname = match_list[0]
            log.info("Proposed WCS for use after Grism/Prism examination: {}".format(grism_wcsname))
            break

    if not grism_wcsname:
        log.error("")
        log.error("None of the preferred WCS names are present in the common set of WCS names for the Grism/Prism images.")
        log.error("There is a problem with this visit.  Deleting all SVM Grism/Prism FLT/FLC files.")
        log.error("")
        try:
            for image_file in tdp.grism_edp_list:
               os.remove(image_file.full_filename)
               log.warning("Deleted Grism/Prism image {}.".format(image_file.full_filename))
            tdp.grism_edp_list = []
            return grism_product_list
        except OSError:
            pass
        sys.exit(1)

        log.info("\nProposed WCS for use after Grism/Prism examination: {}".format(grism_wcsname))

    # If it is now known there is at least one viable Grism/Prism observation, it is desired to make
    # the active Grism/Prism WCS be the same for all of the Grism/Prism as well as the direct images
    # in the visit.  The catch can be if any of the direct images do not have an 'a priori' solution
    # which matches the Grism/Prism images.  In this case, it will be necessary to "fall back" to
    # the pipeline WCS solution.
    log.info("\n***** Direct Image Processing *****")

    # Loop over all the direct images for this detector in the visit to update the WCS
    for edp in tdp.edp_list:
        filename = edp.full_filename
        hdu = fits.open(filename)

        # Make sure the WCS is up-to-date - doing the update for the
        # direct images here so there is no need to modify a
        # pre-existing class which could cause an issue
        drizcorr = hdu[0].header['DRIZCORR']
        if drizcorr == "OMIT":
            updatewcs.updatewcs(filename, use_db=True)

        hdu.close()
        # Insure HDRNAME keywords are properly populated in SCI extensions.
        _verify_sci_hdrname(filename)

    direct_wcs_set, skip_direct_list, d_keyword_wcs_names_dict, direct_dict = collect_wcs_names(tdp.edp_list, 'DIRECT')
    log.info("WCS solutions common to all viable direct images: {}".format(direct_wcs_set))

    # Are the grism_wcs_set and the direct_wcs_set disjoint?  If they are disjoint, there can
    # be no common WCS solution.
    is_disjoint = grism_wcs_set.isdisjoint(direct_wcs_set)
    if not is_disjoint:
        # Generate the intersection between the Grism/Prism and direct images
        grism_wcs_set &= direct_wcs_set
        log.info("Common WCS solutions exist between the Grism/Prism and the direct images.")
        log.info("    The common WCS solutions are: {}".format(grism_wcs_set))

        # There is a preference for the active WCS for the viable images in the visit
        # Check the grism_wcs_set for the preferential solutions
        match_list = []
        final_wcsname = ''
        for wcs_item in wcs_preference:
            match_list = fnmatch.filter(grism_wcs_set, wcs_item)
            if match_list:
                final_wcsname = match_list[0]
                log.info("Final WCS solution to use for all Grism/Prism and direct images: {}".format(final_wcsname))
                break

        if final_wcsname:
            # Finally, if the image is not in a skip list, reset the primary WCS in all the images

            for g_edp in tdp.grism_edp_list:
                filename = g_edp.full_filename
                if filename not in skip_grism_list:
                    log.info("")
                    log.info("Setting the primary WCS for Grism/Prism image {} to {}.".format(filename, final_wcsname))
                    update_active_wcs(filename, final_wcsname)

                    # Add the Grism/Prism images to the manifest as all of the files exist.
                    # Also, add the Grism/Prism trailer filename here.  The file will be created by
                    # by copying the Total Data Product trailer file at the end of processing.
                    grism_product_list.append(filename)
                    grism_product_list.append(g_edp.trl_filename)

            for edp in tdp.edp_list:
                filename = edp.full_filename
                if filename not in skip_direct_list:
                    log.info("")
                    log.info("Setting the primary WCS for direct image {} to {}.".format(filename, final_wcsname))
                    update_active_wcs(filename, final_wcsname)
        else:
            # Do nothing
            pass

    # Disjoint
    else:
        log.info("The Grism/Prism and direct images in this visit have no WCS solution common to both sets of data.")
        log.info("The active WCS solution for each Grism/Prism and direct image is not changed.")

    log.info("Grism product list: {}".format(grism_product_list))

    return grism_product_list

# ------------------------------------------------------------------------------


def collect_wcs_names(edp_list, image_type):
    """
    Utility to collect all the WCS solution names common to the input image list

    Parameters
    ----------
    edp_list: str list
        List containing the SVM FLT/FLC filenames

    image_type: string
        String containing either 'GRISM' or 'DIRECT' to use as a switch for
        output information

    Returns
    -------
    image_wcs_set: set of WCS solutions
        The set contains the WCS solution names common to all of the
        input exposures

    skip_image_list: list
        This is a list of exposures in the input list which should be
        skipped/ignored when updating the active WCS solution

    keyword_wcs_names_dict: dictionary {filename: list}
        The dictionary is used to associate an individual image/filename with
        a list of WCS solution names in the file stored as keywords (not headerlets)

    image_dict: dictionary {filename: list}
        The dictionary is used to associate an individual image/filename with
        a list of *all* WCS solution names in the file

    """

    image_wcs_set = set()
    skip_image_list = []
    exist_image_set = False
    image_dict = {}
    keyword_wcs_names_dict = {}
    # Loop over all the Grism/Prism images for this detector in the visit
    for edp in edp_list:

        filename = edp.full_filename

        # Get all the WCS names which are common to all of the images
        # Note that WCS solutions may be represented as FITS keyword values in the
        # SCI extension and/or as headerlets in the HDRLET extensions.
        # Get the keyword WCS solution names.
        keyword_wcs_names = list(wcsutil.altwcs.wcsnames(filename, ext=1).values())

        # Get the headerlet WCS solution names
        headerlet_wcs_names = wcsutil.headerlet.get_headerlet_kw_names(filename, kw="WCSNAME")
        all_wcs_names = keyword_wcs_names + headerlet_wcs_names
        keyword_wcs_names_dict[filename] = keyword_wcs_names
        image_dict[filename] = all_wcs_names
        if all_wcs_names:
            log.debug("WCS solutions for file {} are {}.".format(filename, all_wcs_names))
            # Initialize a set with wcsnames
            if not exist_image_set:
                image_wcs_set = set(all_wcs_names)
                exist_image_set = True
            # Generate the intersection with the set of existing wcsnames and the wcsnames
            # from the current image
            else:
                image_wcs_set &= set(all_wcs_names)

            # Oops...no common wcsnames
            if not image_wcs_set:
                log.error("There are no common WCS solutions with this image {} and previously processed images".format(filename))
                log.error("    There is a problem with this image/visit.")
                log.error("    Make sure the input data are not *_raw.fits files.")
                sys.exit(1)
        # If there are no WCS solutions, the image could be bad (e.g., EXPTIME=0 or EXPFLAG="TDF-DOWN...")
        else:
            log.warning("There are no WCS solutions in the image {} in this visit.".format(filename))
            skip_image_list.append(filename)
            if image_type == 'GRISM':
                log.warning("    Skip and delete this image.")
                # Delete the SVM FLT/FlC Grism/Prism image as it has no updated WCS
                try:
                    os.remove(filename)
                    log.warning("Deleted Grism/Prism image {}.".format(filename))
                except OSError:
                    pass
            else:
                log.warning("    Skip this image.")

    return image_wcs_set, skip_image_list, keyword_wcs_names_dict, image_dict

# ------------------------------------------------------------------------------


def update_active_wcs(filename, wcsname):
    """
    Utility to update the active/primary WCS solution

    This small utility updates the active/primary WCS solution for the input
    file with the WCS solution indicted by the input parameter "wcsname"

    Parameters
    ----------
    filename : str
        Input/Output SVM FLT/FLC filename - the file is updated in-place

    wcsname : str
        Name of the desired WCS active/primary solution to be set for the filename

    Returns
    -------
    None

    """
    # For exposures with multiple science extensions (multiple chips),
    # generate a combined WCS
    num_sci_ext, extname = util.count_sci_extensions(filename)
    extname_list = [(extname, x + 1) for x in range(num_sci_ext)]

    hdu = fits.open(filename)

    # Check if the desired WCS solution is already the active solution
    # whereupon there is nothing to do
    key = wcsutil.altwcs.getKeyFromName(hdu['SCI', 1].header, wcsname)

    # No need to keep this file handle open anymore
    hdu.close()

    # Case where the desired active solution is not the current active solution
    if key != ' ':
        # Get the distortion model identification of the desired active WCS solution
        tmp_wcsname = wcsname.split('-')[0]
        index = tmp_wcsname.upper().find('IDC_')
        idc_new_string = ''
        if index > -1:
            idc_new_string = tmp_wcsname[index:]

        # Get the headerlet HDRNAMES for comparison to the alternate WCS solutions
        headerlet_hdrnames = wcsutil.headerlet.get_headerlet_kw_names(filename, kw="HDRNAME")

        if headerlet_hdrnames:
            # Examine the alternate WCS solutions to determine if they will be auto-archived due
            # to a distortion model change.  The auto-archiving will happen when the desired WCS
            # solution is installed as the active solution - just deleting duplicates here pro-actively.
            wcs_key_dict = wcsutil.altwcs.wcsnames(filename, ext=1)
            for wkey, wname in wcs_key_dict.items():
                if wkey == ' ':
                    continue

                index = wname.upper().find(idc_new_string.upper())

                # No match so solution will be copied to a headerlet automatically when the new primary is set
                if index == -1 and wkey.upper() != 'O':
                    log.info("Archiving alternate WCS solution as a headerlet as necessary: {}".format(wname))

                    # Now check if the HDRNAME between this solution and a headerlet already exists
                    hdr_keyword = fits.getval(filename, 'HDRNAME{}'.format(wkey.upper()), ext=1)

                    # Solution already exists as a headerlet extension, so just delete it
                    if hdr_keyword in headerlet_hdrnames:
                        wcsutil.altwcs.deleteWCS(filename, extname_list, wcskey=wkey)

            # Get all the WCS solution names
            headerlet_wcsnames = wcsutil.headerlet.get_headerlet_kw_names(filename, kw="WCSNAME")
            keyword_wcs_list = list(wcs_key_dict.values())

            # Prepare to install a new active WCS, but need to do some checking first
            #
            # This returns the first matching instance
            hdrname = headerlet_hdrnames[headerlet_wcsnames.index(wcsname)]
            extensions = []
            extensions = wcsutil.headerlet.find_headerlet_HDUs(filename, hdrname=hdrname)

            # It is possible the hdrname is not unique, so need to delete the dups
            for ext in reversed(extensions[1:]):
                wcsutil.headerlet.delete_headerlet(filename, hdrext=ext)
                log.info("Delete duplicate headerlet extension {} in filename {}.".format(ext, filename))

            log.info("Desired active WCS solution {} has an HDRNAME of {}.".format(wcsname, hdrname))

            # Finally, install the desired WCS as the active WCS solution
            # Is the source of the wcsname for this image from a headerlet extension
            # or from the alternate solutions in the header as the source dictates how
            # the WCS will be made the active WCS. If available, restore a WCS solution
            # from the headerlet extension.
            try:
                wcsutil.headerlet.restore_from_headerlet(filename, hdrname=hdrname, force=True)
            except ValueError as err:
                log.warning("Trapped ValueError - attempting recovery: {}".format(str(err)))
                found_string = [i for i in keyword_wcs_list if wcsname == i]
                if found_string:
                    wcsutil.altwcs.restoreWCS(filename, ext=extname_list, wcsname=found_string[0])
                else:
                    log.warning("Could not restore the common WCS, {}, as the active WCS in this file {}.".format(wcsname, filename))
            except AssertionError:
                _, _, tb = sys.exc_info()
                tb_info = traceback.extract_tb(tb)
                _, _, _, text = tb_info[-1]
                log.warning("Trapped AssertionError: {}.".format(text))
                log.warning("Could not restore the common WCS, {}, as the active WCS in this file {}.".format(wcsname, filename))
        else:
            found_string = [i for i in keyword_wcs_list if wcsname == i]
            if found_string:
                wcsutil.altwcs.restoreWCS(filename, ext=extname_list, wcsname=found_string[0])
            else:
                log.warning("Could not restore the common WCS from alternate WCS solutions, {}, as the active WCS in this file {}.".format(wcsname, filename))
    else:
        log.info("No need to update active WCS solution of {} for {} as it is already the active solution.".format(wcsname, filename))


# ------------------------------------------------------------------------------


def delete_ramp_exposures(obj_list, type_of_list):
    """Delete the Ramp filter objects from the Total Product internal lists

    The Total Data Product (tdp) object is comprised of a list of Filter Product objects,
    as well as a list of Exposure Product objects.  The Ramp filter
    images need to be processed in the same manner as nominal exposures for at least
    some of the processing steps.  Because of this, it was deemed best to keep the
    Ramp exposures in the tdp list attributes, until their final processing
    stage (align_to_gaia), and then delete these objects from the attribute
    lists. This function handles the deletion of Ramp images from the input
    list.

    Parameters
    ----------
    obj_list : list of either FilterProduct ro ExposureProduct objects
        The list is updated in-place.

    type_of_list : str
        String literal indicating if the list contains FILTER or EXPOSURE objects

    Returns
    -------
    ramp_filenames_list : list of str
        List of ramp exposure SVM FLT/FLC filenames and their corresponding trailer filenames

    """
    reported = False
    temp = []
    ramp_filenames_list = []
    while obj_list:
        obj = obj_list.pop()
        if obj.filters.lower().find('fr') == -1:
            temp.append(obj)
        else:
            # Add the Ramp images to the manifest as well as the Ramp trailer filename here.
            # The file will be created by by copying the Total Data Product trailer file at
            # the end of processing.
            if type_of_list == 'EXPOSURE':
                ramp_filenames_list.append(obj.full_filename)
                ramp_filenames_list.append(obj.trl_filename)
                if not reported:
                    log.info('Removing the Ramp images from the Total Data Product exposure list.')
                    log.info('Theses images are not processed beyond the "align to Gaia" stage.')
                    reported = True
    while temp:
        obj_list.append(temp.pop())

    return ramp_filenames_list


def _verify_sci_hdrname(filename):
    """Insures that HDRNAME keyword is populated in SCI extensions.

    This function checks to make sure the HDRNAME keyword in the SCI
    extension of the science image `filename` is populated with a valid
    non-empty string.
    """
    fhdu, closefits = proc_utils._process_input(filename)

    # Find all extensions to be updated
    numext = countExtn(fhdu, extname='SCI')

    for ext in range(1, numext + 1):
        sciext = ('sci', ext)
        scihdr = fhdu[sciext].header
        if 'hdrname' not in scihdr or scihdr['hdrname'].rstrip() == '':
            # We need to create a valid value for the keyword
            # Define new HDRNAME value in case it is needed.
            # Same value for all SCI extensions, so just precompute it and be ready.
            # This code came from 'stwcs.updatewcs.astrometry_utils'
            hdrname = "{}_{}".format(filename.replace('.fits', ''), scihdr['wcsname'])
            # Create full filename for headerlet:
            hfilename = "{}_hlet.fits".format(hdrname)
            # Update the header with the new value, inserting after WCSNAME
            scihdr.set('hdrname', hfilename, 'Name of headerlet file', after='wcsname')

    if closefits:
        fhdu.close()
        del fhdu
