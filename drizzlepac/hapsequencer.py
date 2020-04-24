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
import glob
import logging
import os
import pdb
import pickle
import sys
import traceback

from astropy.table import Table
import numpy as np
import drizzlepac

from drizzlepac.devutils.comparison_tools import compare_sourcelists
from drizzlepac.devutils.comparison_tools.read_hla import read_hla_catalog
from drizzlepac.hlautils import config_utils
from drizzlepac.hlautils import diagnostic_utils
from drizzlepac.hlautils import hla_flag_filter
from drizzlepac.hlautils import poller_utils
from drizzlepac.hlautils import product
from drizzlepac.hlautils import processing_utils as proc_utils
from drizzlepac.hlautils import svm_quality_analysis as svm_qa
from drizzlepac.hlautils.catalog_utils import HAPCatalogs

from stsci.tools import logutil
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

# --------------------------------------------------------------------------------------------------------------


def correct_hla_classic_ra_dec(orig_hla_classic_sl_name, hap_imgname, cattype, log_level):
    """
    This subroutine runs Rick White's read_hla_catalog script to convert the RA and Dec values from a HLA Classic
    sourcelist into same reference frame used by the HAP sourcelists. Additionally, the new RA and Dec values
    are then transformed to produce X and Y coordinates that are in the HAP image frame of reference. A new
    version of the input file with the converted X and Y and RA and Dec values is written to the current
    working directory named <INPUT SOURCELIST NAME>_corrected.txt.

    Parameters
    ----------
    orig_hla_classic_sl_name : string
        name of the HLA Classic sourcelist whose RA and Dec values will be converted.

    hap_imgname : string
        name of HAP image. The WCS info from this image will be used to transform the updated RA and DEC
        values to new X and Y values in the this image's frame of reference

    cattype : string
        HLA Classic catalog type. Either 'sex' (source extractor) or 'dao' (DAOphot).

    log_level : int
        The desired level of verboseness in the log statements displayed on the screen and written to the .log file.

    Returns
    -------
    mod_sl_name : string
        Name of the new version of the input file with the converted RA and Dec values
    """
    try:
        mod_sl_name = os.path.basename(orig_hla_classic_sl_name)

        # Execute read_hla_catalog.read_hla_catalog() to convert RA and Dec values
        dataset = mod_sl_name.replace("_{}phot.txt".format(cattype), "")
        modcat = read_hla_catalog.read_hla_catalog(dataset, cattype=cattype, applyomega=True, multiwave=False,
                                                   verbose=True, trim=False, log_level=log_level)
        # sort catalog with updated RA, DEC values so that ordering is the same as the uncorrected table and everything maps correclty.
        if cattype == "dao":
            sortcoltitle = "ID"
            x_coltitle = "X-Center"
            y_coltitle = "Y-Center"
        if cattype == "sex":
            sortcoltitle = "NUMBER"
            x_coltitle = "X-IMAGE"
            y_coltitle = "Y-IMAGE"

        modcat.sort(sortcoltitle)

        # Identify RA and Dec column names in the new catalog table object
        for ra_col_title in ["ra", "RA", "ALPHA_J2000", "alpha_j2000"]:
            if ra_col_title in modcat.colnames:
                true_ra_col_title = ra_col_title
                log.debug("RA Col_name: {}".format(true_ra_col_title))
                break
        for dec_col_title in ["dec", "DEC", "Dec", "DELTA_J2000", "delta_j2000"]:
            if dec_col_title in modcat.colnames:
                true_dec_col_title = dec_col_title
                log.debug("DEC Col_name: {}".format(true_dec_col_title))
                break

        # transform new RA and Dec values into X and Y values in the HAP reference frame
        ra_dec_values = np.stack((modcat[true_ra_col_title], modcat[true_dec_col_title]), axis=1)
        new_xy_values = hla_flag_filter.rdtoxy(ra_dec_values, hap_imgname, '[sci,1]')

        # get HLA Classic sourcelist data, replace existing RA and Dec column data with the converted RA and Dec column data
        cat = Table.read(orig_hla_classic_sl_name, format='ascii')
        cat['RA'] = modcat[true_ra_col_title]
        cat['DEC'] = modcat[true_dec_col_title]

        # update existing X and Y values with new X and Y values transformed from new RA and Dec values.
        cat[x_coltitle] = new_xy_values[:, 0]
        cat[y_coltitle] = new_xy_values[:, 1]

        # Write updated version of HLA Classic catalog to current working directory
        mod_sl_name = mod_sl_name.replace(".txt", "_corrected.txt")
        log.info("Corrected version of HLA Classic file {} with new X, Y and RA, Dec values written to {}.".format(os.path.basename(orig_hla_classic_sl_name), mod_sl_name))
        cat.write(mod_sl_name, format="ascii.csv")
        return mod_sl_name

    except Exception:
        log.warning("There was a problem converting the RA and Dec values. Using original uncorrected HLA Classic sourcelist instead.")
        log.warning("Comparisons may be of questionable quality.")
        return orig_hla_classic_sl_name


# ----------------------------------------------------------------------------------------------------------------------

def create_catalog_products(total_obj_list, log_level, diagnostic_mode=False, phot_mode='both'):
    """This subroutine utilizes hlautils/catalog_utils module to produce photometric sourcelists for the specified
    total drizzle product and it's associated child filter products.

    Parameters
    ----------
    total_obj_list : drizzlepac.hlautils.Product.TotalProduct
        total drizzle product that will be processed by catalog_utils. catalog_utils will also create photometric
        sourcelists for the child filter products of this total product.

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the .log file.

    diagnostic_mode : bool, optional
        generate ds9 region file counterparts to the photometric sourcelists? Default value is False.

    phot_mode : str, optional
        Which algorithm should be used to generate the sourcelists? 'aperture' for aperture photometry;
        'segment' for segment map photometry; 'both' for both 'segment' and 'aperture'. Default value is 'both'.

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
        # Instantiate filter catalog product object
        total_product_catalogs = HAPCatalogs(total_product_obj.drizzle_filename,
                                             total_product_obj.configobj_pars.get_pars('catalog generation'),
                                             total_product_obj.configobj_pars.get_pars('quality control'),
                                             log_level,
                                             types=input_phot_mode,
                                             diagnostic_mode=diagnostic_mode)

        # Generate an "n" exposure mask which has the image footprint set to the number
        # of exposures which constitute each pixel.
        total_product_obj.generate_footprint_mask()

        # Identify sources in the input image and delay writing the total detection
        # catalog until the photometric measurements have been done on the filter
        # images and some of the measurements can be appended to the total catalog
        total_product_catalogs.identify(mask=total_product_obj.mask)

        # Determine how to continue if "aperture" or "segment" fails to find sources for this total
        # detection product - take into account the initial setting of phot_mode.
        # If no sources were found by either the point or segmentation algorithms, go on to
        # the next total detection product (detector) in the visit with the initially requested
        # phot_mode.  If the point or segmentation algorithms found sources, need to continue
        # processing for that (those) algorithm(s) only.

        # When both algorithms have been requested...
        if input_phot_mode == 'both':
            # If no sources found with either algorithm, skip to the next total detection product
            if total_product_catalogs.catalogs['aperture'].sources is None and total_product_catalogs.catalogs['segment'].sources is None:
                del total_product_catalogs.catalogs['aperture']
                del total_product_catalogs.catalogs['segment']
                continue

            # Only point algorithm found sources, continue to the filter catalogs for just point
            if total_product_catalogs.catalogs['aperture'].sources is not None and total_product_catalogs.catalogs['segment'].sources is None:
                phot_mode = 'aperture'
                del total_product_catalogs.catalogs['segment']

            # Only segment algorithm found sources, continue to the filter catalogs for just segmentation
            if total_product_catalogs.catalogs['aperture'].sources is None and total_product_catalogs.catalogs['segment'].sources is not None:
                phot_mode = 'segment'
                del total_product_catalogs.catalogs['aperture']

        # Only requested the point algorithm
        elif input_phot_mode == 'aperture':
            if total_product_catalogs.catalogs['aperture'].sources is None:
                del total_product_catalogs.catalogs['aperture']
                continue

        # Only requested the segmentation algorithm
        elif input_phot_mode == 'segment':
            if total_product_catalogs.catalogs['segment'].sources is None:
                del total_product_catalogs.catalogs['segment']
                continue

        # Build dictionary of total_product_catalogs.catalogs[*].sources to use for
        # filter photometric catalog generation
        sources_dict = {}
        for cat_type in total_product_catalogs.catalogs.keys():
            sources_dict[cat_type] = {}
            sources_dict[cat_type]['sources'] = total_product_catalogs.catalogs[cat_type].sources
            # For the segmentation source finding, both the segmentation image AND the segmentation catalog
            # computed for the total object need to be provided to the filter objects
            if cat_type == "segment":
                sources_dict['segment']['kernel'] = total_product_catalogs.catalogs['segment'].kernel
                sources_dict['segment']['source_cat'] = total_product_catalogs.catalogs['segment'].source_cat

        log.info("Generating filter product source catalogs")
        for filter_product_obj in total_product_obj.fdp_list:
            # Instantiate filter catalog product object
            filter_product_catalogs = HAPCatalogs(filter_product_obj.drizzle_filename,
                                                  total_product_obj.configobj_pars.get_pars('catalog generation'),
                                                  total_product_obj.configobj_pars.get_pars('quality control'),
                                                  log_level,
                                                  types=phot_mode,
                                                  diagnostic_mode=diagnostic_mode,
                                                  tp_sources=sources_dict)

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

            # write out CI and FWHM values to file (if IRAFStarFinder was used instead of DAOStarFinder) for hla_flag_filter parameter optimization.
            if diagnostic_mode:
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

            log.info("Writing out filter product catalog")
            # Write out photometric (filter) catalog(s)
            filter_product_catalogs.write()

            # Load a dictionary with the filter subset table for each catalog...
            subset_columns_dict = {}
            for cat_type in filter_product_catalogs.catalogs.keys():
                subset_columns_dict[cat_type] = {}
                subset_columns_dict[cat_type]['subset'] = \
                    filter_product_catalogs.catalogs[cat_type].subset_filter_source_cat

            # ...and append the filter columns to the total detection product catalog.
            total_product_catalogs.combine(subset_columns_dict)

            # append filter product catalogs to list
            if phot_mode in ['aperture', 'both']:
                product_list.append(filter_product_obj.point_cat_filename)
            if phot_mode in ['segment', 'both']:
                product_list.append(filter_product_obj.segment_cat_filename)

        log.info("Writing out total product catalog")
        # write out list(s) of identified sources
        total_product_catalogs.write()

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
    total_obj_list: list
        List of TotalProduct objects, one object per instrument/detector combination is
        a visit.  The TotalProduct objects are comprised of FilterProduct and ExposureProduct
        objects.

    RETURNS
    -------
    product_list: list
        A list of output products
    """
    log.info("Processing with astrodrizzle version {}".format(drizzlepac.astrodrizzle.__version__))
    # Get rules files
    for imgname in glob.glob("*fl?.fits"):
        proc_utils.get_rules_file(imgname)

    # Keep track of all the products created for the output manifest
    product_list = []

    # For each detector (as the total detection product are instrument- and detector-specific),
    # create the drizzle-combined filtered image, the drizzled exposure (aka single) images,
    # and finally the drizzle-combined total detection image.
    for total_obj in total_obj_list:
        log.info("~" * 118)
        # Get the common WCS for all images which are part of a total detection product,
        # where the total detection product is detector-dependent.
        meta_wcs = total_obj.generate_metawcs()

        # Create drizzle-combined filter image as well as the single exposure drizzled image
        for filt_obj in total_obj.fdp_list:
            log.info("~" * 118)

            log.info("CREATE DRIZZLE-COMBINED FILTER IMAGE: {}\n".format(filt_obj.drizzle_filename))
            filt_obj.wcs_drizzle_product(meta_wcs)
            product_list.append(filt_obj.drizzle_filename)
            product_list.append(filt_obj.trl_filename)

            # Create individual single drizzled images
            for exposure_obj in filt_obj.edp_list:
                log.info("~" * 118)

                log.info("CREATE SINGLE DRIZZLED IMAGE: {}".format(exposure_obj.drizzle_filename))
                exposure_obj.wcs_drizzle_product(meta_wcs)
                product_list.append(exposure_obj.drizzle_filename)
                product_list.append(exposure_obj.full_filename)
                product_list.append(exposure_obj.headerlet_filename)
                product_list.append(exposure_obj.trl_filename)

        # Create drizzle-combined total detection image after the drizzle-combined filter image and
        # drizzled exposure images in order to take advantage of the cosmic ray flagging.
        log.info("CREATE DRIZZLE-COMBINED TOTAL IMAGE: {}\n".format(total_obj.drizzle_filename))
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
    for rules_filename in glob.glob("*_header_hla.rules"):
        log.info("Removed rules file {}".format(rules_filename))
        os.remove(rules_filename)

    # Add primary header information to all objects
    for total_obj in total_obj_list:
        total_obj = poller_utils.add_primary_fits_header_as_attr(total_obj)
        for filt_obj in total_obj.fdp_list:
            filt_obj = poller_utils.add_primary_fits_header_as_attr(filt_obj)
            for exposure_obj in filt_obj.edp_list:
                exposure_obj = poller_utils.add_primary_fits_header_as_attr(exposure_obj)

    # Return product list for creation of pipeline manifest file
    return product_list

# ----------------------------------------------------------------------------------------------------------------------


def run_hap_processing(input_filename, diagnostic_mode=False, use_defaults_configs=True,
                       input_custom_pars_file=None, output_custom_pars_file=None, phot_mode="both",
                       log_level=logutil.logging.INFO):
    """
    Run the HST Advanced Products (HAP) generation code.  This routine is the sequencer or
    controller which invokes the high-level functionality to process the single visit data.

    Parameters
    ----------
    input_filename: string
        The 'poller file' where each line contains information regarding an exposures taken
        during a single visit.

    diagnostic_mode : bool, optional
        Allows printing of additional diagnostic information to the log.  Also, can turn on
        creation and use of pickled information.

    use_defaults_configs: bool, optional
        If True, use the configuration parameters in the 'default' portion of the configuration
        JSON files.  If False, use the configuration parameters in the "parameters" portion of
        the file.  The default is True.

    input_custom_pars_file: string, optional
        Represents a fully specified input filename of a configuration JSON file which has been
        customized for specialized processing.  This file should contain ALL the input parameters
        necessary for processing.  If there is a filename present for this parameter, the
        'use_defaults_configs' parameter is ignored. The default is None.

    output_custom_pars_file: string, optional
        Fully specified output filename which contains all the configuration parameters
        available during the processing session.  The default is None.

    phot_mode : str, optional
        Which algorithm should be used to generate the sourcelists? 'aperture' for aperture photometry;
        'segment' for segment map photometry; 'both' for both 'segment' and 'aperture'. Default value is 'both'.

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the .log file.
        Default value is 20, or 'info'.


    RETURNS
    -------
    return_value: integer
        A return exit code used by the calling Condor/OWL workflow code: 0 (zero) for success, 1 for error
    """
    # This routine needs to return an exit code, return_value, for use by the calling
    # Condor/OWL workflow code: 0 (zero) for success, 1 for error condition
    return_value = 0
    log.setLevel(log_level)
    # Define trailer file (log file) that will contain the log entries for all processing
    if isinstance(input_filename, str):  # input file is a poller file -- easy case
        logname = input_filename.replace('.out', '.log')
    else:
        logname = 'svm_process.log'
    # Initialize total trailer filename as temp logname
    logging.basicConfig(filename=logname, format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)
    # start processing
    starting_dt = datetime.datetime.now()
    log.info("Run start time: {}".format(str(starting_dt)))
    total_obj_list = []
    product_list = []
    try:
        # Parse the poller file and generate the the obs_info_dict, as well as the total detection
        # product lists which contain the ExposureProduct, FilterProduct, and TotalProduct objects
        # A poller file contains visit data for a single instrument.  The TotalProduct discriminant
        # is the detector.  A TotalProduct object is comprised of FilterProducts and ExposureProducts
        # where its FilterProduct is distinguished by the filter in use, and the ExposureProduct
        # is the atomic exposure data.
        log.info("Parse the poller and determine what exposures need to be combined into separate products.\n")
        obs_info_dict, total_obj_list = poller_utils.interpret_obset_input(input_filename, log_level)

        # Generate the name for the manifest file which is for the entire visit.  It is fine
        # to use only one of the Total Products to generate the manifest name as the name is not
        # dependent on the detector.
        # Example: instrument_programID_obsetID_manifest.txt (e.g.,wfc3_b46_06_manifest.txt)
        manifest_name = total_obj_list[0].manifest_name
        log.info("\nGenerate the manifest name for this visit.")
        log.info("The manifest will contain the names of all the output products.")

        # The product_list is a list of all the output products which will be put into the manifest file
        product_list = []

        # Update all of the product objects with their associated configuration information.
        for total_item in total_obj_list:
            log.info("Preparing configuration parameter values for total product {}".format(total_item.drizzle_filename))
            total_item.configobj_pars = config_utils.HapConfig(total_item,
                                                               log_level=log_level,
                                                               use_defaults=use_defaults_configs,
                                                               input_custom_pars_file=input_custom_pars_file,
                                                               output_custom_pars_file=output_custom_pars_file)
            for filter_item in total_item.fdp_list:
                log.info("Preparing configuration parameter values for filter product {}".format(filter_item.drizzle_filename))
                filter_item.configobj_pars = config_utils.HapConfig(filter_item,
                                                                    log_level=log_level,
                                                                    use_defaults=use_defaults_configs,
                                                                    input_custom_pars_file=input_custom_pars_file,
                                                                    output_custom_pars_file=output_custom_pars_file)
            for expo_item in total_item.edp_list:
                log.info("Preparing configuration parameter values for exposure product {}".format(expo_item.drizzle_filename))
                expo_item.configobj_pars = config_utils.HapConfig(expo_item,
                                                                  log_level=log_level,
                                                                  use_defaults=use_defaults_configs,
                                                                  input_custom_pars_file=input_custom_pars_file,
                                                                  output_custom_pars_file=output_custom_pars_file)
                expo_item = poller_utils.add_primary_fits_header_as_attr(expo_item, log_level)

        log.info("The configuration parameters have been read and applied to the drizzle objects.")

        reference_catalog = run_align_to_gaia(total_obj_list, log_level=log_level, diagnostic_mode=diagnostic_mode)
        if reference_catalog:
            product_list += [reference_catalog]

        # Run AstroDrizzle to produce drizzle-combined products
        log.info("\n{}: Create drizzled imagery products.".format(str(datetime.datetime.now())))
        driz_list = create_drizzle_products(total_obj_list)
        product_list += driz_list

        # Create source catalogs from newly defined products (HLA-204)
        log.info("{}: Create source catalog from newly defined product.\n".format(str(datetime.datetime.now())))
        if "total detection product 00" in obs_info_dict.keys():
            catalog_list = create_catalog_products(total_obj_list, log_level,
                                                   diagnostic_mode=diagnostic_mode,
                                                   phot_mode=phot_mode)
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

            # Compute WCS differences between the primary WCS and the alternates, except for OPUS WCS
            svm_qa.report_wcs(total_obj_list)

        # 9: Compare results to HLA classic counterparts (if possible)
        if diagnostic_mode:
            run_sourcelist_comparision(total_obj_list, diagnostic_mode=diagnostic_mode, log_level=log_level)
        # Write out manifest file listing all products generated during processing
        log.info("Creating manifest file {}.".format(manifest_name))
        log.info("  The manifest contains the names of products generated during processing.")
        with open(manifest_name, mode='w') as catfile:
            [catfile.write("{}\n".format(name)) for name in product_list]
        # 10: Return exit code for use by calling Condor/OWL workflow code: 0 (zero) for success, 1 for error condition
        return_value = 0
    except Exception:
        return_value = 1
        print("\a\a\a")
        exc_type, exc_value, exc_tb = sys.exc_info()
        traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)
        logging.exception("message")

    finally:
        end_dt = datetime.datetime.now()
        log.info('Processing completed at {}'.format(str(end_dt)))
        log.info('Total processing time: {} sec'.format((end_dt - starting_dt).total_seconds()))
        log.info("Return exit code for use by calling Condor/OWL workflow code: 0 (zero) for success, 1 for error ")
        log.info("Return condition {}".format(return_value))
        logging.shutdown()
        # Append total trailer file (from astrodrizzle) to all total log files
        if total_obj_list:
            for tot_obj in total_obj_list:
                proc_utils.append_trl_file(tot_obj.trl_filename, logname, clean=False)
        # Now remove single temp log file
        if os.path.exists(logname):
            os.remove(logname)
        else:
            print("Master log file not found.  Please check logs to locate processing messages.")
        return return_value


# ------------------------------------------------------------------------------------------------------------

def run_align_to_gaia(total_obj_list, log_level=logutil.logging.INFO, diagnostic_mode=False):
    # Run align.py on all input images sorted by overlap with GAIA bandpass
    log.info("\n{}: Align the all filters to GAIA with the same fit".format(str(datetime.datetime.now())))
    gaia_obj = None
    # Start by creating a FilterProduct instance which includes ALL input exposures
    for tot_obj in total_obj_list:
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
        align_table, filt_exposures = gaia_obj.align_to_gaia(output=diagnostic_mode, fit_label='SVM')

        for tot_obj in total_obj_list:
            tot_obj.generate_metawcs()
        log.info("\n{}: Finished aligning gaia_obj to GAIA".format(str(datetime.datetime.now())))

        # Return the name of the alignment catalog
        if align_table is None:
            gaia_obj.refname = None

        return gaia_obj.refname

        #
        # Composite WCS fitting should be done at this point so that all exposures have been fit to GAIA at
        # the same time (on the same frame)
        #

# ----------------------------------------------------------------------------------------------------------------------

def run_sourcelist_comparision(total_list, diagnostic_mode = False, log_level=logutil.logging.INFO):
    """ This subroutine automates execution of drizzlepac/devutils/comparison_tools/compare_sourcelist_flagging.py to
    compare HAP-generated filter catalogs with their HLA classic counterparts.

    NOTE: In order for this subroutine to run, the following environment variables need to be set:
    - HLA_CLASSIC_BASEPATH
    - HLA_BUILD_VER

    Alternatively, if the HLA classic path is unavailable, The comparison can be run using locally stored HLA classic
    files. The relevant HLA classic imagery and sourcelist files must be placed in a subdirectory of the current working
    directory called 'hla_classic'.

    Parameters
    ----------
    total_list: list
        List of TotalProduct objects, one object per instrument/detector combination is
        a visit.  The TotalProduct objects are comprised of FilterProduct and ExposureProduct
        objects.

    diagnostic_mode : Boolean, optional.
        create intermediate diagnostic files? Default value is False.

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the .log file.
        Default value is 20, or 'info'.

    RETURNS
    -------
    Nothing.
    """
    # get HLA classic path details from envroment variables
    hla_classic_basepath = os.getenv('HLA_CLASSIC_BASEPATH')
    hla_build_ver = os.getenv("HLA_BUILD_VER")
    for tot_obj in total_list:
        combo_comp_pdf_list = []
        if hla_classic_basepath and hla_build_ver and os.path.exists(hla_classic_basepath):
            hla_cassic_basepath = os.path.join(hla_classic_basepath, tot_obj.instrument, hla_build_ver)
            hla_classic_path = os.path.join(hla_cassic_basepath, tot_obj.prop_id, tot_obj.prop_id + "_" + tot_obj.obset_id)  # Generate path to HLA classic products
        elif os.path.exists(os.path.join(os.getcwd(), "hla_classic")):  # For local testing
            hla_classic_basepath = os.path.join(os.getcwd(), "hla_classic")
            hla_classic_path = hla_classic_basepath
        else:
            return  # bail out if HLA classic path can't be found.
        for filt_obj in tot_obj.fdp_list:
            hap_imgname = filt_obj.drizzle_filename
            hla_imgname = glob.glob("{}/{}{}_dr*.fits".format(hla_classic_path, filt_obj.basename, filt_obj.filters))[0]
            if not os.path.exists(hap_imgname) or not os.path.exists(hla_imgname):  # Skip filter if one or both of the images can't be found
                continue
            for hap_sourcelist_name in [filt_obj.point_cat_filename, filt_obj.segment_cat_filename]:
                if hap_sourcelist_name.endswith("point-cat.ecsv"):
                    hla_classic_cat_type = "dao"
                    plotfile_prefix = filt_obj.product_basename + "_point"
                else:
                    hla_classic_cat_type = "sex"
                    plotfile_prefix = filt_obj.product_basename + "_segment"
                if hla_classic_basepath and hla_build_ver and os.path.exists(hla_classic_basepath):
                    hla_sourcelist_name = "{}/logs/{}{}_{}phot.txt".format(hla_classic_path, filt_obj.basename,
                                                                           filt_obj.filters, hla_classic_cat_type)
                else:
                    hla_sourcelist_name = "{}/{}{}_{}phot.txt".format(hla_classic_path, filt_obj.basename,
                                                                      filt_obj.filters, hla_classic_cat_type)
                if not os.path.exists(hap_sourcelist_name) or not os.path.exists(hla_sourcelist_name):  # Skip catalog type if one or both of the catalogs can't be found
                    continue

                # convert HLA Classic RA and Dec values to HAP reference frame so the RA and Dec comparisons are correct
                updated_hla_sourcelist_name = correct_hla_classic_ra_dec(hla_sourcelist_name, hap_imgname, hla_classic_cat_type, log_level)
                log.info("HAP image:                   {}".format(os.path.basename(hap_imgname)))
                log.info("HLA Classic image:           {}".format(os.path.basename(hla_imgname)))
                log.info("HAP catalog:                 {}".format(os.path.basename(hap_sourcelist_name)))
                log.info("HLA Classic catalog:         {}".format(os.path.basename(updated_hla_sourcelist_name)))

                # once all file exist checks are passed, execute sourcelist comparision
                return_status = compare_sourcelists.comparesourcelists(slNames=[updated_hla_sourcelist_name,
                                                                       hap_sourcelist_name],
                                                                       imgNames=[hla_imgname, hap_imgname],
                                                                       good_flag_sum=255,
                                                                       plotGen="file",
                                                                       plotfile_prefix=plotfile_prefix,
                                                                       output_json_filename=hap_sourcelist_name.replace(".ecsv", "_compare_sourcelists.json"),
                                                                       verbose=True,
                                                                       log_level=log_level,
                                                                       debugMode=diagnostic_mode)
                combo_comp_pdf_filename = "{}_comparision_plots.pdf".format(plotfile_prefix)
                if os.path.exists(combo_comp_pdf_filename):
                    combo_comp_pdf_list.append(combo_comp_pdf_filename)
        if len(combo_comp_pdf_list) > 0:  # combine all plots generated by compare_sourcelists.py for this total object into a single pdf file
            total_combo_comp_pdf_filename = "{}_comparision_plots.pdf".format(tot_obj.drizzle_filename[:-9].replace("_total", ""))
            compare_sourcelists.pdf_merger(total_combo_comp_pdf_filename, combo_comp_pdf_list)
            log.info("Sourcelist comparison plots saved to file {}.".format(total_combo_comp_pdf_filename))
# ----------------------------------------------------------------------------------------------------------------------


def run_sourcelist_flagging(filter_product_obj, filter_product_catalogs, log_level, diagnostic_mode=False):
    """
    Super-basic and profoundly inelegant interface to hla_flag_filter.py.

    Execute haputils.hla_flag_filter.run_source_list_flaging() to populate the "Flags" column in the catalog tables
    generated by HAPcatalogs.measure().

    Parameters
    ----------
    filter_product_obj : drizzlepac.hlautils.product.FilterProduct object
        object containing all the relevant info for the drizzled filter product

    filter_product_catalogs : drizzlepac.hlautils.catalog_utils.HAPCatalogs object
        drizzled filter product catalog object

    log_level : int
        The desired level of verboseness in the log statements displayed on the screen and written to the .log file.

    diagnostic_mode : Boolean, optional.
        create intermediate diagnostic files? Default value is False.

    Returns
    -------
    filter_product_catalogs : drizzlepac.hlautils.catalog_utils.HAPCatalogs object
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

        filter_product_catalogs.catalogs[cat_type].source_cat = hla_flag_filter.run_source_list_flagging(drizzled_image,
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

    return filter_product_catalogs


def _get_envvar_switch(envvar_name):
    """
    This private routine interprets the environment variable, SVM_QUALITY_TESTING,
    if specified.  NOTE: This is a copy of the routine in runastrodriz.py.  This
    code should be put in a common place.
    """
    if envvar_name in os.environ:
        val = os.environ[envvar_name].lower()
        if val not in envvar_bool_dict:
            msg = "ERROR: invalid value for {}.".format(envvar_name)
            msg += "  \n    Valid Values: on, off, yes, no, true, false"
            raise ValueError(msg)
        switch_val = envvar_bool_dict[val]
    else:
        switch_val = None

    return switch_val
