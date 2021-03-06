#!/usr/bin/env python

""" This script defines the HST Advanced Products (HAP) Multi-visit (MVM) generation
    portion of the calibration pipeline.  This portion of the pipeline produces mosaic
    products. This script provides similar functionality as compared to the Hubble
    Legacy Archive (HLA) pipeline in that it provides the overall sequence of
    the processing.

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

from drizzlepac.haputils import config_utils
from drizzlepac.haputils import poller_utils
from drizzlepac.haputils import product
from drizzlepac.haputils import processing_utils as proc_utils
from drizzlepac.haputils import svm_quality_analysis as svm_qa

from stsci.tools import logutil
from stwcs import wcsutil


__taskname__ = 'hapmultisequencer'
MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)
__version__ = 0.1
__version_date__ = '01-May-2020'

# Environment variable which controls the quality assurance testing
# for the Single Visit Mosaic processing.
envvar_bool_dict = {'off': False, 'on': True, 'no': False, 'yes': True, 'false': False, 'true': True}
envvar_qa_svm = "SVM_QUALITY_TESTING"

# --------------------------------------------------------------------------------------------------------------


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
    # Get rules files
    rules_files = {}

    log.info("Processing with astrodrizzle version {}".format(drizzlepac.astrodrizzle.__version__))
    # Get rules files
    for imgname in glob.glob("*skycell*fl?.fits"):
        rules_files[imgname] = proc_utils.get_rules_file(imgname)

    # Keep track of all the products created for the output manifest
    product_list = []

    # For each detector (as the total detection product are instrument- and detector-specific),
    # create the drizzle-combined filtered image, the drizzled exposure (aka single) images,
    # and finally the drizzle-combined total detection image.
    for filt_obj in total_obj_list:
        filt_obj.rules_file = rules_files[filt_obj.edp_list[0].full_filename]

        log.info("~" * 118)
        # Get the common WCS for all images which are part of a total detection product,
        # where the total detection product is detector-dependent.
        meta_wcs = filt_obj.generate_metawcs()

        log.info("CREATE DRIZZLE-COMBINED FILTER IMAGE: {}\n".format(filt_obj.drizzle_filename))
        filt_obj.wcs_drizzle_product(meta_wcs)
        product_list.append(filt_obj.drizzle_filename)
        product_list.append(filt_obj.trl_filename)

        # Add individual single input images with updated WCS headers to manifest
        for exposure_obj in filt_obj.edp_list:
            product_list.append(exposure_obj.full_filename)

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
    for filt_obj in total_obj_list:
        filt_obj = poller_utils.add_primary_fits_header_as_attr(filt_obj)

    # Return product list for creation of pipeline manifest file
    return product_list

# ----------------------------------------------------------------------------------------------------------------------


def run_mvm_processing(input_filename, diagnostic_mode=False, use_defaults_configs=True,
                       input_custom_pars_file=None, output_custom_pars_file=None, phot_mode="both",
                       log_level=logutil.logging.INFO):
    """
    Run the HST Advanced Products (HAP) generation code.  This routine is the sequencer or
    controller which invokes the high-level functionality to process the multi-visit data.

    Parameters
    ----------
    input_filename: string
        The 'poller file' where each line contains information regarding an exposures considered
        part of the multi-visit.

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
        Which algorithm should be used to generate the sourcelists? 'aperture' for aperture (point) photometry;
        'segment' for isophotal photometry; 'both' for both 'segment' and 'aperture'. Default value is 'both'.

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
        logname = 'mvm_process.log'

    # Initialize total trailer filename as temp logname
    logging.basicConfig(filename=logname, format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)
    # start processing
    starting_dt = datetime.datetime.now()
    log.info("Run start time: {}".format(str(starting_dt)))
    total_obj_list = []
    product_list = []
    try:
        # Parse the MVM poller file and generate the the obs_info_dict, as well as the total detection
        # product lists which contain the ExposureProduct, FilterProduct, and TotalProduct objects
        # A poller file contains visit data for a single instrument.  The TotalProduct discriminant
        # is the detector.  A TotalProduct object is comprised of FilterProducts and ExposureProducts
        # where its FilterProduct is distinguished by the filter in use, and the ExposureProduct
        # is the atomic exposure data.
        log.info("Parse the poller and determine what exposures need to be combined into separate products.\n")
        obs_info_dict, total_obj_list = poller_utils.interpret_mvm_input(input_filename, log_level,
                                                                         layer_method='all')

        # Generate the name for the manifest file which is for the entire multi-visit.  It is fine
        # to use only one of the Total Products to generate the manifest name as the name is not
        # dependent on the detector.
        # Example: instrument_programID_obsetID_manifest.txt (e.g.,wfc3_b46_06_manifest.txt)
        manifest_name = total_obj_list[0].manifest_name
        log.info("\nGenerate the manifest name for this multi-visit.")
        log.info("The manifest will contain the names of all the output products.")

        # The product_list is a list of all the output products which will be put into the manifest file
        product_list = []

        # Update the SkyCellProduct objects with their associated configuration information.
        for filter_item in total_obj_list:
            filter_item.generate_metawcs()
            filter_item.generate_footprint_mask()
            log.info("Preparing configuration parameter values for filter product {}".format(filter_item.drizzle_filename))
            filter_item.configobj_pars = config_utils.HapConfig(filter_item,
                                                                log_level=log_level,
                                                                use_defaults=use_defaults_configs,
                                                                input_custom_pars_file=input_custom_pars_file,
                                                                output_custom_pars_file=output_custom_pars_file)

        log.info("The configuration parameters have been read and applied to the drizzle objects.")

        reference_catalog = run_align_to_gaia(total_obj_list, log_level=log_level, diagnostic_mode=diagnostic_mode)
        if reference_catalog:
            product_list += [reference_catalog]

        # Run AstroDrizzle to produce drizzle-combined products
        log.info("\n{}: Create drizzled imagery products.".format(str(datetime.datetime.now())))
        driz_list = create_drizzle_products(total_obj_list)
        product_list += driz_list

        # Store total_obj_list to a pickle file to speed up development
        if False:
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

            # Number of sources in Point and Segment catalogs
            total_catalog_list = [i for i in catalog_list if 'total' in i]
            fits_list = [i for i in driz_list if 'fits' in i]
            total_drizzle_list = [i for i in fits_list if 'total' in i]
            svm_qa.compare_num_sources(total_catalog_list, total_drizzle_list, log_level=log_level)

            # Get point/segment cross-match RA/Dec statistics
            for filter_obj in total_obj_list:
                svm_qa.compare_ra_dec_crossmatches(filter_obj, log_level=log_level)

            # Identify the number of GAIA sources in final product footprints
            for filter_obj in total_obj_list:
                svm_qa.find_gaia_sources(filter_obj, log_level=log_level)

            # Photometry of cross-matched sources in Point and Segment catalogs for Filter products
            tot_len = len(total_obj_list)
            filter_drizzle_list = []
            temp_list = []
            for tot in total_obj_list:
                temp_list = [x.drizzle_filename for x in tot.fdp_list]
                filter_drizzle_list.extend(temp_list)
            svm_qa.compare_photometry(filter_drizzle_list, log_level=log_level)

        # 9: Compare results to HLA classic counterparts (if possible)
        # if diagnostic_mode:
            # run_sourcelist_comparison(total_obj_list, diagnostic_mode=diagnostic_mode, log_level=log_level)

        # Insure manifest file does not contain duplicate entries
        # Use of numpy.unique preserves the order of the entries in the product list
        product_list = np.unique(product_list).tolist()
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
        catalog_list = [gaia_obj.configobj_pars.pars['alignment'].pars_multidict['all']['run_align']['catalog_list'][0]]  # For now, just pass in a single catalog name as list
        align_table, filt_exposures = gaia_obj.align_to_gaia(catalog_list=catalog_list,
                                                             output=diagnostic_mode,
                                                             fit_label='MVM')

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
