#!/usr/bin/env python

"""This script duplicates the functionality of the HLA pipeline.

"""
import datetime
import glob
import os
import pickle
import sys
import traceback
import shutil

from astropy.io import fits
import drizzlepac
from drizzlepac import alignimages
from drizzlepac import astrodrizzle
from drizzlepac import wcs_functions
from drizzlepac.hlautils import poller_utils
from drizzlepac.hlautils import processing_utils as dpu
from drizzlepac.hlautils import sourcelist_generation
from stsci.tools import logutil

from drizzlepac.hlautils import product

__taskname__ = 'runhlaprocessing'
log = logutil.create_logger('runhlaprocessing', level=logutil.logging.INFO, stream=sys.stdout)

__version__ = 0.1
__version_date__ = '19-Mar-2019'

# ----------------------------------------------------------------------------------------------------------------------
# set up instrument/detector-specific params.
# Parameters imported from the following HLA classic parameter files:
# - https://grit.stsci.edu/HLA/hla/tree/master/software/trunk/HLApipeline/HLApipe/param/parameter_acs_hrc.cfg
# - https://grit.stsci.edu/HLA/hla/tree/master/software/trunk/HLApipeline/HLApipe/param/parameter_acs_sbc.cfg
# - https://grit.stsci.edu/HLA/hla/tree/master/software/trunk/HLApipeline/HLApipe/param/parameter_acs_wfc.cfg
# - https://grit.stsci.edu/HLA/hla/tree/master/software/trunk/HLApipeline/HLApipe/param/parameter_wfc3_ir.cfg
# - https://grit.stsci.edu/HLA/hla/tree/master/software/trunk/HLApipeline/HLApipe/param/parameter_wfc3_uvis.cfg

param_dict = {
    "ACS HRC": {
        "astrodrizzle": {
            "SCALE": 0.025,
            "PIXFRAC": 1.0,
            "KERNEL": "square",
            "OUTNX": None,
            "OUTNY": None,
            "ROT": 0.0,
            "BITS": 256},
        "ci filter": {
            "ci_daolower_limit": 0.9,
            "ci_daoupper_limit": 1.6,
            "ci_selower_limit": 0.9,
            "ci_seupper_limit": 1.6},
        "dao": {
            "TWEAK_FWHMPSF": 0.073,
            "TWEAK_THRESHOLD": 3.0,
            "aperture_1": 0.03,
            "aperture_2": 0.125,
            "bthresh": 5.0},
        "sourcex": {
            "fwhm": 0.073,
            "thresh": 1.4,
            "bthresh": 5.0,
            "source_box": 7},
        "swarm filter": {
            "upper_epp_limit": 70000.,
            "lower_epp_limit": 2000.,
            "eppsky_limit": 1000.,
            "swarm_thresh": 1.,
            "clip_radius_list": [120.0, 100.0, 80.0, 60.0, 40.0, 20.0, 10.0, 5.0, 2.0, 0.0],
            "scale_factor_list": [0.0, 1.778106e-05, 3.821292e-05, 9.017166e-05, 2.725184e-04, 1.269197e-03, 7.007126e-03, 3.839166e-02, 2.553349e-01, 1.000000e+00],
            "proximity_binary": "no"}},
    "ACS SBC": {
        "astrodrizzle": {
            "SCALE": 0.03,
            "PIXFRAC": 1.0,
            "KERNEL": "square",
            "OUTNX": None,
            "OUTNY": None,
            "ROT": 0.0,
            "BITS": 256},
        "ci filter": {
            "ci_daolower_limit": 0.15,
            "ci_daoupper_limit": 0.45,
            "ci_selower_limit": 0.15,
            "ci_seupper_limit": 0.45},
        "dao": {
            "TWEAK_FWHMPSF": 0.065,
            "TWEAK_THRESHOLD": 3.0,
            "aperture_1": 0.07,
            "aperture_2": 0.125,
            "bthresh": 5.0},
        "sourcex": {
            "fwhm": 0.065,
            "thresh": 1.4,
            "bthresh": 5.0,
            "source_box": 7},
        "swarm filter": {
            "upper_epp_limit": 70000.,
            "lower_epp_limit": 2000.,
            "eppsky_limit": 1000.,
            "swarm_thresh": 1.,
            "clip_radius_list": [120.0, 100.0, 80.0, 60.0, 40.0, 20.0, 10.0, 5.0, 2.0, 0.0],
            "scale_factor_list": [0.0, 1.778106e-05, 3.821292e-05, 9.017166e-05, 2.725184e-04, 1.269197e-03, 7.007126e-03, 3.839166e-02, 2.553349e-01, 1.000000e+00],
            "proximity_binary": "no"}},
    "ACS WFC": {
        "astrodrizzle": {
            "SCALE": 0.05,
            "PIXFRAC": 1.0,
            "KERNEL": "square",
            "OUTNX": None,
            "OUTNY": None,
            "ROT": 0.0,
            "BITS": 256},
        "ci filter": {
            "ci_daolower_limit": 0.9,
            "ci_daoupper_limit": 1.23,
            "ci_selower_limit": 0.9,
            "ci_seupper_limit": 1.23},
        "dao": {
            "TWEAK_FWHMPSF": 0.076,
            "TWEAK_THRESHOLD": 3.0,
            "aperture_1": 0.05,  # update from 0.15
            "aperture_2": 0.15,  # update from 0.25
            "bthresh": 5.0},
        "sourcex": {
            "fwhm": 0.13,
            "thresh": -1.2,
            "bthresh": 5.0,
            "source_box": 5},
        "swarm filter": {
            "upper_epp_limit": 70000.,
            "lower_epp_limit": 2000.,
            "eppsky_limit": 1000.,
            "swarm_thresh": 1.,
            "clip_radius_list": [120., 100., 80., 60., 40., 30., 20., 10., 5., 2., 0.],
            "scale_factor_list": [0.0, 0.000000e+00, 6.498530e-06, 3.687270e-05, 1.412972e-04, 3.151877e-04, 1.023391e-03, 3.134859e-03, 2.602436e-02, 1.820539e-01, 1.000000e+00],
            "proximity_binary": "no"}},
    "WFC3 IR": {
        "astrodrizzle": {
            "SCALE": 0.09,
            "PIXFRAC": 1.0,
            "KERNEL": "square",
            "OUTNX": None,
            "OUTNY": None,
            "ROT": 0.0,
            "BITS": 768},
        "ci filter": {
            "ci_daolower_limit": 0.25,
            "ci_daoupper_limit": 0.55,
            "ci_selower_limit": 0.25,
            "ci_seupper_limit": 0.55},
        "dao": {
            "TWEAK_FWHMPSF": 0.14,
            "TWEAK_THRESHOLD": 3.0,
            "aperture_1": 0.15,
            "aperture_2": 0.45,
            "bthresh": 5.0},
        "sourcex": {
            "fwhm": 0.14,
            "thresh": 1.4,
            "bthresh": 5.0,
            "source_box": 7},
        "swarm filter": {
            "upper_epp_limit": 70000.,
            "lower_epp_limit": 2000.,
            "eppsky_limit": 100.,
            "swarm_thresh": 1.,
            "clip_radius_list": [140., 120., 100., 80., 60., 40., 20., 10., 5., 2., 0.],
            #                   x10    x10    x10   x10   x10   x10    x10   x10  x10  x2,
            "scale_factor_list": [1.5e-5, 2.3e-5, 4.e-5, 8.e-5, 2.e-4, 0.0006, 0.015, 0.05, 0.15, 0.9, 1.],
            # "scale_factor_list_orig": [1.5e-5, 2.3e-5, 4.e-5, 8.e-5, 2.e-4, 0.0006, 0.005, 0.05, 0.15, 0.9, 1.],
            "proximity_binary": "yes"}},
    "WFC3 UVIS": {
        "astrodrizzle": {
            "SCALE": 0.04,
            "PIXFRAC": 1.0,
            "KERNEL": "square",
            "OUTNX": None,
            "OUTNY": None,
            "ROT": 0.0,
            "BITS": 256},
        "ci filter": {
            "ci_daolower_limit": 0.75,
            "ci_daoupper_limit": 1.0,
            "ci_selower_limit": 0.75,
            "ci_seupper_limit": 1.0},
        "dao": {
            "TWEAK_FWHMPSF": 0.076,
            "TWEAK_THRESHOLD": 3.0,
            "aperture_1": 0.05,
            "aperture_2": 0.15,
            "bthresh": 5.0},
        "sourcex": {
            "fwhm": 0.076,
            "thresh": 1.4,
            "bthresh": 5.0,
            "source_box": 7},
        "swarm filter": {
            "upper_epp_limit": 70000.,
            "lower_epp_limit": 2000.,
            "eppsky_limit": 1000.,
            "swarm_thresh": 1.,
            "clip_radius_list": [120., 100., 80., 60., 40., 20., 10., 5., 2., 0.],
            "scale_factor_list": [2.3e-6, 4.e-6, 8.e-6, 2.e-5, 0.0005, 0.005, 0.005, 0.015, 0.45, 1.],
            # "scale_factor_list_orig": [2.3e-6, 4.e-6, 8.e-6, 2.e-5, 6.e-5, 0.0005, 0.005, 0.015, 0.45, 1.],
            "proximity_binary": "yes"}}}
# ----------------------------------------------------------------------------------------------------------------------

def create_drizzle_products(obs_info_dict, total_list, filt_list, expo_list, meta_wcs):
    """
    Run astrodrizzle to produce products specified in obs_info_dict.

    Parameters
    ----------
    obs_info_dict : dictionary
        Dictionary containing all relevant information required to process the dataset.

    RETURNS
    -------
    Nothing.
    """
    log.info("Processing with astrodrizzle version {}".format(drizzlepac.astrodrizzle.__version__))

    # Get rules files
    for imgname in glob.glob("*fl?.fits"):
        dpu.get_rules_file(imgname)

    # FIX - This will be updated with the new configuration management.
    cfgfile_path = os.path.join(os.path.dirname(__file__), "pars")
    total_config_obj = '{}{}astrodrizzle_total_hap.cfg'.format(cfgfile_path, os.path.sep)
    filt_config_obj = '{}{}astrodrizzle_filter_hap.cfg'.format(cfgfile_path, os.path.sep)
    expo_config_obj = '{}{}astrodrizzle_single_hap.cfg'.format(cfgfile_path, os.path.sep)
    product_list = []

    # Create drizzle-combined total detection image using the meta_wcs as the reference output 
    # for each instrument/detector combination (single instrument, multiple detectors)
    for total_obj in total_list:
        log.info("~" * 118)

        # Make sure to create the drizzle filename before a call to the drizzle functionality
        # FIX in pipeline_poller_utils
        _ = total_obj.create_drizzle_filename()
        log.info("CREATE TOTAL DRIZZLE-COMBINED IMAGE: {}\n".format(total_obj.drizzle_filename))
        total_obj.drizzle_product(meta_wcs, total_config_obj)
        product_list.append(total_obj.drizzle_filename)

        # Create drizzle-combined filter image using the meta_wcs as the reference output 
        print("length of total FDP list: {}".format(len(total_obj.fdp_list)))
        for filt_obj in total_obj.fdp_list:
                log.info("~" * 118)

                # Make sure to create the drizzle filename before a call to the drizzle functionality
                # FIX in pipeline_poller_utils
                _ = filt_obj.create_drizzle_filename()
                log.info("CREATE DRIZZLE-COMBINED FILTER IMAGE: {}\n".format(filt_obj.drizzle_filename))
                filt_obj.drizzle_product(meta_wcs, filt_config_obj)
                product_list.append(filt_obj.drizzle_filename)

                # Create individual singlely-drizzled images using meta_wcs
                for exposure_obj in filt_obj.edp_list:
                    log.info("~" * 118)
                    log.info("CREATE SINGLY DRIZZLED IMAGE: {}".format(exposure_obj.drizzle_filename))

                    exposure_obj.drizzle_product(meta_wcs, expo_config_obj)
                    product_list.append(exposure_obj.drizzle_filename)

    # Ensure that all drizzled products have headers that are to specification
    log.info("Updating these drizzle products for CAOM compatibility:")
    for filename in product_list:
        log.info("    {}".format(filename))
        dpu.refine_product_headers(filename, obs_info_dict)

    # Remove rules files copied to the current working directory
    for rules_filename in glob.glob("*_header_hla.rules"):
        log.info("Removed rules file {}".format(rules_filename))
        os.remove(rules_filename)

    # Return product list for creation of pipeline manifest file
    return product_list

# ----------------------------------------------------------------------------------------------------------------------


def run_hla_processing(input_filename, result=None, debug=False):
    # This routine needs to return an exit code, return_value, for use by the calling 
    # Condor/OWL workflow code: 0 (zero) for success, 1 for error condition
    return_value = 0

    starting_dt = datetime.datetime.now()
    log.info("Run start time: {}".format(str(starting_dt)))
    product_list = []
    try:
        # 1: Parse the poller file and generate the the obs_info_dict, as well as the single exposure, filter, 
        # and total detection product lists which contain the ExposureProduct, FilterProduct, and 
        # TotalProduct objects
        log.info("1: Parse the poller and determine what exposures need to be combined into separate products")
        obs_info_dict, expo_list, filt_list, total_list = poller_utils.interpret_obset_input(input_filename)

        # 2: Run alignimages.py on images on a filter-by-filter basis.
        # Process each filter object which contains a list of exposure objects/products,
        # regardless of detector.
        log.info("2: Run alignimages.py on images on a filter-by-filter basis.")
        exposure_filenames = []
        for filt_obj in filt_list:
            align_table, filt_exposures = filt_obj.align_to_gaia()

            # Report results and track the output files
            # FIX - Add info here in the case of alignment working on data that should not be aligned
            # as well as outright failure (exception vs msgs)
            if align_table:
                log.info("ALIGN_TABLE: {}".format(align_table))
                # FIX
                # os.remove("alignimages.log")  # FIX This log needs to be included in total product trailer file
                for row in align_table:
                    if row['status'] == 0:
                        log.info("Successfully aligned {} to {} astrometric frame\n".format(row['imageName'], row['catalog']))
                    # Alignment did not work for this particular image
                    # FIX - If alignment did not work for an image, it seems this exposure should
                    # be removed from the exposure lists.  TotalProduct and FilterProduct need
                    # methods to do this.
                    else:
                        log.info("Could not align {} to absolute astrometric frame\n".format(row['imageName']))

                hdrlet_list = align_table['headerletFile'].tolist()
                product_list += hdrlet_list
                exposure_filenames += filt_exposures

            else:
                log.info("{}: Alignimages step skipped.".format(obs_category))

        # 3: Run meta wcs code to get common WCS for all images in this obset_id, regardless
        # of detector.
        # FIX (1) Intended for this to be a method of TotalProduct, but it should be
        # associated with all the exposures really used in the alignment (the "as built") 
        # as is done here. 
        # This function used based upon WH analysis but make sure to set
        # the size of the output image. This comment is related to the previously mentioned issue.
        log.info("3: Run make_mosaic_wcs to create a common WCS for all images aligned in the previous step.")
        log.info("The following images will be used: ")
        for imgname in exposure_filenames:
            log.info("{}".format(imgname))
        if exposure_filenames:
            meta_wcs = wcs_functions.make_mosaic_wcs(exposure_filenames)

        # 4: Run AstroDrizzle to produce drizzle-combined products
        log.info("4: (WIP) Create drizzled imagery products")
        driz_list = create_drizzle_products(obs_info_dict, total_list, filt_list, expo_list, meta_wcs)
        product_list += driz_list

        # MDD ENDED HERE

        """
        # 7: Create source catalogs from newly defined products (HLA-204)
        log.info("7: (WIP) Create source catalog from newly defined product")
        if debug:
            pickle_filename = input_filename.replace(".out", ".pickle")
            if os.path.exists(pickle_filename):
                os.remove(pickle_filename)
            pickle_out = open(pickle_filename, "wb")
            pickle.dump([obs_info_dict, param_dict], pickle_out)
            pickle_out.close()
            print("Wrote obs_info_dict to pickle file {}".format(pickle_filename))

        if 'total detection product 00' in obs_info_dict.keys():
            catalog_list = sourcelist_generation.create_sourcelists(obs_info_dict, param_dict)
            product_list += catalog_list
        else:
            print("Sourcelist generation step skipped.")

        # 8: (OPTIONAL) Determine whether there are any problems with alignment or photometry of product
        log.info("8: (TODO) (OPTIONAL) Determine whether there are any problems with alignment or photometry"
                 "of product")
        # TODO: QUALITY CONTROL SUBROUTINE CALL GOES HERE.

        # 9: Write out manifest file listing all products generated during processing
        log.info("Creating manifest file {}".format(manifest_name))
        log.info("  Manifest contains the names of products generated during processing.")
        with open(manifest_name, mode='w') as catfile:
            [catfile.write("{}\n".format(name)) for name in product_list]

        """
        # 10: Return exit code for use by calling Condor/OWL workflow code: 0 (zero) for success, 1 for error condition
        return_value = 0
    except Exception:
        return_value = 1
        if debug:
            log.info("\a\a\a")
            exc_type, exc_value, exc_tb = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)
    finally:
        log.info('Total processing time: {} sec'.format((datetime.datetime.now() - starting_dt).total_seconds()))
        log.info("9: Return exit code for use by calling Condor/OWL workflow code: 0 (zero) for success, 1 for error "
                 "condition")
        return return_value

