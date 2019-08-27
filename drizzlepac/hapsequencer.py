#!/usr/bin/env python

""" This script defines the HST Advanced Products (HAP) generation portion of the
    calibration pipeline.  This portion of the pipeline produces mosaic and catalog
    products. This script provides similar functionality as compared to the Hubble
    Legacy Archive (HLA) pipeline in that it provides the overall sequence of 
    the processing.
"""
import datetime
import glob
import os
import sys
import traceback

from astropy.io import fits
import drizzlepac
from drizzlepac import alignimages
from drizzlepac import astrodrizzle
from drizzlepac import wcs_functions
from drizzlepac.hlautils.catalog_utils import HAPCatalogs
from drizzlepac.hlautils import config_utils
from drizzlepac.hlautils import poller_utils
from drizzlepac.hlautils import processing_utils as proc_utils
from drizzlepac.hlautils import sourcelist_generation
from stsci.tools import logutil

__taskname__ = 'hapsequencer'
log = logutil.create_logger('hapsequencer', level=logutil.logging.INFO, stream=sys.stdout)

__version__ = 0.1
__version_date__ = '19-Mar-2019'

# ----------------------------------------------------------------------------------------------------------------------
# set up instrument/detector-specific params
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


def create_catalog_products(total_list, debug=False, phot_mode='both'):
    """This subroutine utilizes hlautils/catalog_utils module to produce photometric sourcelists for the specified
    total drizzle product and it's associated child filter products.

    Parameters
    ----------
    total_list : drizzlepac.hlautils.Product.TotalProduct
        total drizzle product that will be processed by catalog_utils. catalog_utils will also create photometric
        sourcelists for the child filter products of this total product.

    debug : bool, optional
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
    for total_product_obj in total_list:
        # determine total product filename
        if os.path.exists(total_product_obj.product_basename+"_drc.fits"):
            total_product_name = total_product_obj.product_basename + "_drc.fits"
        else:
            total_product_name = total_product_obj.product_basename + "_drz.fits"

        # Instantiate filter catalog product object
        total_product_catalogs = HAPCatalogs(total_product_name,
                                             total_product_obj.pars.get_pars('catalog generation'),
                                             types=phot_mode,
                                             debug=debug)

        # Identify sources to be measured by filter photometry step
        total_product_catalogs.identify()

        # write out list(s) of identified sources
        total_product_catalogs.write()

        # append total product catalogs to list
        if phot_mode in ['aperture', 'both']:
            product_list.append(total_product_obj.point_cat_filename)
        if phot_mode in ['segment', 'both']:
            product_list.append(total_product_obj.segment_cat_filename)

        # build dictionary of total_product_catalogs.catalogs[*].sources to use for
        # filter photometric catalog generation
        sources_dict = {}
        for cat_type in total_product_catalogs.catalogs.keys():
            sources_dict[cat_type] = {}
            sources_dict[cat_type]['sources'] = total_product_catalogs.catalogs[cat_type].sources
            if cat_type == "segment":
                sources_dict['segment']['kernel'] = total_product_catalogs.catalogs['segment'].kernel

        for filter_product_obj in total_product_obj.fdp_list:
            # determine filter product filename
            if os.path.exists(filter_product_obj.product_basename + "_drc.fits"):
                filter_product_name = filter_product_obj.product_basename + "_drc.fits"
            else:
                filter_product_name = filter_product_obj.product_basename + "_drz.fits"

            # Instantiate filter catalog product object
            filter_product_catalogs = HAPCatalogs(filter_product_name,
                                                  filter_product_obj.pars.get_pars('catalog generation'),
                                                  types=phot_mode,
                                                  debug=debug,
                                                  tp_sources=sources_dict)
            # Perform photometry
            filter_product_catalogs.measure()

            # Write out photometric catalog(s)
            filter_product_catalogs.write()

            # append filter product catalogs to list
            if phot_mode in ['aperture', 'both']:
                product_list.append(filter_product_obj.point_cat_filename)
            if phot_mode in ['segment', 'both']:
                product_list.append(filter_product_obj.segment_cat_filename)

    return product_list


# ----------------------------------------------------------------------------------------------------------------------

def create_drizzle_products(obs_info_dict, total_list, meta_wcs):
    """
    Run astrodrizzle to produce products specified in the total_list.

    Parameters
    ----------
    obs_info_dict : dictionary
        Dictionary containing all relevant information required to process the dataset.

    total_list: list
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

    # TODO: This will be updated with the new configuration management.
    cfgfile_path = os.path.join(os.path.dirname(__file__), "pars")
    total_config_obj = '{}{}astrodrizzle_total_hap.cfg'.format(cfgfile_path, os.path.sep)
    filt_config_obj = '{}{}astrodrizzle_filter_hap.cfg'.format(cfgfile_path, os.path.sep)
    expo_config_obj = '{}{}astrodrizzle_single_hap.cfg'.format(cfgfile_path, os.path.sep)
    product_list = []

    # Create drizzle-combined total detection image using the meta_wcs as the reference output
    # for each instrument/detector combination (single instrument, multiple detectors)
    for total_obj in total_list:
        log.info("~" * 118)

        # TODO: Use WCS
        # log.info("CREATE TOTAL DRIZZLE-COMBINED IMAGE: {}\n".format(total_obj.drizzle_filename))
        # total_obj.wcs_drizzle_product(meta_wcs, total_config_obj)
        # product_list.append(total_obj.drizzle_filename)

        ref_total_combined_image = total_obj.ref_drizzle_filename
        log.info("CREATE TEMPORARY REFERENCE TOTAL DRIZZLED IMAGE: {}\n".format(ref_total_combined_image))
        adriz_in_list = [element.full_filename for element in total_obj.edp_list]
        log.info("Ref total combined image. {} {}".format(ref_total_combined_image, adriz_in_list))
        astrodrizzle.AstroDrizzle(input=adriz_in_list, output=ref_total_combined_image,
                                  configobj=total_config_obj)
        log.info("Finished creating TEMP REFERENCE TOTAL DRIZZLED IMAGE: {}\n".format(ref_total_combined_image))

        # Extract shape of ref_total_combined_image for explicit use in AstroDrizzle for all other products.
        rtci = fits.open(ref_total_combined_image)
        total_shape = rtci[('sci', 1)].data.shape
        rtci.close()

        # Create drizzle-combined filter image
        for filt_obj in total_obj.fdp_list:
            log.info("~" * 118)

            # TODO: Use WCS
            # filt_obj.wcs_drizzle_product(meta_wcs, filt_config_obj)
            log.info("CREATE DRIZZLE-COMBINED FILTER IMAGE: {}\n".format(filt_obj.drizzle_filename))
            filt_obj.image_drizzle_product(ref_total_combined_image, total_shape, filt_config_obj)
            product_list.append(filt_obj.drizzle_filename)

            # Create individual single drizzled images
            for exposure_obj in filt_obj.edp_list:
                log.info("~" * 118)

                # TODO: Use WCS
                # exposure_obj.wcs_drizzle_product(meta_wcs, expo_config_obj)
                log.info("CREATE SINGLE DRIZZLED IMAGE: {}".format(exposure_obj.drizzle_filename))
                exposure_obj.image_drizzle_product(ref_total_combined_image, total_shape, expo_config_obj)
                product_list.append(exposure_obj.drizzle_filename)

        # TODO: Use WCS
        # total_obj.wcs_drizzle_product(meta_wcs, total_config_obj)

        log.info("CREATE DRIZZLE-COMBINED TOTAL IMAGE: {}\n".format(total_obj.drizzle_filename))
        total_obj.image_drizzle_product(ref_total_combined_image, total_shape, total_config_obj)
        product_list.append(total_obj.drizzle_filename)

    # Ensure that all drizzled products have headers that are to specification
    try:
        log.info("Updating these drizzle products for CAOM compatibility:")
        for filename in product_list:
            log.info("    {}".format(filename))
            proc_utils.refine_product_headers(filename, total_list)
    except Exception:
        print("Trouble updating drizzle products for CAOM.")
        # TODO: Maybe these next two lines are just temporary for debugging
        exc_type, exc_value, exc_tb = sys.exc_info()
        traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)

    # Remove rules files copied to the current working directory
    for rules_filename in glob.glob("*_header_hla.rules"):
        log.info("Removed rules file {}".format(rules_filename))
        os.remove(rules_filename)

    # Return product list for creation of pipeline manifest file
    return product_list

# ----------------------------------------------------------------------------------------------------------------------


def run_hla_processing(input_filename, result=None, debug=False, use_defaults_configs=True,
                       input_custom_pars_file=None, output_custom_pars_file=None, phot_mode='both'):
    # This routine needs to return an exit code, return_value, for use by the calling
    # Condor/OWL workflow code: 0 (zero) for success, 1 for error condition
    return_value = 0

    starting_dt = datetime.datetime.now()
    log.info("Run start time: {}".format(str(starting_dt)))
    product_list = []
    try:
        # Parse the poller file and generate the the obs_info_dict, as well as the total detection 
        # product lists which contain the ExposureProduct, FilterProduct, and TotalProduct objects
        log.info("Parse the poller and determine what exposures need to be combined into separate products")
        obs_info_dict, total_list = poller_utils.interpret_obset_input(input_filename)

        # Generate the name for the manifest file which is for the entire visit.  It is fine
        # to use only one of the Total Products to generate the manifest name as it is not
        # dependent on the detector.
        # instrument_programID_obsetID_manifest.txt (e.g.,wfc3_b46_06_manifest.txt)
        manifest_name = total_list[0].manifest_name

        # Set up all pipeline parameter values that will be used later on in the run.
        for total_item in total_list:
            total_item.pars = config_utils.HapConfig(total_item,
                                                     use_defaults=use_defaults_configs,
                                                     input_custom_pars_file=input_custom_pars_file,
                                                     output_custom_pars_file=output_custom_pars_file)

            for filter_item in total_item.fdp_list:
                filter_item.pars = config_utils.HapConfig(filter_item,
                                                          use_defaults=use_defaults_configs,
                                                          input_custom_pars_file=input_custom_pars_file,
                                                          output_custom_pars_file=output_custom_pars_file)
            for expo_item in total_item.edp_list:
                expo_item.pars = config_utils.HapConfig(expo_item,
                                                        use_defaults=use_defaults_configs,
                                                        input_custom_pars_file=input_custom_pars_file,
                                                        output_custom_pars_file=output_custom_pars_file)

        # Run alignimages.py on images on a filter-by-filter basis.
        # Process each filter object which contains a list of exposure objects/products,
        # regardless of detector.
        log.info("Run alignimages.py on images on a filter-by-filter basis.")
        exposure_filenames = []
        for tot_obj in total_list:
            for filt_obj in tot_obj.fdp_list:
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
                    log.info("Alignimages step skipped.")

        # Run meta wcs code to get common WCS for all images in this obset_id, regardless of detector.
        # FIX (1) Intended for this to be a method of TotalProduct, but it should be
        # associated with all the exposures really used in the alignment (the "as built")
        # as is done here.
        # This function used based upon WH analysis but make sure to set
        # the size of the output image. This comment is related to the previously mentioned issue.
        # This produced incompatible results.  Perhaps accessing wrong dimension information.
        """
        log.info("Run make_mosaic_wcs to create a common WCS for all images aligned in the previous step.")
        log.info("The following images will be used: ")
        for imgname in exposure_filenames:
            log.info("{}".format(imgname))
        if exposure_filenames:
            meta_wcs = wcs_functions.make_mosaic_wcs(exposure_filenames)
        """
        # Not using meta_wcs at this time
        meta_wcs = []

        # Run AstroDrizzle to produce drizzle-combined products
        log.info("Create drizzled imagery products")
        driz_list = create_drizzle_products(obs_info_dict, total_list, meta_wcs)
        product_list += driz_list

        # Create source catalogs from newly defined products (HLA-204)
        log.info("Create source catalog from newly defined product")
        if 'total detection product 00' in obs_info_dict.keys():
            catalog_list = create_catalog_products(total_list, debug=debug, phot_mode=phot_mode)
            product_list += catalog_list
        else:
            print("Sourcelist generation step skipped.")
        """
        # 8: (OPTIONAL) Determine whether there are any problems with alignment or photometry of product
        log.info("8: (TODO) (OPTIONAL) Determine whether there are any problems with alignment or photometry"
                 "of product")
        # TODO: QUALITY CONTROL SUBROUTINE CALL GOES HERE.
        """

        # Write out manifest file listing all products generated during processing
        log.info("Creating manifest file {}".format(manifest_name))
        log.info("  Manifest contains the names of products generated during processing.")
        with open(manifest_name, mode='w') as catfile:
            [catfile.write("{}\n".format(name)) for name in product_list]

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

# ----------------------------------------------------------------------------------------------------------------------
# TODO: COMMENT OUT/REMOVE BELOW COMMAND-LINE INTERFACE PRIOR TO PULL REQUEST

# if __name__ == '__main__':
#     os.system("rm -f *.*")
#     os.system("cp orig/* .")
#     out_pars_file = sys.argv[1].replace(".out", "_config.json")
#
#     x = run_hla_processing(sys.argv[1], debug=True, output_custom_pars_file=out_pars_file, phot_mode='aperture')
#     print("\a")
