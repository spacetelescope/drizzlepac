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

import drizzlepac
from drizzlepac.hlautils.catalog_utils import HAPCatalogs
from drizzlepac.hlautils import config_utils
from drizzlepac.hlautils import hla_flag_filter
from drizzlepac.hlautils import poller_utils
from drizzlepac.hlautils import processing_utils as proc_utils
from stsci.tools import logutil
from stwcs import wcsutil

__taskname__ = 'hapsequencer'
log = logutil.create_logger('hapsequencer', level=logutil.logging.INFO, stream=sys.stdout)

__version__ = 0.1
__version_date__ = '19-Mar-2019'

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
        # Instantiate filter catalog product object
        total_product_catalogs = HAPCatalogs(total_product_obj.drizzle_filename, total_product_obj.configobj_pars.get_pars('catalog generation'), total_product_obj.configobj_pars.get_pars('quality control'), types=phot_mode, debug=debug)

        # Generate an "n" exposure mask which has the image footprint set to the number
        # of exposures which constitute each pixel.
        total_product_obj.generate_footprint_mask()

        # Identify sources in the input image and delay writing the total detection
        # catalog until the photometric measurements have been done on the filter
        # images and some of the measurements can be appended to the total catalog
        total_product_catalogs.identify(mask=total_product_obj.mask)

        # Build dictionary of total_product_catalogs.catalogs[*].sources to use for
        # filter photometric catalog generation
        sources_dict = {}
        for cat_type in total_product_catalogs.catalogs.keys():
            sources_dict[cat_type] = {}
            sources_dict[cat_type]['sources'] = total_product_catalogs.catalogs[cat_type].sources
            # FIX MDD Remove?
            if cat_type == "segment":
                sources_dict['segment']['kernel'] = total_product_catalogs.catalogs['segment'].kernel

        for filter_product_obj in total_product_obj.fdp_list:

            # Instantiate filter catalog product object
            filter_product_catalogs = HAPCatalogs(filter_product_obj.drizzle_filename, total_product_obj.configobj_pars.get_pars('catalog generation'), total_product_obj.configobj_pars.get_pars('quality control'), types=phot_mode, debug=debug, tp_sources=sources_dict)
            # Perform photometry
            filter_name = filter_product_obj.filters
            filter_product_catalogs.measure(filter_name)

            filter_product_catalogs = run_sourcelist_flagging(filter_product_obj, filter_product_catalogs, debug)

            # Replace zero-value total-product catalog 'Flags' column values with meaningful filter-product catalog
            # 'Flags' column values
            filter_product_catalogs.catalogs['aperture'].subset_filter_source_cat[
                'Flags_{}'.format(filter_product_obj.filters)] = \
            filter_product_catalogs.catalogs['aperture'].source_cat['Flags']

            # Write out photometric (filter) catalog(s)
            filter_product_catalogs.write()

            # Load a dictionary with a subset table for each catalog
            subset_columns_dict = {}
            for cat_type in filter_product_catalogs.catalogs.keys():
                subset_columns_dict[cat_type] = {}
                subset_columns_dict[cat_type]['subset'] = filter_product_catalogs.catalogs[cat_type].subset_filter_source_cat

            # ...append the new columns to the total detection project catalog.
            total_product_catalogs.combine(subset_columns_dict)

            # append filter product catalogs to list
            if phot_mode in ['aperture', 'both']:
                product_list.append(filter_product_obj.point_cat_filename)
            if phot_mode in ['segment', 'both']:
                product_list.append(filter_product_obj.segment_cat_filename)

        # write out list(s) of identified sources
        total_product_catalogs.write()

        # append total product catalogs to manifest list
        if phot_mode in ['aperture', 'both']:
            product_list.append(total_product_obj.point_cat_filename)
        if phot_mode in ['segment', 'both']:
            product_list.append(total_product_obj.segment_cat_filename)


    return product_list


# ----------------------------------------------------------------------------------------------------------------------

def create_drizzle_products(total_list):
    """
    Run astrodrizzle to produce products specified in the total_list.

    Parameters
    ----------
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

    # Keep track of all the products created for the output manifest
    product_list = []

    # For each detector (as the total detection product are instrument- and detector-specific),
    # create the drizzle-combined filtered image, the drizzled exposure (aka single) images,
    # and finally the drizzle-combined total detection image.
    for total_obj in total_list:
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
        fits_files = [file for file in product_list if "fits" in file]
        for filename in fits_files:
            log.info("    {}".format(filename))
            proc_utils.refine_product_headers(filename, total_list)
    except Exception:
        print("Trouble updating drizzle products for CAOM.")
        exc_type, exc_value, exc_tb = sys.exc_info()
        traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)

    # Remove rules files copied to the current working directory
    for rules_filename in glob.glob("*_header_hla.rules"):
        log.info("Removed rules file {}".format(rules_filename))
        os.remove(rules_filename)

    # Return product list for creation of pipeline manifest file
    return product_list

# ----------------------------------------------------------------------------------------------------------------------


def run_hap_processing(input_filename, debug=False, use_defaults_configs=True,
                       input_custom_pars_file=None, output_custom_pars_file=None, phot_mode="both"):
    """
    Run the HST Advanced Products (HAP) generation code.  This routine is the sequencer or
    controller which invokes the high-level functionality to process the single visit data.

    Parameters
    ----------
    input_filename: string
        The 'poller file' where each line contains information regarding an exposures taken
        during a single visit.

    debug : bool, optional
        Allows printing of additional diagnostic information to the log.  Also, can turn on
        creation and use of pickled information.

    use_default_configs: bool, optional
        If True, use the configuration parameters in the 'default' portion of the configuration
        JSON files.  If False, use the configuration parameters in the "parameters" portion of
        the file.  The default is True.

    input_custom_pars_file: string, optional
        Represents a fully specified input filename of a configuration JSON file which has been
        customized for specialized processing.  This file should contain ALL the input parameters
        necessary for processing.  If there is a filename present for this parameter, the
        'use_default_configs' parameter is ignored. The default is None.

    output_custom_pars_file: string, optional
        Fully specified output filename which contains all the configuration parameters
        available during the processing session.  The default is None.

    phot_mode : str, optional
        Which algorithm should be used to generate the sourcelists? 'aperture' for aperture photometry;
        'segment' for segment map photometry; 'both' for both 'segment' and 'aperture'. Default value is 'both'.


    RETURNS
    -------
    return_value: integer
        A return exit code used by the calling Condor/OWL workflow code: 0 (zero) for success, 1 for error

    """
    # This routine needs to return an exit code, return_value, for use by the calling
    # Condor/OWL workflow code: 0 (zero) for success, 1 for error condition
    return_value = 0

    starting_dt = datetime.datetime.now()
    log.info("Run start time: {}".format(str(starting_dt)))
    product_list = []
    try:
        # Parse the poller file and generate the the obs_info_dict, as well as the total detection
        # product lists which contain the ExposureProduct, FilterProduct, and TotalProduct objects
        # A poller file contains visit data for a single instrument.  The TotalProduct discriminant
        # is the detector.  A TotalProduct object is comprised of FilterProducts and ExposureProducts
        # where its FilterProduct is distinguished by the filter in use, and the ExposureProduct
        # is the atomic exposure data.
        log.info("Parse the poller and determine what exposures need to be combined into separate products.\n")
        obs_info_dict, total_list = poller_utils.interpret_obset_input(input_filename)

        # Generate the name for the manifest file which is for the entire visit.  It is fine
        # to use only one of the Total Products to generate the manifest name as the name is not
        # dependent on the detector.
        # Example: instrument_programID_obsetID_manifest.txt (e.g.,wfc3_b46_06_manifest.txt)
        manifest_name = total_list[0].manifest_name
        log.info("\nGenerate the manifest name for this visit.  The manifest will contain the names of all the output products.")

        # The product_list is a list of all the output products which will be put into the manifest file
        product_list = []

        # Update all of the product objects with their associated configuration information.
        for total_item in total_list:
            total_item.configobj_pars = config_utils.HapConfig(total_item,
                                                               use_defaults=use_defaults_configs,
                                                               input_custom_pars_file=input_custom_pars_file,
                                                               output_custom_pars_file=output_custom_pars_file)
            for filter_item in total_item.fdp_list:
                filter_item.configobj_pars = config_utils.HapConfig(filter_item,
                                                                    use_defaults=use_defaults_configs,
                                                                    input_custom_pars_file=input_custom_pars_file,
                                                                    output_custom_pars_file=output_custom_pars_file)
            for expo_item in total_item.edp_list:
                expo_item.configobj_pars = config_utils.HapConfig(expo_item,
                                                                  use_defaults=use_defaults_configs,
                                                                  input_custom_pars_file=input_custom_pars_file,
                                                                  output_custom_pars_file=output_custom_pars_file)

        log.info("The configuration parameters have been read and applied to the drizzle objects.")

        # Run alignimages.py on images on a filter-by-filter basis.
        # Process each filter object which contains a list of exposure objects/products.
        log.info("\nAlign the images on a filter-by-filter basis.")
        for tot_obj in total_list:
            for filt_obj in tot_obj.fdp_list:
                align_table, filt_exposures = filt_obj.align_to_gaia(output=debug)

                # Report results and track the output files
                if align_table:
                    log.info("ALIGN_TABLE: {}".format(align_table.filtered_table))
                    # os.remove("alignimages.log")  # FIX This log needs to be included in total product trailer file
                    for row in align_table.filtered_table:
                        log.info(row['status'])
                        if row['status'] == 0:
                            log.info("Successfully aligned {} to {} astrometric frame\n".format(row['imageName'], row['catalog']))

                        # Alignment did not work for this particular image
                        # FIX - If alignment did not work for an image, it seems this exposure should
                        # be removed from the exposure lists.  TotalProduct and FilterProduct need
                        # methods to do this.
                        else:
                            log.info("Could not align {} to absolute astrometric frame\n".format(row['imageName']))

                    hdrlet_list = align_table.filtered_table['headerletFile'].tolist()
                    product_list += hdrlet_list
                    product_list += filt_exposures

                else:
                    log.warning("Step to align the images has failed. No alignment table has been generated.")

        # Run AstroDrizzle to produce drizzle-combined products
        log.info("\nCreate drizzled imagery products.")
        driz_list = create_drizzle_products(total_list)
        product_list += driz_list

        # Create source catalogs from newly defined products (HLA-204)
        log.info("Create source catalog from newly defined product.\n")
        if "total detection product 00" in obs_info_dict.keys():
            catalog_list = create_catalog_products(total_list, debug=debug, phot_mode=phot_mode)
            product_list += catalog_list
        else:
            log.warning("No total detection product has been produced.  The sourcelist generation step has been skipped.")
        """
        # 8: (OPTIONAL) Determine whether there are any problems with alignment or photometry of product
        log.info("8: (TODO) (OPTIONAL) Determine whether there are any problems with alignment or photometry"
                 "of product")
        # TODO: QUALITY CONTROL SUBROUTINE CALL GOES HERE.
        """

        # Write out manifest file listing all products generated during processing
        log.info("Creating manifest file {}.".format(manifest_name))
        log.info("  The manifest contains the names of products generated during processing.")
        with open(manifest_name, mode='w') as catfile:
            [catfile.write("{}\n".format(name)) for name in product_list]

        # 10: Return exit code for use by calling Condor/OWL workflow code: 0 (zero) for success, 1 for error condition
        return_value = 0
    except Exception:
        return_value = 1
        log.info("\a\a\a")
        exc_type, exc_value, exc_tb = sys.exc_info()
        traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)
    finally:
        log.info('Total processing time: {} sec'.format((datetime.datetime.now() - starting_dt).total_seconds()))
        log.info("Return exit code for use by calling Condor/OWL workflow code: 0 (zero) for success, 1 for error "
                 "condition {}".format(return_value))
        return return_value


# ----------------------------------------------------------------------------------------------------------------------

def run_sourcelist_flagging(filter_product_obj, filter_product_catalogs, debug = True):
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

    debug : Boolean, optional.
        create intermediate diagnostic files? Default value is True.

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
    hla_flag_filter.make_mask_file(drizzled_image)
    if filter_product_obj.configobj_pars.use_defaults:
        ci_lookup_file_path = "default_parameters/any"
    else:
        ci_lookup_file_path = "user_parameters/any"
    output_custom_pars_file = filter_product_obj.configobj_pars.output_custom_pars_file
    for cat_type in filter_product_catalogs.catalogs.keys():
        exptime = filter_product_catalogs.catalogs[cat_type].image.imghdu[0].header['exptime'] #TODO: This works for ACS. Make sure that it also works for WFC3. Look at "TEXPTIME"
        catalog_name = filter_product_catalogs.catalogs[cat_type].sourcelist_filename
        catalog_data = filter_product_catalogs.catalogs[cat_type].source_cat
        drz_root_dir = os.getcwd()
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
                                                                                                         ci_lookup_file_path,
                                                                                                         output_custom_pars_file,
                                                                                                         debug)

    return filter_product_catalogs




