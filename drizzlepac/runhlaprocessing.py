#!/usr/bin/env python

"""This script duplicates the functionality of the HLA pipeline.

"""
import argparse
import collections
import datetime
import glob
import os
import pdb
import pickle
import sys
import traceback
import shutil

from astropy.io import fits
import drizzlepac
from drizzlepac import alignimages
from drizzlepac import astrodrizzle
from drizzlepac import util
from drizzlepac import wcs_functions
from drizzlepac.hlautils import pipeline_poller_utils
from drizzlepac.hlautils import processing_utils as dpu
from drizzlepac.hlautils import sourcelist_generation
from stsci.tools import logutil

__taskname__ = 'runhlaprocessing'
log = logutil.create_logger('runhlaprocessing', level=logutil.logging.INFO, stream=sys.stdout)

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
            "thresh": None,
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


def convert_base10_base36(in_number):
    """
    Convert base-10 numbers to base-36ish, in the same style that HST visits are named

    Parameters
    ----------
    in_number : integer
        base 10 value to convert to base 36.

    Returns
    --------
    out_val : string
        converted base 36 value
    """
    if in_number < 100:
        out_val = "{}{}".format("0"*(2-len(str(in_number))), in_number)
    elif (in_number > 99) and (in_number < 360):
        alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        alphadict = {}
        for item in enumerate(list(alphabet)):
            alphadict[item[0]] = item[1]
        c1 = (in_number - 100)//26
        c2 = (in_number - 100) % 26
        out_val = "{}{}".format(c1, alphadict[c2])
    else:

        chars = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'

        sign = '-' if in_number < 0 else ''
        in_number = abs(in_number)
        out_val = ''

        while in_number > 0:
            in_number, remainder = divmod(in_number, 36)
            out_val = chars[remainder] + out_val

    return(out_val)


# ----------------------------------------------------------------------------------------------------------------------

def perform_processing(input_filename, **kwargs):
    """
    Main calling subroutine.

    Parameters
    ----------
    input_filename : string
        Name of the input csv file containing information about the files to
        be processed

    debug : Boolean
        display all tracebacks, and debug information?

    Updates
    -------
    return_value : list
        a simple status value. '0' for a successful run and '1' for a failed
        run
    """
    return_value = []
    run_hla_processing(input_filename, result=return_value, **kwargs)
    return(return_value[0])

# ----------------------------------------------------------------------------------------------------------------------


def rename_subproduct_files(obs_info_dict_item):
    """
    renames subproduct images (single-exposure products)

    Parameters
    ----------
    obs_info_dict_item : dictionary
        obs_info_dict singleton that may contain files to be renamed.

    Returns
    -------
    Nothing.
    """
    # Bail out if there are no subproducts to rename.
    if "subproduct #0 filenames" not in obs_info_dict_item.keys():
        log.info("No subproduct image files to rename.")
        return()
    else:
        for key in obs_info_dict_item.keys():
            log.info("Subproduct image files found.")
            if key.startswith("subproduct"):
                dest_imgname = obs_info_dict_item[key]["image"]
                imgname_root = dest_imgname.split("_")[-2]
                image_rootname = [i for i in obs_info_dict_item['files'] if i.startswith("{}".format(imgname_root))][0].split("_")[0]
                src_imgname = "{}_single_sci.fits".format(image_rootname)

                # rename single_sci.fits image
                os.rename(src_imgname, dest_imgname)
                log.info("RENAME {} ~~> {}".format(src_imgname, dest_imgname))


# ----------------------------------------------------------------------------------------------------------------------


def restructure_obs_info_dict(obs_info_dict):
    """
    restructures obs_info_dict so that single exposure product names become part of their parent filter products or
    total detection products.

    Parameters
    ----------
    obs_info_dict : dictionary
        dictionary to be restructured.

    Returns
    -------
    restructured_dict : ordered dictionary
        reordered and restructured dict

    """
    restructured_dict = collections.OrderedDict()
    single_exposure_dict = {}
    priority_list = ['single exposure product',
                     'filter product',
                     'total detection product',
                     'multivisit mosaic product']
    # 1: reorder dictionary from most complicated product to least complicated product
    for obs_category in priority_list:
        for ctr_b10 in range(0, 10000):
            ctr_b36 = convert_base10_base36(ctr_b10)
            category_to_search = "{} {}".format(obs_category, ctr_b36)
            if category_to_search in obs_info_dict.keys():
                print("{} FOUND.".format(category_to_search))
                restructured_dict[category_to_search] = obs_info_dict[category_to_search]
                if obs_category.startswith("single"):
                    single_exposure_dict[obs_info_dict[category_to_search]['files'][0]] = category_to_search
            else:
                print("{} not found.".format(category_to_search))
                break

    # 2: have the most complicated products 'absorb' the generated product names of associated single exposure files so
    #    they are not generated twice.
    temp_restructured_dict = restructured_dict.copy()

    for category in temp_restructured_dict.keys():
        if not category.startswith('single exposure product'):
            try:
                for subprod_ctr, imgname in zip(range(0, len(temp_restructured_dict[category]['files'])),
                                                temp_restructured_dict[category]['files']):
                    single_exp_dict_key = single_exposure_dict[imgname]
                    restructured_dict[category]["subproduct #{} filenames".format(subprod_ctr)] = \
                        temp_restructured_dict[single_exp_dict_key]['product filenames']
                    del restructured_dict[single_exp_dict_key]
                    del single_exposure_dict[imgname]
            except:
                continue

    # 3: add field "associated filter products"
    for total_driz_product in [x for x in restructured_dict.keys() if x.startswith('total detection product')]:
        restructured_dict[total_driz_product]['associated filter products'] = [y for y in restructured_dict.keys() if
        set(restructured_dict[y]['files']).issubset(restructured_dict[total_driz_product]['files']) and not
                                                                               y.startswith('total detection product')]
    return(restructured_dict)


# ----------------------------------------------------------------------------------------------------------------------


def run_astrodrizzle(obs_info_dict):
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

    # 0: get rules files for step #6.
    for imgname in glob.glob("*fl?.fits"):
        dpu.get_rules_file(imgname)

    cfgfile_path = os.path.join(os.path.dirname(__file__), "pars")
    for tdp_keyname in [oid_key for oid_key in list(obs_info_dict.keys()) if
                        oid_key.startswith('total detection product')]:  # loop over total filtered products
        # 1: Create temp. total drizzled image used to align all subsequent products
        log.info("~" * 118)
        log.info("CREATE TEMP REFERENCE TOTAL DRIZZLED IMAGE\n")
        ref_total_combined_image = "{}ref_{}".format(obs_info_dict[tdp_keyname]['product filenames']['image'][:-8],
                                                     obs_info_dict[tdp_keyname]['product filenames']['image'][-8:])
        adriz_in_list = obs_info_dict[tdp_keyname]['files']
        log.info("Ref total combined image. {} {}".format(ref_total_combined_image,adriz_in_list,ref_total_combined_image))
        astrodrizzle.AstroDrizzle(input=adriz_in_list,output=ref_total_combined_image,
                                  configobj='{}{}astrodrizzle_total_hap.cfg'.format(cfgfile_path, os.path.sep))

        log.info("Finished creating TEMP REFERENCE TOTAL DRIZZLED IMAGE\n")        
        # Extract shape of ref_total_combined_image for explicit use in AstroDrizzle for all other products.
        rtci = fits.open(ref_total_combined_image)
        total_shape = rtci[('sci',1)].data.shape
        rtci.close()
        
        for fp_keyname in obs_info_dict[tdp_keyname]['associated filter products']:
            # 2: Create drizzle-combined filter image using the temp ref image as astrodrizzle param 'final_refimage'
            log.info("~" * 118)
            log.info("CREATE DRIZZLE-COMBINED FILTER IMAGE\n")
            filter_combined_imagename = obs_info_dict[fp_keyname]['product filenames']['image']
            adriz_in_list = obs_info_dict[fp_keyname]['files']
            trlname = '_'.join(filter_combined_imagename.split('_')[:-1]+['trl.log'])
            print("FILTER PRODUCT trailer file: {}".format(trlname))
            log.info("Filter combined image.... {} {}".format(filter_combined_imagename,adriz_in_list))
            astrodrizzle.AstroDrizzle(input=adriz_in_list, output=filter_combined_imagename,
                                      final_refimage=ref_total_combined_image,
                                      final_outnx=total_shape[1],
                                      final_outny=total_shape[0],
                                      runfile=trlname,
                                      configobj='{}{}astrodrizzle_filter_hap.cfg'.format(cfgfile_path, os.path.sep))
            # Rename Astrodrizzle log file as a trailer file
            shutil.move(trlname, trlname.replace('.log','.txt'))

            # 3: Create individual singly-drizzled images using the temp ref image as astrodrizzle param 'final_refimage'
            for sp_name in [sp_key for sp_key in list(obs_info_dict[fp_keyname].keys()) if
                            sp_key.startswith('subproduct #')]:
                log.info("~" * 118)
                log.info("CREATE SINGLY DRIZZLED IMAGE")
                single_drizzled_filename = obs_info_dict[fp_keyname][sp_name]["image"]
                imgname_root = single_drizzled_filename.split("_")[-2]
                trlname = '_'.join(single_drizzled_filename.split('_')[:-1]+['trl.log'])
                adriz_in_file = [i for i in obs_info_dict[fp_keyname]['files'] if i.startswith(imgname_root)][0]
                log.info("Single drizzled image.... {} {}".format(single_drizzled_filename,adriz_in_file))
                astrodrizzle.AstroDrizzle(input=adriz_in_file, output=single_drizzled_filename,
                                          final_refimage=ref_total_combined_image,
                                          final_outnx=total_shape[1],
                                          final_outny=total_shape[0],
                                          runfile=trlname,
                                          configobj='{}{}astrodrizzle_single_hap.cfg'.format(cfgfile_path, os.path.sep))
                # Rename Astrodrizzle log file as a trailer file
                shutil.move(trlname, trlname.replace('.log','.txt'))

        # 4 Create total image using the temp ref image as astrodrizzle param 'final_refimage'
        log.info("~" * 118)
        log.info("CREATE TOTAL DRIZZLE-COMBINED IMAGE\n")
        total_combined_image = obs_info_dict[tdp_keyname]['product filenames']['image']
        adriz_in_list = obs_info_dict[tdp_keyname]['files']
        trlname = '_'.join(total_combined_image.split('_')[:-1]+['trl.log'])
        log.info("Total combined image..... {} {}".format(total_combined_image,adriz_in_list))
        astrodrizzle.AstroDrizzle(input=adriz_in_list,output=total_combined_image,
                                  final_refimage=ref_total_combined_image,
                                  final_outnx=total_shape[1],
                                  final_outny=total_shape[0],
                                  runfile=trlname,
                                  configobj='{}{}astrodrizzle_total_hap.cfg'.format(cfgfile_path, os.path.sep))
        # Rename Astrodrizzle log file as a trailer file
        shutil.move(trlname, trlname.replace('.log','.txt'))

        # 5: remove reference total temp file
        log.info("Removed temp ref file {}".format(ref_total_combined_image))
        os.remove(ref_total_combined_image)

    # 6: Ensure that all drizzled products is have headers that are to spec.
    drcfiles = sorted(glob.glob('*drc.fits'))
    for d in drcfiles:
        dpu.refine_product_headers(d, obs_info_dict)

    # 7: remove rules files copied to the CWD in step #0
    for rules_filename in glob.glob("*_header_hla.rules"):
        log.info("Removed rules file {}".format(rules_filename))
        os.remove(rules_filename)


# ----------------------------------------------------------------------------------------------------------------------


def run_hla_processing(input_filename, result=None, debug=True):
    starting_dt = datetime.datetime.now()
    log.info("Run start time: {}".format(str(starting_dt)))
    try:
        # 1: Apply rules to determine what exposures need to be combined into separate products (HLA-211 or a new
        # ticket if necessary)
        log.info("1: Apply rules to determine what exposures need to be combined into separate products")
        obs_info_dict = pipeline_poller_utils.interpret_obset_input(input_filename)

        # 2: generate an output names for each defined product...
        log.info("2: generate an output names for each defined product")
        for obs_category in obs_info_dict.keys():
            obs_info_dict[obs_category]['product filenames'] = \
                pipeline_poller_utils.run_generator(obs_category, obs_info_dict[obs_category]["info"])
            for key in obs_info_dict[obs_category].keys():
                log.info("{}: {}".format(key, obs_info_dict[obs_category][key]))

        # 3: restructure obs_info_dict so that it's ready for processing.
        log.info("4: restructure obs_info_dict so that it's ready for processing.")
        # obs_info_dict_old = restructure_obs_info_dict(obs_info_dict_old)
        obs_info_dict = restructure_obs_info_dict(obs_info_dict)

        # 4: run alignimages.py on images on a filter-by-filter basis.
        log.info("4: run alignimages.py on images on a filter-by-filter basis for {}".format(obs_category))
        wcs_input_list = []
        for obs_category in obs_info_dict.keys():
            if 'subproduct #0 filenames' in obs_info_dict[obs_category].keys():
                #create dictionary mapping flc/flt.fits file names to their corresponding HAP-compatible headerlet
                # filenames
                headerlet_filenames={}
                for fitsname in obs_info_dict[obs_category]['files']:
                    for dict_item in obs_info_dict[obs_category].keys():
                        if dict_item.startswith('subproduct #'):
                            if obs_info_dict[obs_category][dict_item]['image'].find(fitsname[:-10]) > 0:
                                headerlet_filenames[fitsname] = \
                                    obs_info_dict[obs_category][dict_item]['image'][:-8]+"hlet.fits"

                run_perform_align(obs_info_dict[obs_category]['files'],headerlet_filenames)

                for item in obs_info_dict[obs_category]['files']:
                    wcs_input_list.append(item)
            else:
                log.info("{}: Alignimages step skipped.".format(obs_category))

        # 5: run meta wcs code to get common WCS for all images.
        log.info("5: run make_mosaic_wcs to create a common WCS for all images aligned in the previous step.")
        log.info("The following images will be used: ")
        for imgname in wcs_input_list:
            log.info("{}".format(imgname))
        if wcs_input_list:
            meta_wcs = wcs_functions.make_mosaic_wcs(wcs_input_list)

        # 6: Run AstroDrizzle to produce drizzle-combined products
        log.info("6: (WIP) Create drizzled imagery products")
        run_astrodrizzle(obs_info_dict)

        # 7: Create source catalogs from newly defined products (HLA-204)
        log.info("7: (WIP) Create source catalog from newly defined product")
        if debug:
            pickle_filename = input_filename.replace(".out",".pickle")
            if os.path.exists(pickle_filename):
                os.remove(pickle_filename)
            pickle_out = open(pickle_filename, "wb")
            pickle.dump([obs_info_dict,param_dict], pickle_out)
            pickle_out.close()
            print("Wrote obs_info_dict to pickle file {}".format(pickle_filename))
        if 'total detection product 00' in obs_info_dict.keys():
            sourcelist_generation.create_sourcelists(obs_info_dict, param_dict)
        else:
            print("Sourcelist generation step skipped.")

        # 8: (OPTIONAL) Determine whether there are any problems with alignment or photometry of product
        log.info("8: (TODO) (OPTIONAL) Determine whether there are any problems with alignment or photometry of "
                 "product")
        # TODO: QUALITY CONTROL SUBROUTINE CALL GOES HERE.

    # 9: Return exit code for use by calling Condor/OWL workflow code: 0 (zero) for success, 1 for error condition
        return_value = 0
    except:
        return_value = 1
        if debug:
            log.info("\a\a\a")
            exc_type, exc_value, exc_tb = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)
    finally:
        log.info('Total processing time: {} sec'.format((datetime.datetime.now() - starting_dt).total_seconds()))
        log.info("9: Return exit code for use by calling Condor/OWL workflow code: 0 (zero) for success, 1 for error "
                 "condition")
        result.append(return_value)

# ----------------------------------------------------------------------------------------------------------------------


def run_perform_align(filelist,headerlet_filenames):
    """
    executes drizzlepac.alignimages.perform_align(). If run is successful, and a good fit solution is found, the newly
    created headerlets are applied as the primary WCS in the in flc.fits or flt.fits images.

    Parameters
    ----------
    filelist : list
        List of files to be processed by drizzlepac.alignimages.perform_align().

    headerlet_filenames : dictionary
        dictionary that maps the flt/flc.fits file name to the corresponding custom headerlet filename.

    Returns
    -------
    Nothing.
    """
    try:
        align_table = alignimages.perform_align(filelist, debug=False, runfile='alignimages.log', update_hdr_wcs=True,headerlet_filenames=headerlet_filenames)
        os.remove("alignimages.log")
        for row in align_table:
            if row['status'] == 0:
                log.info("Successfully aligned {} to {} astrometric frame\n".format(row['imageName'], row['catalog']))
            else:
                log.info("Could not align {} to absolute astrometric frame\n".format(row['imageName']))

    except Exception:
        # Something went wrong with alignment to GAIA, so report this
        log.info("EXCEPTION encountered in alignimages...\n")
        exc_type, exc_value, exc_tb = sys.exc_info()
        traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)
        log.info("   No correction to absolute astrometric frame applied!\n")

# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process images, produce drizzled images and sourcelists')
    parser.add_argument('input_filename', help='Name of the input csv file containing information about the files to '
                        'be processed')
    parser.add_argument('-d', '--debug', required=False, action='store_true', help='If this option is turned on, the '
                        'align_images.perform_align() will attempt to use saved sourcelists stored in a pickle file '
                        'generated during a previous run. Using a saved sorucelist instead of generating new '
                        'sourcelists greatly reduces overall run time. If the pickle file does not exist, the program '
                        'will generate new sourcelists and save them in a pickle file named after the first input '
                        'file.')
    ARGS = parser.parse_args()

    rv = perform_processing(ARGS.input_filename, debug=ARGS.debug)
    print("Return Value: ",rv)
