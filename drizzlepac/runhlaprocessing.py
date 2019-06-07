#!/usr/bin/env python

"""This script duplicates the functionality of the HLA pipeline.

"""
import argparse
import collections
import datetime
import os
import pdb
import pickle
import sys
import traceback

import drizzlepac
from drizzlepac import alignimages
from drizzlepac import util
from drizzlepac import wcs_functions
from drizzlepac.hlautils import pipeline_poller_utils
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
            "fwhm": 0.076,
            "thresh": 1.4,
            "bthresh": 5.0,
            "source_box": 7},
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


def generate_test_data(file_name):
    """
    Generates test data for use during development

    Returns
    -------
    obs_info_dict,filelist
    """
    # WFC3 UVIS/IR visit 11665_04, full visit
    if file_name == "ib4604.out":
        obs_info_dict = {'total detection product 00':
                         {'info': '11665 B46 WFC3 UVIS',
                          'files': ['ib4604fmq_flc.fits',
                                    'ib4604fxq_flc.fits',
                                    'ib4604fnq_flc.fits',
                                    'ib4604fuq_flc.fits',
                                    'ib4604fqq_flc.fits',
                                    'ib4604fvq_flc.fits']},
                         'total detection product 01':
                             {'info': '11665 B46 WFC3 IR',
                              'files': ['ib4604g2q_flt.fits',
                                        'ib4604g6q_flt.fits',
                                        'ib4604g3q_flt.fits',
                                        'ib4604g8q_flt.fits']},
                         'filter product 00':
                             {'info': '11665 B46 WFC3 IR F110W',
                              'files': ['ib4604g2q_flt.fits',
                                        'ib4604g6q_flt.fits']},
                         'single exposure product 00':
                             {'info': '11665 B46 WFC3 UVIS F555W IB4604FMQ',
                              'files': ['ib4604fmq_flc.fits']},
                         'single exposure product 01':
                             {'info': '11665 B46 WFC3 UVIS F555W IB4604FXQ',
                              'files': ['ib4604fxq_flc.fits']},
                         'filter product 01':
                             {'info': '11665 B46 WFC3 IR F160W',
                              'files': ['ib4604g3q_flt.fits',
                                        'ib4604g8q_flt.fits']},
                         'single exposure product 02':
                             {'info': '11665 B46 WFC3 UVIS F555W IB4604FNQ',
                              'files': ['ib4604fnq_flc.fits']},
                         'single exposure product 03':
                             {'info': '11665 B46 WFC3 UVIS F555W IB4604FUQ',
                              'files': ['ib4604fuq_flc.fits']},
                         'filter product 02':
                             {'info': '11665 B46 WFC3 UVIS F555W',
                              'files': ['ib4604fqq_flc.fits',
                                        'ib4604fvq_flc.fits']},
                         'single exposure product 04':
                             {'info': '11665 B46 WFC3 UVIS F555W IB4604FQQ',
                              'files': ['ib4604fqq_flc.fits']},
                         'single exposure product 05':
                             {'info': '11665 B46 WFC3 UVIS F555W IB4604FVQ',
                              'files': ['ib4604fvq_flc.fits']},
                         'single exposure product 06':
                             {'info': '11665 B46 WFC3 IR F110W IB4604G2Q',
                              'files': ['ib4604g2q_flt.fits']},
                         'single exposure product 07':
                             {'info': '11665 B46 WFC3 IR F110W IB4604G6Q',
                              'files': ['ib4604g6q_flt.fits']},
                         'single exposure product 08':
                             {'info': '11665 B46 WFC3 IR F160W IB4604G3Q',
                              'files': ['ib4604g3q_flt.fits']},
                         'single exposure product 09':
                             {'info': '11665 B46 WFC3 IR F160W IB4604G8Q',
                              'files': ['ib4604g8q_flt.fits']}}
    if file_name == "ib4604uvis.out":  # WFC3 UVIS/IR visit 11665_04
        obs_info_dict = {'total detection product 00':  # full visit
                         {'info': '11665 B46 WFC3 UVIS',
                          'files': ['ib4604fmq_flc.fits',
                                    'ib4604fxq_flc.fits',
                                    'ib4604fnq_flc.fits',
                                    'ib4604fuq_flc.fits',
                                    'ib4604fqq_flc.fits',
                                    'ib4604fvq_flc.fits']},
                         'single exposure product 00':
                             {'info': '11665 B46 WFC3 UVIS F555W IB4604FMQ',
                              'files': ['ib4604fmq_flc.fits']},
                         'single exposure product 01':
                             {'info': '11665 B46 WFC3 UVIS F555W IB4604FXQ',
                              'files': ['ib4604fxq_flc.fits']},
                         'single exposure product 02':
                             {'info': '11665 B46 WFC3 UVIS F555W IB4604FNQ',
                              'files': ['ib4604fnq_flc.fits']},
                         'single exposure product 03':
                             {'info': '11665 B46 WFC3 UVIS F555W IB4604FUQ',
                              'files': ['ib4604fuq_flc.fits']},
                         'filter product 00':
                             {'info': '11665 B46 WFC3 UVIS F555W',
                              'files': ['ib4604fqq_flc.fits',
                                        'ib4604fvq_flc.fits']},
                         'single exposure product 04':
                             {'info': '11665 B46 WFC3 UVIS F555W IB4604FQQ',
                              'files': ['ib4604fqq_flc.fits']},
                         'single exposure product 05':
                             {'info': '11665 B46 WFC3 UVIS F555W IB4604FVQ',
                              'files': ['ib4604fvq_flc.fits']}}
    if file_name == "ib4604ir.out":
        obs_info_dict = {'total detection product 00':  # Just the WFC3/IR portion of 11665_04
                         {'info': '11665 B46 WFC3 IR',
                          'files': ['ib4604g2q_flt.fits',
                                    'ib4604g6q_flt.fits',
                                    'ib4604g3q_flt.fits',
                                    'ib4604g8q_flt.fits']},
                         'filter product 00':
                             {'info': '11665 B46 WFC3 IR F110W',
                              'files': ['ib4604g2q_flt.fits',
                                        'ib4604g6q_flt.fits']},
                         'filter product 01':
                             {'info': '11665 B46 WFC3 IR F160W',
                              'files': ['ib4604g3q_flt.fits',
                                        'ib4604g8q_flt.fits']},
                         'single exposure product 00':
                             {'info': '11665 B46 WFC3 IR F110W IB4604G2Q',
                              'files': ['ib4604g2q_flt.fits']},
                         'single exposure product 01':
                             {'info': '11665 B46 WFC3 IR F110W IB4604G6Q',
                              'files': ['ib4604g6q_flt.fits']},
                         'single exposure product 02':
                             {'info': '11665 B46 WFC3 IR F160W IB4604G3Q',
                              'files': ['ib4604g3q_flt.fits']},
                         'single exposure product 03':
                             {'info': '11665 B46 WFC3 IR F160W IB4604G8Q',
                              'files': ['ib4604g8q_flt.fits']}}
    if file_name == "j92c01.out":  # obs_info_dict/filelist definition for ACS/WFC visit 10265_01
        obs_info_dict = {"single exposure product 00":
                         {"info": "10265 01S ACS WFC F606W j92c01b4q",
                             "files": ["j92c01b4q_flc.fits"]},
                         "single exposure product 01":
                             {"info": "10265 01S ACS WFC F606W j92c01b5q",
                              "files": ["j92c01b5q_flc.fits"]},
                         "single exposure product 02":
                             {"info": "10265 01S ACS WFC F606W j92c01b7q",
                              "files": ["j92c01b7q_flc.fits"]},
                         "single exposure product 03":
                             {"info": "10265 01S ACS WFC F606W j92c01b9q",
                              "files": ["j92c01b9q_flc.fits"]},
                         "filter product 00":
                             {"info": "10265 01S ACS WFC F606W",
                              "files": ['j92c01b4q_flc.fits',
                                        'j92c01b5q_flc.fits',
                                        'j92c01b7q_flc.fits',
                                        'j92c01b9q_flc.fits']},
                         "total detection product 00":
                             {"info": "10265 01S ACS WFC",
                              "files": ['j92c01b4q_flc.fits',
                                        'j92c01b5q_flc.fits',
                                        'j92c01b7q_flc.fits',
                                        'j92c01b9q_flc.fits']}}
    return(obs_info_dict)

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
                print("\a")
                pdb.set_trace()

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
    # 3: add field "associated 
    for total_driz_product in [x for x in restructured_dict.keys() if x.startswith('total detection product')]:
        restructured_dict[total_driz_product]['associated filter products'] = [y for y in restructured_dict.keys() if
        restructured_dict[y]['info'].startswith(restructured_dict[total_driz_product]['info']) and not
                                                                               y.startswith('total detection product')]
    return(restructured_dict)

# ----------------------------------------------------------------------------------------------------------------------


def run_astrodrizzle(filelist, adriz_param_dict, outfilename, custom_wcs=None):
    """
    Run astrodrizzle on user-specified file(s) with specified parameters.

    Parameters
    ----------
    filelist : list
        List of files to be processed by astrodrizzle.

    adriz_param_dict : dictionary
        Dictionary containing instrument/specific values for astrodrizzle paramters "PIXSCALE", "PIXFRAC":, "KERNEL",
        "OUTNX", "OUTNY", "ROT", and "DRIZ_BITS".

    outfilename : string
        name of the output drizzle-combined image.

    custom_wcs : HSTWCS object
        The composite WCS created by wcs_functions.make_mosaic_wcs()

    RETURNS
    -------
    Nothing.
    """
    log.info("Processing with astrodrizzle version {}".format(drizzlepac.astrodrizzle.__version__))
    # Define parameters which need to be set specifically for
    #    pipeline use of astrodrizzle
    pipeline_pars = {'mdriztab': True,
                     'stepsize': 10,
                     'output': outfilename,
                     'preserve': False,
                     'resetbits': 4096}

    # splice in parameters from instrument/detector-specific astrodrizzle dictionary
    for key in adriz_param_dict.keys():
        if key in ["SCALE", "PIXFRAC", "KERNEL", "OUTNX", "OUTNY", "ROT", "BITS"]:
            pipeline_pars["final_{}".format(key.lower())] = adriz_param_dict[key]
            pipeline_pars["driz_sep_{}".format(key.lower())] = adriz_param_dict[key]
        else:
            pipeline_pars[key] = adriz_param_dict[key]
    # prep custom_wcs values
    if custom_wcs:
        custom_pars = wcs_functions.create_mosaic_pars(custom_wcs)
    # merge custom_pars into pipeline_pars
        log.info("Recombobulating Astrodrizzle input parameters")
        pipeline_keys = pipeline_pars.keys()
        for custom_key in custom_pars.keys():
            if custom_key in pipeline_keys:
                if custom_pars[custom_key] != pipeline_pars[custom_key]:
                    log.info("Updating pipeline_pars value '{}' from {} to {}".format(custom_key,
                                                                                      pipeline_pars[custom_key],
                                                                                      custom_pars[custom_key]))
                    pipeline_pars[custom_key] = custom_pars[custom_key]
            else:
                log.info("Inserting custom_pars value '{}' = {} into pipeline_pars.".format(custom_key,
                                                                                            custom_pars[custom_key]))
                pipeline_pars[custom_key] = custom_pars[custom_key]
        log.info("AstroDrizzle parameter recombobulation successful.")

    # Execute astrodrizzle
    b = drizzlepac.astrodrizzle.AstroDrizzle(input=filelist, runfile="astrodrizzle.log",
                                             configobj='defaults', in_memory=None,
                                             num_cores=None, **pipeline_pars)

# ----------------------------------------------------------------------------------------------------------------------


@util.with_logging
def run_hla_processing(input_filename, result=None, debug=True):
    starting_dt = datetime.datetime.now()
    log.info("Run start time: {}".format(str(starting_dt)))
    try:
        # 1: Interpret input csv file as an astropy table with defined column names (HLA-211)
        log.info("1: (TODO) Interpret input csv file as an astropy table with defined column names")
        # TODO: SUBROUTINE CALL GOES HERE.

        # 2: Apply rules to determine what exposures need to be combined into separate products (HLA-211 or a new
        # ticket if necessary)
        log.info("2: Apply rules to determine what exposures need to be combined into separate products")
        obs_info_dict_old = generate_test_data(input_filename)  # TODO: REMOVE once all previous steps are up and running
        obs_info_dict = pipeline_poller_utils.interpret_obset_input(input_filename)
        # 3: generate an output names for each defined product...
        log.info("3: generate an output names for each defined product")
        for obs_category in obs_info_dict.keys():
            obs_info_dict[obs_category]['product filenames'] = \
                pipeline_poller_utils.run_generator(obs_category, obs_info_dict[obs_category]["info"])
            for key in obs_info_dict[obs_category].keys():
                log.info("{}: {}".format(key, obs_info_dict[obs_category][key]))

        # 4: restructure obs_info_dict so that it's ready for processing.
        log.info("4: restructure obs_info_dict so that it's ready for processing.")
        obs_info_dict = restructure_obs_info_dict(obs_info_dict)

        # 5: run alignimages.py on images on a filter-by-filter basis.
        log.info("5: run alignimages.py on images on a filter-by-filter basis for {}".format(obs_category))
        wcs_input_list = []
        for obs_category in obs_info_dict.keys():
            if 'subproduct #0 filenames' in obs_info_dict[obs_category].keys():

                run_perform_align(obs_info_dict[obs_category]['files'])
                for item in obs_info_dict[obs_category]['files']:
                    wcs_input_list.append(item)
            else:
                log.info("{}: Alignimages step skipped.".format(obs_category))

        # 6: run meta wcs code to get common WCS for all images.
        log.info("6: run make_mosaic_wcs to create a common WCS for all images aligned in the previous step.")
        log.info("The following images will be used: ")
        for imgname in wcs_input_list:
            log.info("{}".format(imgname))
        if wcs_input_list:
            meta_wcs = wcs_functions.make_mosaic_wcs(wcs_input_list)

        # 7: Run AstroDrizzle to produce filter-level products.
        log.info("7: Run AstroDrizzle to produce filter-level products.")
        for obs_category in obs_info_dict.keys():
            if 'subproduct #0 filenames' in obs_info_dict[obs_category].keys():
                adriz_param_dict = {}
                for inst_det in param_dict.keys():
                        if obs_info_dict[obs_category]['info'].find(inst_det) != -1:
                            adriz_param_dict = param_dict[inst_det]['astrodrizzle'].copy()
                            log.info("Using {} AstroDrizzle parameters for {}.".format(inst_det, obs_category))
                            break
                # Turn on astrodrizzle step 7a: Custom WCS for final output
                adriz_param_dict["final_wcs"] = True
                run_astrodrizzle(obs_info_dict[obs_category]['files'],
                                 adriz_param_dict,
                                 obs_info_dict[obs_category]['product filenames']['image'],
                                 custom_wcs=meta_wcs)

                rename_subproduct_files(obs_info_dict[obs_category])
            else:
                log.info("{}: Filter-by-Filter AstroDrizzle step skipped.".format(obs_category))

        # 8: Run AstroDrizzle to produce total detection products
        log.info("8: Run AstroDrizzle to produce total detection products")
        for obs_category in obs_info_dict.keys():
            if obs_category.startswith("total detection product"):
                adriz_param_dict = {}
                for inst_det in param_dict.keys():
                        if obs_info_dict[obs_category]['info'].find(inst_det) != -1:
                            adriz_param_dict = param_dict[inst_det]['astrodrizzle'].copy()
                            log.info("Using {} AstroDrizzle parameters for {}.".format(inst_det, obs_category))
                            break
                # Turn off all astrodrizzle steps EXCEPT steps 7 and 7a.
                adriz_param_dict["static"] = False
                adriz_param_dict["skysub"] = False
                adriz_param_dict["driz_separate"] = False
                adriz_param_dict["driz_sep_wcs"] = False
                adriz_param_dict["median"] = False
                adriz_param_dict["blot"] = False
                adriz_param_dict["driz_combine"] = True
                adriz_param_dict["final_wcs"] = True
                run_astrodrizzle(obs_info_dict[obs_category]['files'],
                                 adriz_param_dict,
                                 obs_info_dict[obs_category]['product filenames']['image'],
                                 custom_wcs=meta_wcs)
            else:
                log.info("{}: Total detection AstroDrizzle step skipped.".format(obs_category))

        # 9: Create source catalogs from newly defined products (HLA-204)
        log.info("9: (WIP) Create source catalog from newly defined product")
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

        # 10: (OPTIONAL) Determine whether there are any problems with alignment or photometry of product
        log.info("10: (TODO) (OPTIONAL) Determine whether there are any problems with alignment or photometry of "
                 "product")
        # TODO: QUALITY CONTROL SUBROUTINE CALL GOES HERE.

        # 11: (OPTIONAL/TBD) Create trailer file for new product to provide information on processing done to generate
        # the new product.

    # 12: Return exit code for use by calling Condor/OWL workflow code: 0 (zero) for success, 1 for error condition
        return_value = 0
    except:
        return_value = 1
        if debug:
            log.info("\a\a\a")
            exc_type, exc_value, exc_tb = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)

    log.info('Total processing time: {} sec'.format((datetime.datetime.now() - starting_dt).total_seconds()))
    log.info("7: Return exit code for use by calling Condor/OWL workflow code: 0 (zero) for success, 1 for error "
             "condition")
    result.append(return_value)

# ----------------------------------------------------------------------------------------------------------------------


def run_perform_align(filelist):
    """
    executes drizzlepac.alignimages.perform_align(). If run is successful, and a good fit solution is found, the newly
    created headerlets are applied as the primary WCS in the in flc.fits or flt.fits images.

    Parameters
    ----------
    filelist : list
        List of files to be processed by drizzlepac.alignimages.perform_align().

    Returns
    -------
    Nothing.
    """
    try:
        align_table = alignimages.perform_align(filelist, debug=True, runfile='alignimages.log', update_hdr_wcs=True)
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