#!/usr/bin/env python

"""This script contains code to support creation of source extractor-like and daophot-like sourcelists.

"""
import os
import pdb
import sys
import traceback

from astropy.io import fits
import numpy
from stsci.tools import logutil

__taskname__ = 'sourcelist_generation'

log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout)

# ----------------------------------------------------------------------------------------------------------------------
# set up instrument/detector-specific params
# Params imported from the following HLA classic parameter files:
# - https://grit.stsci.edu/HLA/hla/tree/master/software/trunk/HLApipeline/HLApipe/param/parameter_acs_hrc.cfg
# - https://grit.stsci.edu/HLA/hla/tree/master/software/trunk/HLApipeline/HLApipe/param/parameter_acs_sbc.cfg
# - https://grit.stsci.edu/HLA/hla/tree/master/software/trunk/HLApipeline/HLApipe/param/parameter_acs_wfc.cfg
# - https://grit.stsci.edu/HLA/hla/tree/master/software/trunk/HLApipeline/HLApipe/param/parameter_wfc3_ir.cfg
# - https://grit.stsci.edu/HLA/hla/tree/master/software/trunk/HLApipeline/HLApipe/param/parameter_wfc3_uvis.cfg

phot_param_dict = {
    "ACS HRC": {
        "dao": {
            "TWEAK_FWHMPSF": 0.073,
            "TWEAK_THRESHOLD": 3.0,
            "aperture_1": 0.03,
            "aperture_2": 0.125,
            "bthresh": 5.0}},
    "ACS SBC": {
        "dao": {
            "TWEAK_FWHMPSF": 0.065,
            "TWEAK_THRESHOLD": 3.0,
            "aperture_1": 0.07,
            "aperture_2": 0.125,
            "bthresh": 5.0}},
    "ACS WFC": {
        "dao":{
            "TWEAK_FWHMPSF": 0.076,
            "TWEAK_THRESHOLD": 3.0,
            "aperture_1": 0.05,  # update from 0.15
            "aperture_2": 0.15,  # update from 0.25
            "bthresh": 5.0}},
    "WFC3 IR": {
        "dao": {
            "TWEAK_FWHMPSF": 0.14,
            "TWEAK_THRESHOLD": 3.0,
            "aperture_1": 0.15,
            "aperture_2": 0.45,
            "bthresh": 5.0}},
    "WFC3 UVIS": {
        "dao": {
            "TWEAK_FWHMPSF": 0.076,
            "TWEAK_THRESHOLD": 3.0,
            "aperture_1": 0.05,
            "aperture_2": 0.15,
            "bthresh": 5.0}}}
# ----------------------------------------------------------------------------------------------------------------------

def create_sourcelists(obs_info_dict):
    """Make sourcelists

    Parameters
    ----------
    obs_info_dict : dictionary
        Dictionary containing all information about the images being processed

    Returns
    -------
    """
    print("----------------------------------------------------------------------------------------------------------------------")
    os.system("clear")
    for key1 in obs_info_dict.keys():
        for key2 in obs_info_dict[key1].keys():
            print(key1,key2,obs_info_dict[key1][key2])   # TODO: REMOVE THIS SECTION BEFORE ACTUAL USE
        print()
    print("----------------------------------------------------------------------------------------------------------------------")
    log.info("SOURCELIST CREATION OCCURS HERE!")

    for tdp_keyname in [oid_key for oid_key in list(obs_info_dict.keys()) if
                        oid_key.startswith('total detection product')]:  # loop over total filtered products
        # 0: Map image filename to correspoinding catalog filename for total detection product and the associated filter products
        totdet_product_cat_dict = {}
        filter_product_cat_dict = {}
        totdet_product_cat_dict[obs_info_dict[tdp_keyname]['product filenames']['image']] = obs_info_dict[tdp_keyname]['product filenames']['source catalog']
        for fp_keyname in obs_info_dict[tdp_keyname]['associated filter products']:
            filter_product_cat_dict[obs_info_dict[fp_keyname]['product filenames']['image']] = obs_info_dict[fp_keyname]['product filenames']['source catalog']

        inst_det = "{} {}".format(obs_info_dict[tdp_keyname]['info'].split()[-2],
                                  obs_info_dict[tdp_keyname]['info'].split()[-1])
        # 1: Generate daophot-like sourcelist(s)
        create_daophot_like_sourcelists(totdet_product_cat_dict,filter_product_cat_dict,inst_det)

        # 2: Generate source extractor-like sourcelist(s)
        create_se_like_sourcelists(obs_info_dict)


# ----------------------------------------------------------------------------------------------------------------------


def create_daophot_like_sourcelists(totdet_product_cat_dict,filter_product_cat_dict,inst_det):
    """Make daophot-like sourcelists

    Parameters
    ----------
    totdet_product_cat_dict : dictionary
        Dictionary mapping image filename to corresponding catalog filename for total detection product

    filter_product_cat_dict : dictionary
        Dictionary mapping image filename to corresponding catalog filename for total detection product

    inst_det : string
        Text string containing space-deliminated instrument name, detector name (upper case) (i.e. WFC3_UVIS)

    Returns
    -------
    """

    log.info("DAOPHOT-LIKE SOURCELIST CREATION OCCURS HERE!")

    totfiltprod_filename_list = filter_product_cat_dict.keys()

    # ### (0) ### Collect applicable parameters
    log.info("### (0) ### Collect applicable parameters")

    fwhm = float(phot_param_dict[inst_det]["dao"]["TWEAK_FWHMPSF"])
    thresh = float(phot_param_dict[inst_det]["dao"]["TWEAK_THRESHOLD"])
    ap_diameter1 = float(phot_param_dict[inst_det]["dao"]["aperture_1"])
    ap_diameter2 = float(phot_param_dict[inst_det]["dao"]["aperture_2"])
    daofind_basic_param = [fwhm, thresh, ap_diameter1, ap_diameter2]

    # ----------------------------------------
    # Calculate mean readnoise value (float):
    # ----------------------------------------
    readnoise_dictionary_drzs = get_readnoise(totfiltprod_filename_list)

    # -----------------------------
    # Get scale arcseconds / pixel
    # -----------------------------
    scale_dict_drzs = stwcs_get_scale(totfiltprod_filename_list)

    for img in readnoise_dictionary_drzs.keys():
        print(img,readnoise_dictionary_drzs[img],scale_dict_drzs[img])


    # ### (1) ###  White-light source list
    # Create source lists (Returns: name of white-light source-list with path (string)):

    # ### (2) ###  Extract sources that fall "close" to 'INDEF' regions.
    # Take out any sources from the white-light source list falling within 'remove_radius' of a flag.

    # ### (3) ### Feed corrected whitelight source lists into daophot with science images

    # ### (4) ### Gather columns and put in nice format (dictated by: "column_keys_phot.cfg")
    # This will convert columns from xy to ra and dec (controlled by: "column_keys_phot.cfg")



# ----------------------------------------------------------------------------------------------------------------------


def create_se_like_sourcelists(obs_info_dict):
    """Make source extractor-like sourcelists

    Parameters
    ----------
    obs_info_dict : dictionary
        Dictionary containing all information about the images being processed

    Returns
    -------
    """

    log.info("SOURCE EXTRACTOR-LIKE SOURCELIST CREATION OCCURS HERE!")

# ----------------------------------------------------------------------------------------------------------------------


#~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-

def extract_name(stringWpath):
    """
    This task will extract just the name of a specific filename that includes the path in the name: 'stringWpath'.

    Tested.

    :param stringWpath: filename with full path
    :type stringWpath: string
    :returns: naked filename stripped of its path
    """
    while "/" == stringWpath[-1]:
        stringWpath = stringWpath[:-1]
    stringname = stringWpath.split("/")[-1]

    return stringname

# ......................................................................................................................

def get_readnoise(listofimages):
    """
    This task will grab the average readnoise for HST
    data from the header.  This task is known to call the
    correct header keys for ACS UVIS data as well as WFC3
    UVIS and IR data, and WFPC2 data.

    Tested.

    :param listofimages: list of images that will be used to get readnoise values
    :type listofimages: list
    :returns: A dictionary of read noise values keyed by image name
    """
    dictionary_output = {}
    for individual_image in listofimages:
        try:
            ref_readnoise = get_mean_readnoise(individual_image)
        except:  # XXX what kind of exception here?
            exc_type, exc_value, exc_tb = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)
            log.info("ALERT: Readnoise could not be gathered from the header for image {}".format(extract_name(individual_image)))
            ref_readnoise = 0.0
        dictionary_output[individual_image] = ref_readnoise
    return dictionary_output

# ......................................................................................................................

def get_mean_readnoise(image):
    """
    This subroutine computes mean readnoise values

    :param image: image filename
    :type image: string
    :return: mean readnoise
    """
    imghdu = fits.open(image)
    header = imghdu[0].header
    filename = image.split('/')[-1]
    if filename.startswith('hst_'):
        inst = filename.split('_')[3]
        detect = filename.split('_')[4]
    else:
        inst = header['instrume'].lower()
        if inst != 'wfpc2': detect = header['DETECTOR']
        if inst == 'wfpc2':
            rn = header['ROOTNAME'].lower()
            rn = rn.split("/")[-1].split("_")[0]
            if rn[len(rn) - 1] != 'x': detect = 'wfpc2'
            if rn[len(rn) - 1] == 'x': detect = 'pc'

    if inst != 'wfpc2':
        readnsea = header['READNSEA']
        readnseb = header['READNSEB']
        readnsec = header['READNSEC']
        readnsed = header['READNSED']
        readnoise = numpy.mean([readnsea, readnseb, readnsec, readnsed])
    else:
        gain_val = header['ATODGAIN']
        # below values taken from WFPC2 instrument handbook, Table 4.2, Page 81
        if gain_val == 15.0:
            readnsea = 7.02  # PC1
            readnseb = 7.84  # WF2
            readnsec = 6.99  # WF3
            readnsed = 8.32  # WF4
        else:
            # gain = 7 (default)
            readnsea = 5.24  # PC1
            readnseb = 5.51  # WF2
            readnsec = 5.22  # WF3
            readnsed = 5.19  # WF4
        if detect == 'pc':
            readnoise = readnsea
        else:
            readnoise = numpy.mean([readnsea, readnseb, readnsec, readnsed])
    return (readnoise)

# ......................................................................................................................

def stwcs_get_scale(listofimages):
    """
    This task will grab the arcsec/pixel scale for HST
    data from the WCS information of the header.

    Note: Assumes science image is extension: "[1]".

    Tested.

    :param listofimages: list of images that will be used to get scale values
    :type listofimages: list
    :returns: A dictionary of scale values keyed by image name
    """
    import stwcs
    from stwcs import wcsutil
    dictionary_output = {}
    for individual_image in listofimages:
        try:
            wcs1 = stwcs.wcsutil.HSTWCS(individual_image + "[1]")
            scale2return = wcs1.pscale
        except:  # XXX what kind of exception here?
            print
            "ALERT: Pixel-scale could not be gathered from the header for image %s." % (extract_name(individual_image))
            scale2return = False
        dictionary_output[individual_image] = scale2return

    return dictionary_output