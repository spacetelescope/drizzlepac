#!/usr/bin/env python

"""This script contains code to support creation of source extractor-like and daophot-like sourcelists.

"""
import pdb
import sys

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
    log.info("SOURCELIST CREATION OCCURS HERE!")

    # 0: Generate daophot-like sourcelists
    create_daophot_like_sourcelists(obs_info_dict)

    # 1: Generate source extractor-like sourcelists
    create_se_like_sourcelists(obs_info_dict)


# ----------------------------------------------------------------------------------------------------------------------


def create_daophot_like_sourcelists(obs_info_dict):
    """Make daophot-like sourcelists

    Parameters
    ----------
    obs_info_dict : dictionary
        Dictionary containing all information about the images being processed

    Returns
    -------
    """

    log.info("DAOPHOT-LIKE SOURCELIST CREATION OCCURS HERE!")

    for item in [foo for foo in obs_info_dict.keys() if foo.startswith('total detection product')]:
        print(item)

    # ### (1) ### Collect applicable parameters

    # ### (2) ###  White-light source list
    # Create source lists (Returns: name of white-light source-list with path (string)):

    # ### (3) ###  Extract sources that fall "close" to 'INDEF' regions.
    # Take out any sources from the white-light source list falling within 'remove_radius' of a flag.

    # ### (4) ### Feed corrected whitelight source lists into daophot with science images

    # ### (5) ### Gather columns and put in nice format (dictated by: "column_keys_phot.cfg")
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