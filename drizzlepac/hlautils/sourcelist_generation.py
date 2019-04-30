#!/usr/bin/env python

"""This script contains code to support creation of source extractor-like and daophot-like sourcelists.

"""
import sys

from stsci.tools import logutil

__taskname__ = 'sourcelist_generation'

log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout)

# ----------------------------------------------------------------------------------------------------------------------

def create_sourcelists(obs_info_dict,param_dict):
    """Make sourcelists

    Parameters
    ----------
    obs_info_dict : dictionary
        Dictionary containing all information about the images being processed

    param_dict : dictionary
        Dictionary containing instrument/detector-specific values for source identification and photometry ops

    Returns
    -------
    """
    log.info("SOURCELIST CREATION OCCURS HERE!")

    # 0: Generate daophot-like sourcelists
    create_daophot_like_sourcelists(obs_info_dict,param_dict)

    # 1: Generate source extractor-like sourcelists
    create_se_like_sourcelists(obs_info_dict,param_dict)


# ----------------------------------------------------------------------------------------------------------------------


def create_daophot_like_sourcelists(obs_info_dict):
    """Make daophot-like sourcelists

    Parameters
    ----------
    obs_info_dict : dictionary
        Dictionary contianing all information about the images being processed

    Returns
    -------
    """

    log.info("DAOPHOT-LIKE SOURCELIST CREATION OCCURS HERE!")


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
        Dictionary contianing all information about the images being processed

    Returns
    -------
    """

    log.info("SOURCE EXTRACTOR-LIKE SOURCELIST CREATION OCCURS HERE!")

# ----------------------------------------------------------------------------------------------------------------------