#!/usr/bin/env python

"""This script contains code to support creation of source extractor-like and daophot-like sourcelists.

"""
import pdb
import sys

from stsci.tools import logutil

from drizzlepac import util

__taskname__ = 'sourcelist_generation'

log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout)

# ----------------------------------------------------------------------------------------------------------------------




def create_dao_like_coordlists():
    """Make daofind-like coordinate lists

    Parameters
    ----------

    Returns
    -------
    """

    log.info("DAOFIND-LIKE COORDINATE LIST CREATION OCCURS HERE!")


# ----------------------------------------------------------------------------------------------------------------------


def create_dao_like_sourcelists():
    """Make DAOphot-like sourcelists

    Parameters
    ----------


    Returns
    -------
    """
    log.info("DAOPHOT-LIKE SOURCELIST CREATION OCCURS HERE!")

# ----------------------------------------------------------------------------------------------------------------------


def create_se_like_coordlists():
    """Make source extractor-like coordinate lists

    Parameters
    ----------

    Returns
    -------
    """

    log.info("SOURCE EXTRACTOR-LIKE COORDINATE LIST CREATION OCCURS HERE!")


# ----------------------------------------------------------------------------------------------------------------------


def create_se_like_sourcelists():
    """Make source extractor-like sourcelists

    Parameters
    ----------

    Returns
    -------
    """

    log.info("SOURCE EXTRACTOR-LIKE SOURCELIST CREATION OCCURS HERE!")


# ----------------------------------------------------------------------------------------------------------------------


def create_sourcelists(obs_info_dict, param_dict):
    """Main calling code. Make sourcelists

    Parameters
    ----------
    obs_info_dict : dictionary
        Dictionary containing all information about the images being processed

    param_dict : dictionary
        Dictionary of instrument/detector - specific drizzle, source finding and photometric parameters

    Returns
    -------
    """
    log.info("----------------------------------------------------------------------------------------------------------------------")
    for key1 in list(obs_info_dict.keys()):
        for key2 in list(obs_info_dict[key1].keys()):
            log.info("obs_info_dict[{}][{}]: {}".format(key1,key2,obs_info_dict[key1][key2]))  # TODO: REMOVE THIS SECTION BEFORE ACTUAL USE

    log.info("----------------------------------------------------------------------------------------------------------------------")
    log.info("SOURCELIST CREATION OCCURS HERE!")

    create_se_like_coordlists()
    create_se_like_sourcelists()

    create_dao_like_coordlists()
    create_dao_like_sourcelists()




# ----------------------------------------------------------------------------------------------------------------------


@util.with_logging
def run_create_sourcelists(obs_info_dict, param_dict):
    """ subroutine to run create_sourcelists and produce log file when not run sourcelist_generation is not run from
    hlaprocessing.py.

    Parameters
    ----------
    obs_info_dict : dictionary
        Dictionary containing all information about the images being processed

    param_dict : dictionary
        Dictionary of instrument/detector - specific drizzle, source finding and photometric parameters

    Returns
    -------
    """

    create_sourcelists(obs_info_dict, param_dict)


