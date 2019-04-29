#!/usr/bin/env python

"""This script contains code to support creation of source extractor-like and daophot-like sourcelists.

"""
import sys

from stsci.tools import logutil

__taskname__ = 'sourcelist_generation'

log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout)

# ----------------------------------------------------------------------------------------------------------------------

def create_sourcelists():
    """Make sourcelists

    Parameters
    ----------

    Returns
    -------
    """

    log.info("SOURCELIST CREATION OCCURS HERE!")

    create_daophot_like_sourcelists()

    create_se_like_sourcelists()
# ----------------------------------------------------------------------------------------------------------------------


def create_daophot_like_sourcelists():
    """Make daophot-like sourcelists

    Parameters
    ----------

    Returns
    -------
    """

    log.info("DAOPHOT-LIKE SOURCELIST CREATION OCCURS HERE!")

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