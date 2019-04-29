#!/usr/bin/env python

"""This script contains code to support creation of source extractor-like and daophot-like sourcelists.

"""
import sys

from stsci.tools import logutil

__taskname__ = 'sourcelist_generation'

log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout)

def create_sourcelists():
    """Make sourcelists

    Parameters
    ----------

    Returns
    -------
    """

    log.info("SOURCELIST CREATION OCCURS HERE!")