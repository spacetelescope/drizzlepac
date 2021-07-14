#!/usr/bin/env python

"""
MAIN DOCSTRING GOES HERE!
"""
# TODO: Add main docstring


import datetime
import glob
import logging
import pdb
import sys


from astropy.table import Table
import numpy as np
import drizzlepac

from drizzlepac.haputils import config_utils
from drizzlepac.haputils import poller_utils
from drizzlepac.haputils import product
from drizzlepac.haputils import processing_utils as proc_utils
from drizzlepac.haputils import svm_quality_analysis as svm_qa

from stsci.tools import logutil
from stwcs import wcsutil

__taskname__ = 'make_custom_mosaic'
MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)

__version__ = 0.1
__version_date__ = '14-July-2021'