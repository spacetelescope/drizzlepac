#!/usr/bin/env python

import os
import pdb
import sys

from drizzlepac.hlautils import astrometric_utils
from drizzlepac.hlautils import diagnostic_utils
from stsci.tools import logutil


__taskname__ = 'svm_catalog_analysis.py'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)
# ======================================================================================================================


if __name__ == "__main__":
    # Testing

