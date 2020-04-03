#!/usr/bin/env python

import glob
import os
import pdb
import sys

from drizzlepac.hlautils import astrometric_utils
from drizzlepac.hlautils import diagnostic_utils
from stsci.tools import logutil


__taskname__ = 'svm_quality_analysis'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)
# ======================================================================================================================


if __name__ == "__main__":
    # Testing
    img_list = glob.glob("*dr?.fits")
    for imgname in img_list:
        print(imgname)

