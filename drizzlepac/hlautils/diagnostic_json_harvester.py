#!/usr/bin/env python

"""This script 'harvests' information stored in the .json files produced by
drizzlepac/hlautils/svm_quality_analysis.py and stores it as a Pandas DataFrame"""

# Standard library imports
import glob
import json

# Related third party imports
import pandas as pd

# Local application imports
import drizzlepac.hlautils.diagnostic_utils as du
from stsci.tools import logutil


__taskname__ = 'diagnostic_json_harvester'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)

def json_harvester(log_level=logutil.logging.INFO):
    """Main calling function

    Parameters
    ----------
    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the .log file.
        Default value is 'INFO'.

    Returns
    -------
    TBD
    """
    log.setLevel(log_level)



# ======================================================================================================================


if __name__ == "__main__":
    # Testing