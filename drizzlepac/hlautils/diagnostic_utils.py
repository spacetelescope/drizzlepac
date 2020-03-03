#!/usr/bin/env python

"""Code related to the creation and modifaction of diagnostic .json files"""

import collections
import json
import pdb
import sys
import numpy as np


from stsci.tools import logutil

__taskname__ = 'diagnostic_utils'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)
# ======================================================================================================================

class HapDiagnostic(object):
    def __init__(self, diag_prod, log_level=logutil.logging.NOTSET):
        """Add insightful docstring content here! # TODO: FLESH OUT!

        Parameters
        ----------
        diag_prod : Not sure yet.
            Definitely something, just not sure what yet. # TODO: FLESH OUT!

        log_level : int, optional
            The desired level of verboseness in the log statements displayed on the screen and written to the .log file.
            Default value is 'NOTSET'.

        Returns
        -------
        Nothing.
        """
        log.setLevel(log_level)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def instantiate_dict(self):
        """Creates a new diagnostic dictionary using the standard format. # TODO: FLESH OUT!

        Parameters
        ----------
        Nothing.

        Returns
        -------
        # TODO: FLESH OUT!
        """

        out_dict = collections.OrderedDict()
        out_dict['header'] = collections.OrderedDict()
        out_dict['data'] = {}

        with open("diag_test.json","w") as json_file:
            json.dump(out_dict, json_file, indent=4)

        pdb.set_trace()

# ======================================================================================================================
if __name__ == "__main__":
    """
    Preliminary list of standardized header values:
    - git revision id
    - creation date/time
    - ippsss
    - proposal
    - visit
    - instrument
    - detector
    - filter
    - catalog type
    - source script (what generated the data?)
    
    Preliminary data section item format:
    Each numpy table column would be converted into a nested dictionary with the fillowing items
    - title (same as the key to this dictionary(?)
    - units
    - dtype
    - data
    - masking information
    """
    HapDiagnostic.instantiate_dict(12)