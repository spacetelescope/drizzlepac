#!/usr/bin/env python
"""script to run hla_flag_filter and compare_sourcelists to test hla_flag_filter parameter updates"""

import pdb
import pickle
import sys

from drizzlepac import hapsequencer
from drizzlepac.devutils.comparison_tools import compare_sourcelists
from drizzlepac.hlautils import hla_flag_filter
from stsci.tools import logutil
from stwcs import wcsutil

def run_stuff(input_pickle_filename):
    hff_inputs = pickle.load(open(input_pickle_filename, "rb"))
    pdb.set_trace()

# =======================================================================================================================
if __name__ == "__main__":
    input_pickle_filename = sys.argv[1]
    run_stuff(input_pickle_filename)