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
# -----------------------------------------------------------------------------------------------------------
def run_hla_flag_fitler(hff_inputs):
    source_list = hla_flag_filter.run_source_list_flagging(hff_inputs["drizzled_image"],
                                                           hff_inputs["flt_list"],
                                                           hff_inputs["param_dict"],
                                                           hff_inputs["exptime"],
                                                           hff_inputs["plate_scale"],
                                                           hff_inputs["median_sky"],
                                                           hff_inputs["catalog_name"],
                                                           hff_inputs["catalog_data"],
                                                           hff_inputs["cat_type"],
                                                           hff_inputs["drz_root_dir"],
                                                           hff_inputs["hla_flag_msk"],
                                                           hff_inputs["ci_lookup_file_path"],
                                                           hff_inputs["output_custom_pars_file"],
                                                           hff_inputs["log_level"],
                                                           hff_inputs["diagnostic_mode"])
    return source_list
# -----------------------------------------------------------------------------------------------------------

def run_stuff(input_pickle_filename):
    hff_inputs = pickle.load(open(input_pickle_filename, "rb"))
    sourcelist = run_hla_flag_fitler(hff_inputs)
    pdb.set_trace()

# =======================================================================================================================
if __name__ == "__main__":
    input_pickle_filename = sys.argv[1]
    run_stuff(input_pickle_filename)