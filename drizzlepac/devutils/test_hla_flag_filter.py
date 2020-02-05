#!/usr/bin/env python
"""script to run hla_flag_filter and compare_sourcelists to test hla_flag_filter parameter updates"""

import os
import pdb
import pickle
import sys

from astropy.table import Table
from drizzlepac import hapsequencer
from drizzlepac.devutils.comparison_tools import compare_sourcelists
from drizzlepac.hlautils import hla_flag_filter
from stsci.tools import logutil
from stwcs import wcsutil

__taskname__ = 'test_hla_flag_filter'
MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)


# -----------------------------------------------------------------------------------------------------------
def run_correct_hla_classic_ra_dec(hff_inputs):
    """set up inputs and execute hapsequencer.correct_hla_classic_ra_dec()

    Parameters
    ----------
    hff_inputs : dict
        dictionary containing all necessary inputs to run hla_flag_filter.run_source_list_flagging().

    Returns
    -------
    Nothing.
    """
    print("PLACEHOLDER")


# -----------------------------------------------------------------------------------------------------------
def run_hla_flag_filter(hff_inputs):
    """execute hla_flag_filter.run_source_list_flagging() and write out the catalog to file

    Parameters
    ----------
    hff_inputs : dict
        dictionary containing all necessary inputs to run hla_flag_filter.run_source_list_flagging().

    Returns
    -------
    Nothing.
    """
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
    log.info("renamed existing catalog {} -> {}".format(hff_inputs["catalog_name"],
                                                        hff_inputs["catalog_name"].replace(".ecsv","_old.ecsv")))
    os.rename(hff_inputs["catalog_name"],hff_inputs["catalog_name"].replace(".ecsv","_old.ecsv"))
    source_list.write(hff_inputs['catalog_name'], format='ascii.ecsv')
    log.info('Wrote new catalog {}.'.format(hff_inputs["catalog_name"]))


# -----------------------------------------------------------------------------------------------------------



def run_stuff(input_pickle_filename):

    hff_inputs = pickle.load(open(input_pickle_filename, "rb"))

    run_hla_flag_filter(hff_inputs)

    run_correct_hla_classic_ra_dec(hff_inputs)







# =======================================================================================================================
if __name__ == "__main__":
    log_level = 10
    log.setLevel(log_level)
    input_pickle_filename = sys.argv[1]
    run_stuff(input_pickle_filename)