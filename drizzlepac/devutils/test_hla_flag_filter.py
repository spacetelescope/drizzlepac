#!/usr/bin/env python
"""script to run hla_flag_filter and compare_sourcelists to test hla_flag_filter parameter updates"""

import glob
import logging
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

# ----------------------------------------------------------------------------------------------------------------------

def run_compare_sourcelists(hff_inputs, log_level):
    """locate HLA classic image and sourcelists, convert HLA classic sorucelist X, Y, RA and Dec to HAP ref frame for
    apples-to-apples comparision, and run comparision code

    Parameters
    ----------
    hff_inputs : dict
        dictionary containing all necessary inputs to run hla_flag_filter.run_source_list_flagging().

    log_level : int
        desired log level

    Returns
    -------
    Nothing.
    """
    # get HLA classic path details from environment variables
    hla_classic_basepath = os.getenv('HLA_CLASSIC_BASEPATH')
    hla_build_ver = os.getenv("HLA_BUILD_VER")
    parse_imgname = hff_inputs['drizzled_image'].split("_")
    instrument = parse_imgname[3]
    prop_id = parse_imgname[1]
    obset_id = parse_imgname[2]
    filter = parse_imgname[5]
    hla_product_basename = "_".join(parse_imgname[0:6])
    hap_product_basename = "_".join(parse_imgname[0:7])
    if hla_classic_basepath and hla_build_ver and os.path.exists(hla_classic_basepath):
        hla_cassic_basepath = os.path.join(hla_classic_basepath, instrument, hla_build_ver)
        hla_classic_path = os.path.join(hla_cassic_basepath, prop_id, prop_id + "_" + obset_id)  # Generate path to HLA classic products
    elif os.path.exists(os.path.join(os.getcwd(), "hla_classic")):  # For local testing
        hla_classic_basepath = os.path.join(os.getcwd(), "hla_classic")
        hla_classic_path = hla_classic_basepath
    else:
        return  # bail out if HLA classic path can't be found.


    hap_imgname = hff_inputs['drizzled_image']
    hla_imgname = glob.glob("{}/{}_dr*.fits".format(hla_classic_path,hla_product_basename))[0]

    if not os.path.exists(hap_imgname) or not os.path.exists(hla_imgname):  # bail if one or both of the images can't be found
        sys.exit("Error: one or both of the images can't be found")
    hap_sourcelist_name = hff_inputs['catalog_name']
    if hap_sourcelist_name.endswith("point-cat.ecsv"):
        hla_classic_cat_type = "dao"
        plotfile_prefix = hap_product_basename + "_point"
    else:
        hla_classic_cat_type = "sex"
        plotfile_prefix = hap_product_basename + "_segment"

    if hla_classic_basepath and hla_build_ver and os.path.exists(hla_classic_basepath):
        hla_sourcelist_name = "{}/logs/{}_{}phot.txt".format(hla_classic_path,
                                                             hla_product_basename,
                                                             hla_classic_cat_type)
    else:
        hla_sourcelist_name = "{}/{}_{}phot.txt".format(hla_classic_path,
                                                        hla_product_basename,
                                                        hla_classic_cat_type)
    plotfile_prefix += "_testing"
    if not os.path.exists(hap_sourcelist_name) or not os.path.exists(hla_sourcelist_name):  # Skip catalog type if one or both of the catalogs can't be found
        sys.exit() # TODO: Circle back to this

    # convert HLA Classic RA and Dec values to HAP reference frame so the RA and Dec comparisons are correct
    updated_hla_sourcelist_name = hapsequencer.correct_hla_classic_ra_dec(hla_sourcelist_name, hap_imgname,
                                                                          hla_classic_cat_type, log_level)
    log.info("HAP image:                   {}".format(os.path.basename(hap_imgname)))
    log.info("HLA Classic image:           {}".format(os.path.basename(hla_imgname)))
    log.info("HAP catalog:                 {}".format(os.path.basename(hap_sourcelist_name)))
    log.info("HLA Classic catalog:         {}".format(os.path.basename(updated_hla_sourcelist_name)))

    # once all file exist checks are passed, execute sourcelist comparision

    return_status = compare_sourcelists.comparesourcelists([updated_hla_sourcelist_name, hap_sourcelist_name],
                                                           [hla_imgname, hap_imgname], good_flag_sum=255,
                                                           plotGen="file", plotfile_prefix=plotfile_prefix,
                                                           verbose=True, log_level=log_level, debugMode=True)
# ----------------------------------------------------------------------------------------------------------------------

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


# ----------------------------------------------------------------------------------------------------------------------


def run_stuff(input_pickle_filename):
    log_level = 10
    log.setLevel(log_level)

    hff_inputs = pickle.load(open(input_pickle_filename, "rb"))

    # run_hla_flag_filter(hff_inputs) # TODO: uncomment once run_correct_hla_classic_ra_dec() up and running

    run_compare_sourcelists(hff_inputs, log_level)


# =======================================================================================================================
if __name__ == "__main__":
    input_pickle_filename = sys.argv[1]
    run_stuff(input_pickle_filename)