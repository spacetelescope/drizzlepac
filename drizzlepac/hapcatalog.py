#!/usr/bin/env python
import argparse
import datetime
import os
import sys

from stsci.tools import logutil

from drizzlepac import util
from drizzlepac.hlautils.catalog_utils import HAPCatalogs
from drizzlepac.hlautils import poller_utils

import pdb

log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout)


@util.with_logging
def run_catalog_utils(total_list, debug, phot_mode):
    """Super simple testing interface for the catalog generation code.

    Parameters
    ----------
    args : argparse.Namespace object
        command-line input arguments

    starting_dt : datetime.datetime object
        start date/time of current run.

    Returns
    -------
    Nothing.
    """
    starting_dt = datetime.datetime.now() # TODO: remove prior to final integration
    log.info("Run start time: {}".format(str(starting_dt))) # TODO: remove prior to final integration

    for total_product_obj in total_list:
        if os.path.exists(total_product_obj.product_basename+"_drc.fits"):
            total_product_name = total_product_obj.product_basename + "_drc.fits"
        else:
            total_product_name = total_product_obj.product_basename + "_drz.fits"
        total_product_catalogs = HAPCatalogs(total_product_name, types=phot_mode, debug=debug)
        total_product_catalogs.identify()
        #total_product_catalogs.measure() # <----DON'T NEED TO DO PHOTOMETRY ON IDENTIFIED SOURCES. REMOVE!
        total_product_catalogs.write()

        # build dictionary of total_product_catalogs.catalogs[*].sources to use for filter photometric catalog generation
        sources_dict = {}
        for cat_type in total_product_catalogs.catalogs.keys():
            sources_dict[cat_type] = {}
            sources_dict[cat_type]['sources'] = total_product_catalogs.catalogs[cat_type].sources
            if cat_type == "segment":
                sources_dict['segment']['kernel'] = total_product_catalogs.catalogs['segment'].kernel
        for filter_product_obj in total_product_obj.fdp_list:
            if os.path.exists(filter_product_obj.product_basename + "_drc.fits"):
                filter_product_name = filter_product_obj.product_basename + "_drc.fits"
            else:
                filter_product_name = filter_product_obj.product_basename + "_drz.fits"
            filter_product_catalogs = HAPCatalogs(filter_product_name, types=phot_mode, debug=debug, tp_sources=sources_dict)
            # filter_product_catalogs.identify() # <------FILTER IMAGES DON'T NEED TO RE-IDENTIFY SORUCES AGAIN. REMOVE!
            filter_product_catalogs.measure()
            filter_product_catalogs.write()


    log.info('Total processing time: {} sec\a'.format((datetime.datetime.now() - starting_dt).total_seconds())) # TODO: remove prior to final integration


# ======================================================================================================================


def main():
    """Super simple testing interface for the catalog_utils code."""
    parser = argparse.ArgumentParser(description='test interface for sourcelist_generation')
    parser.add_argument('input_file', help="input filename (ends with '.out'")
    parser.add_argument('-d', '--debug', required=False, choices=['True', 'False'], default='False', help='debug mode on? (generate region files?)')
    parser.add_argument('-m', '--phot_mode', required=False, choices=['point', 'segment', 'both'], default='both', help="which photometry mode should be run? 'point' for point-soruce only; 'seg' for segment only, and 'both' for both point-source and segment photometry. ")
    args = parser.parse_args()
    if args.debug == "True":
        args.debug = True
    else:
        args.debug = False

    log.info("python {} {} -d {} -m {}".format(os.path.realpath(__file__), args.input_file, args.debug, args.phot_mode))

    obs_info_dict, expo_list, filt_list, total_list = poller_utils.interpret_obset_input(args.input_file)

    run_catalog_utils(total_list, args.debug, args.phot_mode)


# ======================================================================================================================


if __name__ == '__main__':
    main()



