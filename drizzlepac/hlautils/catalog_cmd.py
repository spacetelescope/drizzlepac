#!/usr/bin/env python
import argparse
import datetime
import os
import sys

from stsci.tools import logutil

from .. import util
from .catalog_utils import HAPPointCatalog, HAPSegmentCatalog

log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout)

@util.with_logging
def run_catalog_utils(args, starting_dt):
    """Super simple testing interface for the above code.

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
    log.info("Run start time: {}".format(str(starting_dt)))
    log.info("python {} {} -f {} -d {} -m {}".format(os.path.realpath(__file__),
                                               args.total_product_name,
                                               " ".join(args.filter_product_list),
                                               args.debug, args.phot_mode))


    if args.phot_mode in ['point', 'both']:
        total_point_product = HAPPointCatalog(args.total_product_name)
        total_point_product.ps_source_cat = total_point_product.identify_point_sources()
        total_point_product.write_to(total_point_product.ps_source_cat, write_region_file=args.debug)

    if args.phot_mode in ['seg', 'both']:
        total_seg_product = HAPSegmentCatalog(args.total_product_name)
        total_seg_product.segmap, \
        total_seg_product.kernel, \
        total_seg_product.bkg_dao_rms = \
            total_seg_product.create_sextractor_like_sourcelists(se_debug=args.debug)

    for filter_img_name in args.filter_product_list:

        if args.phot_mode in ['point', 'both']:
            filter_point_product = HAPPointCatalog(filter_img_name)
            filter_point_product.ps_phot_cat = filter_point_product.perform_point_photometry(total_point_product.ps_source_cat)
            filter_point_product.write_to(filter_point_product.ps_phot_cat, write_region_file=args.debug)

        if args.phot_mode in ['seg', 'both']:
            filter_seg_product = HAPSegmentCatalog(filter_img_name)
            filter_seg_product.measure_source_properties(total_seg_product.segmap, total_seg_product.kernel)

    log.info('Total processing time: {} sec\a'.format((datetime.datetime.now() - starting_dt).total_seconds()))


# ======================================================================================================================



if __name__ == '__main__':
    """Super simple testing interface for the catalog_utils code."""

    starting_dt = datetime.datetime.now()

    parser = argparse.ArgumentParser(description='test interface for sourcelist_generation')
    parser.add_argument('total_product_name', help="total product filename")
    parser.add_argument('-f', '--filter_product_list', nargs='+', required=True,
                        help="Space-separated list of one or more total filter products")
    parser.add_argument('-d', '--debug', required=False, choices=['True', 'False'], default='False', help='debug mode on? (generate region files?)')
    parser.add_argument('-m', '--phot_mode', required=False, choices=['point', 'seg', 'both'], default='both', help="which photometry mode should be run? 'point' for point-soruce only; 'seg' for segment only, and 'both' for both point-source and segment photometry. ")
    args = parser.parse_args()
    if args.debug == "True":
        args.debug = True
    else:
        args.debug = False

    run_catalog_utils(args, starting_dt)
