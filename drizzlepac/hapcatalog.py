#!/usr/bin/env python
import argparse
import datetime

import os
import sys

from stsci.tools import logutil

from drizzlepac import util
from drizzlepac.devutils.confirm_execution import confirm_execution
from drizzlepac.hlautils.catalog_utils import HAPCatalogs
from drizzlepac import hapsequencer
from drizzlepac.hlautils import config_utils
from drizzlepac.hlautils import poller_utils

log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout)


@util.with_logging
def run_catalog_utils(total_list, debug=False, phot_mode='both'):
    product_list = hapsequencer.create_catalog_products(total_list, debug=debug, phot_mode=phot_mode)
    return product_list

def main():
    """Super simple testing interface for the catalog_utils code."""
    parser = argparse.ArgumentParser(description='test interface for sourcelist_generation')
    parser.add_argument('input_file', help="input filename (ends with '.out'")
    parser.add_argument('-d', '--debug', required=False, choices=['True', 'False'], default='False', help='debug mode on? (generate region files?)')
    parser.add_argument('-m', '--phot_mode', required=False, choices=['aperture', 'segment', 'both'], default='both', help="which photometry mode should be run? 'aperture' for aperture only; 'seg' for segment only, and 'both' for both aperture and segment photometry.")
    args = parser.parse_args()
    if args.debug == "True":
        args.debug = True
    else:
        args.debug = False

    log.info("python {} {} -d {} -m {}".format(os.path.realpath(__file__), args.input_file, args.debug, args.phot_mode))

    obs_info_dict, total_list = poller_utils.interpret_obset_input(args.input_file)
    out_pars_file = 'pars.json'
    for total_item in total_list:
        total_item.configobj_pars = config_utils.HapConfig(total_item, output_custom_pars_file=out_pars_file,use_defaults=True)
        for filter_item in total_item.fdp_list:
            filter_item.configobj_pars = config_utils.HapConfig(filter_item, output_custom_pars_file=out_pars_file,use_defaults=True)
        for expo_item in total_item.edp_list:
            expo_item.configobj_pars = config_utils.HapConfig(expo_item, output_custom_pars_file=out_pars_file,use_defaults=True)

    starting_dt = datetime.datetime.now()
    log.info("Run start time: {}".format(str(starting_dt)))

    product_list = run_catalog_utils(total_list, args.debug, args.phot_mode)

    log.info('Total processing time: {} sec\a'.format((datetime.datetime.now() - starting_dt).total_seconds()))

    for item in product_list:
        print(item)
# ======================================================================================================================



if __name__ == '__main__':
    print("Current working directory: "+os.getcwd())
    confirm_execution()

    cmd_list = ['rm -f *.*','cp orig/* .']
    for cmd in cmd_list:
        print(cmd)
        os.system(cmd)
    main()
    print("\a\a\a")
