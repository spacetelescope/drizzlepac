#!/usr/bin/env python
import argparse
import datetime
import os
import pdb
import sys

from stsci.tools import logutil

from drizzlepac import util
from drizzlepac.hlautils.catalog_utils import HAPCatalogs
from drizzlepac.hlautils import config_utils
from drizzlepac.hlautils import poller_utils


log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout)


@util.with_logging
def run_catalog_utils(total_list, debug=False, phot_mode='both'):
    """This subroutine utilizes hlautils/catalog_utils module to produce photometric sourcelists for the specified
    total drizzle product and it's associated child filter products.

    Parameters
    ----------
    total_list : drizzlepac.hlautils.Product.TotalProduct
        total drizzle product that will be processed by catalog_utils. catalog_utils will also create photometric
        sourcelists for the child filter products of this total product.

    debug : bool, optional
        generate ds9 region file counterparts to the photometric sourcelists? Default value is False.

    phot_mode : str, optional
        Which algorithm should be used to generate the sourcelists? 'aperture' for aperture photometry;
        'segment' for segment map photometry; 'both' for both 'segment' and 'aperture'. Default value is 'both'.

    Returns
    -------
    Nothing.
    """


    product_list = []
    for total_product_obj in total_list:
        # determine total product filename
        if os.path.exists(total_product_obj.product_basename+"_drc.fits"):
            total_product_name = total_product_obj.product_basename + "_drc.fits"
        else:
            total_product_name = total_product_obj.product_basename + "_drz.fits"

        # Instantiate filter catalog product object
        total_product_catalogs = HAPCatalogs(total_product_name, types=phot_mode, debug=debug)

        # Identify sources to be measured by filter photometry step
        total_product_catalogs.identify()

        #write out list(s) of identified sources
        total_product_catalogs.write()

        #append total product catalogs to list
        if phot_mode in ['aperture', 'both']:
            product_list.append(total_product_obj.point_cat_filename)
        if phot_mode in ['segment', 'both']:
            product_list.append(total_product_obj.segment_cat_filename)

        # build dictionary of total_product_catalogs.catalogs[*].sources to use for
        # filter photometric catalog generation
        sources_dict = {}
        for cat_type in total_product_catalogs.catalogs.keys():
            sources_dict[cat_type] = {}
            sources_dict[cat_type]['sources'] = total_product_catalogs.catalogs[cat_type].sources
            if cat_type == "segment":
                sources_dict['segment']['kernel'] = total_product_catalogs.catalogs['segment'].kernel

        for filter_product_obj in total_product_obj.fdp_list:
            # determine filter product filename
            if os.path.exists(filter_product_obj.product_basename + "_drc.fits"):
                filter_product_name = filter_product_obj.product_basename + "_drc.fits"
            else:
                filter_product_name = filter_product_obj.product_basename + "_drz.fits"

            # Instantiate filter catalog product object
            filter_product_catalogs = HAPCatalogs(filter_product_name, types=phot_mode,
                                                  debug=debug, tp_sources=sources_dict)
            # Perform photometry
            filter_product_catalogs.measure()

            # Write out photometric catalog(s)
            filter_product_catalogs.write()

            # append filter product catalogs to list
            if phot_mode in ['aperture', 'both']:
                product_list.append(filter_product_obj.point_cat_filename)
            if phot_mode in ['segment', 'both']:
                product_list.append(filter_product_obj.segment_cat_filename)
    return product_list
# ======================================================================================================================


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

    for total_item in total_list:
        total_item.pars = config_utils.HapConfig(total_item, use_defaults=True)

        for filter_item in total_item.fdp_list:
            filter_item.pars = config_utils.HapConfig(filter_item, use_defaults=True)

        for expo_item in total_item.edp_list:
            expo_item.pars = config_utils.HapConfig(expo_item, use_defaults=True)

    starting_dt = datetime.datetime.now()  # TODO: remove prior to final integration
    log.info("Run start time: {}".format(str(starting_dt)))  # TODO: remove prior to final integration

    product_list = run_catalog_utils(total_list, args.debug, args.phot_mode)

    log.info('Total processing time: {} sec\a'.format((datetime.datetime.now() - starting_dt).total_seconds()))  # TODO: remove prior to final integration

    for item in product_list:
        print(item)


# ======================================================================================================================


if __name__ == '__main__':
    main()



