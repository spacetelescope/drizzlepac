#!/usr/bin/env python

"""This script duplicates the functionality of the HLA pipeline.

"""
import argparse
from drizzlepac import generate_final_product_filenames
from drizzlepac import runastrodriz
from drizzlepac import util
import pdb
from stsci.tools import logutil
import sys
import traceback

__taskname__ = 'runhlaprocessing'

log = logutil.create_logger('runhlaprocessing', level=logutil.logging.INFO, stream=sys.stdout)

__version__ = 0.1
__version_date__ = '19-Mar-2019'

# ----------------------------------------------------------------------------------------------------------------------

def perform_processing(input_filename, **kwargs):
    """
    Main calling subroutine.

    Parameters
    ----------
    input_filename : string
        Name of the input csv file containing information about the files to
        be processed

    debug : Boolean
        display all tracebacks, and debug information?

    Updates
    -------
    return_value : list
        a simple status value. '0' for a successful run and '1' for a failed
        run
    """
    return_value = []
    run_processing(input_filename,result=return_value,**kwargs)
    return(return_value[0])

# ----------------------------------------------------------------------------------------------------------------------

@util.with_logging
def run_processing(input_filename, result=None, debug=True):
    try:
        # 1: Interpret input csv file as an astropy table with defined column names (HLA-211)

        # 2: Apply rules to determine what exposures need to be combined into separate products (HLA-211 or a new ticket if necessary)
        # PLACEHOLDER TO TEST STEP 3 FUNCTIONALITY UNTIL STEP 2 IS COMPLETED
        obs_info_dict = {}
        obs_info_dict["single exposure product 00"] = "50 A1S WFC3 IR F110W ia1s70jrq" #test proposal_id padding
        obs_info_dict["single exposure product 01"] = "11150 A1S WFC3 UVIS F110W ia1s70jtq"
        obs_info_dict["single exposure product 02"] = "11150 A1S WFC3 IR F110W ia1s70jvq"
        obs_info_dict["single exposure product 03"] = "11150 A1S WFC3 IR F110W ia1s70jwq"
        obs_info_dict["single exposure product 04"] = "11150 A1S WFC3 IR F160W ia1s70jkq"
        obs_info_dict["single exposure product 05"] = "11150 A1S WFC3 IR F160W ia1s70jmq"
        obs_info_dict["single exposure product 06"] = "11150 A1S WFC3 IR F160W ia1s70joq"
        obs_info_dict["single exposure product 07"] = "11150 A1S WFC3 IR F160W ia1s70jpq"
        obs_info_dict["single exposure product 08"] = "10182 A1S ACS HRC PR200LPOL120UV j90za1hyq" #determine maximum generated name length
        obs_info_dict["filter product 00"] = "11150 A1S WFC3 IR F110W"
        obs_info_dict["filter product 01"] = "11150 A1S WFC3 IR F160W"
        obs_info_dict["total detection product 00"] = "11150 A1S WFC3 IR"
        obs_info_dict['multivisit mosaic product 00'] = "1234567 ACS WFC F606W"

        # 3: For each defined product...
        for obs_category in obs_info_dict.keys():
        #   3.1: generate an output name
            product_filename_dict = generate_final_product_filenames.run_generator(obs_category, obs_info_dict[obs_category])
            for key in product_filename_dict.keys():
                log.info("{}: {}".format(key, product_filename_dict[key]))

        #   3.2: Run astrodrizzle on inputs which define the new product using parameters defined by HLA along with the newly defined output name

        #   3.3: Create source catalog from newly defined product (HLA-204)

        #   3.4: (OPTIONAL) Determine whether there are any problems with alignment or photometry of product

        # 4: (OPTIONAL/TBD) Create trailer file for new product to provide information on processing done to generate the new product.

        # 5: Return exit code for use by calling Condor/OWL workflow code: 0 (zero) for success, 1 for error condition
        return_value = 0
    except:
        return_value = 1
        if debug:

            exc_type, exc_value, exc_tb = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)

    result.append(return_value)

# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    PARSER = argparse.ArgumentParser(description='Process images, produce drizzled images and sourcelists')
    PARSER.add_argument('input_filename',help = 'Name of the input csv file containing information about the files to be processed')
    PARSER.add_argument( '-d', '--debug', required=False,choices=['True','False'],default='True',help='display all tracebacks, and debug information? If not otherwise specifed, the default value is "True".')
    ARGS = PARSER.parse_args()

    if ARGS.debug == "True":
        ARGS.debug = True
    else:
        ARGS.debug = False

    rv = perform_processing(ARGS.input_filename,debug=ARGS.debug)
    print(rv)
