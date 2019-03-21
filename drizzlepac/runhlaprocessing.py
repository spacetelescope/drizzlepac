#!/usr/bin/env python

"""This script duplicates the functionality of the HLA pipeline.

"""
import argparse
import drizzlepac
from drizzlepac import generate_final_product_filenames
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
# set up instrument/detector-specific astrodrizzle params
astrodrizzle_param_dict = {}
astrodrizzle_param_dict["ACS HRC"] = {
    "SCALE": 0.025,
    "PIXFRAC": 1.0,
    "KERNEL": "square",
    "OUTNX": None,
    "OUTNY": None,
    "ROT": 0.0,
    "BITS": 256}
astrodrizzle_param_dict["ACS SBC"] = {
    "SCALE": 0.03,
    "PIXFRAC": 1.0,
    "KERNEL": "square",
    "OUTNX": None,
    "OUTNY": None,
    "ROT": 0.0,
    "BITS": 256}
astrodrizzle_param_dict["ACS WFC"] = {
    "SCALE": 0.05,
    "PIXFRAC": 1.0,
    "KERNEL": "square",
    "OUTNX": None,
    "OUTNY": None,
    "ROT": 0.0,
    "BITS": 256}
astrodrizzle_param_dict["WFC3 IR"] = {
    "SCALE": 0.09,
    "PIXFRAC": 1.0,
    "KERNEL": "square",
    "OUTNX": None,
    "OUTNY": None,
    "ROT": 0.0,
    "BITS": 768}
astrodrizzle_param_dict["WFC3 UVIS"] = {
    "SCALE": 0.04,
    "PIXFRAC": 1.0,
    "KERNEL": "square",
    "OUTNX": None,
    "OUTNY": None,
    "ROT": 0.0,
    "BITS": 256}

# ----------------------------------------------------------------------------------------------------------------------

def convert_base10_base36(in_number):
    """
    Convert base-10 numbers to base-36ish, in the same style that HST visits are named

    Parameters
    ----------
    in_number : integer
        base 10 value to convert to base 36.

    Returns
    --------
    out_val : string
        converted base 36 value
    """
    if in_number < 100:
        out_val = "{}{}".format("0"*(2-len(str(in_number))),in_number)
    elif (in_number > 99) and (in_number<360) :
        alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        alphadict = {}
        for item in enumerate(list(alphabet)):
            alphadict[item[0]] = item[1]
        c1 = (in_number - 100)//26
        c2 = (in_number - 100)%26
        out_val = "{}{}".format(c1,alphadict[c2])
    else:

        chars = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'

        sign = '-' if in_number < 0 else ''
        in_number = abs(in_number)
        out_val = ''

        while in_number > 0:
            in_number, remainder = divmod(in_number, 36)
            out_val = chars[remainder] + out_val

    return(out_val)

# ----------------------------------------------------------------------------------------------------------------------

def generate_test_data():
    """
    Generates test data for use during development

    Returns
    -------
    obs_info_dict,filelist
    """
    obs_info_dict = {}
    # obs_info_dict definition #1
    # obs_info_dict["single exposure product 00"] = "50 A1S WFC3 IR F110W ia1s70jrq"  # test proposal_id padding
    # obs_info_dict["single exposure product 01"] = "11150 A1S WFC3 UVIS F110W ia1s70jtq"
    # obs_info_dict["single exposure product 02"] = "11150 A1S WFC3 IR F110W ia1s70jvq"
    # obs_info_dict["single exposure product 03"] = "11150 A1S WFC3 IR F110W ia1s70jwq"
    # obs_info_dict["single exposure product 04"] = "11150 A1S WFC3 IR F160W ia1s70jkq"
    # obs_info_dict["single exposure product 05"] = "11150 A1S WFC3 IR F160W ia1s70jmq"
    # obs_info_dict["single exposure product 06"] = "11150 A1S WFC3 IR F160W ia1s70joq"
    # obs_info_dict["single exposure product 07"] = "11150 A1S WFC3 IR F160W ia1s70jpq"
    # obs_info_dict["single exposure product 08"] = "10182 A1S ACS HRC PR200LPOL120UV j90za1hyq"  # determine maximum generated name length
    # obs_info_dict["filter product 00"] = "11150 A1S WFC3 IR F110W"
    # obs_info_dict["filter product 01"] = "11150 A1S WFC3 IR F160W"
    # obs_info_dict["total detection product 00"] = "11150 A1S WFC3 IR"
    # obs_info_dict['multivisit mosaic product 00'] = "1234567 ACS WFC F606W"

    # obs_info_dict/filelist definition for ACS/WFC visit 10265_01
    # obs_info_dict["single exposure product 00"] = "10265 01S ACS WFC F606W j92c01b4q"
    # obs_info_dict["single exposure product 01"] = "10265 01S ACS WFC F606W j92c01b5q"
    # obs_info_dict["single exposure product 02"] = "10265 01S ACS WFC F606W j92c01b7q"
    # obs_info_dict["single exposure product 03"] = "10265 01S ACS WFC F606W j92c01b9q"
    obs_info_dict["filter product 00"] = "10265 01S ACS WFC F606W"

    file_list = ['j92c01b4q_flc.fits','j92c01b5q_flc.fits','j92c01b7q_flc.fits','j92c01b9q_flc.fits']

    return(obs_info_dict,file_list)

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
        # TODO: SUBROUTINE CALL GOES HERE.

        # 2: Apply rules to determine what exposures need to be combined into separate products (HLA-211 or a new ticket if necessary)
        # TODO: SUBROUTINE CALL GOES HERE.
        obs_info_dict,file_list = generate_test_data() #TODO: REMOVE once all previous steps are up and running

        # 3: For each defined product...
        for obs_category in obs_info_dict.keys():
        #   3.1: generate an output name
            product_filename_dict = generate_final_product_filenames.run_generator(obs_category, obs_info_dict[obs_category])
            for key in product_filename_dict.keys():
                log.info("{}: {}".format(key, product_filename_dict[key]))

        #   3.2: align images with alignimages.perform_align() (I THINK)
        # TODO: SUBROUTINE CALL GOES HERE.

        #   3.3: Run astrodrizzle on inputs which define the new product using parameters defined by HLA along with the
        #        newly defined output name
            for inst_det in astrodrizzle_param_dict.keys():
                if obs_info_dict[obs_category].find(inst_det) != -1:
                    adriz_param_dict=astrodrizzle_param_dict[inst_det]
                    break
            run_astrodrizzle(file_list,adriz_param_dict,product_filename_dict['image'])

        #   3.4: Create source catalog from newly defined product (HLA-204)
        # TODO: SUBROUTINE CALL GOES HERE.

        #   3.5: (OPTIONAL) Determine whether there are any problems with alignment or photometry of product

        # 4: (OPTIONAL/TBD) Create trailer file for new product to provide information on processing done to generate the new product.

        # 5: Return exit code for use by calling Condor/OWL workflow code: 0 (zero) for success, 1 for error condition

        return_value = 0
    except:
        return_value = 1
        if debug:
            log.info("\a\a\a")
            exc_type, exc_value, exc_tb = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)

    result.append(return_value)

# ----------------------------------------------------------------------------------------------------------------------

def run_astrodrizzle(filelist,adriz_param_dict,outfilename):
    """
    Run astrodrizzle on user-specified file(s) with specified parameters.

    Parameters
    ----------
    filelist : list
        List of files to be processed by astrodrizzle.

    adriz_param_dict : dictionary
        Dictionary containing instrument/specific values for astrodrizzle paramters "PIXSCALE", "PIXFRAC":, "KERNEL",
        "OUTNX", "OUTNY", "ROT", and "DRIZ_BITS".

    outfilename: name of the output drizzle-combined image.

    RETURNS
    -------
    Nothing.
    """
    log.info("Processing with astrodrizzle version {}".format(drizzlepac.astrodrizzle.__version__))
    # Define parameters which need to be set specifically for
    #    pipeline use of astrodrizzle
    pipeline_pars = {'mdriztab': True,
                     'stepsize': 10,
                     'output': outfilename,
                     'preserve': False,
                     'resetbits': 4096,
                     'final_wcs': True,
                     }

    # splice in parameters from instrument/detector-specific astrodrizzle dictionary
    for key in adriz_param_dict.keys():
        pipeline_pars["final_{}".format(key.lower())] = adriz_param_dict[key]
        pipeline_pars["driz_sep_{}".format(key.lower())] = adriz_param_dict[key]

    # Execute astrodrizzle
    b = drizzlepac.astrodrizzle.AstroDrizzle(input=filelist, runfile="astrodrizzle.log",
                                             configobj='defaults', in_memory=None,
                                             num_cores=None, **pipeline_pars)


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