#!/usr/bin/env python

"""This script duplicates the functionality of the HLA pipeline.

"""
import argparse
import collections
import datetime
import drizzlepac
from drizzlepac import alignimages
from drizzlepac import generate_final_product_filenames
from drizzlepac import util
from drizzlepac import wcs_functions
import os
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
    obs_info_dict["single exposure product 00"] = {"info": "10265 01S ACS WFC F606W j92c01b4q", "files":["j92c01b4q_flc.fits"]}
    obs_info_dict["single exposure product 01"] = {"info": "10265 01S ACS WFC F606W j92c01b5q", "files":["j92c01b5q_flc.fits"]}
    obs_info_dict["single exposure product 02"] = {"info": "10265 01S ACS WFC F606W j92c01b7q", "files":["j92c01b7q_flc.fits"]}
    obs_info_dict["single exposure product 03"] = {"info": "10265 01S ACS WFC F606W j92c01b9q", "files":["j92c01b9q_flc.fits"]}
    obs_info_dict["filter product 00"] = {"info": "10265 01S ACS WFC F606W", "files":['j92c01b4q_flc.fits',
                                                                                      'j92c01b5q_flc.fits',
                                                                                      'j92c01b7q_flc.fits',
                                                                                      'j92c01b9q_flc.fits']}
    obs_info_dict["total detection product 00"] = {"info": "10265 01S ACS WFC", "files":['j92c01b4q_flc.fits',
                                                                                               'j92c01b5q_flc.fits',
                                                                                               'j92c01b7q_flc.fits',
                                                                                               'j92c01b9q_flc.fits']}
    return(obs_info_dict)

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
    run_hla_processing(input_filename,result=return_value,**kwargs)
    return(return_value[0])

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
def rename_subproduct_files(obs_info_dict_item):
    """
    renames subproduct images (single-exposure products)

    Parameters
    ----------
    obs_info_dict_item : dictionary
        obs_info_dict singleton that may contain files to be renamed.

    Returns
    -------
    Nothing.
    """
    # Bail out if there are no subproducts to rename.
    if "subproduct #0 filenames" not in obs_info_dict_item.keys():
        log.info("No subproduct image files to rename.")
        return()
    else:
        for key in obs_info_dict_item.keys():
            log.info("Subproduct image files found.")
            if key.startswith("subproduct"):
                dest_imgname = obs_info_dict_item[key]["image"]
                imgname_root = dest_imgname.split("_")[-2]
                src_imgname = "{}_single_sci.fits".format(imgname_root)

                #rename single_sci.fits image
                os.rename(src_imgname,dest_imgname)
                log.info("RENAME {} ~~> {}".format(src_imgname,dest_imgname))

# ----------------------------------------------------------------------------------------------------------------------

def restructure_obs_info_dict(obs_info_dict):
    """
    restructures obs_info_dict so that single exposure product names become part of their parent filter products or
    total detection products.

    Parameters
    ----------
    obs_info_dict : dictionary
        dictionary to be restructured.

    Returns
    -------
    restructured_dict : ordered dictionary
        reordered and restructured dict

    """
    restructured_dict = collections.OrderedDict()
    single_exposure_dict = {}
    priority_list = ['single exposure product',
                     'filter product',
                     'total detection product',
                     'multivisit mosaic product']
    # 1: reorder dictionary from most complicated product to least complicated product
    for obs_category in priority_list:
        for ctr_b10 in range(0,10000):
            ctr_b36=convert_base10_base36(ctr_b10)
            category_to_search = "{} {}".format(obs_category,ctr_b36)
            if category_to_search in obs_info_dict.keys():
                print("{} FOUND.".format(category_to_search))
                restructured_dict[category_to_search] = obs_info_dict[category_to_search]
                if obs_category.startswith("single"):
                    single_exposure_dict[obs_info_dict[category_to_search]['files'][0]]=category_to_search
            else:
                print("{} not found.".format(category_to_search))
                break

    # 2: have the most complicated products 'absorb' the generated product names of associated single exposure files so
    #    they are not generated twice.
    temp_restructured_dict = restructured_dict.copy()

    for category in temp_restructured_dict.keys():
        if not category.startswith('single exposure product'):
            try:
                for subprod_ctr,imgname in zip(range(0,len(temp_restructured_dict[category]['files'])),temp_restructured_dict[category]['files']):
                    single_exp_dict_key=single_exposure_dict[imgname]
                    restructured_dict[category]["subproduct #{} filenames".format(subprod_ctr)]=temp_restructured_dict[single_exp_dict_key]['product filenames']
                    del restructured_dict[single_exp_dict_key]
                    del single_exposure_dict[imgname]
            except:
                continue
    return(restructured_dict)

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
@util.with_logging
def run_hla_processing(input_filename, result=None, debug=True):
    startingDT = datetime.datetime.now()
    log.info("Run start time: {}".format(str(startingDT)))
    try:
        # 1: Interpret input csv file as an astropy table with defined column names (HLA-211)
        log.info("1: (TODO) Interpret input csv file as an astropy table with defined column names")
        # TODO: SUBROUTINE CALL GOES HERE.

        # 2: Apply rules to determine what exposures need to be combined into separate products (HLA-211 or a new ticket if necessary)
        log.info("2: Apply rules to determine what exposures need to be combined into separate products")
        # TODO: SUBROUTINE CALL GOES HERE.
        obs_info_dict = generate_test_data() #TODO: REMOVE once all previous steps are up and running

        # 3: generate an output names for each defined product...
        log.info("3: generate an output names for each defined product")
        for obs_category in obs_info_dict.keys():
            obs_info_dict[obs_category]['product filenames'] = generate_final_product_filenames.run_generator(obs_category, obs_info_dict[obs_category]["info"])
            for key in obs_info_dict[obs_category].keys():
                log.info("{}: {}".format(key, obs_info_dict[obs_category][key]))

        # 4: restructure obs_info_dict so that it's ready for processing.
        log.info("4: restructure obs_info_dict so that it's ready for processing.")
        obs_info_dict = restructure_obs_info_dict(obs_info_dict)

        # 5: run alignimages.py on images on a filter-by-filter basis.
        log.info("5: run alignimages.py on images on a filter-by-filter basis for {}".format(obs_category))
        wcs_input_list=[]
        for obs_category in obs_info_dict.keys():
            if 'subproduct #0 filenames' in obs_info_dict[obs_category].keys():

                run_perform_align(obs_info_dict[obs_category]['files'])
                wcs_input_list += obs_info_dict[obs_category]['files']
            else:
                log.info("{}: Alignimages step skipped.".format(obs_category))

        #6: run meta wcs code to get common WCS for all images.
        log.info("run make_mosaic_wcs to create a common WCS for all images aligned in the previous step.")
        log.info("The following images will be used: ")
        for imgname in wcs_input_list:log.info("{}".format(imgname))
        if wcs_input_list: meta_wcs = wcs_functions.make_mosaic_wcs(wcs_input_list)

        # # 7: Run AstroDrizzle to produce filter-level products.
        # for obs_category in obs_info_dict.keys():
        #     if 'subproduct #0 filenames' in obs_info_dict[obs_category].keys():
        #
        #
        # pdb.set_trace()
        # # 5: For each defined product...
        # for obs_category in obs_info_dict.keys():
        # #   5.1: align images with alignimages.perform_align()
        #     log.info("5.1: align images with alignimages.perform_align()")
        #     if "subproduct #0 filenames" in obs_info_dict[obs_category].keys():
        #         run_perform_align(obs_info_dict[obs_category]['files'])
        #     else:
        #         log.info("{}: Align_images step skipped.".format(obs_category))
        #
        # #   5.2: Run astrodrizzle on inputs which define the new product using parameters defined by HLA along with the
        # #        newly defined output name
        #     log.info("5.2: Run AstroDrizzle")
        #     for inst_det in astrodrizzle_param_dict.keys():
        #         if obs_info_dict[obs_category]['info'].find(inst_det) != -1:
        #             adriz_param_dict=astrodrizzle_param_dict[inst_det]
        #             break
        #     run_astrodrizzle(obs_info_dict[obs_category]['files'],adriz_param_dict,obs_info_dict[obs_category]['product filenames']['image'])
        #     rename_subproduct_files(obs_info_dict[obs_category])
        #
        # #   5.3: Create source catalog from newly defined product (HLA-204)
        #     log.info("5.3: (TODO) Create source catalog from newly defined product")
        # # TODO: SUBROUTINE CALL GOES HERE.
        #
        # #   5.4: (OPTIONAL) Determine whether there are any problems with alignment or photometry of product
        #     log.info("5.4: (TODO) (OPTIONAL) Determine whether there are any problems with alignment or photometry of product")
        # # 6: (OPTIONAL/TBD) Create trailer file for new product to provide information on processing done to generate the new product.

        # 7: Return exit code for use by calling Condor/OWL workflow code: 0 (zero) for success, 1 for error condition
        return_value = 0
    except:
        return_value = 1
        if debug:
            log.info("\a\a\a")
            exc_type, exc_value, exc_tb = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)

    log.info('Total processing time: {} sec'.format((datetime.datetime.now() - startingDT).total_seconds()))
    log.info("7: Return exit code for use by calling Condor/OWL workflow code: 0 (zero) for success, 1 for error condition")
    result.append(return_value)

# ----------------------------------------------------------------------------------------------------------------------

def run_perform_align(filelist):
    """
    executes drizzlepac.alignimages.perform_align(). If run is successful, and a good fit solution is found, the newly
    created headerlets are applied as the primary WCS in the in flc.fits or flt.fits images.

    Parameters
    ----------
    filelist : list
        List of files to be processed by drizzlepac.alignimages.perform_align().

    Returns
    -------
    Nothing.
    """
    try:
        align_table = alignimages.perform_align(filelist,debug=True,runfile='alignimages.log',update_hdr_wcs=True)
        for row in align_table:
            if row['status'] == 0:
                log.info("Successfully aligned {} to {} astrometric frame\n".format(row['imageName'], row['catalog']))
            else:
                log.info("Could not align {} to absolute astrometric frame\n".format(row['imageName']))

    except Exception:
        # Something went wrong with alignment to GAIA, so report this
        log.info("EXCEPTION encountered in alignimages...\n")
        exc_type, exc_value, exc_tb = sys.exc_info()
        traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)
        log.info("   No correction to absolute astrometric frame applied!\n")
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