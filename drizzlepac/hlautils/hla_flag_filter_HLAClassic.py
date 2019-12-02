#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4 ai :
"""Identify and flag sources as either stellar sources, extended sources or anomalous sources

Anomalous sources fall into several categories:

- Saturated sources: the pixel values in the cores of these sources are maxed out at the detector maximum.
- Hot pixels: These sources have intensity profiles that rise much more sharply than a typical gaussian profile for
    that intensity should; Usually indicative of cosmic ray events or hot pixels in the detector.
- Swam detection: These are false detections that occur in areas surrounding bright sources where structures associated
    with an actual bright source are mistaken for actual sources (e.g. diffraction spikes)
- Edge detections: These false detections are from regions where there are a low (or a null) number of contributing
    exposures


**Flag Value Identification Information**

========== ============
Flag Value Flag Meaning
---------- ------------
0          Stellar Source
1          Extended Source (Concentration Index > CI Upper Limit)
4          Saturated Source
16         Concentration Index < CI Lower Limit (i.e. Hot Pixels)
32         Swarm Detection
64         Nexp filtered detections, i.e. edge detections
========== ============

Where the concentration index (CI) = mag1 - mag2


Dependencies
------------
* ci_table.py
"""
import os
import string
import sys
import pdb
import glob
import time
import shutil
import math
import numpy
import datetime
import itertools
import scipy

import astropy.io.fits as pyfits
from astropy.table import Table
new_table = pyfits.BinTableHDU.from_columns
if not hasattr(pyfits,'__version__'):
    pyfits.__version__ = '1.3'

getheader = pyfits.getheader
getdata = pyfits.getdata

# import pywcs
from stwcs import wcsutil
from scipy import spatial

from drizzlepac.hlautils import ci_table
from drizzlepac import util
from stsci.tools import logutil

__taskname__ = 'hla_flag_filter'

log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout)

# config = configparser.ConfigParser()
# Configs=util.toolbox.Configs()
# Rename=util.toolbox.Rename()
# Headers=util.toolbox.Headers()
# Util = util.toolbox.Util()

x_limit = 4096.
y_limit = 2051.


@util.with_logging
def run_source_list_flaging(all_drizzled_filelist, filter_sorted_flt_dict,param_dict, exp_dictionary_scis,
                            dict_newTAB_matched2drz, phot_table_matched2drz, proc_type, drz_root_dir, debug=True):
    """Simple calling subroutine that executes the other flagging subroutines.
    
    Parameters
    ----------
    all_drizzled_filelist : list
        List of drizzled images to process

    filter_sorted_flt_dict : dictionary
        dictionary containing lists of calibrated images sorted (also keyed) by filter name.
    
    param_dict : dictionary
        Dictionary of instrument/detector - specific drizzle, source finding and photometric parameters

    exp_dictionary_scis : dictionary
        dictionary of exposure time values keyed by drizzled image name.

    dict_newTAB_matched2drz : dictionary
        dictionary of source lists keyed by drizzled image name.

    phot_table_matched2drz : dictionary
        dictionary of source lists tables (already read into memory) keyed by drizzled image name.

    proc_type : string
        sourcelist generation type.

    drz_root_dir : string
        Root directory of drizzled images.

    debug : bool
        write intermediate files?
    
    Returns
    -------
    catalog_data : astropy.Table object
        drizzled filter product catalog data with updated flag values
    """
    # -----------------------
    # FLAG FILTER PROCESSING
    # -----------------------
    log.info("************************** * * * HLA_FLAG_FILTER * * * **************************")

    drizzled_image = all_drizzled_filelist[0] # TODO: remove once all code is dictinary-independant
    catalog_name = dict_newTAB_matched2drz[drizzled_image] # TODO: remove once all code is dictinary-independant
    catalog_data = phot_table_matched2drz[drizzled_image] # TODO: remove once all code is dictinary-independant
    for filt_key in filter_sorted_flt_dict.keys(): flt_list = filter_sorted_flt_dict[filt_key] # TODO: remove once all code is dictinary-independant
    exptime = exp_dictionary_scis[drizzled_image] # TODO: remove once all code is dictinary-independant

    # -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
    # Flag sources based on concentration index.
    log.info("ci_filter({} {} {} {} {} {})".format(drizzled_image, catalog_name, "<CATALOG DATA>", proc_type, param_dict,
                                                debug))
    catalog_data = ci_filter(drizzled_image, catalog_name, catalog_data, proc_type, param_dict, debug)

    # -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
    # Flag saturated sources
    log.info("HLASaturationFlags({} {} {} {} {} {} {})".format(drizzled_image, flt_list, catalog_name, "<Catalog Data>",
                                                            proc_type, param_dict, debug))
    catalog_data = HLASaturationFlags(drizzled_image, flt_list, catalog_name, catalog_data, proc_type, param_dict, debug)

    # -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
    # Flag swarm sources
    log.info("HLASwarmFlags({} {} {} {} {} {} {})".format(drizzled_image, catalog_name, "<Catalog Data>", exptime,
                                                       proc_type, param_dict, debug))
    catalog_data = HLASwarmFlags(drizzled_image, catalog_name, catalog_data, exptime, proc_type, param_dict, debug)



    # -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
    # Flag sources from regions where there are a low (or a null) number of contributing exposures
    log.info("HLANexpFlags({} {} {} {} {} {} {})".format(drizzled_image, flt_list, param_dict, catalog_name,
                                                         "<Catalog Data>", drz_root_dir, debug))
    catalog_data = HLANexpFlags(drizzled_image, flt_list, param_dict, catalog_name, catalog_data, drz_root_dir, debug)

    return catalog_data

# ======================================================================================================================

def ci_filter(drizzled_image, catalog_name, catalog_data, proc_type, param_dict, debug):
    """This subroutine flags sources based on concentration index.  Sources below the minimum CI value are
    flagged as hot pixels/CRs (flag=16). Sources above the maximum (for stars) are flagged as extended (flag=1).
    It also flags sources below the detection limit in mag_aper2 (flag=8).

    Parameters
    ----------
    drizzled_image : string
        drizzled filter product image filename

    catalog_name : string
        drizzled filter product catalog filename

    catalog_data : astropy.Table object
        drizzled filter product catalog data

    proc_type : string
        Sourcelist generation type

    param_dict : dictionary
        Dictionary of instrument/detector - specific drizzle, source finding and photometric parameters

    debug : bool
        write intermediate files?

    Returns
    -------
    catalog_data : astropy.Table object
        drizzled filter product catalog data with updated flag values
    """


    # column indices for SE and DAO catalogs
    if proc_type == 'sexphot':
        imag1 = 5
        imag2 = 6
        imerr1 = 7
        imerr2 = 8
    elif proc_type == 'daophot':
        imag1 = 5
        imag2 = 7
        imerr1 = 6
        imerr2 = 8
    else:
        raise ValueError("Unknown proc_type '%s', must be 'sexphot' or 'daophot'" % (proc_type,))

    catalog_name_root = catalog_name.split('.')[0]
    if proc_type == 'sexphot':
        ci_lower_limit = float(param_dict['quality control']['ci filter']['ci_selower_limit'])
        ci_upper_limit = float(param_dict['quality control']['ci filter']['ci_seupper_limit'])
        snr = float(param_dict['quality control']['ci filter']['sourcex_bthresh']) # TODO: Figure out where the best place for bthresh to be

    if proc_type == 'daophot':
        ci_lower_limit = float(param_dict['quality control']['ci filter']['ci_daolower_limit'])
        ci_upper_limit = float(param_dict['quality control']['ci filter']['ci_daoupper_limit'])
        snr = float(param_dict['quality control']['ci filter']['dao_bthresh']) # TODO: Figure out where the best place for bthresh to be

    # replace CI limits with values from table if possible
    cidict = ci_table.get_ci_from_file(drizzled_image, ci_lower=ci_lower_limit, ci_upper=ci_upper_limit)
    ci_lower_limit = cidict['ci_lower_limit']
    ci_upper_limit = cidict['ci_upper_limit']

    log.info(' ')
    log.info('ci limits for {}'.format(drizzled_image))
    log.info('ci_lower_limit = {}'.format(ci_lower_limit))
    log.info('ci_upper_limit = {}'.format(ci_upper_limit))
    log.info(' ')


    if debug:
        failed_index_list=[]
    for i, table_row in enumerate(catalog_data):
        try:
            table_row[-1] = int(table_row[-1])
        except ValueError:
            table_row[-1] = 0

        ci_value = table_row[-2]
        if ci_value:
            ci_value = float(ci_value)
        merr1 = table_row[imerr1]
        if not merr1:
            merr1 = numpy.nan
        else:
            merr1 = float(merr1)
        merr2 = table_row[imerr2]
        if not merr2:
            merr2 = numpy.nan
        else:
            merr2 = float(merr2)
        good_snr = merr2 <= 2.5 / (snr * numpy.log(10))
        ci_err = numpy.sqrt(merr1 ** 2 + merr2 ** 2)

        if not good_snr:
            table_row[-1] |= 8

        if not ci_value or (not numpy.isfinite(ci_err)) or ci_value < ci_lower_limit - ci_err:
            table_row[-1] |= 16


        if not ci_value or ci_value > ci_upper_limit:
            table_row[-1] |= 1

        if not ci_value and debug:
            failed_index_list.append(i)

    if debug:
        # Write out list of ONLY failed rows to to file
        catalog_name_failed = catalog_name_root + '_Failed-CI.txt'
        catalog_data_failed = catalog_data.copy()
        all_indicies = range(0, len(catalog_data))
        rows_to_remove = [z for z in all_indicies if z not in failed_index_list]
        catalog_data_failed.remove_rows(rows_to_remove)
        catalog_data_failed.write(catalog_name_failed,delimiter=",",format='ascii')

        # Write out intermediate catalog with updated flags
        catalog_name = catalog_name_root + 'CIFILT.txt'
        catalog_data.write(catalog_name, delimiter=",", format='ascii')


    return catalog_data

# ======================================================================================================================

def HLASaturationFlags(drizzled_image, flt_list, catalog_name, catalog_data, proc_type, param_dict, debug):
    """Identifies and flags saturated sources.

    Parameters
    ----------
    drizzled_image : string
        drizzled filter product image filename

    flt_list : list
        list of calibrated images that were drizzle-combined to produce image specified by input parameter
        'drizzled_image'

    catalog_name : string
        drizzled filter product catalog filename to process

    catalog_data : astropy.Table object
        drizzled filter product catalog data to process

    proc_type : string
        sourcelist generation type.

    param_dict : dictionary
        Dictionary of instrument/detector - specific drizzle, source finding and photometric parameters

    debug : bool
        write intermediate files?

    Returns
    -------
    phot_table_rows : astropy.Table object
        drizzled filter product catalog data with updated flag values
    """
    image_split = drizzled_image.split('/')[-1]
    channel = drizzled_image.split("_")[-3].upper()  # TODO: May need to be refactored to adjust for new names, and fact that ACS has two filters

    if channel == 'IR': #TODO: Test and IR case just to make sure that IR shouldn't be skipped
        return catalog_data

    phot_table = catalog_name
    phot_table_root = phot_table.split('.')[0]

    #        for flt_image in flt_images:
    #            os.system('cp '+flt_image+' .')

    # -------------------------------------------------------------------
    # STEP THROUGH EACH APPLICABLE FLT IMAGE, DETERMINE THE COORDINATES
    # FOR ALL SATURATION FLAGGED PIXELS, AND TRANSFORM THESE COORDINATES
    # INTO THE DRIZZLED IMAGE REFERENCE FRAME.
    # -------------------------------------------------------------------
    main_drizzled_filelist = drizzled_image
    main_drizzled_filelist_orig = drizzled_image

    drz_filter = drizzled_image.split("_")[5]  # TODO: REFACTOR FOR HAP. this is just a short-term hack to get things working for HLA
    num_flts_in_main_driz = len(flt_list)
    flt_list.sort()

    log.info(' ')
    log.info("Current Working Directory: {}".format(os.getcwd()))
    log.info(' ')
    log.info('LIST OF FLTS IN {}: {}'.format(drizzled_image.split('/')[-1], flt_list))
    log.info(' ')
    log.info('NUMBER OF FLTS IN {}: {}'.format(drizzled_image.split('/')[-1], num_flts_in_main_driz))
    log.info(' ')

    # ----------------------------------------------------
    # EXTRACT DQ DATA FROM FLT IMAGE AND CREATE A LIST
    # OF "ALL" PIXEL COORDINATES WITH A FLAG VALUE OF 256
    # ----------------------------------------------------
    if ((channel.lower() != 'wfpc2') and (channel.lower() != 'pc')):
        image_ext_list = ["[sci,1]", "[sci,2]"]
        dq_sat_bit = 256
    if channel.lower() == 'wfpc2':
        image_ext_list = ["[sci,1]", "[sci,2]", "[sci,3]", "[sci,4]"]
        dq_sat_bit = 8
    if channel.lower() == 'pc':
        image_ext_list = ["[sci,1]"]
        dq_sat_bit = 8

    # build list of arrays
    drz_sat_xy_coords_list = []

    for flt_cnt, flt_image in enumerate(flt_list):
        for ext_cnt, image_ext in enumerate(image_ext_list):
            ext_part = image_ext.split(',')[1].split(']')[0]
            try:
                if ((channel.lower() != 'wfpc2') and (channel.lower() != 'pc')): flt_data = getdata(flt_image, 'DQ',
                                                                                                    int(ext_part))
                if ((channel.lower() == 'wfpc2') or (channel.lower() == 'pc')): flt_data = getdata(
                    flt_image.replace("_c0m", "_c1m"), 'SCI', int(ext_part))
            except KeyError:
                log.info(' ')
                log.info('WARNING: There is only one set of file extensions in {}'.format(flt_image))
                log.info(' ')

                continue

            # ----------------------------------------------------
            # DETERMINE IF ANY OF THE PIXELS LOCATED IN THE GRID
            # HAVE A BIT VALUE OF 256, I.E. FULL WELL SATURATION.
            # ----------------------------------------------------
            # NOTE: NUMPY ARRAYS REPORT Y COORD VALUES FIRST AND
            #       X COORD VALUES SECOND AS FOLLOWS:
            #
            #       --> numpy.shape(flt_data)
            #       (2051, 4096)
            #
            #       WHERE 2051 IS THE NUMBER OF PIXELS IN THE Y
            #       DIRECTION, AND 4096 IS THE NUMBER OF PIXELS
            #       IN THE X DIRECTION.
            # ----------------------------------------------------
            bit_flt_data = dq_sat_bit & flt_data
            complete_sat_coords = numpy.where(bit_flt_data == dq_sat_bit)

            if len(complete_sat_coords[0]) == 0:
                continue

            # -------------------------------------------------
            # RESTRUCTURE THE LIST OF X AND Y COORDINATES FROM
            # THE FLT FILE THAT HAVE BEEN FLAGGED AS SATURATED
            # -------------------------------------------------
            nsat = len(complete_sat_coords[0])
            x_y_array = numpy.empty((nsat, 2), dtype=int)
            x_y_array[:, 0] = complete_sat_coords[1]
            x_y_array[:, 1] = complete_sat_coords[0]

            # ---------------------------------------------------
            # WRITE FLT COORDS TO A FILE FOR DIAGNOSTIC PURPOSES
            # ---------------------------------------------------
            flt_xy_coord_out = flt_image.split('/')[-1].split('.')[0] + '_sci' + str(ext_cnt + 1) + '.txt'
            outfile = open(flt_xy_coord_out, 'w')
            for flt_xy_coord in x_y_array:
                x = flt_xy_coord[0]
                y = flt_xy_coord[1]
                outfile.write(str(x) + '     ' + str(y) + '\n')
            outfile.close()

            # ----------------------------------------------------
            # CONVERT SATURATION FLAGGED X AND Y COORDINATES FROM
            # THE FLT IMAGE INTO RA AND DEC
            # ----------------------------------------------------
            flt_ra_dec_coords = xytord(x_y_array, flt_image, image_ext)

            # -------------------------------------------------
            # CONVERT RA & DEC VALUES FROM FLT REFERENCE FRAME
            # TO THAT OF THE DRIZZLED IMAGE REFERENCE FRAME
            # -------------------------------------------------
            drz_sat_xy_coords_list.append(rdtoxy(flt_ra_dec_coords, drizzled_image, "[sci,1]"))

            log.info(' ')
            log.info('FLT IMAGE = {}'.format(flt_image.split('/')[-1]))
            log.info('IMAGE EXT = {}'.format(image_ext))
            log.info(' ')

    # ----------------------------------------------------------------
    # IF NO SATURATION FLAGS EXIST IN ANY OF THE FLT FILES, THEN SKIP
    # ----------------------------------------------------------------
    if len(drz_sat_xy_coords_list) == 0:
        log.info(' ')
        log.info('*******************************************************************************************')
        log.info('NO SATURATION FLAGGED PIXELS EXIST IN ANY OF THE FLT FILES FOR:')
        log.info('     --> {}'.format(drizzled_image.split('/')[-1]))
        log.info('*******************************************************************************************')
        log.info(' ')

        return catalog_data

    # ------------------------------
    # now concatenate all the arrays
    # ------------------------------
    full_satList = numpy.concatenate(drz_sat_xy_coords_list)

    # --------------------------------------------
    # WRITE RA & DEC FLT CONVERTED X & Y DRIZZLED
    # IMAGE COORDINATES TO A TEXT FILE
    # --------------------------------------------
    drz_coord_file = drizzled_image.split('/')[-1].split('.')[0] + '_ALL_FLT_SAT_FLAG_PIX.txt'
    drz_coord_out = open(drz_coord_file, 'w')
    for coord in full_satList:
        drz_coord_out.write(str(coord[0]) + '     ' + str(coord[1]) + '\n')
    drz_coord_out.close()

    # --------------------------------------------------------
    # READ IN FULL DRIZZLED IMAGE-BASED CATALOG AND SAVE
    # X AND Y COORDINATE VALUES TO LISTS FOR LATER COMPARISON
    # --------------------------------------------------------
    # full_drz_cat = dict_newTAB_matched2drz[drizzled_image]
    # inputfile = open(full_drz_cat, 'r')
    # all_detections = inputfile.readlines()
    # inputfile.close()
    all_detections = catalog_data

    nrows = len(all_detections)
    full_coordList = numpy.empty((nrows, 2), dtype=numpy.float)
    for row_count, detection in enumerate(all_detections):
        full_coordList[row_count, 0] = float(detection[0])
        full_coordList[row_count, 1] = float(detection[1])


    # ----------------------------------------------------
    # CREATE SUB-GROUPS OF SATURATION-FLAGGED COORDINATES
    # ----------------------------------------------------
    proc_time1 = time.ctime()
    log.info(' ')
    log.info('PROC_TIME_1: {}'.format(proc_time1))
    log.info(' ')

    # ----------------------------------
    # Convert aperture radius to pixels
    # ----------------------------------
    ap2 = param_dict['catalog generation']['aperture_2']
    plate_scale = param_dict['catalog generation']['scale'] #TODO: Need to move scale value out of 'catalog generation' > 'dao' to somewhere more sensable
    if proc_type == 'daophot':
        radius = round((ap2/plate_scale) + 0.5) * 2. # TODO: WHY DOES DAOPHOT RADIUS VALUE MULTIPLIED BY 2 BUT SEXPHOT VALUE IS NOT?

    if proc_type == 'sexphot':
        radius = round((ap2 / plate_scale) + 0.5) # TODO: WHY DOES DAOPHOT RADIUS VALUE MULTIPLIED BY 2 BUT SEXPHOT VALUE IS NOT?

    log.info(' ')
    log.info('THE RADIAL DISTANCE BEING USED IS {} PIXELS'.format(str(radius)))
    log.info(' ')

    # do the cross-match using xymatch
    log.info('Matching {} saturated pixels with {} catalog sources'.format(len(full_satList), len(full_coordList)))
    psat, pfull = xymatch(full_satList, full_coordList, radius, multiple=True, verbose=False)
    log.info('Found cross-matches (including duplicates)'.format(len(psat)))
    saturation_flag = numpy.zeros(len(full_coordList), dtype=bool)
    saturation_flag[pfull] = True

    proc_time2 = time.ctime()
    log.info(' ')
    log.info('PROC_TIME_2: {}'.format(proc_time2))
    log.info(' ')

    # ------------------------------------------------------------------
    # REMOVE DUPLICATE DETECTIONS FROM THE LIST, "group", CREATTED FROM
    # MATCHING SATURATION FLAGGED FLT PIXELS TO FINAL SOURCE DETECTIONS
    # ------------------------------------------------------------------

    nsaturated = saturation_flag.sum()
    if nsaturated == 0:
        log.info(' ')
        log.info('**************************************************************************************')
        log.info('NOTE: NO SATURATED SOURCES WERE FOUND FOR: {}'.format(image_split))
        log.info('**************************************************************************************')
        log.info(' ')

        return catalog_data

    else:
        log.info(' ')
        log.info('FLAGGED {} SOURCES'.format(nsaturated))
        log.info(' ')

        if debug:
            sat_coord_file = drizzled_image.split('/')[-1].split('.')[0] + '_INTERMEDIATE.txt'
            sat_coord_out = open(sat_coord_file, 'w')
            for sat_coord in full_coordList[saturation_flag, :]:
                sat_coord_out.write(str(sat_coord[0]) + '     ' + str(sat_coord[1]) + '\n')
            sat_coord_out.close()

        # --------------------------------------------------------------------------
        # WRITE SAT FLAGS TO OUTPUT PHOT TABLE BASED ON flag_src_central_pixel_list
        # --------------------------------------------------------------------------
        phot_table = catalog_name
        phot_table_root = phot_table.split('.')[0]

        phot_table_rows = catalog_data
        for i, table_row in enumerate(phot_table_rows):
            if saturation_flag[i]:
                table_row[-1] = int(table_row[-1]) | 4

        phot_table_rows = HLA_flag4and8_hunter_killer(phot_table_rows)

        if debug:
            phot_table_temp = phot_table_root + '_SATFILT.txt'
            phot_table_rows.write(phot_table_temp, delimiter=",", format='ascii')
        return phot_table_rows

# ======================================================================================================================

def HLASwarmFlags(drizzled_image, catalog_name, catalog_data, exptime, proc_type, param_dict, debug):

    """Identifies and flags swarm sources.

    Parameters
    ----------
    drizzled_image : string
        Name of drizzled image to process
    
    catalog_name : string
        drizzled filter product catalog filename to process

    catalog_data : astropy.Table object
        drizzled filter product catalog data to process

    exptime : float
        exposure of the specified drizzled image

    proc_type : string
        sourcelist generation type.

    param_dict : dictionary
        Dictionary of instrument/detector - specific drizzle, source finding and photometric parameters

    debug : bool
        write intermediate files?
    
    Returns
    -------
    catalog_data : astropy.Table object
        drizzled filter product catalog data with updated flag values
    """

    drz_imgPath_split = drizzled_image.split('/')
    drz_img_split = drz_imgPath_split[-1].split('_')# TODO: May need to be refactored to adjust for new names,
    data_type = drz_img_split[4]
    drz_filter = drz_img_split[5]

    # ============================= SKY COMPUTATION ==============================
    # ----------------------------------------------------------------------------
    # In the pipeline, AstroDrizzle sky computation is based on the statistical
    # distribution of pixel values in an input image. It is performed by
    # iterative sigma-clipping, starting with the full range of pixel intensity
    # values, to calculate the standard deviation. By default, pixel values
    # deviating from the median value (specified by parameter skystat) by over
    # four sigma are rejected, and this operation is repeated for a total of five
    # iterations. The median value of the final distribution is used as the sky
    # value.
    # ----------------------------------------------------------------------------
    # ============================= SKY COMPUTATION ==============================

    median_sky = get_median_sky(drizzled_image) #TODO: maybe get value from filter product object
   # single_rms = rms_dict[drizzled_image]
   # rms_array = pyfits.getdata(single_rms,0)
   # rms_subarray = rms_array[rms_array > 0.0]
   # median_sky = Util.binmode(rms_subarray[rms_subarray < 5000.])[0]

    log.info(' ')
    log.info('MEDIAN SKY VALUE = {}'.format(median_sky))
    log.info(' ')

    # ==========================================
    # ------------------------------------------
    # ------------------------------------------
    # CREATE LIST OF POTENTIAL SWARM DETECTIONS
    # ------------------------------------------
    # ------------------------------------------
    # ==========================================
    phot_table_root = catalog_name.split('/')[-1].split('.')[0]

    image_split = drizzled_image.split('/')[-1]
    channel = drizzled_image.split("_")[-3].lower() # TODO: May need to be refactored to adjust for new names, and fact that ACS has two filters

    ap2 = param_dict['catalog generation']['aperture_2']
    if proc_type not in ('sexphot', 'daophot'):
        raise ValueError("Unknown catalog type '%s'" % proc_type)

    # ----------------------------------
    # Convert aperture radius to pixels
    # ----------------------------------
    radius = ap2 / float(param_dict['catalog generation']['scale']) # TODO: this value should be probably be somewhere else
    log.info(' ')
    log.info('Aperture Size = {}'.format(ap2))
    log.info('Pixel Scale = {} arcsec per pixel'.format(float(param_dict['catalog generation']['scale']))) # TODO: this value should be probably be somewhere else
    log.info(' ')
    area = math.pi * radius**2

    nrows = len(catalog_data)

    complete_src_list = numpy.empty((nrows,6), dtype=numpy.float)

    for row_num,row in enumerate(catalog_data[0:]):
        x_val = float(row[0])
        y_val = float(row[1])

        if proc_type == 'sexphot':
            # mag = row_split[6]
            flux = row[10]
            sky = row[13]
        elif proc_type == 'daophot':
            # mag = row_split[7]
            flux = row[11]
            sky = row[9]

        if not flux:
            flux = 0.0

        if not sky:
            sky = 0.0

        electronpp = flux / area * exptime
        eppsky = electronpp / median_sky
        complete_src_list[row_num,:] = [x_val,y_val,flux,electronpp,sky,eppsky]

    if len(complete_src_list) == 0:
        return catalog_data

    # view into the complete_src_list array for convenience
    swarm_epp_listA = complete_src_list[:,3]

    # swarm flag array
    # this will get set as candidates to flag are found
    swarm_flag = numpy.zeros(nrows, dtype=bool)

    # ------------------------------------------------------------
    # WRITE SUBSET SOURCE LIST TO AN OUTPUT FILE FOR VERIFICATION
    # ------------------------------------------------------------
    final_complete_source_file = open(phot_table_root+'_SWFILT_COMPLETE_SOURCE_FILE.txt','w')
    final_complete_source_file.write("# ------------------------------------------------------------------------------------------------\n")
    final_complete_source_file.write("# X-Center   Y-Center     Flux        ElectronPP          Sky         EPPSKY_Ratio \n")
    final_complete_source_file.write("# ------------------------------------------------------------------------------------------------\n")
    for i, complete_src_value in enumerate(complete_src_list):
        final_complete_source_file.write(str(complete_src_value[0])+'     '+
                                         str(complete_src_value[1])+'     '+
                                         str(complete_src_value[2])+'     '+
                                         str(complete_src_value[3])+'     '+
                                         str(complete_src_value[4])+'     '+
                                         str(complete_src_value[5])+'\n')

    final_complete_source_file.close()

    # ======================================================================
    # ----------------------------------------------------------------------
    # ----------------------------------------------------------------------
    # Introduce 2 thresholds:
    # -----------------------
    # A minimum electronpp, and a minimum electronpp/sky.
    # The thresholds should have different values for IR and UVIS.
    #
    # For IR, sources that have electronpp > 100k, OR
    # ((electronpp > 100*sky) AND (electronpp > 10k)), should be considered.
    #
    # For UVIS, I would set the thresholds at (electronpp > 100k) OR
    # ((electronpp > 1000*sky) AND (electronpp > 10k)).
    # ----------------------------------------------------------------------
    # ----------------------------------------------------------------------
    # ======================================================================
    upper_epp_limit = float(param_dict["quality control"]["swarm filter"]["upper_epp_limit"])
    lower_epp_limit = float(param_dict["quality control"]["swarm filter"]["lower_epp_limit"])
    eppsky_limit_cfg = float(param_dict["quality control"]["swarm filter"]["eppsky_limit"])
    selfradius = float(param_dict["quality control"]["swarm filter"]["selfradius"]) # TODO: optimize selfradius values for ACS/HRC, ACS/SBC in quality control param files

    eppsky_limit = eppsky_limit_cfg * median_sky

    # ----------------------------------------------------------
    # UVIS --> EPP > 100000. OR (EPP > 1000*sky AND EPP > 10000)
    # IR   --> EPP > 100000. OR (EPP > 100*sky AND EPP > 10000)
    # ----------------------------------------------------------

    initial_central_pixel_positions = numpy.where(numpy.logical_or(swarm_epp_listA > upper_epp_limit,
                                                  numpy.logical_and(swarm_epp_listA > eppsky_limit,
                                                                    swarm_epp_listA > lower_epp_limit)))[0]
    initial_central_pixel_list = complete_src_list[initial_central_pixel_positions,:]
    if len(initial_central_pixel_positions) == 0:
        # no bright objects
        # copy empty lists so output file is created anyway
        final_central_pixel_positions = initial_central_pixel_positions
        final_flag_src_central_pixel_list = initial_central_pixel_list
    else:
        # ---------------------------------------------------------
        # Remove duplicate central pixel position swarm candidates
        # Keep only objects that are the brightest within 20 pixels
        # ---------------------------------------------------------

        # -------------------------------------------
        # match initial central pixel list to itself
        # -------------------------------------------
        if data_type.upper() == 'IR':

            # -------------------------------------------------------
            # Define EPP cut values for filtering multiple detections
            # from a given central positions for a swarm candidate
            # -------------------------------------------------------
            cuts = param_dict["quality control"]["swarm filter"]["cuts_list"]
            cuts = list(map(float, cuts))

            selfradii = param_dict["quality control"]["swarm filter"]["selfradii_list"]
            selfradii = list(map(float, selfradii))

            p1 = []
            p2 = []
            for cut_cnt,cut in enumerate(cuts):

                # --------------------------------------------------------------------
                # Extract indices of detections that are within the set EPP cut range
                # --------------------------------------------------------------------
                if cut_cnt == 0:
                    cut_value_positions = numpy.where(initial_central_pixel_list[:,3:4] > cut)[0]
                else:
                    cut_value_positions = numpy.where(numpy.logical_and(initial_central_pixel_list[:,3:4] >= cut,
                                                                        initial_central_pixel_list[:,3:4] <= cuts[cut_cnt-1]))[0]

                # -----------------------------------------------
                # If no detections exist for the specified EPP
                # cut range, then continue to the next cut range
                # -----------------------------------------------
                if len(cut_value_positions) == 0:
                    continue

                # -----------------------------------------------------------------------
                # Determine all matches for detections in "cut_value_positions"
                # within the radius value identified for the cut range being implemented
                # -----------------------------------------------------------------------
                p1_sub, p2_sub = xymatch(initial_central_pixel_list[cut_value_positions,:][:,0:2],
                                         initial_central_pixel_list[:,0:2], selfradii[cut_cnt],
                                         multiple=True, stack=False, verbose=False)

                # ------------------------------------------
                # For each cut range, add the corresponding
                # matches to each detection to a final list
                # ------------------------------------------
                for p1_arr in p1_sub:
                    p1.append(p1_arr)

                for p2_arr in p2_sub:
                    p2.append(p2_arr)

                # Not sure if this is still needed???
                # ------------------------------------
                if cut_cnt == len(cuts) - 1:
                    if len(p1) == 0 and len(p2) == 0:
                        p1, p2 = xymatch(initial_central_pixel_list[:,0:2], initial_central_pixel_list[:,0:2],
                                         selfradius, multiple=True, stack=False, verbose=False)

            # ---------------------------------------------------------------------
            # each object is guaranteed to have at least one match (itself)
            # get brightest of each group of matches by building a list of indices
            # ---------------------------------------------------------------------
            exclude_index = None
            for i1, i2 in zip(p1,p2):
                flux2 = initial_central_pixel_list[i2,2]

                # -------------------------------------------------------------
                # Verify that there is more than one detection in a given group
                # otherwise no detection is added to exclude index because
                # there is only one detection for the source being evaluated
                # -------------------------------------------------------------
                if len(i2[numpy.where(flux2 < numpy.max(flux2))]) > 0:

                    # ----------------------------------------------------------
                    # Add all detections in grouping with a flux value less than
                    # that of the maximum flux value to an array to be excluded
                    # ----------------------------------------------------------
                    if exclude_index is None:
                        exclude_index = i2[numpy.where(flux2 < numpy.max(flux2))]
                    else:
                        exclude_index = numpy.concatenate((exclude_index,i2[numpy.where(flux2 < numpy.max(flux2))]),axis=0)

                    exclude_index = exclude_index.astype(numpy.int32)

            # -----------------------------------------------------------
            # exclude_index can have multiple copies of the same index
            # use exclude_bool array to get a list of the unique indices
            # -----------------------------------------------------------
            exclude_bool = numpy.ones(len(initial_central_pixel_list),dtype=bool)
            if not (exclude_index is None):
                exclude_bool[exclude_index] = False
            out_values = numpy.where(exclude_bool)[0]

            # -------------------------------------------------------------------------------
            # Create final source list based on where the excluded detection indices are not
            # -------------------------------------------------------------------------------
            final_central_pixel_positions = initial_central_pixel_positions[out_values]
            final_flag_src_central_pixel_list = initial_central_pixel_list[out_values,:]

        else:

            p1, p2 = xymatch(initial_central_pixel_list[:,0:2], initial_central_pixel_list[:,0:2], selfradius,
                             multiple=True, stack=False, verbose=False)

            # ---------------------------------------------------------------------
            # each object is guaranteed to have at least one match (itself)
            # get brightest of each group of matches by building a list of indices
            # ---------------------------------------------------------------------
            keep_index = numpy.arange(len(initial_central_pixel_list),dtype=int)
            for i1, i2 in zip(p1,p2):
                flux2 = initial_central_pixel_list[i2,2]
                keep_index[i1] = i2[flux2.argmax()]

            # --------------------------------------------------------
            # keep_index can have multiple copies of the same index
            # use keep_bool array to get a list of the unique indices
            # --------------------------------------------------------
            keep_bool = numpy.zeros(len(initial_central_pixel_list),dtype=bool)
            keep_bool[keep_index] = True
            in_values = numpy.where(keep_bool)[0]
            final_central_pixel_positions = initial_central_pixel_positions[in_values]
            final_flag_src_central_pixel_list = initial_central_pixel_list[in_values,:]


    # ---------------------------------------------------
    # WRITE CENTRAL PIXEL POSITIONS FOR SWARMS TO A FILE
    # ---------------------------------------------------
    if debug:
        cetrl_pix_pos_file = phot_table_root+'_SWFILT_CENTRAL-PIX-POS.txt'
        drz_coord_out = open(cetrl_pix_pos_file,'w')
        for i in range(len(final_flag_src_central_pixel_list)):
            drz_coord_out.write(str(final_flag_src_central_pixel_list[i,0])+'     '+
                                str(final_flag_src_central_pixel_list[i,1])+'     '+
                                str(final_flag_src_central_pixel_list[i,2])+'     '+
                                str(final_flag_src_central_pixel_list[i,3])+'     '+
                                str(final_flag_src_central_pixel_list[i,4])+'     '+
                                str(final_flag_src_central_pixel_list[i,5])+'\n')
        drz_coord_out.close()

    # ==========================================================================
    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # EXTRACT THE CENTRAL PIXEL POSITIONS IN final_flag_src_central_pixel_list,
    # FROM swarm_xListB AND swarm_yListB
    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # ==========================================================================

    swarm_thresh = float(param_dict["quality control"]["swarm filter"]["swarm_thresh"])
    clip_radius_list = param_dict["quality control"]["swarm filter"]["clip_radius_list"]
    clip_radius_list = list(map(float, clip_radius_list))
    scale_factor_list = param_dict["quality control"]["swarm filter"]["scale_factor_list"]
    scale_factor_list = list(map(float, scale_factor_list))
    log.info('SWARM FILTER CLIP_RADIUS_LIST: {}'.format(clip_radius_list))
    log.info('SWARM FILTER SCALE_FACTOR_LIST: {}'.format(scale_factor_list))

    # get list of objects not in the central pixel list
    keep = numpy.ones(nrows, dtype=bool)
    keep[final_central_pixel_positions] = False
    notcentral_index = numpy.where(keep)[0]
    swarm_listB = complete_src_list[notcentral_index,:]

    # views into the swarm_listB array for convenience
    swarm_xListB = swarm_listB[:,0]
    swarm_yListB = swarm_listB[:,1]

    # ---------------------------------------------------------------------
    # ITERATIVELY CLIP SOURCES CONTAINED WITHIN RINGS AT SPECIFIED RADIUS
    # VALUES, PROGRESSIVELY MOVING CLOSER TO THE CENTRAL SOURCE
    # ---------------------------------------------------------------------

    # do the cross-match using xymatch
    log.info('Matching {} swarm centers with {} catalog sources'.format(len(final_flag_src_central_pixel_list),len(swarm_listB)))
    pcentral, pfull = xymatch(final_flag_src_central_pixel_list[:,0:2], swarm_listB[:,0:2], clip_radius_list[0], multiple=True, stack=False, verbose=False)

    #XXX RLW: the ring list is needed only for testing, get rid of it when code works

    if debug:
        ring_index_list = []
        ring_refepp_list = []
        ring_thresh_list = []
        ring_count = []

    for pindex, ii in enumerate(pcentral):

        central_pixel_value = final_flag_src_central_pixel_list[ii,:]
        log.info(' ')
        log.info('CENTRAL PIXEL VALUE: {}'.format(central_pixel_value))
        log.info(' ')

        base_epp = central_pixel_value[3]
        coords = central_pixel_value[0:2]

        allmatches = pfull[pindex]

        if len(allmatches) == 0:
            # (this should not happen using xymatch)
            log.info(' ')
            log.info('------------------------------------------')
            log.info('NOTE: NO SWARM CANDIDATES FOR THIS SOURCE ')
            log.info('------------------------------------------')
            log.info(' ')
            continue

        distsq = (swarm_xListB[allmatches]-coords[0])**2 + (swarm_yListB[allmatches]-coords[1])**2
        sind = distsq.argsort()
        allmatches = allmatches[sind]
        distsq = distsq[sind]
        rcut = distsq.searchsorted(numpy.array(clip_radius_list)**2)
        for radius_cnt in range(1,len(clip_radius_list)):

            # -------------------------------------------
            # ISOLATE THE DETECTIONS WITHIN A GIVEN RING
            # -------------------------------------------

            matches = allmatches[rcut[radius_cnt]:rcut[radius_cnt-1]]

            if len(matches) == 0:
                log.info(' ')
                log.info('------------------------------------------')
                log.info('NOTE: ALL MATCHES/DETECTIONS IN THIS RING ')
                log.info('      HAVE PREVIOUSLY BEEN ACCOUNTED FOR  ')
                log.info('------------------------------------------')
                log.info(' ')

                continue

            # -----------------------------------------------------------
            # CALCULATE THE MEDIAN SKY VALUE FOR THE GROUP OF DETECTIONS
            # CONTAINED WITHIN THE SPECIFIED RING BEING PROCESSED
            # -----------------------------------------------------------
            ref_epp = base_epp * scale_factor_list[radius_cnt-1]

            # -----------------------------------------------------------------------------------
            # DIFFERENTIATE BETWEEN GOOD DETECTIONS AND SWARM DETECTIONS WITHIN SPECIFIED RINGS
            # -----------------------------------------------------------------------------------
            ring = swarm_listB[matches,:]
            w = numpy.where(ring[:,3]/ref_epp < swarm_thresh)
            if len(w) > 0:
                swarm_flag[notcentral_index[matches[w]]] = True


            #XXX RLW: following needed only for testing, get rid of it when code works
            if debug:
                ring_index_list.append(matches)
                ring_count.append(len(matches))
                ring_refepp_list.append(ring[:,3]/ref_epp)
                ring_thresh_list.append(swarm_thresh)

    #XXX RLW: following needed only for testing, get rid of it when code works
    if debug:
        # -----------------------------------------------------------------------------------------
        # WRITE CLIPPED SOURCES CONTAINED WITHIN RINGS TO AN OUTPUT FILE FOR INTERMEDIATE ANALYSIS
        # -----------------------------------------------------------------------------------------
        ring_source_file = phot_table_root+'_SWFILT_RING-SOURCE-INFO.txt'
        ring_src_outfile = open(ring_source_file,'w')
        ring_src_outfile.write("# ------------------------------------------------------------------------------------------------\n")
        ring_src_outfile.write("# X-Center   Y-Center     Flux        ElectronPP          Sky        SrcEPP/RefEPP   Swarm Thresh \n")
        ring_src_outfile.write("# ------------------------------------------------------------------------------------------------\n")

        if ring_index_list:
            ring_index_list = numpy.concatenate(ring_index_list)

            # select just the lowest value of refepp/swarm threshold for each source
            # create array with extra columns
            ring_source_list = numpy.empty((len(ring_index_list),9), dtype=numpy.float)
            ring_source_list[:,0:6] = swarm_listB[ring_index_list,:]
            ring_source_list[:,6] = numpy.concatenate(ring_refepp_list)
            ring_source_list[:,7] = numpy.repeat(ring_thresh_list,ring_count)
            ring_source_list[:,8] = ring_source_list[:,6] / ring_source_list[:,7]

            # sort by x, y, and refepp
            # tricky here: get a view with named columns, then specify names as sort items
            ring_source_list.view(','.join(['f8']*9)).sort(order=['f0','f1','f8'],axis=0)

            # keep just first entry when the same source appears more than once
            keep = numpy.ones(len(ring_index_list),dtype=bool)
            # keep[1:] = numpy.any(ring_source_list[1:,0:2]!=ring_source_list[:-1,0:2], axis=1)
            keep[1:] = numpy.logical_or(ring_source_list[1:,0]!=ring_source_list[:-1,0],
                                        ring_source_list[1:,1]!=ring_source_list[:-1,1])
            ring_source_list = ring_source_list[keep,:]

#                numpy.savetxt(ring_src_outfile,ring_source_list,delimiter='     ')

            for ring_source in ring_source_list:
                ring_src_outfile.write(str(ring_source[0])+'     '+
                                       str(ring_source[1])+'     '+
                                       str(ring_source[2])+'     '+
                                       str(ring_source[3])+'     '+
                                       str(ring_source[4])+'     '+
                                       str(ring_source[5])+'     '+
                                       str(ring_source[6])+'     '+
                                       str(ring_source[7])+'\n')

        ring_src_outfile.close()
        #XXX RLW: end of testing code


    # ===================================================================================
    # -----------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------
    #                          ----- PROXIMITY FILTER -----
    # EXTRACT ADDITIONAL SWARM DETECTIONS BASED ON THE SWARM CANDIDATE CENTRAL POSITIONS,
    # DEFINING THE REMOVAL RADIUS AROUND EACH POSITION BASED ON THAT SOURCE'S EPP
    # -----------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------
    # ===================================================================================

    proximity_flag = numpy.zeros(nrows, dtype=bool)

    proximity_choice = param_dict["quality control"]["swarm filter"]["proximity_binary"]

    if proximity_choice:
        if len(final_flag_src_central_pixel_list) > 0:
            ctrList_radiusList = param_dict["quality control"]["swarm filter"]["ctrList_radiusList"] # TODO: optimize ctrList_radiusList for ACS wfc, hrc, sbc in quality control config files
            ctrList_radiusList = list(map(int, ctrList_radiusList))

            ctrList_thresholdList = param_dict["quality control"]["swarm filter"]["ctrList_thresholdList"] # TODO: optimize ctrList_thresholdList for ACS wfc, hrc, sbc in quality control config files
            ctrList_thresholdList = list(map(int, ctrList_thresholdList))

            for ctrList_cnt,(threshold,radius) in enumerate(zip(ctrList_thresholdList, ctrList_radiusList)):

                if ctrList_cnt == 0:
                    ctr_list_cut = final_flag_src_central_pixel_list[:,3] > threshold
                else:
                    ctr_list_cut = numpy.logical_and(final_flag_src_central_pixel_list[:,3] > threshold,
                                                     final_flag_src_central_pixel_list[:,3] <= ctrList_thresholdList[ctrList_cnt-1])

                ctr_list_cut1 = final_flag_src_central_pixel_list[ctr_list_cut, :]
                pcentral, pfull = xymatch(ctr_list_cut1[:,0:2], swarm_listB[:,0:2], radius, multiple=True, verbose=False)
                proximity_flag[notcentral_index[pfull]] = True

        log.info("Proximity filter flagged {} sources".format(proximity_flag.sum()))

        # --------------------------------------------------------------------------
        # WRITE NEAR CENTRAL POSITION SWARM LIST TO AN OUTPUT FILE FOR VERIFICATION
        # --------------------------------------------------------------------------
        if debug:
            near_swmList = complete_src_list[proximity_flag, :]
            final_near_swarm_file = open(phot_table_root+'_SWFILT_NEAR_SWARM_FILE.txt','w')
            for swarm_value in near_swmList:
                final_near_swarm_file.write(str(swarm_value[0])+'     '+
                                       str(swarm_value[1])+'     '+
                                       str(swarm_value[2])+'     '+
                                       str(swarm_value[3])+'     '+
                                       str(swarm_value[4])+'\n')
            final_near_swarm_file.close()

    # -------------------------------------------------------------------------
    # EXTRACT DETECTIONS FROM THE complete_src_list THAT ARE NOT FLAGGED
    # -------------------------------------------------------------------------

    combined_flag = numpy.logical_or(swarm_flag,proximity_flag)
    final_swarm_list = complete_src_list[combined_flag, :]
    final_source_list = complete_src_list[numpy.logical_not(combined_flag), :]

    log.info(' ')
    log.info('************************************************')
    log.info('INITIAL LENGTH OF complete_src_list = {}'.format(len(complete_src_list)))
    log.info(' ')
    log.info('LENGTH OF final_source_list = {}'.format(len(final_source_list)))
    log.info('LENGTH OF final_swarm_list = {}'.format(len(final_swarm_list)))
    log.info('TOTAL LENGTH = {}'.format(len(final_source_list) + len(final_swarm_list)))
    log.info(' ')
    log.info('MEDIAN SKY VALUE = {}'.format(median_sky))
    log.info('************************************************')
    log.info(' ')

    # ----------------------------------------------------
    # WRITE SWARM LIST TO AN OUTPUT FILE FOR VERIFICATION
    # ----------------------------------------------------
    if debug:
        final_swarm_file = open(phot_table_root+'_SWFILT_SWARM_FILE.txt','w')
        for swarm_value in final_swarm_list:
            final_swarm_file.write(str(swarm_value[0])+'     '+
                                   str(swarm_value[1])+'     '+
                                   str(swarm_value[2])+'     '+
                                   str(swarm_value[3])+'     '+
                                   str(swarm_value[4])+'\n')
        final_swarm_file.close()

    # ----------------------------------------------------
    # WRITE SOURCE LIST TO AN OUTPUT FILE FOR VERIFICATION
    # ----------------------------------------------------
    if debug:
        final_source_file = open(phot_table_root+'_SWFILT_SOURCE_FILE.txt','w')
        for source_value in final_source_list:
            final_source_file.write(str(source_value[0])+'     '+
                                    str(source_value[1])+'     '+
                                    str(source_value[2])+'     '+
                                    str(source_value[3])+'     '+
                                    str(source_value[4])+'\n')
        final_source_file.close()



    # Update catalog_data flag values
    for i,table_row in enumerate(catalog_data[0:]):
        if combined_flag[i]:
            table_row[-1] |= 32

    if debug:
        # Write out intermediate catalog with updated flags
        phot_table_temp = phot_table_root + '_SWFILT.txt'
        catalog_data.write(phot_table_temp, delimiter=",",format='ascii')

    return catalog_data

# ======================================================================================================================

def HLANexpFlags(drizzled_image, flt_list, param_dict, catalog_name, catalog_data, drz_root_dir, debug):
    """flags out sources from regions where there are a low (or a null) number of contributing exposures
   
    drizzled_image : string
        Name of drizzled image to process

    flt_list : list
        list of calibrated images that were drizzle-combined to produce image specified by input parameter
        'drizzled_image'

    param_dict : dictionary
        Dictionary of instrument/detector - specific drizzle, source finding and photometric parameters.

    catalog_name : string
        drizzled filter product catalog filename to process

    catalog_data : astropy.Table object
        drizzled filter product catalog data to process

    drz_root_dir : dictionary of source lists keyed by drizzled image name.

    debug : bool
        write intermediate files?

    Returns
    -------
    catalog_data : astropy.Table object
        drizzled filter product catalog data with updated flag values
    
    """
    # ------------------
    # CREATE NEXP IMAGE
    # ------------------
    image_split = drizzled_image.split('/')[-1]
    channel = drizzled_image.split("_")[-3].upper() # TODO: May need to be refactored to adjust for new names, and fact that ACS has two filters

    #if channel == 'IR':
    #    continue

    # ---------
    # METHOD 1:
    # ---------
    ctx = getdata(drizzled_image, 3)

    if channel in ['UVIS','IR','WFC','HRC']: ncombine = getheader(drizzled_image,1)['NCOMBINE']
    if channel in ['WFPC2','PC']:
        ndrizim=getheader(drizzled_image,0)['NDRIZIM']
        ncombine=ndrizim/4
    if channel == 'SBC':
        ndrizim=getheader(drizzled_image,0)['NDRIZIM']
        ncombine = ndrizim
    ctxarray = arrayfy_ctx(ctx, ncombine)

    nexp_array_ctx = ctxarray.sum(axis=-1)
    nexp_image_ctx = drizzled_image.split('.')[0]+'_NCTX.fits'
    if not os.path.isfile(nexp_image_ctx):
        hdr = getheader(drizzled_image, 1)
        pyfits.writeto(nexp_image_ctx, numpy.float32(nexp_array_ctx), hdr)

    # ---------
    # METHOD 2:
    # ---------
    drz_data = getdata(drizzled_image, 1)

    ## this bit is added to get the mask integrated into the exp map
    maskfile = drizzled_image.replace('_drz.fits','_msk.fits')
    if os.path.isfile(maskfile):
        mask_data = getdata(maskfile)
        mask_array = (mask_data==0.0).astype(numpy.int32)

    component_drz_img_list = get_component_drz_list(drizzled_image, drz_root_dir, flt_list)# TODO: This might be a problem for HAP
    nx = drz_data.shape[0]
    ny = drz_data.shape[1]
    nexp_array = numpy.zeros((nx, ny), dtype = numpy.int32)

    for comp_drz_img in component_drz_img_list:
        comp_drz_data = (getdata(comp_drz_img) != 0).astype(numpy.int32)
        try:
            nexp_array += comp_drz_data
        except ValueError:
            log.info("WARNING: Astrodrizzle added an extra-row/column...")
            nexp_array += comp_drz_data[0:nx,0:ny]

    if os.path.isfile(maskfile):
        nexp_array = nexp_array * mask_array
    else:
        log.info("something's wrong: maskfile {} is not a file".format(maskfile))
        sys.exit()
    nexp_image = drizzled_image.split('.')[0]+'_NEXP.fits'
    if not os.path.isfile(nexp_image):
        hdr = getheader(drizzled_image, 1)
        pyfits.writeto(nexp_image, numpy.float32(nexp_array), hdr)

    # -------------------------------------------------------
    # EXTRACT FLUX/NEXP INFORMATION FROM NEXP IMAGE BASED ON
    # THE SOURCE DETECTION POSITIONS PREVIOUSLY ESTABLISHED
    # -------------------------------------------------------

    phot_table_root = catalog_name.split('/')[-1].split('.')[0]

    nrows = len(catalog_data)
    cat_coords = numpy.empty((nrows,2),dtype=float)
    for line_cnt,phot_table_line in enumerate(catalog_data):
        x_coord = phot_table_line[0]
        y_coord = phot_table_line[1]
        cat_coords[line_cnt, :] = [x_coord, y_coord]
    # ----------------------------------
    # Convert aperture radius to pixels
    # ----------------------------------

    ap2 = param_dict['catalog generation']['aperture_2']
    plate_scale = param_dict['catalog generation']['scale'] #TODO: Need to move scale value out of 'catalog generation' > 'dao' to somewhere more sensable
    radius = ap2/plate_scale

    num_exp = round(numpy.max(nexp_array))
    if num_exp<=1 or channel in ('IR','SBC'):
        # Keep everything that has an exposure for detectors without CRs or
        # when there is only one exposure
        artifact_filt = 0.5
    elif num_exp > 5:
        # Flag sources with <= 2 exposures when there are > 5 total
        # We are always using the 'imedian' combination in that case, and it
        # does not do a very good job of rejecting CRs with only 2 available
        # exposures
        artifact_filt = 2.5
    else:
        artifact_filt = 1.5

    icoords = (cat_coords+0.5).astype(int)
    # note x & y are swapped so they can index the numpy array nexp_array
    # catalog x is second subscript, catalog y is first subscript
    ix = icoords[:,1]
    iy = icoords[:,0]

    # get list of neighboring pixels that are within radius
    iradius = int(radius+1)
    idiam = iradius*2+1
    gx, gy = numpy.mgrid[0:idiam,0:idiam] - iradius
    gx = gx.ravel()
    gy = gy.ravel()
    w = numpy.where(gx**2+gy**2 <= radius**2)[0]
    gx = gx[w]
    gy = gy[w]

    # check the pixel values for low nexp

    # this version uses numpy broadcasting sum gx+ix is [len(gx),nrows]
    gx = (gx[:,numpy.newaxis] + ix).clip(0,nexp_array.shape[0]-1)
    gy = (gy[:,numpy.newaxis] + iy).clip(0,nexp_array.shape[1]-1)
    artifact_flag = nexp_array[gx,gy].min(axis=0) < artifact_filt

    log.info('FLAGGING {} OF {} SOURCES'.format(artifact_flag.sum(),nrows))

    # -------------------------------------------------------------------
    # WRITE NEXP FLAGS TO OUTPUT PHOT TABLE BASED ON nexp_phot_data_list
    # -------------------------------------------------------------------
    #nexp_outfile_good = open('temp_outfile_NEXP_GOOD.txt','w')
    #nexp_outfile_bad = open('temp_outfile_NEXP_BAD.txt','w')

    # Add flag bit to appropriate sources
    for i, table_row in enumerate(catalog_data):
        if artifact_flag[i]:
            table_row[-1] |= 64

    if debug:
        # Write out intermediate catalog with updated flags
        phot_table_temp = phot_table_root + '_NEXPFILT.txt'
        catalog_data.write(phot_table_temp, delimiter=",", format='ascii')


    return catalog_data


# +++++++++++++++++++++++++++++++++++++++ OLD VERSIONS +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# These subroutines are the older versions of the subroutines that did not use persistant in-memory storage of the catalogs.
def run_source_list_flaging_OLD(all_drizzled_filelist, working_hla_red, filter_sorted_flt_dict, param_dict,
                                readnoise_dictionary_drzs, scale_dict_drzs, zero_point_AB_dict, exp_dictionary_scis,
                                detection_image, dict_newTAB_matched2drz, phot_table_matched2drz, proc_type,
                                drz_root_dir,rms_dict):
    """Simple calling subroutine that executes the other flagging subroutines. OLDER FILE-BASED VERSIONS!

    Parameters
    ----------
    all_drizzled_filelist : list
        List of drizzled images to process

    working_hla_red : string
        full path to working directory

    filter_sorted_flt_dict : dictionary
        dictionary containing lists of calibrated images sorted (also keyed) by filter name.

    param_dict : dictionary
        Dictionary of instrument/detector - specific drizzle, source finding and photometric parameters

    readnoise_dictionary_drzs : dictionary
        dictionary of readnoise values keyed by drizzled image.

    scale_dict_drzs : dictionary
        dictionary of scale values keyed by drizzled image name.

    zero_point_AB_dict : dictionary
        ***UNUSED*** dictionary of zeropoint values keyed by drizzled image name.

    exp_dictionary_scis : dictionary
        dictionary of exposure time values keyed by drizzled image name.

    detection_image : string
        name of drizzled image.

    dict_newTAB_matched2drz : dictionary
        dictionary of source lists keyed by drizzled image name.

    phot_table_matched2drz : dictionary
        dictionary of source lists tables (already read into memory) keyed by drizzled image name.

    proc_type : string
        sourcelist generation type.

    drz_root_dir : string
        Root directory of drizzled images.

    rms_dict : dictionary
        dictionary of RMS image counterparts to drizzled images. Keyed by drizzled image name.

    Returns
    -------
    Nothing!
    """
    ci_filter_OLD(all_drizzled_filelist, dict_newTAB_matched2drz,working_hla_red, proc_type, param_dict)

    HLASaturationFlags_OLD(all_drizzled_filelist, working_hla_red, filter_sorted_flt_dict, readnoise_dictionary_drzs,
                           scale_dict_drzs, exp_dictionary_scis, dict_newTAB_matched2drz, proc_type, param_dict)

    HLASwarmFlags_OLD(all_drizzled_filelist, dict_newTAB_matched2drz, working_hla_red, exp_dictionary_scis,
                      filter_sorted_flt_dict, detection_image, proc_type, rms_dict, param_dict)

    HLANexpFlags_OLD(all_drizzled_filelist, working_hla_red, filter_sorted_flt_dict, param_dict, readnoise_dictionary_drzs,
                 scale_dict_drzs, exp_dictionary_scis, dict_newTAB_matched2drz, drz_root_dir)

# ======================================================================================================================

def ci_filter_OLD(all_drizzled_filelist, dict_newTAB_matched2drz, working_hla_red, proc_type, param_dict):
    """This subroutine flags sources based on concentration index.  Sources below the minimum CI value are
    flagged as hot pixels/CRs (flag=16). Sources above the maximum (for stars) are flagged as extended (flag=1).
    It also flags sources below the detection limit in mag_aper2 (flag=8).

    Parameters
    ----------
    all_drizzled_filelist : list
        list of drizzled images

    dict_newTAB_matched2drz : dictionary
        dictionary of source lists keyed by drizzled image name.

    working_hla_red : string
        ***UNUSED*** full path of working directory.

    proc_type : string
        Sourcelist generation type

    param_dict : dictionary
        Dictionary of instrument/detector - specific drizzle, source finding and photometric parameters

    Returns
    -------
    Nothing!
    """

    # column indices for SE and DAO catalogs
    if proc_type == 'sexphot':
        imag1 = 5
        imag2 = 6
        imerr1 = 7
        imerr2 = 8
    elif proc_type == 'daophot':
        imag1 = 5
        imag2 = 7
        imerr1 = 6
        imerr2 = 8
    else:
        raise ValueError("Unknown proc_type '%s', must be 'sexphot' or 'daophot'" % (proc_type,))

    for drizzled_image in all_drizzled_filelist:
        phot_table = dict_newTAB_matched2drz[drizzled_image]
        phot_table_root = phot_table.split('.')[0]
        if proc_type == 'sexphot':
            ci_lower_limit = float(param_dict['quality control']['ci filter']['ci_selower_limit'])
            ci_upper_limit = float(param_dict['quality control']['ci filter']['ci_seupper_limit'])
            snr = float(param_dict['catalog generation']['sourcex']['bthresh'])

        if proc_type == 'daophot':
            ci_lower_limit = float(param_dict['quality control']['ci filter']['ci_daolower_limit'])
            ci_upper_limit = float(param_dict['quality control']['ci filter']['ci_daoupper_limit'])
            snr = float(param_dict['catalog generation']['dao']['bthresh'])


        # replace CI limits with values from table if possible
        cidict = ci_table.get_ci_from_file(drizzled_image, ci_lower=ci_lower_limit, ci_upper=ci_upper_limit)
        ci_lower_limit = cidict['ci_lower_limit']
        ci_upper_limit = cidict['ci_upper_limit']

        log.info(' ')
        log.info('ci limits for {}'.format(drizzled_image))
        log.info('ci_lower_limit = {}'.format(ci_lower_limit))
        log.info('ci_upper_limit = {}'.format(ci_upper_limit))
        log.info(' ')

        phot_table = dict_newTAB_matched2drz[drizzled_image]
        phot_table_root = phot_table.split('.')[0]

        phot_table_in = open(phot_table,'r')
        phot_table_rows = phot_table_in.readlines()
        phot_table_in.close()

        phot_table_temp = phot_table_root+'_temp.txt'
        phot_table_out = open(phot_table_temp,'w')
        failed_ci_table_out = open(phot_table_root+'_Failed-CI.txt','w')

        for i,table_row in enumerate(phot_table_rows):

            if i == 0:
                phot_table_out.write(table_row)
                continue
            else:
                row_split = table_row.split(',')
                try:
                    flag_value = int(row_split[-1])
                except ValueError:
                    flag_value = 0

                ci_value = row_split[-2]
                if ci_value != '':
                    ci_value = float(ci_value)
                merr1 = row_split[imerr1]
                if merr1 == '':
                    merr1 = numpy.nan
                else:
                    merr1 = float(merr1)
                merr2 = row_split[imerr2]
                if merr2 == '':
                    merr2 = numpy.nan
                else:
                    merr2 = float(merr2)
                good_snr = merr2 <= 2.5/(snr*numpy.log(10))
                ci_err = numpy.sqrt(merr1**2+merr2**2)

                if not good_snr:
                    flag_value |= 8

                if ci_value == '' or (not numpy.isfinite(ci_err)) or ci_value < ci_lower_limit - ci_err:
                    flag_value |= 16
                    # print(i,ci_value)
                    # input()

                if ci_value != '':
                    if ci_value > ci_upper_limit:
                        flag_value |= 1


                row_split[-1] = '%d\n' % flag_value
                table_row = ','.join(row_split)
                phot_table_out.write(table_row)
                if ci_value == '':
                    failed_ci_table_out.write(table_row)

        phot_table_out.close()
        failed_ci_table_out.close()

        os.system('mv '+phot_table+' '+phot_table+'.PreCIFilt')
        os.system('mv '+phot_table_temp+' '+phot_table)

# ======================================================================================================================

def HLASaturationFlags_OLD(all_drizzled_filelist, working_hla_red, filter_sorted_flt_dict, readnoise_dictionary_drzs,
                           scale_dict_drzs, exp_dictionary_scis, dict_newTAB_matched2drz, proc_type, param_dict):

    """Identifies and flags saturated sources.

    Parameters
    ----------
    all_drizzled_filelist : list
        List of drizzled images to process.

    working_hla_red : string
        ***UNUSED*** full path to working directory

    filter_sorted_flt_dict : dictionary
        dictionary containing lists of calibrated images sorted (also keyed) by filter name.

    readnoise_dictionary_drzs : dictionary
        ***UNUSED*** dictionary of readnoise values keyed by drizzled image.

    scale_dict_drzs : dictionary
        ***UNUSED*** dictionary of scale values keyed by drizzled image.

    exp_dictionary_scis : dictionary
        ***UNUSED*** dictionary of exposure time values keyed by drizzled image.

    dict_newTAB_matched2drz : dictionary
        dictionary of source lists keyed by drizzled image name.

    proc_type : string
        sourcelist generation type.

    param_dict : dictionary
        Dictionary of instrument/detector - specific drizzle, source finding and photometric parameters

    Returns
    -------
    Nothing!
    """
    for drizzled_image in all_drizzled_filelist:
        image_split = drizzled_image.split('/')[-1]
        channel = drizzled_image.split("_")[-3].upper()

        if channel == 'IR':
            continue

        phot_table = dict_newTAB_matched2drz[drizzled_image]
        phot_table_root = phot_table.split('.')[0]

#        for flt_image in flt_images:
#            os.system('cp '+flt_image+' .')

        # -------------------------------------------------------------------
        # STEP THROUGH EACH APPLICABLE FLT IMAGE, DETERMINE THE COORDINATES
        # FOR ALL SATURATION FLAGGED PIXELS, AND TRANSFORM THESE COORDINATES
        # INTO THE DRIZZLED IMAGE REFERENCE FRAME.
        # -------------------------------------------------------------------
        main_drizzled_filelist = [drizzled_image]
        main_drizzled_filelist_orig = [drizzled_image]

        drz_filter = drizzled_image.split("_")[5]
        list_of_flts_in_main_driz = filter_sorted_flt_dict[drz_filter.lower()]
        num_flts_in_main_driz = len(list_of_flts_in_main_driz)
        list_of_flts_in_main_driz.sort()

        log.info(' ')
        log.info("Current Working Directory: {}".format(os.getcwd()))
        log.info(' ')
        log.info('LIST OF FLTS IN {}: {}'.format(drizzled_image.split('/')[-1],list_of_flts_in_main_driz))
        log.info(' ')
        log.info('NUMBER OF FLTS IN {}: {}'.format(drizzled_image.split('/')[-1],num_flts_in_main_driz))
        log.info(' ')

        # ----------------------------------------------------
        # EXTRACT DQ DATA FROM FLT IMAGE AND CREATE A LIST
        # OF "ALL" PIXEL COORDINATES WITH A FLAG VALUE OF 256
        # ----------------------------------------------------
        if ((channel.lower() != 'wfpc2') and (channel.lower() != 'pc')):
            image_ext_list = ["[sci,1]","[sci,2]"]
            dq_sat_bit=256
        if channel.lower() == 'wfpc2':
            image_ext_list = ["[sci,1]","[sci,2]","[sci,3]","[sci,4]"]
            dq_sat_bit = 8
        if channel.lower() == 'pc':
            image_ext_list = ["[sci,1]"]
            dq_sat_bit = 8

        # build list of arrays
        drz_sat_xy_coords_list = []

        for flt_cnt,flt_image in enumerate(list_of_flts_in_main_driz):
            for ext_cnt,image_ext in enumerate(image_ext_list):
                ext_part = image_ext.split(',')[1].split(']')[0]
                try:
                    if ((channel.lower() != 'wfpc2') and (channel.lower() != 'pc')): flt_data = getdata(flt_image,'DQ',int(ext_part))
                    if ((channel.lower() == 'wfpc2') or (channel.lower() == 'pc')): flt_data = getdata(flt_image.replace("_c0m","_c1m"),'SCI',int(ext_part))
                except KeyError:
                    log.info(' ')
                    log.info('WARNING: There is only one set of file extensions in {}'.format(flt_image))
                    log.info(' ')

                    continue

                # ----------------------------------------------------
                # DETERMINE IF ANY OF THE PIXELS LOCATED IN THE GRID
                # HAVE A BIT VALUE OF 256, I.E. FULL WELL SATURATION.
                # ----------------------------------------------------
                # NOTE: NUMPY ARRAYS REPORT Y COORD VALUES FIRST AND
                #       X COORD VALUES SECOND AS FOLLOWS:
                #
                #       --> numpy.shape(flt_data)
                #       (2051, 4096)
                #
                #       WHERE 2051 IS THE NUMBER OF PIXELS IN THE Y
                #       DIRECTION, AND 4096 IS THE NUMBER OF PIXELS
                #       IN THE X DIRECTION.
                # ----------------------------------------------------
                bit_flt_data = dq_sat_bit & flt_data
                complete_sat_coords = numpy.where(bit_flt_data == dq_sat_bit)

                if len(complete_sat_coords[0]) == 0:
                    continue

                # -------------------------------------------------
                # RESTRUCTURE THE LIST OF X AND Y COORDINATES FROM
                # THE FLT FILE THAT HAVE BEEN FLAGGED AS SATURATED
                # -------------------------------------------------
                nsat = len(complete_sat_coords[0])
                x_y_array = numpy.empty((nsat,2),dtype=int)
                x_y_array[:,0] = complete_sat_coords[1]
                x_y_array[:,1] = complete_sat_coords[0]

                # ---------------------------------------------------
                # WRITE FLT COORDS TO A FILE FOR DIAGNOSTIC PURPOSES
                # ---------------------------------------------------
                flt_xy_coord_out = flt_image.split('/')[-1].split('.')[0]+'_sci'+str(ext_cnt+1)+'.txt'
                outfile = open(flt_xy_coord_out,'w')
                for flt_xy_coord in x_y_array:
                    x = flt_xy_coord[0]
                    y = flt_xy_coord[1]
                    outfile.write(str(x)+'     '+str(y)+'\n')
                outfile.close()

                # ----------------------------------------------------
                # CONVERT SATURATION FLAGGED X AND Y COORDINATES FROM
                # THE FLT IMAGE INTO RA AND DEC
                # ----------------------------------------------------
                flt_ra_dec_coords = xytord(x_y_array, flt_image, image_ext)

                # -------------------------------------------------
                # CONVERT RA & DEC VALUES FROM FLT REFERENCE FRAME
                # TO THAT OF THE DRIZZLED IMAGE REFERENCE FRAME
                # -------------------------------------------------
                drz_sat_xy_coords_list.append(rdtoxy(flt_ra_dec_coords, drizzled_image, "[sci,1]"))

                log.info(' ')
                log.info('FLT IMAGE = {}'.format(flt_image.split('/')[-1]))
                log.info('IMAGE EXT = {}'.format(image_ext))
                log.info(' ')

        # ----------------------------------------------------------------
        # IF NO SATURATION FLAGS EXIST IN ANY OF THE FLT FILES, THEN SKIP
        # ----------------------------------------------------------------
        if len(drz_sat_xy_coords_list) == 0:
            log.info(' ')
            log.info('*******************************************************************************************')
            log.info('NO SATURATION FLAGGED PIXELS EXIST IN ANY OF THE FLT FILES FOR:')
            log.info('     --> {}'.format(drizzled_image.split('/')[-1]))
            log.info('*******************************************************************************************')
            log.info(' ')

            continue

        # ------------------------------
        # now concatenate all the arrays
        # ------------------------------
        full_satList = numpy.concatenate(drz_sat_xy_coords_list)

        # --------------------------------------------
        # WRITE RA & DEC FLT CONVERTED X & Y DRIZZLED
        # IMAGE COORDINATES TO A TEXT FILE
        # --------------------------------------------
        drz_coord_file = drizzled_image.split('/')[-1].split('.')[0]+'_ALL_FLT_SAT_FLAG_PIX.txt'
        drz_coord_out = open(drz_coord_file,'w')
        for coord in full_satList:
            drz_coord_out.write(str(coord[0])+'     '+str(coord[1])+'\n')
        drz_coord_out.close()

        # --------------------------------------------------------
        # READ IN FULL DRIZZLED IMAGE-BASED CATALOG AND SAVE
        # X AND Y COORDINATE VALUES TO LISTS FOR LATER COMPARISON
        # --------------------------------------------------------
        full_drz_cat = dict_newTAB_matched2drz[drizzled_image]
        inputfile=open(full_drz_cat,'r')
        all_detections=inputfile.readlines()
        inputfile.close()

        nrows = len(all_detections)-1
        full_coordList = numpy.empty((nrows,2), dtype=numpy.float)
        for row_count,detection in enumerate(all_detections[1:]):
            ss = detection.split(',')
            full_coordList[row_count,0] = float(ss[0])
            full_coordList[row_count,1] = float(ss[1])

        # ----------------------------------------------------
        # CREATE SUB-GROUPS OF SATURATION-FLAGGED COORDINATES
        # ----------------------------------------------------
        proc_time1=time.ctime()
        log.info(' ')
        log.info('PROC_TIME_1: {}'.format(proc_time1))
        log.info(' ')

        # ----------------------------------
        # Convert aperture radius to pixels
        # ----------------------------------
        ap2 = param_dict['catalog generation']['aperture_2']
        if proc_type == 'daophot':
            if channel == 'IR':
                radius = round((ap2 / 0.09) + 0.5) * 2.
            if channel == 'UVIS':
                radius = round((ap2 / 0.04) + 0.5) * 2.
            if channel == 'WFC':
                radius = round((ap2 / 0.05) + 0.5) * 2.
            if channel == 'HRC':
                radius = round((ap2 / 0.027) + 0.5) * 2.
            if channel == 'WFPC2':
                radius = round((ap2 / 0.1) + 0.5) * 2.
            if channel == 'PC':
                radius = round((ap2 / 0.046) + 0.5) * 2.


        if proc_type == 'sexphot':
            if channel == 'IR':
                radius = round((ap2 / 0.09) + 0.5)
            if channel == 'UVIS':
                radius = round((ap2 / 0.04) + 0.5)
            if channel == 'WFC':
                radius = round((ap2 / 0.05) + 0.5)
            if channel == 'HRC':
                radius = round((ap2 / 0.027) + 0.5)
            if channel == 'WFPC2':
                radius = round((ap2 / 0.1) + 0.5)
            if channel == 'PC':
                radius = round((ap2 / 0.046) + 0.5) * 2.


        log.info(' ')
        log.info('THE RADIAL DISTANCE BEING USED IS {} PIXELS'.format(str(radius)))
        log.info(' ')

        # do the cross-match using xymatch
        log.info('Matching {} saturated pixels with {} catalog sources'.format(len(full_satList),len(full_coordList)))
        psat, pfull = xymatch(full_satList, full_coordList, radius, multiple=True, verbose=False)
        log.info('Found cross-matches (including duplicates)'.format(len(psat)))
        saturation_flag = numpy.zeros(len(full_coordList),dtype=bool)
        saturation_flag[pfull] = True

        proc_time2=time.ctime()
        log.info(' ')
        log.info('PROC_TIME_2: {}'.format(proc_time2))
        log.info(' ')

        # ------------------------------------------------------------------
        # REMOVE DUPLICATE DETECTIONS FROM THE LIST, "group", CREATTED FROM
        # MATCHING SATURATION FLAGGED FLT PIXELS TO FINAL SOURCE DETECTIONS
        # ------------------------------------------------------------------

        nsaturated = saturation_flag.sum()
        if nsaturated == 0:
            log.info(' ')
            log.info('**************************************************************************************')
            log.info('NOTE: NO SATURATED SOURCES WERE FOUND FOR: {}'.format(image_split))
            log.info('**************************************************************************************')
            log.info(' ')

            continue

        else:
            log.info(' ')
            log.info('FLAGGED {} SOURCES'.format(nsaturated))
            log.info(' ')

            sat_coord_file = drizzled_image.split('/')[-1].split('.')[0]+'_INTERMEDIATE.txt'
            sat_coord_out = open(sat_coord_file,'w')
            for sat_coord in full_coordList[saturation_flag,:]:
                sat_coord_out.write(str(sat_coord[0])+'     '+str(sat_coord[1])+'\n')
            sat_coord_out.close()

            # --------------------------------------------------------------------------
            # WRITE SAT FLAGS TO OUTPUT PHOT TABLE BASED ON flag_src_central_pixel_list
            # --------------------------------------------------------------------------
            phot_table = dict_newTAB_matched2drz[drizzled_image]
            phot_table_root = phot_table.split('.')[0]

            phot_table_in = open(phot_table,'r')
            phot_table_rows = phot_table_in.readlines()
            phot_table_in.close()

            phot_table_temp = phot_table_root+'_SATFILT.txt'
            phot_table_out = open(phot_table_temp,'w')

            phot_table_out.write(phot_table_rows[0])
            for i,table_row in enumerate(phot_table_rows[1:]):
                if saturation_flag[i]:
                    row_split = table_row.split(',')
                    sat_flag = int(row_split[-1]) | 4
                    row_split[-1] = str(sat_flag)+'\n'
                    table_row = ','.join(row_split)
                phot_table_out.write(table_row)

            phot_table_out.close()

            os.system('mv '+phot_table+' '+phot_table+'.PreSatFilt')
            os.system('mv '+phot_table_temp+' '+phot_table)

            log.info(' ')
            log.info('FINAL SAT-FILT PHOT_TABLE: {}'.format(phot_table))
            log.info(' ')
            HLA_flag4and8_hunter_killer_OLD(phot_table)

# ======================================================================================================================

def HLASwarmFlags_OLD(all_drizzled_filelist, dict_newTAB_matched2drz, working_hla_red, exp_dictionary_scis,
                      filter_sorted_flt_dict, detection_image, proc_type, rms_dict, param_dict):
    """Identifies and flags swarm sources.

    Parameters
    ----------
    all_drizzled_filelist : list
        List of drizzled images to process

    dict_newTAB_matched2drz : dictionary
        dictionary of source lists keyed by drizzled image name.

    working_hla_red : string
        full path to working directory.

    exp_dictionary_scis : dictionary
        dictionary of exposure time values keyed by drizzled image.

    filter_sorted_flt_dict : dictionary
        dictionary containing lists of calibrated images sorted (also keyed) by filter name.

    detection_image : string
        Name of multi-filter composite 'detection' image

    proc_type : string
        sourcelist generation type.

    rms_dict : dictionary
        dictionary of RMS image counterparts to drizzled images. Keyed by drizzled image name.

    param_dict : dictionary
        Dictionary of instrument/detector - specific drizzle, source finding and photometric parameters

    Returns
    -------
    Nothing!
    """
    for drizzled_image in all_drizzled_filelist:

        drz_imgPath_split = drizzled_image.split('/')
        drz_img_split = drz_imgPath_split[-1].split('_')
        data_type = drz_img_split[4]
        drz_filter = drz_img_split[5]

        # ============================= SKY COMPUTATION ==============================
        # ----------------------------------------------------------------------------
        # In the pipeline, AstroDrizzle sky computation is based on the statistical
        # distribution of pixel values in an input image. It is performed by
        # iterative sigma-clipping, starting with the full range of pixel intensity
        # values, to calculate the standard deviation. By default, pixel values
        # deviating from the median value (specified by parameter skystat) by over
        # four sigma are rejected, and this operation is repeated for a total of five
        # iterations. The median value of the final distribution is used as the sky
        # value.
        # ----------------------------------------------------------------------------
        # ============================= SKY COMPUTATION ==============================

        median_sky = get_median_sky(drizzled_image)
        # single_rms = rms_dict[drizzled_image]
        # rms_array = pyfits.getdata(single_rms,0)
        # rms_subarray = rms_array[rms_array > 0.0]
        # median_sky = Util.binmode(rms_subarray[rms_subarray < 5000.])[0]

        log.info(' ')
        log.info('MEDIAN SKY VALUE = {}'.format(median_sky))
        log.info(' ')

        # ==========================================
        # ------------------------------------------
        # ------------------------------------------
        # CREATE LIST OF POTENTIAL SWARM DETECTIONS
        # ------------------------------------------
        # ------------------------------------------
        # ==========================================
        phot_table = dict_newTAB_matched2drz[drizzled_image]
        phot_table_root = phot_table.split('/')[-1].split('.')[0]

        image_split = drizzled_image.split('/')[-1]
        channel = drizzled_image.split("_")[
            -3].lower()

        ap2 = param_dict['catalog generation']['aperture_2']
        if proc_type not in ('sexphot', 'daophot'):
            raise ValueError("Unknown catalog type '%s'" % proc_type)

        # ----------------------------------
        # Convert aperture radius to pixels
        # ----------------------------------
        radius = ap2 / float(
            param_dict['catalog generation']['scale'])
        log.info(' ')
        log.info('Aperture Size = {}'.format(ap2))
        log.info('Pixel Scale = {} arcsec per pixel'.format(float(
            param_dict['catalog generation']['scale'])))
        log.info(' ')
        area = math.pi * radius ** 2
        exptime = exp_dictionary_scis[drizzled_image]

        log.info('Reading catalog from{}'.format(phot_table))
        phot_table_in = open(phot_table, 'r')
        phot_table_rows = phot_table_in.readlines()
        phot_table_in.close()

        nrows = len(phot_table_rows) - 1

        complete_src_list = numpy.empty((nrows, 6), dtype=numpy.float)

        for row_num, row in enumerate(phot_table_rows[1:]):

            row_split = row.split(',')
            x_val = float(row_split[0])
            y_val = float(row_split[1])

            if proc_type == 'sexphot':
                # mag = row_split[6]
                flux = row_split[10]
                sky = row_split[13]
            elif proc_type == 'daophot':
                # mag = row_split[7]
                flux = row_split[11]
                sky = row_split[9]

            if flux.strip():
                flux = float(flux)
            else:
                flux = 0.0

            if sky.strip():
                sky = float(sky)
            else:
                sky = 0.0

            electronpp = flux / area * exptime
            eppsky = electronpp / median_sky
            complete_src_list[row_num, :] = [x_val, y_val, flux, electronpp, sky, eppsky]

        if len(complete_src_list) == 0:
            continue

        # view into the complete_src_list array for convenience
        swarm_epp_listA = complete_src_list[:, 3]

        # swarm flag array
        # this will get set as candidates to flag are found
        swarm_flag = numpy.zeros(nrows, dtype=bool)

        # ------------------------------------------------------------
        # WRITE SUBSET SOURCE LIST TO AN OUTPUT FILE FOR VERIFICATION
        # ------------------------------------------------------------
        final_complete_source_file = open(phot_table_root + '_SWFILT_COMPLETE_SOURCE_FILE.txt', 'w')
        final_complete_source_file.write(
            "# ------------------------------------------------------------------------------------------------\n")
        final_complete_source_file.write(
            "# X-Center   Y-Center     Flux        ElectronPP          Sky         EPPSKY_Ratio \n")
        final_complete_source_file.write(
            "# ------------------------------------------------------------------------------------------------\n")
        for i, complete_src_value in enumerate(complete_src_list):
            final_complete_source_file.write(
                str(complete_src_value[0]) + '     ' + str(complete_src_value[1]) + '     ' + str(
                    complete_src_value[2]) + '     ' + str(complete_src_value[3]) + '     ' + str(
                    complete_src_value[4]) + '     ' + str(complete_src_value[5]) + '\n')

        final_complete_source_file.close()

        # ======================================================================
        # ----------------------------------------------------------------------
        # ----------------------------------------------------------------------
        # Introduce 2 thresholds:
        # -----------------------
        # A minimum electronpp, and a minimum electronpp/sky.
        # The thresholds should have different values for IR and UVIS.
        #
        # For IR, sources that have electronpp > 100k, OR
        # ((electronpp > 100*sky) AND (electronpp > 10k)), should be considered.
        #
        # For UVIS, I would set the thresholds at (electronpp > 100k) OR
        # ((electronpp > 1000*sky) AND (electronpp > 10k)).
        # ----------------------------------------------------------------------
        # ----------------------------------------------------------------------
        # ======================================================================
        upper_epp_limit = float(param_dict["quality control"]["swarm filter"]["upper_epp_limit"])
        lower_epp_limit = float(param_dict["quality control"]["swarm filter"]["lower_epp_limit"])
        eppsky_limit_cfg = float(param_dict["quality control"]["swarm filter"]["eppsky_limit"])

        if data_type.upper() == 'UVIS':
            eppsky_limit = eppsky_limit_cfg * median_sky
            selfradius = 20.0

        if data_type.upper() == 'IR':
            eppsky_limit = eppsky_limit_cfg * median_sky
            selfradius = 10.0

        if data_type.upper() == 'WFC':  # acs/wfc
            eppsky_limit = eppsky_limit_cfg * median_sky
            selfradius = 20.0

        if data_type.upper() == 'HRC':
            eppsky_limit = eppsky_limit_cfg * median_sky
            selfradius = 20.0  # JUST USING ACS/WFC VALUE. PROBABLY NEEDS TO BE OPTIMIZED FOR HRC.

        if data_type.upper() == 'WFPC2':
            eppsky_limit = eppsky_limit_cfg * median_sky
            selfradius = 20.0  # JUST USING ACS/WFC VALUE PROBABLY NEEDS TO BE OPTIMIZED FOR WFPC2.

        if data_type.upper() == 'PC':
            eppsky_limit = eppsky_limit_cfg * median_sky
            selfradius = 20.0  # JUST USING ACS/WFC VALUE PROBABLY NEEDS TO BE OPTIMIZED FOR PC.

        # ----------------------------------------------------------
        # UVIS --> EPP > 100000. OR (EPP > 1000*sky AND EPP > 10000)
        # IR   --> EPP > 100000. OR (EPP > 100*sky AND EPP > 10000)
        # ----------------------------------------------------------

        initial_central_pixel_positions = numpy.where(numpy.logical_or(swarm_epp_listA > upper_epp_limit,
                                                                       numpy.logical_and(swarm_epp_listA > eppsky_limit,
                                                                                         swarm_epp_listA > lower_epp_limit)))[
            0]
        initial_central_pixel_list = complete_src_list[initial_central_pixel_positions, :]
        if len(initial_central_pixel_positions) == 0:
            # no bright objects
            # copy empty lists so output file is created anyway
            final_central_pixel_positions = initial_central_pixel_positions
            final_flag_src_central_pixel_list = initial_central_pixel_list
        else:
            # ---------------------------------------------------------
            # Remove duplicate central pixel position swarm candidates
            # Keep only objects that are the brightest within 20 pixels
            # ---------------------------------------------------------

            # -------------------------------------------
            # match initial central pixel list to itself
            # -------------------------------------------
            if data_type.upper() == 'IR':

                # -------------------------------------------------------
                # Define EPP cut values for filtering multiple detections
                # from a given central positions for a swarm candidate
                # -------------------------------------------------------
                cuts = [2000000., 1800000., 1000000., 500000., 70000., 20000., 0.]
                selfradii = [25., 100., 35., 30., 20., 15., 10.]

                p1 = []
                p2 = []
                for cut_cnt, cut in enumerate(cuts):

                    # --------------------------------------------------------------------
                    # Extract indices of detections that are within the set EPP cut range
                    # --------------------------------------------------------------------
                    if cut_cnt == 0:
                        cut_value_positions = numpy.where(initial_central_pixel_list[:, 3:4] > cut)[0]
                    else:
                        cut_value_positions = numpy.where(numpy.logical_and(initial_central_pixel_list[:, 3:4] >= cut,
                                                                            initial_central_pixel_list[:, 3:4] <= cuts[
                                                                                cut_cnt - 1]))[0]

                    # -----------------------------------------------
                    # If no detections exist for the specified EPP
                    # cut range, then continue to the next cut range
                    # -----------------------------------------------
                    if len(cut_value_positions) == 0:
                        continue

                    # -----------------------------------------------------------------------
                    # Determine all matches for detections in "cut_value_positions"
                    # within the radius value identified for the cut range being implemented
                    # -----------------------------------------------------------------------
                    p1_sub, p2_sub = xymatch(initial_central_pixel_list[cut_value_positions, :][:, 0:2],
                                             initial_central_pixel_list[:, 0:2], selfradii[cut_cnt], multiple=True,
                                             stack=False, verbose=False)

                    # ------------------------------------------
                    # For each cut range, add the corresponding
                    # matches to each detection to a final list
                    # ------------------------------------------
                    for p1_arr in p1_sub:
                        p1.append(p1_arr)

                    for p2_arr in p2_sub:
                        p2.append(p2_arr)

                    # Not sure if this is still needed???
                    # ------------------------------------
                    if cut_cnt == len(cuts) - 1:
                        if len(p1) == 0 and len(p2) == 0:
                            p1, p2 = xymatch(initial_central_pixel_list[:, 0:2], initial_central_pixel_list[:, 0:2],
                                             selfradius, multiple=True, stack=False, verbose=False)

                # ---------------------------------------------------------------------
                # each object is guaranteed to have at least one match (itself)
                # get brightest of each group of matches by building a list of indices
                # ---------------------------------------------------------------------
                exclude_index = None
                for i1, i2 in zip(p1, p2):
                    flux2 = initial_central_pixel_list[i2, 2]

                    # -------------------------------------------------------------
                    # Verify that there is more than one detection in a given group
                    # otherwise no detection is added to exclude index because
                    # there is only one detection for the source being evaluated
                    # -------------------------------------------------------------
                    if len(i2[numpy.where(flux2 < numpy.max(flux2))]) > 0:

                        # ----------------------------------------------------------
                        # Add all detections in grouping with a flux value less than
                        # that of the maximum flux value to an array to be excluded
                        # ----------------------------------------------------------
                        if exclude_index is None:
                            exclude_index = i2[numpy.where(flux2 < numpy.max(flux2))]
                        else:
                            exclude_index = numpy.concatenate(
                                (exclude_index, i2[numpy.where(flux2 < numpy.max(flux2))]), axis=0)

                        exclude_index = exclude_index.astype(numpy.int32)

                # -----------------------------------------------------------
                # exclude_index can have multiple copies of the same index
                # use exclude_bool array to get a list of the unique indices
                # -----------------------------------------------------------
                exclude_bool = numpy.ones(len(initial_central_pixel_list), dtype=bool)
                if not (exclude_index is None):
                    exclude_bool[exclude_index] = False
                out_values = numpy.where(exclude_bool)[0]

                # -------------------------------------------------------------------------------
                # Create final source list based on where the excluded detection indices are not
                # -------------------------------------------------------------------------------
                final_central_pixel_positions = initial_central_pixel_positions[out_values]
                final_flag_src_central_pixel_list = initial_central_pixel_list[out_values, :]

            else:

                p1, p2 = xymatch(initial_central_pixel_list[:, 0:2], initial_central_pixel_list[:, 0:2], selfradius,
                                 multiple=True, stack=False, verbose=False)

                # ---------------------------------------------------------------------
                # each object is guaranteed to have at least one match (itself)
                # get brightest of each group of matches by building a list of indices
                # ---------------------------------------------------------------------
                keep_index = numpy.arange(len(initial_central_pixel_list), dtype=int)
                for i1, i2 in zip(p1, p2):
                    flux2 = initial_central_pixel_list[i2, 2]
                    keep_index[i1] = i2[flux2.argmax()]

                # --------------------------------------------------------
                # keep_index can have multiple copies of the same index
                # use keep_bool array to get a list of the unique indices
                # --------------------------------------------------------
                keep_bool = numpy.zeros(len(initial_central_pixel_list), dtype=bool)
                keep_bool[keep_index] = True
                in_values = numpy.where(keep_bool)[0]
                final_central_pixel_positions = initial_central_pixel_positions[in_values]
                final_flag_src_central_pixel_list = initial_central_pixel_list[in_values, :]

        # ---------------------------------------------------
        # WRITE CENTRAL PIXEL POSITIONS FOR SWARMS TO A FILE
        # ---------------------------------------------------
        cetrl_pix_pos_file = phot_table_root + '_SWFILT_CENTRAL-PIX-POS.txt'
        drz_coord_out = open(cetrl_pix_pos_file, 'w')
        for i in range(len(final_flag_src_central_pixel_list)):
            drz_coord_out.write(str(final_flag_src_central_pixel_list[i, 0]) + '     ' + str(
                final_flag_src_central_pixel_list[i, 1]) + '     ' + str(
                final_flag_src_central_pixel_list[i, 2]) + '     ' + str(
                final_flag_src_central_pixel_list[i, 3]) + '     ' + str(
                final_flag_src_central_pixel_list[i, 4]) + '     ' + str(
                final_flag_src_central_pixel_list[i, 5]) + '\n')
        drz_coord_out.close()

        # ==========================================================================
        # --------------------------------------------------------------------------
        # --------------------------------------------------------------------------
        # EXTRACT THE CENTRAL PIXEL POSITIONS IN final_flag_src_central_pixel_list,
        # FROM swarm_xListB AND swarm_yListB
        # --------------------------------------------------------------------------
        # --------------------------------------------------------------------------
        # ==========================================================================

        swarm_thresh = float(param_dict["quality control"]["swarm filter"]["swarm_thresh"])
        clip_radius_list = param_dict["quality control"]["swarm filter"]["clip_radius_list"]
        clip_radius_list = list(map(float, clip_radius_list))
        scale_factor_list = param_dict["quality control"]["swarm filter"]["scale_factor_list"]
        scale_factor_list = list(map(float, scale_factor_list))
        log.info('SWARM FILTER CLIP_RADIUS_LIST: {}'.format(clip_radius_list))
        log.info('SWARM FILTER SCALE_FACTOR_LIST: {}'.format(scale_factor_list))

        # get list of objects not in the central pixel list
        keep = numpy.ones(nrows, dtype=bool)
        keep[final_central_pixel_positions] = False
        notcentral_index = numpy.where(keep)[0]
        swarm_listB = complete_src_list[notcentral_index, :]

        # views into the swarm_listB array for convenience
        swarm_xListB = swarm_listB[:, 0]
        swarm_yListB = swarm_listB[:, 1]

        # ---------------------------------------------------------------------
        # ITERATIVELY CLIP SOURCES CONTAINED WITHIN RINGS AT SPECIFIED RADIUS
        # VALUES, PROGRESSIVELY MOVING CLOSER TO THE CENTRAL SOURCE
        # ---------------------------------------------------------------------

        # do the cross-match using xymatch
        log.info('Matching {} swarm centers with {} catalog sources'.format(len(final_flag_src_central_pixel_list),
                                                                            len(swarm_listB)))
        pcentral, pfull = xymatch(final_flag_src_central_pixel_list[:, 0:2], swarm_listB[:, 0:2], clip_radius_list[0],
                                  multiple=True, stack=False, verbose=False)

        # XXX RLW: the ring list is needed only for testing, get rid of it when code works
        testing = True
        if testing:
            ring_index_list = []
            ring_refepp_list = []
            ring_thresh_list = []
            ring_count = []

        for pindex, ii in enumerate(pcentral):

            central_pixel_value = final_flag_src_central_pixel_list[ii, :]
            log.info(' ')
            log.info('CENTRAL PIXEL VALUE: {}'.format(central_pixel_value))
            log.info(' ')

            base_epp = central_pixel_value[3]
            coords = central_pixel_value[0:2]

            allmatches = pfull[pindex]

            if len(allmatches) == 0:
                # (this should not happen using xymatch)
                log.info(' ')
                log.info('------------------------------------------')
                log.info('NOTE: NO SWARM CANDIDATES FOR THIS SOURCE ')
                log.info('------------------------------------------')
                log.info(' ')
                continue

            distsq = (swarm_xListB[allmatches] - coords[0]) ** 2 + (swarm_yListB[allmatches] - coords[1]) ** 2
            sind = distsq.argsort()
            allmatches = allmatches[sind]
            distsq = distsq[sind]
            rcut = distsq.searchsorted(numpy.array(clip_radius_list) ** 2)
            for radius_cnt in range(1, len(clip_radius_list)):

                # -------------------------------------------
                # ISOLATE THE DETECTIONS WITHIN A GIVEN RING
                # -------------------------------------------

                matches = allmatches[rcut[radius_cnt]:rcut[radius_cnt - 1]]

                if len(matches) == 0:
                    log.info(' ')
                    log.info('------------------------------------------')
                    log.info('NOTE: ALL MATCHES/DETECTIONS IN THIS RING ')
                    log.info('      HAVE PREVIOUSLY BEEN ACCOUNTED FOR  ')
                    log.info('------------------------------------------')
                    log.info(' ')

                    continue

                # -----------------------------------------------------------
                # CALCULATE THE MEDIAN SKY VALUE FOR THE GROUP OF DETECTIONS
                # CONTAINED WITHIN THE SPECIFIED RING BEING PROCESSED
                # -----------------------------------------------------------
                ref_epp = base_epp * scale_factor_list[radius_cnt - 1]

                # -----------------------------------------------------------------------------------
                # DIFFERENTIATE BETWEEN GOOD DETECTIONS AND SWARM DETECTIONS WITHIN SPECIFIED RINGS
                # -----------------------------------------------------------------------------------
                ring = swarm_listB[matches, :]
                w = numpy.where(ring[:, 3] / ref_epp < swarm_thresh)
                if len(w) > 0:
                    swarm_flag[notcentral_index[matches[w]]] = True

                # XXX RLW: following needed only for testing, get rid of it when code works
                if testing:
                    ring_index_list.append(matches)
                    ring_count.append(len(matches))
                    ring_refepp_list.append(ring[:, 3] / ref_epp)
                    ring_thresh_list.append(swarm_thresh)

        # XXX RLW: following needed only for testing, get rid of it when code works
        if testing:
            # -----------------------------------------------------------------------------------------
            # WRITE CLIPPED SOURCES CONTAINED WITHIN RINGS TO AN OUTPUT FILE FOR INTERMEDIATE ANALYSIS
            # -----------------------------------------------------------------------------------------
            ring_source_file = phot_table_root + '_SWFILT_RING-SOURCE-INFO.txt'
            ring_src_outfile = open(ring_source_file, 'w')
            ring_src_outfile.write(
                "# ------------------------------------------------------------------------------------------------\n")
            ring_src_outfile.write(
                "# X-Center   Y-Center     Flux        ElectronPP          Sky        SrcEPP/RefEPP   Swarm Thresh \n")
            ring_src_outfile.write(
                "# ------------------------------------------------------------------------------------------------\n")

            if ring_index_list:
                ring_index_list = numpy.concatenate(ring_index_list)

                # select just the lowest value of refepp/swarm threshold for each source
                # create array with extra columns
                ring_source_list = numpy.empty((len(ring_index_list), 9), dtype=numpy.float)
                ring_source_list[:, 0:6] = swarm_listB[ring_index_list, :]
                ring_source_list[:, 6] = numpy.concatenate(ring_refepp_list)
                ring_source_list[:, 7] = numpy.repeat(ring_thresh_list, ring_count)
                ring_source_list[:, 8] = ring_source_list[:, 6] / ring_source_list[:, 7]

                # sort by x, y, and refepp
                # tricky here: get a view with named columns, then specify names as sort items
                ring_source_list.view(','.join(['f8'] * 9)).sort(order=['f0', 'f1', 'f8'], axis=0)

                # keep just first entry when the same source appears more than once
                keep = numpy.ones(len(ring_index_list), dtype=bool)
                # keep[1:] = numpy.any(ring_source_list[1:,0:2]!=ring_source_list[:-1,0:2], axis=1)
                keep[1:] = numpy.logical_or(ring_source_list[1:, 0] != ring_source_list[:-1, 0],
                                            ring_source_list[1:, 1] != ring_source_list[:-1, 1])
                ring_source_list = ring_source_list[keep, :]

                #                numpy.savetxt(ring_src_outfile,ring_source_list,delimiter='     ')

                for ring_source in ring_source_list:
                    ring_src_outfile.write(str(ring_source[0]) + '     ' + str(ring_source[1]) + '     ' + str(
                        ring_source[2]) + '     ' + str(ring_source[3]) + '     ' + str(ring_source[4]) + '     ' + str(
                        ring_source[5]) + '     ' + str(ring_source[6]) + '     ' + str(ring_source[7]) + '\n')

            ring_src_outfile.close()  # XXX RLW: end of testing code

        # ===================================================================================
        # -----------------------------------------------------------------------------------
        # -----------------------------------------------------------------------------------
        #                          ----- PROXIMITY FILTER -----
        # EXTRACT ADDITIONAL SWARM DETECTIONS BASED ON THE SWARM CANDIDATE CENTRAL POSITIONS,
        # DEFINING THE REMOVAL RADIUS AROUND EACH POSITION BASED ON THAT SOURCE'S EPP
        # -----------------------------------------------------------------------------------
        # -----------------------------------------------------------------------------------
        # ===================================================================================

        proximity_flag = numpy.zeros(nrows, dtype=bool)

        proximity_choice = param_dict["quality control"]["swarm filter"]["proximity_binary"]

        if proximity_choice == 'yes':
            if len(final_flag_src_central_pixel_list) > 0:
                # XXX these ought to come from config files
                if data_type == 'ir':
                    # ctrList_radiusList = [25,100,80,30,50,20,15,10]
                    # ctrList_thresholdList = [2000000,1800000,500000,250000,100000,40000,20000,2000]
                    ctrList_radiusList = [125, 100, 80, 30, 50, 20, 15]
                    ctrList_thresholdList = [2000000, 1800000, 500000, 250000, 100000, 40000, 20000]
                if data_type == 'uvis':
                    ctrList_radiusList = [40, 35, 20, 15, 10]
                    ctrList_thresholdList = [100000, 70000, 50000, 10000, 2000]
                if data_type == 'wfc':  # acs/wfc
                    ctrList_radiusList = [40, 35, 20, 15,
                                          10]  # just copied UVIS values. These values need to be optimized for ACS/WFC
                    ctrList_thresholdList = [100000, 70000, 50000, 10000,
                                             2000]  # just copied UVIS values. These values need to be optimized for ACS/WFC
                if data_type == 'hrc':
                    ctrList_radiusList = [40, 35, 20, 15,
                                          10]  # just copied UVIS values. These values need to be optimized for ACS/HRC
                    ctrList_thresholdList = [100000, 70000, 50000, 10000,
                                             2000]  # just copied UVIS values. These values need to be optimized for ACS/HRC
                if data_type == 'wfpc2':  # just copied UVIS values. These values need to be optimized for WFPC2
                    ctrList_radiusList = [40, 35, 20, 15, 10]
                    ctrList_thresholdList = [100000, 70000, 50000, 10000, 2000]
                if data_type == 'pc':  # just copied UVIS values. These values need to be optimized for PC
                    ctrList_radiusList = [40, 35, 20, 15, 10]
                    ctrList_thresholdList = [100000, 70000, 50000, 10000, 2000]

                for ctrList_cnt, (threshold, radius) in enumerate(zip(ctrList_thresholdList, ctrList_radiusList)):

                    if ctrList_cnt == 0:
                        ctr_list_cut = final_flag_src_central_pixel_list[:, 3] > threshold
                    else:
                        ctr_list_cut = numpy.logical_and(final_flag_src_central_pixel_list[:, 3] > threshold,
                                                         final_flag_src_central_pixel_list[:, 3] <=
                                                         ctrList_thresholdList[ctrList_cnt - 1])

                    ctr_list_cut1 = final_flag_src_central_pixel_list[ctr_list_cut, :]
                    pcentral, pfull = xymatch(ctr_list_cut1[:, 0:2], swarm_listB[:, 0:2], radius, multiple=True,
                                              verbose=False)
                    proximity_flag[notcentral_index[pfull]] = True

            log.info("Proximity filter flagged {} sources".format(proximity_flag.sum()))

            # --------------------------------------------------------------------------
            # WRITE NEAR CENTRAL POSITION SWARM LIST TO AN OUTPUT FILE FOR VERIFICATION
            # --------------------------------------------------------------------------
            near_swmList = complete_src_list[proximity_flag, :]
            final_near_swarm_file = open(phot_table_root + '_SWFILT_NEAR_SWARM_FILE.txt', 'w')
            for swarm_value in near_swmList:
                final_near_swarm_file.write(
                    str(swarm_value[0]) + '     ' + str(swarm_value[1]) + '     ' + str(swarm_value[2]) + '     ' + str(
                        swarm_value[3]) + '     ' + str(swarm_value[4]) + '\n')
            final_near_swarm_file.close()

        # -------------------------------------------------------------------------
        # EXTRACT DETECTIONS FROM THE complete_src_list THAT ARE NOT FLAGGED
        # -------------------------------------------------------------------------

        combined_flag = numpy.logical_or(swarm_flag, proximity_flag)
        final_swarm_list = complete_src_list[combined_flag, :]
        final_source_list = complete_src_list[numpy.logical_not(combined_flag), :]

        log.info(' ')
        log.info('************************************************')
        log.info('INITIAL LENGTH OF complete_src_list = {}'.format(len(complete_src_list)))
        log.info(' ')
        log.info('LENGTH OF final_source_list = {}'.format(len(final_source_list)))
        log.info('LENGTH OF final_swarm_list = {}'.format(len(final_swarm_list)))
        log.info('TOTAL LENGTH = {}'.format(len(final_source_list) + len(final_swarm_list)))
        log.info(' ')
        log.info('MEDIAN SKY VALUE = {}'.format(median_sky))
        log.info('************************************************')
        log.info(' ')

        # ----------------------------------------------------
        # WRITE SWARM LIST TO AN OUTPUT FILE FOR VERIFICATION
        # ----------------------------------------------------
        final_swarm_file = open(phot_table_root + '_SWFILT_SWARM_FILE.txt', 'w')
        for swarm_value in final_swarm_list:
            final_swarm_file.write(
                str(swarm_value[0]) + '     ' + str(swarm_value[1]) + '     ' + str(swarm_value[2]) + '     ' + str(
                    swarm_value[3]) + '     ' + str(swarm_value[4]) + '\n')
        final_swarm_file.close()

        # ----------------------------------------------------
        # WRITE SOURCE LIST TO AN OUTPUT FILE FOR VERIFICATION
        # ----------------------------------------------------
        final_source_file = open(phot_table_root + '_SWFILT_SOURCE_FILE.txt', 'w')
        for source_value in final_source_list:
            final_source_file.write(
                str(source_value[0]) + '     ' + str(source_value[1]) + '     ' + str(source_value[2]) + '     ' + str(
                    source_value[3]) + '     ' + str(source_value[4]) + '\n')
        final_source_file.close()

        # =================================================================
        # -----------------------------------------------------------------
        # -----------------------------------------------------------------
        # WRITE SWARM FLAGS TO OUTPUT PHOT TABLE BASED ON final_swarm_list
        # -----------------------------------------------------------------
        # -----------------------------------------------------------------
        # =================================================================
        phot_table_temp = phot_table_root + '_SWFILT.txt'
        phot_table_out = open(phot_table_temp, 'w')

        phot_table_in = open(phot_table, 'r')
        phot_table_rows = phot_table_in.readlines()
        phot_table_in.close()

        phot_table_out.write(phot_table_rows[0])
        for i, table_row in enumerate(phot_table_rows[1:]):
            if combined_flag[i]:
                row_split = table_row.split(',')
                sat_flag = int(row_split[-1]) | 32
                row_split[-1] = str(sat_flag) + '\n'
                table_row = ','.join(row_split)
            phot_table_out.write(table_row)

        phot_table_out.close()

        os.system('mv ' + phot_table + ' ' + phot_table + '.PreSwarmFilt')
        os.system('mv ' + phot_table_temp + ' ' + phot_table)

        log.info(' ')
        log.info('FINAL SWAR-FILT PHOT_TABLE: {}'.format(phot_table))
        log.info(' ')

# ======================================================================================================================

def HLANexpFlags_OLD(all_drizzled_filelist, working_hla_red, filter_sorted_flt_dict, param_dict,
                     readnoise_dictionary_drzs,scale_dict_drzs, exp_dictionary_scis, dict_newTAB_matched2drz,
                     drz_root_dir):
    """flags out sources from regions where there are a low (or a null) number of contributing exposures

    all_drizzled_filelist : list
        List of drizzled images to process

    working_hla_red : string
        ***UNUSED*** full path to working directory.

    filter_sorted_flt_dict : dictionary
        dictionary containing lists of calibrated images sorted (also keyed) by filter name.

    param_dict : dictionary
        Dictionary of instrument/detector - specific drizzle, source finding and photometric parameters

    readnoise_dictionary_drzs : dictionary
        ***UNUSED*** dictionary of readnoise values keyed by drizzled image.

    scale_dict_drzs : dictionary
        ***UNUSED*** dictionary of scale values keyed by drizzled image.

    exp_dictionary_scis :  dictionary
        ***UNUSED*** dictionary of exposure time values keyed by drizzled image.

    dict_newTAB_matched2drz : dictionary
        dictionary of source lists keyed by drizzled image name.

    drz_root_dir : dictionary of source lists keyed by drizzled image name.

    Returns
    -------
    nothing!

    """
    # ------------------
    # CREATE NEXP IMAGE
    # ------------------
    for drizzled_image in all_drizzled_filelist:
        image_split = drizzled_image.split('/')[-1]
        channel = drizzled_image.split("_")[
            -3].upper()

        # if channel == 'IR':
        #    continue

        # ---------
        # METHOD 1:
        # ---------
        ctx = getdata(drizzled_image, 3)

        if channel in ['UVIS', 'IR', 'WFC', 'HRC']: ncombine = getheader(drizzled_image, 1)['NCOMBINE']
        if channel in ['WFPC2', 'PC']:
            ndrizim = getheader(drizzled_image, 0)['NDRIZIM']
            ncombine = ndrizim / 4
        if channel == 'SBC':
            ndrizim = getheader(drizzled_image, 0)['NDRIZIM']
            ncombine = ndrizim
        ctxarray = arrayfy_ctx(ctx, ncombine)

        nexp_array_ctx = ctxarray.sum(axis=-1)
        nexp_image_ctx = drizzled_image.split('.')[0] + '_NCTX.fits'
        if not os.path.isfile(nexp_image_ctx):
            hdr = getheader(drizzled_image, 1)
            pyfits.writeto(nexp_image_ctx, numpy.float32(nexp_array_ctx), hdr)

        # ---------
        # METHOD 2:
        # ---------
        drz_data = getdata(drizzled_image, 1)

        ## this bit is added to get the mask integrated into the exp map
        maskfile = drizzled_image.replace('_drz.fits', '_msk.fits')
        if os.path.isfile(maskfile):
            mask_data = getdata(maskfile)
            mask_array = (mask_data == 0.0).astype(numpy.int32)

        component_drz_img_list = get_component_drz_list(drizzled_image, drz_root_dir, filter_sorted_flt_dict)

        nx = drz_data.shape[0]
        ny = drz_data.shape[1]
        nexp_array = numpy.zeros((nx, ny), dtype=numpy.int32)

        for comp_drz_img in component_drz_img_list:
            comp_drz_data = (getdata(comp_drz_img) != 0).astype(numpy.int32)
            try:
                nexp_array += comp_drz_data
            except ValueError:
                log.info("WARNING: Astrodrizzle added an extra-row/column...")
                nexp_array += comp_drz_data[0:nx, 0:ny]

        if os.path.isfile(maskfile):
            nexp_array = nexp_array * mask_array
        else:
            log.info("something's wrong: maskfile {} is not a file".format(maskfile))
            sys.exit()
        nexp_image = drizzled_image.split('.')[0] + '_NEXP.fits'
        if not os.path.isfile(nexp_image):
            hdr = getheader(drizzled_image, 1)
            pyfits.writeto(nexp_image, numpy.float32(nexp_array), hdr)

        # -------------------------------------------------------
        # EXTRACT FLUX/NEXP INFORMATION FROM NEXP IMAGE BASED ON
        # THE SOURCE DETECTION POSITIONS PREVIOUSLY ESTABLISHED
        # -------------------------------------------------------
        phot_table = dict_newTAB_matched2drz[drizzled_image]
        phot_table_root = phot_table.split('/')[-1].split('.')[0]

        phot_table_in = open(phot_table, 'r')
        phot_table_lines = phot_table_in.readlines()
        phot_table_in.close()

        nrows = len(phot_table_lines) - 1
        cat_coords = numpy.empty((nrows, 2), dtype=float)
        for line_cnt, phot_table_line in enumerate(phot_table_lines):

            if line_cnt == 0:
                continue

            phot_table_line_split = phot_table_line.split(',')
            x_coord = float(phot_table_line_split[0])
            y_coord = float(phot_table_line_split[1])
            cat_coords[line_cnt - 1, :] = [x_coord, y_coord]

        # ----------------------------------
        # Convert aperture radius to pixels
        # ----------------------------------

        ap2 = param_dict['catalog generation']['aperture_2']

        if channel == 'IR':
            radius = ap2 / 0.09
        if channel == 'UVIS':
            radius = ap2 / 0.04
        if channel == 'WFC':  # ACS/WFC
            radius = ap2 / 0.05
        if channel == 'HRC':
            radius = ap2 / 0.025
        if channel == 'SBC':
            radius = ap2 / 0.03
        if channel == 'WFPC2':
            radius = ap2 / 0.1
        if channel == 'PC':
            radius = ap2 / 0.05

        num_exp = round(numpy.max(nexp_array))
        if num_exp <= 1 or channel in ('IR', 'SBC'):
            # Keep everything that has an exposure for detectors without CRs or
            # when there is only one exposure
            artifact_filt = 0.5
        elif num_exp > 5:
            # Flag sources with <= 2 exposures when there are > 5 total
            # We are always using the 'imedian' combination in that case, and it
            # does not do a very good job of rejecting CRs with only 2 available
            # exposures
            artifact_filt = 2.5
        else:
            artifact_filt = 1.5

        icoords = (cat_coords + 0.5).astype(int)
        # note x & y are swapped so they can index the numpy array nexp_array
        # catalog x is second subscript, catalog y is first subscript
        ix = icoords[:, 1]
        iy = icoords[:, 0]

        # get list of neighboring pixels that are within radius
        iradius = int(radius + 1)
        idiam = iradius * 2 + 1
        gx, gy = numpy.mgrid[0:idiam, 0:idiam] - iradius
        gx = gx.ravel()
        gy = gy.ravel()
        w = numpy.where(gx ** 2 + gy ** 2 <= radius ** 2)[0]
        gx = gx[w]
        gy = gy[w]

        # check the pixel values for low nexp

        # this version uses numpy broadcasting sum gx+ix is [len(gx),nrows]
        gx = (gx[:, numpy.newaxis] + ix).clip(0, nexp_array.shape[0] - 1)
        gy = (gy[:, numpy.newaxis] + iy).clip(0, nexp_array.shape[1] - 1)
        artifact_flag = nexp_array[gx, gy].min(axis=0) < artifact_filt

        log.info('FLAGGING {} OF {} SOURCES'.format(artifact_flag.sum(), nrows))

        # -------------------------------------------------------------------
        # WRITE NEXP FLAGS TO OUTPUT PHOT TABLE BASED ON nexp_phot_data_list
        # -------------------------------------------------------------------
        # nexp_outfile_good = open('temp_outfile_NEXP_GOOD.txt','w')
        # nexp_outfile_bad = open('temp_outfile_NEXP_BAD.txt','w')

        phot_table_temp = phot_table_root + '_NEXPFILT.txt'
        phot_table_out = open(phot_table_temp, 'w')

        phot_table_in = open(phot_table, 'r')
        phot_table_rows = phot_table_in.readlines()
        phot_table_in.close()

        phot_table_out.write(phot_table_rows[0])
        for i, table_row in enumerate(phot_table_rows[1:]):
            if artifact_flag[i]:
                row_split = table_row.split(',')
                nexp_flag = int(row_split[-1]) | 64
                row_split[-1] = str(nexp_flag) + '\n'
                table_row = ','.join(row_split)
            phot_table_out.write(table_row)

        phot_table_out.close()
        # nexp_outfile_good.close()
        # nexp_outfile_bad.close()

        os.system('mv ' + phot_table + ' ' + phot_table + '.PreNexpFilt')
        os.system('mv ' + phot_table_temp + ' ' + phot_table)

        log.info('Created new version of {}'.format(phot_table))

# ======================================================================================================================

def HLA_flag4and8_hunter_killer_OLD(photfilename):
    """This function searches through photometry catalogs for sources whose flags contain
    both bits 4 (multi-pixel saturation), and 8 (faint magnitude limit).
    If found, the subroutine removes the "8" bit value from the set of flags for that source.

    Parameters
    ----------
    photfilename : string
        name of sourcelist to process

    Returns
    -------
    nothing!
    """

    # for flag_value in

    inf=open(photfilename)
    phot_lines=inf.readlines()
    inf.close()
    fout=open(photfilename,'w')
    conf_ctr=0
    log.info("Searching {} for flag 4 + flag 8 conflicts....".format(photfilename))
    for phot_line in phot_lines:
        phot_line=phot_line.strip()
        parse_pl=phot_line.split(',')
        x=parse_pl[0]
        if (x[0] == 'X'):out_line=phot_line
        else:
            flagval=int(parse_pl[-1])
            if ((flagval & 4 >0) and (flagval & 8 >0)):
                conf_ctr+=1
                parse_pl[-1]=str(int(parse_pl[-1])-8)
            out_line = ""
            for item in parse_pl:out_line=out_line+"{},".format(item)
            out_line=out_line[:-1]
        fout.write("%s\n"%(out_line))
    fout.close()
    if conf_ctr == 0: log.info("No conflicts found.")
    if conf_ctr == 1: log.info("{} conflict fixed.".format(conf_ctr))
    if conf_ctr > 1:  log.info("{} conflicts fixed.".format(conf_ctr))

# +++++++++++++++++++++++++++++++++++++++++ END OLD VERSIONS +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def get_component_drz_list(drizzled_image, drz_root_dir, flt_file_names):

    """Get a list of the drizzled exposure images associated with this combined drizzled image

    Usually this can just use glob to get a list of all the drizzled exposures for this
    filter, but it also handles the case where some exposures were not used (e.g., for
    scan mode images).

    drizzled_image : string
        Name of combined (level 2) drizzled image

    drz_root_dir : string
        Location of drizzled exposures

    filter_sorted_flt_dict : dictionary
        dictionary containing lists of calibrated images sorted (also keyed) by filter name.

    Returns
    -------
    rv : list
        a list of drizzled exposure images associated with the specified combined drizzled image
    """

    drz_img_split = drizzled_image.split('/')[-1].split('drz')
    component_drz_img_list = glob.glob(os.path.join(drz_root_dir,drz_img_split[0])+'*_drz.fits')
    component_drz_img_list.sort()

    drz_filter = drizzled_image.split("_")[5]  # TODO: REFACTOR FOR HAP. this is just a short-term hack to get things working for HLA

    if type(flt_file_names).__name__ == 'dict':
        list_of_flts = flt_file_names[drz_filter.lower()]
    if type(flt_file_names).__name__ == 'list':
        list_of_flts = flt_file_names

    if len(list_of_flts) == len(component_drz_img_list):
        # length match means we use them all
        return component_drz_img_list
    elif len(list_of_flts) > len(component_drz_img_list):
        # this must be a bug?
        log.info("ERROR: too few drizzled exposures for {}".format(drz_filter))
        log.info("Drizzled exposure list: {}".format("\n".join(component_drz_img_list)))
        log.info("flt exposure list: {}".format("\n".join(list_of_flts)))
        log.info("Plowing ahead with the full drizzled list")
        return component_drz_img_list
    # check the drz headers to see which ipppssoots are included
    ipdict = {}
    for ipname in list_of_flts:
        fname = os.path.split(ipname)[-1]
        fname = fname.split('_')[0].lower()
        ipdict[fname] = 1
    rv = []
    for drzfile in component_drz_img_list:
        fh = pyfits.open(drzfile)
        rootname = fh[0].header.get('rootname','')
        fh.close()
        fname = os.path.split(rootname)[-1]
        fname = fname.split('_')[0].lower()
        if fname in ipdict:
            rv.append(drzfile)
    if len(list_of_flts) != len(rv):
        # this must be a bug?
        log.info("ERROR: mismatch after filtering in exposure lists for {}".format(drz_filter))
        log.info("Filtered drizzled exposure list: {}".format("\n".join(rv)))
        log.info("flt exposure list: {}".format("\n".join(list_of_flts_in_main_driz)))
        log.info("Plowing ahead with the filtered drizzled list")
    return rv


# =============================================================================
# -----------------------------------------------------------------------------
# ------------------------ AUXILLIARY FUNCTIONS BELOW -------------------------
# -----------------------------------------------------------------------------
# =============================================================================


def xymatch(cat1, cat2, sep, multiple=False, stack=True, verbose=True):
    """Routine to match two lists of objects by position using 2-D Cartesian distances.

    Matches positions in cat1 with positions in cat2, for matches within separation (sep).
    If more than one match is found, the nearest is returned.
    Setting multiple=True returns all matching pairs rather than just the closest.

    Input catalogs need not be sorted. They should be 2-element arrays [:,2] with
    cat1[:,0] = x1 and cat1[:,1] = y1.

    Returns an array of indices for cat2 that correspond to the closest match
    (within sep) in cat2 for each object in cat1, so x1[i] matches x2[return_value[i]].
    Note that objects in cat1 with no match within sep have indices
    -N-1 with N the length of cat2, so that IndexErrors will be
    raised if trying to assign these indices.

    If multiple is true, returns a tuple (p1,p2) such that cat1[p1]
    and cat2[p2] are within sep.  p1 and p2 may include multiple pointers to
    the same objects in cat1 or cat2.  In this case objects that don't match
    are simply omitted from the lists.

    The stack parameter applies only when multiple is true.
    If stack is true (the default), the returned matching pointers are stacked
    into a single big array, so both p1 and p2 are 1-D arrays of
    length nmatches.

    If stack is false, p1 is a list of indices into cat1, and p2 is
    a list of *array* indices into cat2.  So p2[k] is all the sources
    that match p1[k].  This version maybe be more useful if you need
    to look at the groups of sources associated with cat1 objects.
    Only cat1 objects that have a match are included in these lists.

    Set verbose to False to run without status info.

    Marcel Haas, 2012-06-29, after IDL routine xymatch.pro by Rick White
    With some tweaks by Rick White
    
    Parameters
    ----------
    cat1 : numpy.ndarray
        list of x,y source coords to match.
    
    cat2 : numpy.ndarray
        list of x,y source coords to match.
    
    sep : float
        maximum separation (in pixels) allowed for source matching.
    
    multiple : Boolean
        If multiple is true, returns a tuple (p1,p2) such that cat1[p1] and cat2[p2] are within sep.  p1 and p2 may include multiple pointers to the same objects in cat1 or cat2.  In this case objects that don't match are simply omitted from the lists. Default value is 'False'.
    
    stack : Boolean
        If stack is true, the returned matching pointers are stacked into a single big array, so both p1 and p2 are 1-D arrays of length nmatches. Default value is 'True'.
    
    verbose : Boolean
        print verbose output? Default value is 'True'.

    Returns
    -------
    Varies; Depending on inputs, either just 'p2', or 'p1' and 'p2'. p1 and p2 are lists of matched indicies
    """
    if not (isinstance(cat1, numpy.ndarray) and len(cat1.shape)==2 and cat1.shape[1]==2):
        raise ValueError("cat1 must be a [N,2] array")
    if not (isinstance(cat2, numpy.ndarray) and len(cat2.shape)==2 and cat2.shape[1]==2):
        raise ValueError("cat2 must be a [N,2] array")

    x1 = cat1[:,0]
    y1 = cat1[:,1]
    x2 = cat2[:,0]
    y2 = cat2[:,1]

    # Sort the arrays by increasing y-coordinate
    is1 = y1.argsort()
    x1 = x1[is1]
    y1 = y1[is1]
    is2 = y2.argsort()
    x2 = x2[is2]
    y2 = y2[is2]
    # find search limits in y2 for each object in y1
    kvlo = y2.searchsorted(y1-sep,'left').clip(0, len(y2))
    kvhi = y2.searchsorted(y1+sep,'right').clip(kvlo, len(y2))

    nnomatch = 0
    n1 = len(x1)
    if multiple:
        # build lists of array segments for matches
        p1 = []
        p2 = []
    else:
        p2 = numpy.zeros(n1, dtype='int') - len(x2) - 1
    t0 = time.time()
    sepsq = sep**2
    for i in range(n1):
        y = y1[i]
        x = x1[i]
        klo = kvlo[i]
        khi = kvhi[i]
        dx = numpy.abs(x2[klo:khi] - x)
        w = (dx <= sep).nonzero()[0]
        if len(w) == 0:
            # Nothing matched
            nnomatch += 1
        else:
            distsq = (x - x2[klo+w])**2 + (y - y2[klo+w])**2

            if multiple:
                ww = (distsq <= sepsq).nonzero()[0]
                if len(ww) == 0:
                    nnomatch += 1
                else:
                    if stack:
                        p1.append(numpy.zeros(len(ww),dtype='int')+is1[i])
                    else:
                        p1.append(is1[i])
                    p2.append(is2[klo + w[ww]])

            else:
                if dist.min() <= sep:
                    p2[is1[i]] = is2[klo+w[dist.argmin()]]
                else:
                    nnomatch += 1

        if verbose and (i+1) % 10000 == 0:
            log.info("%.1f s: Finished %d of %d (%d unmatched)" % (time.time()-t0, i+1, n1, nnomatch))

    if verbose:
        log.info("%.1f s: Finished %d (%d unmatched)" % (time.time()-t0, n1, nnomatch))
    if multiple:
        if stack:
            if len(p1) == 0:
                # no matches found
                # return empty integer arrays that are still usable as indices
                return numpy.array([],dtype=int), numpy.array([],dtype=int)
            else:
                return numpy.concatenate(p1),numpy.concatenate(p2)
        else:
            return (p1, p2)
    else:
        return p2

# ======================================================================================================================

def sorted_median(a):
    """Compute the median for a 1-D numpy array that is already sorted

    Parameters
    ----------
    a : numpy.ndarray
        list of values what will be processed.
     
    Returns
    -------
    med : numpy.float32
        median value
    """

    ll = len(a)
    if (ll % 2) == 1:
        med = a[ll//2]
    else:
        med = 0.5*(a[ll//2]+a[ll//2 - 1])
    return med

# ======================================================================================================================

def get_median_sky(drizzled_image):
    """Read drizzled image from FITS file and compute sky
    
    Parameters
    ----------
    drizzled_image : string
        name of image that will be used to compute median sky value.

    Returns
    -------
    median_sky : float
        median sky value for image specified by 'drizzled_image'. 
    """

    driz_img_data = getdata(drizzled_image,1).ravel()
    # change the datatype to swap pixels if necessary
    # this makes searchsorted much faster
    dstring = str(driz_img_data.dtype)
    if dstring[0] in '<>':
        driz_img_data = driz_img_data.astype(dstring[1:])
    # ignore zero pixels (missing data)
    driz_img_data = driz_img_data[driz_img_data != 0]
    if len(driz_img_data) == 0:
        # all zero pixels so sky is zero too
        median_sky = 0.0
    else:
        driz_img_data.sort()
        for clip_cnt in range(5):
            driz_data_std = numpy.std(driz_img_data)
            # Computing median (easy since array is sorted)
            driz_data_med = sorted_median(driz_img_data)
            upper_limit = driz_data_med + (4. * driz_data_std)
            lower_limit = driz_data_med - (4. * driz_data_std)
            subIndexUpper = driz_img_data.searchsorted(upper_limit,'right')
            subIndexLower = driz_img_data.searchsorted(lower_limit,'left')
            driz_img_data = driz_img_data[subIndexLower:subIndexUpper]
        median_sky = sorted_median(driz_img_data)
        if median_sky == 0.0:
            median_sky = numpy.mean(driz_img_data)
    del driz_img_data # free memory
    return median_sky

# ======================================================================================================================

def rdtoxy(rd_coord_array, image, image_ext):
    """converts RA and dec to x,y image coords.

    rd_coord_array : numpy.ndarray
        array containing RA and dec values to convert.
    
    image : string
        drizzled image whose WCS info will be used in the coordinate conversion. 
    
    image_ext : string
        fits image extension to be used in the conversion.
    
    Returns
    xy_arr: array
        array of converted x,y coordinate value pairs
    """

    scifile = image + image_ext
    wcs = wcsutil.HSTWCS(scifile)
    try:
        xy_arr = wcs.wcs_sky2pix(rd_coord_array,1)
    except AttributeError:
        xy_arr = wcs.wcs_world2pix(rd_coord_array,1)
    return (xy_arr)

# ======================================================================================================================

def xytord(xy_coord_array, image, image_ext):
    """converts x,y image coords to RA and dec.

    xy_coord_array : numpy.ndarray
        array containing image x,y coord values to convert.

    image : string
        drizzled image whose WCS info will be used in the coordinate conversion.

    image_ext : string
        fits image extension to be used in the conversion.

    Returns
    -------
    rd_arr : array
        an array of converted RA and dec value pairs
    """

    scifile = image + image_ext
    wcs = wcsutil.HSTWCS(scifile)
    try:
        rd_arr = wcs.all_pix2sky(xy_coord_array,1)
    except AttributeError:
        rd_arr = wcs.all_pix2world(xy_coord_array,1)
    return (rd_arr)

# ========================================================================================================
# ========================================= NEW FUNCTIONS 131002 =========================================
# ========================================================================================================
def extract_name(stringWpath):
    """This task will extract just the name of  specific filename that includes the path in the name: 'stringWpath'.
    
    Tested.
    
    stringWpath : string
        input path to be processed

    Returns
    -------
    stringname : string
        name of string
    """
    while "/" == stringWpath[-1]:
        stringWpath = stringWpath[:-1]
    stringname = stringWpath.split("/")[-1]
    return stringname

# ======================================================================================================================
def arrayfy_ctx(ctx, maxindex):

    """Function to turn the context array returned by AstroDrizzle
    into a bit datacube with the third dimension equal to maxindex.
    Requires care since vanilla python seems to be a bit loose with
    data types, while arrays in numpy are strictly typed.
    Comments indicate the why of certain operations.
    Current version requires maxindex to be specified.  In upgrade, could
    estimate maxindex from the highest non-zero bit set in ctx.
    (In principle maxindex can be greater, if the last images contribute
    to no pixels, or smaller, if the calling routine chooses to ignore
    the contribution of further images.)

    Per AstroDrizzle specifications, ctx can be a 2-dimensional or 3-dimensional
    array of 32-bit integers.  It will be 2-dimensional if fewer than 32 images
    are combined; 3-dimensional if more images are combined, in which case the
    ctx[:,:,0] contains the bit values for images 1-32, ctx[:,:,1] for images
    33-64, and so forth.

    Parameters
    ----------
    ctx : numpy.ndarray
        input context array to be converted

    maxindex : int
        maximum index value to process

    Returns
    -------
    ctxarray : numpy.ndarray
        ctxarray, The input context array converted to datacube form.
    """
    nx = ctx.shape[0]
    ny = ctx.shape[1]
    nz = 1
    if ctx.ndim > 2:
        nz = ctx.shape[2]
    n3 = maxindex
    # Need to find out how to handle the case in which maxindex is not specified

    ctxarray = numpy.zeros ( (nx, ny, n3), dtype="bool" )
    comparison = numpy.zeros ( (nx, ny), dtype="int32")
    
    for i in range(n3):
        ilayer = int(i/32)
        ibit = i - 32*ilayer
        cc = comparison + 2**ibit
        if nz > 1:
            ctxarray [:,:,i] = numpy.bitwise_and (ctx[:,:,ilayer], cc)
        else:
            ctxarray [:,:,i] = numpy.bitwise_and (ctx[:,:], cc)
    return ctxarray

# ======================================================================================================================

def HLA_flag4and8_hunter_killer(catalog_data):
    """This function searches through photometry catalogs for sources whose flags contain
    both bits 4 (multi-pixel saturation), and 8 (faint magnitude limit).
    If found, the subroutine removes the "8" bit value from the set of flags for that source.

    Parameters
    ----------
    catalog_data : astropy Table object
        catalog data to process

    Returns
    -------
    catalog_data : astropy Table object
        input catalog data with updated flags
    """
    conf_ctr=0
    log.info("Searching for flag 4 + flag 8 conflicts....")
    if 'FLAGS' in catalog_data.keys():
        flag_col_title = 'FLAGS'
    elif 'Flags' in catalog_data.keys():
        flag_col_title = 'Flags'
    else:
        sys.exit("ERROR! Unrecognized catalog format!")
    for catalog_line in catalog_data:
        if ((catalog_line[flag_col_title] & 4 >0) and (catalog_line[flag_col_title] & 8 >0)):
            conf_ctr+=1
            catalog_line[flag_col_title]=int(catalog_line[flag_col_title])-8
    if conf_ctr == 0: log.info("No conflicts found.")
    if conf_ctr == 1: log.info("{} conflict fixed.".format(conf_ctr))
    if conf_ctr > 1:  log.info("{} conflicts fixed.".format(conf_ctr))

    return catalog_data