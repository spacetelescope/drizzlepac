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
import glob
import json
import math
import os
import sys
import time

from astropy.io import fits as fits
from astropy.table import Table
import numpy
import scipy
import scipy.ndimage

from drizzlepac.hlautils import ci_table
from stsci.tools import logutil
from stwcs import wcsutil


__taskname__ = 'hla_flag_filter'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)


def run_source_list_flagging(drizzled_image, flt_list, param_dict, exptime, plate_scale, median_sky,
                             catalog_name, catalog_data, proc_type, drz_root_dir, hla_flag_msk, ci_lookup_file_path,
                             output_custom_pars_file, log_level, diagnostic_mode):

    """Simple calling subroutine that executes the other flagging subroutines.

    Parameters
    ----------
    drizzled_image : string
        drizzled filter product image filename

    flt_list : list
        list of calibrated images that were drizzle-combined to produce image specified by input parameter
        'drizzled_image'

    param_dict : dictionary
        Dictionary of instrument/detector - specific drizzle, source finding and photometric parameters

    exptime : float
        drizzled filter product exposure time in seconds

    plate_scale : float
        plate scale, in arcseconds/pixel

    median_sky : float
        median sky value

    catalog_name : string
        drizzled filter product catalog filename

    catalog_data : astropy.Table object
        drizzled filter product catalog data

    proc_type : string
        sourcelist generation type.

    drz_root_dir : string
        Root directory of drizzled images.

    hla_flag_msk : numpy.ndarray object
        mask array used by hla_nexp_flags().

    ci_lookup_file_path : string
        final path elements of the concentration index lookup file

    output_custom_pars_file : string
        name of the output config file

    log_level : int
        The desired level of verboseness in the log statements displayed on the screen and written to the .log file.

    diagnostic_mode : bool
        write intermediate files?

    Returns
    -------
    catalog_data : astropy.Table object
        drizzled filter product catalog data with updated flag values
    """
    # set logging level to user-specified level
    log.setLevel(log_level)

    # Relevant equivalent column titles for aperture and segment catalogs
    all_column_titles = {
        "aperture": {
            "x_coltitle": "X-Center",
            "y_coltitle": "Y-Center",
        },
        "segment": {
            "x_coltitle": "X-Centroid",
            "y_coltitle": "Y-Centroid",
        }
    }
    if proc_type not in all_column_titles.keys():
        log.error("Unknown proc_type '{}', must be 'aperture' or 'segment'".format(proc_type))
        raise ValueError("Unknown proc_type '{}', must be 'aperture' or 'segment'".format(proc_type))
    column_titles = all_column_titles[proc_type]
    # -----------------------
    # FLAG FILTER PROCESSING
    # -----------------------
    log.info("************************** * * * HLA_FLAG_FILTER * * * **************************")
    # -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
    # Flag sources based on concentration index.
    log.info("Determining concentration indices for sources.")
    log.debug("ci_filter({} {} {} {} {} {} {} {} {} {})".format(drizzled_image, catalog_name, "<CATALOG DATA>",
                                                                proc_type, param_dict, ci_lookup_file_path,
                                                                output_custom_pars_file, column_titles, log_level,
                                                                diagnostic_mode))
    catalog_data = ci_filter(drizzled_image, catalog_name, catalog_data, proc_type, param_dict, ci_lookup_file_path,
                             output_custom_pars_file, column_titles, log_level, diagnostic_mode)

    # -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
    # Flag saturated sources
    log.info("Flagging saturated sources in the catalogs.")
    log.debug("hla_saturation_flags({} {} {} {} {} {} {} {} {})".format(drizzled_image, flt_list, catalog_name,
                                                                        "<Catalog Data>", proc_type, param_dict,
                                                                        plate_scale, column_titles, diagnostic_mode))
    catalog_data = hla_saturation_flags(drizzled_image, flt_list, catalog_name, catalog_data, proc_type, param_dict,
                                        plate_scale, column_titles, diagnostic_mode)

    # -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
    # Flag swarm sources
    log.info("Flagging possible swarm features in catalogs")
    log.debug("hla_swarm_flags({} {} {} {} {} {} {} {} {} {})".format(drizzled_image, catalog_name, "<Catalog Data>",
                                                                      exptime, plate_scale, median_sky, proc_type,
                                                                      param_dict, column_titles, diagnostic_mode))
    catalog_data = hla_swarm_flags(drizzled_image, catalog_name, catalog_data, exptime, plate_scale, median_sky,
                                   proc_type, param_dict, column_titles, diagnostic_mode)

    # -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
    # Flag sources from regions where there are a low (or a null) number of contributing exposures
    log.info("Flagging sources from regions observed with only a small number of exposures.")
    log.debug("hla_nexp_flags({} {} {} {} {} {} {} {} {} {})".format(drizzled_image, flt_list, param_dict, plate_scale,
                                                                     catalog_name, "<Catalog Data>", drz_root_dir,
                                                                     "<MASK_ARRAY>", column_titles, diagnostic_mode))
    catalog_data = hla_nexp_flags(drizzled_image, flt_list, param_dict, plate_scale, catalog_name, catalog_data,
                                  drz_root_dir, hla_flag_msk, column_titles, diagnostic_mode)

    display_catalog_bit_populations(catalog_data['Flags'])
    return catalog_data

# ======================================================================================================================


def ci_filter(drizzled_image, catalog_name, catalog_data, proc_type, param_dict, ci_lookup_file_path,
              output_custom_pars_file, column_titles, log_level, diagnostic_mode):
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

    ci_lookup_file_path : string
        final path elements of the concentration index lookup file

    output_custom_pars_file : string
        name of the output config file

    column_titles : dictionary
        Relevant column titles

    log_level : int
        The desired level of verboseness in the log statements displayed on the screen and written to the .log file.

    diagnostic_mode : bool
        write intermediate files?

    Returns
    -------
    catalog_data : astropy.Table object
        drizzled filter product catalog data with updated flag values
    """
    catalog_name_root = catalog_name.split('.')[0]
    ci_lower_limit = float(param_dict['quality control']['ci filter'][proc_type]['ci_lower_limit'])
    ci_upper_limit = float(param_dict['quality control']['ci filter'][proc_type]['ci_upper_limit'])
    snr = float(param_dict['quality control']['ci filter'][proc_type]['bthresh'])

    # replace CI limits with values from table if possible
    cidict = ci_table.get_ci_from_file(drizzled_image, ci_lookup_file_path, log_level,
                                       diagnostic_mode=diagnostic_mode, ci_lower=ci_lower_limit,
                                       ci_upper=ci_upper_limit)  # TODO: add values for ACS/SBC
    ci_lower_limit = cidict['ci_lower_limit']
    ci_upper_limit = cidict['ci_upper_limit']

    # if an output custom param file was created and the CI values were updated by ci_table.get_ci_from_file,
    # update output custom param file with new CI values
    if output_custom_pars_file:
        if ci_upper_limit != float(param_dict['quality control']['ci filter'][proc_type]['ci_lower_limit']) or \
                ci_upper_limit != float(param_dict['quality control']['ci filter'][proc_type]['ci_upper_limit']):
            log.info("CI limits updated.")
            with open(output_custom_pars_file) as f:
                json_data = json.load(f)
            if ci_lookup_file_path.startswith("default"):
                param_set = "default_values"
            else:
                param_set = "parameters"

            if ci_lower_limit != float(param_dict['quality control']['ci filter'][proc_type]['ci_lower_limit']):
                json_data[drizzled_image[:-9]][param_set]["quality control"]["ci filter"][proc_type]["ci_lower_limit"]\
                    = ci_lower_limit

            if ci_upper_limit != float(param_dict['quality control']['ci filter'][proc_type]['ci_upper_limit']):
                json_data[drizzled_image[:-9]][param_set]["quality control"]["ci filter"][proc_type]["ci_upper_limit"]\
                    = ci_upper_limit

            with open(output_custom_pars_file, 'w') as f:
                json.dump(json_data, f, indent=4)
            log.info("Updated custom pars file {}".format(output_custom_pars_file))

    log.info(' ')
    log.info('ci limits for {}'.format(drizzled_image))
    log.info('ci_lower_limit = {}'.format(ci_lower_limit))
    log.info('ci_upper_limit = {}'.format(ci_upper_limit))
    log.info(' ')

    failed_index_list = []
    for i, table_row in enumerate(catalog_data):
        try:
            table_row["Flags"] = int(table_row["Flags"])
        except ValueError:
            table_row["Flags"] = 0

        ci_value = table_row["CI"]
        if ci_value:
            ci_value = float(ci_value)
        merr1 = table_row["MagErrAp1"]
        if not merr1:
            merr1 = numpy.nan
        else:
            merr1 = float(merr1)
        merr2 = table_row["MagErrAp2"]
        if not merr2:
            merr2 = numpy.nan
        else:
            merr2 = float(merr2)
        good_snr = merr2 <= 2.5 / (snr * numpy.log(10))
        ci_err = numpy.sqrt(merr1 ** 2 + merr2 ** 2)

        if not good_snr:
            table_row["Flags"] |= 8

        if not ci_value or (not numpy.isfinite(ci_err)) or ci_value < ci_lower_limit - ci_err:
            table_row["Flags"] |= 16

        if not ci_value or ci_value > ci_upper_limit:
            table_row["Flags"] |= 1

        if not ci_value and diagnostic_mode:
            failed_index_list.append(i)

    if diagnostic_mode:
        # Write out list of ONLY failed rows to to file
        catalog_name_failed = catalog_name_root + '_Failed-CI.txt'
        catalog_data_failed = catalog_data.copy()
        all_indicies = range(0, len(catalog_data))
        rows_to_remove = [z for z in all_indicies if z not in failed_index_list]
        catalog_data_failed.remove_rows(rows_to_remove)
        catalog_data_failed.write(catalog_name_failed, delimiter=",", format='ascii')

        # Write out intermediate catalog with updated flags
        catalog_name = catalog_name_root + 'CIFILT.txt'
        catalog_data.write(catalog_name, delimiter=",", format='ascii')
    return catalog_data

# ======================================================================================================================


def hla_saturation_flags(drizzled_image, flt_list, catalog_name, catalog_data, proc_type, param_dict, plate_scale,
                         column_titles, diagnostic_mode):
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

    plate_scale : float
        plate scale, in arcseconds/pixel

    column_titles : dictionary
        Relevant column titles

    diagnostic_mode : bool
        write intermediate files?

    Returns
    -------
    phot_table_rows : astropy.Table object
        drizzled filter product catalog data with updated flag values
    """
    image_split = drizzled_image.split('/')[-1]
    channel = drizzled_image.split("_")[4].upper()

    if channel == 'IR':  # TODO: Test and IR case just to make sure that IR shouldn't be skipped
        return catalog_data

    # -------------------------------------------------------------------
    # STEP THROUGH EACH APPLICABLE FLT IMAGE, DETERMINE THE COORDINATES
    # FOR ALL SATURATION FLAGGED PIXELS, AND TRANSFORM THESE COORDINATES
    # INTO THE DRIZZLED IMAGE REFERENCE FRAME.
    # -------------------------------------------------------------------
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
        if channel.lower() in ['wfc', 'uvis']:
            image_ext_list = ["[sci,1]", "[sci,2]"]
        if channel.lower() in ['sbc', 'hrc']:
            image_ext_list = ["[sci,1]"]
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
                if ((channel.lower() != 'wfpc2') and (channel.lower() != 'pc')):
                    flt_data = fits.getdata(flt_image, 'DQ', int(ext_part))
                if ((channel.lower() == 'wfpc2') or (channel.lower() == 'pc')):
                    flt_data = fits.getdata(flt_image.replace("_c0m", "_c1m"), 'SCI', int(ext_part))
            except KeyError:
                log.info(' ')
                log.info('WARNING: There is only one set of file extensions in {}'.format(flt_image))
                log.info(' ')

                continue

            # TODO: Should we also look for pixels flagged with DQ value 2048 (A to D saturation) for ACS data?

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
            if diagnostic_mode:
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
    full_sat_list = numpy.concatenate(drz_sat_xy_coords_list)

    # --------------------------------------------
    # WRITE RA & DEC FLT CONVERTED X & Y DRIZZLED
    # IMAGE COORDINATES TO A TEXT FILE
    # --------------------------------------------
    if diagnostic_mode:
        drz_coord_file = drizzled_image.split('/')[-1].split('.')[0] + '_ALL_FLT_SAT_FLAG_PIX.txt'
        drz_coord_out = open(drz_coord_file, 'w')
        for coord in full_sat_list:
            drz_coord_out.write(str(coord[0]) + '     ' + str(coord[1]) + '\n')
        drz_coord_out.close()

    # ----------------------------------------------------
    # GET SOURCELIST X AND Y VALUES
    # ----------------------------------------------------
    all_detections = catalog_data

    nrows = len(all_detections)
    full_coord_list = numpy.empty((nrows, 2), dtype=numpy.float)
    for row_count, detection in enumerate(all_detections):
        full_coord_list[row_count, 0] = float(detection[column_titles["x_coltitle"]])
        full_coord_list[row_count, 1] = float(detection[column_titles["y_coltitle"]])

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

    radius = round((ap2/plate_scale) + 0.5) * 2.

    log.info(' ')
    log.info('THE RADIAL DISTANCE BEING USED IS {} PIXELS'.format(str(radius)))
    log.info(' ')

    # do the cross-match using xymatch
    log.info('Matching {} saturated pixels with {} catalog sources'.format(len(full_sat_list), len(full_coord_list)))
    psat, pfull = xymatch(full_sat_list, full_coord_list, radius, multiple=True, verbose=False)
    log.info('Found cross-matches (including duplicates)'.format(len(psat)))
    saturation_flag = numpy.zeros(len(full_coord_list), dtype=bool)
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

        if diagnostic_mode:
            sat_coord_file = drizzled_image.split('/')[-1].split('.')[0] + '_INTERMEDIATE.txt'
            sat_coord_out = open(sat_coord_file, 'w')
            for sat_coord in full_coord_list[saturation_flag, :]:
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
                table_row["Flags"] = int(table_row["Flags"]) | 4

        phot_table_rows = flag4and8_hunter_killer(phot_table_rows, column_titles)

        if diagnostic_mode:
            phot_table_temp = phot_table_root + '_SATFILT.txt'
            phot_table_rows.write(phot_table_temp, delimiter=",", format='ascii')
        return phot_table_rows

# ======================================================================================================================


def hla_swarm_flags(drizzled_image, catalog_name, catalog_data, exptime, plate_scale, median_sky, proc_type, param_dict,
                    column_titles, diagnostic_mode):

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

    plate_scale : float
        plate scale, in arcseconds/pixel

    median_sky : float
        median sky value

    proc_type : string
        sourcelist generation type.

    param_dict : dictionary
        Dictionary of instrument/detector - specific drizzle, source finding and photometric parameters

    column_titles : dictionary
        Relevant column titles

    diagnostic_mode : bool
        write intermediate files?

    Returns
    -------
    catalog_data : astropy.Table object
        drizzled filter product catalog data with updated flag values
    """
    drz_img_path_split = drizzled_image.split('/')
    drz_img_split = drz_img_path_split[-1].split('_')
    data_type = drz_img_split[4]

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

    ap2 = param_dict['catalog generation']['aperture_2']
    if proc_type not in ('segment', 'aperture'):
        log.error("Unknown catalog type '{}', must be 'aperture' or 'segment'".format(proc_type))
        raise ValueError("Unknown catalog type '%s'" % proc_type)

    # ----------------------------------
    # Convert aperture radius to pixels
    # ----------------------------------
    radius = ap2 / plate_scale
    log.info(' ')
    log.info('Aperture Size = {}'.format(ap2))
    log.info('Pixel Scale = {} arcsec per pixel'.format(plate_scale))
    log.info(' ')
    area = math.pi * radius**2

    nrows = len(catalog_data)

    complete_src_list = numpy.empty((nrows, 6), dtype=numpy.float)

    for row_num, row in enumerate(catalog_data[0:]):
        x_val = float(row[column_titles["x_coltitle"]])
        y_val = float(row[column_titles["y_coltitle"]])
        flux = row["FluxAp2"]
        sky = row["MSkyAp2"]

        if not flux:
            flux = 0.0

        if not sky:
            sky = 0.0

        electronpp = flux / area * exptime
        eppsky = electronpp / median_sky
        complete_src_list[row_num, :] = [x_val, y_val, flux, electronpp, sky, eppsky]

    if len(complete_src_list) == 0:
        return catalog_data

    # view into the complete_src_list array for convenience
    swarm_epp_list_a = complete_src_list[:, 3]

    # swarm flag array
    # this will get set as candidates to flag are found
    swarm_flag = numpy.zeros(nrows, dtype=bool)

    # ------------------------------------------------------------
    # WRITE SUBSET SOURCE LIST TO AN OUTPUT FILE FOR VERIFICATION
    # ------------------------------------------------------------
    if diagnostic_mode:
        final_complete_source_file = open(phot_table_root+'_SWFILT_COMPLETE_SOURCE_FILE.txt', 'w')
        final_complete_source_file.write("# {}\n".format("-"*96))
        swfilt_table_header = "# X-Center   Y-Center     Flux        ElectronPP          Sky         EPPSKY_Ratio \n"
        final_complete_source_file.write(swfilt_table_header)
        final_complete_source_file.write("# {}\n".format("-"*96))
        for i, complete_src_value in enumerate(complete_src_list):
            final_complete_source_file.write(str(complete_src_value[0]) + '     ' +
                                             str(complete_src_value[1]) + '     ' +
                                             str(complete_src_value[2]) + '     ' +
                                             str(complete_src_value[3]) + '     ' +
                                             str(complete_src_value[4]) + '     ' +
                                             str(complete_src_value[5]) + '\n')

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
    selfradius = float(param_dict["quality control"]["swarm filter"]["selfradius"])  # TODO: optimize selfradius values for ACS/HRC, ACS/SBC in quality control param files

    eppsky_limit = eppsky_limit_cfg * median_sky

    # ----------------------------------------------------------
    # UVIS --> EPP > 100000. OR (EPP > 1000*sky AND EPP > 10000)
    # IR   --> EPP > 100000. OR (EPP > 100*sky AND EPP > 10000)
    # ----------------------------------------------------------

    initial_central_pixel_positions = numpy.where(numpy.logical_or(swarm_epp_list_a > upper_epp_limit,
                                                  numpy.logical_and(swarm_epp_list_a > eppsky_limit,
                                                                    swarm_epp_list_a > lower_epp_limit)))[0]
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
            cuts = param_dict["quality control"]["swarm filter"]["cuts_list"]
            cuts = list(map(float, cuts))

            selfradii = param_dict["quality control"]["swarm filter"]["selfradii_list"]
            selfradii = list(map(float, selfradii))

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
                                                                        initial_central_pixel_list[:, 3:4] <=
                                                                        cuts[cut_cnt-1]))[0]

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
                                         initial_central_pixel_list[:, 0:2], selfradii[cut_cnt],
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
                        exclude_index = numpy.concatenate((exclude_index,
                                                           i2[numpy.where(flux2 < numpy.max(flux2))]), axis=0)

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
    if diagnostic_mode:
        cetrl_pix_pos_file = phot_table_root + '_SWFILT_CENTRAL-PIX-POS.txt'
        drz_coord_out = open(cetrl_pix_pos_file, 'w')
        for i in range(len(final_flag_src_central_pixel_list)):
            drz_coord_out.write(str(final_flag_src_central_pixel_list[i, 0]) + '     ' +
                                str(final_flag_src_central_pixel_list[i, 1]) + '     ' +
                                str(final_flag_src_central_pixel_list[i, 2]) + '     ' +
                                str(final_flag_src_central_pixel_list[i, 3]) + '     ' +
                                str(final_flag_src_central_pixel_list[i, 4]) + '     ' +
                                str(final_flag_src_central_pixel_list[i, 5]) + '\n')
        drz_coord_out.close()

    # ==========================================================================
    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # EXTRACT THE CENTRAL PIXEL POSITIONS IN final_flag_src_central_pixel_list,
    # FROM swarm_x_list_b AND swarm_y_list_b
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
    swarm_list_b = complete_src_list[notcentral_index, :]

    # views into the swarm_list_b array for convenience
    swarm_x_list_b = swarm_list_b[:, 0]
    swarm_y_list_b = swarm_list_b[:, 1]

    # ---------------------------------------------------------------------
    # ITERATIVELY CLIP SOURCES CONTAINED WITHIN RINGS AT SPECIFIED RADIUS
    # VALUES, PROGRESSIVELY MOVING CLOSER TO THE CENTRAL SOURCE
    # ---------------------------------------------------------------------

    # do the cross-match using xymatch
    log.info('Matching {} swarm centers with {} catalog sources'.format(len(final_flag_src_central_pixel_list),
                                                                        len(swarm_list_b)))
    pcentral, pfull = xymatch(final_flag_src_central_pixel_list[:, 0:2], swarm_list_b[:, 0:2],
                              clip_radius_list[0], multiple=True, stack=False, verbose=False)

    # TODO: RLW: the ring list is needed only for testing, get rid of it when code works

    if diagnostic_mode:
        ring_index_list = []
        ring_refepp_list = []
        ring_thresh_list = []
        ring_count = []

    for pindex, ii in enumerate(pcentral):

        central_pixel_value = final_flag_src_central_pixel_list[ii, :]
        log.debug(' ')
        log.debug('CENTRAL PIXEL VALUE: {}'.format(central_pixel_value))
        log.debug(' ')

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

        distsq = (swarm_x_list_b[allmatches]-coords[0])**2 + (swarm_y_list_b[allmatches]-coords[1])**2
        sind = distsq.argsort()
        allmatches = allmatches[sind]
        distsq = distsq[sind]
        rcut = distsq.searchsorted(numpy.array(clip_radius_list)**2)
        for radius_cnt in range(1, len(clip_radius_list)):

            # -------------------------------------------
            # ISOLATE THE DETECTIONS WITHIN A GIVEN RING
            # -------------------------------------------

            matches = allmatches[rcut[radius_cnt]:rcut[radius_cnt-1]]

            if len(matches) == 0:
                log.debug(' ')
                log.debug('------------------------------------------')
                log.debug('NOTE: ALL MATCHES/DETECTIONS IN THIS RING ')
                log.debug('      HAVE PREVIOUSLY BEEN ACCOUNTED FOR  ')
                log.debug('------------------------------------------')
                log.debug(' ')

                continue

            # -----------------------------------------------------------
            # CALCULATE THE MEDIAN SKY VALUE FOR THE GROUP OF DETECTIONS
            # CONTAINED WITHIN THE SPECIFIED RING BEING PROCESSED
            # -----------------------------------------------------------
            ref_epp = base_epp * scale_factor_list[radius_cnt-1]

            # -----------------------------------------------------------------------------------
            # DIFFERENTIATE BETWEEN GOOD DETECTIONS AND SWARM DETECTIONS WITHIN SPECIFIED RINGS
            # -----------------------------------------------------------------------------------
            ring = swarm_list_b[matches, :]
            w = numpy.where(ring[:, 3]/ref_epp < swarm_thresh)
            if len(w) > 0:
                swarm_flag[notcentral_index[matches[w]]] = True

            # TODO: RLW: following needed only for testing, get rid of it when code works
            if diagnostic_mode:
                ring_index_list.append(matches)
                ring_count.append(len(matches))
                ring_refepp_list.append(ring[:, 3]/ref_epp)
                ring_thresh_list.append(swarm_thresh)

    # TODO: RLW: following needed only for testing, get rid of it when code works
    if diagnostic_mode:
        # -----------------------------------------------------------------------------------------
        # WRITE CLIPPED SOURCES CONTAINED WITHIN RINGS TO AN OUTPUT FILE FOR INTERMEDIATE ANALYSIS
        # -----------------------------------------------------------------------------------------
        ring_source_file = phot_table_root+'_SWFILT_RING-SOURCE-INFO.txt'
        ring_src_outfile = open(ring_source_file, 'w')
        ring_src_outfile.write("# {}\n".format("-"*96))
        swfilt_ring_file_header = "# X-Center   Y-Center     Flux        ElectronPP"
        swfilt_ring_file_header += "          Sky        SrcEPP/RefEPP   Swarm Thresh \n"
        ring_src_outfile.write(swfilt_ring_file_header)
        ring_src_outfile.write("# {}\n".format("-"*96))

        if ring_index_list:
            ring_index_list = numpy.concatenate(ring_index_list)

            # select just the lowest value of refepp/swarm threshold for each source
            # create array with extra columns
            ring_source_list = numpy.empty((len(ring_index_list), 9), dtype=numpy.float)
            ring_source_list[:, 0:6] = swarm_list_b[ring_index_list, :]
            ring_source_list[:, 6] = numpy.concatenate(ring_refepp_list)
            ring_source_list[:, 7] = numpy.repeat(ring_thresh_list, ring_count)
            ring_source_list[:, 8] = ring_source_list[:, 6] / ring_source_list[:, 7]

            # sort by x, y, and refepp
            # tricky here: get a view with named columns, then specify names as sort items
            ring_source_list.view(','.join(['f8']*9)).sort(order=['f0', 'f1', 'f8'], axis=0)

            # keep just first entry when the same source appears more than once
            keep = numpy.ones(len(ring_index_list), dtype=bool)
            keep[1:] = numpy.logical_or(ring_source_list[1:, 0] != ring_source_list[:-1, 0],
                                        ring_source_list[1:, 1] != ring_source_list[:-1, 1])
            ring_source_list = ring_source_list[keep, :]

            for ring_source in ring_source_list:
                ring_src_outfile.write(str(ring_source[0]) + '     ' +
                                       str(ring_source[1]) + '     ' +
                                       str(ring_source[2]) + '     ' +
                                       str(ring_source[3]) + '     ' +
                                       str(ring_source[4]) + '     ' +
                                       str(ring_source[5]) + '     ' +
                                       str(ring_source[6]) + '     ' +
                                       str(ring_source[7]) + '\n')
        ring_src_outfile.close()
        # XXX RLW: end of testing code

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
            ctr_list_radius_list = param_dict["quality control"]["swarm filter"]["ctrList_radiusList"]  # TODO: optimize ctr_list_radius_list for ACS wfc, hrc, sbc in quality control config files
            ctr_list_radius_list = list(map(int, ctr_list_radius_list))

            ctr_list_threshold_list = param_dict["quality control"]["swarm filter"]["ctrList_thresholdList"]  # TODO: optimize ctr_list_threshold_list for ACS wfc, hrc, sbc in quality control config files
            ctr_list_threshold_list = list(map(int, ctr_list_threshold_list))

            for ctr_list_cnt, (threshold, radius) in enumerate(zip(ctr_list_threshold_list, ctr_list_radius_list)):

                if ctr_list_cnt == 0:
                    ctr_list_cut = final_flag_src_central_pixel_list[:, 3] > threshold
                else:
                    ctr_list_cut = numpy.logical_and(final_flag_src_central_pixel_list[:, 3] > threshold,
                                                     final_flag_src_central_pixel_list[:, 3] <=
                                                     ctr_list_threshold_list[ctr_list_cnt-1])

                ctr_list_cut1 = final_flag_src_central_pixel_list[ctr_list_cut, :]
                pcentral, pfull = xymatch(ctr_list_cut1[:, 0:2], swarm_list_b[:, 0:2],
                                          radius, multiple=True, verbose=False)
                proximity_flag[notcentral_index[pfull]] = True

        log.info("Proximity filter flagged {} sources".format(proximity_flag.sum()))

        # --------------------------------------------------------------------------
        # WRITE NEAR CENTRAL POSITION SWARM LIST TO AN OUTPUT FILE FOR VERIFICATION
        # --------------------------------------------------------------------------
        if diagnostic_mode:
            near_swm_list = complete_src_list[proximity_flag, :]
            final_near_swarm_file = open(phot_table_root+'_SWFILT_NEAR_SWARM_FILE.txt', 'w')
            for swarm_value in near_swm_list:
                final_near_swarm_file.write(str(swarm_value[0]) + '     ' +
                                            str(swarm_value[1]) + '     ' +
                                            str(swarm_value[2]) + '     ' +
                                            str(swarm_value[3]) + '     ' +
                                            str(swarm_value[4]) + '\n')
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
    if diagnostic_mode:
        final_swarm_file = open(phot_table_root+'_SWFILT_SWARM_FILE.txt', 'w')
        for swarm_value in final_swarm_list:
            final_swarm_file.write(str(swarm_value[0]) + '     ' +
                                   str(swarm_value[1]) + '     ' +
                                   str(swarm_value[2]) + '     ' +
                                   str(swarm_value[3]) + '     ' +
                                   str(swarm_value[4]) + '\n')
        final_swarm_file.close()

    # ----------------------------------------------------
    # WRITE SOURCE LIST TO AN OUTPUT FILE FOR VERIFICATION
    # ----------------------------------------------------
    if diagnostic_mode:
        final_source_file = open(phot_table_root+'_SWFILT_SOURCE_FILE.txt', 'w')
        for source_value in final_source_list:
            final_source_file.write(str(source_value[0]) + '     ' +
                                    str(source_value[1]) + '     ' +
                                    str(source_value[2]) + '     ' +
                                    str(source_value[3]) + '     ' +
                                    str(source_value[4]) + '\n')
        final_source_file.close()

    # Update catalog_data flag values
    for i, table_row in enumerate(catalog_data[0:]):
        if combined_flag[i]:
            table_row["Flags"] |= 32

    if diagnostic_mode:
        # Write out intermediate catalog with updated flags
        phot_table_temp = phot_table_root + '_SWFILT.txt'
        catalog_data.write(phot_table_temp, delimiter=",", format='ascii')

    return catalog_data

# ======================================================================================================================


def hla_nexp_flags(drizzled_image, flt_list, param_dict, plate_scale, catalog_name, catalog_data, drz_root_dir,
                   mask_data, column_titles, diagnostic_mode):
    """flags out sources from regions where there are a low (or a null) number of contributing exposures

    drizzled_image : string
        Name of drizzled image to process

    flt_list : list
        list of calibrated images that were drizzle-combined to produce image specified by input parameter
        'drizzled_image'

    param_dict : dictionary
        Dictionary of instrument/detector - specific drizzle, source finding and photometric parameters.

    plate_scale : float
        plate scale, in arcseconds/pixel

    catalog_name : string
        drizzled filter product catalog filename to process

    catalog_data : astropy.Table object
        drizzled filter product catalog data to process

    drz_root_dir :
        dictionary of source lists keyed by drizzled image name.

    mask_data : numpy.ndarray object
        mask array used by hla_nexp_flags().

    column_titles : dictionary
        Relevant column titles

    diagnostic_mode : bool
        write intermediate files?

    Returns
    -------
    catalog_data : astropy.Table object
        drizzled filter product catalog data with updated flag values
    """
    # ------------------
    # CREATE NEXP IMAGE
    # ------------------
    channel = drizzled_image.split("_")[4].upper()

    # if channel == 'IR':  # TODO: This was commented out in the HLA classic era, prior to adaption to the HAP pipeline. Ask Rick about it.
    #    return catalog_data

    drz_data = fits.getdata(drizzled_image, 1)

    component_drz_img_list = get_component_drz_list(drizzled_image, drz_root_dir, flt_list)
    nx = drz_data.shape[0]
    ny = drz_data.shape[1]
    nexp_array = numpy.zeros((nx, ny), dtype=numpy.int32)

    for comp_drz_img in component_drz_img_list:
        comp_drz_data = (fits.getdata(comp_drz_img) != 0).astype(numpy.int32)
        try:
            nexp_array += comp_drz_data
        except ValueError:
            log.info("WARNING: Astrodrizzle added an extra-row/column...")
            nexp_array += comp_drz_data[0:nx, 0:ny]

    # this bit is added to get the mask integrated into the exp map
    mask_array = (mask_data == 0.0).astype(numpy.int32)
    nexp_array = nexp_array * mask_array
    # -------------------------------------------------------
    # EXTRACT FLUX/NEXP INFORMATION FROM NEXP IMAGE BASED ON
    # THE SOURCE DETECTION POSITIONS PREVIOUSLY ESTABLISHED
    # -------------------------------------------------------

    phot_table_root = catalog_name.split('/')[-1].split('.')[0]

    nrows = len(catalog_data)
    cat_coords = numpy.empty((nrows, 2), dtype=float)
    for line_cnt, phot_table_line in enumerate(catalog_data):
        x_coord = phot_table_line[column_titles["x_coltitle"]]
        y_coord = phot_table_line[column_titles["y_coltitle"]]
        cat_coords[line_cnt, :] = [x_coord, y_coord]
    # ----------------------------------
    # Convert aperture radius to pixels
    # ----------------------------------

    ap2 = param_dict['catalog generation']['aperture_2']
    radius = ap2/plate_scale

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

    icoords = (cat_coords+0.5).astype(int)
    # note x & y are swapped so they can index the numpy array nexp_array
    # catalog x is second subscript, catalog y is first subscript
    ix = icoords[:, 1]
    iy = icoords[:, 0]

    # get list of neighboring pixels that are within radius
    iradius = int(radius+1)
    idiam = iradius*2+1
    gx, gy = numpy.mgrid[0:idiam, 0:idiam] - iradius
    gx = gx.ravel()
    gy = gy.ravel()
    w = numpy.where(gx**2+gy**2 <= radius**2)[0]
    gx = gx[w]
    gy = gy[w]

    # check the pixel values for low nexp

    # this version uses numpy broadcasting sum gx+ix is [len(gx), nrows]
    gx = (gx[:, numpy.newaxis] + ix).clip(0, nexp_array.shape[0]-1)
    gy = (gy[:, numpy.newaxis] + iy).clip(0, nexp_array.shape[1]-1)
    artifact_flag = nexp_array[gx, gy].min(axis=0) < artifact_filt

    log.info('FLAGGING {} OF {} SOURCES'.format(artifact_flag.sum(), nrows))

    # Add flag bit to appropriate sources
    for i, table_row in enumerate(catalog_data):
        if artifact_flag[i]:
            table_row["Flags"] |= 64

    if diagnostic_mode:
        # Write out intermediate catalog with updated flags
        phot_table_temp = phot_table_root + '_NEXPFILT.txt'
        catalog_data.write(phot_table_temp, delimiter=",", format='ascii')

    return catalog_data

# ======================================================================================================================


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
    drizzle_file_suffex = drizzled_image[-8:-5]
    drz_img_split = drizzled_image.split('/')[-1].split("_"+drizzle_file_suffex)
    component_drz_img_list = glob.glob(os.path.join(drz_root_dir,
                                                    drz_img_split[0])+'*_{}.fits'.format(drizzle_file_suffex))
    component_drz_img_list.sort()
    for item in component_drz_img_list:
        if item.endswith(drizzled_image):
            component_drz_img_list.remove(item)
    drz_filter = drizzled_image.split("_")[5]

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
        fh = fits.open(drzfile)
        rootname = fh[0].header.get('rootname', '')
        fh.close()
        fname = os.path.split(rootname)[-1]
        fname = fname.split('_')[0].lower()
        if fname in ipdict:
            rv.append(drzfile)
    if len(list_of_flts) != len(rv):
        # this must be a bug?
        log.info("ERROR: mismatch after filtering in exposure lists for {}".format(drz_filter))
        log.info("Filtered drizzled exposure list: {}".format("\n".join(rv)))
        log.info("flt exposure list: {}".format("\n".join(list_of_flts)))
        log.info("Plowing ahead with the filtered drizzled list")
    return rv

# =============================================================================
# -----------------------------------------------------------------------------
# ------------------------ AUXILIARY FUNCTIONS BELOW --------------------------
# -----------------------------------------------------------------------------
# =============================================================================


def xymatch(cat1, cat2, sep, multiple=False, stack=True, verbose=True):
    """Routine to match two lists of objects by position using 2-D Cartesian distances.

    Matches positions in cat1 with positions in cat2, for matches within separation (sep).
    If more than one match is found, the nearest is returned.
    Setting multiple=True returns all matching pairs rather than just the closest.

    Input catalogs need not be sorted. They should be 2-element arrays [:, 2] with
    cat1[:, 0] = x1 and cat1[:, 1] = y1.

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
        If multiple is true, returns a tuple (p1,p2) such that cat1[p1] and cat2[p2] are within sep.
        p1 and p2 may include multiple pointers to the same objects in cat1 or cat2.  In this case objects that don't
        match are simply omitted from the lists. Default value is 'False'.

    stack : Boolean
        If stack is true, the returned matching pointers are stacked into a single big array, so both p1 and p2 are 1-D
        arrays of length nmatches. Default value is 'True'.

    verbose : Boolean
        print verbose output? Default value is 'True'.

    Returns
    -------
    Varies; Depending on inputs, either just 'p2', or 'p1' and 'p2'. p1 and p2 are lists of matched indices
    """
    if not (isinstance(cat1, numpy.ndarray) and len(cat1.shape) == 2 and cat1.shape[1] == 2):
        log.error("catalog 1 must be a [N, 2] array")
        raise ValueError("cat1 must be a [N, 2] array")
    if not (isinstance(cat2, numpy.ndarray) and len(cat2.shape) == 2 and cat2.shape[1] == 2):
        log.error("catalog 2 must be a [N, 2] array")
        raise ValueError("cat2 must be a [N, 2] array")

    x1 = cat1[:, 0]
    y1 = cat1[:, 1]
    x2 = cat2[:, 0]
    y2 = cat2[:, 1]

    # Sort the arrays by increasing y-coordinate
    is1 = y1.argsort()
    x1 = x1[is1]
    y1 = y1[is1]
    is2 = y2.argsort()
    x2 = x2[is2]
    y2 = y2[is2]
    # find search limits in y2 for each object in y1
    kvlo = y2.searchsorted(y1-sep, 'left').clip(0, len(y2))
    kvhi = y2.searchsorted(y1+sep, 'right').clip(kvlo, len(y2))

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
                        p1.append(numpy.zeros(len(ww), dtype='int')+is1[i])
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
                return numpy.array([], dtype=int), numpy.array([], dtype=int)
            else:
                return numpy.concatenate(p1), numpy.concatenate(p2)
        else:
            return (p1, p2)
    else:
        return p2

# ======================================================================================================================


def rdtoxy(rd_coord_array, image, image_ext):
    """converts RA and dec to x, y image coords.

    rd_coord_array : numpy.ndarray
        array containing RA and dec values to convert.

    image : string
        drizzled image whose WCS info will be used in the coordinate conversion.

    image_ext : string
        fits image extension to be used in the conversion.

    Returns
    xy_arr: array
        array of converted x, y coordinate value pairs
    """

    scifile = image + image_ext
    wcs = wcsutil.HSTWCS(scifile)
    try:
        xy_arr = wcs.wcs_sky2pix(rd_coord_array, 1)
    except AttributeError:
        xy_arr = wcs.wcs_world2pix(rd_coord_array, 1)
    return (xy_arr)

# ======================================================================================================================


def xytord(xy_coord_array, image, image_ext):
    """converts x, y image coords to RA and dec.

    xy_coord_array : numpy.ndarray
        array containing image x, y coord values to convert.

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
        rd_arr = wcs.all_pix2sky(xy_coord_array, 1)
    except AttributeError:
        rd_arr = wcs.all_pix2world(xy_coord_array, 1)
    return (rd_arr)

# ======================================================================================================================


def flag4and8_hunter_killer(catalog_data, column_titles):
    """This function searches through photometry catalogs for sources whose flags contain
    both bits 4 (multi-pixel saturation), and 8 (faint magnitude limit).
    If found, the subroutine removes the "8" bit value from the set of flags for that source.

    Parameters
    ----------
    catalog_data : astropy Table object
        catalog data to process

    column_titles : dictionary
        Relevant column titles

    Returns
    -------
    catalog_data : astropy Table object
        input catalog data with updated flags
    """
    conf_ctr = 0
    log.info("Searching for flag 4 + flag 8 conflicts....")
    for catalog_line in catalog_data:
        if ((catalog_line["Flags"] & 4 > 0) and (catalog_line["Flags"] & 8 > 0)):
            conf_ctr += 1
            catalog_line["Flags"] = int(catalog_line["Flags"]) - 8
    if conf_ctr == 0:
        log.info("No conflicts found.")
    if conf_ctr == 1:
        log.info("{} conflict fixed.".format(conf_ctr))
    if conf_ctr > 1:
        log.info("{} conflicts fixed.".format(conf_ctr))

    return catalog_data

# ======================================================================================================================


def make_mask_array(drz_image):
    """
    Creates _msk.fits mask file that contains pixel values of 1 outside the drizzled image footprint and pixel values
    of 0 inside the footprint. This file is used by subroutine hla_nexp_flags().

    Parameters
    ----------
    drz_image : string
        drizzled image filename

    Returns
    -------
    mask : numpy.ndarray object
        mask array
    """
    mask = fits.open(drz_image)[1].data != 0
    dilate = scipy.ndimage.morphology.binary_dilation
    erode = scipy.ndimage.morphology.binary_erosion
    kernel1 = numpy.ones((25, 25), dtype=int)
    kernel2 = numpy.ones((31, 31), dtype=int)
    # add padding around the edge so pixels close to image boundary are correct
    padding = 13
    bigmask = numpy.pad(mask, padding, 'constant')
    # strip the padding back off after creating mask
    mask = (erode(dilate(bigmask, kernel1), kernel2) == 0)[padding:-padding, padding:-padding]
    mask = mask.astype(numpy.int16)
    return mask


# ======================================================================================================================

def deconstruct_flag(flagval):
    """Breaks down an integer flag value into individual component bit values.

    Parameters
    ----------
    flagval : int
        Flag value to deconstruct

    Returns
    -------
    out_idx_list : list
        a 9-element numpy array of 0s and 1s. Each element of the array represents the presence of a particular
        bit value (element 0 = bit 0, element 1 = bit 1, ..., element 3 = bit 4 and so on...)
    """
    bitlist = [1, 2, 4, 8, 16, 32, 64, 128]
    flagval = int(flagval)
    # out_bit_list = []
    out_idx_list = numpy.zeros(9, dtype=int)
    if flagval == 0:
        # out_bit_list = [0]
        out_idx_list[0] = 1
    if flagval > 0:
        idx = 1
        for bit in bitlist:
            if flagval & bit > 0:
                # out_bit_list.append(bit)
                out_idx_list[idx] = 1
            if bit > flagval:
                break
            idx += 1
    return out_idx_list


# ======================================================================================================================

def display_catalog_bit_populations(flag_data):
    """Breaks all input flag values down into their constituent bit values and displays a bit-by-bit population summary

    Parameters
    ----------
    flag_data : astropy.table.column.Column object
        'Flags' column of a given sourcelist to analyze

    Returns
    -------
    Nothing.
    """
    bit_list = [0, 1, 2, 4, 8, 16, 32, 64, 128]
    flag_meanings = ['Point Source',
                     'Extended Source',
                     'Single-Pixel Saturation',
                     'Multi-Pixel Saturation',
                     'Faint Magnitude Limit',
                     'Hot Pixel',
                     'Swarm Detection',
                     'Edge and Chip Gap',
                     'Bleeding and Cosmic Rays']
    flag_counts = numpy.zeros(9, dtype=int)
    n_sources = len(flag_data)
    for flagval in flag_data:
        flag_counts += deconstruct_flag(flagval)
    max_length = 5
    for bitval in flag_counts:
        max_length = max([max_length, len(str(bitval))])
    log.info("{}".format("-"*60))
    log.info("{}FLAG BIT VALUE POPULATION SUMMARY".format(" "*13))
    log.info("Bit   Meaning{}Count Percentage".format(" "*20))
    fill_char = " "
    for ctr in range(0, len(bit_list)):
        bit_val = bit_list[ctr]
        pct_val = 100.0*(float(flag_counts[ctr])/float(n_sources))
        padding1 = 6 - len(str(bit_val))
        padding2 = 27 - len(flag_meanings[ctr])
        padding3 = max_length-len(str(flag_counts[ctr]))
        if pct_val == 100.:
            padding4 = 3
        elif pct_val >= 10.:
            padding4 = 4
        else:
            padding4 = 5
        log.info("{}{}{}{}{}{}{}{:.3f}%".format(bit_val, fill_char*padding1, flag_meanings[ctr], padding2*fill_char,
                                                fill_char*padding3, flag_counts[ctr], fill_char*padding4, pct_val))
    log.info("{}".format(" -- " * 15))
    log.info("NOTE: As the flag value for a given source can be composed ")
    log.info("of multiple bits, the above percentage values need not add")
    log.info("up to 100%.")
    log.info("{}\n".format("-" * 60))
