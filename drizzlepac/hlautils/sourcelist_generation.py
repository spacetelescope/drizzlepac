#!/usr/bin/env python

"""This script contains code to support creation of source extractor-like and daophot-like sourcelists.

"""
import datetime
import os
import pdb
import sys
import traceback

from astropy.io import fits, ascii
from astropy.stats import sigma_clipped_stats
from astropy.table import Table, Column, MaskedColumn
import numpy
from photutils import aperture_photometry, Background2D, CircularAperture, CircularAnnulus, detection, findstars
from photutils import MedianBackground, SExtractorBackground, StdBackgroundRMS
import scipy
from stsci.tools import logutil

from drizzlepac import util
import hla_flag_filter
from photometry_tools import iraf_style_photometry

__taskname__ = 'sourcelist_generation'

log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout)

# ----------------------------------------------------------------------------------------------------------------------


def add_header_phot_tab(phot_table, drz_image, param_dict):
    """
    This task will make a header for a DAOPhot table based on information contained in the header of the 'drz_image'.

    Tested.

    Parameters
    ----------
    phot_table : string
        DAOPhot sourcelist filename.

    drz_image : string
        filename of drizzled image whose header information will be used to populate the DAOPhot table header.

    param_dict : dictionary
        dictionary of drizzle, source finding and photometric parameters

    Returns
    -------
    nothing.
    """
    loaded_fits = fits.open(drz_image)

    # -------------------------------
    # FITS-compliant UTC time string
    # -------------------------------
    fits_string_date = datetime.datetime.utcnow().strftime("%A-%m-%dT%H:%M:%S")

    # ---------------------------------
    # Human-friendly local time string
    # ---------------------------------
    string_date = datetime.datetime.now().strftime("%A %B %d %H:%M:%S %Y")
    fname = extract_name(drz_image)

    # ---------------------------------------------------
    # Get available header keywords from drizzled image:
    # ("not available" is returned if not found)
    # ---------------------------------------------------
    prop = get_head_val_opened_fits(loaded_fits, "proposid")
    tname = get_head_val_opened_fits(loaded_fits, "targname")
    inst = get_head_val_opened_fits(loaded_fits, "instrume")
    detect = drz_image.split("_")[-2].lower() # TODO: May need to be refactored to adjust for new names, and fact that ACS has two filters
    filt = drz_image.split("_")[-1].replace(".fits","").lower() # TODO: May need to be refactored to adjust for new names, and fact that ACS has two filters
    im_ra = get_head_val_opened_fits(loaded_fits, "crval1", ext=1)
    im_dec = get_head_val_opened_fits(loaded_fits, "crval2", ext=1)
    orient = get_head_val_opened_fits(loaded_fits, "pa_aper", ext=1)

    fwhm = param_dict['dao']['TWEAK_FWHMPSF']
    thresh = param_dict['dao']['TWEAK_THRESHOLD']
    scale = param_dict['astrodrizzle']['SCALE']

    # --------------------
    # Fill in the header:
    # --------------------
    head_outfile = "#DAOPhot Catalog\n#  Processed by the HLA pipeline\n#\n"
    head_outfile = head_outfile + "#---------------------------------------------------#\n"
    head_outfile = head_outfile + "# Data Release Version: %s (DAOPHOT COMPLETE CATALOG)\n" % (
        os.getenv("HLA_BUILD_VER"))
    head_outfile = head_outfile + "# Processed On:         %s\n" % string_date
    head_outfile = head_outfile + "# Proposal ID:          %s\n" % prop
    head_outfile = head_outfile + "# File Name:            %s\n" % fname
    head_outfile = head_outfile + "# Target Name:          %s\n" % tname
    head_outfile = head_outfile + "# Instrument:           %s\n" % inst
    head_outfile = head_outfile + "# Detector:             %s\n" % detect
    head_outfile = head_outfile + "# Image Center RA:      %s\n" % im_ra
    head_outfile = head_outfile + "# Image Center Dec:     %s\n" % im_dec
    head_outfile = head_outfile + "# Filter:               %s\n" % filt
    head_outfile = head_outfile + "# Orientation:          %s\n" % orient
    head_outfile = head_outfile + "# Aperture Correction:  Not Available\n"
    head_outfile = head_outfile + "# CTE_Corr:             No\n"
    head_outfile = head_outfile + "# CTE_Date:             Not Applicable\n"
    head_outfile = head_outfile + "# FWHM:                 %f\n" % fwhm
    head_outfile = head_outfile + "# Threshold:            %f\n" % thresh
    head_outfile = head_outfile + "# Scale:                %f\n" % scale
    head_outfile = head_outfile + "#---------------------------------------------------#\n#\n"

    # -----------------
    # Append to table:
    # -----------------
    dirname = os.path.dirname(os.path.abspath(phot_table))
    rootname = extract_name(phot_table).split("_daophot")[0]
    phot_table_new = find_unique_name(rootname + "_daohead.cat", dirname, 'no', FITS_file_rootname=rootname,
                                             Suffix="daohead.cat")
    #    phot_table_new = Rename.unique_name(rootname + "_daohead.cat",suffix = "daohead.cat")
    phot_table_new = os.path.join(dirname, phot_table_new)
    readlines = [head_outfile]

    # ----------------------------------
    #    save(phot_table_new, readlines)
    newfile = open(phot_table_new, "w")
    rows = len(readlines)
    for row in range(0, rows):
        newfile.write(readlines[row])
    newfile.close()
    # ----------------------------------


# ----------------------------------------------------------------------------------------------------------------------


def average_values_from_dict(Dictionary):
    """Average all values within a dictionary.  This task is used for the 'hla_reduction.py' source listing procedure
    in estimating DAOFind parameters for the white-light image from the parameters for each of  the composite images of
    the white-light.

    Tested.

    Parameters
    ----------
    Dictionary : dictionary
        Input dictionary containing values to average

    Returns
    -------
    final : float
        The average value of the input dictionary
    """
    all_vL = list(Dictionary.values())
    try:
        num_tot = 0.0
        for item in all_vL:
            num_tot = num_tot + item
        final = float(num_tot / len(all_vL))
    except:  # XXX what kind of exception here?
        log.info("ALERT: Cannot average dictionary, passing first value instead.")
        final = all_vL[0]

    return final


# ----------------------------------------------------------------------------------------------------------------------


def binmode(data, bins=None):
    """Compute statistical mode of values in input data array

    Parameters
    ----------
    data : numpy array
        input values used to compute statistical mode

    bins : int
        number of bins to use for histogram generation. If not explicitly specified, the default value is 'None'.

    Returns
    -------
    mbin : float
        midpoint of the histogram bin with the largest value

    mbins : array
        array of bin edge values
    """
    from scipy import array, where, isinf, isnan, sort, zeros
    from scipy import arange, histogram, stats
    data = array(data)
    mdx = where(~isnan(data))
    data = data[mdx]
    mmdx = where(~isinf(data))
    data = data[mmdx]
    if bins != None:
        m, mbins = histogram(data, bins=bins)
    else:
        step = 1 / 100.
        splits = arange(0, 1 + step, step)
        bin_edges = stats.mstats.mquantiles(data, splits)
        bins = sort(list(set(bin_edges)))
        rebins = arange(min(bins[1:]), max(bins[:-1]), (max(bins[:-1]) - min(bins[1:])) * step)
        if len(rebins) > 0:
            m, mbins = histogram(data, bins=rebins)
        else:
            m = zeros(1);
            mbins = zeros(2)
    mdx = where(m == max(m))
    mbin = 0.5 * (mbins[mdx[0][0] + 1] + mbins[mdx[0][0]])
    return (mbin, mbins)


# ----------------------------------------------------------------------------------------------------------------------


def compute_abmag_zeropoint(imgname,inst_det):
    """Compute photometric AB mag zeropoint

    Parameters
    ----------
    imgname : string
        Name of the image to use in calculations

    inst_det : string
        Space-separated text string containing instrument name, detector name (upper case) of the products being
        processed.(i.e. WFC3_UVIS)

    Returns
    -------
    zpt_value : float
        AB magnitude photometric zeropoint value
    """

    if inst_det.lower().startswith('acs'):
        exten=1
    if inst_det.lower().startswith('wfc3'):
        exten=0
    if inst_det.lower().startswith('wfpc2'):
        exten=1
    try:
        photFlam = float(fits.getheader(ingname, exten)['PHOTFLAM'])
        photPlam = float(fits.getheader(ingname, exten)['PHOTPLAM'])
    except:
        exten = 1
        photFlam = float(fits.getheader(imgname, exten)['PHOTFLAM'])
        photPlam = float(fits.getheader(imgname, exten)['PHOTPLAM'])

    ref_stmag_zpt = -2.5 * (numpy.log10(photFlam)) - 21.10
    zpt_value = ref_stmag_zpt - (5 * numpy.log10(photPlam)) + 18.6921

    return(zpt_value)


# ----------------------------------------------------------------------------------------------------------------------


def compute_background (image, threshold=None):
    bkg_estimator = SExtractorBackground()
    bkgrms_estimator = StdBackgroundRMS()
    bkg = None
    bkg_dao_rms = None

    exclude_percentiles = [10, 25, 50, 75]
    for percentile in exclude_percentiles:
        log.info("Percentile in use: {}".format(percentile))
        try:
            bkg = Background2D(image, (50, 50),
                               filter_size=(3, 3),
                               bkg_estimator=bkg_estimator,
                               bkgrms_estimator=bkgrms_estimator,
                               exclude_percentile=percentile,
                               edge_method='pad')
            print('bkg: ', bkg)
        except Exception:
            bkg = None
            continue

        if bkg is not None:
            # If it succeeds, stop and use that value
            bkg_rms = (5. * bkg.background_rms)
            bkg_rms_mean = bkg.background.mean() + 5. * bkg_rms.std()
            default_threshold = bkg.background + bkg_rms
            bkg_dao_rms = bkg.background_rms
            if threshold is None:
                threshold = default_threshold
            elif threshold < 0:
                threshold = -1 * threshold * default_threshold
                log.info("{} based on {}".format(threshold.max(), default_threshold.max()))
                bkg_rms_mean = threshold.max()
            else:
                bkg_rms_mean = 3. * threshold

            if bkg_rms_mean < 0:
                bkg_rms_mean = 0.
            break

    # If Background2D does not work at all, define default scalar values for
    # the background to be used in source identification
    if bkg is None:
        bkg_rms_mean = max(0.01, image.min())
        bkg_rms = bkg_rms_mean * 5
        bkg_dao_rms = bkg_rms_mean

    return bkg, bkg_rms, bkg_dao_rms, threshold


# ----------------------------------------------------------------------------------------------------------------------



def conv_nan_zero(img_arr, replace_val=0.0, reverse=False):
    """Replace NaNs in an image arr with zeros

    Parameters
    ----------
    img_arr : numpy array of floating-point values
        image array to process

    replace_val : float
        replacement value. If not explicitly specified, the default value is '0.0'.

    reverse : Boolean
        perform the reverse operation instead (replace zeros with NaNs) (True/False)? If not explicitly specified, the
        default value is 'False'.

    Returns
    -------
    out_arr : numpy array
        de-NaNed (or de-zeroed if reverse = True) version of input array.
    """
    import numpy

    if reverse:
        Zeros = numpy.where(img_arr == replace_val)
        out_arr = numpy.copy(img_arr)
        out_arr[Zeros] = numpy.nan
    else:
        NaNs = numpy.where(img_arr != img_arr)
        out_arr = numpy.copy(img_arr)
        out_arr[NaNs] = replace_val
    return (out_arr)

# ----------------------------------------------------------------------------------------------------------------------


def create_dao_like_coordlists(totdet_product_cat_dict,filter_product_cat_dict,inst_det,param_dict,rms_dict):
    """Make daophot-like coordinate lists

    Parameters
    ----------
    totdet_product_cat_dict : dictionary
        Dictionary mapping image filename to corresponding catalog filename for total detection product

    filter_product_cat_dict : dictionary
        Dictionary mapping image filename to corresponding catalog filename for total detection product

    inst_det : string
        Space-separated text string containing instrument name, detector name (upper case) of the products being
        processed.(i.e. WFC3_UVIS)

    param_dict : dictionary
        dictionary of drizzle, source finding and photometric parameters

    rms_dict : dictionary
        dictionary containing 2d rms frames keyed by image name

    Returns
    -------
    daofind_white_sources : string
        Coordinate list filename
    """
    log.info("DAOPHOT-LIKE SOURCELIST CREATION OCCURS HERE!")

    totfiltprod_filename_list = list(filter_product_cat_dict.keys())
    tdp_imagename = list(totdet_product_cat_dict.keys())[0]
    tdp_catname = list(totdet_product_cat_dict.values())[0]
    # ### (1) ### Collect applicable parameters
    log.info("### (1) ### Collect applicable parameters")


    # ----------------------------------------
    # Calculate mean readnoise value (float):
    # ----------------------------------------
    readnoise_dictionary_drzs = get_readnoise(totfiltprod_filename_list)

    # -----------------------------
    # Get scale arcseconds / pixel
    # -----------------------------
    scale_dict_drzs = stwcs_get_scale(totfiltprod_filename_list)

    for img in list(readnoise_dictionary_drzs.keys()):
        log.info("{} {} {}".format(img,readnoise_dictionary_drzs[img],scale_dict_drzs[img]))


    # ### (2) ###  Use daostarfinder to create a sourcelist from the total detection image.

    log.info('### (2) ###  Use daostarfinder to create a sourcelist from the total detection image {}'
             .format(tdp_imagename))

    exp_dictionary_scis = {tdp_imagename: fits.getval(tdp_imagename,keyword='exptime')} # TODO: Quick and dirty hack for testing. FIND BETTER WAY TO DO THIS BEFORE DEPLOYMENT

    daofind_white_sources = run_daofind(param_dict,
                                        whitelightimage = tdp_imagename,
                                        rms_array = rms_dict[tdp_imagename],
                                        readnoise_dictionary_drzs = readnoise_dictionary_drzs,
                                        scale_dict_drzs = scale_dict_drzs,
                                        exp_dictionary_scis = exp_dictionary_scis,
                                        working_dir = os.getcwd())

    # ~~~~~~~~~~~~~Bail out if daofind can't locate a single source~~~~~~~~~~~~~~
    with open(daofind_white_sources) as f:
        sl_lines = f.readlines()
        if len(sl_lines) <= 4:
            log.info("*** WARNING: DAOFIND was unable to locate any sources in the detection image. No _daophot.txt sourcelist will be produced. ***")
            return(None)

    # ### (3) ###  Extract sources that fall "close" to 'INDEF' regions.
    # Take out any sources from the white-light source list falling within 'remove_radius' of a flag.
    log.info("\n(3) phot")
    flag_white = replace_NaN_w_flag_whitelight(tdp_imagename, rms_dict[tdp_imagename],os.getcwd())
    dict_source_lists_filtered = {}
    for drzimage in totfiltprod_filename_list:
        print("DRZIMAGE: ",drzimage)
        drzrootname = drzimage.replace(".fits", "")
        flag_name = drzrootname+"_01_sexflag.fits"
        flag_image,rms_data = replace_NaN_w_flag_whitelight(drzimage,
                                                   rms_dict[drzimage],
                                                   os.getcwd(),sci_ind=1,
                                                   flag_name = flag_name,
                                                   root_n = drzrootname + "_")
        rms_dict[drzimage] = rms_data
        daofind_white_open = open(daofind_white_sources)
        daofind_white_lines = daofind_white_open.readlines()
        daofind_white_open.close()
        #sci_sources = find_unique_name(extract_name(flag_image) + ".coo", os.path.dirname(flag_white), 'yes')
        sci_sources = unique_name(extract_name(flag_name) + "A.coo", suffix = ".coo")

        sci_sources = os.path.join(os.path.dirname(flag_name), sci_sources)
        # -------------------------------------
        #        save(sci_sources, daofind_white_lines)
        newfile = open(sci_sources, "w")
        rows = len(daofind_white_lines)
        for row in range(0, rows):
            newfile.write(daofind_white_lines[row])
        newfile.close()
        # -------------------------------------

        dict_source_lists_filtered[drzimage] = sci_sources


    return(dict_source_lists_filtered)


# ----------------------------------------------------------------------------------------------------------------------


def create_dao_like_sourcelists(dict_source_lists_filtered,filter_sorted_flt_dict,img_name,detection_image,inst_det,
                                param_dict,rms_dict):
    """Make source extractor-like sourcelists

    Parameters
    ----------
    dict_source_lists_filtered : dictionary
        Dictionary containing source extractor sourcelist file path keyed by total_drz.fits image

    filter_sorted_flt_dict : dictionary
        Dictionary of input flc/flt images used to create total drizzle-combined filter image 'img_name' sorted by
        filter name (lowercase).

    img_name : string
        drizzled total filter image to be processed

    detection_image : string
        total detection product filename

    inst_det : string
        Space-separated text string containing instrument name, detector name (upper case) of the products being
        processed.(i.e. WFC3_UVIS)

    param_dict : dictionary
        Dictionary of instrument/detector - specific drizzle, source finding and photometric parameters

    rms_dict : dictionary
        dictionary of rms images (numpy.ndarrays) keyed by the imagename listed in img_list

    Returns
    -------
    """
    # ### 4a ### compute daophot_process inputs
    log.info("DAOPHOT-LIKE SOURCELIST CREATION OCCURS HERE!")

    readnoise_dict = get_readnoise([img_name])
    log.info("{} readnoise: {}".format(img_name, readnoise_dict[img_name]))

    scale_dict = {}
    scale_dict[img_name] = param_dict['astrodrizzle']['SCALE']
    log.info("{} Scale: {}".format(img_name, scale_dict[img_name]))

    abmag_zpt_dict = {}
    abmag_zpt_dict[img_name] = compute_abmag_zeropoint(img_name, inst_det)
    log.info("{} AB magnitude zeropoint: {}".format(img_name, abmag_zpt_dict[img_name]))

    exptime_dict = {}
    exptime_dict[img_name] = fits.getval(img_name, keyword='exptime')
    log.info("{} Exposure time: {}".format(img_name, exptime_dict[img_name]))


    # ### (4) ### Feed corrected whitelight source lists into daophot with science images
    dict_newTAB_matched2drz = daophot_process([img_name],
                                 dict_source_lists_filtered,
                                 param_dict,
                                 readnoise_dict,
                                 scale_dict,
                                 abmag_zpt_dict,
                                 exptime_dict,
                                 os.getcwd(),
                                 rms_dict)

    # ### (5) ### Gather columns and put in nice format (dictated by: "column_keys_phot.cfg")
    # This will convert columns from xy to ra and dec (controlled by: "column_keys_phot.cfg")
    log.info("dict_newTAB_matched2drz: {}".format(dict_newTAB_matched2drz))

    hla_flag_filter.run_source_list_flaging([img_name], os.getcwd(), filter_sorted_flt_dict, param_dict,
                                            readnoise_dict, scale_dict, abmag_zpt_dict, exptime_dict, detection_image,
                                            dict_newTAB_matched2drz, 'daophot', os.getcwd(), rms_dict)
# ----------------------------------------------------------------------------------------------------------------------

def Create_MedDivImage(whitelightimage):
    """Computes a median-divided image.

    Parameters
    ----------

    whitelightimage : string
        Name of the drizzled image to process

    Returns
    -------
    medDivImg : basestring
        Name of the newly created median-divided image

    wht_data : numpy array
        WHT extension of the whitelightimage image *medDivImg*
    """
    from scipy import signal

    # ----------------------------
    # Create Median-Divided Image
    # ----------------------------
    medImg = whitelightimage + '_med.fits'
    medDivImg = whitelightimage + '_med_div.fits'

    white_light_hdu = fits.open(whitelightimage)

    wl_data = white_light_hdu[1].data
    wht_data = white_light_hdu[2].data

    wht_data = conv_nan_zero(wht_data)
    wl_sm = signal.medfilt2d(wl_data, 13)

    wl_med_val = numpy.median(wl_data[wl_data > 0.0])

    wl_sm_norm = wl_sm / wl_med_val

    wl_sm_norm = numpy.maximum(wl_sm_norm, 1.0)

    med_div = wl_data / wl_sm_norm

    med_div[numpy.where(wht_data <= 0.)] = -1

    pri_hdr = fits.getheader(whitelightimage, 0)
    fits.append(medDivImg, numpy.float32(med_div), pri_hdr)

    return medDivImg, wht_data


# ----------------------------------------------------------------------------------------------------------------------


def create_rms_image(img_list,thresh):
    """Creates RMS image for use by source-finding code.

    Parameters
    ----------
    img_list : list
        list of images to process

    thresh : float
        source extractor 'thresh' value

    Returns
    -------
    rms_dict : dictionary
        dictionary of rms images (numpy.ndarrays) keyed by the imagename listed in img_list
    """
    rms_dict = {}

    for imgname in img_list:
        img_data = fits.getdata(imgname, 1)
        bkg, bkg_rms, rms_array, threshold = compute_background(img_data,thresh)
        rms_dict[imgname] = rms_array

    return(rms_dict)


# ----------------------------------------------------------------------------------------------------------------------


def create_se_like_coordlists():
    """Make source extractor-like coordinate lists

    Parameters
    ----------

    Returns
    -------
    """

    log.info("SOURCE EXTRACTOR-LIKE COORDINATE LIST CREATION OCCURS HERE!")


# ----------------------------------------------------------------------------------------------------------------------


def create_se_like_sourcelists():
    """Make source extractor-like sourcelists

    Parameters
    ----------

    Returns
    -------
    """

    log.info("SOURCE EXTRACTOR-LIKE SOURCELIST CREATION OCCURS HERE!")


# ----------------------------------------------------------------------------------------------------------------------


def create_sourcelists(obs_info_dict, param_dict):
    """Main calling code. Make sourcelists

    Parameters
    ----------
    obs_info_dict : dictionary
        Dictionary containing all information about the images being processed

    param_dict : dictionary
        Dictionary of instrument/detector - specific drizzle, source finding and photometric parameters

    Returns
    -------
    """
    log.info("----------------------------------------------------------------------------------------------------------------------")
    for key1 in list(obs_info_dict.keys()):
        for key2 in list(obs_info_dict[key1].keys()):
            log.info("obs_info_dict[{}][{}]: {}".format(key1,key2,obs_info_dict[key1][key2]))  # TODO: REMOVE THIS SECTION BEFORE ACTUAL USE

    log.info("----------------------------------------------------------------------------------------------------------------------")
    log.info("SOURCELIST CREATION OCCURS HERE!")

    for tdp_keyname in [oid_key for oid_key in list(obs_info_dict.keys()) if
                        oid_key.startswith('total detection product')]:  # loop over total filtered products
        log.info("=====> {} <======".format(tdp_keyname))
        # 0: Map image filename to correspoinding catalog filename for total detection product and the associated
        # filter products
        totdet_product_cat_dict = {}
        filter_product_cat_dict = {}
        filter_component_dict = {}
        totdet_product_cat_dict[obs_info_dict[tdp_keyname]['product filenames']['image']] = \
            obs_info_dict[tdp_keyname]['product filenames']['source catalog']
        for fp_keyname in obs_info_dict[tdp_keyname]['associated filter products']:
            filter_product_cat_dict[obs_info_dict[fp_keyname]['product filenames']['image']] = \
                obs_info_dict[fp_keyname]['product filenames']['source catalog']
            filter_component_dict[obs_info_dict['filter product 00']['info'].split()[-1].lower()] = \
                obs_info_dict[fp_keyname]['files']

        inst_det = "{} {}".format(obs_info_dict[tdp_keyname]['info'].split()[-2],
                                  obs_info_dict[tdp_keyname]['info'].split()[-1])
        # 1: Generate source extractor-like sourcelist(s)
        create_se_like_coordlists()

        # 2: Generate daophot-like sourcelist(s)
        drz_img_list = list(totdet_product_cat_dict.keys())+list(filter_product_cat_dict.keys())
        rms_dict = create_rms_image(drz_img_list,param_dict[inst_det]['sourcex']['thresh'])
        dict_source_lists_filtered = create_dao_like_coordlists(totdet_product_cat_dict,filter_product_cat_dict,inst_det,param_dict[inst_det],rms_dict)

        # 3: Generate daophot-like and source extractor-like sourcelists from coordinate lists for each filter
        # assocatied with the current total detection product
        for img_name in filter_product_cat_dict.keys():
            sourcelist_name = filter_product_cat_dict[img_name]


            create_se_like_sourcelists()

            if dict_source_lists_filtered != None:
                create_dao_like_sourcelists(dict_source_lists_filtered,filter_component_dict,img_name,
                                            obs_info_dict[tdp_keyname]['product filenames']['image'],inst_det,
                                            param_dict[inst_det],rms_dict)
            else:
                log.info("Empty coordinate file. DAO sourcelist {} NOT created.".format(sourcelist_name))


# ----------------------------------------------------------------------------------------------------------------------


def daophot_process(all_drizzled_filelist, dict_source_lists_filtered, param_dict, readnoise_dictionary_drzs,
                    scale_dict_drzs, zero_point_AB_dict, exp_dictionary_scis, working_dir, rms_dict, ext=1,
                    Verbose=True, WHT=False):
    """
    This task will run the photometric calculations on a list of drizzled images. The coordinates for the sources are 
    found in the 'dict_source_lists_filtered' dictionary. This dictionary has keys that match the images inside '
    all_drizzled_filelist', and values corresponding to the source list.

    Tested.

    Parameters
    ----------
    all_drizzled_filelist : list
        List of images to process
        
    dict_source_lists_filtered : dictionary
        Dictionary containing source extractor sourcelist file path keyed by total_drz.fits image
        
    param_dict : dictionary
        dictionary of drizzle, source finding and photometric parameters
        
    readnoise_dictionary_drzs : dictionary
        dictionary of readnoise values matched to the science images.
        
    scale_dict_drzs : dictionary
        dictionary of pixel scale arcsecs/pixel values matched to the science images.
    
    zero_point_AB_dict : dictionary
        dictionary of AB magnitude zeropoint as values matched to the science images.
    
    exp_dictionary_scis : dictionary
        dictionary of exposure time values matched to the science images.
        
    working_dir : string
        where to store temporary products.
        
    rms_dict : dictionary
        dictionary of RMS values

    config_file : string? 
        instrument-specific param file
    
    ext : int
        extension of the science image inside the drizzled stack of each science image inside 'all_drizzled_filelist'. 
        Default value is '1'.
    
    Verbose : Boolean
        Generate verbose output? Default value is 'True'.
        
    WHT : Boolean
        Process weight image? Default value is 'False'.
        
    Returns
    -------
    return_dict : dictionary
        dictionary with the drizzled images as values matched to the photometric catalog generated by 
        daophot_style_photometry().
    """

    log.info(' ')
    log.info('**********************************')
    log.info('DAOPHOT FUNCTION INPUT PARAMETERS:')
    log.info('**********************************')
    log.info('ALL DRIZZLED FILE LIST: {}'.format(all_drizzled_filelist))

    # Gather basics:
    fwhm = param_dict['dao']['TWEAK_FWHMPSF']
    log.info('FWHM = {}'.format(fwhm))
    thresh = param_dict['dao']['TWEAK_THRESHOLD']
    log.info('THRESH = {}'.format(thresh))
    aps_like = "%s,%s" % (param_dict['dao']['aperture_1'], param_dict['dao']['aperture_2'])
    log.info('APS_LIKE = {}'.format(aps_like))
    apertures_sort = [param_dict['dao']['aperture_1'], param_dict['dao']['aperture_2']]
    log.info('APERTURES_SORT = {}'.format(apertures_sort))
    apertures_sort.sort()  # Sort smallest to largest.
    log.info('APERTURES_SORT = {}'.format(apertures_sort))
    largest_ap = apertures_sort[-1]
    log.info('LARGEST_AP = {}'.format(largest_ap))

    # Format image specific parameters:
    return_dict = {}
    for image_drz in all_drizzled_filelist:
        log.info('IMAGE_DRZ = {}'.format(image_drz))
        readnoise = readnoise_dictionary_drzs[image_drz]
        log.info('READNOISE = {}'.format(readnoise))
        scale = scale_dict_drzs[image_drz]
        log.info('SCALE = {}'.format(scale))
        annulus = scale * 5.
        log.info('ANNULUS = {}'.format(annulus))
        dannulus = scale * 5.
        log.info('DANNULUS = {}'.format(dannulus))
        zeropt = zero_point_AB_dict[image_drz]
        log.info('ZEROPT = {}'.format(zeropt))
        exptime = exp_dictionary_scis[image_drz]
        log.info('EXPTIME = {}'.format(exptime))

        rms_array = rms_dict[image_drz]
        #rms_array = pyfits.getdata(single_rms, 0)
        rms_image_median = binmode(rms_array[numpy.isfinite(rms_array) & (rms_array > 0.0)])[0]
        log.info(' ')
        # log.info('single rms image = {}'.format(single_rms))
        log.info("Median from RMS image = {}".format(rms_image_median))
        log.info(' ')

        # Define image input
        image2run = image_drz + "[%s]" % (ext)
        log.info('IMAGE2RUN = {}'.format(image2run))
        coordinates_filtered = dict_source_lists_filtered[image_drz]
        log.info('COORDINATES_FILTERED = {}'.format(coordinates_filtered))

        if Verbose:
            log.info("     Performing daophot_style_photometry on sources in {}\n".format(extract_name(image2run)))

        # ------------
        # Run DAOPhot
        # ------------
        name_daoOUT = extract_name(image_drz)  # (Want name to look like: "HST_10048_a1_ACS_HRC_F344N_daophot_tmp.txt".)
        log.info('NAME_DAOOUT = {}'.format(name_daoOUT))
        name_daoOUT = name_daoOUT.replace(".fits", "")
        log.info('NAME_DAOOUT = {}'.format(name_daoOUT))
        name_daoOUT = name_daoOUT.replace("_drz", "")
        log.info('NAME_DAOOUT = {}'.format(name_daoOUT))

        name_daoOUT = find_unique_name(name_daoOUT + "_daophot_tmp.txt", working_dir, 'no',
                                              FITS_file_rootname=name_daoOUT, Suffix="daophot_tmp.txt")
        log.info('NAME_DAOOUT = {}'.format(name_daoOUT))

        output_dao = os.path.join(working_dir, name_daoOUT)
        log.info('OUTPUT_DAO = {}'.format(output_dao))
        log.info(' ')

        # -----------------------------------------
        # PLATE SCALES & APERTURES
        # -------------------------
        # UVIS --> 0.040 arcsec/pixel
        # APERTURES --> 0.05 arcsec & 0.15 arcsec
        #               1.25 pixels & 3.75 pixels
        #
        # IR --> 0.13 arcsec/pixel
        # APERTURES --> 0.15 arcsec & 0.45 arcsec
        #               1.15 pixels & 3.46 pixels
        #
        # CI = mag1(0.15 arcsec)-mag2(0.45 arcsec)
        # -----------------------------------------

        if WHT:
            wht_image2run = image_drz + "[2]"
            wht_name_daoOUT = string.split(name_daoOUT, '.')[0] + '_WHT.txt'
            wht_output_dao = os.path.join(working_dir, wht_name_daoOUT)
            log.info("WEIGHT EXT")
            log.info('image = {}'.format(wht_image2run))
            log.info('coords = {}'.format(coordinates_filtered))
            log.info('output = {}'.format(wht_output_dao))
            daophot_style_photometry(wht_image2run, None, coordinates_filtered, wht_output_dao, scale, apertures_sort,
                                     annulus, dannulus, 'mode', exptime, zeropt, param_dict)
            return_dict[image_drz] = wht_output_dao
        else:
            log.info("SCI EXT")
            log.info('image = {}'.format(image2run))
            log.info('coords = {}'.format(coordinates_filtered))
            log.info('output = {}'.format(output_dao))
            output_dao = output_dao.replace("daophot_tmp.txt", "daophot.txt")
            daophot_style_photometry(image2run, None, coordinates_filtered, output_dao, scale, apertures_sort, annulus,
                                     dannulus, 'mode', exptime, zeropt, param_dict)
            return_dict[image_drz] = output_dao
        log.info('return_dict: {}'.format(return_dict))
    return return_dict


# ----------------------------------------------------------------------------------------------------------------------


def daophot_style_photometry(imgFile, errFile, cooFile, outFile, platescale, radiiArcsec, skyAnnulus, dSkyAnnulus,
                             salgorithm, gain, zeroPoint, param_dict):
    """generates iraf.daophot-like photometric sourcelist using spacetelescope.wfc3_photometry and astropy.photutils

    Parameters
    ----------
    imgFile : string
        name of the fits file that was used to generate the x,y coordinates stored in cooFile.

    errFile : string
        name of the fits file that contains the corresponding RMS map for imgFile

    cooFile : string
        name of the file that contains the x,y coordinates for photometry

    outFile : string
        name of the output sourcelist file

    platescale : float
        instrument platescale in arcseconds per pixel.

    radiiArcsec : list of floats
        list of photometric aperture radii (in units of arcseconds) to use

    skyAnnulus : float
        inner radius (in arcseconds) of the annulus that will be used to compute the background sky value.

    dskyannulus : float
        width (in arcseconds) of the annulus that will be used to compute the background sky value.

    salgorithm : string
        the statistic used to calculate the background. Choices are 'mean', 'median', or 'mode'. All measurements are
        sigma clipped. NOTE: From DAOPHOT, mode = 3 * median - 2 * mean.

    gain : float
        gain in electrons per adu (only use if image units aren't e-).

    zeroPoint : float
        photometric zeropoint used to compute magnitude values from flux values

    param_dict : dictionary
        dictionary of drizzle, source finding and photometric parameters

    Returns
    -------
    Nothing!
    """
    # convert input values whose units are arcseconds to pixles
    radii = []
    for ctr in range(0, len(radiiArcsec)):
        radii.append(radiiArcsec[ctr] / platescale)  # radii.append(round(radiiArcsec[ctr] / platescale))
    skyAnnulus /= platescale
    dSkyAnnulus /= platescale

    verbose = True
    if verbose:
        log.info("SUMMARY OF INPUT PARAMETERS")
        log.info("imgFile:          {}".format(imgFile))
        log.info("errFile:          {}".format(errFile))
        log.info("cooFile:          {}".format(cooFile))
        log.info("outFile:          {}".format(outFile))
        log.info("platescale:       {}".format(platescale))
        log.info("radii (pixels):   {}".format(radii))
        log.info("radii (arcsec):   {}".format(radiiArcsec))
        log.info("annulus:          {}".format(skyAnnulus))
        log.info("dSkyAnnulus:      {}".format(dSkyAnnulus))
        log.info("salgorithm:       {}".format(salgorithm))
        log.info("gain:             {}".format(gain))
        log.info("zeropoint:        {}".format(zeroPoint))

    # read in data from  image, err and coo files
    # read in image data
    parse_imgFile = imgFile.split("[")
    imgFile = parse_imgFile[0]
    if len(parse_imgFile) == 2:
        fitsExt = int(parse_imgFile[1].replace("]", ""))
        imgData = fits.getdata(imgFile, ext=fitsExt)
    if len(parse_imgFile) == 2:
        imgData = fits.getdata(imgFile)

    # read in rms data
    if errFile:
        parse_errFile = errFile.split("[")
        errFile = parse_errFile[0]
        if len(parse_errFile) == 2:
            fitsExt = int(parse_errFile[1].replace("]", ""))
            errData = fits.getdata(errFile, ext=fitsExt)
        if len(parse_errFile) == 1:
            errData = fits.getdata(errFile)
    if not errFile:
        errData = None

    # read in coordinate data
    x, y = numpy.loadtxt(cooFile, skiprows=4, usecols=(0, 1), unpack=True)

    # Do photometry
    # adjust coods for calculations that assume origin value of 0, rather than 1.
    x = x - 1.
    y = y - 1.

    bgAps = CircularAnnulus((x, y), r_in=skyAnnulus, r_out=skyAnnulus + dSkyAnnulus)  # compute background

    photAps = [CircularAperture((x, y), r=r) for r in radii]
    photometry_tbl = iraf_style_photometry(photAps, bgAps, data=imgData, platescale=platescale, error_array=errData,
                                           bg_method=salgorithm, epadu=gain, zero_point=zeroPoint)

    # convert coords back to origin value = 1 rather than 0
    photometry_tbl["XCENTER"] = photometry_tbl["XCENTER"] + 1.
    photometry_tbl["YCENTER"] = photometry_tbl["YCENTER"] + 1.

    # calculate and add RA and DEC columns to table
    ra, dec = Transform_list_xy_to_RA_Dec(photometry_tbl["XCENTER"], photometry_tbl["YCENTER"], imgFile)
    raCol = Column(name="RA", data=ra, dtype=numpy.float64)
    DecCol = Column(name="DEC", data=dec, dtype=numpy.float64)
    photometry_tbl.add_column(raCol, index=2)
    photometry_tbl.add_column(DecCol, index=3)

    # TODO: Execute Daophot_ap_correct() here!
    photometry_tbl.write("photometry_tbl.csv", format='ascii.csv', overwrite=True)
    log.info("WROTE photometry_tbl.csv !")
    # Calculate and add concentration index (CI) column to table
    ci_data = photometry_tbl["MAG_{}".format(radiiArcsec[0])].data - photometry_tbl[
        "MAG_{}".format(radiiArcsec[1])].data
    ciMask = numpy.logical_and(numpy.abs(ci_data) > 0.0, numpy.abs(ci_data) < 1.0e-30)
    bigBadIndex = numpy.where(abs(ci_data) > 1.0e20)
    ciMask[bigBadIndex] = True
    ciCol = MaskedColumn(name="CI", data=ci_data, dtype=numpy.float64, mask=ciMask)
    photometry_tbl.add_column(ciCol)

    # Add zero-value "Flags" column in preparation for source flagging
    flagCol = Column(name="Flags", data=numpy.zeros_like(photometry_tbl['ID']), dtype=numpy.int64)
    photometry_tbl.add_column(flagCol)

    # Add null-value "TotMag(<outer radiiArc>)" and "TotMag(<outer radiiArc>)" columns
    emptyTotMag = MaskedColumn(name="TotMag({})".format(radiiArcsec[1]), fill_value=None, mask=True,
                               length=len(photometry_tbl["XCENTER"].data), dtype=numpy.int64)
    emptyTotMagErr = MaskedColumn(name="TotMagErr({})".format(radiiArcsec[1]), fill_value=None, mask=True,
                                  length=len(photometry_tbl["XCENTER"].data), dtype=numpy.int64)
    photometry_tbl.add_column(emptyTotMag)
    photometry_tbl.add_column(emptyTotMagErr)

    # build final output table
    finalColOrder = ["XCENTER", "YCENTER", "RA", "DEC", "ID", "MAG_{}".format(radiiArcsec[0]),
                     "MAG_{}".format(radiiArcsec[1]), "MERR_{}".format(radiiArcsec[0]),
                     "MERR_{}".format(radiiArcsec[1]), "MSKY", "STDEV", "FLUX_{}".format(radiiArcsec[1]),
                     "TotMag({})".format(radiiArcsec[1]), "TotMagErr({})".format(radiiArcsec[1]), "CI", "Flags"]
    outputPhotometryTable = photometry_tbl[finalColOrder]

    # format output table columns
    finalColFormat = {"RA": "13.10f", "DEC": "13.10f", "MAG_{}".format(radiiArcsec[0]): '6.3f',
                      "MAG_{}".format(radiiArcsec[1]): '6.3f', "MERR_{}".format(radiiArcsec[0]): '6.3f',
                      "MERR_{}".format(radiiArcsec[1]): '6.3f', "MSKY": '10.8f', "STDEV": '10.8f',
                      "FLUX_{}".format(radiiArcsec[1]): '10.8f', "CI": "7.3f"}
    for fcf_key in list(finalColFormat.keys()):
        outputPhotometryTable[fcf_key].format = finalColFormat[fcf_key]

    # change some column titles to match old daophot.txt files
    rename_dict = {"XCENTER": "X-Center", "YCENTER": "Y-Center",
                   "MAG_{}".format(radiiArcsec[0]): "MagAp({})".format(radiiArcsec[0]),
                   "MAG_{}".format(radiiArcsec[1]): "MagAp({})".format(radiiArcsec[1]),
                   "MERR_{}".format(radiiArcsec[0]): "MagErr({})".format(radiiArcsec[0]),
                   "MERR_{}".format(radiiArcsec[1]): "MagErr({})".format(radiiArcsec[1]),
                   "MSKY": "MSky({})".format(radiiArcsec[1]), "STDEV": "Stdev({})".format(radiiArcsec[1]),
                   "FLUX_{}".format(radiiArcsec[1]): "Flux({})".format(radiiArcsec[1])}
    for oldColTitle in rename_dict:
        outputPhotometryTable.rename_column(oldColTitle, rename_dict[oldColTitle])
        log.info("Column '{}' renamed '{}'".format(oldColTitle, rename_dict[oldColTitle]))

    outputPhotometryTable.write(outFile, format='ascii.csv', overwrite=True)

    log.info("Wrote {}".format(outFile))

    add_header_phot_tab(outFile, imgFile, param_dict)


# ----------------------------------------------------------------------------------------------------------------------


def extract_name(stringWpath):
    """
    This task will extract just the name of a specific filename that includes the path in the name: 'stringWpath'.

    Tested.

    Parameters
    ----------
    stringWpath : string
        filename with full path

    Returns
    -------
    stringname : string
        naked filename stripped of its path
    """
    while "/" == stringWpath[-1]:
        stringWpath = stringWpath[:-1]
    stringname = stringWpath.split("/")[-1]

    return stringname


# ----------------------------------------------------------------------------------------------------------------------


def filter_daolist(infile, outfile, mask_image, edgemask=0):
    """Read DAOfind source list and filter out images beyond edge of mask. Returns new count of sources in file

    Parameters
    ----------
    infile : string
        Name of .coo input file to read

    outfile : string
        Name of .coo filtered output file to write (may be same as infile to replace infile)

    mask_image : numpy.ndarray
        boolean image with true values in pixels to reject

    edgemask : int
        expand masked region by edgemask pixels

    Returns
    -------
    count : int
        Updated count of the number of sources in file
    """

    ny, nx = mask_image.shape
    if edgemask:
        # use binary dilation to expand mask
        ky, kx = numpy.ogrid[-edgemask:edgemask + 1, -edgemask:edgemask + 1]
        kernel = kx * kx + ky * ky <= edgemask * edgemask
        mask_image = scipy.ndimage.morphology.binary_dilation(mask_image, structure=kernel)

    src_infile = open(infile, 'r')
    src_lines = src_infile.readlines()
    src_infile.close()

    if outfile == infile:
        mod_outfile = infile + ".mod"
    else:
        mod_outfile = outfile
    src_outfile = open(mod_outfile, 'w')

    count = 0
    for i, src_line in enumerate(src_lines):
        if src_line.startswith('#'):
            src_outfile.write(src_line)
        else:
            src_line_split = src_line.strip().split()
            if len(src_line_split) > 3:
                mag_value = src_line_split[2]
                if mag_value != "INDEF":
                    x_cen = int(float(src_line_split[0])) - 1
                    y_cen = int(float(src_line_split[1])) - 1
                    if x_cen >= 0 and x_cen < nx and y_cen >= 0 and y_cen < ny and not mask_image[y_cen, x_cen]:
                        src_outfile.write(src_line)
                        count += 1
    src_outfile.close()
    if outfile != mod_outfile:
        os.rename(mod_outfile, outfile)
    return count


# ----------------------------------------------------------------------------------------------------------------------


def find_unique_name(beginname, dir, isdir, FITS_file_rootname=False, Suffix=False):
    """Generates a unique filename based on the input rootname/directory name/filename and path.

    Parameters
    ----------
    beginname: string
        input rootname/directory name/filename

    dir : string
        path used for the unique filename search

    isdir : string
        is **beginname** a directory (yes/no)?

    FITS_file_rootname : string/boolean
        If specified, this value will interpreted as a fits file rootname. This rootname will be used as basis for the
        output filename. Default value is Boolean 'False'.

    Suffix : string/boolean
        Output filename suffix. Default value is Boolean 'False'.

    Returns
    -------
    filename : string
        unique filename
    """
    directory4file = os.listdir(dir)

    if FITS_file_rootname:
        filename = beginname
        p = 1
        while filename in directory4file:
            p = p + 1
            k = str(p)
            if p < 10:
                filename = FITS_file_rootname + "0" + k + "_" + Suffix
            else:
                filename = FITS_file_rootname + k + "_" + Suffix
    elif isdir == "yes":
        filename = beginname + "_1"
        p = 1
        while filename in directory4file:
            p = p + 1
            k = "_" + str(p)
            filename = beginname + k
    else:
        filename = beginname + "_1.txt"
        p = 1
        while filename in directory4file:
            p = p + 1
            k = "_" + str(p) + ".txt"
            filename = beginname + k

    return filename

# ----------------------------------------------------------------------------------------------------------------------


def get_mean_readnoise(image):
    """This subroutine computes mean readnoise values

    Parameters
    ----------
    image : string
        image filename

    Returns
    -------
    readnoise : float
        average readnoise value
    """
    imghdu = fits.open(image)
    header = imghdu[0].header
    filename = image.split('/')[-1]
    if filename.startswith('hst_'):
        inst = filename.split('_')[3]
        detect = filename.split('_')[4]
    else:
        inst = header['instrume'].lower()
        if inst != 'wfpc2': detect = header['DETECTOR']
        if inst == 'wfpc2':
            rn = header['ROOTNAME'].lower()
            rn = rn.split("/")[-1].split("_")[0]
            if rn[len(rn) - 1] != 'x': detect = 'wfpc2'
            if rn[len(rn) - 1] == 'x': detect = 'pc'

    if inst != 'wfpc2':
        readnsea = header['READNSEA']
        readnseb = header['READNSEB']
        readnsec = header['READNSEC']
        readnsed = header['READNSED']
        readnoise = numpy.mean([readnsea, readnseb, readnsec, readnsed])
    else:
        gain_val = header['ATODGAIN']
        # below values taken from WFPC2 instrument handbook, Table 4.2, Page 81
        if gain_val == 15.0:
            readnsea = 7.02  # PC1
            readnseb = 7.84  # WF2
            readnsec = 6.99  # WF3
            readnsed = 8.32  # WF4
        else:
            # gain = 7 (default)
            readnsea = 5.24  # PC1
            readnseb = 5.51  # WF2
            readnsec = 5.22  # WF3
            readnsed = 5.19  # WF4
        if detect == 'pc':
            readnoise = readnsea
        else:
            readnoise = numpy.mean([readnsea, readnseb, readnsec, readnsed])
    return (readnoise)


# ----------------------------------------------------------------------------------------------------------------------


def get_head_val_opened_fits(loaded_fits, key_word, ext=0):
    """
    Return a specific header value for an image loaded already with pyfits.
    Used in the header creation of the DAOPhot catalog.

    Tested.


    Parameters
    ----------
    loaded_fits : HDUList object
        A previously opened fits file already in memory

    key_word : string
        Name of the header field to be returned

    ext : int
        extension of **loaded_fits** to query. If not specified, the default value is 0, the primary header.

    Returns
    -------
    hd_value : type varies
        If found, the header value specified by input parameters will be returned. If the header keyword specified
        can't be found, text string 'not available' will be returned.
    """
    # if loaded_fits[ext].header.has_key(key_word): .has_key depricated
    #     hd_value = loaded_fits[ext].header[key_word]
    # else:
    #     hd_value = "not available"

    hval_check = key_word in loaded_fits[ext].header
    if hval_check:
        hd_value = loaded_fits[ext].header[key_word]
    else:
        hd_value = "not available"

    return hd_value


# ----------------------------------------------------------------------------------------------------------------------


def get_readnoise(listofimages):
    """
    This task will grab the average readnoise for HST data from the header.  This task is known to call the correct
    header keys for ACS UVIS data as well as WFC3 UVIS and IR data, and WFPC2 data.

    Tested.

    Parameters
    ----------
    listofimages : list
        list of images that will be used to get readnoise values

    Returns
    -------
    dictionary_output : dictionary
        A dictionary of readnoise values keyed by image name
    """
    dictionary_output = {}
    for individual_image in listofimages:
        try:
            ref_readnoise = get_mean_readnoise(individual_image)
        except:  # XXX what kind of exception here?
            exc_type, exc_value, exc_tb = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)
            log.info("ALERT: Readnoise could not be gathered from the header for image {}".format(extract_name(individual_image)))
            ref_readnoise = 0.0
        dictionary_output[individual_image] = ref_readnoise
    return dictionary_output


# ----------------------------------------------------------------------------------------------------------------------


def replace_NaN_w_flag_whitelight(whitelight_im, white_rms_data, working_dir, sci_ind=1, rms_ind=0,
                                  flag_name="white_sexflag.fits", root_n="white"):
    """Take the white light image and the RMS image and replace all NaN values with 'replace_val', and output flag-map
    for Source Extractor built off of the 'OR' comparison of the flag-map for the white light image and the flag-map
    for the white-RMS:

    "External flags come from 'flag-maps': these are images with the same size as the one where objects are detected,
    where integer numbers can be used to flag some pixels" - Source Extractor manual, section 9.

    Tested.

    Parameters
    ----------
    whitelight_im : string
        White-light image filename

    white_rms_data : numpy.ndarray
        White-light RMS image filename

    working_dir : string
        Path where products are saved.

    sci_ind : int
        White-light image fits file extension to operate on. Default value is 0.

    rms_ind : int
        RMS image fits file extension to operate on.Default value is 0.

    flag_name : string
        Name of the flag-map output image produced by this subroutine. Default value is 'white_sexflag.fits'.

    root_n : string
        Fits file rootname to be used as basis for output filename. Default value is 'white'.

    Returns
    -------
    flag_imageWpath: string
        Full filename (with path) of the 'de-NaNed' image.
    """
    log.info(" Working on white light - replace_NaN_w_flag_whitelight 1 ")


    whitearray = fits.open(whitelight_im, mode='update')


    whitearray_data = whitearray[sci_ind].data


    # ----------------
    # Make flag array
    # ----------------
    new_flag_array = numpy.zeros(numpy.shape(whitearray_data), dtype=int)
    where_NaNs_r = numpy.where(whitearray_data != whitearray_data)
    new_flag_array[where_NaNs_r] = 1
    where_NaNs_r_in_rms = numpy.where(white_rms_data != white_rms_data)
    new_flag_array[where_NaNs_r_in_rms] = 1

    flag_name = find_unique_name(flag_name, working_dir, 'no', FITS_file_rootname=root_n, Suffix="sexflag.fits")
    #    flag_name = Rename.unique_name(flag_name, suffix = "sexflag.fits")

    flag_imageWpath = os.path.join(working_dir, flag_name)
    fits.writeto(flag_imageWpath, numpy.int32(new_flag_array))

    log.info(" Working on white light - replace_NaN_w_flag_whitelight 2 ")

    # -----------------------------
    # Replace NaN's in white-light
    # -----------------------------
    whitearray_data = conv_nan_zero(whitearray_data)
    whitearray[sci_ind].data = whitearray_data
    whitearray.flush()

    log.info(" Working on white light - replace_NaN_w_flag_whitelight 3 ")

    return flag_imageWpath,white_rms_data


# ----------------------------------------------------------------------------------------------------------------------


def Transform_list_xy_to_RA_Dec(list_of_x,list_of_y, drizzled_image):
    """Transform lists of X and Y coordinates to lists of RA and Dec coordinates

    Tested.

    list_of_x : list
        list of x coordinates to convert

    list_of_y : list
        list of y coordinates to convert

    drizzled_image : string
        Name of the image that corresponds to the table from DAOPhot. This image is used to re-write x and y coordinates in RA and Dec.

    Returns
    -------
    RA : list
        A list of right ascension values

    Dec : list
        A list declination values.
    """
    import stwcs

    wcs1_drz = stwcs.wcsutil.HSTWCS(drizzled_image + "[1]")
    origin = 1
    # *origin* is the coordinate in the upper left corner of the
    # image.  In FITS and Fortran standards, this is 1.  In Numpy and C
    # standards this is 0.
    try:
        skyposish = wcs1_drz.all_pix2sky(list_of_x, list_of_y, origin)
    except AttributeError:
        skyposish = wcs1_drz.all_pix2world(list_of_x, list_of_y, origin)
    RA = skyposish[0]
    Dec = skyposish[1]

    return RA,Dec


# ----------------------------------------------------------------------------------------------------------------------


@util.with_logging
def run_create_sourcelists(obs_info_dict, param_dict):
    """ subroutine to run create_sourcelists and produce log file when not run sourcelist_generation is not run from
    hlaprocessing.py.

    Parameters
    ----------
    obs_info_dict : dictionary
        Dictionary containing all information about the images being processed

    param_dict : dictionary
        Dictionary of instrument/detector - specific drizzle, source finding and photometric parameters

    Returns
    -------
    """

    create_sourcelists(obs_info_dict, param_dict)


# ----------------------------------------------------------------------------------------------------------------------


def run_daofind(param_dict, filelist=None, source_match=50000., verbose=True,whitelightimage=None, rms_array=None,
                readnoise_dictionary_drzs=None, scale_dict_drzs=None,exp_dictionary_scis=None, working_dir=None,
                sl_ext = 0,sharphi=None, sharplo=None,edgemask=5):
    """Generates sourcelists using DAOfind.

    Parameters
    ----------
    param_dict : dictionary
        dictionary of drizzle, source finding and photometric parameters

    X filelist : list
        Default value is 'None'.

    X source_match : float
        Source matching list length threshold. Default value is '50000.'.

    X verbose : Boolean
        Print verbose output? Default value is 'True'.

    whitelightimage : string
        Name of the multi-filter composite image produced by hla_reduction.py. Default value is 'None'.

    rms_array : numpy.ndarray
        Multi-filter composite RMS image. Default value is 'None'.

    readnoise_dictionary_drzs : Dictionary
        Dictionary containing readnoise values keyed by filter-specific drizzled image filename. Default value is
        'None'.

    X scale_dict_drzs : Dictionary
        **UNUSED** Dictionary containing scale values keyed by filter-specific drizzled image filename. Default value
        is 'None'.

    exp_dictionary_scis : Dictionary
        Dictionary containing exposure time values keyed by filter-specific drizzled image filename. Default value is '
        None'.

    working_dir : string
        Working directory. Default value is 'None'.

    sl_ext : integer
        FITS image extension to perform DAOfind on. Default value is '0'.

    sharphi : float
        Upper limit on SHARPNESS parameter. Default value is 'None'.

    sharplo : float
        Lower limit on SHARPNESS parameter. Default value is 'None'.

    edgemask : integer
        Remove sources within edgemask pixels of image edge. Default value is '5.'.

    Returns
    -------
    output_dao : string
        daofind output coordinate filename.
    """


    # set up default inputs for run_DAOStarFinder()
    daoParams_default = {}
    daoParams_default["fwhm"] = 2.5
    daoParams_default["scale"] = 1.0
    daoParams_default["threshold"] = 4.0
    daoParams_default["sigma"] = 0.0
    daoParams_default["ratio"] = 1.0
    daoParams_default["theta"] = 0.0
    daoParams_default["sharplo"] = 0.0
    daoParams_default["sharphi"] = 1.0
    daoParams_default["roundlo"] = -1.0
    daoParams_default["roundhi"] = 1.0

    daoParams=daoParams_default #set up dictionary and get default values.

    ### get daofind parameters
    fwhm = float(param_dict["dao"]["TWEAK_FWHMPSF"])
    thresh = float(param_dict["dao"]["TWEAK_THRESHOLD"])
    ap_diameter1 = float(param_dict["dao"]["aperture_1"])
    ap_diameter2 = float(param_dict["dao"]["aperture_2"])
    scale = float(param_dict['astrodrizzle']['SCALE'])

    daoParams["fwhm"] = fwhm
    daoParams["threshold"] = thresh
    daoParams["scale"] = scale
    daoParams["ratio"] = 0.8
    log.info(' ')
    log.info('run_daofind INPUT PARAMETERS:')
    log.info('-----------------------------')
    log.info('fwhm = {}'.format(fwhm))
    log.info('thresh = {}'.format(thresh))
    log.info('scale = {}'.format(scale))
    if sharphi:
        log.info('sharphi = {}'.format(sharphi))
        daoParams["sharphi"] = sharphi

    if sharplo:
        log.info('sharplo = ',sharplo)
        daoParams["sharplo"] = sharplo

    readnoise = average_values_from_dict(readnoise_dictionary_drzs)
    exptime = average_values_from_dict(exp_dictionary_scis)

    # ----------------------------
    # Create Median-Divided Image
    # ----------------------------
    medDivImg,wht_data = Create_MedDivImage(whitelightimage)



    # rms_array = fits.getdata(whitelightrms,0)
    rms_image_median = binmode(rms_array[numpy.isfinite(rms_array) & (rms_array > 0.0)])[0]
    log.info("Median from RMS image = {}".format(rms_image_median))
    daoParams["sigma"] = rms_image_median

    # log.info('white light rms image = {}'.format(whitelightrms))
    log.info('sigma = {}'.format(daoParams["sigma"]))
    log.info('readnoise = {}'.format(readnoise))
    log.info('exptime = {}'.format(exptime))
    log.info(' ')

    name_daoOUT = extract_name(whitelightimage)
    name_daoOUT = unique_name(name_daoOUT + ".coo", suffix=".coo")
    output_dao = os.path.join(working_dir, name_daoOUT)

    try:
        log.info("image = {}".format(medDivImg))
        log.info("output = {}".format(output_dao))

        run_DAOStarFinder(medDivImg,output_dao,daoParams,debug=False)
    except: #XXX what kind of exception here?
        log.info(' ')
        log.info('****************************************************************')
        log.info('WARNING: THE MEDIAN-DIVIDED IMAGE CONTAINS MULTIPLE EXTENSIONS; ')
        log.info('                  STARTING FLAG AND FILTER PROCESSING.          ')
        log.info('****************************************************************')
        log.info(' ')
        log.info("image = {}".format(medDivImg+"[1]"))
        log.info("output = {}".format(output_dao))
        run_DAOStarFinder(medDivImg + "[1]", output_dao, daoParams, debug=False)
    reject_image = numpy.zeros(wht_data.shape,dtype=numpy.int16)
    reject_image[numpy.where(wht_data <= 0.)] = 1

    rej_img = whitelightimage+'_rej_img.fits'
    fits.append(rej_img, numpy.float32(reject_image))

    mod_output_dao = output_dao+".mod"
    filter_daolist(output_dao, mod_output_dao, reject_image, edgemask)

    os.rename(output_dao, output_dao+".OLD ")
    os.rename(mod_output_dao, output_dao)

    return output_dao


# ----------------------------------------------------------------------------------------------------------------------


def run_DAOStarFinder(imgName,outFilename,daoParams,debug=False):
    """runs astropy.photutils.DAOStarFinder() to identify sources in *imgName*.

    Writes the following information to file *outFileName* for each detected source:
        * X centroid
        * Y centroid
        * Magnitude
        * Sharpness
        * S-Round
        * G-Round
        * Source ID number

    Parameters
    ----------
    imgName : string
        Input image name and (optionally) the extension  in square brackets. Can be in the form
        [<EXT_NAME>,<GROUP_NUM>] or simply [<EXT_NUM>].

    outFilename : string
        name of the file that information on the detected sources will be written to.

    daoParams : dictionary
        dictionary of floating point values passed to DAOStarFinder. Values are as follows:

        * fwhm: The full-width at half-maximum of the point spread function in scale units. Default value = 2.5
        * scale: Image scale in units per pixel. Used to convert *in_fwhm* value from arcseconds to pixels.
        Default value = 1.0
        * threshold: Threshold in sigma for feature detection. Default value = 4.0
        * sigma: Standard deviation of background in counts. Used to convert *in_threshold* to counts.
        Default value = 0.0
        * ratio: Ratio of minor to major axis of Gaussian kernel. Default value = 1.0
        * theta: Position angle of major axis of Gaussian kernel. Default value = 0.0
        * sharplo: Minimum bound on sharpness for feature detection. Default value = 0.0
        * sharphi: Maximum bound on sharpness for feature detection. Default value = 1.0
        * roundlo: Minimum bound on roundness for feature detection. Default value = -1.0
        * roundhi: Minimum bound on roundness for feature detection. Default value = 1.0

    debug : Boolean
        write out 'xcentroid' and 'ycentroid' values to separate file for troubleshooting purposes (True/False)?
        Default value is 'False'.

    Returns
    -------
    Nothing.
    """
    if imgName.endswith("]"):
        parse_imgname=imgName.split("fits[")
        imgName=parse_imgname[0]+"fits"
        imgHDU = fits.open(imgName)
        rawext=parse_imgname[1][:-1]
        parse_rawext=rawext.split(",")
        if len(parse_rawext) == 1:
            extnum = int(parse_rawext[0])
        else:
            extName=parse_rawext[0]
            extGroupNum=int(parse_rawext[1])
            extname_ctr=0
            for extctr in range(0, len(imgHDU)):
                if imgHDU[extctr].name.upper() == extName.upper(): extname_ctr += 1
                if ((imgHDU[extctr].name.upper() == extName.upper()) and (extname_ctr == extGroupNum)):
                    extnum = extctr
                    break
    else:
        imgHDU = fits.open(imgName)
        extnum=0

    imgData = imgHDU[extnum].data

    mean, median, std = sigma_clipped_stats(imgData, sigma=3.0, iters=5)

    daofind = findstars.DAOStarFinder(fwhm=daoParams["fwhm"] / daoParams["scale"], threshold=daoParams["threshold"] * daoParams["sigma"], ratio=daoParams["ratio"], theta=daoParams["theta"], sharplo=daoParams["sharplo"], sharphi=daoParams["sharphi"], roundlo=daoParams["roundlo"], roundhi=daoParams["roundhi"])
    sources = daofind(imgData - median)

    if os.path.exists(outFilename):
        cmd="rm -f "+outFilename
        log.info(cmd)
        os.system(cmd)
    fout = open(outFilename, 'w')
    fout.write("#N XCENTER   YCENTER   MAG      SHARPNESS   SROUND      GROUND      ID         \ \n")
    fout.write("#U pixels    pixels    #        #           #           #           #          \ \n")
    fout.write("#F %-13.3f   %-10.3f   %-9.3f   %-12.3f     %-12.3f     %-12.3f     %-6d       \ \n")
    fout.write("#\n")
    for line in sources:
        fout.write("   %-13.3f%-10.3f%-9.3f%-12.3f%-12.3f%-12.3f%-6d\n"%(line["xcentroid"]+1.0, line["ycentroid"]+1.0,line["mag"],line["sharpness"],line["roundness1"],line["roundness2"],line["id"]))
    fout.close()
    log.info("Wrote {}".format(outFilename))

    if debug:
        outfilename=outFilename+".xy"
        if os.path.exists(outfilename):
            cmd="rm -f "+outfilename
            log.info(cmd)
            os.system(cmd)
        fout = open(outfilename, 'w')
        for line in sources:
            fout.write("%f %f\n"%(line["xcentroid"]+1.0, line["ycentroid"]+1.0)) #Write only X and Y coords to file.
        fout.close()
        log.info("Wrote {}".format(outFilename))
    return()


# ----------------------------------------------------------------------------------------------------------------------


def stwcs_get_scale(listofimages):
    """
    This task will grab the arcsec/pixel scale for HST data from the WCS information of the header.

    Note: Assumes science image is extension: "[1]".

    Tested.

    Parameters
    ----------
    listofimages : list
        list of images that will be used to get scale values

    Returns
    -------
    dictionary_output : dictionary
        A dictionary of scale values keyed by image name
    """
    import stwcs
    from stwcs import wcsutil
    dictionary_output = {}
    for individual_image in listofimages:
        try:
            wcs1 = stwcs.wcsutil.HSTWCS(individual_image + "[1]")
            scale2return = wcs1.pscale
        except:  # XXX what kind of exception here?
            log.info("ALERT: Pixel-scale could not be gathered from the header for image {}."
                     .format(extract_name(individual_image)))
            scale2return = False
        dictionary_output[individual_image] = scale2return

    return dictionary_output


# ----------------------------------------------------------------------------------------------------------------------


def unique_name(path, suffix='', dbg=False):
    """If path name exists, append with count, or incremented count

    Parameters
    ----------
    path : string
        file path

    suffix : string
        file suffix

    dbg : Boolean
        Turn debug mode on? (True/False)

    Returns
    -------
    outstr : string
        a unique filename
    """
    check = True
    k = 0
    namestr = path.split('_')
    try:
        k = int(namestr[-1])
    except:
        k = int(1)
    if dbg: pdb.set_trace()
    while check:
        outstr = ''
        for ll in namestr[:-1]:
            outstr += ll + '_'
        outstr += str(k) + suffix
        if not os.path.exists(outstr):
            check = False
        else:
            k += 1
    if dbg: pdb.set_trace()
    return (outstr)

