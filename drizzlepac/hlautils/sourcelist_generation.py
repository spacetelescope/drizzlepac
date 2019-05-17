#!/usr/bin/env python

"""This script contains code to support creation of source extractor-like and daophot-like sourcelists.

"""
import os
import pdb
import sys
import traceback

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
import numpy
from photutils import detection, findstars
from stsci.tools import logutil

__taskname__ = 'sourcelist_generation'

log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout)

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
    import numpy
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


def create_daophot_like_sourcelists(totdet_product_cat_dict,filter_product_cat_dict,inst_det,param_dict):
    """Make daophot-like sourcelists

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

    Returns
    -------
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

    whitelightrms_image = create_rms_image()

    exp_dictionary_scis = {tdp_imagename: fits.getval(tdp_imagename,keyword='exptime')} # TODO: Quick and dirty hack for testing. FIND BETTER WAY TO DO THIS BEFORE DEPLOYMENT
    daofind_white_sources = run_daofind(param_dict,
                                        whitelightimage = tdp_imagename,
                                        whitelightrms = whitelightrms_image,
                                        readnoise_dictionary_drzs = readnoise_dictionary_drzs,
                                        scale_dict_drzs = scale_dict_drzs,
                                        exp_dictionary_scis = exp_dictionary_scis,
                                        working_dir = os.getcwd())
    # ### (3) ###  Extract sources that fall "close" to 'INDEF' regions.
    # Take out any sources from the white-light source list falling within 'remove_radius' of a flag.

    # ### (4) ### Feed corrected whitelight source lists into daophot with science images

    # ### (5) ### Gather columns and put in nice format (dictated by: "column_keys_phot.cfg")
    # This will convert columns from xy to ra and dec (controlled by: "column_keys_phot.cfg")


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


def create_rms_image():
    """Creates RMS image for use by source-finding code.

    Parameters
    ----------


    Returns
    -------
    rms_imgname : string
        filename of freshly computed RMS image
    """
    rms_img_filename = "foo_rms.fits"

    return(rms_img_filename)


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
    os.system("clear")
    for key1 in list(obs_info_dict.keys()):
        for key2 in list(obs_info_dict[key1].keys()):
            log.info("obs_info_dict[{}][{}]: {}".format(key1,key2,obs_info_dict[key1][key2]))  # TODO: REMOVE THIS SECTION BEFORE ACTUAL USE

    log.info("----------------------------------------------------------------------------------------------------------------------")
    log.info("SOURCELIST CREATION OCCURS HERE!")

    for tdp_keyname in [oid_key for oid_key in list(obs_info_dict.keys()) if
                        oid_key.startswith('total detection product')]:  # loop over total filtered products
        # 0: Map image filename to correspoinding catalog filename for total detection product and the associated
        # filter products
        totdet_product_cat_dict = {}
        filter_product_cat_dict = {}
        totdet_product_cat_dict[obs_info_dict[tdp_keyname]['product filenames']['image']] = \
            obs_info_dict[tdp_keyname]['product filenames']['source catalog']
        for fp_keyname in obs_info_dict[tdp_keyname]['associated filter products']:
            filter_product_cat_dict[obs_info_dict[fp_keyname]['product filenames']['image']] = \
                obs_info_dict[fp_keyname]['product filenames']['source catalog']

        inst_det = "{} {}".format(obs_info_dict[tdp_keyname]['info'].split()[-2],
                                  obs_info_dict[tdp_keyname]['info'].split()[-1])

        # 1: Generate source extractor-like sourcelist(s)
        create_se_like_sourcelists()

        # 2: Generate daophot-like sourcelist(s)
        create_daophot_like_sourcelists(totdet_product_cat_dict,filter_product_cat_dict,inst_det,param_dict[inst_det])


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


def run_daofind(param_dict, filelist=None, source_match=50000., verbose=True,whitelightimage=None, whitelightrms=None,
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

    whitelightrms : string
        Name of the multi-filter composite RMS image produced by hla_reduction.py.Default value is 'None'.

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
    # medDivImg,wht_data = Create_MedDivImage(whitelightimage)
    # rms_array = fits.getdata(whitelightrms,0)
    # rms_image_median = Util.binmode(rms_array[numpy.isfinite(rms_array) & (rms_array > 0.0)])[0]
    # #rms_image_median = numpy.median(rms_array[numpy.isfinite(rms_array) & (rms_array > 0.0)])
    # log.info("Median from RMS image = {}".format(rms_image_median))
    #
    # daoParams["sigma"] = rms_image_median

    whitelightdata = fits.getdata(whitelightimage, 1)
    mean, median, std = sigma_clipped_stats(whitelightdata, sigma=3.0) # TODO: Quick and dirty estimation of background RMS to move things along. ** REFINE PRIOR TO DEPLOYMENT **
    daoParams["sigma"] = std

    log.info('white light rms image = {}'.format(whitelightrms))
    log.info('sigma = {}'.format(std))
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

