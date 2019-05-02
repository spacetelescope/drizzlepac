#!/usr/bin/env python

"""This script contains code to support creation of source extractor-like and daophot-like sourcelists.

"""
import os
import pdb
import sys
import traceback

from astropy.io import fits
import numpy
from stsci.tools import logutil

__taskname__ = 'sourcelist_generation'

log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout)

# ----------------------------------------------------------------------------------------------------------------------
# set up instrument/detector-specific params
# Params imported from the following HLA classic parameter files:
# - https://grit.stsci.edu/HLA/hla/tree/master/software/trunk/HLApipeline/HLApipe/param/parameter_acs_hrc.cfg
# - https://grit.stsci.edu/HLA/hla/tree/master/software/trunk/HLApipeline/HLApipe/param/parameter_acs_sbc.cfg
# - https://grit.stsci.edu/HLA/hla/tree/master/software/trunk/HLApipeline/HLApipe/param/parameter_acs_wfc.cfg
# - https://grit.stsci.edu/HLA/hla/tree/master/software/trunk/HLApipeline/HLApipe/param/parameter_wfc3_ir.cfg
# - https://grit.stsci.edu/HLA/hla/tree/master/software/trunk/HLApipeline/HLApipe/param/parameter_wfc3_uvis.cfg

phot_param_dict = {
    "ACS HRC": {
        "dao": {
            "TWEAK_FWHMPSF": 0.073,
            "TWEAK_THRESHOLD": 3.0,
            "aperture_1": 0.03,
            "aperture_2": 0.125,
            "bthresh": 5.0}},
    "ACS SBC": {
        "dao": {
            "TWEAK_FWHMPSF": 0.065,
            "TWEAK_THRESHOLD": 3.0,
            "aperture_1": 0.07,
            "aperture_2": 0.125,
            "bthresh": 5.0}},
    "ACS WFC": {
        "dao":{
            "TWEAK_FWHMPSF": 0.076,
            "TWEAK_THRESHOLD": 3.0,
            "aperture_1": 0.05,  # update from 0.15
            "aperture_2": 0.15,  # update from 0.25
            "bthresh": 5.0}},
    "WFC3 IR": {
        "dao": {
            "TWEAK_FWHMPSF": 0.14,
            "TWEAK_THRESHOLD": 3.0,
            "aperture_1": 0.15,
            "aperture_2": 0.45,
            "bthresh": 5.0}},
    "WFC3 UVIS": {
        "dao": {
            "TWEAK_FWHMPSF": 0.076,
            "TWEAK_THRESHOLD": 3.0,
            "aperture_1": 0.05,
            "aperture_2": 0.15,
            "bthresh": 5.0}}}

# ----------------------------------------------------------------------------------------------------------------------


def create_daophot_like_sourcelists(totdet_product_cat_dict,filter_product_cat_dict,inst_det):
    """Make daophot-like sourcelists

    Parameters
    ----------
    totdet_product_cat_dict : dictionary
        Dictionary mapping image filename to corresponding catalog filename for total detection product

    filter_product_cat_dict : dictionary
        Dictionary mapping image filename to corresponding catalog filename for total detection product

    inst_det : string
        Text string containing space-deliminated instrument name, detector name (upper case) (i.e. WFC3_UVIS)

    Returns
    -------
    """

    log.info("DAOPHOT-LIKE SOURCELIST CREATION OCCURS HERE!")

    totfiltprod_filename_list = list(filter_product_cat_dict.keys())
    tdp_imagename = list(totdet_product_cat_dict.keys())[0]
    tdp_catname = list(totdet_product_cat_dict.values())[0]
    # ### (1) ### Collect applicable parameters
    log.info("### (1) ### Collect applicable parameters")

    fwhm = float(phot_param_dict[inst_det]["dao"]["TWEAK_FWHMPSF"])
    thresh = float(phot_param_dict[inst_det]["dao"]["TWEAK_THRESHOLD"])
    ap_diameter1 = float(phot_param_dict[inst_det]["dao"]["aperture_1"])
    ap_diameter2 = float(phot_param_dict[inst_det]["dao"]["aperture_2"])
    daofind_basic_param = [fwhm, thresh, ap_diameter1, ap_diameter2]

    # ----------------------------------------
    # Calculate mean readnoise value (float):
    # ----------------------------------------
    readnoise_dictionary_drzs = get_readnoise(totfiltprod_filename_list)

    # -----------------------------
    # Get scale arcseconds / pixel
    # -----------------------------
    scale_dict_drzs = stwcs_get_scale(totfiltprod_filename_list)

    for img in readnoise_dictionary_drzs.keys():
        log.info("{} {} {}".format(img,readnoise_dictionary_drzs[img],scale_dict_drzs[img]))


    # ### (2) ###  Use daostarfinder to create a sourcelist from the total detection image.

    log.info('### (2) ###  Use daostarfinder to create a sourcelist from the total detection image {}'.format(tdp_imagename))
    daofind_white_sources = run_daofind(config_file,
                                        sourcelist_create = True,
                                        whitelightimage = whitelightimage_string,
                                        whitelightrms = whitelightrms_string,
                                        daofind_basic_param = daofind_basic_param,
                                        readnoise_dictionary_drzs = readnoise_dictionary_drzs,
                                        scale_dict_drzs = scale_dict_drzs,
                                        exp_dictionary_scis = exp_dictionary_scis,
                                        working_dir = working_hla_red,
                                        detector = detector)
    # ### (3) ###  Extract sources that fall "close" to 'INDEF' regions.
    # Take out any sources from the white-light source list falling within 'remove_radius' of a flag.

    # ### (4) ### Feed corrected whitelight source lists into daophot with science images

    # ### (5) ### Gather columns and put in nice format (dictated by: "column_keys_phot.cfg")
    # This will convert columns from xy to ra and dec (controlled by: "column_keys_phot.cfg")


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


def create_sourcelists(obs_info_dict):
    """Make sourcelists

    Parameters
    ----------
    obs_info_dict : dictionary
        Dictionary containing all information about the images being processed

    Returns
    -------
    """
    print("----------------------------------------------------------------------------------------------------------------------")
    os.system("clear")
    for key1 in obs_info_dict.keys():
        for key2 in obs_info_dict[key1].keys():
            print(key1,key2,obs_info_dict[key1][key2])   # TODO: REMOVE THIS SECTION BEFORE ACTUAL USE
        print()
    print("----------------------------------------------------------------------------------------------------------------------")
    log.info("SOURCELIST CREATION OCCURS HERE!")

    for tdp_keyname in [oid_key for oid_key in list(obs_info_dict.keys()) if
                        oid_key.startswith('total detection product')]:  # loop over total filtered products
        # 0: Map image filename to correspoinding catalog filename for total detection product and the associated filter products
        totdet_product_cat_dict = {}
        filter_product_cat_dict = {}
        totdet_product_cat_dict[obs_info_dict[tdp_keyname]['product filenames']['image']] = obs_info_dict[tdp_keyname]['product filenames']['source catalog']
        for fp_keyname in obs_info_dict[tdp_keyname]['associated filter products']:
            filter_product_cat_dict[obs_info_dict[fp_keyname]['product filenames']['image']] = obs_info_dict[fp_keyname]['product filenames']['source catalog']

        inst_det = "{} {}".format(obs_info_dict[tdp_keyname]['info'].split()[-2],
                                  obs_info_dict[tdp_keyname]['info'].split()[-1])
        # 1: Generate daophot-like sourcelist(s)
        create_daophot_like_sourcelists(totdet_product_cat_dict,filter_product_cat_dict,inst_det)

        # 2: Generate source extractor-like sourcelist(s)
        create_se_like_sourcelists()


# ----------------------------------------------------------------------------------------------------------------------


def extract_name(stringWpath):
    """
    This task will extract just the name of a specific filename that includes the path in the name: 'stringWpath'.

    Tested.

    :param stringWpath: filename with full path
    :type stringWpath: string
    :returns: naked filename stripped of its path
    """
    while "/" == stringWpath[-1]:
        stringWpath = stringWpath[:-1]
    stringname = stringWpath.split("/")[-1]

    return stringname


# ----------------------------------------------------------------------------------------------------------------------


def get_mean_readnoise(image):
    """
    This subroutine computes mean readnoise values

    :param image: image filename
    :type image: string
    :return: mean readnoise
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
    This task will grab the average readnoise for HST
    data from the header.  This task is known to call the
    correct header keys for ACS UVIS data as well as WFC3
    UVIS and IR data, and WFPC2 data.

    Tested.

    :param listofimages: list of images that will be used to get readnoise values
    :type listofimages: list
    :returns: A dictionary of read noise values keyed by image name
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


def stwcs_get_scale(listofimages):
    """
    This task will grab the arcsec/pixel scale for HST
    data from the WCS information of the header.

    Note: Assumes science image is extension: "[1]".

    Tested.

    :param listofimages: list of images that will be used to get scale values
    :type listofimages: list
    :returns: A dictionary of scale values keyed by image name
    """
    import stwcs
    from stwcs import wcsutil
    dictionary_output = {}
    for individual_image in listofimages:
        try:
            wcs1 = stwcs.wcsutil.HSTWCS(individual_image + "[1]")
            scale2return = wcs1.pscale
        except:  # XXX what kind of exception here?
            print
            "ALERT: Pixel-scale could not be gathered from the header for image %s." % (extract_name(individual_image))
            scale2return = False
        dictionary_output[individual_image] = scale2return

    return dictionary_output


# ----------------------------------------------------------------------------------------------------------------------


def run_daofind(config_file, filelist=None, source_match=50000., verbose=True,
                sourcelist_create=False, whitelightimage=None, whitelightrms=None,
                daofind_basic_param=None, readnoise_dictionary_drzs=None, scale_dict_drzs=None,
                exp_dictionary_scis=None, working_dir=None, detector=None, sl_ext = 0,
                sharphi=None, sharplo=None, edgemask=5):
    """
    Generates sourcelists using DAOfind.

    :param config_file:  name of the detector-specific configuration file to use.
    :param filelist: Default value is 'None'.
    :param source_match: Source matching list length threshold. Default value is '50000.'.
    :param verbose: Print verbose output? Default value is 'True'.
    :param sourcelist_create: Create DAOphot sourcelists? Default value is 'False'.
    :param whitelightimage: Name of the multi-filter composite image produced by hla_reduction.py. Default value is 'None'.
    :param whitelightrms: Name of the multi-filter composite RMS image produced by hla_reduction.py.Default value is 'None'.
    :param daofind_basic_param: **UNUSED** List of values that will be used for daofind parameters 'fwhm', 'thresh', 'ap_diameter1', 'ap_diameter2'. Default value is 'None'.
    :param readnoise_dictionary_drzs: Dictionary containing readnoise values keyed by filter-specific drizzled image filename. Default value is 'None'.
    :param scale_dict_drzs: **UNUSED** Dictionary containing scale values keyed by filter-specific drizzled image filename. Default value is 'None'.
    :param exp_dictionary_scis: Dictionary containing exposure time values keyed by filter-specific drizzled image filename. Default value is 'None'.
    :param working_dir: Working directory. Default value is 'None'.
    :param detector: **UNUSED** Detector name. Default value is 'None'.
    :param sl_ext: FITS image extension to perform DAOfind on. Default value is '0'.
    :param sharphi: Upper limit on SHARPNESS parameter. Default value is 'None'.
    :param sharplo: Lower limit on SHARPNESS parameter. Default value is 'None'.
    :param edgemask: Remove sources within edgemask pixels of image edge. Default value is '5.'.
    :type config_file: string
    :type filelist: list
    :type source_match: float
    :type verbose: Boolean
    :type sourcelist_create: Boolean
    :type whitelightimage: string
    :type whitelightrms: string
    :type daofind_basic_param: List
    :type readnoise_dictionary_drzs: Dictionary
    :type scale_dict_drzs: Dictionary
    :type exp_dictionary_scis: Dictionary
    :type working_dir: string
    :type detector: string
    :type sl_ext: integer
    :type sharphi: float
    :type sharplo: float
    :type edgemask: integer
    :returns: if *sourcelist_create* is True, the name of daofind output coordinate file. If False, a dictionary of coo filenames, keyed by image name.
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
    config.read(software+'/param/'+config_file)
    fwhm = float(Configs.loadcfgs(config,'DAOFIND PARAMETERS','TWEAK_FWHMPSF'))
    thresh = float(Configs.loadcfgs(config,'DAOFIND PARAMETERS','TWEAK_THRESHOLD'))
    scale = float(Configs.loadcfgs(config,'ASTRODRIZZLE PARAMETERS','PIXSCALE'))


    daoParams["fwhm"] = fwhm
    daoParams["threshold"] = thresh
    daoParams["scale"] = scale
    daoParams["ratio"] = 0.8
    print ' '
    print 'run_daofind INPUT PARAMETERS:'
    print '-----------------------------'
    print 'config_file = ',software+'/param/'+config_file
    print 'fwhm = ',fwhm
    print 'thresh = ',thresh
    print 'scale = ',scale
    if sharphi:
        print 'sharphi = ',sharphi
        daoParams["sharphi"] = sharphi

    if sharplo:
        print 'sharplo = ',sharplo
        daoParams["sharplo"] = sharplo
    if sourcelist_create:

        readnoise = average_values_from_dict(readnoise_dictionary_drzs)
        exptime = average_values_from_dict(exp_dictionary_scis)

        # ----------------------------
        # Create Median-Divided Image
        # ----------------------------
        medDivImg,wht_data = Create_MedDivImage(whitelightimage)
        rms_array = pyfits.getdata(whitelightrms,0)
        rms_image_median = Util.binmode(rms_array[numpy.isfinite(rms_array) & (rms_array > 0.0)])[0]
        #rms_image_median = numpy.median(rms_array[numpy.isfinite(rms_array) & (rms_array > 0.0)])
        print "Median from RMS image = ",rms_image_median

        daoParams["sigma"] = rms_image_median

        print 'white light rms image = ',whitelightrms
        print 'sigma = ',rms_image_median
        print 'readnoise = ',readnoise
        print 'exptime = ',exptime
        print ' '

        name_daoOUT = extract_name(whitelightimage)
        name_daoOUT = Rename.unique_name(name_daoOUT + ".coo", suffix=".coo")
        output_dao = os.path.join(working_dir, name_daoOUT)

        try:
            print "image = ",medDivImg
            print "output = ",output_dao

            run_DAOStarFinder(medDivImg,output_dao,daoParams,debug=False)
        except: #XXX what kind of exception here?
            print ' '
            print '****************************************************************'
            print 'WARNING: THE MEDIAN-DIVIDED IMAGE CONTAINS MULTIPLE EXTENSIONS; '
            print '                  STARTING FLAG AND FILTER PROCESSING.          '
            print '****************************************************************'
            print ' '
            print "image = ",medDivImg+"[1]"
            print "output = ",output_dao
            run_DAOStarFinder(medDivImg + "[1]", output_dao, daoParams, debug=False)
        reject_image = numpy.zeros(wht_data.shape,dtype=numpy.int16)
        reject_image[numpy.where(wht_data <= 0.)] = 1

        rej_img = whitelightimage+'_rej_img.fits'
        pyfits.append(rej_img, numpy.float32(reject_image))

        mod_output_dao = output_dao+".mod"
        filter_daolist(output_dao, mod_output_dao, reject_image, edgemask)

        os.rename(output_dao, output_dao+".OLD ")
        os.rename(mod_output_dao, output_dao)

        return output_dao

    # ---------------------------------------------------
    # Run DAOFIND For The Creation of AstroDrizzle Input
    # ---------------------------------------------------
    coo_dict={}
    orig_thresh = thresh
    for image in filelist:
        try:
            readnoise = Headers.get_mean_readnoise(image)
        except (KeyError, IndexError):
            print "ALERT: Readnoise could not be gathered from the header for image %s."%(image)
            readnoise = 0.0

        exptime = pyfits.getheader(image)['EXPTIME']

        # estimate rms directly from image
        rms_array = pyfits.getdata(image,0)
        rms_image_median = rms_from_image(rms_array)
        print "Median from RMS MAD estimate = ",rms_image_median

        daoParams["sigma"] = rms_image_median
        daoParams["threshold"] = thresh
        daoParams["ratio"] = 0.8

        print 'sigma = ',rms_image_median
        print 'readnoise = ',readnoise
        print 'exptime = ',exptime
        print 'image = ',image

        if verbose: print "Finding sources in %s" %(image.split('/')[-1])
        outcoo=image.split('/')[-1]+'.coo'

        thresh = orig_thresh
        thresh1 = ns1 = None  # previous threshold and source count
        while True:
            daoParams["threshold"] = thresh
            print "image = ",image+'[%d]' %sl_ext
            print "output = ",outcoo
            print "verify = no"
            run_DAOStarFinder(image+'[%d]' %sl_ext, outcoo, daoParams, debug=False)

            infile = open(outcoo,'r')
            image_coo = infile.readlines()
            infile.close()

            if len(image_coo) <= source_match:
                if edgemask:
                    # create a mask that is true for bad pixels
                    if verbose:
                        print "Masking within",edgemask,"pixels of edge"
                    fh = pyfits.open(image)
                    mask = fh[sl_ext].data == 0
                    fh.close()
                    count = filter_daolist(outcoo, outcoo, mask, edgemask)
                    infile = open(outcoo,'r')
                    image_coo = infile.readlines()
                    infile.close()
                break

            # if number of sources > source_match, increase threshold
            thresh0 = thresh1
            ns0 = ns1
            thresh1 = thresh
            ns1 = len(image_coo)
            if thresh0 is None:
                thresh = thresh*1.1
            else:
                # linear projection from last 2 values with a little acceleration to get better threshold
                # shoot for 5% below the maximum threshold
                thresh = thresh0 + 1.3 * (thresh1-thresh0)/(ns1-ns0) * (0.95*source_match-ns0)
            print "# of stars %d > source match %d, upping threshold from %3.2f to %3.2f" %(ns1,
                                                                                            source_match,
                                                                                            thresh1,
                                                                                            thresh)
            os.remove(outcoo)

        if verbose:
            print '%d sources above %5.2f-sigma added to %s\n' %(len(image_coo),float(thresh), outcoo)
        coo_dict[outcoo] = len(image_coo)
    return(coo_dict)


# ----------------------------------------------------------------------------------------------------------------------
