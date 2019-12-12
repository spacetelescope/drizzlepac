#!/usr/bin/env python
"""Utility to create the Sextractor-like source catalog using PhotUtils

The function, create_sextractor_like_sourcelists, processes the provided
input file (presumptively a drizzled, white-light image) looking for sources.
A segmentation map and the total detection product catalog are the critical
outputs of this routine.

The function, measure_source_properties, uses the segmentation map generated
by create_sextractor_like_sourcelists, and tries to measure various properties
of the objects at the location of the found sources.  The measurements are
done in the drizzled, specific filter images of the visit. The filter product
is the output of the routine.

"""
import sys
# from distutils.version import LooseVersion

# import numpy as np
try:
    from matplotlib import pyplot as plt
except Exception:
    plt = None

import astropy.units as u
# from astropy.io import ascii
from astropy.io import fits as fits
from astropy.table import Column, Table
from astropy.convolution import Gaussian2DKernel, MexicanHat2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
# import photutils
# from photutils import detect_sources, source_properties, deblend_sources
# from photutils import Background2D, MedianBackground, SExtractorBackground, StdBackgroundRMS
from stwcs.wcsutil import HSTWCS
from stsci.tools import logutil


__taskname__ = "se_source_generation"

log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout)

_all__ = ["create_sextractor_like_sourcelists", "measure_source_properties"]

def create_sextractor_like_sourcelists(source_filename, catalog_filename, param_dict, se_debug=False):
    """Use photutils to find sources in image based on segmentation.

    Parameters
    ----------
    source_filename : string
        Filename of the "white light" drizzled image (aka the total detection product) which
        is used for the detection of sources

    catalog_filename : string
        Name of the output source catalog for the total detection product

    param_dict : dictionary
        dictionary of drizzle, source finding, and photometric parameters

    se_debug : bool, optional
        Specify whether or not to plot the image and segmentation image for
        visualization and debugging purposes

    Returns
    -------
    segm : `photutils.segmentation.SegmentationImage`
        Two-dimensional segmentation image where found source regions are labeled with
        unique, non-zero positive integers.

    kernel :

    bkg : `~photutils.background.Background2D` or None
        A background map based upon the `~photutils.background.SExtractorBackground`
        estimator

    bkg_rms_mean : float
        Mean bkg.background FIX

    """

    # Open the "white light" image and get the SCI image data
    imghdu = fits.open(source_filename)
    imgarr = imghdu['sci',1].data

    # Get the HSTWCS object from the first extension
    imgwcs = HSTWCS(imghdu, 1)

    # Get header information to annotate the output catalogs
    keyword_dict = _get_header_data(imghdu)

    # Get the instrument/detector-specific values from the param_dict
    fwhm = param_dict["sourcex"]["fwhm"]
    size_source_box = param_dict["sourcex"]["source_box"]
    threshold_flag = param_dict["sourcex"]["thresh"]

    # Report configuration values to log
    log.info("{}".format("=" * 80))
    log.info("")
    log.info("SExtractor-like source finding settings for Photutils segmentation")
    log.info("Total Detection Product - Input Parameters")
    log.info("FWHM: {}".format(fwhm))
    log.info("size_source_box: {}".format(size_source_box))
    log.info("threshold_flag: {}".format(threshold_flag))
    log.info("")
    log.info("{}".format("=" * 80))

    # Only use a single kernel for now
    kernel_list = [Gaussian2DKernel, MexicanHat2DKernel]
    kernel_in_use = kernel_list[0]

    bkg, bkg_dao_rms, threshold = _compute_background(imgarr, nsigma=5., threshold_flag=threshold_flag)

    # FIX imgarr should be background subtracted, sextractor uses the filtered_data image
    imgarr_bkgsub = imgarr - bkg.background

    # *** FIX: should size_source_box size be used in all these places? ***
    # Create a 2D filter kernel - this will be used to smooth the input
    # image prior to thresholding in detect_sources().
    sigma = fwhm * gaussian_fwhm_to_sigma
    kernel = kernel_in_use(sigma, x_size=size_source_box, y_size=size_source_box)
    kernel.normalize()

    # Source segmentation/extraction
    # If the threshold includes the background level, then the input image
    # should NOT be background subtracted.
    # Note: SExtractor has "connectivity=8" which is the default for this function
    segm = detect_sources(imgarr, threshold, npixels=size_source_box,
                          filter_kernel=kernel)

    # For debugging purposes...
    if se_debug:
        # Write out a catalog which can be used as an overlay for image in ds9
        cat = source_properties(imgarr_bkgsub, segm, background=bkg.background,
                                filter_kernel=kernel, wcs=imgwcs)
        table = cat.to_table()

        # Copy out only the X and Y coordinates to a "debug table" and
        # cast as an Astropy Table
        tbl = Table(table["xcentroid", "ycentroid"])

        # Construct the debug output filename and write the catalog
        indx = catalog_filename.find("ecsv")
        outname = catalog_filename[0:indx] + "reg"

        tbl["xcentroid"].info.format = ".10f"  # optional format
        tbl["ycentroid"].info.format = ".10f"

        # Add one to the X and Y table values to put the data onto a one-based system,
        # particularly for display with DS9
        tbl["xcentroid"] = tbl["xcentroid"] + 1
        tbl["ycentroid"] = tbl["ycentroid"] + 1
        tbl.write(outname, format="ascii.commented_header")
        log.info("Wrote debug source catalog: {}".format(outname))

        """
        # Generate a graphic of the image and the segmented image
        norm = ImageNormalize(stretch=SqrtStretch())
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))
        ax1.imshow(imgarr, origin="lower", cmap="Greys_r", norm=norm)
        ax1.set_title("Data")
        ax2.imshow(segm, origin="lower", cmap=segm.cmap(random_state=12345))
        ax2.set_title("Segmentation Image")
        plt.show()
        """

    # TROUBLESOME at this time
    # Deblending is a combination of multi-thresholding and watershed
    # segmentation. Sextractor uses a multi-thresholding technique.
    # npixels = number of connected pixels in source
    # npixels and filter_kernel should match those used by detect_sources()
    # Note: SExtractor has "connectivity=8" which is the default for this function
    """
    segm = deblend_sources(imgarr, segm, npixels=size_source_box,
                           filter_kernel=kernel, nlevels=32,
                           contrast=0.005)
    print("after deblend. ", segm)
    """

    """
    if se_debug:
        # Generate a graphic of the image and the segmented image
        norm = ImageNormalize(stretch=SqrtStretch())
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))
        ax1.imshow(imgarr, origin="lower", cmap="Greys_r", norm=norm)
        ax1.set_title("Data")
        ax2.imshow(segm, origin="lower", cmap=segm.cmap(random_state=12345))
        ax2.set_title("Segmentation Image")
        plt.show()
    """

    # Regenerate the source catalog with presumably now only good sources
    seg_cat = source_properties(imgarr_bkgsub, segm, background=bkg.background,
                                filter_kernel=kernel, wcs=imgwcs)

    _write_catalog(seg_cat, keyword_dict, catalog_filename)

    return segm, kernel, bkg_dao_rms


def measure_source_properties(segm, kernel, source_filename, catalog_filename, param_dict):
    """Use the positions of the sources identified in the white light image to
    measure properties of these sources in the filter images

    An instrument/detector combination may have multiple filter-level products.
    This routine is called for each filter image which is then measured to generate
    a filter-level source catalog based on object positions measured in the total
    detection product image.

    Parameters
    ----------
    segm : `~astropy.photutils.segmentation` Segmentation image
        Two-dimensional image of labeled source regions based on the "white light" drizzed product

    kernel : `~astropy.convolution`
        Two dimensional function of a specified FWHM used to smooth the image and
        used in the detection of sources as well as for the determination of the
        source properties (this routine)

    source_filename : string
        Filename of the filter drizzled image (aka the filter detection product) which
        is used for the measurement of properties of the previously found sources

    catalog_filename : string
        Name of the output source catalog for the filter detection product

    param_dict : dictionary
        dictionary of drizzle, source finding, and photometric parameters

    Returns
    -------

    """

    # Open the filter-level image
    imghdu = fits.open(source_filename)
    imgarr = imghdu['sci', 1].data

    # Get the HSTWCS object from the first extension
    imgwcs = HSTWCS(imghdu, 1)

    # Get header information to annotate the output catalogs
    keyword_dict = _get_header_data(imghdu, product="fdp")

    # Get the instrument/detector-specific values from the param_dict
    fwhm = param_dict["sourcex"]["fwhm"]
    size_source_box = param_dict["sourcex"]["source_box"]
    threshold_flag = param_dict["sourcex"]["thresh"]

    # Report configuration values to log
    log.info("{}".format("=" * 80))
    log.info("")
    log.info("SExtractor-like source property measurements based on Photutils segmentation")
    log.info("Filter Level Product - Input Parameters")
    log.info("FWHM: {}".format(fwhm))
    log.info("size_source_box: {}".format(size_source_box))
    log.info("threshold_flag: {}".format(threshold_flag))
    log.info("")
    log.info("{}".format("=" * 80))

    # The data needs to be background subtracted when computing the source properties
    bkg, _, _ = _compute_background(imgarr, nsigma=5., threshold_flag=threshold_flag)

    imgarr_bkgsub = imgarr - bkg.background

    # Compute source properties...
    seg_cat = source_properties(imgarr_bkgsub, segm, background=bkg.background,
                                filter_kernel=kernel, wcs=imgwcs)
    print(Table(seg_cat.to_table()).colnames)

    # Write the source catalog
    _write_catalog(seg_cat, keyword_dict, catalog_filename, product="fdp")


def _compute_background(image, box_size=50, win_size=3, nsigma=5., threshold_flag=None):
    """Use Background2D to determine the background of the input image.

    Parameters
    ----------
    image : ndarray
        Numpy array of the science extension from the observations FITS file.

    box_size : int
        Size of box along each axis

    win_size : int
        Size of 2D filter to apply to the background image

    nsigma : float
        Number of sigma above background

    threshold_flag : float or None
        Value from the image which serves as the limit for determining sources.
        If None, compute a default value of (background+5*rms(background)).
        If threshold < 0.0, use absolute value as scaling factor for default value.


    Returns
    -------
    bkg : 2D ndarray
        Background image

    bkg_dao_rms : ndarry
        Background RMS image

    threshold : ndarray
        Numpy array representing the background plus RMS

    """
    # Report configuration values to log
    log.info("")
    log.info("Computation of white light image background - Input Parameters")
    log.info("Box size: {}".format(box_size))
    log.info("Window size: {}".format(win_size))
    log.info("NSigma: {}".format(nsigma))

    # SExtractorBackground ans StdBackgroundRMS are the defaults
    bkg_estimator = SExtractorBackground()
    bkgrms_estimator = StdBackgroundRMS()
    bkg = None
    bkg_dao_rms = None

    exclude_percentiles = [10, 25, 50, 75]
    for percentile in exclude_percentiles:
        log.info("")
        log.info("Percentile in use: {}".format(percentile))
        try:
            bkg = Background2D(image, box_size,
                               filter_size=win_size,
                               bkg_estimator=bkg_estimator,
                               bkgrms_estimator=bkgrms_estimator,
                               exclude_percentile=percentile,
                               edge_method="pad")
        except Exception:
            bkg = None
            continue

        if bkg is not None:
            # Set the bkg_rms at "nsigma" sigma above background
            bkg_rms = nsigma * bkg.background_rms
            default_threshold = bkg.background + bkg_rms
            bkg_rms_mean = bkg.background.mean() + nsigma * bkg_rms.std()
            bkg_dao_rms = bkg.background_rms
            if threshold_flag is None:
                threshold = default_threshold
            elif threshold_flag < 0:
                threshold = -1 * threshold_flag * default_threshold
                log.info("Background threshold set to {} based on {}".format(threshold.max(), default_threshold.max()))
                bkg_rms_mean = threshold.max()
            else:
                bkg_rms_mean = 3. * threshold_flag
                threshold = bkg_rms_mean

            if bkg_rms_mean < 0:
                bkg_rms_mean = 0.
            break

    # If Background2D does not work at all, define default scalar values for
    # the background to be used in source identification
    if bkg is None:
        bkg_rms_mean = max(0.01, image.min())
        bkg_rms = nsigma * bkg_rms_mean
        bkg_dao_rms = bkg_rms_mean
        threshold = bkg_rms_mean + bkg_rms

    # *** FIX: Need to do something for bkg if bkg is None ***

    # Report other useful quantities
    log.info("")
    log.info("Mean background: {}".format(bkg.background.mean()))
    log.info("Mean threshold: {}".format(bkg_rms_mean))
    log.info("")
    log.info("{}".format("=" * 80))

    return bkg, bkg_dao_rms, threshold


def _write_catalog(seg_cat, keyword_dict, catalog_filename, product="tdp"):
    """Actually write the specified source catalog out to disk

    Parameters
    ----------
    seg_cat : list of `~photutils.SourceProperties` objects
        List of SourceProperties objects, one for each source found in the
        specified detection product

    keyword_dict : dict
        Dictionary containing FITS keyword values pertaining to the data product

    catalog_filename : str
        Official generated name for the output catalog

    product : str, optional
        Identification string for the catalog product being written.  This
        controls the data being put into the catalog product
    """

    # Convert the list of SourceProperties objects to a QTable and
    # document in column metadata Photutils columns which map to SExtractor columns
    seg_table = Table(seg_cat.to_table())
    radec_data = seg_table["sky_centroid_icrs"]
    ra_icrs = radec_data.ra.degree
    dec_icrs = radec_data.dec.degree

    num_sources = len(seg_table)

    # If the output is for the total detection product, then only
    # a subset of the full catalog is needed.
    if product.lower() == "tdp":

        # [x|y]centroid are in pixels, physical data coordinates
        seg_subset_table = seg_table["xcentroid", "ycentroid"]

        # Add metadata to the output subset table
        seg_subset_table = _annotate_table(seg_subset_table, keyword_dict, num_sources, product=product)

        seg_subset_table["xcentroid"].description = "SExtractor Column x_image"
        seg_subset_table["ycentroid"].description = "SExtractor Column y_image"
        seg_subset_table["RA_icrs"] = ra_icrs
        seg_subset_table["Dec_icrs"] = dec_icrs
        seg_subset_table["RA_icrs"].description = "SExtractor Column RA"
        seg_subset_table["Dec_icrs"].description = "SExtractor Column Dec"
        seg_subset_table["RA_icrs"].unit = u.deg
        seg_subset_table["Dec_icrs"].unit = u.deg

        # Write out the official total detection product source catalog
        seg_subset_table["xcentroid"].info.format = ".10f"
        seg_subset_table["ycentroid"].info.format = ".10f"
        seg_subset_table["RA_icrs"].info.format = ".10f"
        seg_subset_table["Dec_icrs"].info.format = ".10f"
        print("seg_subset_table (white light image): ", seg_subset_table)

        seg_subset_table.write(catalog_filename, format="ascii.ecsv")
        log.info("Wrote source catalog: {}".format(catalog_filename))

    # else the product is the "filter detection product"
    else:

        seg_table = _annotate_table(seg_table, keyword_dict, num_sources, product=product)

        # Rework the current table for output
        del seg_table["id"]
        del seg_table["sky_centroid"]
        del seg_table["sky_centroid_icrs"]
        rr = Column(ra_icrs, name="RA_icrs", description="SExtractor Column RA", unit=u.deg)
        dd = Column(dec_icrs, name="Dec_icrs", description="SExtractor Column Dec", unit=u.deg)
        seg_table.add_column(dd, index=2)
        seg_table.add_column(rr, index=2)

        # Add a description for columns which map to SExtractor catalog columns
        seg_table["xcentroid"].description = "SExtractor Column x_image"
        seg_table["ycentroid"].description = "SExtractor Column y_image"
        seg_table["background_at_centroid"].description = "SExtractor Column background"
        seg_table["source_sum"].description = "SExtractor Column flux_iso"
        seg_table["source_sum_err"].description = "SExtractor Column fluxerr_iso"
        # FIX: is mapping to _image or _world?  _image
        seg_table["cxx"].description = "SExtractor Column cxx_image, ellipse parameter"
        seg_table["cyy"].description = "SExtractor Column cyy_image, ellipse parameter"
        seg_table["cxy"].description = "SExtractor Column cxy_image, ellipse parameter"
        # FIX: is the mapping to _image or _world?
        seg_table["covar_sigx2"].description = "SExtractor Column x2_image, (0,0) element of covariance matrix"
        seg_table["covar_sigy2"].description = "SExtractor Column y2_image, (1,1) element of covariance matrix"
        seg_table["covar_sigxy"].description = "SExtractor Column xy_image, (0,1) and (1,0) elements of covariance matrix"

        if 'xmin' in seg_table:
            cols = ['xmin', 'xmax', 'ymin', 'ymax']
        else:
            cols = ['bbox_xmin', 'bbox_xmax', 'bbox_ymin', 'bbox_ymax']
        seg_table[cols[0]].description = "SExtractor Column xmin_image"
        seg_table[cols[1]].description = "SExtractor Column xmax_image"
        seg_table[cols[2]].description = "SExtractor Column ymin_image"
        seg_table[cols[3]].description = "SExtractor Column ymax_image"

        # Write out the official filter detection product source catalog
        seg_table["xcentroid"].info.format = ".10f"
        seg_table["ycentroid"].info.format = ".10f"
        seg_table["RA_icrs"].info.format = ".10f"
        seg_table["Dec_icrs"].info.format = ".10f"
        print("seg_table (filter): {}".format(seg_table))

        seg_table.write(catalog_filename, format="ascii.ecsv")
        log.info("Wrote filter source catalog: {}".format(catalog_filename))

def _get_header_data(imghdu, product="tdp"):
    """Read FITS keywords from the primary or extension header and store the
    information in a dictionary

    Parameters
    ----------
    imghdu : HDUList
        An HDU object pertaining to the FITS file being processed
    """

    keyword_dict = {}

    keyword_dict["proposal_id"] = imghdu[0].header["PROPOSID"]
    keyword_dict["image_file_name"] = imghdu[0].header['FILENAME'].upper()
    keyword_dict["target_name"] = imghdu[0].header["TARGNAME"].upper()
    keyword_dict["date_obs"] = imghdu[0].header["DATE-OBS"]
    keyword_dict["instrument"] = imghdu[0].header["INSTRUME"].upper()
    keyword_dict["detector"] = imghdu[0].header["DETECTOR"].upper()
    keyword_dict["target_ra"] = imghdu[0].header["RA_TARG"]
    keyword_dict["target_dec"] = imghdu[0].header["DEC_TARG"]
    keyword_dict["expo_start"] = imghdu[0].header["EXPSTART"]
    keyword_dict["texpo_time"] = imghdu[0].header["TEXPTIME"]
    keyword_dict["ccd_gain"] = imghdu[0].header["CCDGAIN"]
    keyword_dict["aperture_pa"] = imghdu[0].header["PA_V3"]

    # The total detection product has the FILTER keyword in
    # the primary header - read it for any instrument.
    #
    # For the filter detection product:
    # WFC3 only has FILTER, but ACS has FILTER1 and FILTER2
    # in the primary header.
    if product.lower() == "tdp":
        keyword_dict["filter"] = imghdu[0].header["FILTER"]
    # The filter detection product...
    else:
        if keyword_dict["instrument"] == "ACS":
            keyword_dict["filter1"] = imghdu[0].header["FILTER1"]
            keyword_dict["filter2"] = imghdu[0].header["FILTER2"]
        else:
            keyword_dict["filter1"] = imghdu[0].header["FILTER"]
            keyword_dict["filter2"] = ""

    # Get the HSTWCS object from the first extension
    keyword_dict["wcs_name"] = imghdu[1].header["WCSNAME"]
    keyword_dict["wcs_type"] = imghdu[1].header["WCSTYPE"]
    print('WCSTYPE: {}'.format(keyword_dict["wcs_type"]))
    keyword_dict["orientation"] = imghdu[1].header["ORIENTAT"]
    keyword_dict["aperture_ra"] = imghdu[1].header["RA_APER"]
    keyword_dict["aperture_dec"] = imghdu[1].header["DEC_APER"]

    return(keyword_dict)

def _annotate_table(data_table, keyword_dict, num_sources, product="tdp"):
    """Add state metadata to the output source catalog

    Parameters
    ----------
    data_table : QTable
        Table of source properties

    keyword_dict : dict
        Dictionary containing FITS keyword values pertaining to the data product

    num_sources : int
        Number of sources (items) in table

    product : str, optional
        Identification string for the catalog product being written.  This
        controls the data being put into the catalog product

    Returns
    -------
    data_table : QTable
        Table of source properties updatd to contain state metadata

    """

    data_table.meta["WCSNAME"] = keyword_dict["wcs_name"]
    data_table.meta["WCSTYPE"] = keyword_dict["wcs_type"]
    data_table.meta["Proposal ID"] = keyword_dict["proposal_id"]
    data_table.meta["Image File Name"] = keyword_dict['image_file_name']
    data_table.meta["Target Name"] = keyword_dict["target_name"]
    data_table.meta["Date Observed"] = keyword_dict["date_obs"]
    # FIX
    if product.lower() == "tdp":
        data_table.meta["Time Observed"] = " "
        data_table.meta["Filter"] = keyword_dict["filter"]
    else:
        data_table.meta["Time Observed"] = "FIX ME"
        data_table.meta["Filter 1"] = keyword_dict["filter1"]
        data_table.meta["Filter 2"] = keyword_dict["filter2"]
    data_table.meta["Instrument"] = keyword_dict["instrument"]
    data_table.meta["Detector"] = keyword_dict["detector"]
    data_table.meta["Target RA"] = keyword_dict["target_ra"]
    data_table.meta["Target DEC"] = keyword_dict["target_dec"]
    data_table.meta["Orientation"] = keyword_dict["orientation"]
    data_table.meta["Aperture RA"] = keyword_dict["aperture_ra"]
    data_table.meta["Aperture DEC"] = keyword_dict["aperture_dec"]
    data_table.meta["Aperture PA"] = keyword_dict["aperture_pa"]
    data_table.meta["Exposure Start"] = keyword_dict["expo_start"]
    data_table.meta["Total Exposure Time"] = keyword_dict["texpo_time"]
    data_table.meta["CCD Gain"] = keyword_dict["ccd_gain"]
    data_table.meta["Number of sources"] = num_sources
    data_table.meta[""] = " "
    data_table.meta[""] = "Absolute coordinates are in a zero-based coordinate system."

    return(data_table)


def run_photutils():

    """
    param_dict = {"sourcex" : {"fwhm" : 0.25, "source_box" : 5, "thresh" : None}}
    white_light_filename = "hst_11150_70_wfc3_ir_total_ia1s70_drz.fits"
    tdp_catalog_filename = "hst_11150_70_wfc3_ir_total_ia1s70_segment-cat.ecsv"
    fp_filename_1 = "hst_11150_70_wfc3_ir_f110w_ia1s70_drz.fits"
    fp_catalog_filename_1 = "hst_11150_70_wfc3_ir_f110w_ia1s70_segment-cat.ecsv"
    fp_filename_2 = "hst_11150_70_wfc3_ir_f160w_ia1s70_drz.fits"
    fp_catalog_filename_2 = "hst_11150_70_wfc3_ir_f160w_ia1s70_segment-cat.ecsv"

    param_dict = {"sourcex" : {"fwhm" : 0.13, "source_box" : 5, "thresh" : None}}
    white_light_filename = "hst_10595_06_acs_wfc_total_j9es06_drc.fits"
    tdp_catalog_filename = "hst_10595_06_acs_wfc_total_j9es06_segment-cat.ecsv"

    fp_filename_1 = "hst_10595_06_acs_wfc_f435w_j9es06_drc.fits"
    fp_catalog_filename_1 = "hst_10595_06_acs_wfc_f435w_j9es06_segment-cat.ecsv"

    fp_filename_2 = "hst_10595_06_acs_wfc_f606w_j9es06_drc.fits"
    fp_catalog_filename_2 = "hst_10595_06_acs_wfc_f606w_j9es06_segment-cat.ecsv"

    fp_filename_3 = "hst_10595_06_acs_wfc_f814w_j9es06_drc.fits"
    fp_catalog_filename_3 = "hst_10595_06_acs_wfc_f814w_j9es06_segment-cat.ecsv"
    """

    param_dict = {"sourcex" : {"fwhm" : 0.13, "source_box" : 5, "thresh" : None}}
    white_light_filename = "hst_10265_01_acs_wfc_total_j92c01_drc.fits"
    tdp_catalog_filename = "hst_10265_01_acs_wfc_total_j92c01_segment-cat.ecsv"

    fp_filename_1 = "hst_10265_01_acs_wfc_f606w_j92c01_drc.fits"
    fp_catalog_filename_1 = "hst_10265_01_acs_wfc_f606w_j92c01_segment-cat.ecsv"


    segmap, kernel, bkg_dao_rms = create_sextractor_like_sourcelists(white_light_filename,
                                                                     tdp_catalog_filename,
                                                                     param_dict,
                                                                     se_debug=False)
    measure_source_properties(segmap, kernel, fp_filename_1, fp_catalog_filename_1, param_dict)
    print("Measured filter 1")

    """
    measure_source_properties(segmap, kernel, fp_filename_2, fp_catalog_filename_2, param_dict)
    print("Measured filter 2")

    measure_source_properties(segmap, kernel, fp_filename_3, fp_catalog_filename_3, param_dict)
    print("Measured filter 3")
    """

    return
