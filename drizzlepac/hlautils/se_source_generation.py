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
from distutils.version import LooseVersion

import numpy as np
try:
    from matplotlib import pyplot as plt
except Exception:
    plt = None

import astropy.units as u
from astropy.io import ascii
from astropy.io import fits as fits
from astropy.table import Column
from astropy.convolution import Gaussian2DKernel, MexicanHat2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
import photutils
from photutils import detect_sources, source_properties, deblend_sources
from photutils import Background2D, MedianBackground, SExtractorBackground, StdBackgroundRMS
from stwcs.wcsutil import HSTWCS
from stsci.tools import logutil


__taskname__ = "se_source_generation"

log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout)

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

    # Open the "white light" image and get the image data
    imghdu = fits.open(source_filename)
    imgarr = imghdu[1].data

    # Get the HSTWCS object from the first extension
    imgwcs = HSTWCS(source_filename, 1)

    # Get the instrument/detector-specific values from the param_dict 
    fwhm = param_dict["sourcex"]["fwhm"]
    size_source_box = param_dict["sourcex"]["size_source_box"]
    threshold = param_dict["sourcex"]["threshold"]

    # Report configuration values to log
    log.info("{}".format("=" * 80))
    log.info("")
    log.info("SExtractor-like source finding settings for Photutils segmentation")
    log.info("Total Detection Product - Input Parameters")
    log.info("FWHM: {}".format(fwhm))
    log.info("size_source_box: {}".format(size_source_box))
    log.info("threshold: {}".format(threshold))
    log.info("")
    log.info("{}".format("=" * 80))

    # Only use a single kernel for now
    kernel_list = [Gaussian2DKernel, MexicanHat2DKernel]
    kernel_in_use = kernel_list[0]

    bkg, bkg_dao_rms, threshold = compute_background(imgarr, nsigma=5., threshold=threshold)

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
    segm = detect_sources(imgarr, threshold, npixels=size_source_box,
                          filter_kernel=kernel)

    # For debugging purposes...
    if se_debug:
        # Write out a catalog which can be used as an overlay for image in ds9
        cat = source_properties(imgarr_bkgsub, segm, filter_kernel=kernel, wcs=imgwcs)
        table = cat.to_table()

        # Copy out only the X and Y coordinates to a "debug table" 
        tbl = table["xcentroid", "ycentroid"]

        # Construct the debug output filename and write the catalog
        indx = catalog_filename.find("ecsv")
        outname = catalog_filename[0:indx] + "reg"

        tbl["xcentroid"].info.format = ".10f"  # optional format
        tbl["ycentroid"].info.format = ".10f"
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
    seg_cat = source_properties(imgarr_bkgsub, segm, filter_kernel=kernel, wcs=imgwcs)

    write_catalog(seg_cat, catalog_filename)

    # FIX: All of these may not be needed so clean up
    return segm, kernel, bkg, bkg_dao_rms, imgarr_bkgsub


def measure_source_properties(segm, imgarr_bkgsub, kernel, source_filename, catalog_filename, param_dict):
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
        Filename of the filter drizzled image (aka the filter data product) which
        is used for the measurement of properties of the previously found sources

    catalog_filename : string
        Name of the output source catalog for the filter data product

    param_dict : dictionary
        dictionary of drizzle, source finding, and photometric parameters

    Returns
    -------

    """

    # Open the filter-level image
    imghdu = fits.open(source_filename)
    imgarr = imghdu[1].data

    # Get the HSTWCS object from the first extension
    imgwcs = HSTWCS(source_filename, 1)

    # Get the instrument/detector-specific values from the param_dict 
    fwhm = param_dict["sourcex"]["fwhm"]
    size_source_box = param_dict["sourcex"]["size_source_box"]
    threshold = param_dict["sourcex"]["threshold"]

    # Report configuration values to log
    log.info("{}".format("=" * 80))
    log.info("")
    log.info("SExtractor-like source property measurements based on Photutils segmentation")
    log.info("Filter Level Product - Input Parameters")
    log.info("FWHM: {}".format(fwhm))
    log.info("size_source_box: {}".format(size_source_box))
    log.info("threshold: {}".format(threshold))
    log.info("")
    log.info("{}".format("=" * 80))

    # The data needs to be background subtracted when computing the source properties
    bkg, _, _ = compute_background(imgarr, nsigma=5., threshold=threshold)
    imgarr_bkgsubX = imgarr - bkg.background

    # Compute source properties...
    seg_cat = source_properties(imgarr_bkgsubX, segm, filter_kernel=kernel, wcs=imgwcs)

    # Write the source catalog
    write_catalog(seg_cat, catalog_filename, product="fdp")


def compute_background(image, box_size=50, win_size=3, nsigma=5., threshold=None):
    """Use photutils to find sources in image based on segmentation.

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
        bkg_rms = nsigma * bkg_rms_mean
        bkg_dao_rms = bkg_rms_mean
        threshold = bkg_rms_mean + bkg_rms

    # Report other useful quantities
    log.info("")
    log.info("Mean background: {}".format(bkg.background.mean()))
    log.info("Mean threshold: {}".format(bkg_rms_mean))
    log.info("")
    log.info("{}".format("=" * 80))

    return bkg, bkg_dao_rms, threshold

def write_catalog(seg_cat, catalog_filename, product="tdp"):

    # Convert the list of SourceProperties objects to a QTable and
    # manipulate and rename columns to match SExtractor output
    seg_table = seg_cat.to_table()
    radec_data = seg_table["sky_centroid_icrs"]
    ra_icrs = radec_data.ra.degree
    dec_icrs = radec_data.dec.degree

    # If the output is for the total detection product, then only
    # a subset of the full catalog is needed.
    if product.lower() == "tdp":

        # [x|y]centroid are in pixels, physical data coordinates
        seg_subset_table = seg_table["xcentroid", "ycentroid"]
        seg_subset_table.rename_column("xcentroid", "x_image")
        seg_subset_table.rename_column("ycentroid", "y_image")
        seg_subset_table["RA_icrs"] = ra_icrs
        seg_subset_table["Dec_icrs"] = dec_icrs

        # Write out the official total detection product source catalog
        seg_subset_table["x_image"].info.format = ".10f"
        seg_subset_table["y_image"].info.format = ".10f"
        seg_subset_table["RA_icrs"].info.format = ".10f"
        seg_subset_table["Dec_icrs"].info.format = ".10f"
        seg_subset_table["RA_icrs"].unit = u.deg
        seg_subset_table["Dec_icrs"].unit = u.deg
        print("seg_subset_table (white): ", seg_subset_table)

        seg_subset_table.write(catalog_filename, format="ascii.ecsv")
        log.info("Wrote source catalog: {}".format(catalog_filename))

    else:

        # Rework the current table for output
        del seg_table["id"]
        del seg_table["sky_centroid"]
        del seg_table["sky_centroid_icrs"]
        rr = Column(ra_icrs, name="RA_icrs")
        dd = Column(dec_icrs, name="Dec_icrs")
        seg_table.add_column(dd, index=2)
        seg_table.add_column(rr, index=2)

        # Rename column names to be more consistent with Sextractor column names
        seg_table.rename_column("xcentroid", "x_image")
        seg_table.rename_column("ycentroid", "y_image")
        seg_table.rename_column("background_at_centroid", "background")
        seg_table.rename_column("source_sum", "flux")
        seg_table.rename_column("source_sum_err", "flux_err")
        seg_table.rename_column("cxx", "cxx_image")
        seg_table.rename_column("cxy", "cxy_image")
        seg_table.rename_column("cyy", "cyy_image")
        seg_table.rename_column("covar_sigx2", "x2_image")
        seg_table.rename_column("covar_sigy2", "y2_image")
        seg_table.rename_column("covar_sigxy", "xy_image")

        # Write out the official total detection product source catalog
        seg_table["x_image"].info.format = ".10f"
        seg_table["y_image"].info.format = ".10f"
        seg_table["RA_icrs"].info.format = ".10f"
        seg_table["Dec_icrs"].info.format = ".10f"
        seg_table["RA_icrs"].unit = u.deg
        seg_table["Dec_icrs"].unit =  u.deg
        print("seg_table (filter): {}".format(seg_table))

        seg_table.write(catalog_filename, format="ascii.ecsv")
        log.info("Wrote filter source catalog: {}".format(catalog_filename))

def run_photutils():

    # white_light_filename = "j92c01010_drc.fits"
    """
    param_dict = {"sourcex" : {"fwhm" : 0.25, "size_source_box" : 5, "threshold" : None}}
    white_light_filename = "hst_11150_70_wfc3_ir_total_ia1s70_drz.fits"
    tdp_catalog_filename = "hst_11150_70_wfc3_ir_total_ia1s70_segment-cat.ecsv"
    fp_filename_1 = "hst_11150_70_wfc3_ir_f110w_ia1s70_drz.fits"
    fp_catalog_filename_1 = "hst_11150_70_wfc3_ir_f110w_ia1s70_segment-cat.ecsv"
    fp_filename_2 = "hst_11150_70_wfc3_ir_f160w_ia1s70_drz.fits"
    fp_catalog_filename_2 = "hst_11150_70_wfc3_ir_f160w_ia1s70_segment-cat.ecsv"

    """
    param_dict = {"sourcex" : {"fwhm" : 0.13, "size_source_box" : 5, "threshold" : None}}
    white_light_filename = "hst_10595_06_acs_wfc_total_j9es06_drc.fits"
    tdp_catalog_filename = "hst_10595_06_acs_wfc_total_j9es06_segment-cat.ecsv"

    fp_filename_1 = "hst_10595_06_acs_wfc_f435w_j9es06_drc.fits"
    fp_catalog_filename_1 = "hst_10595_06_acs_wfc_f435w_j9es06_segment-cat.ecsv"

    fp_filename_2 = "hst_10595_06_acs_wfc_f606w_j9es06_drc.fits"
    fp_catalog_filename_2 = "hst_10595_06_acs_wfc_f606w_j9es06_segment-cat.ecsv"

    fp_filename_3 = "hst_10595_06_acs_wfc_f814w_j9es06_drc.fits"
    fp_catalog_filename_3 = "hst_10595_06_acs_wfc_f814w_j9es06_segment-cat.ecsv"

    segmap, kernel, bkg, bkg_dao_rms, bkgsub  = create_sextractor_like_sourcelists(white_light_filename,
                                                                                            tdp_catalog_filename,
                                                                                            param_dict,
                                                                                            se_debug=True)
    measure_source_properties(segmap, bkgsub, kernel, fp_filename_1, fp_catalog_filename_1, param_dict)
    print("Measured filter 1")

    measure_source_properties(segmap, bkgsub, kernel, fp_filename_2, fp_catalog_filename_2, param_dict)
    print("Measured filter 2")

    measure_source_properties(segmap, bkgsub, kernel, fp_filename_3, fp_catalog_filename_3, param_dict)
    print("Measured filter 3")

    return
