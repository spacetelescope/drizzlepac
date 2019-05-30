#!/usr/bin/env python
"""Utility to create the Sextractor-like source catalog using PhotUtils

The function, create_astrometric_catalog, allows the user to query an
astrometric catalog online to generate a catalog of astrometric sources that
should fall within the field-of-view of all the input images.

This module relies on the definition of an environment variable to specify
the URL of the astrometric catalog to use for generating this
reference catalog. ::

    ASTROMETRIC_CATALOG_URL  -- URL of web service that can be queried to
                                obtain listing of astrometric sources,
                                sky coordinates, and magnitudes.

"""
import sys
from distutils.version import LooseVersion

import numpy as np
try:
    from matplotlib import pyplot as plt
except Exception:
    plt = None

from astropy.io import fits as fits
from astropy.convolution import Gaussian2DKernel, MexicanHat2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
import photutils
from photutils import detect_sources, source_properties, deblend_sources
from photutils import Background2D, MedianBackground, SExtractorBackground, StdBackgroundRMS
from stwcs.wcsutil import HSTWCS
from stsci.tools import logutil


__taskname__ = 'se_source_generation'


log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout)


# _all__ = ['run_photutils', 'create_sextractor_like_sourcelists', classify_sources']

# def create_sextractor_like_sourcelists(source_filename, catalog_filename, param_dict, se_debug=False):
def create_sextractor_like_sourcelists(source_filename, catalog_filename, se_debug=False):
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

    bkg_rms : float
        N times the bkg.background_rms where N = 5 FIX

    bkg_rms_mean : float
        Mean bkg.background FIX

    """

    # Open the "white light" image and get the image data
    imghdu = fits.open(source_filename)
    imgarr = imghdu[1].data

    # Get the HSTWCS object from the first extension
    imgwcs = HSTWCS(source_filename, 1)

    # Get the instrument/detector-specific values from the
    # param_dict for now - will eventually come from a
    # configuration file.
    # fwhm = float(param_dict['sourcex']['fwhm'])
    # size_source_box = float(param_dict['sourcex']['source_box'])
    # threshold = float(param_dict['sourcex']['threshold'])
    fwhm = 0.073
    size_source_box = 7
    threshold = 1.4

    # Report configuration values to log
    log.info('====================')
    log.info('')
    log.info('SExtractor-like source finding settings for Photutils segmentation')
    log.info('Total Detection Product')
    log.info('FWHM: {}'.format(fwhm))
    log.info('size_source_box: {}'.format(size_source_box))
    log.info('threshold: {}'.format(threshold))
    log.info('')
    log.info('====================')

    # Only use a single kernel for now
    kernel_list = [Gaussian2DKernel, MexicanHat2DKernel]
    kernel_in_use = kernel_list[0]

    bkg, bkg_rms, bkg_dao_rms, threshold = compute_background(imgarr, threshold=threshold)

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
    # FIX: This currently generates a bad: detect.py:132: RuntimeWarning: invalid value
    # encountered in greater check_normalization=True) > threshold)
    print('Threshold: {}'.format(threshold))
    segm = detect_sources(imgarr, threshold, npixels=size_source_box,
                          filter_kernel=kernel)

    # For debugging purposes...
    if se_debug:
        # Write out a catalog which can be used as an overlay for image in ds9
        cat = source_properties(imgarr_bkgsub, segm, filter_kernel=kernel, wcs=imgwcs)
        table = cat.to_table()
        cnames = table.colnames

        # Move `id` column from the first to the last column for ds9
        cnames.append(cnames[0])
        del cnames[0]
        tbl = table[cnames[0:2]]

        # Construct the debug output filename and write the catalog
        indx = catalog_filename.find('ecsv')
        outname = catalog_filename[0:indx] + 'reg'

        tbl['xcentroid'].info.format = '.10f'  # optional format
        tbl['ycentroid'].info.format = '.10f'
        tbl.write(outname, format='ascii.commented_header')
        log.info("Wrote debug source catalog: {}".format(outname))

        """
        # Generate a graphic of the image and the segmented image
        norm = ImageNormalize(stretch=SqrtStretch())
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))
        ax1.imshow(imgarr, origin='lower', cmap='Greys_r', norm=norm)
        ax1.set_title('Data')
        ax2.imshow(segm, origin='lower', cmap=segm.cmap(random_state=12345))
        ax2.set_title('Segmentation Image')
        plt.show()
        """

    # TROUBLESOME at this time
    # Deblending is a combination of multi-thresholding and watershed
    # segmentation. Sextractor uses a multi-thresholding technique.
    # npixels = number of connected pixels in source
    # npixels and filter_kernel should match those used by detect_sources()
    # segm = deblend_sources(imgarr, segm, npixels=size_source_box,
    #                       filter_kernel=kernel, nlevels=16,
    #                       contrast=0.01)
    # print('after deblend. ', segm)

    # Modify the segmentation map to clean out possible cosmic rays
    # Should imgarr be imgarr_bkg here?
    segm = modify_segmentation_map(imgarr_bkgsub, segm, kernel)

    # Regenerate the source catalog with presumably now only good sources
    seg_cat = source_properties(imgarr_bkgsub, segm, filter_kernel=kernel, wcs=imgwcs)
    seg_table = seg_cat.to_table()
    radec_data = seg_table['sky_centroid_icrs']
    ra_icrs = radec_data.ra.degree
    dec_icrs = radec_data.dec.degree
    print('seg_table (white): ', seg_table)

    # Construct a table with every value in its own column, rename columns to
    # map to SExtractor output, and possibly only write out a subset of all
    # the possible columns
    #
    # [x|y]centroid are in pixels, physical data coordinates
    seg_subset_table = seg_table['xcentroid', 'ycentroid']
    # RA and Dec are decimal degrees ICRS
    seg_subset_table['RA_icrs'] = ra_icrs
    seg_subset_table['Dec_icrs'] = dec_icrs

    # Write out the official total detection product source catalog
    seg_subset_table['xcentroid'].info.format = '.10f'
    seg_subset_table['ycentroid'].info.format = '.10f'
    seg_subset_table['RA_icrs'].info.format = '.10f'
    seg_subset_table['Dec_icrs'].info.format = '.10f'

    seg_subset_table.write(catalog_filename, format='ascii.commented_header')
    log.info("Wrote source catalog: {}".format(catalog_filename))

    # FIX: All of these may not be needed so clean up
    return segm, kernel, bkg, bkg_rms, bkg_dao_rms


# def measure_source_properties(segm, kernel, source_filename, catalog_filename, param_dict):
def measure_source_properties(segm, kernel, source_filename, catalog_filename):
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

    # Get the instrument/detector-specific values from the
    # param_dict for now - will eventually come from a
    # configuration file.
    # fwhm = float(param_dict['sourcex']['fwhm'])
    # size_source_box = float(param_dict['sourcex']['source_box'])
    # threshold = float(param_dict['sourcex']['threshold'])
    fwhm = 0.073
    size_source_box = 7
    threshold = 1.4

    # Report configuration values to log
    log.info('====================')
    log.info('')
    log.info('SExtractor-like source property measurements based on Photutils segmentation')
    log.info('Filter Level Product')
    log.info('FWHM: {}'.format(fwhm))
    log.info('size_source_box: {}'.format(size_source_box))
    log.info('threshold: {}'.format(threshold))
    log.info('')
    log.info('====================')

    # The data needs to be background subtracted when computing the source properties
    bkg, _, _, _ = compute_background(imgarr, threshold=threshold)
    imgarr_bkgsub = imgarr - bkg.background

    # Compute source properties...
    seg_cat = source_properties(imgarr_bkgsub, segm, filter_kernel=kernel, wcs=imgwcs)

    # ...and convert the output to a table
    seg_table = seg_cat.to_table()
    radec_data = seg_table['sky_centroid_icrs']
    ra_icrs = radec_data.ra.degree
    dec_icrs = radec_data.dec.degree

    # RA and Dec are decimal degrees ICRS
    del seg_table['sky_centroid']
    del seg_table['sky_centroid_icrs']
    seg_table['RA_icrs'] = ra_icrs
    seg_table['Dec_icrs'] = dec_icrs

    # Rename column names to be more consistent with Sextractor column names
    seg_table.rename_column('xcentroid', 'x_image')
    seg_table.rename_column('ycentroid', 'y_image')
    seg_table.rename_column('background_at_centroid', 'background')
    seg_table.rename_column('source_sum', 'flux')
    seg_table.rename_column('source_sum_err', 'flux_err')
    seg_table.rename_column('cxx', 'cxx_image')
    seg_table.rename_column('cxy', 'cxy_image')
    seg_table.rename_column('cyy', 'cyy_image')
    seg_table.rename_column('covar_sigx2', 'x2_image')
    seg_table.rename_column('covar_sigy2', 'y2_image')
    seg_table.rename_column('covar_sigxy', 'xy_image')

    # Write out the official total detection product source catalog
    seg_table['x_image'].info.format = '.10f'
    seg_table['y_image'].info.format = '.10f'
    seg_table['RA_icrs'].info.format = '.10f'
    seg_table['Dec_icrs'].info.format = '.10f'
    print('seg_table (filter): {}'.format(seg_table))

    seg_table.write(catalog_filename, format='ascii.commented_header')
    log.info("Wrote filter source catalog: {}".format(catalog_filename))


def compute_background(image, threshold=None):
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

def modify_segmentation_map(image, segm, kernel):
    cat = source_properties(image, segm, filter_kernel=kernel)
    if len(cat) > 0:
        print('Original length of catalog: ', segm.nlabels)
        # Remove likely cosmic-rays based on central_moments classification
        bad_srcs = np.where(classify_sources(cat) == 0)[0] + 1

        if LooseVersion(photutils.__version__) >= '0.7':
            segm.remove_labels(bad_srcs)
        else:
            # this is the photutils >= 0.7 fast code for removing labels
            segm.check_labels(bad_srcs)
            bad_srcs = np.atleast_1d(bad_srcs)
            if len(bad_srcs) != 0:
                idx = np.zeros(segm.max_label + 1, dtype=int)
                idx[segm.labels] = segm.labels
                idx[bad_srcs] = 0
                segm.data = idx[segm.data]

        print('Modified length of catalog: ', segm.nlabels)
    return segm

def classify_sources(catalog, sources=None):
    """ Convert moments_central attribute for source catalog into star/cr flag.

    This algorithm interprets the central_moments from the source_properties
    generated for the sources as more-likely a star or a cosmic-ray.  It is not
    intended or expected to be precise, merely a means of making a first cut at
    removing likely cosmic-rays or other artifacts.

    Parameters
    ----------
    catalog : `~photutils.SourceCatalog`
        The photutils catalog for the image/chip.

    sources : tuple
        Range of objects from catalog to process as a tuple of (min, max).
        If None (default) all sources are processed.

    Returns
    -------
    srctype : ndarray
        An ndarray where a value of 1 indicates a likely valid, non-cosmic-ray
        source, and a value of 0 indicates a likely cosmic-ray.
    """
    moments = catalog.moments_central
    if sources is None:
        sources = (0, len(moments))
    num_sources = sources[1] - sources[0]
    srctype = np.zeros((num_sources,), np.int32)
    for src in range(sources[0], sources[1]):
        # Protect against spurious detections
        src_x = catalog[src].xcentroid
        src_y = catalog[src].ycentroid
        if np.isnan(src_x) or np.isnan(src_y):
            continue
        x, y = np.where(moments[src] == moments[src].max())
        if (x[0] > 1) and (y[0] > 1):
            srctype[src] = 1

    return srctype

def run_photutils():

    white_light_filename = "j92c01010_drc.fits"
    tdp_catalog_filename = "hst_10265_01s_acs_wfc_total_j92c01_segment-cat.ecsv"
    fp_filename_1 = "hst_10265_01s_acs_wfc_f606w_j92c01_drc.fits"
    fp_catalog_filename_1 = "hst_10265_01s_acs_wfc_f606w_j92c01_segment-cat.ecsv"
    fp_filename_2 = "hst_10265_01s_acs_wfc_f888w_j92c01_drc.fits"
    fp_catalog_filename_2 = "hst_10265_01s_acs_wfc_f888w_j92c01_segment-cat.ecsv"

    segmap, kernel, bkg, bkg_rms, bkg_dao_rms = create_sextractor_like_sourcelists(white_light_filename,
                                                                                   tdp_catalog_filename,
                                                                                   se_debug=True)

    # measure_source_properties(segmap, source_filename, catalog_filename, param_dict)
    measure_source_properties(segmap, kernel, fp_filename_1, fp_catalog_filename_1)
    print("Measured filter 1")

    measure_source_properties(segmap, kernel, fp_filename_2, fp_catalog_filename_2)
    print("Measured filter 2")

    return
