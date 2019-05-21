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
import os
import sys
from distutils.version import LooseVersion

import numpy as np
from scipy import ndimage
from lxml import etree
try:
    from matplotlib import pyplot as plt
except Exception:
    plt = None

from astropy import units as u
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
from astropy.io import fits as fits
from astropy.io import ascii
from astropy.convolution import Gaussian2DKernel, MexicanHat2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.nddata.bitmask import bitfield_to_boolean_mask
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize

import photutils
from photutils import detect_sources, source_properties, deblend_sources
from photutils import Background2D, MedianBackground, SExtractorBackground
from tweakwcs import FITSWCS
from stwcs.distortion import utils
from stwcs import wcsutil
from stwcs.wcsutil import headerlet, HSTWCS
from stsci.tools import fileutil as fu
from stsci.tools import parseinput
from stsci.tools import logutil
from stsci.tools.fileutil import countExtn


__taskname__ = 'se_source_generation'

# Module-level dictionary contains instrument/detector-specific parameters used later on in the script.
detector_specific_params = {"acs": {"hrc": {"fwhmpsf": 0.152,  # 0.073
                                            "classify": True,
                                            "threshold": None},
                                    "sbc": {"fwhmpsf": 0.13,  # 0.065
                                            "classify": False,
                                            "threshold": 2.0},
                                    "wfc": {"fwhmpsf": 0.13,  # 0.076,
                                            "classify": True,
                                            "threshold": -1.1}},
                            "wfc3": {"ir": {"fwhmpsf": 0.25,  # 0.14
                                            "classify": False,
                                            "threshold": None},
                                     "uvis": {"fwhmpsf": 0.152,  # 0.076
                                              "classify": True,
                                              "threshold": None}}}

log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout)


#__all__ = ['run_photutils', 'create_sextractor_like_sourcelists', 'generate_se_catalog', 
           'classify_sources']

"""

Primary function for creating an astrometric reference catalog.
name: hst_<propid>_<obsetid>_<instr>_<detector>_total_<ipppss>_cat.txt
where is DAOPHOT or SEXTRACTOR?

"""

#def create_sextractor_like_sourcelists(totdet_product_cat_dict, param_dict):
def create_sextractor_like_sourcelists(imgarr, vmax=None, plt_debug=False):
    """Use photutils to find sources in image based on segmentation.

    Parameters
    ----------
    totdet_product_cat_dict : dictionary
        Dictionary which maps the image filename for a total detection product to the associate catalog filename

    param_dict : dictionary
        dictionary of drizzle, source finding, and photometric parameters

    Ximgarr : ndarray
        Numpy array of the science extension from the observations FITS file.

    vmax : float, optional
        If plotting the sources, scale the image to this maximum value.

    plt_debug : bool, optional
        Specify whether or not to plot the image and segmentation image for
        visualization and debugging purposes

    Returns
    -------
    return src_table, segm, bkg, bkg_rms, bkg_rms_mean
    src_table : `~astropy.table.QTable`
        Table where each row represents a source
  
    segm : `photutils.segmentation.SegmentationImage`
        Two-dimensional segmentation image where found source regions are labeled with 
        unique, non-zero positive integers.

    bkg : `~photutils.background.Background2D` or None 
        A background map based upon the `~photutils.background.SExtractorBackground` 
        estimator

    bkg_rms : float
        N times the bkg.background_rms where N = 5 FIX
        
    bkg_rms_mean : float
        Mean bkg.background FIX

    """

    # Get the image data
    filename = list(totdet_product_cat_dict.keys())[0]
    imgarr = fits.open(filename)
    imgwcs = HSTWCS(imgarr, 1)

    # Default parameters
    seParams_default = {}
    seParams_default['fwhm'] = 3.0
    seParams_default['size_source_box'] = 3.0
    seParams_default['threshold'] = None
    seParams = seParams_default
    
    fwhm = float(param_dict['se']['TWEAK_FWHMPSF'])
    size_source_box = float(param_dict['se']['TWEAK_SOURCE_BOX'])
    threshold = float(param_dict['se']['TWEAK_THRESHOLD'])

    seParams['fwhm'] = fwhm
    seParams['size_source_box'] = size_source_box
    seParams['threshold'] = threshold

    # Report configuration values to log
    log.info('====================')
    log.info('')
    log.info('SExtractor-like source finding settings for Photutils segmentation')
    log.info('FWHM: {}'.format(fwhm))
    log.info('size_source_box: {}'.format(size_source_box))
    log.info('threshold: {}'.format(threshold))
    log.info('')
    log.info('====================')

    # Only use a single kernel for now
    kernel_list = [Gaussian2DKernel, MexicanHat2DKernel]
    kernel_in_use = kernel_list[0]

    bkg, bkg_rms, bkg_rms_mean, threshold = compute_background(imgarr, threshold=threshold)

    # FIX imgarr should be background subtracted, sextractor uses the filtered_data image
    # Can input wcs to source_properties to have sky coords
    # See Photutils docs
    imgarr_bkgsub = imgarr - bkg.background

    # *** FIX: should size_source_box size be used in all these places? ***
    # Create a 2D filter kernel
    sigma = fwhm * gaussian_fwhm_to_sigma
    kernel = kernel_in_use(sigma, x_size=size_source_box, y_size=size_source_box)
    kernel.normalize()
  
    # Source segmentation/extraction
    # If the threshold includes the background level, then the input image
    # should NOT be background subtracted.
    print('Threshold: {}'.format(threshold))
    segm = detect_sources(imgarr, threshold, npixels=size_source_box,
                          filter_kernel=kernel)
    print(segm)

    # TESTING
    # Move 'id' column from first to last position
    # Makes it consistent for remainder of code
    # Overlay on image looks good - yeah.
    cat = source_properties(imgarr_bkgsub, segm, filter_kernel=kernel)
    table = cat.to_table()
    cnames = table.colnames
    cnames.append(cnames[0])
    del cnames[0]
    tbl = table[cnames[0:2]]

    outname = 'md_segm.reg'
    tbl['xcentroid'].info.format = '.10f'  # optional format
    tbl['ycentroid'].info.format = '.10f'
    tbl.write(outname, format='ascii.commented_header')
    log.info("Wrote source catalog: {}".format(outname))
    # END TESTING

    if plt_debug:
        norm = ImageNormalize(stretch=SqrtStretch())
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))
        ax1.imshow(imgarr, origin='lower', cmap='Greys_r', norm=norm)
        ax1.set_title('Data')
        ax2.imshow(segm, origin='lower', cmap=segm.cmap(random_state=12345))
        ax2.set_title('Segmentation Image')
        plt.show()

    # Deblending is a combination of multi-thresholding and watershed
    # segmentation. Sextractor uses a multi-thresholding technique.
    # npixels = number of connected pixels in source
    # npixels and filter_kernel should match those used by detect_sources()
    #segm = deblend_sources(imgarr, segm, npixels=size_source_box,
    #                       filter_kernel=kernel, nlevels=16,
    #                       contrast=0.01)
    #print('after deblend. ', segm)

    # If classify is turned on, it should modify the segmentation map
    # Should imgarr be imgarr_bkg here?  I think so!
    segm = modify_segmentation_map(imgarr_bkgsub, segm, kernel)

    # Regenerate the source catalog with presumably now only good sources
    src_cat = source_properties(imgarr_bkgsub, segm, filter_kernel=kernel)
    src_table = src_cat.to_table()
    print("segm: ", segm)
    print("tbl: ", src_table)

    return src_table, segm, bkg, bkg_rms, bkg_rms_mean


"""
#def measure_source_properties(totdet_product_cat_dict, filter_product_cat_dict,
#                                       inst_det, param_dict):
def measure_source_properties(imgarr, fwhm=3.0, threshold=None, size_source_box=7,
                       classify=True, vmax=None, deblend=True, plt_debug=False):

    # Regenerate the source catalog with presumably now good sources
    src_cat = source_properties(imgarr_bkgsub, segm, filter_kernel=kernel)

    # Collect properties for each identified source
    # FIX Make this a subroutine
    if len(src_cat) < 1:
        log.info("No detected sources!")
        return None
    else:
        src_table = src_cat.to_table()
        log.info("Total Number of detected sources via Photutils segmentation: {}".format(len(src_table)))
        #src_table['cxx'].info.format =  pixels**(-2)
        #src_table['cyy'].info.format =  pixels**(-2)
        #src_table['cxy'].info.format =  pixels**(-2)

        # Rename column names to be more consistent with Sextractor column names
        src_table.rename_column('source_sum', 'flux')
        src_table.rename_column('source_sum_err', 'flux_err')

        print(src_table)

    return src_table, segm
"""

def compute_background (image, threshold=None):
    #Previous: bkg_estimator = MedianBackground()
    bkg_estimator = SExtractorBackground()
    print('bkg estimator: ', bkg_estimator)
    #bkgrms_estimator = StdBackgroundRMS()
    bkg = None

    exclude_percentiles = [10, 25, 50, 75]
    for percentile in exclude_percentiles:
        log.info("Percentile in use: {}".format(percentile))
        try:
            bkg = Background2D(image, (50, 50), filter_size=(3, 3),
                               bkg_estimator=bkg_estimator,
                               exclude_percentile=percentile)
            print('bkg: ', bkg)
        except Exception:
            bkg = None
            continue

        if bkg is not None:
            # If it succeeds, stop and use that value
            bkg_rms = (5. * bkg.background_rms)
            bkg_rms_mean = bkg.background.mean() + 5. * bkg_rms.std()
            default_threshold = bkg.background + bkg_rms
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

    return bkg, bkg_rms, bkg_rms_mean, threshold

# do not get confused between SourceProperties the class and source_properties the method.
# cat = SourceProperties(imgarr_bkgsub, segm, filtered_data=?) for Sextractor
# FIX - also input wcs to get sky coords
def modify_segmentation_map(image, segm, kernel):
    cat = source_properties(image, segm, filter_kernel=kernel)
    if len(cat) > 0:
        print('modify. len(cat): ',len(cat))
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

# Main entry point to the SExtractor-like analysis and source catalog 
# generation
def generate_se_catalogs(image, dbg_output=False, **detector_pars):
    """ Build SExtractor-like source catalogs using photutils.

    This function is essentially a high-level controller for the generation
    of the SExtractor-like catalogs.  A "white light" image (*drz.fits or 
    *drc.fits) will be analyzed first in order to determine the positions of
    sources.  These positions will then be used in the analysis of the 
    "instrument/detector/filter" combined image to produce the corresponding
    source catalogs.

    The catalog returned by this function includes sources found in all chips
    of the input image with the positions translated to the coordinate frame
    defined by the reference WCS `refwcs`.  The sources will be
    - identified using photutils segmentation-based source finding code
    - ignore any input pixel which has been flagged as 'bad' in the DQ
    array, should a DQ array be found in the input HDUList.
    - classified as probable cosmic-rays (if enabled) using central_moments
    properties of each source, with these sources being removed from the
    catalog.

    Parameters
    ----------
    image : `~astropy.io.fits.HDUList`
        Input image as an astropy.io.fits HDUList.
    fwhm : float
        Full-width half-maximum (fwhm) of the PSF in pixels.
    dbg_output : bool, optional
        Specify whether or not to generate a "regions" catalog file which is
        compatible with being ingested by DS9 for overlay purposes.
    plot : bool, optional
        Specify whether or not to create a plot of the sources on a view of the image.

    Returns
    -------
    source_cats : dict
        Dict of astropy Tables identified by chip number with
        each table containing sources from image extension ``('sci', chip)``.

    """
    if not isinstance(image, fits.HDUList):
        raise ValueError("Input {} not fits.HDUList object".format(image))

    # Build source catalog for entire image
    source_cats = {}
    outname = None

    seg_table, segmap, bkg, bkg_rms, bkg_rms_mean = create_sextractor_like_sourcelists(imgarr, 
                                                                                           vmax=None, 
                                                                                           plt_debug=False)
        
        """
        # If requested, write out a catalog which can be used as an overlay region in DS9
        if dbg_output:
            # Move 'id' column from first to last position
            # Makes it consistent for remainder of code
            cnames = seg_table.colnames
            cnames.append(cnames[0])
            del cnames[0]
            tbl = seg_table[cnames[0:2]]

            outname = 'hst_{}_segm.reg'.format(image[0].header['rootname'])
            tbl['xcentroid'].info.format = '.10f'  # optional format
            tbl['ycentroid'].info.format = '.10f'
            #tbl['flux'].info.format = '.10f'
            tbl.write(outname, format='ascii.commented_header')
            log.info("Wrote source catalog: {}".format(outname))

        if plot and plt is not None:
            norm = len(segm.labels)
            if vmax is None:
                norm = ImageNormalize(stretch=SqrtStretch())
            fig, ax = plt.subplots(2, 2, figsize=(8, 8))
            ax[0][0].imshow(imgarr, origin='lower', cmap='Greys_r', norm=norm, vmax=vmax)
            ax[0][1].imshow(segm, origin='lower', cmap=segm.cmap(random_state=12345))
            ax[0][1].set_title('Segmentation Map')
            ax[1][0].imshow(bkg.background, origin='lower')
            if not isinstance(threshold, float):
                ax[1][1].imshow(threshold, origin='lower')

        # Write out the official source catalog
        seg_table['xcentroid'].info.format = '.10f'  # optional format
        seg_table['ycentroid'].info.format = '.10f'
        seg_table['flux'].info.format = '.10f'
        outroot = 'hst'
        if not outroot.endswith('_segment-cat.ecsv'):
            outroot += '_segment-cat.ecsv'
        seg_table.write(outroot, format='ascii.commented_header')
        log.info("Wrote source catalog: {}".format(outroot))

        seg_tab_phot = seg_table

        source_cats[chip] = seg_tab_phot

        #catalog_name = hst_<propid>_<obsetid>_<instr>_<detector>_total_<ipppss>_cat.ecsv 
        """

    return source_cats

# Driver with functions or actions cribbed from alignimages/generate_source_catalogs or 
# astrometric_utils/generate_source_catalog (generate_se_catalog in this file)
def run_photutils(imgname):
    imghdu = fits.open(imgname)

    imgprimaryheader = imghdu[0].header
    instrument = imgprimaryheader['INSTRUME'].lower()
    detector = imgprimaryheader['DETECTOR'].lower()

    sourcecatalogdict = {}
    sourcecatalogdict[imgname] = {}

    # get instrument/detector-specific image alignment parameters
    if instrument in detector_specific_params.keys():
        if detector in detector_specific_params[instrument].keys():
            detector_pars = detector_specific_params[instrument][detector]
            # to allow generate_source_catalog to get detector specific parameters
            #detector_pars.update(pars)
            sourcecatalogdict[imgname]["params"] = detector_pars
        else:
            sys.error("ERROR! Unrecognized detector '{}'. Exiting...".format(detector))
            log.exit("ERROR! Unrecognized detector '{}'. Exiting...".format(detector))
    else:
        sys.error("ERROR! Unrecognized instrument '{}'. Exiting...".format(instrument))
        log.exit("ERROR! Unrecognized instrument '{}'. Exiting...".format(instrument))

    # Identify sources in image, convert coords from chip x, y form to reference WCS sky RA, Dec form.
    imgwcs = HSTWCS(imghdu, 1)
    # Convert fwhmpsf from arsec to pixels
    fwhmpsf_pix = sourcecatalogdict[imgname]["params"]['fwhmpsf'] / imgwcs.pscale

    sourcecatalogdict[imgname]["catalog_table"] = \
        generate_se_catalogs(imghdu, fwhm=fwhmpsf_pix, dbg_output=True, **detector_pars)

    return
    
