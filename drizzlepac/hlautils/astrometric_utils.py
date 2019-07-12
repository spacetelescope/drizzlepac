"""Utilities to support creation of astrometrically accurate reference catalogs

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
from io import BytesIO
import csv
import requests
import inspect
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
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.nddata.bitmask import bitfield_to_boolean_mask
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize

import photutils  # needed to check version
from photutils import detect_sources, source_properties, deblend_sources
from photutils import Background2D, MedianBackground
from photutils import DAOStarFinder
from tweakwcs import FITSWCS
from stwcs.distortion import utils
from stwcs import wcsutil
from stsci.tools import fileutil as fu
from stsci.tools import parseinput
from stsci.tools import logutil
from stsci.tools.fileutil import countExtn
import pysynphot as S

from ..tweakutils import build_xy_zeropoint

__taskname__ = 'astrometric_utils'

log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout)


ASTROMETRIC_CAT_ENVVAR = "ASTROMETRIC_CATALOG_URL"
DEF_CAT_URL = 'http://gsss.stsci.edu/webservices'

if ASTROMETRIC_CAT_ENVVAR in os.environ:
    SERVICELOCATION = os.environ[ASTROMETRIC_CAT_ENVVAR]
else:
    SERVICELOCATION = DEF_CAT_URL

MODULE_PATH = os.path.dirname(inspect.getfile(inspect.currentframe()))
VEGASPEC = os.path.join(os.path.dirname(MODULE_PATH),
                        'data', 'alpha_lyr_stis_008.fits')

__all__ = ['build_reference_wcs', 'create_astrometric_catalog', 'compute_radius',
           'find_gsc_offset', 'get_catalog',
           'extract_sources', 'find_hist2d_offset', 'generate_source_catalog',
           'classify_sources']

"""

Primary function for creating an astrometric reference catalog.

"""


def create_astrometric_catalog(inputs, catalog="GAIADR2", output="ref_cat.ecsv",
                               gaia_only=False, table_format="ascii.ecsv",
                               existing_wcs=None):
    """Create an astrometric catalog that covers the inputs' field-of-view.

    Parameters
    ----------
    input : str, list
        Filenames of images to be aligned to astrometric catalog

    catalog : str, optional
        Name of catalog to extract astrometric positions for sources in the
        input images' field-of-view. Default: GAIADR2. Options available are
        documented on the catalog web page.

    output : str, optional
        Filename to give to the astrometric catalog read in from the master
        catalog web service.  If None, no file will be written out.

    gaia_only : bool, optional
        Specify whether or not to only use sources from GAIA in output catalog

    existing_wcs : `~stwcs.wcsutil.HSTWCS`
        existing WCS object specified by the user

    Notes
    -----
    This function will point to astrometric catalog web service defined
    through the use of the ASTROMETRIC_CATALOG_URL environment variable.

    Returns
    -------
    ref_table : `~astropy.table.Table`
        Astropy Table object of the catalog

    """
    inputs, _ = parseinput.parseinput(inputs)

    # start by creating a composite field-of-view for all inputs
    # This default output WCS will have the same plate-scale and orientation
    # as the first chip in the list, which for WFPC2 data means the PC.
    # Fortunately, for alignment, this doesn't matter since no resampling of
    # data will be performed
    if existing_wcs is not None:
        outwcs = existing_wcs
    else:
        outwcs = build_reference_wcs(inputs)
    radius = compute_radius(outwcs)
    ra, dec = outwcs.wcs.crval

    # perform query for this field-of-view
    ref_dict = get_catalog(ra, dec, sr=radius, catalog=catalog)
    colnames = ('ra', 'dec', 'mag', 'objID', 'GaiaID')
    col_types = ('f8', 'f8', 'f4', 'U25', 'U25')
    ref_table = Table(names=colnames, dtype=col_types)

    # Add catalog name as meta data
    ref_table.meta['catalog'] = catalog
    ref_table.meta['gaia_only'] = gaia_only

    # rename coordinate columns to be consistent with tweakwcs
    ref_table.rename_column('ra', 'RA')
    ref_table.rename_column('dec', 'DEC')

    # extract just the columns we want...
    num_sources = 0
    for source in ref_dict:
        if 'GAIAsourceID' in source:
            g = source['GAIAsourceID']
            if gaia_only and g.strip() == '':
                continue
        else:
            g = "-1"  # indicator for no source ID extracted
        r = float(source['ra'])
        d = float(source['dec'])
        m = -999.9  # float(source['mag'])
        o = source['objID']
        num_sources += 1
        ref_table.add_row((r, d, m, o, g))

    # Write out table to a file, if specified
    if output:
        ref_table.write(output, format=table_format)
        log.info("Created catalog '{}' with {} sources".format(output, num_sources))

    return ref_table


def build_reference_wcs(inputs, sciname='sci'):
    """Create the reference WCS based on all the inputs for a field"""
    # start by creating a composite field-of-view for all inputs
    wcslist = []
    for img in inputs:
        nsci = countExtn(img)
        for num in range(nsci):
            extname = (sciname, num + 1)
            if sciname == 'sci':
                extwcs = wcsutil.HSTWCS(img, ext=extname)
            else:
                # Working with HDRLET as input and do the best we can...
                extwcs = read_hlet_wcs(img, ext=extname)

            wcslist.append(extwcs)

    # This default output WCS will have the same plate-scale and orientation
    # as the first chip in the list, which for WFPC2 data means the PC.
    # Fortunately, for alignment, this doesn't matter since no resampling of
    # data will be performed
    outwcs = utils.output_wcs(wcslist)

    return outwcs


def get_catalog(ra, dec, sr=0.1, catalog='GSC241'):
    """ Extract catalog from VO web service.

    Parameters
    ----------
    ra : float
        Right Ascension (RA) of center of field-of-view (in decimal degrees)

    dec : float
        Declination (Dec) of center of field-of-view (in decimal degrees)

    sr : float, optional
        Search radius (in decimal degrees) from field-of-view center to use
        for sources from catalog.  Default: 0.1 degrees

    catalog : str, optional
        Name of catalog to query, as defined by web-service.  Default: 'GSC241'

    Returns
    -------
    csv : CSV object
        CSV object of returned sources with all columns as provided by catalog

    """
    serviceType = 'vo/CatalogSearch.aspx'
    spec_str = 'RA={}&DEC={}&SR={}&FORMAT={}&CAT={}&MINDET=5'
    headers = {'Content-Type': 'text/csv'}
    fmt = 'CSV'

    spec = spec_str.format(ra, dec, sr, fmt, catalog)
    serviceUrl = '{}/{}?{}'.format(SERVICELOCATION, serviceType, spec)
    rawcat = requests.get(serviceUrl, headers=headers)
    r_contents = rawcat.content.decode()  # convert from bytes to a String
    rstr = r_contents.split('\r\n')
    # remove initial line describing the number of sources returned
    # CRITICAL to proper interpretation of CSV data
    del rstr[0]
    r_csv = csv.DictReader(rstr)

    return r_csv


def compute_radius(wcs):
    """Compute the radius from the center to the furthest edge of the WCS."""

    ra, dec = wcs.wcs.crval
    img_center = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)
    wcs_foot = wcs.calc_footprint()
    img_corners = SkyCoord(ra=wcs_foot[:, 0] * u.degree,
                           dec=wcs_foot[:, 1] * u.degree)
    radius = img_center.separation(img_corners).max().value

    return radius


def find_gsc_offset(image, input_catalog='GSC1', output_catalog='GAIA'):
    """Find the GSC to GAIA offset based on guide star coordinates

    Parameters
    ----------
    image : str
        Filename of image to be processed.

    Returns
    -------
    delta_ra, delta_dec : tuple of floats
        Offset in decimal degrees of image based on correction to guide star
        coordinates relative to GAIA.
    """
    serviceType = "GSCConvert/GSCconvert.aspx"
    spec_str = "TRANSFORM={}-{}&IPPPSSOOT={}"

    if 'rootname' in fits.getheader(image):
        ippssoot = fits.getval(image, 'rootname').upper()
    else:
        ippssoot = fu.buildNewRootname(image).upper()

    spec = spec_str.format(input_catalog, output_catalog, ippssoot)
    serviceUrl = "{}/{}?{}".format(SERVICELOCATION, serviceType, spec)
    rawcat = requests.get(serviceUrl)
    if not rawcat.ok:
        log.info("Problem accessing service with:\n{{}".format(serviceUrl))
        raise ValueError

    delta_ra = delta_dec = None
    tree = BytesIO(rawcat.content)
    for _, element in etree.iterparse(tree):
        if element.tag == 'deltaRA':
            delta_ra = float(element.text)
        elif element.tag == 'deltaDEC':
            delta_dec = float(element.text)

    return delta_ra, delta_dec


def build_auto_kernel(imgarr, whtarr, fwhm=3.0, threshold=None, source_box=7,
                      isolation_size=50, saturation_limit=70000.):
    """Build kernel for use in source detection based on image PSF
    This algorithm looks for an isolated point-source that is non-saturated to use as a template
    for the source detection kernel.  Failing to find any suitable sources, it will return a
    Gaussian2DKernel based on the provided FWHM as a default.
    Parameters
    ----------
    imgarr : ndarray
        Image array (ndarray object) with sources to be identified
    fwhm : float
        Value of FWHM to use for creating a Gaussian2DKernel object in case no suitable source
        can be identified in the image.
    threshold : float
        Value from the image which serves as the limit for determining sources.
        If None, compute a default value of (background+5*rms(background)).
        If threshold < 0.0, use absolute value as scaling factor for default value.
    source_box : int
        Size of box (in pixels) which defines the minimum size of a valid source.
    isolation_size : int
        Separation (in pixels) to use to identify sources that are isolated from any other sources
        in the image.
    saturation_limit : float
        Flux in the image that represents the onset of saturation for a pixel.
    Notes
    ------
    Ideally, it would be best to determine the saturation_limit value from the data itself,
    perhaps by looking at the pixels flagged (in the DQ array) as saturated and selecting
    the value less than the minimum flux of all those pixels, or maximum pixel value in the
    image if non-were flagged as saturated (in the DQ array).
    """

    # Try to use PSF derived from image as detection kernel
    # Kernel must be derived from well-isolated sources not near the edge of the image
    kern_img = imgarr.copy()
    edge = source_box * 2
    kern_img[:edge, :] = 0.0
    kern_img[-edge:, :] = 0.0
    kern_img[:, :edge] = 0.0
    kern_img[:, -edge:] = 0.0
    peaks = photutils.detection.find_peaks(kern_img, threshold=threshold, box_size=isolation_size)
    if len(peaks['peak_value'][peaks['peak_value'] > saturation_limit]):
        # Make sure only peaks less than saturation limit are evaluated
        peaks['peak_value'][peaks['peak_value'] >= saturation_limit] = 0.
    # Sort based on peak_value to identify brightest sources for use as a kernel
    peaks.sort('peak_value')

    wht_box = 2 # Weight image cutout box size is 2 x wht_box + 1 pixels on a side

    # Identify position of brightest, non-saturated peak (in numpy index order)
    for peak_ctr in range(-1, -1 * len(peaks) - 1, -1):
        log.info(peak_ctr)
        kernel_pos = [peaks['y_peak'][peak_ctr], peaks['x_peak'][peak_ctr]]

        kernel = imgarr[kernel_pos[0] - source_box:kernel_pos[0] + source_box + 1,
                        kernel_pos[1] - source_box:kernel_pos[1] + source_box + 1]

        kernel_wht = whtarr[kernel_pos[0] - wht_box:kernel_pos[0] + wht_box + 1,
                        kernel_pos[1] - wht_box:kernel_pos[1] + wht_box + 1]

        # search square cut-out (of size 2 x wht_box + 1 pixels on a side) of weight image centered on peak coords for
        # zero-value pixels. Reject peak if any are found.
        if len(np.where(kernel_wht == 0.)[0]) == 0:
            log.info("Kernel source PSF located at [{},{}]".format(kernel_pos[1], kernel_pos[0]))
            break
        else:
            kernel[:] = 0.0

    if kernel.sum() > 0.0:
        kernel /= kernel.sum() # Normalize the new kernel to a total flux of 1.0
    else:
        # Generate a default kernel using a simple 2D Gaussian
        sigma = fwhm * gaussian_fwhm_to_sigma
        kernel = Gaussian2DKernel(sigma, x_size=source_box, y_size=source_box)
        kernel.normalize()

    return kernel


def extract_sources(img, dqmask=None, fwhm=3.0, threshold=None, source_box=7,
                    classify=True, centering_mode="starfind", nlargest=None,
                    outroot=None, plot=False, vmax=None, deblend=False):
    """Use photutils to find sources in image based on segmentation.

    Parameters
    ----------
    img : ndarray
        Numpy array of the science extension from the observations FITS file.
    dqmask : ndarray
        Bitmask which identifies whether a pixel should be used (1) in source
        identification or not(0). If provided, this mask will be applied to the
        input array prior to source identification.
    fwhm : float
        Full-width half-maximum (fwhm) of the PSF in pixels.
    threshold : float or None
        Value from the image which serves as the limit for determining sources.
        If None, compute a default value of (background+5*rms(background)).
        If threshold < 0.0, use absolute value as scaling factor for default value.
    source_box : int
        Size of box (in pixels) which defines the minimum size of a valid source.
    classify : bool
        Specify whether or not to apply classification based on invarient moments
        of each source to determine whether or not a source is likely to be a
        cosmic-ray, and not include those sources in the final catalog.
    centering_mode : str
        "segmentaton" or "starfind"
        Algorithm to use when computing the positions of the detected sources.
        Centering will only take place after `threshold` has been determined, and
        sources are identified using segmentation.  Centering using `segmentation`
        will rely on `photutils.segmentation.source_properties` to generate the
        properties for the source catalog.  Centering using `starfind` will use
        `photutils.IRAFStarFinder` to characterize each source in the catalog.
    nlargest : int, None
        Number of largest (brightest) sources in each chip/array to measure
        when using 'starfind' mode.
    outroot : str, optional
        If specified, write out the catalog of sources to the file with this name rootname.
    plot : bool, optional
        Specify whether or not to create a plot of the sources on a view of the image.
    vmax : float, optional
        If plotting the sources, scale the image to this maximum value.
    deblend : bool, optional
        Specify whether or not to apply photutils deblending algorithm when
        evaluating each of the identified segments (sources) from the chip.
    """
    # apply any provided dqmask for segmentation only
    if dqmask is not None:
        imgarr = img.copy()
        imgarr[dqmask] = 0
    else:
        imgarr = img

    bkg_estimator = MedianBackground()
    bkg = None

    exclude_percentiles = [10, 25, 50, 75]
    for percentile in exclude_percentiles:
        try:
            bkg = Background2D(imgarr, (50, 50), filter_size=(3, 3),
                               bkg_estimator=bkg_estimator,
                               exclude_percentile=percentile)
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
        bkg_rms_mean = max(0.01, imgarr.min())
        bkg_rms = bkg_rms_mean * 5

    sigma = fwhm * gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma, x_size=source_box, y_size=source_box)
    kernel.normalize()
    segm = detect_sources(imgarr, threshold, npixels=source_box,
                          filter_kernel=kernel)
    # photutils >= 0.7: segm=None; photutils < 0.7: segm.nlabels=0
    if segm is None or segm.nlabels == 0:
        log.info("No detected sources!")
        return None, None

    if deblend:
        segm = deblend_sources(imgarr, segm, npixels=5,
                               filter_kernel=kernel, nlevels=16,
                               contrast=0.01)
    # If classify is turned on, it should modify the segmentation map
    if classify:
        cat = source_properties(imgarr, segm)
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

    # convert segm to mask for daofind
    if centering_mode == 'starfind':
        src_table = None
        # daofind = IRAFStarFinder(fwhm=fwhm, threshold=5.*bkg.background_rms_median)
        log.info("Setting up DAOStarFinder with: \n    fwhm={}  threshold={}".format(fwhm, bkg_rms_mean))
        daofind = DAOStarFinder(fwhm=fwhm, threshold=bkg_rms_mean)

        # Identify nbrightest/largest sources
        if nlargest is not None:
            nlargest = min(nlargest, len(segm.labels))
            if LooseVersion(photutils.__version__) >= '0.7':
                large_labels = segm.labels[
                    np.flip(np.argsort(segm.areas))[: nlargest]]
            else:
                # for photutils < 0.7
                areas = np.array([area for area in np.bincount(segm.data.ravel())[1:] if area != 0])
                large_labels = segm.labels[
                    np.flip(np.argsort(areas))[: nlargest]]

        log.info("Looking for sources in {} segments".format(len(segm.labels)))

        for segment in segm.segments:
            # check needed for photutils <= 0.6; it can be removed when
            # the drizzlepac depends on photutils >= 0.7
            if segment is None:
                continue
            if nlargest is not None and segment.label not in large_labels:
                continue  # Move on to the next segment
            # Get slice definition for the segment with this label
            seg_slice = segment.slices
            seg_yoffset = seg_slice[0].start
            seg_xoffset = seg_slice[1].start

            # Define raw data from this slice
            detection_img = img[seg_slice]
            # zero out any pixels which do not have this segments label
            detection_img[segm.data[seg_slice] == 0] = 0

            # Detect sources in this specific segment
            seg_table = daofind.find_stars(detection_img)

            # Pick out brightest source only
            if src_table is None and seg_table:
                # Initialize final master source list catalog
                src_table = Table(names=seg_table.colnames,
                                  dtype=[dt[1] for dt in seg_table.dtype.descr])

            if seg_table and seg_table['peak'].max() == detection_img.max():
                max_row = np.where(seg_table['peak'] == seg_table['peak'].max())[0][0]
                # Add row for detected source to master catalog
                # apply offset to slice to convert positions into full-frame coordinates
                seg_table['xcentroid'] += seg_xoffset
                seg_table['ycentroid'] += seg_yoffset
                src_table.add_row(seg_table[max_row])

    else:
        cat = source_properties(img, segm)
        src_table = cat.to_table()
        # Make column names consistent with IRAFStarFinder column names
        src_table.rename_column('source_sum', 'flux')
        src_table.rename_column('source_sum_err', 'flux_err')

    if src_table is not None:
        log.info("Total Number of detected sources: {}".format(len(src_table)))
    else:
        log.info("No detected sources!")
        return None, None

    # Move 'id' column from first to last position
    # Makes it consistent for remainder of code
    cnames = src_table.colnames
    cnames.append(cnames[0])
    del cnames[0]
    tbl = src_table[cnames]

    if outroot:
        tbl['xcentroid'].info.format = '.10f'  # optional format
        tbl['ycentroid'].info.format = '.10f'
        tbl['flux'].info.format = '.10f'
        if not outroot.endswith('.cat'):
            outroot += '.cat'
        tbl.write(outroot, format='ascii.commented_header')
        log.info("Wrote source catalog: {}".format(outroot))

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
    return tbl, segm


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


def generate_source_catalog(image, dqname="DQ", output=False, fwhm=3.0, **detector_pars):
    """ Build source catalogs for each chip using photutils.

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
    dqname : str
        EXTNAME for the DQ array, if present, in
        the input image HDUList.
    output : bool
        Specify whether or not to write out a separate catalog file for all the
        sources found in each chip.
    fwhm : float
        Full-width half-maximum (fwhm) of the PSF in pixels.

    Returns
    -------
    source_cats : dict
        Dict of astropy Tables identified by chip number with
        each table containing sources from image extension ``('sci', chip)``.

    """
    if not isinstance(image, fits.HDUList):
        raise ValueError("Input {} not fits.HDUList object".format(image))

    # remove parameters that are not needed by subsequent functions
    del detector_pars['fwhmpsf']

    # Build source catalog for entire image
    source_cats = {}
    numSci = countExtn(image, extname='SCI')
    outroot = None

    for chip in range(numSci):
        chip += 1
        # find sources in image
        if output:
            rootname = image[0].header['rootname']
            outroot = '{}_sci{}_src'.format(rootname, chip)

        imgarr = image['sci', chip].data

        # apply any DQ array, if available
        dqmask = None
        if image.index_of(dqname):
            dqarr = image[dqname, chip].data

            # "grow out" regions in DQ mask flagged as saturated by several
            # pixels in every direction to prevent the
            # source match algorithm from trying to match multiple sources
            # from one image to a single source in the
            # other or vice-versa.
            # Create temp DQ mask containing all pixels flagged with any value EXCEPT 256
            non_sat_mask = bitfield_to_boolean_mask(dqarr, ignore_flags=256)

            # Create temp DQ mask containing saturated pixels ONLY
            sat_mask = bitfield_to_boolean_mask(dqarr, ignore_flags=~256)

            # Grow out saturated pixels by a few pixels in every direction
            grown_sat_mask = ndimage.binary_dilation(sat_mask, iterations=5)

            # combine the two temporary DQ masks into a single composite DQ mask.
            dqmask = np.bitwise_or(non_sat_mask, grown_sat_mask)

            # dqmask = bitfield_to_boolean_mask(dqarr, good_mask_value=False)
            # TODO: <---Remove this old no-sat bit grow line once this
            # thing works

        seg_tab, segmap = extract_sources(imgarr, dqmask=dqmask, outroot=outroot, fwhm=fwhm, **detector_pars)
        seg_tab_phot = seg_tab

        source_cats[chip] = seg_tab_phot

    return source_cats


def generate_sky_catalog(image, refwcs, dqname="DQ", output=False):
    """Build source catalog from input image using photutils.

    This script borrows heavily from build_source_catalog.

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
        Input image.
    refwcs : `~stwcs.wcsutil.HSTWCS`
        Definition of the reference frame WCS.
    dqname : str, optional
        EXTNAME for the DQ array, if present, in the input image.
    output : bool, optional
        Specify whether or not to write out a separate catalog file for all the
        sources found in each chip.

    Returns
    --------
    master_cat : `~astropy.table.Table`
        Source catalog for all 'valid' sources identified from all chips of the
        input image with positions translated to the reference WCS coordinate
        frame.

    """
    # Extract source catalogs for each chip
    source_cats = generate_source_catalog(image, dqname=dqname, output=output)

    # Build source catalog for entire image
    master_cat = None
    numSci = countExtn(image, extname='SCI')
    # if no refwcs specified, build one now...
    if refwcs is None:
        refwcs = build_reference_wcs([image])
    for chip in range(numSci):
        chip += 1
        # work with sources identified from this specific chip
        seg_tab_phot = source_cats[chip]
        if seg_tab_phot is None:
            continue
        # Convert pixel coordinates from this chip to sky coordinates
        chip_wcs = wcsutil.HSTWCS(image, ext=('sci', chip))
        seg_ra, seg_dec = chip_wcs.all_pix2world(seg_tab_phot['xcentroid'], seg_tab_phot['ycentroid'], 1)
        # Convert sky positions to pixel positions in the reference WCS frame
        seg_xy_out = refwcs.all_world2pix(seg_ra, seg_dec, 1)
        seg_tab_phot['xcentroid'] = seg_xy_out[0]
        seg_tab_phot['ycentroid'] = seg_xy_out[1]
        if master_cat is None:
            master_cat = seg_tab_phot
        else:
            master_cat = vstack([master_cat, seg_tab_phot])

    return master_cat


def compute_photometry(catalog, photmode):
    """ Compute magnitudes for sources from catalog based on observations photmode.

    Parameters
    ----------
    catalog : `~astropy.table.Table`
        Astropy Table with 'source_sum' column for the measured flux for each source.

    photmode : str
        Specification of the observation filter configuration used for the exposure
        as reported by the 'PHOTMODE' keyword from the PRIMARY header.

    Returns
    -------
    phot_cat : `~astropy.table.Table`
        Astropy Table object of input source catalog with added column for
        VEGAMAG photometry (in magnitudes).
    """
    # Determine VEGAMAG zero-point using pysynphot for this photmode
    photmode = photmode.replace(' ', ', ')
    vega = S.FileSpectrum(VEGASPEC)
    bp = S.ObsBandpass(photmode)
    vegauvis = S.Observation(vega, bp)
    vegazpt = 2.5 * np.log10(vegauvis.countrate())

    # Use zero-point to convert flux values from catalog into magnitudes
    # source_phot = vegazpt - 2.5*np.log10(catalog['source_sum'])
    source_phot = vegazpt - 2.5 * np.log10(catalog['flux'])
    source_phot.name = 'vegamag'
    # Now add this new column to the catalog table
    catalog.add_column(source_phot)

    return catalog


def filter_catalog(catalog, bright_limit=1.0, max_bright=None, min_bright=20, colname="vegamag"):
    """ Create a new catalog selected from input based on photometry.

    Parameters
    ----------
    catalog : `~astropy.table.Table`
        Table containing the full set of identified sources.
    bright_limit : float
        Fraction of catalog based on brightness that should be retained.
        Value of 1.00 means full catalog.
    max_bright : int
        Maximum number of sources to keep regardless of `bright_limit`.
    min_bright : int
        Minimum number of sources to keep regardless of `bright_limit`.
    colname : str
        Name of column to use for selection/sorting.

    Returns
    -------
    new_catalog : `~astropy.table.Table`
        New table which only has the sources that meet the selection criteria.
    """

    # sort by magnitude
    phot_column = catalog[colname]
    num_sources = len(phot_column)
    sort_indx = np.argsort(phot_column)
    if max_bright is None:
        max_bright = num_sources

    # apply limits, insuring no more than full catalog gets selected
    limit_num = max(int(num_sources * bright_limit), min_bright)
    limit_num = min(max_bright, limit_num, num_sources)

    # Extract sources identified by selection
    new_catalog = catalog[sort_indx[:limit_num]]

    return new_catalog


def build_self_reference(filename, clean_wcs=False):
    """ This function creates a reference, undistorted WCS that can be used to
    apply a correction to the WCS of the input file.

    Parameters
    ----------
    filename : str
        Filename of image which will be corrected, and which will form the basis
        of the undistorted WCS.

    clean_wcs : bool
        Specify whether or not to return the WCS object without any distortion
        information, or any history of the original input image.  This converts
        the output from `utils.output_wcs()` into a pristine `~stwcs.wcsutil.HSTWCS` object.

    Returns
    -------
    customwcs : `~stwcs.wcsutil.HSTWCS`
        HSTWCS object which contains the undistorted WCS representing the entire
        field-of-view for the input image.

    Examples
    --------
    This function can be used with the following syntax to apply a shift/rot/scale
    change to the same image:

    >>> import buildref
    >>> from drizzlepac import updatehdr
    >>> filename = "jce501erq_flc.fits"
    >>> wcslin = buildref.build_self_reference(filename)
    >>> updatehdr.updatewcs_with_shift(filename, wcslin, xsh=49.5694,
    ... ysh=19.2203, rot = 359.998, scale = 0.9999964)

    """
    if 'sipwcs' in filename:
        sciname = 'sipwcs'
    else:
        sciname = 'sci'

    wcslin = build_reference_wcs([filename], sciname=sciname)

    if clean_wcs:
        wcsbase = wcslin.wcs
        customwcs = build_hstwcs(wcsbase.crval, wcsbase.crpix,
                                 wcslin.naxis1, wcslin.naxis2,
                                 wcslin.pscale, wcslin.orientat)
    else:
        customwcs = wcslin
    return customwcs


def read_hlet_wcs(filename, ext):
    """Insure `~stwcs.wcsutil.HSTWCS` includes all attributes of a full image WCS.

    For headerlets, the WCS does not contain information about the size of the
    image, as the image array is not present in the headerlet.
    """
    hstwcs = wcsutil.HSTWCS(filename, ext=ext)
    if hstwcs.naxis1 is None:
        hstwcs.naxis1 = int(hstwcs.wcs.crpix[0] * 2.)  # Assume crpix is center of chip
        hstwcs.naxis2 = int(hstwcs.wcs.crpix[1] * 2.)

    return hstwcs


def build_hstwcs(crval, crpix, naxis1, naxis2, pscale, orientat):
    """ Create an `~stwcs.wcsutil.HSTWCS` object for a default instrument without
    distortion based on user provided parameter values.
    """
    wcsout = wcsutil.HSTWCS()
    wcsout.wcs.crval = crval.copy()
    wcsout.wcs.crpix = crpix.copy()
    wcsout.naxis1 = naxis1
    wcsout.naxis2 = naxis2
    wcsout.wcs.cd = fu.buildRotMatrix(orientat) * [-1, 1] * pscale / 3600.0
    # Synchronize updates with astropy.wcs objects
    wcsout.wcs.set()
    wcsout.setPscale()
    wcsout.setOrient()
    wcsout.wcs.ctype = ['RA---TAN', 'DEC--TAN']

    return wcsout


def within_footprint(img, wcsobj, x, y):
    """Determine whether input x, y fall in the science area of the image.

    Parameters
    ----------
    img : ndarray
        ndarray of image where non-science areas are marked with value of NaN.

    wcsobj : `~stwcs.wcsutil.HSTWCS`
        HSTWCS or WCS object with naxis terms defined.

    x, y : ndarray
        arrays of x, y positions for sources to be checked.

    Returns
    -------
    x, y : ndarray
        New arrays which have been trimmed of all sources that fall outside
        the science areas of the image

    """
    # start with limits of WCS shape

    sky = wcsobj.pixel_to_world(x, y, 1)
    inmask = wcsobj.footprint_contains(sky)
    x = x[inmask]
    y = y[inmask]
    return x, y


def find_hist2d_offset(filename, reference, refwcs=None, refnames=['ra', 'dec'],
                       match_tolerance=5., chip_catalog=True, search_radius=15.0,
                       min_match=10, classify=True):
    """Iteratively look for the best cross-match between the catalog and ref.

    Parameters
    ----------
    filename : `~astropy.io.fits.HDUList` or str
        Single image to extract sources for matching to
        the external astrometric catalog.

    reference : str or `~astropy.table.Table`
        Reference catalog, either as a filename or ``astropy.Table``
        containing astrometrically accurate sky coordinates for astrometric
        standard sources.

    refwcs : `~stwcs.wcsutil.HSTWCS`
        This WCS will define the coordinate frame which will
        be used to determine the offset. If None is specified, use the
        WCS from the input image `filename` to build this WCS using
        `build_self_reference()`.

    refnames : list
        List of table column names for sky coordinates of astrometric
        standard sources from reference catalog.

    match_tolerance : float
        Tolerance (in pixels) for recognizing that a source position matches
        an astrometric catalog position.  Larger values allow for lower
        accuracy source positions to be compared to astrometric catalog

    chip_catalog : bool
        Specify whether or not to write out individual source catalog for
        each chip in the image.

    search_radius : float
        Maximum separation (in arcseconds) from source positions to look
        for valid cross-matches with reference source positions.

    min_match : int
        Minimum number of cross-matches for an acceptable determination of
        the offset.

    classify : bool
        Specify whether or not to use central_moments classification to
        ignore likely cosmic-rays/bad-pixels when generating the source
        catalog.

    Returns
    -------
        best_offset : tuple
            Offset in input image pixels between image source positions and
            astrometric catalog positions that results in largest number of
            matches of astrometric sources with image sources

        seg_xy, ref_xy : `~astropy.table.Table`
            Source catalog and reference catalog, respectively, used for
            determining the offset.  Each catalog includes sources for the entire
            field-of-view, not just a single chip.
    """
    # Interpret input image to generate initial source catalog and WCS
    if isinstance(filename, str):
        image = fits.open(filename)
        rootname = filename.split("_")[0]
    else:
        image = filename
        rootname = image[0].header['rootname']

    # check to see whether reference catalog can be found
    if not os.path.exists(reference):
        log.info("Could not find input reference catalog: {}".format(reference))
        raise FileNotFoundError

    # Extract reference WCS from image
    if refwcs is None:
        refwcs = build_self_reference(image, clean_wcs=True)
    log.info("Computing offset for field-of-view defined by:")
    log.info(refwcs)

    # read in reference catalog
    if isinstance(reference, str):
        refcat = ascii.read(reference)
    else:
        refcat = reference
    log.info("\nRead in reference catalog with {} sources.".format(len(refcat)))

    ref_ra = refcat[refnames[0]]
    ref_dec = refcat[refnames[1]]

    # Build source catalog for entire image
    img_cat = generate_source_catalog(image, refwcs, output=chip_catalog, classify=classify)
    img_cat.write(filename.replace(".fits", "_xy.cat"), format='ascii.no_header',
                  overwrite=True)

    # Retrieve source XY positions in reference frame
    seg_xy = np.column_stack((img_cat['xcentroid'], img_cat['ycentroid']))
    seg_xy = seg_xy[~np.isnan(seg_xy[:, 0])]

    # Translate reference catalog positions into input image coordinate frame
    xref, yref = refwcs.all_world2pix(ref_ra, ref_dec, 1)

    # look for only sources within the viewable area of the exposure to
    # determine the offset
    xref, yref = within_footprint(image, refwcs, xref, yref)
    ref_xy = np.column_stack((xref, yref))
    log.info("\nWorking with {} astrometric sources for this field".format(len(ref_xy)))

    # write out astrometric reference catalog that was actually used
    ref_ra_img, ref_dec_img = refwcs.all_pix2world(xref, yref, 1)
    ref_tab = Table([ref_ra_img, ref_dec_img, xref, yref], names=['ra', 'dec', 'x', 'y'])
    ref_tab.write(reference.replace('.cat', '_{}.cat'.format(rootname)),
                  format='ascii.fast_commented_header', overwrite=True)
    searchrad = search_radius / refwcs.pscale

    # Use 2d-Histogram builder from drizzlepac.tweakreg -- for demo only...
    xp, yp, nmatches, zpqual = build_xy_zeropoint(seg_xy, ref_xy,
                                                  searchrad=searchrad,
                                                  histplot=False, figure_id=1,
                                                  plotname=None, interactive=False)
    hist2d_offset = (xp, yp)
    log.info('best offset {} based on {} cross-matches'.format(hist2d_offset, nmatches))

    return hist2d_offset, seg_xy, ref_xy


##############################
#
# Functions to support working with Tweakwcs
#
##############################
def build_wcscat(image, group_id, source_catalog):
    """ Return a list of `~tweakwcs.tpwcs.FITSWCS` objects for all chips in an image.

    Parameters
    ----------
    image : str, `~astropy.io.fits.HDUList`
        Either filename or HDUList of a single HST observation.

    group_id : int
        Integer ID for group this image should be associated with; primarily
        used when separate chips are in separate files to treat them all as one
        exposure.

    source_catalog : dict
        If provided, these catalogs will be attached as `catalog`
        entries in each chip's ``FITSWCS`` object.  It should be provided as a
        dict of astropy Tables identified by chip number with
        each table containing sources from image extension ``('sci', chip)`` as
        generated by `generate_source_catalog()`.

    Returns
    -------
    wcs_catalogs : list of `~tweakwcs.tpwcs.FITSWCS`
        List of `~tweakwcs.tpwcs.FITSWCS` objects defined for all chips in input image.

    """
    open_file = False
    if isinstance(image, str):
        hdulist = fits.open(image)
        open_file = True
    elif isinstance(image, fits.HDUList):
        hdulist = image
    else:
        log.info("Wrong type of input, {}, for build_wcscat...".format(type(image)))
        raise ValueError

    wcs_catalogs = []
    numsci = countExtn(hdulist)
    for chip in range(1, numsci + 1):
        w = wcsutil.HSTWCS(hdulist, ('SCI', chip))

        imcat = source_catalog[chip]
        # rename xcentroid/ycentroid columns, if necessary, to be consistent with tweakwcs
        if 'xcentroid' in imcat.colnames:
            imcat.rename_column('xcentroid', 'x')
            imcat.rename_column('ycentroid', 'y')

        wcscat = FITSWCS(
            w,
            meta={
                'chip': chip,
                'group_id': group_id,
                'filename': image,
                'catalog': imcat,
                'name': image
            }
        )

        wcs_catalogs.append(wcscat)

    if open_file:
        hdulist.close()

    return wcs_catalogs
