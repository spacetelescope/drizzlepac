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
import requests
import inspect
import sys
import time
from distutils.version import LooseVersion
import copy

import numpy as np
import scipy.stats as st
from scipy import ndimage
from scipy.stats import pearsonr
from lxml import etree
try:
    from matplotlib import pyplot as plt
except Exception:
    plt = None

from astropy import units as u
from astropy.table import Table, vstack, Column
from astropy.coordinates import SkyCoord
from astropy.io import fits as fits
from astropy.io import ascii
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma, gaussian_sigma_to_fwhm, sigma_clipped_stats
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.time import Time
from astropy.utils.decorators import deprecated

import photutils  # needed to check version
if LooseVersion(photutils.__version__) < '1.1.0':
    OLD_PHOTUTILS = True
    from photutils.segmentation import (detect_sources,
                                        deblend_sources, make_source_mask)
    from photutils.segmentation import source_properties as SourceCatalog
else:
    OLD_PHOTUTILS = False
    from photutils.segmentation import (detect_sources, SourceCatalog,
                                        deblend_sources, make_source_mask,
                                        SegmentationImage)

from photutils.detection import DAOStarFinder, find_peaks

from photutils.background import (Background2D, MMMBackground,
                                  SExtractorBackground, StdBackgroundRMS)
from photutils.psf import (IntegratedGaussianPRF, DAOGroup,
                           IterativelySubtractedPSFPhotometry)

from tweakwcs import FITSWCS
from stwcs.distortion import utils
from stwcs import wcsutil
from stwcs.updatewcs import astrometry_utils
from stsci.tools import fileutil as fu
from stsci.tools import parseinput
from stsci.tools import logutil
from stsci.tools.fileutil import countExtn
from stsci.tools.bitmask import bitfield_to_boolean_mask

import stsci.tools

from ..tweakutils import build_xy_zeropoint, ndfind

__taskname__ = 'astrometric_utils'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)

ASTROMETRIC_CAT_ENVVAR = "ASTROMETRIC_CATALOG_URL"
DEF_CAT_URL = 'http://gsss.stsci.edu/webservices'

if ASTROMETRIC_CAT_ENVVAR in os.environ:
    SERVICELOCATION = os.environ[ASTROMETRIC_CAT_ENVVAR]
else:
    SERVICELOCATION = DEF_CAT_URL

MODULE_PATH = os.path.dirname(inspect.getfile(inspect.currentframe()))

__all__ = ['create_astrometric_catalog', 'compute_radius',
           'build_auto_kernel', 'find_fwhm',
           'get_catalog', 'get_catalog_from_footprint',
           'extract_sources', 'find_hist2d_offset', 'generate_source_catalog',
           'classify_sources', 'within_footprint',
           'compute_similarity', 'determine_focus_index', 'max_overlap_diff']

FOCUS_DICT = {'exp': [], 'prod': [], 'stats': {},
              'exp_pos': None, 'prod_pos': None,
              'alignment_verified': False, 'alignment_quality': -1,
              'expnames': "", 'prodname': ""}

EXP_LIMIT = 0.05  # hard-limit of exptime weighting for comparing images
EXP_RATIO = 0.2

# A radius of 25 pixels (1" in WFC3, 1.25" in ACS) corresponds to >=95% total flux for a point-source
# Any source larger than this, would either be saturated or blended with other sources
# in either case, going to a smaller kernel will help with source identification.
MAX_AREA_LIMIT = 1964

#
# Dictionary of supported astrometric catalogs with column name translations
# These are the minimum columns necessary for alignment to work properly while
# taking into account all available information from the catalogs.
# Blank column names indicate the catalog does not have that column,
# and will set to None in the output table.
# Output table using these common column names will be produced by
# `create_astrometric_catalog`.
#
SUPPORTED_CATALOGS = {
    'GAIAEDR3': {'RA': 'ra', 'RA_error': 'ra_error',
                 'DEC': 'dec', 'DEC_error': 'dec_error',
                 'pmra': 'pmra', 'pmra_error': 'pmra_error',
                 'pmdec': 'pmdec', 'pmdec_error': 'pmdec_error',
                 'mag': 'mag', 'objID': 'objID', 'epoch': 'epoch'},
    'GAIADR2': {'RA': 'ra', 'RA_error': 'ra_error',
                'DEC': 'dec', 'DEC_error': 'dec_error',
                'pmra': 'pmra', 'pmra_error': 'pmra_error',
                'pmdec': 'pmdec', 'pmdec_error': 'pmdec_error',
                'mag': 'mag', 'objID': 'objID', 'epoch': 'epoch'},
    'GAIADR1': {'RA': 'ra', 'RA_error': 'ra_error',
                'DEC': 'dec', 'DEC_error': 'dec_error',
                'pmra': 'pmra', 'pmra_error': 'pmra_error',
                'pmdec': 'pmdec', 'pmdec_error': 'pmdec_error',
                'mag': 'mag', 'objID': 'objID', 'epoch': 'epoch'},
    '2MASS': {'RA': 'ra', 'RA_error': '',
              'DEC': 'dec', 'DEC_error': '',
              'pmra': '', 'pmra_error': '',
              'pmdec': '', 'pmdec_error': '',
              'mag': 'mag', 'objID': 'objid', 'epoch': 'jdate'},
    'PS1DR2': {'RA': 'ra', 'RA_error': 'raMeanErr',
               'DEC': 'dec', 'DEC_error': 'decMeanErr',
               'pmra': '', 'pmra_error': '',
               'pmdec': '', 'pmdec_error': '',
               'mag': 'mag', 'objID': 'objid', 'epoch': 'epochMean'},
    'PS1BEST': {'RA': 'ra', 'RA_error': 'raMeanErr',
               'DEC': 'dec', 'DEC_error': 'decMeanErr',
               'pmra': '', 'pmra_error': '',
               'pmdec': '', 'pmdec_error': '',
               'mag': 'mag', 'objID': 'objid', 'epoch': 'epochMean'},
    'GSC242': {'RA': 'ra', 'RA_error': 'raErr',
               'DEC': 'dec', 'DEC_error': 'decErr',
               'pmra': 'rapm', 'pmra_error': 'rapmErr',
               'pmdec': 'decpm', 'pmdec_error': 'decpmErr',
               'mag': 'mag', 'objID': 'objID', 'epoch': 'epoch'},
    'GSC241': {'RA': 'ra', 'RA_error': 'raErr',
               'DEC': 'dec', 'DEC_error': 'decErr',
               'pmra': 'rapm', 'pmra_error': 'rapmErr',
               'pmdec': 'decpm', 'pmdec_error': 'decpmErr',
               'mag': 'mag', 'objID': 'objID', 'epoch': 'epoch'},
    }


# CRBIT definitions
CRBIT = 4096
"""

Primary function for creating an astrometric reference catalog.

"""


def create_astrometric_catalog(inputs, catalog="GAIAedr3", output="ref_cat.ecsv",
                               gaia_only=False, table_format="ascii.ecsv",
                               existing_wcs=None, num_sources=None,
                               use_footprint=False, full_catalog=False,
                               user_epoch='match', log_level=logutil.logging.NOTSET):
    """Create an astrometric catalog that covers the inputs' field-of-view.

    This function will return a table containing sources derived from the
    selected astrometric catalog specified by the user for the field-of-view
    covered by the list of input images.

    .. note::  The output table generated by this function will have a
        pre-defined set of columns with common names for catalogs included
        in `SUPPORTED_CATALOGS`.  Columns which don't exist in the catalog
        will be added with null (None) values.  The common column names are:

          * **RA, DEC** : RA, DEC of source
          * **RA_error, DEC_error** : error in RA, DEC position
          * **objID** : object ID from catalog
          * **pmra, pmdec** : Proper motion in RA and DEC
          * **pmra_error, pmdec_error** : errors in the proper motion
          * **mag** : reported magnitude for the source
          * **epoch** : epoch of RA, DEC after applying proper motions from catalog (if any)

    All tables where columns have been renamed have '.meta["converted"]' set to 'True'.

    Parameters
    ----------
    input : str, list
        Filenames of images to be aligned to astrometric catalog

    catalog : str, optional
        Name of catalog to extract astrometric positions for sources in the
        input images' field-of-view. Default: GAIAeDR3. Only those catalogs
        listed in `SUPPORTED_CATALOGS` will be returned with standardized
        column names.  Catalog names are case-insensitive.

    output : str, optional
        Filename to give to the astrometric catalog read in from the master
        catalog web service.  If None, no file will be written out.

    gaia_only : bool, optional
        Specify whether or not to only use sources from GAIA in output catalog

    existing_wcs : `~stwcs.wcsutil.HSTWCS`
        existing WCS object specified by the user

    num_sources : int
        Maximum number of brightest/faintest sources to return in catalog.
        If `num_sources` is negative, return that number of the faintest
        sources.  By default, all sources are returned.

    use_footprint : bool, optional
        Use the image footprint as the source identification bounds insted of using a circle centered on a
        given RA and Dec? By default, use_footprint = False.

    full_catalog : bool, optional
        Return the full set of columns provided by the web service.

    user_epoch : str or float, optional
        Epoch of returned astrometric catalog.  If 'match', match the
        epoch of the specified observation.  If None, return the values
        at the epoch of the catalog.  If a decimal year is specified,
        returned values will be for that epoch.

    Notes
    -----
    This function will point to astrometric catalog web service defined
    through the use of the ASTROMETRIC_CATALOG_URL environment variable.

    Returns
    -------
    ref_table : `~astropy.table.Table`
        Astropy Table object of the catalog

    """
    # Initialize logging for this user-callable function
    log.setLevel(log_level)

    inputs, _ = parseinput.parseinput(inputs)

    # start by creating a composite field-of-view for all inputs
    # This default output WCS will have the same plate-scale and orientation
    # as the first chip in the list, which for WFPC2 data means the PC.
    # Fortunately, for alignment, this doesn't matter since no resampling of
    # data will be performed
    if existing_wcs is not None:
        outwcs = existing_wcs
    else:
        outwcs = build_reference_wcs(inputs)  # TODO: replace with stwcs.updatewcs.astrometry_utils version

    if use_footprint:
        footprint = outwcs.calc_footprint()
        with open('footprint.reg', 'w') as reg_out:
            for item in footprint:
                reg_out.write("{}, {}\n".format(item[0], item[1]))
    else:
        radius = compute_radius(outwcs)
        ra, dec = outwcs.wcs.crval

    if user_epoch == 'match':
        # Get the observation date
        epoch = Time(fits.getval(inputs[0], 'date-obs')).decimalyear
    else:
        epoch = user_epoch

    # perform query for this field-of-view
    if use_footprint:
        ref_table = get_catalog_from_footprint(footprint, epoch=epoch, catalog=catalog)
    else:
        ref_table = get_catalog(ra, dec, sr=radius, epoch=epoch, catalog=catalog)

    if not ref_table:
        return ref_table

    # Add catalog name as meta data
    ref_table.meta['catalog'] = catalog
    ref_table.meta['gaia_only'] = gaia_only
    ref_table.meta['epoch'] = epoch
    ref_table.meta['converted'] = False

    # Convert common set of columns into standardized column names
    convert_astrometric_table(ref_table, catalog)

    if not full_catalog:
        ref_table = ref_table['RA', 'DEC', 'mag', 'objID']

    sources = len(ref_table)

    # sort table by magnitude, fainter to brightest
    ref_table.sort('mag', reverse=True)

    if num_sources is not None:
        indx = -1 * num_sources
        ref_table = ref_table[:indx] if num_sources < 0 else ref_table[indx:]
        sources_type = "faintest" if num_sources < 0 else "brightest"
        sources = abs(num_sources)
    else:
        sources_type = ""

    # Write out table to a file, if specified
    if output:
        ref_table.write(output, format=table_format, overwrite=True)
        log.info("Created catalog '{}' with {} {} sources".format(
                  output, sources, sources_type))

    return ref_table


def convert_astrometric_table(table, catalog_name):
    """Convert a table with varying column names into a more standardized table"""
    cat_name = catalog_name.upper()
    if cat_name in SUPPORTED_CATALOGS:
        cat_dict = SUPPORTED_CATALOGS[cat_name]
        table.meta['converted'] = True
        log.debug('Converting supported catalog: {}'.format(catalog_name))
    else:
        table.meta['converted'] = False
        log.debug('NOT converting catalog: {}'.format(catalog_name))
        return table

    # Now for each column specified in this dict,
    # translate the old column name into the new
    # standardized name
    for new, old in cat_dict.items():
        if new != 'epoch':
            rename_colname(table, old, new)
        elif old != '':
            # Add decimalyear epoch column as new column for existing 'epoch'-like column
            # determine format of data in old time/date column
            val = table[old][0]
            if val > 2400000:
                dfmt = 'jd'
            elif val > 10000:
                dfmt = 'mjd'
            else:
                dfmt = None

            if dfmt:
                # Insure 'epoch' is decimal year and add it as a new column
                newtimes = Time(table[old], format=dfmt).decimalyear
                timecol = Column(data=newtimes, name=new)
                table.add_column(timecol)
            else:
                # Otherwise, no format conversion needed, so simply rename it
                # since we already know the old column exists
                table.rename_column(old, new)
        else:
            # Insure at least an empty column is provided for 'epoch'
            rename_colname(table, old, new)

    return table


def rename_colname(table, colname, newcol):
    """Convert column name in table to user-specified name"""
    # If table is missing column, add a column with values of None
    if colname != '':
        table.rename_column(colname, newcol)
    else:
        empty_column = Column(data=[None]*len(table), name=newcol, dtype=np.float64)
        table.add_column(empty_column)


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
    if len(wcslist) == 1 and not wcslist[0].has_distortion and wcslist[0].idcscale is None:
        outwcs = wcslist[0]
    else:
        outwcs = utils.output_wcs(wcslist)

    return outwcs


def get_catalog(ra, dec, sr=0.1, epoch=None, catalog='GSC241'):
    """ Extract reference catalog from VO web service.

    Queries the catalog available at the ``SERVICELOCATION`` specified
    for this module to get any available astrometric source catalog entries
    around the specified position in the sky based on a cone-search.

    Parameters
    ----------
    ra : float
        Right Ascension (RA) of center of field-of-view (in decimal degrees)

    dec : float
        Declination (Dec) of center of field-of-view (in decimal degrees)

    sr : float, optional
        Search radius (in decimal degrees) from field-of-view center to use
        for sources from catalog.  Default: 0.1 degrees

    epoch : float, optional
        Catalog positions returned for this field-of-view will have their
        proper motions applied to represent their positions at this date, if
        a value is specified at all, for catalogs with proper motions.

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
    epoch_str = '&EPOCH={:.3f}'

    base_spec = spec_str.format(ra, dec, sr, fmt, catalog)
    spec = base_spec + epoch_str.format(epoch) if epoch else base_spec

    serviceUrl = '{}/{}?{}'.format(SERVICELOCATION, serviceType, spec)
    log.debug("Getting catalog using: \n    {}".format(serviceUrl))
    rawcat = requests.get(serviceUrl, headers=headers)
    r_contents = rawcat.content.decode()  # convert from bytes to a String
    rstr = r_contents.split('\r\n')
    # remove initial line describing the number of sources returned
    # CRITICAL to proper interpretation of CSV data
    if rstr[0].startswith('Error'):
        # Try again without EPOCH
        serviceUrl = '{}/{}?{}'.format(SERVICELOCATION, serviceType, base_spec)
        log.debug("Getting catalog using: \n    {}".format(serviceUrl))
        rawcat = requests.get(serviceUrl, headers=headers)
        r_contents = rawcat.content.decode()  # convert from bytes to a String
        rstr = r_contents.split('\r\n')

    # If we still have an error returned by the web-service, report the exact error
    if rstr[0].startswith('Error'):
        log.warning("Astrometric catalog generation FAILED with: \n{}".format(rstr))

    del rstr[0]
    r_csv = Table.read(rstr, format='ascii.csv')
    return r_csv


def get_catalog_from_footprint(footprint, epoch=None, catalog='GSC241'):
    """ Extract catalog from VO web service based on the specified footprint

    Queries the catalog available at the ``SERVICELOCATION`` specified
    for this module to get any available astrometric source catalog entries
    around the specified position in the sky based on the input footprint.

    Parameters
    ----------
    footprint : numpy.ndarray
        Array of RA, Dec points that describe the footprint polygon

    epoch : float, optional
        Catalog positions returned for this field-of-view will have their
        proper motions applied to represent their positions at this date, if
        a value is specified at all, for catalogs with proper motions.

    catalog : str, optional
        Name of catalog to query, as defined by web-service.  Default: 'GSC241'

    Returns
    -------
    csv : CSV object
        CSV object of returned sources with all columns as provided by catalog

    """
    serviceType = 'vo/CatalogSearch.aspx'
    spec_str = 'STCS=polygon{}&FORMAT={}&CAT={}&MINDET=5'
    headers = {'Content-Type': 'text/csv'}
    fmt = 'CSV'
    epoch_str = '&EPOCH={:.3f}'

    footprint_string = ""
    for item in footprint:
        footprint_string += "%20{}%20{}".format(item[0],item[1])

    base_spec = spec_str.format(footprint_string, fmt, catalog)
    spec = base_spec + epoch_str.format(epoch) if epoch else base_spec

    serviceUrl = '{}/{}?{}'.format(SERVICELOCATION, serviceType, spec)
    log.debug("Getting catalog using: \n    {}".format(serviceUrl))

    rawcat = requests.get(serviceUrl, headers=headers)
    r_contents = rawcat.content.decode()  # convert from bytes to a String
    rstr = r_contents.split('\r\n')

    # Check for any error conditions, and try without the epoch in case that was the problem
    if rstr[0].startswith('Error'):
        # Try again without EPOCH
        serviceUrl = '{}/{}?{}'.format(SERVICELOCATION, serviceType, base_spec)
        log.debug("Getting catalog using: \n    {}".format(serviceUrl))
        rawcat = requests.get(serviceUrl, headers=headers)
        r_contents = rawcat.content.decode()  # convert from bytes to a String
        rstr = r_contents.split('\r\n')

    # If we still have an error returned by the web-service, report the exact error
    if rstr[0].startswith('Error'):
        log.warning("Astrometric catalog generation FAILED with: \n{}".format(rstr))

    # remove initial line describing the number of sources returned
    # CRITICAL to proper interpretation of CSV data
    del rstr[0]

    r_csv = Table.read(rstr, format='ascii.csv')
    return r_csv


def compute_radius(wcs):
    """Compute the radius from the center to the furthest edge of the WCS.

    Parameters
    -----------
    wcs : `~astropy.wcsutil.HSTWCS` object
        HSTWCS object for an exposure including a valid footprint

    """

    ra, dec = wcs.wcs.crval
    img_center = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)
    wcs_foot = wcs.calc_footprint()
    img_corners = SkyCoord(ra=wcs_foot[:, 0] * u.degree,
                           dec=wcs_foot[:, 1] * u.degree)
    radius = img_center.separation(img_corners).max().value

    return radius

@deprecated(since="3.2.1",
            message="Replaced by stwcs.updatewcs.astrometry_utils.find_gsc_offset",
            name='drizzlepac.haputils.astrometric_utils.find_gsc_offset',
            alternative='stwcs.updatewcs.astrometry_utils.find_gsc_offset')
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


def compute_2d_background(imgarr, box_size, win_size,
                          bkg_estimator=SExtractorBackground,
                          rms_estimator=StdBackgroundRMS):
    """Compute a 2D background for the input array.

    This function uses `~photutils.background.Background2D` to determine
    an adaptive background that takes into account variations in flux
    across the image.

    Parameters
    ==========
    imgarr : ndarray
        NDarray of science data for which the background needs to be computed

    box_size : integer
        The box_size along each axis for Background2D to use.

    win_size : integer
        The window size of the 2D median filter to apply to the low-resolution map as the
        `filter_size` parameter in Background2D.

    bkg_estimator : function
        The name of the function to use as the estimator of the background.

    rms_estimator : function
        The name of the function to use for estimating the RMS in the background.

    Returns
    =======
    bkg_background : ndarray
        NDarray the same shape as the input image array which contains the determined
        background across the array.  If Background2D fails for any reason, a simpler
        sigma-clipped single-valued array will be computed instead.

    bkg_median : float
        The median value (or single sigma-clipped value) of the computed background.

    bkg_rms : ndarray
        NDarray the same shape as the input image array which contains the RMS of the
        background across the array.  If Background2D fails for any reason, a simpler
        sigma-clipped single-valued array will be computed instead.

    bkg_rms_median : float
        The median value (or single sigma-clipped value) of the RMS of the computed
        background.
    """

    # SExtractorBackground and StdBackgroundRMS are the defaults
    bkg = None

    exclude_percentiles = [10, 25, 50, 75]
    for percentile in exclude_percentiles:
        log.info("")
        log.info("Percentile in use: {}".format(percentile))
        try:
            bkg = Background2D(imgarr, (box_size, box_size), filter_size=(win_size, win_size),
                               bkg_estimator=bkg_estimator(),
                               bkgrms_estimator=rms_estimator(),
                               exclude_percentile=percentile, edge_method="pad")

        except Exception:
            bkg = None
            continue

        if bkg is not None:
            bkg_background = bkg.background
            bkg_median = bkg.background_median
            bkg_rms = bkg.background_rms
            bkg_rms_median = bkg.background_rms_median
            break

    # If Background2D does not work at all, define default scalar values for
    # the background to be used in source identification
    if bkg is None:
        log.info("Background2D failure detected. Using alternative background calculation instead....")
        mask = make_source_mask(imgarr, nsigma=2, npixels=5, dilate_size=11)
        sigcl_mean, sigcl_median, sigcl_std = sigma_clipped_stats(imgarr, sigma=3.0, mask=mask, maxiters=9)
        bkg_median = max(0.0, sigcl_median)
        bkg_rms_median = sigcl_std
        # create background frame shaped like imgarr populated with sigma-clipped median value
        bkg_background = np.full_like(imgarr, bkg_median)
        # create background frame shaped like imgarr populated with sigma-clipped standard deviation value
        bkg_rms = np.full_like(imgarr, sigcl_std)

    return bkg_background, bkg_median, bkg_rms, bkg_rms_median

def build_auto_kernel(imgarr, whtarr, fwhm=3.0, threshold=None, source_box=7,
                      good_fwhm=[1.0, 4.0], num_fwhm=30,
                      isolation_size=11, saturation_limit=70000., log_level=logutil.logging.NOTSET):
    """Build kernel for use in source detection based on image PSF
    This algorithm looks for an isolated point-source that is non-saturated to use as a template
    for the source detection kernel.  Failing to find any suitable sources, it will return a
    Gaussian2DKernel based on the provided FWHM as a default.

    Parameters
    ----------
    imgarr : ndarray
        Image array (ndarray object) with sources to be identified

    whtarr : ndarray
        Image array (ndarray object) for `imgarr` that provides an estimate of the
        weight or errors for the processed image.

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
    # Initialize logging for this user-callable function
    log.setLevel(log_level)

    # Try to use PSF derived from image as detection kernel
    # Kernel must be derived from well-isolated sources not near the edge of the image
    kern_img = imgarr.copy()
    edge = source_box * 2
    kern_img[:edge, :] = 0.0
    kern_img[-edge:, :] = 0.0
    kern_img[:, :edge] = 0.0
    kern_img[:, -edge:] = 0.0
    kernel_psf = False

    peaks = find_peaks(kern_img, threshold=threshold * 5, box_size=isolation_size)
    if peaks is None or (peaks is not None and len(peaks) == 0):
        tmean = threshold.mean() if isinstance(threshold, np.ndarray) else threshold
        if tmean > kern_img.mean():
            kern_stats = sigma_clipped_stats(kern_img)
            threshold = kern_stats[2]
        peaks = find_peaks(kern_img, threshold=threshold, box_size=isolation_size)

    if peaks is not None:
        # Sort based on peak_value to identify brightest sources for use as a kernel
        peaks.sort('peak_value', reverse=True)

        if saturation_limit:
            sat_peaks = np.where(peaks['peak_value'] > saturation_limit)[0]
            sat_index = sat_peaks[-1] + 1 if len(sat_peaks) > 0 else 0
            peaks['peak_value'][:sat_index] = 0.

        wht_box = 2  # Weight image cutout box size is 2 x wht_box + 1 pixels on a side
        fwhm_attempts = 0
        # Identify position of brightest, non-saturated peak (in numpy index order)
        for peak_ctr in range(len(peaks)):
            kernel_pos = [peaks['y_peak'][peak_ctr], peaks['x_peak'][peak_ctr]]

            kernel = imgarr[kernel_pos[0] - source_box:kernel_pos[0] + source_box + 1,
                            kernel_pos[1] - source_box:kernel_pos[1] + source_box + 1].copy()

            kernel_wht = whtarr[kernel_pos[0] - wht_box:kernel_pos[0] + wht_box + 1,
                         kernel_pos[1] - wht_box:kernel_pos[1] + wht_box + 1].copy()

            minsize = min(kernel.shape)
            # search square cut-out (of size 2 x wht_box + 1 pixels on a side) of weight image centered on peak coords for
            # zero-value pixels. Reject peak if any are found.
            if len(np.where(kernel_wht == 0.)[0]) == 0 and minsize > 11:
                log.debug("Kernel source PSF located at [{},{}]".format(kernel_pos[1], kernel_pos[0]))
            else:
                kernel = None
            if kernel is not None:
                kernel = np.clip(kernel, 0, None)  # insure background subtracted kernel has no negative pixels
                if kernel.sum() > 0.0:
                    kernel /= kernel.sum()  # Normalize the new kernel to a total flux of 1.0
                    kernel_fwhm = find_fwhm(kernel, fwhm)
                    fwhm_attempts += 1
                    if kernel_fwhm is None:
                        kernel = None
                    else:
                        log.debug("Determined FWHM from sample PSF of {:.2f}".format(kernel_fwhm))
                        log.debug("  based on good range of FWHM:  {:.1f} to {:.1f}".format(good_fwhm[0], good_fwhm[1]))
                        if good_fwhm[1] > kernel_fwhm > good_fwhm[0]:  # This makes it hard to work with sub-sampled data (WFPC2?)
                            fwhm = kernel_fwhm
                            kernel_psf = True
                            break
                        else:
                            kernel = None
                if fwhm_attempts == num_fwhm:
                    break
    else:
        kernel = None

    if kernel is None:
        num_peaks = len(peaks) if peaks else 0
        log.warning("Did not find a suitable PSF out of {} possible sources...".format(num_peaks))
        log.warning("Using a Gaussian 2D Kernel for source detection.")
        # Generate a default kernel using a simple 2D Gaussian
        kernel_fwhm = fwhm
        sigma = fwhm * gaussian_fwhm_to_sigma
        k = Gaussian2DKernel(sigma, x_size=source_box, y_size=source_box)
        k.normalize()
        kernel = k.array

    return (kernel, kernel_psf), kernel_fwhm


def find_fwhm(psf, default_fwhm):
    """Determine FWHM for auto-kernel PSF

    This function iteratively fits a Gaussian model to the extracted PSF
    using `photutils.psf.IterativelySubtractedPSFPhotometry` to determine
    the FWHM of the PSF.

    Parameters
    -----------
    psf : ndarray
        Array (preferably a slice) containing the PSF to be measured.

    default_fwhm : float
        Starting guess for the FWHM

    Returns
    --------
    fwhm : float
        Value of the computed Gaussian FWHM for the PSF

    """
    daogroup = DAOGroup(crit_separation=8)
    mmm_bkg = MMMBackground()
    iraffind = DAOStarFinder(threshold=2.5 * mmm_bkg(psf), fwhm=default_fwhm)
    fitter = LevMarLSQFitter()
    sigma_psf = gaussian_fwhm_to_sigma * default_fwhm
    gaussian_prf = IntegratedGaussianPRF(sigma=sigma_psf)
    gaussian_prf.sigma.fixed = False
    itr_phot_obj = IterativelySubtractedPSFPhotometry(finder=iraffind,
                                                      group_maker=daogroup,
                                                      bkg_estimator=mmm_bkg,
                                                      psf_model=gaussian_prf,
                                                      fitter=fitter,
                                                      fitshape=(11, 11),
                                                      niters=2)
    phot_results = itr_phot_obj(psf)
    # Insure none of the fluxes determined by photutils is np.nan
    phot_results['flux_fit'] = np.nan_to_num(phot_results['flux_fit'].data, nan=0)

    if len(phot_results['flux_fit']) == 0:
        return None
    psf_row = np.where(phot_results['flux_fit'] == phot_results['flux_fit'].max())[0][0]
    sigma_fit = phot_results['sigma_fit'][psf_row]
    fwhm = gaussian_sigma_to_fwhm * sigma_fit
    log.debug("Found FWHM: {}".format(fwhm))

    return fwhm


def extract_point_sources(img, dqmask=None, fwhm=3.0, kernel=None,
                          nbright=1000,
                          threshold=200.0, sigma=3.0, source_box=7, log_level=logutil.logging.NOTSET):
    """Use photutils to replicate the IRAF point-source catalogs"""
    # Initialize logging for this user-callable function
    log.setLevel(log_level)

    # Detect threshold using a relatively fast method and
    # subtract off that background.
    nsigma = 5.0
    bkg_thresh, bkg = sigma_clipped_bkg(img, sigma=sigma, nsigma=nsigma)

    sigma = np.sqrt(2.0 * np.abs(bkg[1]))
    x, y, flux, src_id, sharp, round1, round2 = ndfind(img,
                                                       sigma * threshold,
                                                       fwhm, bkg[1],
                                                       nbright=nbright,
                                                       use_sharp_round=True)
    srcs = Table([x, y, flux, src_id], names=['xcentroid', 'ycentroid', 'flux', 'id'])

    """
    # Now, use IRAFStarFinder to identify sources across chip
    starfind = IRAFStarFinder(threshold=bkg[2]*nsigma, fwhm=fwhm)
    srcs = starfind.find_stars(img, mask=dqmask)
    if srcs is not None and high_sn is not None and len(srcs) > high_sn:
        # sort by flux, return high_sn srcs only...
        indx = np.argsort(srcs['flux'])[:high_sn]
        srcs = srcs[indx]
    """
    num_srcs = len(srcs) if srcs is not None else 0
    log.info("Found {} sources".format(num_srcs))

    return srcs


def extract_sources(img, dqmask=None, fwhm=3.0, kernel=None, photmode=None,
                    segment_threshold=None, dao_threshold=None,
                    dao_nsigma=3.0, source_box=7,
                    classify=True, centering_mode="starfind", nlargest=None,
                    outroot=None, plot=False, vmax=None, deblend=False, log_level=logutil.logging.NOTSET):
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
        will rely on `photutils.segmentation.SourceCatalog` to generate the
        properties for the source catalog.  Centering using `starfind` will use
        `photutils.detection.IRAFStarFinder` to characterize each source in the catalog.
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
    # Initialize logging for this user-callable function
    log.setLevel(log_level)

    # apply any provided dqmask for segmentation only
    imgarr = img.copy()
    if dqmask is not None:
        imgarr[dqmask] = 0

    if segment_threshold is None:
        dao_threshold, bkg = sigma_clipped_bkg(imgarr, sigma=4.0, nsigma=dao_nsigma)
        segment_threshold = np.ones(imgarr.shape, imgarr.dtype) * dao_threshold

    segm = detect_sources(imgarr, segment_threshold, npixels=source_box,
                          filter_kernel=kernel, connectivity=4)

    # photutils >= 0.7: segm=None; photutils < 0.7: segm.nlabels=0
    if segm is None or segm.nlabels == 0:
        log.info("No detected sources!")
        return None, None, None

    log.debug("Creating segmentation map for {} ".format(outroot))
    if kernel is not None:
        kernel_area = ((kernel.shape[0] // 2) ** 2) * np.pi
        log.debug("   based on kernel shape of {}".format(kernel.shape))
    else:
        kernel_area = ((source_box // 2) ** 2) * np.pi
        log.debug("   based on a default kernel.")

    num_brightest = 10 if len(segm.areas) > 10 else len(segm.areas)
    mean_area = np.mean(segm.areas)
    max_area = np.sort(segm.areas)[-1 * num_brightest:].mean()
    # This section looks for crowded fields where segments run into each other
    # By reducing the size of the kernel used for segment detection, this can be minimized
    # in crowded fields.  Also, mean area is used to try to avoid this logic for fields with
    # several large extended sources in an otherwise empty field.
    if max_area > MAX_AREA_LIMIT and mean_area > (kernel_area / 2):  # largest > 25-pix radius source
        # reset kernel to only use the central 1/4 area and redefine the segment map
        kcenter = (kernel.shape[0] - 1) // 2
        koffset = (kcenter - 1) // 2
        kernel = kernel[kcenter - koffset: kcenter + koffset + 1,
                        kcenter - koffset: kcenter + koffset + 1].copy()
        kernel /= kernel.sum()  # normalize to total sum == 1
        log.info("Looking for crowded sources using smaller kernel with shape: {}".format(kernel.shape))
        segm = detect_sources(imgarr, segment_threshold, npixels=source_box, filter_kernel=kernel)

    if deblend:
        segm = deblend_sources(imgarr, segm, npixels=5,
                               filter_kernel=kernel, nlevels=32,
                               contrast=0.01)

    # If classify is turned on, it should modify the segmentation map
    dqmap = None
    if classify:
        log.info('Removing suspected cosmic-rays from source catalog')
        cat = SourceCatalog(imgarr, segm)
        # Remove likely cosmic-rays based on central_moments classification
        bad_srcs = np.where(classify_sources(cat, fwhm) == 0)[0] + 1
        # Convert this bad_srcs into a segmap that can be used to update a DQ array
        dqmap = np.zeros_like(segm.data)
        for src in bad_srcs:
            dqmap[segm.data == src] += CRBIT

        segm.remove_labels(bad_srcs)

    if OLD_PHOTUTILS:
        flux_colname = 'source_sum'
        ferr_colname = 'source_sum_err'
    else:
        flux_colname = 'segment_flux'
        ferr_colname = 'segment_fluxerr'

    # convert segm to mask for daofind
    if centering_mode == 'starfind':
        src_table = None

        # Identify nbrightest/largest sources
        if nlargest is not None:
            nlargest = min(nlargest, len(segm.labels))

            # Look for brightest sources by flux...
            src_fluxes = np.array([imgarr[src].max() for src in segm.slices])
            src_labels = np.array([label for label in segm.labels])
            src_brightest = np.flip(np.argsort(src_fluxes))
            large_labels = src_labels[src_brightest]
            log.debug("Brightest sources in segments: \n{}".format(large_labels))
        else:
            src_brightest = np.arange(len(segm.labels))

        log.info("Looking for sources in {} segments".format(len(segm.labels)))

        for indx in src_brightest:
            segment = segm.segments[indx]

        # for segment in segm.segments:
            # check needed for photutils <= 0.6; it can be removed when
            # the drizzlepac depends on photutils >= 0.7
            if segment is None:
                continue

            # Get slice definition for the segment with this label
            seg_slice = segment.slices
            seg_yoffset = seg_slice[0].start
            seg_xoffset = seg_slice[1].start

            dao_threshold = segment_threshold[seg_slice].mean()
            daofind = DAOStarFinder(fwhm=fwhm, threshold=dao_threshold)
            log.debug("Setting up DAOStarFinder with: \n    fwhm={}  threshold={}".format(fwhm, dao_threshold))

            # Define raw data from this slice
            detection_img = img[seg_slice]
            # zero out any pixels which do not have this segments label
            detection_img[segm.data[seg_slice] == 0] = 0

            # Detect sources in this specific segment
            seg_table = daofind.find_stars(detection_img)

            # Pick out brightest source only
            if src_table is None and seg_table:
                # Initialize final master source list catalog
                log.debug("Defining initial src_table based on: {}".format(seg_table.colnames))
                src_table = Table(names=seg_table.colnames,
                                  dtype=[dt[1] for dt in seg_table.dtype.descr])

            if seg_table:
                # This logic will eliminate saturated sources, where the max pixel value is not
                # the center of the PSF (saturated and streaked along the Y axis)
                max_row = np.where(seg_table['peak'] == seg_table['peak'].max())[0][0]

                # Add logic to remove sources which have more than 3 pixels
                # within 10% of the max value in the source segment, a situation
                # which would indicate the presence of a saturated source
                if (detection_img > detection_img.max() * 0.9).sum() > 3:

                    # Revert to segmentation photometry for sat. source posns
                    if OLD_PHOTUTILS:
                        segment_properties = SourceCatalog(detection_img, segment.data)
                    else:
                        segimg = SegmentationImage(segment.data)
                        segment_properties = SourceCatalog(detection_img, segimg)

                    sat_table = segment_properties.to_table()
                    seg_table['flux'][max_row] = sat_table[flux_colname][0]
                    seg_table['peak'][max_row] = sat_table['max_value'][0]
                    if OLD_PHOTUTILS:
                        xcentroid = sat_table['xcentroid'][0].value
                        ycentroid = sat_table['ycentroid'][0].value
                        sky = sat_table['background_mean'][0].value
                    else:
                        xcentroid = sat_table['xcentroid'][0]
                        ycentroid = sat_table['ycentroid'][0]
                        sky = sat_table['local_background'][0]
                    seg_table['xcentroid'][max_row] = xcentroid
                    seg_table['ycentroid'][max_row] = ycentroid
                    seg_table['npix'][max_row] = sat_table['area'][0].value
                    seg_table['sky'][max_row] = sky if sky is not None and not np.isnan(sky) else 0.0
                    seg_table['mag'][max_row] = -2.5 * np.log10(sat_table[flux_colname][0])

                # Add row for detected source to master catalog
                # apply offset to slice to convert positions into full-frame coordinates
                seg_table['xcentroid'] += seg_xoffset
                seg_table['ycentroid'] += seg_yoffset
                src_table.add_row(seg_table[max_row])

            # If we have accumulated the desired number of sources, stop looking for more...
            if nlargest is not None and src_table is not None and len(src_table) == nlargest:
                break
    else:
        log.debug("Determining source properties as src_table...")
        cat = SourceCatalog(img, segm)
        src_table = cat.to_table()
        # Make column names consistent with IRAFStarFinder column names
        src_table.rename_column(flux_colname, 'flux')
        src_table.rename_column(ferr_colname, 'flux_err')
        src_table.rename_column('max_value', 'peak')

    if src_table is not None:
        log.info("Total Number of detected sources: {}".format(len(src_table)))
    else:
        log.info("No detected sources!")
        return None, None, None

    # Include magnitudes for each source for use in verification of alignment through
    # comparison with GAIA magnitudes
    tbl = compute_photometry(src_table, photmode)

    # Insure all IDs are sequential and unique (at least in this catalog)
    tbl['cat_id'] = np.arange(1, len(tbl) + 1)

    if outroot:
        tbl['xcentroid'].info.format = '.10f'  # optional format
        tbl['ycentroid'].info.format = '.10f'
        tbl['flux'].info.format = '.10f'
        if not outroot.endswith('.cat'):
            outroot += '.cat'
        tbl.write(outroot, format='ascii.commented_header', overwrite=True)
        log.info("Wrote source catalog: {}".format(outroot))

    if plot and plt is not None:
        norm = len(segm.labels)
        if vmax is None:
            norm = ImageNormalize(stretch=SqrtStretch())
        fig, ax = plt.subplots(2, 2, figsize=(8, 8))
        ax[0][0].imshow(imgarr, origin='lower', cmap='Greys_r', norm=norm, vmax=vmax)
        ax[0][1].imshow(segm, origin='lower', cmap=segm.cmap(random_state=12345))
        ax[0][1].set_title('Segmentation Map')
        if not isinstance(segment_threshold, float):
            ax[1][1].imshow(segment_threshold, origin='lower')
    return tbl, segm, dqmap


def crclean_image(imgarr, segment_threshold, kernel, fwhm,
                  source_box=7, deblend=False, background=0.0):
    """ Apply moments-based single-image CR clean algorithm to image

    Returns
    -------
    imgarr : ndarray
        Input array with pixels flagged as CRs set to 0.0
    """
    # Perform basic image segmentation to look for sources in image
    segm = detect_sources(imgarr, segment_threshold, npixels=source_box,
                          filter_kernel=kernel, connectivity=4)

    log.debug("Detected {} sources of all types ".format(segm.nlabels))
    # photutils >= 0.7: segm=None; photutils < 0.7: segm.nlabels=0
    if segm is None or segm.nlabels == 0:
        log.info("No sources to identify as cosmic-rays!")
        return imgarr

    if deblend:
        segm = deblend_sources(imgarr, segm, npixels=5,
                               filter_kernel=kernel, nlevels=32,
                               contrast=0.01)

    # Measure segment sources to compute centralized moments for each source
    cat = SourceCatalog(imgarr, segm)
    # Remove likely cosmic-rays based on central_moments classification
    bad_srcs = classify_sources(cat, fwhm)
    # Get the label IDs for sources flagged as CRs, IDs are 1-based not 0-based
    bad_ids = np.where(bad_srcs == 0)[0] + 1
    log.debug("Recognized {} sources as cosmic-rays".format(len(bad_ids)))

    # Instead of ignoring them by removing them from the catalog
    # zero all these sources out from the image itself.
    # Create a mask from the segmentation image with only the flagged CRs
    cr_segm = copy.deepcopy(segm)
    cr_segm.keep_labels(bad_ids)
    cr_mask = np.where(cr_segm.data > 0)
    # Use the mask to replace all pixels associated with the CRs
    # with values specified as the background
    if not isinstance(background, np.ndarray):
        background = np.ones_like(imgarr, dtype=np.float32) * background

    imgarr[cr_mask] = background[cr_mask]
    del cr_segm

    log.debug("Finshed zeroing out cosmic-rays")

    return imgarr

def classify_sources(catalog, fwhm, sources=None):
    """ Convert moments_central attribute for source catalog into star/cr flag.

    This algorithm interprets the central_moments from the SourceCatalog
    generated for the sources as more-likely a star or a cosmic-ray.  It is not
    intended or expected to be precise, merely a means of making a first cut at
    removing likely cosmic-rays or other artifacts.

    Parameters
    ----------
    catalog : `~photutils.segmentation.SourceCatalog`
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
    # All these return ndarray objects
    moments = catalog.moments_central  # (num_sources, 4, 4)
    semiminor_axis = catalog.semiminor_sigma.to_value()  # (num_sources,)
    elon = catalog.elongation.to_value()  # (num_sources,)
    # Determine how many sources we have.
    if sources is None:
        sources = (0, len(moments))
    num_sources = sources[1] - sources[0]
    srctype = np.zeros((num_sources,), np.int32)

    # Find which moment is the peak moment for all sources
    mx = np.array([np.where(moments[src] == moments[src].max())[0][0] for src in range(num_sources)], np.int32)
    my = np.array([np.where(moments[src] == moments[src].max())[1][0] for src in range(num_sources)], np.int32)
    # Now avoid processing sources which had NaN as either/both of their X,Y positions (rare, but does happen)
    src_x = catalog.xcentroid  # (num_sources, )
    src_y = catalog.ycentroid  # (num_sources, )
    valid_xy = ~np.bitwise_or(np.isnan(src_x), np.isnan(src_y))

    # Define descriptors of CRs
    # First: look for sources where the moment in X and Y is the first moment
    valid_src = np.bitwise_and(mx > 1, my > 1)[valid_xy]
    # Second: look for sources where the width of the source is less than the PSF FWHM
    valid_width = (semiminor_axis < (0.75 * fwhm))[valid_xy]
    # Third: identify sources which are elongated at least 2:1
    valid_elon = (elon > 2)[valid_xy]
    # Now combine descriptors into a single value
    valid_streak = ~np.bitwise_and(valid_width, valid_elon)
    src_cr = np.bitwise_and(valid_src, valid_streak)
    # Flag identified sources as CRs with a value of 1
    srctype[src_cr] = 1

    return srctype


def generate_source_catalog(image, dqname="DQ", output=False, fwhm=3.0,
                            **detector_pars):
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
    def_fwhmpsf = detector_pars.get('fwhmpsf', 0.13) / 2.0
    if 'fwhmpsf' in detector_pars:
        del detector_pars['fwhmpsf']
    source_box = detector_pars.get('source_box', 7)
    isolation_size = detector_pars.get('isolation_size', 11)
    saturation_limit = detector_pars.get('saturation_limit', 70000.0)
    if 'threshold' in detector_pars:
        del detector_pars['threshold']
    box_size = detector_pars.get('bkg_box_size', 27)
    win_size = detector_pars.get('bkg_filter_size', 3)
    nsigma = detector_pars.get('nsigma', 5)
    sat_flags = detector_pars.get('detector_pars', 256)
    if 'sat_flags' in detector_pars:
        del detector_pars['sat_flags']

    # Build source catalog for entire image
    source_cats = {}
    numSci = countExtn(image, extname='SCI')
    numWht = countExtn(image, extname='WHT')
    outroot = None
    img_inst = image[0].header['instrume']
    img_det = image[0].header['detector']

    for chip in range(numSci):
        chip += 1
        # find sources in image
        if output:
            rootname = image[0].header['rootname']
            outroot = '{}_sci{}_src'.format(rootname, chip)
        try:
            sci_ext = 0 if "{}/{}".format(img_inst, img_det) == "WFC3/IR" else ('sci', chip)
            photmode = {'photflam': image[sci_ext].header['photflam'],
                        'photplam': image[sci_ext].header['photplam']}
        except KeyError:
            photmode = None

        imgarr = image['sci', chip].data
        wcs = wcsutil.HSTWCS(image, ext=('sci', chip))
        def_fwhm = def_fwhmpsf / wcs.pscale

        # apply any DQ array, if available
        dqmask = None
        if image.index_of(dqname):
            dqarr = image[dqname, chip].data.astype(np.int16)

            # "grow out" regions in DQ mask flagged as saturated by several
            # pixels in every direction to prevent the
            # source match algorithm from trying to match multiple sources
            # from one image to a single source in the
            # other or vice-versa.
            # Create temp DQ mask containing all pixels flagged with any value EXCEPT 256
            non_sat_mask = bitfield_to_boolean_mask(dqarr, ignore_flags=sat_flags)

            # Create temp DQ mask containing saturated pixels ONLY
            sat_mask = bitfield_to_boolean_mask(dqarr, ignore_flags=~sat_flags)

            # Grow out saturated pixels by a few pixels in every direction
            grown_sat_mask = ndimage.binary_dilation(sat_mask, iterations=2)

            # combine the two temporary DQ masks into a single composite DQ mask.
            dqmask = np.bitwise_or(non_sat_mask, grown_sat_mask)
            # astropy's code returned the opposite bitmask from what was originally
            # defined by stsci.tools own bitmask code.
            if LooseVersion(stsci.tools.__version__) >= '4.0.0':
                dqmask = ~dqmask

        if numWht > 0:
            whtarr = image['wht', chip].data
        else:
            errarr = image['err', chip].data
            whtarr = errarr.max() / errarr
            whtarr[dqmask] = 0

        bkg_ra, bkg_median, bkg_rms_ra, bkg_rms_median = compute_2d_background(imgarr, box_size, win_size)
        threshold = nsigma * bkg_rms_ra
        dao_threshold = nsigma * bkg_rms_median

        (kernel, kernel_psf), kernel_fwhm = build_auto_kernel(imgarr - bkg_ra, whtarr,
                                                              threshold=threshold,
                                                              fwhm=def_fwhm,
                                                              source_box=source_box,
                                                              isolation_size=isolation_size,
                                                              saturation_limit=saturation_limit)
        log.debug("Built kernel with FWHM = {}".format(kernel_fwhm))

        seg_tab, segmap, crmap = extract_sources(imgarr - bkg_ra, dqmask=dqmask,
                                                 outroot=outroot, kernel=kernel,
                                                 photmode=photmode,
                                                 segment_threshold=threshold,
                                                 dao_threshold=dao_threshold,
                                                 fwhm=kernel_fwhm, **detector_pars)
        del crmap
        source_cats[chip] = seg_tab

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
        refwcs = astrometry_utils.build_reference_wcs([image])
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


def compute_photometry(catalog, photvals):
    """ Compute magnitudes for sources from catalog based on observations photmode.

    Magnitudes will be AB mag values.

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
        ABMAG photometry (in magnitudes).
    """
    if photvals is None:
        source_phot = np.array([-99.99] * len(catalog['flux']))
    else:
        ab_zpt = -2.5 * np.log10(photvals['photflam']) - 21.10 - 5 * np.log10(photvals['photplam']) + 18.692
        source_phot = ab_zpt - 2.5 * np.log10(catalog['flux'])

    # Label the new column
    phot_col = Column(data=source_phot, name='abmag')
    # Now add this new column to the catalog table
    catalog.add_column(phot_col)

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

    wcslin = astrometry_utils.build_reference_wcs([filename], sciname=sciname)

    if clean_wcs:
        wcsbase = wcslin.wcs
        customwcs = build_hstwcs(wcsbase.crval, wcsbase.crpix,
                                 wcslin.naxis1, wcslin.naxis2,
                                 wcslin.pscale, wcslin.orientat)
    else:
        customwcs = wcslin
    return customwcs

@deprecated(since="3.2.1",
            name='drizzlepac.haputils.astrometric_utils.read_hlet_wcs',
            message='Replaced by stwcs.updatewcs.astrometry_utils.read_hlet_wcs',
            alternative='stwcs.updatewcs.astrometry_utils.read_hlet_wcs')
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
    mask : ndarray
        Boolean array of same length as x,y arrays where sources that fall
        within the footprint are True.

    """
    # start with limits of WCS shape

    sky = wcsobj.pixel_to_world(x, y)
    wcsmask = wcsobj.footprint_contains(sky)
    xint = x[wcsmask].astype(np.int32)
    yint = y[wcsmask].astype(np.int32)

    fprint = ~np.isnan(img)
    xymask = np.zeros(img.shape, dtype=bool)
    xymask[yint, xint] = True

    skymask = np.bitwise_and(fprint, xymask)

    mask = [True if skymask[yx]==True else False for yx in zip(yint, xint)]
    mask = np.array(mask).astype(bool)

    # NOTE: There is probably a way to use list comprehension to do this,
    # but for now, this works as intended.
    onimg = wcsmask.copy()
    indx=0
    for i,w in enumerate(onimg):
        if w == True:
            if mask[indx] == False:
                onimg[i] = False
            indx += 1

    return onimg


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
    mask = within_footprint(image, refwcs, xref, yref)
    xref = xref[mask]
    yref = yref[mask]
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
    log.debug('best offset {} based on {} cross-matches'.format(hist2d_offset, nmatches))

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
        fname = image
    elif isinstance(image, fits.HDUList):
        hdulist = image
        fname = image.filename()
    else:
        log.info("Wrong type of input, {}, for build_wcscat...".format(type(image)))
        raise ValueError

    # interpret date for exposure
    epoch = get_epoch(hdulist)
    wcs_catalogs = []
    numsci = countExtn(hdulist)
    for chip in range(1, numsci + 1):
        w = wcsutil.HSTWCS(hdulist, ('SCI', chip))
        imcat = None
        if source_catalog:
            imcat = source_catalog[chip]

        # rename xcentroid/ycentroid columns, if necessary, to be consistent with tweakwcs
        if imcat is None:
            imcat = Table(names=['xcentroid', 'ycentroid', 'mag'])
        if isinstance(imcat, str):
            imcat = Table.read(imcat, format='ascii.fast_commented_header',
                                names=['x', 'y'])
            if 'mag' not in imcat.colnames:
                imcat['mag'] = [-999.9] * len(imcat['x'])

        if 'xcentroid' in imcat.colnames:
            imcat.rename_column('xcentroid', 'x')
            imcat.rename_column('ycentroid', 'y')

        wcscat = FITSWCS(
            w,
            meta={
                'chip': chip,
                'group_id': group_id,
                'filename': fname,
                'rootname': "_".join(fname.split("_")[:-1]),
                'catalog': imcat,
                'epoch': epoch,
                'name': fname
            }
        )

        wcs_catalogs.append(wcscat)

    if open_file:
        hdulist.close()

    return wcs_catalogs


def get_epoch(hdulist):
    """Interpret observation date 'date-obs' as epoch

    Parameters
    ----------
    hdulist : `fits.HDUList`
        Input HDUList object containing the 'date-obs' keyword in
        the Primary header in 'YYYY-MM-DD' format.

    Returns
    -------
    decimal_year : float
        Floating value for the year of the observation

    """
    d = hdulist[0].header['date-obs']
    atime = Time(d)
    return atime.decimalyear


def since_epoch(date): # returns seconds since epoch
    return time.mktime(date.timetuple())


# -------------------------------------------------------------------------------------------------------------
#
#  Utilities and supporting functions for verifying alignment
#
#
def check_mag_corr(imglist, threshold=0.5):
    """Check the correlation between input magnitudes and matched ref magnitudes."""
    mag_checks = []
    for image in imglist:
        input_mags = image.meta['fit_info']['input_mag']
        ref_mags = image.meta['fit_info']['ref_mag']
        if input_mags is not None and len(input_mags) > 1:
            mag_corr, mag_corr_std = pearsonr(input_mags, ref_mags)
            log.info("{} Magnitude correlation: {}".format(image.meta['name'], mag_corr))
            cross_match_check = True if abs(mag_corr) > threshold else False
        else:
            cross_match_check = False
        mag_checks.append(cross_match_check)

    return mag_checks


def rebin(arr, new_shape):
    """Rebin 2D array arr to shape new_shape by summing."""
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).sum(-1).sum(1)


def maxBit(int_val):
    """Return power of 2 for highest bit set for integer"""
    length = 0
    count = 0
    while (int_val):
        count += (int_val & 1)
        length += 1
        int_val >>= 1

    return length - 1


def compute_similarity(image, reference):
    """Compute a similarity index for an image compared to a reference image.

    Similarity index is based on a the general algorithm used in the AmphiIndex algorithm.

        - identify slice of image that is a factor of 256 in size
        - rebin image slice down to a (256,256) image
        - rebin same slice from reference down to a (256,256) image
        - sum the differences of the rebinned slices
        - divide absolute value of difference scaled by reference slice sum

    .. note:: This index will typically return values < 0.1 for similar images, and
              values > 1 for dis-similar images.

    Parameters
    ----------
    image : ndarray
        Image (as ndarray) to measure

    reference : ndarray
        Image which serves as the 'truth' or comparison image.

    Returns
    -------
    similarity_index : float
        Value of similarity index for `image`

    """

    # Insure NaNs are replaced with 0
    image = np.nan_to_num(image[:], nan=0)
    reference = np.nan_to_num(reference[:], nan=0)

    imgshape = (min(image.shape[0], reference.shape[0]),
                min(image.shape[1], reference.shape[1]))
    minsize = min(imgshape[0], imgshape[1])

    # determine largest slice that is a power of 2 in size
    window_bit = maxBit(minsize)
    window = 2**window_bit

    # Define how big the rebinned image should be for computing the sim index
    # Insure a minimum rebinned size of 64x64
    sim_bit = (window_bit - 2) if (window_bit - 2) > 6 else window_bit
    sim_size = 2**sim_bit

    # rebin image and reference
    img = rebin(image[:window, :window], (sim_size, sim_size))
    ref = rebin(reference[:window, :window], (sim_size, sim_size))

    # Compute index
    diffs = np.abs((img - ref).sum())
    sim_indx = diffs / img.sum()
    return sim_indx


def compute_prob(val, mean, sigma):
    """Return z-score for val relative to a distribution

       If abs(z_score) > 1, `val` is most likely not from the
       specified distribution.

    """
    p = st.norm.cdf(x=val, loc=mean, scale=sigma)
    z_score = st.norm.ppf(p)

    return z_score


def determine_focus_index(img, sigma=1.5):
    """Determine blurriness indicator for an image

       This returns a single value that serves as an indication of the
       sharpness of the image based on the max pixel value from the image
       after applying a Laplacian-of-Gaussian filter with sigma.

       Implementation based on discussion from:
       https://stackoverflow.com/questions/7765810/is-there-a-way-to-detect-if-an-image-is-blurry

       as supported by:
       Pertuz, S., et.al., 2013, Pattern Recognition, 46:1415–1432

       This index needs to be based on 'max' value in order to avoid
       field-dependent biases, since the 'max' value will correspond
       to point-source-like sources regardless of the field
       (nebula, galaxy, ...).  Care must be taken, though, to ignore
       cosmic-rays as much as possible as they will mimic real sources
       without providing information on the actual focus through the optics.
       Similarly, saturation regions must also be ignored as they also
       only indicate a detector feature, not the focus through the optics or
       alignment of actual sources.

    """

    img_log = np.abs(ndimage.gaussian_laplace(img, sigma))
    focus_val = img_log.max()
    focus_pos = np.where(img_log == focus_val)

    return focus_val, focus_pos


def compute_zero_mask(imgarr, iterations=8, ext=0):
    """Find section from image with no masked out pixels and max total flux"""
    if isinstance(imgarr, str):
        img_mask = fits.getdata(imgarr, ext=0)
    else:
        img_mask = imgarr.copy()

    img_mask[img_mask > 0] = 1
    img_mask = ndimage.binary_erosion(img_mask, iterations=iterations)

    return img_mask


def build_focus_dict(singlefiles, prodfile, sigma=2.0):

    focus_dict = FOCUS_DICT.copy()
    focus_dict['expnames'] = singlefiles
    focus_dict['prodname'] = prodfile

    # Start by creating the full saturation mask from all single_sci images
    full_sat_mask = None
    for f in singlefiles:
        sat_mask = compute_zero_mask(f)
        if full_sat_mask is None:
            full_sat_mask = sat_mask
        else:
            full_sat_mask = np.bitwise_and(full_sat_mask, sat_mask)

    # Now apply full saturation mask to each single_sci image and compute focus
    for f in singlefiles:
        imgarr = fits.getdata(f)
        imgarr[~full_sat_mask] = 0
        focus_val, focus_pos = determine_focus_index(imgarr, sigma=sigma)
        focus_dict['exp'].append(float(focus_val))
        focus_dict['exp_pos'] = (int(focus_pos[0][0]), int(focus_pos[1][0]))

    # Generate results for drizzle product(s)
    prodarr = fits.getdata(prodfile)
    prodarr[~full_sat_mask] = 0
    # Insure output values are JSON-compliant
    focus_val, focus_pos = determine_focus_index(prodarr, sigma=sigma)
    focus_dict['prod'].append(float(focus_val))
    focus_dict['prod_pos'] = (int(focus_pos[0][0]), int(focus_pos[1][0]))

    # Determine statistics for evalaution
    exparr = np.array(focus_dict['exp'])
    focus_dict['stats'] = {'mean': exparr.mean(), 'std': exparr.std(),
                           'min': exparr.min(), 'max': exparr.max()}
    log.debug("Focus results for {}: \n{}".format(prodfile, focus_dict))
    log.info("Mean Focus computed for {}: {}".format(prodfile, focus_dict['stats']['mean']))

    return focus_dict


def evaluate_focus(focus_dict, tolerance=0.8):
    if focus_dict is None:
        return True

    s = focus_dict['stats']
    min_3sig = min(s['mean'] - 3.0 * s['std'], tolerance * s['mean'])
    max_3sig = s['mean'] + 3.0 * s['std']
    min_prob = compute_prob(min_3sig, s['mean'], s['std'])
    max_prob = compute_prob(max_3sig, s['mean'], s['std'])
    drz_prob = np.array([compute_prob(d, s['mean'], s['std']) for d in focus_dict['prod']])

    if (drz_prob < min_prob).any() or (drz_prob > max_prob).any() or s['std'] > s['min']:
        alignment_verified = False
    else:
        alignment_verified = True

    return alignment_verified


def get_align_fwhm(focus_dict, default_fwhm, src_size=11):
    """Determine FWHM based on position of sharpest focus in the product"""
    pimg = fits.open(focus_dict['prodname'])
    posy, posx = focus_dict['prod_pos']

    prod = pimg[1].data if len(pimg) > 1 else pimg[0].data

    src = prod[posy - src_size:posy + src_size, posx - src_size:posx + src_size]

    # For sources near the edge of the image data, insure that any NaN's are converted to 0
    # This is necessary in order to allow FWHM to be determined
    src = np.nan_to_num(src, nan=0)

    # Normalize to total flux of 1 for FWHM determination
    kernel = src / src.sum()

    fwhm = find_fwhm(kernel, default_fwhm)
    # Be nice and close the FITS image
    pimg.close()

    return fwhm


def max_overlap_diff(total_mask, singlefiles, prodfile, sigma=2.0, scale=1, lsigma=3.0):
    """Determines the difference in the region of max overlap for all drizzled products

    Parameters
    -----------
    total_mask : ndarray
        Mask (array) showing where each input exposure contributes to the final
        drizzle product `prodfile`.  This could be created using
        `cell_utils.SkyFootprint`.

    singlefiles : list
        List of filenames for each single input exposure drizzled onto the same WCS
        as the final drizzle product `prodfile`

    prodfile : str
        Filename for the final drizzle product

    scale : int, optional
        Factor to use in downsizing (resizing smaller) the images to be evaluated.  The larger
        the value, the less sensitive this measurement becomes.

    sigma : float, optional
        Size of default kernel (in pixels) to use for determining the focus.

    Returns
    ---------
    diff_dict : dictionary
        Dictionary of difference scores for each input exposure drizzle product
        (from `singlefiles`) calculated for the region of maximum overlap with the final
        drizzle product `prodfile`.  Entries for each singlefile includes:

        ===================   =============================================
        distance              Hamming distance of singlefile from prodfile
        focus                 focus index of singlefile
        focus_pos             position for best focus in singlefile
        product_focus         focus index for prodfile
        product_focus_pos     position for best focus in prodfile
        ===================   =============================================

    """

    drz = fits.getdata(prodfile, ext=("SCI", 1))
    # Verify that the total_mask has the same dimensions as the drz/single_file images
    if drz.shape != total_mask.shape:
        log.error("Total mask shape {} needs to be the same as input image's shape \
                    {}".format(total_mask.shape, drz.shape))
        raise ValueError

    # Determine regions of overlap in total mask
    min_overlap = total_mask > 1
    max_overlap = total_mask == total_mask.max()

    exptimes = np.array([fits.getval(s, 'exptime') for s in singlefiles])
    exp_weights = exptimes / exptimes.max()
    log.info("Computing diffs for: {}".format(singlefiles))

    diff_dict = {}
    for sfile, exp_weight in zip(singlefiles, exp_weights):
        # start by seeing whether this product overlaps the region of max_overlap
        sdata = fits.getdata(sfile)

        # Create exposure mask corresponding to pixels with drizzled data
        smask = sdata > 0
        if smask.sum() == 0:
            # In some error cases (e.g., jcx552010), blank images get to this point, so treat them as blank
            log.info("Overlap difference for {}: (No valid data)".format(sfile))
            diff_dict[sfile] = {"distance": -1, "xslice": None, "yslice": None}
            diff_dict[sfile]['product_num_sources'] = 0
            diff_dict[sfile]['num_sources'] = 0
            diff_dict[sfile]['focus'] = 0.0
            diff_dict[sfile]['focus_pos'] = (None, None)
            diff_dict[sfile]['product_focus'] = 0.0
            diff_dict[sfile]['product_focus_pos'] = (None, None)
            log.debug("Overlap differences for {} found to be: \n{}".format(sfile, diff_dict[sfile]))
            continue

        # Trim mask down to only include region where the most exposures overlap
        soverlap = smask * max_overlap

        # If, for some reason, the exposure does not overlap the region of
        # max overlap (for example, in a large mosaic with little overlap)
        # resort to using area where single exposure overlaps at least 1 other
        # exposure instead...
        if soverlap.sum() == 0:
            # Use this for computing the difference index
            soverlap = smask * min_overlap

        # get same region from each drizzle product
        drz_region = drz * soverlap
        sfile_region = sdata * soverlap

        # Insure all np.nan's are converted to zeros
        drz_region = np.nan_to_num(drz_region, nan=0)
        sfile_region = np.nan_to_num(sfile_region, nan=0)

        # Also compute focus index for the same region of the single drizzle file
        focus_val, focus_pos = determine_focus_index(sfile_region, sigma=sigma)
        pfocus_val, pfocus_pos = determine_focus_index(drz_region, sigma=sigma)

        # Limit our analysis only to those pixels within the masked region
        #  (modulo slicing limits)
        yr, xr = np.where(soverlap > 0)
        yslice = slice(yr.min(), yr.max(), 1)
        xslice = slice(xr.min(), xr.max(), 1)
        log.debug("overlap region: xslice {}, yslice {}".format(xslice, yslice))

        drz_arr = drz_region[yslice, xslice]
        sfile_arr = sfile_region[yslice, xslice]

        # The number of sources detected is subject to crowding/blending of sources
        # as well as noise from the background (if too low
        #  a background value is used)
        drzlabels, drznum = detect_point_sources(drz_arr, scale=scale, log_sigma=lsigma)
        slabels, snum = detect_point_sources(sfile_arr, scale=scale, exp_weight=exp_weight, log_sigma=lsigma)

        drzsrcs = np.clip(drzlabels, 0, 1).astype(np.int16)
        sfilesrcs = np.clip(slabels, 0, 1).astype(np.int16)

        # Determine number of nonzero pixels being measured in 'truth'/single image
        sfile_num = np.nonzero(sfilesrcs)[0].size

        # Compute distance between difference scores for 'truth' and 'product'
        # and weight by fraction of nonzero pixels in region
        # This produces the HAMMING distance for the two arrays
        # dist = (np.abs(drzsrcs - sfilesrcs).sum() / drz_arr.size) * weight
        dist = (np.abs(drzsrcs - sfilesrcs).sum() / sfile_num) * exp_weight

        # Compute similarity_index for these overlapping regions to compare with hamming dist
        sim = compute_similarity(sfile_arr, drz_arr)
        # Take the min since the Hamming distance can be skewed by bad-pixels and noise
        # more than the similarity_index.
        # similarity_index is scaled by 2 to be scaled the same as the Hamming distance
        dist = min(dist, sim * 2.0)

        # Record results for each exposure compared to the combined drizzle product
        # Number of sources in drz and sfile can include artifacts such as CRs
        # As a result, care must be taken in any comparisons using these values.
        log.info("Overlap difference for {}: {:0.4f}".format(sfile, dist))
        diff_dict[sfile] = {"distance": dist, "xslice": xslice, "yslice": yslice}
        diff_dict[sfile]['product_num_sources'] = drznum
        diff_dict[sfile]['num_sources'] = snum
        diff_dict[sfile]['focus'] = float(focus_val)
        diff_dict[sfile]['focus_pos'] = (int(focus_pos[0][0]), int(focus_pos[1][0]))
        diff_dict[sfile]['product_focus'] = float(pfocus_val)
        diff_dict[sfile]['product_focus_pos'] = (int(pfocus_pos[0][0]), int(pfocus_pos[1][0]))
        log.debug("Overlap differences for {} found to be: \n{}".format(sfile, diff_dict[sfile]))

    return diff_dict


def sigma_clipped_bkg(arr, sigma=3.0, nsigma=4, maxiters=None):
    # Account for input being blank
    if arr.max() == 0:
        return 0.0, [0.0, 0.0, 0.0]
    if maxiters is None:
        maxiters = int(np.log10(arr.max() / 2) + 0.5)

    # Use simple constant background to avoid problems with nebulosity
    bkg = sigma_clipped_stats(arr, sigma=sigma, maxiters=maxiters)
    bkg_total = bkg[0] + nsigma * bkg[2]  # mean + 4 * sigma

    return bkg_total, bkg


def reduce_diff_region(arr, scale=1, background=None, nsigma=4,
                       sigma=3.0, exp_weight=None):
    """Convert the image into a background-removed array"""
    # Provide option to rebin to a smaller image size to minimize
    # impact from high-frequency (pixel-to-pixel) differences in low S/N data
    if scale > 1:
        yend = arr.shape[0] % scale
        xend = arr.shape[1] % scale
        yend = -1 * yend if yend > 0 else None
        xend = -1 * xend if xend > 0 else None
        new_shape = (arr.shape[0] // scale, arr.shape[1] // scale)

        rebin_arr = rebin(arr[:yend, :xend].copy(), new_shape)
    else:
        rebin_arr = arr.copy()

    if background is None:
        """
        if exp_weight is not None and 0.2 >= exp_weight >= EXP_LIMIT:
            sigma = 2.0
            maxiters = 10.
        """
        if exp_weight is not None and exp_weight < EXP_RATIO:
            if EXP_RATIO >= exp_weight >= EXP_LIMIT:
                sigma = 3.0
            elif exp_weight < EXP_LIMIT:
                sigma = 3.
            else:
                pass

        bkg_total, bkg = sigma_clipped_bkg(rebin_arr, sigma=sigma, nsigma=nsigma)
        log.debug("sigma clipped background value: {}".format(bkg_total))
        blank_image = True if (bkg[1] < bkg[2] and bkg[1] < 1.0) else False

    elif isinstance(background, Background2D):
        bkg_total = background.background + nsigma * background.background_rms
        log.debug("background: max={}, mean={}".format(bkg_total.max(), bkg_total.mean()))
        blank_image = True if (background.median < background.median_rms and
                               background.median < 1.0) else False

    if blank_image:
        # median filter image to limit noise-induced variations into overlap differences
        rebin_arr = ndimage.median_filter(rebin_arr, size=5)

    rebin_arr -= bkg_total
    rebin_arr = np.clip(rebin_arr, 0, rebin_arr.max())

    return rebin_arr


def detect_point_sources(arr, background=None, nsigma=4, log_sigma=3.0, scale=1,
                         sigma=3.0, exp_weight=None):
    # Remove background entirely from input array (clip at 0)
    src_arr = reduce_diff_region(arr, background=background, nsigma=nsigma, scale=scale,
                                 sigma=sigma, exp_weight=exp_weight)

    # Compute distance between images using labeled sources
    srclog = -1 * ndimage.gaussian_laplace(src_arr, sigma=log_sigma)

    # zero out wings of sources, only leaving the detected cores/edges...
    srclog[srclog < 0] = 0

    # label sources
    slabels, snum = ndimage.label(srclog)

    return slabels, snum


def diff_score(arr):

    # Convert arrays into 1D arrays along rows and columns, respectively
    rows = arr.flatten()
    cols = arr.flatten("F")
    # Compute pixel-to-pixel differences along row and columns, respectively
    # and convert to boolean result (delta > 0 is True/0, delta < 0 is False/1)
    diff_row = np.diff(rows) > 0
    diff_col = np.diff(cols) > 0
    # Stack row and column 1D array as a single concatenated result
    return np.hstack((diff_row, diff_col)).flatten()


def evaluate_overlap_diffs(diff_dict, limit=1.0):
    """Evaluate whether overlap diffs indicate good alignment or not. """

    max_diff = max([d['distance'] for d in diff_dict.values()])
    verified = max_diff <= limit
    log.info("Maximum overlap difference: {:0.4f}".format(max_diff))
    if verified:
        log.info("Alignment verified based on overlap...")
    else:
        log.info("Alignment NOT verified based on overlap...")

    return verified, max_diff
