"""This script contains code to support creation of photometric sourcelists using two techniques:
aperture photometry and segmentation-map based photometry."""

import os
import sys
import shutil
import warnings

import numpy as np

from astropy.io import fits as fits
from photutils.detection import findstars

from stsci.tools import fileutil as fu

from .. import astrodrizzle

#
# Original imports
#
import copy
import pickle  # FIX Remove

import astropy.units as u
from astropy.stats import sigma_clipped_stats, gaussian_fwhm_to_sigma
from astropy.table import Column, MaskedColumn, Table, join, vstack
from astropy.convolution import RickerWavelet2DKernel
from astropy.coordinates import SkyCoord
from scipy import ndimage

from photutils import CircularAperture, CircularAnnulus, DAOStarFinder, IRAFStarFinder
from photutils import Background2D, SExtractorBackground, StdBackgroundRMS
from photutils import detect_sources, source_properties, deblend_sources
from photutils import make_source_mask
from photutils.utils import calc_total_error
from stsci.tools import logutil
from stwcs.wcsutil import HSTWCS

from . import astrometric_utils
from . import photometry_tools

try:
    from matplotlib import pyplot as plt
except Exception:
    plt = None

CATALOG_TYPES = ['aperture', 'segment']

PSF_PATH = ['pars', 'psfs']

__taskname__ = 'deconvolve_utils'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)


# ======================================================================================================================

def fft_deconv_img(img, psf, freq_limit=0.95):
    """ FFT image deconvolution

    PARAMETERS
    -----------
    img : ndarray
        Numpy array of image to be deconvolved

    psf : ndarray
        Numpy array of PSF to be used for deconvolution.  This array has to have the same shape as
        the input image and should be normalized to a total flux (sum) approximately equal to the
        brightest non-saturated source in the science image being deconvolved.

    freq_limit : float
        Compute regularization parameter (frequency limit) to be used with PSF in deconvolution.
        This value should result in selecting a frequency
        one to two orders of magnitude below the largest spectral component of the point-spread function.

    RETURNS
    -------
    deconv : ndarray
        Numpy array of deconvolved image

    .. note::
        Based on 2017 implementation by Juha Vierinen
        http://www.radio-science.net/2017/09/deconvolution-in-frequency-domain-with.html
    """
    # Insure PSF has no NaN values and is centred in output array, to avoid shifting deconvolved product
    # This also has the advantage of insuring we are only working on a copy of the input PSF
    psf = np.nan_to_num(psf, 0.0)

    # Compute alpha scaling based on PSF frequency limit
    P2 = np.abs(np.copy(psf).flatten())
    P2 = np.sort(P2)[::-1]
    index = int(P2.shape[0] * freq_limit)
    alpha = P2[index]
    del P2

    psf_y, psf_x = np.where(psf == psf.max())[0]
    psf_center = [psf.shape[0] // 2, psf.shape[1] // 2]
    if int(psf_y) != psf_center[0] or int(psf_x) != psf_center[1]:
        # We need to recenter the PSF
        psf_section = np.where(psf != 0.0)
        psf_xr = [psf_section[1].min(), psf_section[1].max()]
        psf_yr = [psf_section[0].min(), psf_section[0].max()]
        psf_size = [(psf_yr[1] - psf_yr[0]) // 2, (psf_xr[1] - psf_xr[0]) // 2]
        # If psf_size is not at least 2 * psf_max//2, increase to that value
        # This will insure that the final position of the PSF in the output is at the exact center
        psf_max = np.where(psf[psf_yr[0]:psf_yr[1], psf_xr[0]:psf_yr[1]] == psf.max())[0]
        psf_len = max(max(psf_max[0] // 2, psf_size[0]), max(psf_max[1] // 2, psf_size[1]))

        centered_psf = psf * 0.0
        centered_psf[psf_center[0] - psf_len: psf_center[0] + psf_len,
                     psf_center[1] - psf_len: psf_center[1] + psf_len] = psf[psf_y - psf_len: psf_y + psf_len,
                                                                             psf_x - psf_len: psf_x + psf_len]
        psf = centered_psf

    # Insure input image also has no NaN values
    img = np.nan_to_num(img, 0.0)

    # FFT point spread function (first column of theory matrix G)
    P = np.fft.fft2(psf)

    # FFT2 measurement
    # Use image in husky_conv.png
    # U^H d
    D = np.fft.fft2(img)

    # -dampped spectral components,
    # -also known as Wiener filtering
    # (conj(S)/(|S|^2 + alpha^2)) U^H d
    M = (np.conj(P) / (np.abs(P)**2.0 + alpha**2.0)) * D

    # maximum a posteriori estimate of deconvolved image
    # m_map = U (conj(S)/(|S|^2 + alpha^2)) U^H d
    m_map = (D.shape[1] * D.shape[0]) * np.fft.fftshift(np.fft.ifft2(M).real)

    return m_map

# ======================================================================================================================

# Functions to manage PSF library for deconvolution

def find_psf(imgname, instr=None, detector=None, filter=None, path_root=None,
             instr_kws=['INSTRUME', 'DETECTOR', 'FILTER']):
    """Pull PSF from library based on unique combination of intrument/detector/filter.

    Parameters
    ===========
    imgname : str
        Image name of science image to be deconvolved.  If the header contains instrument, detector, and filters,
        then those additional parameters do not need to be specified.

    instr : str, optional
        If not included in `imgname`, specify instrument for PSF.

    detector : str, optional
        If not included in `imgname`, specify detector for PSF.

    filter : list, optional
        If not included in `imgname`, specify detector for PSF.

    path_root : str, optional
        Full path to parent directory of PSF library IF not the default path for package.

    instr_kws : list
        List of keywords from input image header where the instrument, detector and filter
        names are specified.

    Returns
    ========
    psfname : str
        Full name, with path, for PSF from library.

    """
    # look for instrument, detector, and filter in input image name
    # Start by looking in the input image header for these values, so we know what to look
    kw_vals = [fits.getval(imgname, kw).lower() for kw in instr_kws]

    if len(kw_vals) < len(instr_kws):
        kw_vals = [instr, detector] + filter

    if kw_vals[0] is None:
        log.error("No valid keywords found.")
        raise ValueError

    if path_root is None:
        path_root = os.path.split(os.path.dirname(__file__))[0]
        for psf_path in PSF_PATH:
            path_root = os.path.join(path_root, psf_path)

    path_root = os.path.join(path_root, kw_vals[0], kw_vals[1])

    psf_name = os.path.join(path_root, "_".join(kw_vals + ['psf.fits']))

    log.debug('Looking for Library PSF {}'.format(psf_name))

    if not os.path.exists(psf_name):
        log.error('No PSF found for keywords {} \n   with values of {}'.format(instr_kws, kw_vals))
        raise ValueError

    return psf_name


def convert_library_psf(calimg, drzimg, psf, scaling, smoothing=1.0):
    """Drizzle library PSF to match science image. """

    # Create copy of input science image based on input psf filename
    psf_root = os.path.basename(psf)
    lib_psf_arr = fits.getdata(psf)
    lib_psf_arr *= scaling

    lib_size = [lib_psf_arr.shape[0] // 2, lib_psf_arr.shape[1] // 2]

    # This will be the name of the new file containing the library PSF that will be drizzled to
    # match the input iamge `drzimg`
    psf_flt_name = psf_root.replace('psf.fits', 'psf_flt.fits')
    psf_drz_name = psf_flt_name.replace('flt.fits', '')  # astrodrizzle will add suffix

    # create version of PSF that will be drizzled
    psf_base = fits.getdata(calimg, ext=1) * 0.0
    # Copy library PSF into this array
    out_cen = [psf_base.shape[0] // 2, psf_base.shape[1] // 2]
    psf_base[out_cen[0] - lib_size[0]: out_cen[0] + lib_size[0],
             out_cen[1] - lib_size[1]: out_cen[1] + lib_size[1]] = lib_psf_arr

    # Write out library PSF FLT file now
    psf_flt = shutil.copy(calimg, psf_flt_name)

    # Update file with library PSF
    flt_hdu = fits.open(psf_flt, mode='update')
    flt_hdu[('sci', 1)].data = psf_base
    num_sci = fu.countExtn(calimg)
    # Also zero out all other science data in this 'PSF' file.
    if num_sci > 1:
        for extn in range(2, num_sci + 1):
            flt_hdu[('sci', extn)].data *= 0.0
    flt_hdu.close()
    del flt_hdu

    # Insure only final drizzle step is run
    drizzle_pars = {}
    drizzle_pars["build"] = True
    drizzle_pars['context'] = False
    drizzle_pars['clean'] = True
    drizzle_pars['in_memory'] = True
    drizzle_pars["resetbits"] = 0

    drizzle_pars["static"] = False
    drizzle_pars["skysub"] = False
    drizzle_pars["driz_separate"] = False
    drizzle_pars["median"] = False
    drizzle_pars["blot"] = False
    drizzle_pars["driz_cr"] = False
    drizzle_pars["driz_combine"] = True
    drizzle_pars['final_wcs'] = True
    drizzle_pars['final_pixfrace'] = smoothing
    drizzle_pars["final_refimage"] = "{}[1]".format(drzimg)

    # Drizzle PSF FLT file to match orientation and plate scale of drizzled science (total detection) image
    astrodrizzle.AstroDrizzle(input=psf_flt_name,
                              output=psf_drz_name,
                              **drizzle_pars)

    return psf_flt_name.replace("flt.fits", "drz.fits")


def get_cutouts(data, star_list, kernel, threshold_eff, exclude_border=False):

    coords = [(row[1], row[0]) for row in star_list]
    convolved_data = findstars._filter_data(data, kernel.data, mode='constant',
                                            fill_value=0.0, check_normalization=False)

    star_cutouts = []
    for (ypeak, xpeak) in coords:
        # now extract the object from the data, centered on the peak
        # pixel in the convolved image, with the same size as the kernel
        x0 = xpeak - kernel.xradius
        x1 = xpeak + kernel.xradius + 1
        y0 = ypeak - kernel.yradius
        y1 = ypeak + kernel.yradius + 1

        if x0 < 0 or x1 > data.shape[1]:
            continue  # pragma: no cover
        if y0 < 0 or y1 > data.shape[0]:
            continue  # pragma: no cover

        slices = (slice(y0, y1), slice(x0, x1))
        data_cutout = data[slices]
        convdata_cutout = convolved_data[slices]

        # correct pixel values for the previous image padding
        if not exclude_border:
            x0 -= kernel.xradius
            x1 -= kernel.xradius
            y0 -= kernel.yradius
            y1 -= kernel.yradius
            xpeak -= kernel.xradius
            ypeak -= kernel.yradius
            slices = (slice(y0, y1), slice(x0, x1))

        star_cutouts.append(findstars._StarCutout(data_cutout, convdata_cutout, slices,
                                        xpeak, ypeak, kernel, threshold_eff))

        return star_cutouts


class UserStarFinder(findstars.StarFinderBase):
    """
    Measure stars in an image using the DAOFIND (`Stetson 1987
    <https://ui.adsabs.harvard.edu/abs/1987PASP...99..191S/abstract>`_)
    algorithm.  Stars measured using DAOFIND can be identified using
    any algorithm defined by the user, with the results passed in as a
    simple list of coords.

    DAOFIND (`Stetson 1987; PASP 99, 191
    <https://ui.adsabs.harvard.edu/abs/1987PASP...99..191S/abstract>`_)
    searches images for local density maxima that have a peak amplitude
    greater than ``threshold`` (approximately; ``threshold`` is applied
    to a convolved image) and have a size and shape similar to the
    defined 2D Gaussian kernel.  The Gaussian kernel is defined by the
    ``fwhm``, ``ratio``, ``theta``, and ``sigma_radius`` input
    parameters.

    ``DAOStarFinder`` finds the object centroid by fitting the marginal x
    and y 1D distributions of the Gaussian kernel to the marginal x and
    y distributions of the input (unconvolved) ``data`` image.

    ``DAOStarFinder`` calculates the object roundness using two methods. The
    ``roundlo`` and ``roundhi`` bounds are applied to both measures of
    roundness.  The first method (``roundness1``; called ``SROUND`` in
    `DAOFIND`_) is based on the source symmetry and is the ratio of a
    measure of the object's bilateral (2-fold) to four-fold symmetry.
    The second roundness statistic (``roundness2``; called ``GROUND`` in
    `DAOFIND`_) measures the ratio of the difference in the height of
    the best fitting Gaussian function in x minus the best fitting
    Gaussian function in y, divided by the average of the best fitting
    Gaussian functions in x and y.  A circular source will have a zero
    roundness.  A source extended in x or y will have a negative or
    positive roundness, respectively.

    The sharpness statistic measures the ratio of the difference between
    the height of the central pixel and the mean of the surrounding
    non-bad pixels in the convolved image, to the height of the best
    fitting Gaussian function at that point.

    Parameters
    ----------
    threshold : float
        The absolute image value above which to select sources.

    fwhm : float
        The full-width half-maximum (FWHM) of the major axis of the
        Gaussian kernel in units of pixels.

    ratio : float, optional
        The ratio of the minor to major axis standard deviations of the
        Gaussian kernel.  ``ratio`` must be strictly positive and less
        than or equal to 1.0.  The default is 1.0 (i.e., a circular
        Gaussian kernel).

    theta : float, optional
        The position angle (in degrees) of the major axis of the
        Gaussian kernel measured counter-clockwise from the positive x
        axis.

    sigma_radius : float, optional
        The truncation radius of the Gaussian kernel in units of sigma
        (standard deviation) [``1 sigma = FWHM /
        (2.0*sqrt(2.0*log(2.0)))``].

    sharplo : float, optional
        The lower bound on sharpness for object detection.

    sharphi : float, optional
        The upper bound on sharpness for object detection.

    roundlo : float, optional
        The lower bound on roundness for object detection.

    roundhi : float, optional
        The upper bound on roundness for object detection.

    sky : float, optional
        The background sky level of the image.  Setting ``sky`` affects
        only the output values of the object ``peak``, ``flux``, and
        ``mag`` values.  The default is 0.0, which should be used to
        replicate the results from `DAOFIND`_.

    exclude_border : bool, optional
        Set to `True` to exclude sources found within half the size of
        the convolution kernel from the image borders.  The default is
        `False`, which is the mode used by `DAOFIND`_.

    brightest : int, None, optional
        Number of brightest objects to keep after sorting the full object list.
        If ``brightest`` is set to `None`, all objects will be selected.

    peakmax : float, None, optional
        Maximum peak pixel value in an object. Only objects whose peak pixel
        values are *strictly smaller* than ``peakmax`` will be selected.
        This may be used to exclude saturated sources. By default, when
        ``peakmax`` is set to `None`, all objects will be selected.

        .. warning::
            `DAOStarFinder` automatically excludes objects whose peak
            pixel values are negative. Therefore, setting ``peakmax`` to a
            non-positive value would result in exclusion of all objects.

    See Also
    --------
    IRAFStarFinder

    Notes
    -----
    For the convolution step, this routine sets pixels beyond the image
    borders to 0.0.  The equivalent parameters in `DAOFIND`_ are
    ``boundary='constant'`` and ``constant=0.0``.

    The main differences between `~photutils.detection.DAOStarFinder`
    and `~photutils.detection.IRAFStarFinder` are:

    * `~photutils.detection.IRAFStarFinder` always uses a 2D
      circular Gaussian kernel, while
      `~photutils.detection.DAOStarFinder` can use an elliptical
      Gaussian kernel.

    * `~photutils.detection.IRAFStarFinder` calculates the objects'
      centroid, roundness, and sharpness using image moments.

    References
    ----------
    .. [1] Stetson, P. 1987; PASP 99, 191
           (https://ui.adsabs.harvard.edu/abs/1987PASP...99..191S/abstract)
    .. [2] https://iraf.net/irafhelp.php?val=daofind

    .. _DAOFIND: https://iraf.net/irafhelp.php?val=daofind
    """

    def __init__(self, threshold, fwhm, ratio=1.0, theta=0.0,
                 sigma_radius=1.5, sharplo=0.2, sharphi=1.0, roundlo=-1.0,
                 roundhi=1.0, sky=0.0, exclude_border=False,
                 brightest=None, peakmax=None):

        if not np.isscalar(threshold):
            raise TypeError('threshold must be a scalar value.')
        self.threshold = threshold

        if not np.isscalar(fwhm):
            raise TypeError('fwhm must be a scalar value.')
        self.fwhm = fwhm

        self.ratio = ratio
        self.theta = theta
        self.sigma_radius = sigma_radius
        self.sharplo = sharplo
        self.sharphi = sharphi
        self.roundlo = roundlo
        self.roundhi = roundhi
        self.sky = sky
        self.exclude_border = exclude_border

        self.kernel = findstars._StarFinderKernel(self.fwhm, self.ratio, self.theta,
                                        self.sigma_radius)
        self.threshold_eff = self.threshold * self.kernel.relerr
        self.brightest = brightest
        self.peakmax = peakmax
        self._star_cutouts = None

    def find_stars(self, data, coords=None, mask=None):
        """
        Find stars in an astronomical image.

        Parameters
        ----------
        data : 2D array_like
            The 2D image array.

        mask : 2D bool array, optional
            A boolean mask with the same shape as ``data``, where a
            `True` value indicates the corresponding element of ``data``
            is masked.  Masked pixels are ignored when searching for
            stars.

        coords : `~astropy.table.Table` or `None`
            A table, such as returned by `find_peaks`, with approximate X,Y positions of identified sources
            If not provided, the DAOFind algorithm will be used to find sources.

        Returns
        -------
        table : `~astropy.table.Table` or `None`
            A table of found stars with the following parameters:

            * ``id``: unique object identification number.
            * ``xcentroid, ycentroid``: object centroid.
            * ``sharpness``: object sharpness.
            * ``roundness1``: object roundness based on symmetry.
            * ``roundness2``: object roundness based on marginal Gaussian
              fits.
            * ``npix``: the total number of pixels in the Gaussian kernel
              array.
            * ``sky``: the input ``sky`` parameter.
            * ``peak``: the peak, sky-subtracted, pixel value of the object.
            * ``flux``: the object flux calculated as the peak density in
              the convolved image divided by the detection threshold.  This
              derivation matches that of `DAOFIND`_ if ``sky`` is 0.0.
            * ``mag``: the object instrumental magnitude calculated as
              ``-2.5 * log10(flux)``.  The derivation matches that of
              `DAOFIND`_ if ``sky`` is 0.0.

            `None` is returned if no stars are found.

        """

        if coords:
            star_cutouts = get_cutouts(data, coords, kernel=self.kernel)
        else:
            star_cutouts = findstars._find_stars(data, self.kernel, self.threshold_eff,
                                                 mask=mask,
                                                 exclude_border=self.exclude_border)

        if star_cutouts is None:
            warnings.warn('No sources were found.', findstars.NoDetectionsWarning)
            return None

        self._star_cutouts = star_cutouts

        star_props = []
        for star_cutout in star_cutouts:
            props = findstars._DAOFindProperties(star_cutout, self.kernel, self.sky)

            if np.isnan(props.dx_hx).any() or np.isnan(props.dy_hy).any():
                continue

            if (props.sharpness <= self.sharplo or
                    props.sharpness >= self.sharphi):
                continue

            if (props.roundness1 <= self.roundlo or
                    props.roundness1 >= self.roundhi):
                continue

            if (props.roundness2 <= self.roundlo or
                    props.roundness2 >= self.roundhi):
                continue

            if self.peakmax is not None and props.peak >= self.peakmax:
                continue

            star_props.append(props)

        nstars = len(star_props)
        if nstars == 0:
            warnings.warn('Sources were found, but none pass the sharpness '
                          'and roundness criteria.', findstars.NoDetectionsWarning)
            return None

        if self.brightest is not None:
            fluxes = [props.flux for props in star_props]
            idx = sorted(np.argsort(fluxes)[-self.brightest:].tolist())
            star_props = [star_props[k] for k in idx]
            nstars = len(star_props)

        table = Table()
        table['id'] = np.arange(nstars) + 1
        columns = ('xcentroid', 'ycentroid', 'sharpness', 'roundness1',
                   'roundness2', 'npix', 'sky', 'peak', 'flux', 'mag')
        for column in columns:
            table[column] = [getattr(props, column) for props in star_props]

        return table
