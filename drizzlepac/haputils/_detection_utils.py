# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module implements common utility functions and classes for the star
finders. This is a copy of photutils private functions, which have been
refactored. Use of these tools can likely be removed after drizzlepac
requires photutils >= 1.2.0.
"""

import math
import warnings

from astropy.convolution import Kernel2D
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.units import Quantity
from astropy.utils import lazyproperty
from astropy.utils.exceptions import AstropyUserWarning
import numpy as np
from photutils.detection import find_peaks
from photutils.utils.exceptions import NoDetectionsWarning


def _filter_data(data, kernel, mode='constant', fill_value=0.0,
                 check_normalization=False):
    """
    Convolve a 2D image with a 2D kernel.

    The kernel may either be a 2D `~numpy.ndarray` or a
    `~astropy.convolution.Kernel2D` object.

    Parameters
    ----------
    data : array_like
        The 2D array of the image.

    kernel : array-like (2D) or `~astropy.convolution.Kernel2D`
        The 2D kernel used to filter the input ``data``. Filtering the
        ``data`` will smooth the noise and maximize detectability of
        objects with a shape similar to the kernel.

    mode : {'constant', 'reflect', 'nearest', 'mirror', 'wrap'}, optional
        The ``mode`` determines how the array borders are handled.  For
        the ``'constant'`` mode, values outside the array borders are
        set to ``fill_value``.  The default is ``'constant'``.

    fill_value : scalar, optional
        Value to fill data values beyond the array borders if ``mode``
        is ``'constant'``.  The default is ``0.0``.

    check_normalization : bool, optional
        If `True` then a warning will be issued if the kernel is not
        normalized to 1.
    """
    from scipy import ndimage

    if kernel is not None:
        if isinstance(kernel, Kernel2D):
            kernel_array = kernel.array
        else:
            kernel_array = kernel

        if check_normalization:
            if not np.allclose(np.sum(kernel_array), 1.0):
                warnings.warn('The kernel is not normalized.',
                              AstropyUserWarning)

        # scipy.ndimage.convolve currently strips units, but be explicit
        # in case that behavior changes
        unit = None
        if isinstance(data, Quantity):
            unit = data.unit
            data = data.value

        # NOTE:  astropy.convolution.convolve fails with zero-sum
        # kernels (used in findstars) (cf. astropy #1647)
        # NOTE: if data is int and kernel is float, ndimage.convolve
        # will return an int image - here we make the data float so
        # that a float image is always returned
        result = ndimage.convolve(data.astype(float), kernel_array,
                                  mode=mode, cval=fill_value)

        if unit is not None:
            result = result * unit  # can't use *= with older astropy

        return result
    else:
        return data


class _StarFinderKernel:
    """
    Class to calculate a 2D Gaussian density enhancement kernel.

    The kernel has negative wings and sums to zero.  It is used by both
    `DAOStarFinder` and `IRAFStarFinder`.

    Parameters
    ----------
    fwhm : float
        The full-width half-maximum (FWHM) of the major axis of the
        Gaussian kernel in units of pixels.

    ratio : float, optional
        The ratio of the minor and major axis standard deviations of the
        Gaussian kernel.  ``ratio`` must be strictly positive and less
        than or equal to 1.0.  The default is 1.0 (i.e., a circular
        Gaussian kernel).

    theta : float, optional
        The position angle (in degrees) of the major axis of the
        Gaussian kernel, measured counter-clockwise from the positive x
        axis.

    sigma_radius : float, optional
        The truncation radius of the Gaussian kernel in units of sigma
        (standard deviation) [``1 sigma = FWHM /
        2.0*sqrt(2.0*log(2.0))``].  The default is 1.5.

    normalize_zerosum : bool, optional
        Whether to normalize the Gaussian kernel to have zero sum, The
        default is `True`, which generates a density-enhancement kernel.

    Notes
    -----
    The class attributes include the dimensions of the elliptical kernel
    and the coefficients of a 2D elliptical Gaussian function expressed
    as:

        ``f(x,y) = A * exp(-g(x,y))``

        where

        ``g(x,y) = a*(x-x0)**2 + 2*b*(x-x0)*(y-y0) + c*(y-y0)**2``

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Gaussian_function
    """

    def __init__(self, fwhm, ratio=1.0, theta=0.0, sigma_radius=1.5,
                 normalize_zerosum=True):

        if fwhm < 0:
            raise ValueError('fwhm must be positive.')

        if ratio <= 0 or ratio > 1:
            raise ValueError('ratio must be positive and less or equal '
                             'than 1.')

        if sigma_radius <= 0:
            raise ValueError('sigma_radius must be positive.')

        self.fwhm = fwhm
        self.ratio = ratio
        self.theta = theta
        self.sigma_radius = sigma_radius
        self.xsigma = self.fwhm * gaussian_fwhm_to_sigma
        self.ysigma = self.xsigma * self.ratio

        theta_radians = np.deg2rad(self.theta)
        cost = np.cos(theta_radians)
        sint = np.sin(theta_radians)
        xsigma2 = self.xsigma**2
        ysigma2 = self.ysigma**2

        self.a = (cost**2 / (2.0 * xsigma2)) + (sint**2 / (2.0 * ysigma2))
        # CCW
        self.b = 0.5 * cost * sint * ((1.0 / xsigma2) - (1.0 / ysigma2))
        self.c = (sint**2 / (2.0 * xsigma2)) + (cost**2 / (2.0 * ysigma2))

        # find the extent of an ellipse with radius = sigma_radius*sigma;
        # solve for the horizontal and vertical tangents of an ellipse
        # defined by g(x,y) = f
        self.f = self.sigma_radius**2 / 2.0
        denom = (self.a * self.c) - self.b**2

        # nx and ny are always odd
        self.nx = 2 * int(max(2, math.sqrt(self.c * self.f / denom))) + 1
        self.ny = 2 * int(max(2, math.sqrt(self.a * self.f / denom))) + 1

        self.xc = self.xradius = self.nx // 2
        self.yc = self.yradius = self.ny // 2

        # define the kernel on a 2D grid
        yy, xx = np.mgrid[0:self.ny, 0:self.nx]
        self.circular_radius = np.sqrt((xx - self.xc)**2 + (yy - self.yc)**2)
        self.elliptical_radius = (self.a * (xx - self.xc)**2
                                  + 2.0 * self.b * (xx - self.xc)
                                  * (yy - self.yc)
                                  + self.c * (yy - self.yc)**2)

        self.mask = np.where(
            (self.elliptical_radius <= self.f)
            | (self.circular_radius <= 2.0), 1, 0).astype(int)
        self.npixels = self.mask.sum()

        # NOTE: the central (peak) pixel of gaussian_kernel has a value of 1.
        self.gaussian_kernel_unmasked = np.exp(-self.elliptical_radius)
        self.gaussian_kernel = self.gaussian_kernel_unmasked * self.mask

        # denom = variance * npixels
        denom = ((self.gaussian_kernel**2).sum()
                 - (self.gaussian_kernel.sum()**2 / self.npixels))
        self.relerr = 1.0 / np.sqrt(denom)

        # normalize the kernel to zero sum
        if normalize_zerosum:
            self.data = ((self.gaussian_kernel
                          - (self.gaussian_kernel.sum() / self.npixels))
                         / denom) * self.mask
        else:
            self.data = self.gaussian_kernel

        self.shape = self.data.shape


class _StarCutout:
    """
    Class to hold a 2D image cutout of a single star for the star finder
    classes.

    Parameters
    ----------
    data : 2D array_like
        The cutout 2D image from the input unconvolved 2D image.

    convdata : 2D array_like
        The cutout 2D image from the convolved 2D image.

    slices : tuple of two slices
        A tuple of two slices representing the minimal box of the cutout
        from the original image.

    xpeak, ypeak : float
        The (x, y) pixel coordinates of the peak pixel.

    kernel : `_StarFinderKernel`
        The convolution kernel.  The shape of the kernel must match that
        of the input ``data``.

    threshold_eff : float
        The absolute image value above which to select sources.  This
        threshold should be the threshold value input to the star finder
        class multiplied by the kernel relerr.
    """

    def __init__(self, data, convdata, slices, xpeak, ypeak, kernel,
                 threshold_eff):

        self.data = data
        self.convdata = convdata
        self.slices = slices
        self.xpeak = xpeak
        self.ypeak = ypeak
        self.kernel = kernel
        self.threshold_eff = threshold_eff

        self.shape = data.shape
        self.nx = self.shape[1]  # always odd
        self.ny = self.shape[0]  # always odd
        self.cutout_xcenter = int(self.nx // 2)
        self.cutout_ycenter = int(self.ny // 2)

        self.xorigin = self.slices[1].start  # in original image
        self.yorigin = self.slices[0].start  # in original image

        self.mask = kernel.mask  # kernel mask
        self.npixels = kernel.npixels  # unmasked pixels
        self.data_masked = self.data * self.mask


def _find_stars(data, kernel, threshold_eff, min_separation=None,
                mask=None, exclude_border=False):
    """
    Find stars in an image.

    Parameters
    ----------
    data : 2D array_like
        The 2D array of the image.

    kernel : `_StarFinderKernel`
        The convolution kernel.

    threshold_eff : float
        The absolute image value above which to select sources.  This
        threshold should be the threshold input to the star finder class
        multiplied by the kernel relerr.

    min_separation : float, optional
        The minimum separation for detected objects in pixels.

    mask : 2D bool array, optional
        A boolean mask with the same shape as ``data``, where a `True`
        value indicates the corresponding element of ``data`` is masked.
        Masked pixels are ignored when searching for stars.

    exclude_border : bool, optional
        Set to `True` to exclude sources found within half the size of
        the convolution kernel from the image borders.  The default is
        `False`, which is the mode used by IRAF's `DAOFIND`_ and
        `starfind`_ tasks.

    Returns
    -------
    objects : list of `_StarCutout`
        A list of `_StarCutout` objects containing the image cutout for
        each source.

    .. _DAOFIND: https://iraf.net/irafhelp.php?val=daofind

    .. _starfind: https://iraf.net/irafhelp.php?val=starfind
    """
    convolved_data = _filter_data(data, kernel.data, mode='constant',
                                  fill_value=0.0, check_normalization=False)

    # define a local footprint for the peak finder
    if min_separation is None:  # daofind
        footprint = kernel.mask.astype(bool)
    else:
        # define a circular footprint
        idx = np.arange(-min_separation, min_separation + 1)
        xx, yy = np.meshgrid(idx, idx)
        footprint = np.array((xx**2 + yy**2) <= min_separation**2, dtype=int)

    # pad the data and convolved image by the kernel x/y radius to allow
    # for detections near the edges
    if not exclude_border:
        ypad = kernel.yradius
        xpad = kernel.xradius
        pad = ((ypad, ypad), (xpad, xpad))
        pad_mode = 'constant'
        data = np.pad(data, pad, mode=pad_mode, constant_values=[0.])
        if mask is not None:
            mask = np.pad(mask, pad, mode=pad_mode, constant_values=[0.])
        convolved_data = np.pad(convolved_data, pad, mode=pad_mode,
                                constant_values=[0.])

    # find local peaks in the convolved data
    with warnings.catch_warnings():
        # suppress any NoDetectionsWarning from find_peaks
        warnings.filterwarnings('ignore', category=NoDetectionsWarning)
        tbl = find_peaks(convolved_data, threshold_eff, footprint=footprint,
                         mask=mask)

    if tbl is None:
        return None

    coords = np.transpose([tbl['y_peak'], tbl['x_peak']])

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

        star_cutouts.append(_StarCutout(data_cutout, convdata_cutout, slices,
                                        xpeak, ypeak, kernel, threshold_eff))

    return star_cutouts


class _DAOFindProperties:
    """
    Class to calculate the properties of each detected star, as defined
    by `DAOFIND`_.

    Parameters
    ----------
    star_cutout : `_StarCutout`
        A `_StarCutout` object containing the image cutout for the star.

    kernel : `_StarFinderKernel`
        The convolution kernel.  The shape of the kernel must match that
        of the input ``star_cutout``.

    sky : float, optional
        The local sky level around the source.  ``sky`` is used only to
        calculate the source peak value, flux, and magnitude.  The
        default is 0.

    .. _DAOFIND: https://iraf.net/irafhelp.php?val=daofind
    """

    def __init__(self, star_cutout, kernel, sky=0.):
        if not isinstance(star_cutout, _StarCutout):
            raise ValueError('data must be an _StarCutout object')

        if star_cutout.data.shape != kernel.shape:
            raise ValueError('cutout and kernel must have the same shape')

        self.cutout = star_cutout
        self.kernel = kernel
        self.sky = sky  # DAOFIND has no sky input -> same as sky=0.

        self.data = star_cutout.data
        self.data_masked = star_cutout.data_masked
        self.npixels = star_cutout.npixels  # unmasked pixels
        self.nx = star_cutout.nx
        self.ny = star_cutout.ny
        self.xcenter = star_cutout.cutout_xcenter
        self.ycenter = star_cutout.cutout_ycenter

    @lazyproperty
    def data_peak(self):
        return self.data[self.ycenter, self.xcenter]

    @lazyproperty
    def conv_peak(self):
        return self.cutout.convdata[self.ycenter, self.xcenter]

    @lazyproperty
    def roundness1(self):
        # set the central (peak) pixel to zero
        cutout_conv = self.cutout.convdata.copy()
        cutout_conv[self.ycenter, self.xcenter] = 0.0  # for sum4

        # calculate the four roundness quadrants.
        # the cutout size always matches the kernel size, which have odd
        # dimensions.
        # quad1 = bottom right
        # quad2 = bottom left
        # quad3 = top left
        # quad4 = top right
        # 3 3 4 4 4
        # 3 3 4 4 4
        # 3 3 x 1 1
        # 2 2 2 1 1
        # 2 2 2 1 1
        quad1 = cutout_conv[0:self.ycenter + 1, self.xcenter + 1:]
        quad2 = cutout_conv[0:self.ycenter, 0:self.xcenter + 1]
        quad3 = cutout_conv[self.ycenter:, 0:self.xcenter]
        quad4 = cutout_conv[self.ycenter + 1:, self.xcenter:]

        sum2 = -quad1.sum() + quad2.sum() - quad3.sum() + quad4.sum()
        if sum2 == 0:
            return 0.

        sum4 = np.abs(cutout_conv).sum()
        if sum4 <= 0:
            return None

        return 2.0 * sum2 / sum4

    @lazyproperty
    def sharpness(self):
        npixels = self.npixels - 1  # exclude the peak pixel
        data_mean = (np.sum(self.data_masked) - self.data_peak) / npixels

        return (self.data_peak - data_mean) / self.conv_peak

    def daofind_marginal_fit(self, axis=0):
        """
        Fit 1D Gaussians, defined from the marginal x/y kernel
        distributions, to the marginal x/y distributions of the original
        (unconvolved) image.

        These fits are used calculate the star centroid and roundness
        ("GROUND") properties.

        Parameters
        ----------
        axis : {0, 1}, optional
            The axis for which the marginal fit is performed:

            * 0: for the x axis
            * 1: for the y axis

        Returns
        -------
        dx : float
            The fractional shift in x or y (depending on ``axis`` value)
            of the image centroid relative to the maximum pixel.

        hx : float
            The height of the best-fitting Gaussian to the marginal x or
            y (depending on ``axis`` value) distribution of the
            unconvolved source data.
        """

        # define triangular weighting functions along each axis, peaked
        # in the middle and equal to one at the edge
        x = self.xcenter - np.abs(np.arange(self.nx) - self.xcenter) + 1
        y = self.ycenter - np.abs(np.arange(self.ny) - self.ycenter) + 1
        xwt, ywt = np.meshgrid(x, y)

        if axis == 0:  # marginal distributions along x axis
            wt = xwt[0]  # 1D
            wts = ywt  # 2D
            size = self.nx
            center = self.xcenter
            sigma = self.kernel.xsigma
            dxx = center - np.arange(size)
        elif axis == 1:  # marginal distributions along y axis
            wt = np.transpose(ywt)[0]  # 1D
            wts = xwt  # 2D
            size = self.ny
            center = self.ycenter
            sigma = self.kernel.ysigma
            dxx = np.arange(size) - center

        # compute marginal sums for given axis
        wt_sum = np.sum(wt)
        dx = center - np.arange(size)

        # weighted marginal sums
        kern_sum_1d = np.sum(self.kernel.gaussian_kernel_unmasked * wts,
                             axis=axis)
        kern_sum = np.sum(kern_sum_1d * wt)
        kern2_sum = np.sum(kern_sum_1d**2 * wt)

        dkern_dx = kern_sum_1d * dx
        dkern_dx_sum = np.sum(dkern_dx * wt)
        dkern_dx2_sum = np.sum(dkern_dx**2 * wt)
        kern_dkern_dx_sum = np.sum(kern_sum_1d * dkern_dx * wt)

        data_sum_1d = np.sum(self.data * wts, axis=axis)
        data_sum = np.sum(data_sum_1d * wt)
        data_kern_sum = np.sum(data_sum_1d * kern_sum_1d * wt)
        data_dkern_dx_sum = np.sum(data_sum_1d * dkern_dx * wt)
        data_dx_sum = np.sum(data_sum_1d * dxx * wt)

        # perform linear least-squares fit (where data = sky + hx*kernel)
        # to find the amplitude (hx)
        # reject the star if the fit amplitude is not positive
        hx_numer = data_kern_sum - (data_sum * kern_sum) / wt_sum
        if hx_numer <= 0.:
            return np.nan, np.nan

        hx_denom = kern2_sum - (kern_sum**2 / wt_sum)
        if hx_denom <= 0.:
            return np.nan, np.nan

        # compute fit amplitude
        hx = hx_numer / hx_denom
        # sky = (data_sum - (hx * kern_sum)) / wt_sum

        # compute centroid shift
        dx = ((kern_dkern_dx_sum
               - (data_dkern_dx_sum - dkern_dx_sum * data_sum))
              / (hx * dkern_dx2_sum / sigma**2))

        hsize = size / 2.
        if abs(dx) > hsize:
            if data_sum == 0.:
                dx = 0.0
            else:
                dx = data_dx_sum / data_sum
                if abs(dx) > hsize:
                    dx = 0.0

        return dx, hx

    @lazyproperty
    def dx_hx(self):
        return self.daofind_marginal_fit(axis=0)

    @lazyproperty
    def dy_hy(self):
        return self.daofind_marginal_fit(axis=1)

    @lazyproperty
    def dx(self):
        return self.dx_hx[0]

    @lazyproperty
    def dy(self):
        return self.dy_hy[0]

    @lazyproperty
    def xcentroid(self):
        return self.cutout.xpeak + self.dx

    @lazyproperty
    def ycentroid(self):
        return self.cutout.ypeak + self.dy

    @lazyproperty
    def hx(self):
        return self.dx_hx[1]

    @lazyproperty
    def hy(self):
        return self.dy_hy[1]

    @lazyproperty
    def roundness2(self):
        """
        The star roundness.

        This roundness parameter represents the ratio of the difference
        in the height of the best fitting Gaussian function in x minus
        the best fitting Gaussian function in y, divided by the average
        of the best fitting Gaussian functions in x and y.  A circular
        source will have a zero roundness.  A source extended in x or y
        will have a negative or positive roundness, respectively.
        """

        if np.isnan(self.hx) or np.isnan(self.hy):
            return np.nan
        else:
            return 2.0 * (self.hx - self.hy) / (self.hx + self.hy)

    @lazyproperty
    def peak(self):
        return self.data_peak - self.sky

    @lazyproperty
    def npix(self):
        """
        The total number of pixels in the rectangular cutout image.
        """

        return self.data.size

    @lazyproperty
    def flux(self):
        return ((self.conv_peak / self.cutout.threshold_eff)
                - (self.sky * self.npix))

    @lazyproperty
    def mag(self):
        if self.flux <= 0:
            return np.nan
        else:
            return -2.5 * np.log10(self.flux)
