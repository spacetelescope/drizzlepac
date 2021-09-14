"""This script contains code to support creation of photometric sourcelists using two techniques:
aperture photometry and segmentation-map based photometry."""

import os
import sys
import shutil
import warnings

import numpy as np

import skimage
from astropy.io import fits as fits

from photutils.detection import StarFinderBase, find_peaks
from photutils.utils import NoDetectionsWarning

from stsci.tools import fileutil as fu

from . import astrometric_utils as amutils
from ._detection_utils import (_DAOFindProperties, _StarCutout,
                               _StarFinderKernel, _find_stars)
from .. import astrodrizzle

#
# Original imports
#

from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from scipy import ndimage
import scipy.signal as ss

from stsci.tools import logutil

try:
    from matplotlib import pyplot as plt
except Exception:
    plt = None

CATALOG_TYPES = ['aperture', 'segment']
INSTR_KWS = ['INSTRUME', 'DETECTOR']
FILTER_KW = "FILTER*"

PSF_PATH = ['pars', 'psfs']

__taskname__ = 'deconvolve_utils'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)


# ======================================================================================================================

def fft_deconv_img(img, psf, freq_limit=0.95, block_size=(1024, 1024)):
    """ FFT image deconvolution

    This function performs a simple 1-step FFT-based deconvolution of the
    input image using the specified PSF.  The input image get transformed
    by FFT section-by-section in blocks defined by the `block_size` parameter
    in order to minimize the memory use for what can be extremely large images.

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

    block_size : tuple, optional
        This specifies how much of the input image will be transformed using an FFT at one time.

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

    psf = _center_psf(psf, block_size)

    # Compute alpha scaling based on PSF frequency limit
    P2 = np.abs(np.copy(psf).flatten())
    P2 = np.sort(P2)[::-1]
    index = int(P2.shape[0] * freq_limit)
    alpha = P2[index]
    del P2

    # Insure input image also has no NaN values, only if necessary
    # This should be less memory-intensive than np.isnan()
    if len(np.where(img == np.nan)[0]) > 0:
        img = np.nan_to_num(img, copy=True, nan=0.0)
    img_shape = img.shape

    # Break up image into blocks that will be deconvolved separately, then
    # pieced back together again.
    # Returns: {'slices':[], 'new_shape':(y, x), 'blocks':[N, M, y, x]}
    block_dict = create_blocks(img, block_size=block_size)
    img_blocks = block_dict['blocks']
    del img

    # FFT point spread function (first column of theory matrix G)
    # In order to avoid memory errors, this code:
    #   - Splits the input into (non-overlapping) blocks using:
    #        slices = skimage.util.view_as_blocks(arr, block_shape=(1024,1024))
    #   -  perform the FFT on each block separately (as if separate exposures)
    #        spectrum = np.fft.fft2(slices)
    m_maps = img_blocks * 0.
    for a in range(img_blocks.shape[0]):
        for b in range(img_blocks.shape[1]):
            block = img_blocks[a, b, :, :]
            m_map = _perform_deconv(block, psf, alpha)
            m_maps[a, b, :, :] = m_map

    # Re-constitute deconvolved image blocks into single image
    # rebuild_arr(block_arr, slices, new_shape, output_shape)
    deconv_img = rebuild_arr(m_maps,
                             block_dict['slices'],
                             block_dict['new_shape'],
                             img_shape)
    return deconv_img

# ======================================================================================================================
# Functions to manage PSF library for deconvolution
#
def _perform_deconv(img_block, psf, alpha):
    P = np.fft.fft2(psf)

    # FFT2 measurement
    # Use image in husky_conv.png
    # U^H d
    D = np.fft.fft2(img_block)

    # -dampped spectral components,
    # -also known as Wiener filtering
    # (conj(S)/(|S|^2 + alpha^2)) U^H d
    M = (np.conj(P) / (np.abs(P)**2.0 + alpha**2.0)) * D

    # maximum a posteriori estimate of deconvolved image
    # m_map = U (conj(S)/(|S|^2 + alpha^2)) U^H d
    m_map = (D.shape[1] * D.shape[0]) * np.fft.fftshift(np.fft.ifft2(M).real)

    zero_mask = (img_block > 0).astype(np.int16)
    m_map *= zero_mask

    return m_map


def _center_psf(psf, img_block_shape):
    """Create centered PSF image with same shape as img_block"""
    psf_y, psf_x = np.where(psf == psf.max())
    psf_y = psf_y[0]
    psf_x = psf_x[0]
    psf_center = [(img_block_shape[0] // 2), (img_block_shape[1] // 2)]
    centered_psf = np.zeros(img_block_shape, dtype=psf.dtype)

    # We need to recenter the PSF
    psf_section = np.where(psf != 0.0)
    psf_xr = [psf_section[1].min(), psf_section[1].max()]
    psf_yr = [psf_section[0].min(), psf_section[0].max()]
    psf_size = [(psf_yr[1] - psf_yr[0]) // 2, (psf_xr[1] - psf_xr[0]) // 2]
    # If psf_size is not at least 2 * psf_max//2, increase to that value
    # This will insure that the final position of the PSF in the output is at the exact center
    psf_max = np.where(psf[psf_yr[0]:psf_yr[1], psf_xr[0]:psf_xr[1]] == psf.max())
    psf_len = max(max(psf_max[0][0] // 2, psf_size[0]), max(psf_max[1][0] // 2, psf_size[1]))

    centered_psf[psf_center[0] - psf_len: psf_center[0] + psf_len,
                 psf_center[1] - psf_len: psf_center[1] + psf_len] = psf[psf_y - psf_len: psf_y + psf_len,
                                                                         psf_x - psf_len: psf_x + psf_len]
    return centered_psf


def pad_arr(arr, block=(1024, 1024)):
    """ Zero-pad an input array up to an integer number of blocks in each dimension

    Parameters
    ----------
    arr : `numpy.ndarray`
        Original input array to be padded to the new size

    block : `tuple` or `int`
        Size of blocks which should be used to define the output size so that the
        output image is an integer number of blocks with this size.  If only an
        integer is specified, then a block size of (block, block) will be used.

    Returns
    -------
    new_arr : `numpy.ndarray`
        Resized output array of size (n*block[0], m*block[1]).

    """
    if isinstance(block, int):
        block = (block, block)
    new_shape = (arr.shape[0] + (block[0] - (arr.shape[0] % block[0])),
                 arr.shape[1] + (block[1] - (arr.shape[1] % block[1])))
    new_arr = np.zeros(new_shape, dtype=arr.dtype)
    new_arr[:arr.shape[0], :arr.shape[1]] = arr
    return new_arr


def create_blocks(arr, block_size=(1024, 1024)):
    """Split input array into uniformly-sized blocks

    This function will split the input array into uniformly-sized blocks
    with the size of each block specified as `block_size`.  The input array
    will be zero-padded in either or both axes in order to expand the array
    to an integer number of blocks in both dimensions.

    Parameters
    ----------
    arr : `numpy.ndarray`
        2-D Input image of size N x M

    block_size : `tuple`, optional
        Tuple specifying the block size n x m

    Returns
    --------
    block_dict : `dict`
        This dictionary contains all the information describing all the blocks
        created from the input array consisting of
        {'slices':[], 'new_shape': (N`, M`), `blocks`: []} where:

        ``"slices"``:
            List of indices [y, x, n, m] describing each block, where,
            (y,x) is index of the block in the original array and
            (n,m) is the number of pixels in the block.
        ``"new_shape"``:
            Full size of the padded array that was cut into blocks
        ``"blocks"``:
            Actual blocks as returned by `skimage.util.view_as_blocks`.

    """
    if arr.shape[0] % block_size[0] != 0 or arr.shape[1] % block_size[1] != 0:
        new_arr = pad_arr(arr, block=block_size)
    else:
        new_arr = arr

    new_shape = new_arr.shape

    # Create blocks from image of size block_size
    # Output set of blocks will have shape: [N, M, block_size[0], block_size[1]]
    blocks = skimage.util.view_as_blocks(new_arr, block_size)

    slices = []
    for a in range(blocks.shape[0]):
        for b in range(blocks.shape[1]):
            slices.append([a, b, slice(block_size[0] * a, block_size[0] * (a + 1)),
                                 slice(block_size[1] * b, block_size[1] * (b + 1))])

    del new_arr
    return {'slices': slices, 'new_shape': new_shape, 'blocks': blocks}


def rebuild_arr(block_arr, slices, new_shape, output_shape):
    """Convert convolved blocks into full array that matches original input array

    Parameters
    -----------
    block_arr : `numpy.ndarray`
        List of image blocks which need to be pieced back together into a single
        array.

    slices : `list`
        List of `slices` from `create_blocks` specifying the location ajd size of each
        block from the original input array.

    new_shape : `tuple`
        Full size of the padded image used to create the uniformly sized blocks.

    output_shape : `tuple`
        The original shape of the array before padding and creating the blocks.

    Returns
    --------
    out_arr : `numpy.ndarray`
        Single array of same size as original input array before padding and splitting into blocks.

    """
    out_arr = np.zeros(new_shape, dtype=block_arr.dtype)
    for s in slices:
        out_arr[(s[2], s[3])] = block_arr[s[0], s[1], :, :]
    return out_arr[:output_shape[0], :output_shape[1]]


def find_psf(imgname, path_root=None):
    """Pull PSF from library based on unique combination of intrument/detector/filter.

    Parameters
    ===========
    imgname : str
        Image name of science image to be deconvolved.  If the header contains
        instrument, detector, and filters, then those additional parameters do not need to be specified.

    path_root : str, optional
        Full path to parent directory of PSF library IF not the default path for package.

    Returns
    ========
    psfnames : list
        List of the Full filenames, with path, for all PSFs from the library
        that apply to the input image.

    """
    # We need to create an average PSF made up of all the filters
    # used to create the image.
    # Start by looking in the input image header for
    # the values for the instrument, detector and filter* keywords,
    # so we know what PSF(s) to look for.
    total_hdu = fits.open(imgname)
    # get list of all input exposures
    input_files = [f.split('[')[0] for f in total_hdu[0].header.get('d*data').values()]
    # If there were (for any reason) no D???DATA keywords,
    #  only use the input image to define the filters for the PSF
    if len(input_files) == 0:
        input_files = [imgname]
    total_hdu.close()
    del total_hdu

    # Reduce the list down to only unique filenames
    input_files = list(dict.fromkeys(input_files))
    # get filter names from each input file
    filter_list = [get_filter_names(f) for f in input_files]
    kw_vals = [filter_list[0][0].lower(), filter_list[0][1].lower()]

    # Set up path to PSF library installed with code based on tree:
    #     drizzlepac/pars/psfs/<instrument>/<detector>
    if path_root is None:
        path_root = os.path.split(os.path.dirname(__file__))[0]
        for psf_path in PSF_PATH:
            path_root = os.path.join(path_root, psf_path)
    path_root = os.path.join(path_root, kw_vals[0], kw_vals[1])

    # Now look for filename associated with selected filter
    psf_names = [os.path.join(path_root, "{}.fits".format("_".join(kw_vals))) for kw_vals in filter_list]
    # Again, remove duplicate entries
    psf_names = list(dict.fromkeys(psf_names))

    log.debug('Looking for Library PSFs:\n  {}'.format(psf_names))

    psfs_exist = [os.path.exists(fname) for fname in psf_names]
    if not all(psfs_exist):
        log.error('Some PSF NOT found for keywords {} \n   with values of {}'.format(INSTR_KWS, filter_list))
        if not any(psfs_exist):
            log.error('NO PSF(s) found for keywords {} \n   with values of {}'.format(INSTR_KWS, filter_list))
            raise ValueError

    log.info("Using Library PSF(s):\n    {}".format([os.path.basename(name) for name in psf_names]))

    return psf_names


def get_filter_names(imgname):
    """Interpret photmode from image

    Parameters
    -----------
    imgname : str
        Filename of observation to extract filter names from

    Returns
    --------
    kw_vals : list
        List containing instrument, detector, and filter names in that order.
    """

    # look for instrument, detector, and filter in input image name
    # Start by looking in the input image header for these values, so we know what to look
    hdu = fits.open(imgname)
    photmode = hdu[('sci', 1)].header.get('photmode')
    if photmode is None:
        photmode = hdu[0].header.get('photmode')

    detector = hdu[0].header.get('detector')
    hdu.close()
    del hdu

    # Remove duplicate entries from photmode, if any (like 2 CLEAR filters)
    filter_list = list(dict.fromkeys(photmode.split(' ')))
    # Key off of instrument and detector (first 2 members of photmode)
    kw_vals = [filter_list[0].lower(), detector.lower()]

    # This will remove all polarizer, grism and prism filter entries as well.
    excl_filters = ['CAL', 'MJD', 'POL', 'GR', 'PR']
    for e in excl_filters:
        indx = [filter_list.index(i) for i in filter_list if i.startswith(e)]
        indx.reverse()
        for i in indx:
            del filter_list[i]

    # Now select the widest-band filter used for this observation
    # This accounts for cross-filter usage.
    found_filter = False
    if len(filter_list[2:]) >= 1:
        # Loop over types of bandpass based on filter names
        # going from widest band to narrowest bands
        bandpass = ['lp', 'w', 'm', 'n']
        for bp in bandpass:
            for f in filter_list[2:]:
                if f.lower().endswith(bp):
                    kw_vals += [f.lower()]
                    found_filter = True
                    break

    if not found_filter:
        kw_vals += ['clear']
    # The result at this point will be:
    #  kw_vals = [instrument, detector, selected filter]
    return kw_vals

def convert_library_psf(calimg, drzimg, psfs,
                        total_flux=100000.0,
                        pixfrac=1.0,
                        clean_psfs=True):
    """Drizzle library PSFs to match science image. """

    psf_flt_names = [_create_input_psf(psfname, calimg, total_flux) for psfname in psfs]

    # Insure only final drizzle step is run
    drizzle_pars = {}
    drizzle_pars["build"] = True
    drizzle_pars['context'] = False
    drizzle_pars['preserve'] = False
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
    drizzle_pars['final_fillval'] = 0.0
    drizzle_pars['final_pixfrac'] = pixfrac
    drizzle_pars["final_refimage"] = "{}[1]".format(drzimg)

    psf_drz_name = psf_flt_names[0].replace('_flt.fits', '')  # astrodrizzle will add suffix
    psf_drz_output = "{}_drz.fits".format(psf_drz_name)

    # Drizzle PSF FLT file to match orientation and plate scale of drizzled science (total detection) image
    astrodrizzle.AstroDrizzle(input=psf_flt_names,
                              output=psf_drz_name,
                              **drizzle_pars)

    if clean_psfs:
        # clean up intermediate files
        for psf_name in psf_flt_names:
            os.remove(psf_name)

    return psf_drz_output


def _create_input_psf(psf_name, calimg, total_flux):

    # Create copy of input science image based on input psf filename
    psf_root = os.path.basename(psf_name)
    lib_psf_arr = fits.getdata(psf_name)
    lib_psf_arr *= total_flux

    lib_size = [lib_psf_arr.shape[0] // 2, lib_psf_arr.shape[1] // 2]

    # create hamming 2d filter to avoid edge effects
    h = ss.hamming(lib_psf_arr.shape[0])
    h2d = np.sqrt(np.outer(h, h))
    lib_psf_arr *= h2d

    # This will be the name of the new file containing the library PSF that will be drizzled to
    # match the input image `drzimg`
    psf_flt_name = psf_root.replace('.fits', '_psf_flt.fits')

    # create version of PSF that will be drizzled
    psf_base = fits.getdata(calimg, ext=1) * 0.0
    # Copy library PSF into this array
    out_cen = [psf_base.shape[0] // 2, psf_base.shape[1] // 2]
    edge = (lib_psf_arr.shape[0] % 2, lib_psf_arr.shape[1] % 2)
    psf_base[out_cen[0] - lib_size[0]: out_cen[0] + lib_size[0] + edge[0],
             out_cen[1] - lib_size[1]: out_cen[1] + lib_size[1] + edge[1]] = lib_psf_arr

    # Write out library PSF FLT file now
    psf_flt = shutil.copy(calimg, psf_flt_name)

    # Update file with library PSF
    flt_hdu = fits.open(psf_flt, mode='update')
    flt_hdu[('sci', 1)].data = psf_base
    flt_hdu[('sci', 1)].header['psf_nx'] = psf_base.shape[1]
    flt_hdu[('sci', 1)].header['psf_ny'] = psf_base.shape[0]
    num_sci = fu.countExtn(calimg)
    # Also zero out all other science data in this 'PSF' file.
    if num_sci > 1:
        for extn in range(2, num_sci + 1):
            flt_hdu[('sci', extn)].data *= 0.0
    flt_hdu.close()
    del flt_hdu, lib_psf_arr

    return psf_flt_name

def get_cutouts(data, star_list, kernel, threshold_eff, exclude_border=False):

    coords = [(row[1], row[0]) for row in star_list]
    convolved_data = data

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
        # Skip slices which include pixels with a value of NaN
        if np.isnan(data_cutout).any():
            continue
        convdata_cutout = convolved_data[slices]

        # correct pixel values for the previous image padding
        if exclude_border:
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


class UserStarFinder(StarFinderBase):
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

    coords : `~astropy.table.Table` or `None`
        A table, such as returned by `find_peaks`, with approximate X,Y positions
        of identified sources.
        If not provided, the DAOFind algorithm will be used to find sources.

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
                 coords=None,
                 brightest=None, peakmax=None):

        if not np.isscalar(threshold):
            raise TypeError('threshold must be a scalar value.')
        self.threshold = threshold

        if not np.isscalar(fwhm):
            raise TypeError('fwhm must be a scalar value.')
        self.fwhm = fwhm

        self.coords = coords
        self.ratio = ratio
        self.theta = theta
        self.sigma_radius = sigma_radius
        self.sharplo = sharplo
        self.sharphi = sharphi
        self.roundlo = roundlo
        self.roundhi = roundhi
        self.sky = sky
        self.exclude_border = exclude_border

        self.kernel = _StarFinderKernel(self.fwhm, self.ratio, self.theta,
                                        self.sigma_radius)
        self.threshold_eff = self.threshold * self.kernel.relerr
        self.brightest = brightest
        self.peakmax = peakmax
        self._star_cutouts = None

    def find_stars(self, data, mask=None):
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
        if self.coords:
            star_cutouts = get_cutouts(data, self.coords,
                                        self.kernel,
                                        self.threshold_eff,
                                        exclude_border=self.exclude_border)
        else:
            star_cutouts = _find_stars(data, self.kernel, self.threshold_eff,
                                       mask=mask,
                                       exclude_border=self.exclude_border)

        if star_cutouts is None:
            warnings.warn('No sources were found.', NoDetectionsWarning)
            return None

        self._star_cutouts = star_cutouts

        star_props = []
        for star_cutout in star_cutouts:
            props = _DAOFindProperties(star_cutout, self.kernel, self.sky)

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
                          'and roundness criteria.', NoDetectionsWarning)
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


# -----------------------------------------------------------------------------
#
# Main user interface
#
# -----------------------------------------------------------------------------
def find_point_sources(drzname, data=None, mask=None,
                       def_fwhm=2.0,
                       box_size=11, block_size=(1024, 1024),
                       diagnostic_mode=False):
    """ Identify point sources most similar to TinyTim PSFs

    Primary user-interface to identifying point-sources in the
    drizzle product image most similar to the TinyTim PSF for the
    filter-combination closest to that found in the drizzled image.
    The PSFs are pulled, by default, from those installed with the
    code as created using the TinyTim PSF modelling software for
    every direct image filter used by the ACS and WFC3 cameras on HST.

    .. note: Sources identified by this function will only have integer pixel
    positions.

    Parameters
    -----------
    drzname : `str`
        Filename of the drizzled image which should be used to find
        point sources.  This will provide the information on the filters
        used on the all the input exposures.

    data : `numpy.ndarray`, optional
        If provided, will be used as the image to be evaluated instead
        of opening the file specified in `drzname`.

    mask : `numpy.ndarray`, optional
        If provided, this mask will be used to eliminate regions in the
        input array from being searched for point sources.  Pixels with
        a value of 0 in the mask indicate what pixels should be ignored.

    def_fwhm : `float`, optional
        Default FWHM to use in case the model PSF can not be accurately
        measured by `photutils`.

    box_size : `int`, optional
        Size of the box used to recognize each point source.

    block_size : `tuple`, optional
        (Y, X) size of the block used by the FFT to process the drizzled image.

    diagnostic_mode : `bool`, optional
        Specify whether or not to provide additional diagnostic messages
        and output while processing.

    Returns
    -------
    peaks : `astropy.table.Table`
        Output from `photutils.detection.find_peaks` for all identified sources
        with columns `x_peak`, `y_peak` and `peak_value`.

    psf_fwhm : `float`
        FWHM (in pixels) of PSF used to identify the sources.

    """
    # determine the name of at least 1 input exposure
    calname = determine_input_image(drzname)
    sep = box_size // 2

    if not isinstance(block_size, tuple):
        block_size = tuple(block_size)

    if data is None:
        # load image
        drzhdu = fits.open(drzname)

        sciext = 0 if len(drzhdu) == 1 else ('sci', 1)
        drz = drzhdu[sciext].data.copy()
        drzhdr = drzhdu[sciext].header.copy()
        drzhdu.close()
        del drzhdu

        if mask is not None:
            # Apply any user-specified mask
            drz *= mask
    else:
        drz = data
        drzhdr = None

    if mask is not None:
        # invert the mask
        invmask = np.invert(mask)
    else:
        invmask = None

    # Identify PSF for image
    psfnames = find_psf(drzname)
    # Load PSF and convert to be consistent (orientation) with image
    clean_psfs = True if not diagnostic_mode else False

    drzpsfname = convert_library_psf(calname, drzname, psfnames,
                                     pixfrac=1.5,
                                     clean_psfs=clean_psfs)
    drzpsf = fits.getdata(drzpsfname)
    # try to measure just the core of the PSF
    # This will be a lot less likely to result in invalid/impossible FWHM values
    max_y, max_x = np.where(drzpsf == drzpsf.max())
    xc = max_x[0]
    yc = max_y[0]
    psf_core = drzpsf[yc - box_size: yc + box_size, xc - box_size: xc + box_size]
    psf_fwhm = amutils.find_fwhm(psf_core, def_fwhm)

    # check value
    if psf_fwhm < 0 or psf_fwhm > 2.0 * def_fwhm:
        # Try a different starting guess for the FWHM
        psf_fwhm = amutils.find_fwhm(psf_core, def_fwhm + 1)

        if psf_fwhm < 0 or psf_fwhm > 2.0 * def_fwhm:
            log.debug("FWHM computed as {}.  Reverting to using default FWHM of {}".format(psf_fwhm, def_fwhm))
            psf_fwhm = def_fwhm

    log.info("Library PSF FWHM computed as {}.".format(psf_fwhm))

    # deconvolve the image with the PSF
    decdrz = fft_deconv_img(drz, drzpsf,
                            block_size=block_size)

    if mask is not None:
        decmask = ndimage.binary_erosion(mask, iterations=box_size)
        decdrz *= decmask

    if diagnostic_mode:
        fits.PrimaryHDU(data=decdrz,
                        header=drzhdr).writeto(drzname.replace('.fits', '_deconv.fits'),
                                               overwrite=True)
        if mask is not None:
            fits.PrimaryHDU(data=decmask.astype(np.uint16)).writeto(drzname.replace('.fits', '_deconv_mask.fits'),
                                                                overwrite=True)
    # find sources in deconvolved image
    dec_peaks = find_peaks(decdrz, threshold=0.0,
                       mask=invmask, box_size=box_size)

    # Use these positions as an initial guess for the final position
    peak_mask = (drz * 0.).astype(np.uint8)
    # Do this by creating a mask for the original input that only
    # includes those pixels with 2 pixels of each peak from the
    # deconvolved image.
    for peak in dec_peaks:
        x = peak['x_peak']
        y = peak['y_peak']
        peak_mask[y - sep: y + sep + 1, x - sep: x + sep + 1] = 1
    drz *= peak_mask
    if diagnostic_mode:
        fits.PrimaryHDU(data=drz).writeto(drzname.replace('.fits', '_peak_mask.fits'), overwrite=True)

    # Use this new mask to find the actual peaks in the original input
    # but only to integer pixel precision.
    peaks = find_peaks(drz, threshold=0., box_size=box_size // 2)
    if len(peaks) == 0:
        peaks = None

    # Remove PSF used, unless running in diagnostic_mode
    if not diagnostic_mode:
        if os.path.exists(drzpsfname):
            os.remove(drzpsfname)
    del peak_mask

    return peaks, psf_fwhm

def determine_input_image(image):
    """Determine the name of an input exposure for the given drizzle product"""
    calimg = None
    with fits.open(image) as hdu:
        calimg = hdu[0].header['d001data']
    if calimg:
        calimg = calimg.split('[')[0]
    else:
        log.warn('No input image found in "D001DATA" keyword for {}'.format(image))

    return calimg
