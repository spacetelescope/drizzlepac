"""This script contains code to support creation of photometric sourcelists using two techniques:
aperture photometry and segmentation-map based photometry."""

import copy
import pickle  # FIX Remove
import sys

import astropy.units as u
from astropy.io import fits as fits
from astropy.stats import sigma_clipped_stats, gaussian_fwhm_to_sigma
from astropy.table import Column, MaskedColumn, Table, join, vstack
from astropy.convolution import RickerWavelet2DKernel
from astropy.coordinates import SkyCoord
import numpy as np
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

__taskname__ = 'deconvolve_utils'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)


# ======================================================================================================================

def fft_deconv_img(img, psf, alpha=500.0):
    #
    # Image deconvolution using FFT
    # (c) 2017 Juha Vierinen
    # http://www.radio-science.net/2017/09/deconvolution-in-frequency-domain-with.html
    #

    # regularization parameter
    # (should be one to two orders of magnitude below the largest spectral component of point-spread function)
    # alpha = 500.0

    # FFT point spread function (first column of theory matrix G)
    P = numpy.fft.fft2(psf)

    # FFT2 measurement
    # Use image in husky_conv.png
    # U^H d
    D = numpy.fft.fft2(img)

    # -dampped spectral components,
    # -also known as Wiener filtering
    # (conj(S)/(|S|^2 + alpha^2)) U^H d
    M = (numpy.conj(P) / (numpy.abs(P)**2.0 + alpha**2.0)) * D

    # maximum a posteriori estimate of deconvolved image
    # m_map = U (conj(S)/(|S|^2 + alpha^2)) U^H d
    m_map = (D.shape[1] * D.shape[0]) * numpy.fft.fftshift(numpy.fft.ifft2(M).real)

    return m_map
