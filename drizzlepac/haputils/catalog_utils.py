"""This script contains code to support creation of photometric sourcelists using two techniques:
aperture photometry and segmentation-map based photometry."""
import copy
import math
import sys
from collections import OrderedDict
from packaging.version import Version

from astropy.io import fits as fits
from astropy.stats import sigma_clipped_stats, gaussian_fwhm_to_sigma
from astropy.table import Column, MaskedColumn, Table, join, vstack
from astropy.convolution import RickerWavelet2DKernel, convolve
from astropy.coordinates import SkyCoord
import numpy as np
from scipy import ndimage, stats

from photutils.segmentation import (detect_sources, SourceCatalog,
                                    deblend_sources)
from photutils.aperture import CircularAperture, CircularAnnulus
from photutils.background import (Background2D, SExtractorBackground,
                                  StdBackgroundRMS)
from photutils.detection import DAOStarFinder, IRAFStarFinder
from photutils.utils import calc_total_error

from stsci.tools import logutil
from stwcs.wcsutil import HSTWCS

from . import astrometric_utils
from . import constants
from . import photometry_tools
from . import deconvolve_utils as decutils
from . import processing_utils as proc_utils

try:
    from matplotlib import pyplot as plt
except Exception:
    plt = None

CATALOG_TYPES = ['aperture', 'segment']

# B1950, J2000 and other valid RADESYS values do not transform natively to ICRS
# Only these options support conversion to/from ICRS in WCSLIB.
RADESYS_OPTIONS = ['FK4', 'FK5', 'ICRS']

id_colname = 'label'
flux_colname = 'segment_flux'
ferr_colname = 'segment_fluxerr'
bac_colname = 'background_centroid'

__taskname__ = 'catalog_utils'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)


class CatalogImage:
    def __init__(self, filename, num_images_mask, log_level):
        # set logging level to user-specified level
        log.setLevel(log_level)

        if isinstance(filename, str):
            self.imghdu = fits.open(filename)
            self.imgname = filename
        else:
            self.imghdu = filename
            self.imgname = filename.filename()

        # This is the "footprint_mask" of the total product object which indicates
        # the number of images which comprise each individual pixel
        self.num_images_mask = num_images_mask

        # Get header information to annotate the output catalogs
        if "total" in self.imgname:
            self.ghd_product = "tdp"
        else:
            self.ghd_product = "fdp"

        # Fits file read
        self.data = self.imghdu[('SCI', 1)].data
        self.wht_image = self.imghdu['WHT'].data.copy()

        # Add blank/empty flag
        self.blank = False
        data_vals = np.nan_to_num(self.data, 0)
        if data_vals.min() == 0 and data_vals.max() == 0:
            self.blank = True
        del data_vals

        # Get the HSTWCS object from the first extension
        self.imgwcs = HSTWCS(self.imghdu, 1)

        self.keyword_dict = self._get_header_data()

        # Populated by self.compute_background()
        self.bkg_background_ra = None
        self.bkg_rms_ra = None
        self.bkg_rms_median = None
        self.bkg_median = None
        self.footprint_mask = None
        self.inv_footprint_mask = None
        self.bkg_type = ""

        # Populated by self.build_kernel()
        self.kernel = None
        self.kernel_fwhm = None
        self.kernel_psf = False

    def close(self):
        self.imghdu.close()
        self.bkg_background_ra = None
        self.bkg_rms_ra = None
        self.bkg_rms_median = None
        # Finished with wht_image, clean up memory immediately...
        del self.wht_image
        self.wht_image = None

    def build_kernel(self, box_size, win_size, fwhmpsf,
                     simple_bkg=False,
                     bkg_skew_threshold=0.5,
                     zero_percent=25.0,
                     negative_percent=15.0,
                     nsigma_clip=3.0,
                     maxiters=3,
                     good_fwhm=[1.5, 3.5]):

        if self.blank:
            return
        if self.bkg_background_ra is None:
            self.compute_background(box_size, win_size,
                                    simple_bkg=simple_bkg,
                                    bkg_skew_threshold=bkg_skew_threshold,
                                    zero_percent=zero_percent,
                                    negative_percent=negative_percent,
                                    nsigma_clip=nsigma_clip,
                                    maxiters=maxiters)

        # Build a Gaussian2DKernel - the use of a custom kernel based upon a PSF
        # of the input image is no longer used here as it could be asymmetric and
        # lead to smeared segments (i.e., centroid coordinates offset from actual source).
        # The value of 11 for the Gaussian kernel is big enough for all detectors. This
        # value is hard-coded here in contrast to being in the configuration files like 
        # the RickerWavelet2DKernel.
        self.kernel_fwhm = fwhmpsf / self.imgwcs.pscale
        self.kernel = astrometric_utils.build_gaussian_kernel(self.kernel_fwhm, npixels=11)

    def compute_background(self, box_size, win_size,
                           bkg_estimator=SExtractorBackground, rms_estimator=StdBackgroundRMS,
                           simple_bkg=False,
                           bkg_skew_threshold=0.5,
                           zero_percent=25.0,
                           negative_percent=15.0,
                           nsigma_clip=3.0,
                           maxiters=3):
        """Use a sigma-clipped algorithm or Background2D to determine the background of the input image.

        There is also a special case background assessment reserved for SBC (MAMA detector) images.

        Parameters
        ----------
        image : ndarray
            Numpy array of the science extension from the observations FITS file.

        box_size : int
            Size of box along each axis

        win_size : int
            Size of 2D filter to apply to the background image

        bkg_estimator : subroutine
            background estimation algorithm

        rms_estimator : subroutine
            RMS estimation algorithm

        simple_bkg : bool, optional
            Forces use of the sigma_clipped_stats algorithm

        bkg_skew_threshold : float, optional
            Discriminator on the skewness computation - below this limit the Background2D algorithm
            will be computed for potential use for the background determination, otherwise
            the sigma_clipped_stats algorithm is used.

        zero_percent : float, optional
            Discriminator on the input image.  The percentage of zero values in the illuminated portion
            of the input image is determined - if there are more zero values than this lower limit, then
            the background is set to an image of constant value zero and the background rms is computed
            based on the pixels which are non-zero in the illuminated portion of the input image.

        negative_percent : float, optional
            Discriminator on the background-subtracted image.  The percentage of negative values in the
            background-subtracted image is determined - below this limit the Background2D algorithm stays in play,
            otherwise the sigma_clipped_stats algorithm is used.

        nsigma_clip : float, optional
            Parameter for the sigma_clipped_stats algorithm - number of standard deviations to use for both
            the lower and upper clipping limit.

        maxiters : float, optional
            Parameter for the sigma_clipped_stats algorithm - number of sigma-clipping iterations to perform

        Attributes
        ----------
        self.bkg_background_ra : 2D ndarray
            Background array

        self.bkg_rms_ra : 2D ndarray
            RMS map array

        self.bkg_median : float
            background median value over entire 2D array

        self.bkg_rms_median : float
            background rms value over entire 2D array

        self.footprint_mask :  bool 2Dndarry
            Footprint of input image set to True for the illuminated portion and False for
            the non-illuminated portion

        self.inv_footprint_mask :  bool 2Dndarry
            Inverse of the footprint_mask

        """
        if self.blank:
            self.bkg_type = 'empty'
            log.info("")
            log.info("*** Background image EMPTY. ***")
            return

        # Negative allowance in sigma
        negative_sigma = -1.0

        # Report configuration values to log
        log.info("")
        log.info("Background Computation")
        log.info("File: {}".format(self.imgname))
        log.info("Zero threshold: {}".format(zero_percent))
        log.info("Sigma-clipped Background Configuration Variables")
        log.info("  Negative percent threshold: {}".format(negative_percent))
        log.info("  Negative sigma: {}".format(negative_sigma))
        log.info("  Nsigma: {}".format(nsigma_clip))
        log.info("  Number of iterations: {}".format(maxiters))
        log.info("Background2D Configuration Variables")
        log.info("  Box size: {}".format(box_size))
        log.info("  Window size: {}".format(win_size))
        log.info("Background discriminant - skew threshold: {}".format(bkg_skew_threshold))

        # SExtractorBackground and StdBackgroundRMS are the defaults
        bkg = None
        is_zero_background_defined = False

        # Make a local copy of the data(image) being processed
        imgdata = copy.deepcopy(self.data)

        # In order to compute the proper statistics on the input data, need to use the footprint
        # mask to get the actual data - illuminated portion (True), non-illuminated (False) based
        # on images which actually contributed to a pixel, even if the contribution was the
        # value of zero. Again, the num_images_mask here is not True/False, but the *number of images*
        # which contributed to the value of a pixel.
        #
        # Note that if there are large shifts between images, the pixels in the gaps between the
        # images have values of 0.  The self.num_images_mask is the mask created by 
        # generate_footprint_mask() in product.py.
        footprint_mask = self.num_images_mask > 0
        self.footprint_mask = ndimage.binary_erosion(footprint_mask, iterations=10)
        self.inv_footprint_mask = np.invert(self.footprint_mask)

        # Get the statistics on the number of illuminated pixels
        num_of_illuminated_pixels = self.footprint_mask.sum()

        # BACKGROUND COMPUTATION 1 (unusual case which should only apply to ACS/SBC data)
        # Get the number of pixels which are zero in the actual footprint mask.  In
        # essence ignore the gaps between images which may also be set to zero and only do
        # this for SBC.
        if (self.imghdu[0].header['DETECTOR'].upper() == "SBC"):
            num_of_zeros = np.count_nonzero(imgdata[self.footprint_mask] == 0)

            # If there are too many background zeros in the image
            # (> number_of_zeros_in_background_threshold), set the background median to
            # zero and the background rms to the real rms of the non-zero values in the image.
            if num_of_zeros / float(num_of_illuminated_pixels) * 100.0 > zero_percent:
                self.bkg_median = 0.0
                self.bkg_rms_median = stats.tstd(imgdata[self.footprint_mask], limits=[0, None], inclusive=[False, True])
                self.bkg_background_ra = np.full_like(imgdata, 0.0)
                self.bkg_rms_ra = np.full_like(imgdata, self.bkg_rms_median)
                self.bkg_type = 'zero_background'

                is_zero_background_defined = True
                log.info(f"Input image contains excessive zero values in the background. Median: {self.bkg_median:.6f} RMS: {self.bkg_rms_median:.6f}")

                # Compute a minimum rms value based upon information directly from the data
                minimum_rms = 1.0 / self.keyword_dict['texpo_time']
                if np.nan_to_num(self.bkg_rms_median) < minimum_rms:
                    self.bkg_rms_median = minimum_rms
                    log.info("")
                    log.info(f"Minimum RMS of input based upon the total exposure time: {minimum_rms:.6f}")
                    log.info(f"Median RMS has been updated - Median: {self.bkg_median:.6f} RMS: {self.bkg_rms_median:.6f}")
                    log.info("")

        # BACKGROUND COMPUTATION 2 (sigma_clipped_stats)
        # If the input data is not the unusual case of SBC "excessive zero background", compute
        # a sigma-clipped background which returns only single values for mean,
        # median, and standard deviation (not 2D images)
        #
        # Computation 1 for rms
        if not is_zero_background_defined:
            log.info("")
            log.info("Computing the background using sigma-clipped statistics algorithm.")
            bkg_mean_full, bkg_median_full, bkg_rms_full = sigma_clipped_stats(imgdata,
                                                                self.inv_footprint_mask,
                                                                sigma=nsigma_clip,
                                                                cenfunc='median',
                                                                maxiters=maxiters)

            # guard against median being negative (can happen for mostly nebulous fields)
            if bkg_median_full < 0.0:
                # Recompute after adjusting input image data so that entire image is positive
                # This corrects for any gross over-subtraction of the background from the image
                imgdata -= (bkg_median_full - bkg_rms_full)
                bkg_mean_full, bkg_median_full, bkg_rms_full = sigma_clipped_stats(imgdata,
                                                                    self.inv_footprint_mask,
                                                                    sigma=nsigma_clip,
                                                                    cenfunc='median',
                                                                    maxiters=maxiters)

            # Compute Pearson’s second coefficient of skewness - this is a criterion
            # for possibly computing a two-dimensional background fit
            # Use the "raw" values generated by sigma_clipped_stats()
            bkg_skew = np.abs(3.0 * (bkg_mean_full - bkg_median_full) / bkg_rms_full)

            # Refine background to better compute the median value
            imgnz = imgdata * self.footprint_mask
            imgnz = imgnz[imgnz > 0.0]  # only want non-negative values in the illuminated region

            # Remove outlier values which are greater than the previously computed
            # (median + bit_of_rms)
            imgvals = imgnz[imgnz < (bkg_median_full + (bkg_rms_full * 0.1))]
            bkg_mean, bkg_median, bkg_rms = sigma_clipped_stats(imgvals,
                                                                None,
                                                                sigma=nsigma_clip,
                                                                cenfunc='median')

            log.info("Computation 1 of RMS - sigma clipped statistics")
            log.info(f"    Sigma-clipped Statistics - Background mean: {bkg_mean:.6f}  median: {bkg_median:.6f}  rms: {bkg_rms:.6f}")

            # Ensure the computed values are not negative
            if bkg_mean < 0.0 or bkg_median < 0.0:
                bkg_mean = max(0, bkg_mean)
                bkg_median = max(0, bkg_median)
                log.info("Updated Sigma-clipped Statistics for non-negative values.")
                log.info(f"    Sigma-clipped Statistics - Background mean: {bkg_mean:.6f}  median: {bkg_median:.6f}  rms: {bkg_rms:.6f}")

            # Computation 2 for rms: Compute a minimum rms value based upon the data
            if self.keyword_dict['detector'].upper() != 'SBC':
                minimum_rms = np.sqrt(self.keyword_dict['numexp']) * self.keyword_dict['readnse'] / self.keyword_dict['texpo_time']
            else:
                minimum_rms = np.sqrt(np.clip(bkg_median * self.keyword_dict['texpo_time'], a_min=1.0, a_max=None)) / self.keyword_dict['texpo_time']
            log.info("Computation 2 of RMS - image characteristics of readnoise, no. of exposures, total exposure time")
            log.info(f"    Minimum RMS: {minimum_rms:.6f}")

            # Computation 3 for rms: MAD (median absolute deviation)
            #
            # Create a version of the input image with zeros and NaN values masked
            masked_imgdata = np.ma.masked_array(imgdata, mask=(imgdata == 0) | np.isnan(imgdata))
            diffs = 2*masked_imgdata[:,2:-2] - masked_imgdata[:,:-4] - masked_imgdata[:,4:]
           
            # Median for flattened array yields the rms estimate
            # Ref: Stoehr et al. 2007, See https://ui.adsabs.harvard.edu/abs/2008ASPC..394..505S
            #
            # 1.482602 = inverse of the expected median value of a normal Gaussian distribution
            # sqrt(6.0) = comes from the propagation of errors for the given pixel combination.
            #    Ref: Rick White - If you have a weighted sum of 3 random numbers, s = a*r1 + b*r2 + c*r3,
            #    then the noise in s is sqrt(a**2+b**2+c**2). In the DER_SNR equation (see reference
            #    paper), the three coefficients are 2, 1 and 1, so the combination is 
            #    sqrt(4 + 1 + 1) = sqrt(6). 
            mad_rms = np.ma.median(np.ma.abs(diffs)) * 1.482602 / np.sqrt(6.0)
            log.info("Computation 3 of RMS - median absolute algorithm")
            log.info(f"    MAD rms: {mad_rms:.6f}")

            # Compare the computed rms, the minimum rms based upon detector characteristics, and
            # the MAD computation of the rms - use the largest of the three values
            bkg_rms = max(bkg_rms, minimum_rms, mad_rms)
            log.info(f"Final Background statistics - mean: {bkg_mean:.6f}  median: {bkg_median:.6f}  rms: {bkg_rms:.6f}")
            log.info("")

            # Generate two-dimensional background and rms images with the attributes of
            # the input data, but the content based on the sigma-clipped statistics.
            # bkg_median ==> background and bkg_rms ==> background rms
            self.bkg_background_ra = np.full_like(imgdata, bkg_median)
            self.bkg_rms_ra = np.full_like(imgdata, bkg_rms)
            self.bkg_median = bkg_median
            self.bkg_rms_median = bkg_rms
            self.bkg_type = 'sigma_clipped_background'
            negative_threshold = negative_sigma * bkg_rms

        # BACKGROUND COMPUTATION 3 (Background2D)
        # The simple_bkg = True is the way to force the background to be computed with the
        # sigma-clipped algorithm, regardless of any other criterion. If simple_bkg == True,
        # the compute_background() is done, otherwise try to use Background2D to compute the background.
        if not simple_bkg and not is_zero_background_defined:

            # If the sigma-clipped background image skew is greater than the threshold,
            # compute a two-dimensional background fit.  A larger skew implies
            # more sources in the field, which requires a more complex background.
            if bkg_skew > bkg_skew_threshold:
                log.info("Computing the background using the Background2D algorithm.")

                exclude_percentiles = [10, 25, 50, 75]
                for percentile in exclude_percentiles:
                    log.info("Percentile in use: {}".format(percentile))
                    try:
                        bkg = Background2D(imgdata, (box_size, box_size), filter_size=(win_size, win_size),
                                           bkg_estimator=bkg_estimator(),
                                           bkgrms_estimator=rms_estimator(),
                                           exclude_percentile=percentile,
                                           coverage_mask=self.inv_footprint_mask)

                    except Exception:
                        bkg = None
                        continue

                    if bkg is not None:
                        bkg_background_ra = bkg.background
                        bkg_rms_ra = bkg.background_rms
                        bkg_rms_median = bkg.background_rms_median
                        bkg_median = bkg.background_median
                        negative_threshold = negative_sigma * bkg.background_rms_median
                        break

                log.info(f"Background2D - Background median: {bkg_median:.6f}  rms: {bkg_rms_median:.6f}")
                log.info("")

                # If computation of a two-dimensional background image were successful, compute the
                # background-subtracted image and evaluate it for the number of negative values.
                #
                # If bkg is None, use the sigma-clipped statistics for the background.
                # If bkg is not None, but the background-subtracted image is too negative, use the
                # sigma-clipped computation for the background.
                if bkg is not None:
                    imgdata_bkgsub = imgdata - bkg_background_ra

                    # Determine how much of the illuminated portion of the background subtracted
                    # image is negative
                    num_negative = np.count_nonzero(imgdata_bkgsub[self.footprint_mask] < negative_threshold)
                    negative_ratio = num_negative / num_of_illuminated_pixels
                    del imgdata_bkgsub

                    # Report this information so the relative percentage and the threshold are known
                    log.info(f"Percentage of negative values in the background subtracted image {100.0*negative_ratio:.2f} vs low threshold of {negative_percent:.2f}.")

                    # If the background subtracted image has too many negative values which may be
                    # indicative of large negative regions, the two-dimensional computed background
                    # fit image should NOT be used.  Use the sigma-clipped data instead.
                    if negative_ratio * 100.0 > negative_percent:
                        log.info(f"Percentage of negative values {100.0*negative_ratio:.2f} in the background subtracted image exceeds the threshold of {negative_percent:.2f}.")
                        log.info("")
                        log.info("*** Use the background image determined from the sigma_clip algorithm. ***")

                    # If the rms computed for the Background2D is larger than the rms computed
                    # for the sigma-clipped algorithm, use the sigma-clipped algorithm instead.
                    # The sigma-clipped values have already been set in BACKGROUND COMPUTATION 2
                    # as the class attributes.
                    elif bkg_rms_median > self.bkg_rms_median:
                        log.info(f"Background2D rms, {bkg_rms_median:.6f}, is greater than sigma-clipped rms, {self.bkg_rms_median:.6f}.")
                        log.info("*** Use the background image determined from the sigma_clip algorithm. ***")

                    # Update the class variables with the background fit data
                    else:
                        self.bkg_background_ra = bkg_background_ra.copy()
                        self.bkg_rms_ra = bkg_rms_ra.copy()
                        self.bkg_rms_median = bkg_rms_median
                        self.bkg_median = bkg_median
                        self.bkg_type = 'twod_background'
                        log.info("")
                        log.info("*** Use the background image determined from the Background2D. ***")

                    del bkg_background_ra, bkg_rms_ra

            # Skewness of sigma_clipped background exceeds threshold
            else:
                log.info("*** Use the background image determined from the sigma_clip algorithm based upon skewness. ***")

        # User requested simple background == sigma_clip algorithm
        else:
            log.info("*** User requested the sigma_clip algorithm to determine the background image. ***")

        log.info("")
        log.info("Computation of image background complete")
        log.info("Found: ")
        log.info(f"    Median background: {self.bkg_median:.6f}")
        log.info(f"    Median RMS background: {self.bkg_rms_median:.6f}")
        log.info("")

        del bkg, imgdata

    def _get_header_data(self):
        """Read FITS keywords from the primary or extension header and store the
        information in a dictionary

        Returns
        -------
        keyword_dict : dictionary
            dictionary of keyword values
        """

        keyword_dict = {}

        keyword_dict["proposal_id"] = self.imghdu[0].header["PROPOSID"]
        keyword_dict["image_file_name"] = self.imghdu[0].header['FILENAME'].upper()
        keyword_dict["target_name"] = self.imghdu[0].header["TARGNAME"].upper()
        keyword_dict["date_obs"] = self.imghdu[0].header["DATE-OBS"]
        keyword_dict["time_obs"] = self.imghdu[0].header["TIME-OBS"]
        keyword_dict["instrument"] = self.imghdu[0].header["INSTRUME"].upper()
        keyword_dict["detector"] = self.imghdu[0].header["DETECTOR"].upper()
        keyword_dict["target_ra"] = self.imghdu[0].header["RA_TARG"]
        keyword_dict["target_dec"] = self.imghdu[0].header["DEC_TARG"]
        keyword_dict["expo_start"] = self.imghdu[0].header["EXPSTART"]
        keyword_dict["texpo_time"] = self.imghdu[0].header["TEXPTIME"]
        keyword_dict["exptime"] = self.imghdu[0].header["EXPTIME"]
        keyword_dict["ndrizim"] = self.imghdu[0].header["NDRIZIM"]
        keyword_dict["numexp"] = self.imghdu[0].header["NUMEXP"]
        if keyword_dict["detector"].upper() != "SBC":
            if keyword_dict["instrument"].upper() == 'WFPC2':
                atodgn = self._get_max_key_value(self.imghdu[0].header, 'ATODGAIN')
                keyword_dict["ccd_gain"] = atodgn
                keyword_dict["atodgn"] = atodgn
                keyword_dict["readnse"] = 5.0
                keyword_dict["gain_keys"] = [self.imghdu[0].header[k[:8]] for k in self.imghdu[0].header["ATODGAI*"]]
            else:
                keyword_dict["ccd_gain"] = self.imghdu[0].header["CCDGAIN"]
                keyword_dict["readnse"] = self._get_max_key_value(self.imghdu[0].header, 'READNSE')
                keyword_dict["atodgn"] = self._get_max_key_value(self.imghdu[0].header, 'ATODGN')
                keyword_dict["gain_keys"] = [self.imghdu[0].header[k[:8]] for k in self.imghdu[0].header["ATODGN*"]]

        keyword_dict["aperture_pa"] = self.imghdu[0].header["PA_V3"]

        # The total detection product has the FILTER keyword in
        # the primary header - read it for any instrument.
        #
        # For the filter detection product:
        # WFC3 only has FILTER, but ACS has FILTER1 and FILTER2
        # in the primary header.
        if self.ghd_product.lower() == "tdp":
            keyword_dict["filter1"] = self.imghdu[0].header["FILTER"]
        # The filter detection product...
        else:
            if keyword_dict["instrument"] == "ACS":
                keyword_dict["filter1"] = self.imghdu[0].header["FILTER1"]
                keyword_dict["filter2"] = self.imghdu[0].header["FILTER2"]
            elif keyword_dict["instrument"] == "WFPC2":
                keyword_dict["filter1"] = self.imghdu[0].header["FILTNAM1"]
                keyword_dict["filter2"] = self.imghdu[0].header["FILTNAM2"]
            else:
                keyword_dict["filter1"] = self.imghdu[0].header["FILTER"]
                keyword_dict["filter2"] = ""

        if keyword_dict["instrument"] == "ACS":
            keyword_dict["aperture_ra"] = self.imghdu[1].header["RA_APER"]
            keyword_dict["aperture_dec"] = self.imghdu[1].header["DEC_APER"]
        elif keyword_dict["instrument"] == "WFPC2":
            keyword_dict["aperture_ra"] = self.imghdu[0].header["RA_TARG"]
            keyword_dict["aperture_dec"] = self.imghdu[0].header["DEC_TARG"]
        else:
            keyword_dict["aperture_ra"] = self.imghdu[1].header["RA_APER"]
            keyword_dict["aperture_dec"] = self.imghdu[1].header["DEC_APER"]


        # Get the HSTWCS object from the first extension
        keyword_dict["wcs_name"] = self.imghdu[1].header["WCSNAME"]
        keyword_dict["wcs_type"] = self.imghdu[1].header["WCSTYPE"]
        keyword_dict["orientation"] = self.imghdu[1].header["ORIENTAT"]
        keyword_dict["photflam"] = proc_utils.find_flt_keyword(self.imghdu, "PHOTFLAM")
        keyword_dict["photplam"] = proc_utils.find_flt_keyword(self.imghdu, "PHOTPLAM")

        # Include WCS keywords as well for use in updating the output catalog metadata
        drzwcs = HSTWCS(self.imghdu, ext=1)
        wcshdr = drzwcs.wcs2header()  # create FITS keywords with CD matrix
        # add them to the keyword_dict preserving
        keyword_dict['wcs'] = {}
        keyword_dict['wcs'].update(wcshdr)

        return keyword_dict


    def _get_max_key_value(self, header, root_of_keyword):
        """Read FITS keywords with the same prefix from primary header and return the maximum value

        Parameters
        ----------
        header : hdu
            The header of a FITS hdu

        root_of_keyword : str
            The common root portion of a FITS keyword  (e.g., READNSE for READNSE[A-D])

        Returns
        -------
        max_value : float
            The maximum value or 1.0 of the keywords examined
        """

        max_value = max(header[root_of_keyword + "*"].values(), default=1.0)

        return max_value


class HAPCatalogs:
    """Generate photometric sourcelist for specified TOTAL or FILTER product image.
    """

    # This reference for the crfactor is Rick White.  It is based upon the physical
    # area of the detector.  From the instrument handbooks:
    # ACS/WFC: 4096**2 pixels, pixel size (15 um)**2
    # ACS/HRC: 1024**2 pixels, pixel size (21 um)**2
    # WFC3/UVIS: 4096**2 pixels, pixel size (15 um)**2
    # WFPC2/PC: 1600**2 pixels, pixel size (15um)**2
    # acs/wfc-> crfactor = {'aperture': 300, 'segment': 150}  # CRs / hr / 4kx4k pixels
    crfactor = {'WFC': {'aperture': 300, 'segment': 150}, 'HRC': {'aperture': 37, 'segment': 18.5},
                'UVIS': {'aperture': 300, 'segment': 150}, 'PC': {'aperture': 46, 'segment': 23}}

    def __init__(self, fitsfile, param_dict, param_dict_qc, num_images_mask, log_level, diagnostic_mode=False, types=None,
                 tp_sources=None):
        # set logging level to user-specified level
        log.setLevel(log_level)

        self.label = "HAPCatalogs"
        self.description = "A class used to generate photometric sourcelists using aperture photometry"

        self.imgname = fitsfile
        self.param_dict = param_dict
        self.param_dict_qc = param_dict_qc
        self.diagnostic_mode = diagnostic_mode
        self.tp_sources = tp_sources  # <---total product catalogs.catalogs[*].sources
        # Determine what types of catalogs have been requested
        if not isinstance(types, list) and types in [None, 'both']:
            types = CATALOG_TYPES

        elif types == 'aperture' or types == 'segment':
            types = [types]
        else:
            if any([t not in CATALOG_TYPES for t in types]):
                log.error("Catalog types {} not supported. Only {} are valid.".format(types, CATALOG_TYPES))
                raise ValueError

        self.types = types

        # Get various configuration variables needed for the background computation
        # Compute the background for this image
        self.image = CatalogImage(fitsfile, num_images_mask, log_level)
        self.image.compute_background(self.param_dict['bkg_box_size'],
                                      self.param_dict['bkg_filter_size'],
                                      simple_bkg=self.param_dict['simple_bkg'],
                                      bkg_skew_threshold=self.param_dict['bkg_skew_threshold'],
                                      zero_percent=self.param_dict['zero_percent'],
                                      negative_percent=self.param_dict['negative_percent'],
                                      nsigma_clip=self.param_dict['nsigma_clip'],
                                      maxiters=self.param_dict['maxiters'])

        self.image.build_kernel(self.param_dict['bkg_box_size'], self.param_dict['bkg_filter_size'],
                                self.param_dict['TWEAK_FWHMPSF'],
                                self.param_dict['simple_bkg'],
                                self.param_dict['bkg_skew_threshold'],
                                self.param_dict['zero_percent'],
                                self.param_dict['negative_percent'],
                                self.param_dict['nsigma_clip'],
                                self.param_dict['maxiters'],
                                self.param_dict['good_fwhm'])

        # Initialize all catalog types here...
        # This does NOT identify or measure sources to create the catalogs at this point...
        # The syntax here is EXTREMELY cludgy, but until a more compact way to do this is found,
        #  it will have to do...
        self.reject_cats = {}
        self.catalogs = {}
        if 'segment' in self.types:
            self.catalogs['segment'] = HAPSegmentCatalog(self.image, self.param_dict, self.param_dict_qc,
                                                         self.diagnostic_mode, tp_sources=tp_sources)
            self.reject_cats['segment'] = False

        if 'aperture' in self.types:
            self.catalogs['aperture'] = HAPPointCatalog(self.image, self.param_dict, self.param_dict_qc,
                                                        self.diagnostic_mode, tp_sources=tp_sources)
            self.reject_cats['aperture'] = False

        self.filters = {}

    def identify(self, **pars):
        """Build catalogs for this image.

        Parameters
        ----------
        types : list
            List of catalog types to be generated.  If None, build all available catalogs.
            Supported types of catalogs include: 'aperture', 'segment'.
        """
        # Support user-input value of 'None' which will trigger generation of all catalog types
        for catalog in self.catalogs:
            log.info("")
            log.info("Identifying {} sources".format(catalog))
            self.catalogs[catalog].identify_sources(**pars)

    def verify_crthresh(self, n1_exposure_time):
        """Verify whether catalogs meet cosmic-ray threshold limits.

        ... note : This algorithm has been modified from the original implementation
                   where if either catalog failed the cosmic ray criterion test, then
                   both catalogs would be rejected.

                   Instead...
                   an empty catalog (n_sources = 0) should not get rejected by the CR
                   contamination test (since it is empty, it clearly is not contaminated
                   by cosmic rays!) Also, if a catalog is empty because it is rejected
                   for some other reason, it should never trigger the rejection of the
                   other type of catalog.

                   thresh = crfactor * n1_exposure_time**2 / texptime

                   Since catalog output (*.ecsv) files are always written, rejection
                   means that all of the rows of measurements are deleted from the
                   output ECSV file.

        """
        for cat_type in self.catalogs:
            crthresh_mask = None
            source_cat = self.catalogs[cat_type].sources if cat_type == 'aperture' else self.catalogs[cat_type].source_cat

            # Collect up the names of all the columns which start with "Flag" - there will be
            # a flag column for each filter in the visit - even the filters which did not end up
            # contributing to the total detection image.
            flag_cols = [colname for colname in source_cat.colnames if colname.startswith('Flag')]
            for colname in flag_cols:
                catalog_crmask = source_cat[colname] < 2
                if crthresh_mask is None:
                    crthresh_mask = catalog_crmask
                else:
                    # Combine masks for all filters for this catalog type
                    crthresh_mask = np.bitwise_or(crthresh_mask, catalog_crmask)
            source_cat.sources_num_good = len(np.where(crthresh_mask)[0])

        log.info("Determining whether point and/or segment catalogs meet cosmic-ray threshold")
        log.info("  based on EXPTIME = {}sec for the n=1 filters".format(n1_exposure_time))

        detector = self.image.keyword_dict["detector"].upper()
        if detector not in ['IR', 'SBC']:
            for cat_type in self.catalogs:
                source_cat = self.catalogs[cat_type]
                if source_cat.sources:
                    thresh = self.crfactor[detector][cat_type] * n1_exposure_time**2 / self.image.keyword_dict['texpo_time']
                    source_cat = source_cat.sources if cat_type == 'aperture' else source_cat.source_cat
                    n_sources = source_cat.sources_num_good
                    all_sources = len(source_cat)
                    log.info("{} catalog with {} good sources out of {} total sources :  CR threshold = {}".format(cat_type, n_sources, all_sources, thresh))
                    if 0 < n_sources < thresh:
                        self.reject_cats[cat_type] = True
                        log.info("{} catalog FAILED CR threshold.".format(cat_type))

        # Ensure if any catalog is rejected, the remaining catalogs are also rejected
        if any(self.reject_cats.values()):
            for k in self.reject_cats.keys():
                self.reject_cats[k] = True

        return self.reject_cats

    def measure(self, filter_name, **pars):
        """Perform photometry and other measurements on sources for this image.

        Parameters
        ----------
        types : list
            List of catalog types to be generated.  If None, build all available catalogs.
            Supported types of catalogs include: 'aperture', 'segment'.
        """
        # Make sure we at least have a default 2D background computed
        for catalog in self.catalogs.values():
            if catalog.sources is None:
                catalog.identify_sources(**pars)

        for catalog in self.catalogs.values():
            catalog.measure_sources(filter_name, **pars)

        for catalog in self.catalogs.values():
            catalog.image.close()

    def write(self, reject_catalogs, **pars):
        """Write catalogs for this image to output files.

        Parameters
        ----------
        reject_catalogs : dictionary
            Dictionary where the keys are the types of catalogs, and the values are
            bools indicating whether or not the catalogs (*.ecsv) should be written
            with their full contents or as empty catalogs.

        types : list
            List of catalog types to be generated.  If None, build all available catalogs.
            Supported types of catalogs include: 'aperture', 'segment'.
        """
        # Make sure we at least have a default 2D background computed
        for catalog in self.catalogs.values():
            if catalog.source_cat is None:
                catalog.source_cat = catalog.sources
            catalog.write_catalog(reject_catalogs)

    def combine(self, subset_dict):
        """Combine subset columns from the filter catalog with the total detection catalog.

        Parameters
        ----------
        subset_dict : dictionary
           Dictionary where the keys are the types of catalogs, and the values are
           the catalog objects.

        """
        for k, v in self.catalogs.items():
            v.combine_tables(subset_dict[k]['subset'])


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
class HAPCatalogBase:
    """Virtual class used to define API for all catalogs"""
    catalog_suffix = ".ecsv"
    catalog_region_suffix = ".reg"
    catalog_format = "ascii.ecsv"
    catalog_type = None

    def __init__(self, image, param_dict, param_dict_qc, diagnostic_mode, tp_sources):
        self.image = image
        self.imgname = image.imgname
        self.param_dict = param_dict
        self.param_dict_qc = param_dict_qc
        self.diagnostic_mode = diagnostic_mode

        self.sourcelist_filename = self.imgname.replace(self.imgname[-9:], self.catalog_suffix)

        # Set the gain for ACS/SBC and WFC3/IR to 1.0
        if self.image.keyword_dict["detector"].upper() in ["IR", "SBC"]:
            self.gain = 1.0
            gain_values = 1.0
        else:
            # Compute average gain - there will always be at least one gain value in the primary header
            gain_keys = self.image.keyword_dict['gain_keys']
            gain_values = [g for g in gain_keys if g > 0.0]
            self.gain = self.image.keyword_dict['exptime'] * np.mean(gain_values)

        # Convert photometric aperture radii from arcsec to pixels
        self.aper_radius_arcsec = [self.param_dict['aperture_1'], self.param_dict['aperture_2']]
        self.aper_radius_list_pixels = []
        for aper_radius in self.aper_radius_arcsec:
            self.aper_radius_list_pixels.append(aper_radius / self.image.imgwcs.pscale)

        # Photometric information
        if not tp_sources:
            log.info("Average gain of {} for input image {}".format(np.mean(gain_values), self.imgname))
            log.info("{}".format("=" * 80))
            log.info("")
            log.info("")
            log.info("SUMMARY OF INPUT PARAMETERS FOR PHOTOMETRY")
            log.info("image name:       {}".format(self.imgname))
            log.info("platescale:       {}".format(self.image.imgwcs.pscale))
            log.info("radii (pixels):   {}".format(self.aper_radius_list_pixels))
            log.info("radii (arcsec):   {}".format(self.aper_radius_arcsec))
            log.info("annulus:          {}".format(self.param_dict['skyannulus_arcsec']))
            log.info("dSkyAnnulus:      {}".format(self.param_dict['dskyannulus_arcsec']))
            log.info("salgorithm:       {}".format(self.param_dict['salgorithm']))
            log.info("gain:             {}".format(self.gain))
            # log.info("ab_zeropoint:     {}".format(self.ab_zeropoint))
            log.info(" ")
            log.info("{}".format("=" * 80))
            log.info("")

        # Initialize attributes which are computed by class methods later
        self.sources = None  # list of identified source positions
        self.source_cat = None  # catalog of sources and their properties
        self.tp_sources = tp_sources

        # Determine what regions we have for source identification
        # Regions are defined as sections of the image which has the same
        # max WHT within a factor of 2.0 (or so).
        # make_wht_masks(whtarr, maskarr, scale=1.5, sensitivity=0.95, kernel=(11,11))
        self_scale = (self.image.keyword_dict['ndrizim'] - 1) / 2
        scale = max(self.param_dict['scale'], self_scale)
        self.tp_masks = None
        if not self.image.blank:
            self.tp_masks = make_wht_masks(
                self.image.wht_image, 
                self.image.inv_footprint_mask,
                scale=scale,
                sensitivity=self.param_dict['sensitivity'],
                kernel=(self.param_dict['region_size'],
                self.param_dict['region_size'])
                )

    def identify_sources(self, **pars):
        pass

    def measure_sources(self, filter_name, **pars):
        pass

    def write_catalog(self, reject_catalogs, **pars):
        pass

    def combine_tables(self, subset_dict):
        pass

    def annotate_table(self, data_table, param_dict_qc, proc_type="aperture", product="tdp"):
        """Add state metadata to the top of the output source catalog.

        Parameters
        ----------
        data_table : QTable
            Table of source properties

        param_dict_qc : dictionary
            Configuration values for quality control step based upon input JSON files (used to build catalog header)

        proc_type : str, optional
            Identification of catalog type: aperture (aka point) or segment

        product : str, optional
            Identification string for the catalog product being written.  This
            controls the data being put into the catalog product

        Returns
        -------
        data_table : QTable
            Table of source properties updatd to contain state metadata

        """
        data_table.meta["h00"] = [" #=================================================================================================="]
        data_table.meta["h01"] = [" # All refereed publications based on data obtained from the HAP must carry the following footnote: "]
        data_table.meta["h02"] = [" #                                                                                                  "]
        data_table.meta["h03"] = [" #     Based on observations made with the NASA/ESA Hubble Space Telescope                          "]
        data_table.meta["h04"] = [" #     and obtained from the Hubble Advanced Products collection generated                          "]
        data_table.meta["h05"] = [" #     by the Space Telescope Science Institute (STScI/NASA).                                       "]
        data_table.meta["h06"] = [" #                                                                                                  "]
        data_table.meta["h07"] = [" # One copy of each paper resulting from data obtained from the HAP should be sent to the STScI.    "]
        data_table.meta["h08"] = [" #=================================================================================================="]

        data_table.meta["WCSNAME"] = self.image.keyword_dict["wcs_name"]
        data_table.meta["WCSTYPE"] = self.image.keyword_dict["wcs_type"]
        data_table.meta["Proposal ID"] = self.image.keyword_dict["proposal_id"]
        data_table.meta["Image File Name"] = self.image.keyword_dict['image_file_name']
        data_table.meta["Target Name"] = self.image.keyword_dict["target_name"]
        data_table.meta["Date Observed"] = self.image.keyword_dict["date_obs"]
        data_table.meta["Time Observed"] = self.image.keyword_dict["time_obs"]
        data_table.meta["Instrument"] = self.image.keyword_dict["instrument"]
        data_table.meta["Detector"] = self.image.keyword_dict["detector"]
        data_table.meta["Target RA"] = self.image.keyword_dict["target_ra"]
        data_table.meta["Target DEC"] = self.image.keyword_dict["target_dec"]
        data_table.meta["Orientation"] = self.image.keyword_dict["orientation"]
        data_table.meta["Aperture RA"] = self.image.keyword_dict["aperture_ra"]
        data_table.meta["Aperture DEC"] = self.image.keyword_dict["aperture_dec"]
        data_table.meta["Aperture PA"] = self.image.keyword_dict["aperture_pa"]

        data_table.meta["Aper1 (arcsec)"] = self.aper_radius_arcsec[0]
        data_table.meta["Aper2 (arcsec)"] = self.aper_radius_arcsec[1]
        if proc_type == "segment":
             data_table.meta["Threshold (sigma)"] = self.final_nsigma
        else:
             data_table.meta["Threshold (sigma)"] = self.param_dict['nsigma']

        data_table.meta["Exposure Start"] = self.image.keyword_dict["expo_start"]
        data_table.meta["Total Exposure Time"] = self.image.keyword_dict["texpo_time"]
        if self.image.keyword_dict["detector"].upper() != "SBC":
            data_table.meta["CCD Gain"] = self.image.keyword_dict["ccd_gain"]
        if product.lower() == "tdp" or self.image.keyword_dict["instrument"].upper() == "WFC3":
            data_table.meta["Filter 1"] = self.image.keyword_dict["filter1"]
            data_table.meta["Filter 2"] = ""
        else:
            data_table.meta["Filter 1"] = self.image.keyword_dict["filter1"]
            data_table.meta["Filter 2"] = self.image.keyword_dict["filter2"]

        # Insert WCS keywords into metadata
        # This relies on the fact that the dicts are always ordered.
        data_table.meta.update(self.image.keyword_dict['wcs'])

        num_sources = len(data_table)
        data_table.meta["Number of sources"] = num_sources

        proc_type = proc_type.lower()
        ci_lower = float(param_dict_qc['ci filter'][proc_type]['ci_lower_limit'])
        ci_upper = float(param_dict_qc['ci filter'][proc_type]['ci_upper_limit'])

        data_table.meta["h09"] = ["#================================================================================================="]
        data_table.meta["h10"] = ["IMPORTANT NOTES"]

        data_table.meta["h10.1"] = ["These catalogs provide aperture photometry in the ABMAG system and are calibrated with"]
        data_table.meta["h10.2"] = ["photometric zeropoints corresponding to an infinite aperture. To convert to total"]
        data_table.meta["h10.3"] = ["magnitudes, aperture corrections must be applied to account for flux falling outside of"]
        data_table.meta["h10.4"] = ["the selected aperture. For details, see Whitmore et al., 2016 AJ, 151, 134W."]
        data_table.meta["h10.5"] = [" "]

        data_table.meta["h11"] = ["The X and Y coordinates in this table are 0-indexed (i.e. the origin is (0,0))."]
        data_table.meta["h12"] = ["RA and Dec values in this table are in sky coordinates (i.e. coordinates at the epoch of observation"]
        data_table.meta["h12.1"] = ["and an {}).".format(self.image.keyword_dict["wcs_type"])]
        data_table.meta["h13"] = ["Magnitude values in this table are in the ABMAG system."]
        data_table.meta["h14"] = ["Column titles in this table ending with Ap1 refer to the inner photometric aperture "]
        data_table.meta["h14.1"] = ["(radius = {} pixels, {} arcsec.".format(self.aper_radius_list_pixels[0],
                                                                             self.aper_radius_arcsec[0])]
        data_table.meta["h15"] = ["Column titles in this table ending with Ap2 refer to the outer photometric aperture "]
        data_table.meta["h15.1"] = ["(radius = {} pixels, {} arcsec.".format(self.aper_radius_list_pixels[1],
                                                                             self.aper_radius_arcsec[1])]
        data_table.meta["h16"] = ["CI = Concentration Index (CI) = MagAp1 - MagAp2."]
        data_table.meta["h17"] = ["Flag Value Identification:"]
        data_table.meta["h17.1"] = ["    0 - Stellar Source ({} < CI < {})".format(ci_lower, ci_upper)]
        data_table.meta["h17.2"] = ["    1 - Extended Source (CI > {})".format(ci_upper)]
        data_table.meta["h17.3"] = ["    2 - Questionable Photometry (Single-Pixel Saturation)"]
        data_table.meta["h17.4"] = ["    4 - Questionable Photometry (Multi-Pixel Saturation)"]
        data_table.meta["h17.5"] = ["    8 - Faint Detection Limit"]
        data_table.meta["h17.6"] = ["   16 - Hot pixels (CI < {})".format(ci_lower)]
        data_table.meta["h17.7"] = ["   32 - False Detection Swarm Around Saturated Source"]
        data_table.meta["h17.8"] = ["   64 - False Detections Near Image Edge"]
        data_table.meta["h17.9"] = ["  128 - Bleeding and Cosmic Rays"]
        data_table.meta["h18"] = ["#================================================================================================="]

        if proc_type == "segment":
            if self.is_big_island:
                data_table.meta["h19"] = ["WARNING: Segmentation catalog is considered to be of poor quality due to a crowded field or large segments."]

        if type(data_table.meta) != OrderedDict:
            data_table.meta = OrderedDict(data_table.meta)

        return (data_table)


# --------------------------------------------------------------------------------------------------------

class HAPPointCatalog(HAPCatalogBase):
    """Generate photometric sourcelist(s) for specified image(s) using aperture photometry of point sources.
    """
    catalog_suffix = "_point-cat.ecsv"
    catalog_type = 'aperture'

    def __init__(self, image, param_dict, param_dict_qc, diagnostic_mode, tp_sources):
        super().__init__(image, param_dict, param_dict_qc, diagnostic_mode, tp_sources)

        # Defined in measure_sources
        self.subset_filter_source_cat = None

    def identify_sources(self, **pars):
        """Create a master coordinate list of sources identified in the specified total detection product image
        """
        # Recognize empty/blank images and treat gracefully
        if self.image.blank:
            self._define_empty_table()
            log.info('Total Detection Product - Point-source catalog is EMPTY')
            return
        source_fwhm = self.image.kernel_fwhm
        # read in sci, wht extensions of drizzled product
        image = np.nan_to_num(self.image.data, copy=True, nan=0.0)

        # Create the background-subtracted image
        image -= self.image.bkg_background_ra
        image = np.clip(image, 0, image.max())  # Insure there are no neg pixels to trip up StarFinder

        if 'drz.fits' in self.image.imgname:
            reg_suffix = 'drz.fits'
        else:
            reg_suffix = 'drc.fits'

        if not self.tp_sources:
            # Report configuration values to log
            log.info("{}".format("=" * 80))
            log.info("")
            log.info("Point-source finding settings")
            log.info("Total Detection Product - Input Parameters")
            log.info("INPUT PARAMETERS")
            log.info("image name: {}".format(self.imgname))
            log.info("{}: {}".format("self.param_dict['dao']['bkgsig_sf']", self.param_dict["dao"]["bkgsig_sf"]))
            log.info("{}: {}".format("self.param_dict['dao']['kernel_sd_aspect_ratio']",
                                     self.param_dict['dao']['kernel_sd_aspect_ratio']))
            log.info("{}: {}".format("self.param_dict['simple_bkg']", self.param_dict['simple_bkg']))
            log.info("{}: {}".format("self.param_dict['nsigma']", self.param_dict['nsigma']))
            log.info("{}: {}".format("self.image.bkg_rms_median", self.image.bkg_rms_median))
            log.info("DERIVED PARAMETERS")
            log.info("{}: {}".format("source_fwhm", source_fwhm))
            log.info("{}: {}".format("threshold", self.param_dict['nsigma'] * self.image.bkg_rms_median))
            log.info("")
            log.info("{}".format("=" * 80))

            sources = None
            last_used_id = 0
            for masknum, mask in enumerate(self.tp_masks):
                # apply mask for each separate range of WHT values
                region = image * mask['mask']

                # Compute separate threshold for each 'region'
                reg_rms = self.image.bkg_rms_ra * np.sqrt(mask['mask'] / mask['rel_weight'].max())
                reg_rms_median = np.nanmedian(reg_rms[reg_rms > 0])
                log.info("Mask {}: rel = {}".format(mask['wht_limit'], mask['rel_weight'].max()))

                # find ALL the sources!!!
                if self.param_dict["starfinder_algorithm"] == "dao":
                    log.info("DAOStarFinder(fwhm={}, threshold={}*{})".format(source_fwhm, self.param_dict['nsigma'],
                                                                              reg_rms_median))
                    daofind = DAOStarFinder(fwhm=source_fwhm,
                                            threshold=self.param_dict['nsigma'] * reg_rms_median)
                    reg_sources = daofind(region, mask=self.image.inv_footprint_mask)
                elif self.param_dict["starfinder_algorithm"] == "iraf":
                    log.info("IRAFStarFinder(fwhm={}, threshold={}*{})".format(source_fwhm, self.param_dict['nsigma'],
                                                                               reg_rms_median))
                    isf = IRAFStarFinder(fwhm=source_fwhm, threshold=self.param_dict['nsigma'] * reg_rms_median)
                    reg_sources = isf(region, mask=self.image.inv_footprint_mask)
                elif self.param_dict["starfinder_algorithm"] == "psf":
                    log.info("UserStarFinder(fwhm={}, threshold={}*{})".format(source_fwhm, self.param_dict['nsigma'],
                                                                              reg_rms_median))
                    # Perform manual detection of sources using theoretical PSFs
                    # Initial test data: ictj65
                    try:
                        # Subtract the detection threshold image so that detection is anything > 0
                        region -= (reg_rms * self.param_dict['nsigma'])
                        # insure no negative values for deconvolution
                        region = np.clip(region, 0., region.max())
                        user_peaks, source_fwhm = decutils.find_point_sources(self.image.imgname,
                                                                 data=region,
                                                                 def_fwhm=source_fwhm,
                                                                 box_size=self.param_dict['region_size'],
                                                                 mask=self.image.footprint_mask,
                                                                 block_size=self.param_dict['block_size'],
                                                                 diagnostic_mode=self.diagnostic_mode)
                    except Exception:
                        # In case we run into out-of-memory error, or any other exception with
                        # PSF use (like CTE or horribly mismatched PSFs), fail-over to using
                        # DAOFind mode instead
                        log.warning("Exception thrown when trying to use PSFs to find sources with UserStarFinder.")
                        user_peaks = None

                    if user_peaks is not None and len(user_peaks) > 0:

                        log.info("UserStarFinder identified {} sources".format(len(user_peaks)))
                        if self.diagnostic_mode:
                            peak_name = "{}_peaks{}.reg".format(self.image.imgname.split('.')[0], masknum)
                            peak_reg = user_peaks['x_peak', 'y_peak']
                            peak_reg['x_peak'] += 1
                            peak_reg['y_peak'] += 1
                            peak_reg.write(peak_name,
                                           format='ascii.fast_no_header',
                                           overwrite=True)

                        # NOTE: It is necessary to use a non-zero threshold to ensure the values
                        # computed in the low levels of Photutils are finite (in v1.13.0) or no
                        # table will be returned in the "src_table = daofind()" line below.
                        daofind = decutils.UserStarFinder(fwhm=source_fwhm,
                                                coords=user_peaks,
                                                threshold=1.0e-10,
                                                sharphi=0.9, sharplo=0.4)

                        _region_name = self.image.imgname.replace(reg_suffix, 'region{}.fits'.format(masknum))
                        if self.diagnostic_mode:
                            fits.PrimaryHDU(data=region).writeto(_region_name, overwrite=True)

                        reg_sources = daofind(region,
                                              mask=self.image.inv_footprint_mask)

                        _region_name = self.image.imgname.replace(reg_suffix, 'starfind_sources{}.ecsv'.format(masknum))
                        if self.diagnostic_mode:
                            reg_sources.write(_region_name, format='ascii.ecsv', overwrite=True)

                    else:
                        # No sources found to match the PSF model, perhaps due to CTE.
                        # Try standard daofind instead
                        log.info("Reverting to DAOStarFinder(fwhm={}, threshold={}*{})".format(source_fwhm, self.param_dict['nsigma'],
                                                                                  reg_rms_median))
                        daofind = DAOStarFinder(fwhm=source_fwhm,
                                                threshold=self.param_dict['nsigma'] * reg_rms_median)
                        reg_sources = daofind(region, mask=self.image.inv_footprint_mask)
                else:
                    err_msg = "'{}' is not a valid 'starfinder_algorithm' parameter input in the catalog_generation parameters json file. Valid options are 'dao' for photutils.detection.DAOStarFinder() or 'iraf' for photutils.detection.IRAFStarFinder().".format(self.param_dict["starfinder_algorithm"])
                    log.error(err_msg)
                    raise ValueError(err_msg)
                log.info("{}".format("=" * 80))
                # Concatenate sources found in each region.
                # NOTE: The unique "id" column must be updated manually so when the sources
                # and reg_sources are "vstacked", all the rows stil have a unique id.
                if reg_sources is not None:
                    # Convert the QTable to an Astropy Table
                    reg_sources = Table(reg_sources)

                    # Remove the extra column added in Photutils v2.0.  A "daofind_mag"
                    # column was added for comparison to the original IRAF DAOFIND algorithm.
                    if "daofind_mag" in reg_sources.colnames:
                        reg_sources.remove_column("daofind_mag")

                    if sources is None:
                        last_used_id = reg_sources['id'][-1]
                        sources = reg_sources
                    else:
                        reg_sources['id'][:] +=  last_used_id
                        sources = vstack([sources, reg_sources])
                        last_used_id = sources['id'][-1]

            # If there are no detectable sources in the total detection image, return as there is nothing more to do.
            if not sources:
                log.warning("No point sources were found in Total Detection Product, {}.".format(self.imgname))
                log.warning("Processing for point source catalogs for this product is ending.")
                self._define_empty_table()
                return

            log.info("Measured {} sources in {}".format(len(sources), self.image.imgname))
            log.info("   colnames: {}".format(sources.colnames))
            # insure centroid columns are not masked
            sources['xcentroid'] = Column(sources['xcentroid'])
            sources['ycentroid'] = Column(sources['ycentroid'])
            # calculate and add RA and DEC columns to table
            ra, dec = self.transform_list_xy_to_ra_dec(sources["xcentroid"], sources["ycentroid"], self.imgname)
            ra_col = Column(name="RA", data=ra, dtype=np.float64)
            dec_col = Column(name="DEC", data=dec, dtype=np.float64)
            sources.add_column(ra_col, index=3)
            sources.add_column(dec_col, index=4)

            for col in sources.colnames:
                if col != 'id':
                    sources[col].info.format = '.8g'  # for consistent table output

            # format output table columns
            final_col_format = {"xcentroid": "10.3f", "ycentroid": "10.3f", "RA": "13.7f", "DEC": "13.7f", "id": "7d"}
            for fcf_key in final_col_format.keys():
                sources[fcf_key].format = final_col_format[fcf_key]

            # descriptions
            final_col_descrip = {"xcentroid": "Pixel Coordinate", "ycentroid": "Pixel Coordinate",
                                 "RA": "Sky coordinate at epoch of observation",
                                 "DEC": "Sky coordinate at epoch of observation",
                                 "id": "Catalog Object Identification Number"}
            for fcd_key in final_col_descrip.keys():
                sources[fcd_key].description = final_col_descrip[fcd_key]

            # add units to columns
            final_col_units = {"xcentroid": "pixel", "ycentroid": "pixel", "RA": "degree", "DEC": "degree",
                               "id": ""}
            for col_title in final_col_units:
                sources[col_title].unit = final_col_units[col_title]

            if self.diagnostic_mode:
                sources.write(self.image.imgname.replace(reg_suffix,'raw-point-cat.ecsv'), format='ascii.ecsv', overwrite=True)
            self.sources = sources

        # if processing filter product, use sources identified by parent total drizzle product identify_sources() run
        if self.tp_sources:
            self.sources = self.tp_sources['aperture']['sources']

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def _define_empty_table(self):
        """Create basic empty table based on total product table to signify no valid source were detected"""
        final_col_format = {"xcentroid": "10.3f", "ycentroid": "10.3f",
                            "RA": "13.7f", "DEC": "13.7f",
                            "id": "7d", "Flags": "5d"}

        final_col_descrip = {"xcentroid": "Pixel Coordinate", "ycentroid": "Pixel Coordinate",
                             "RA": "Sky coordinate at epoch of observation",
                             "DEC": "Sky coordinate at epoch of observation",
                             "id": "Catalog Object Identification Number",
                             "Flags": "Numeric encoding for conditions on detected sources"}

        final_col_units = {"xcentroid": "pixel", "ycentroid": "pixel", "RA": "degree", "DEC": "degree",
                           "id": "", "Flags": ""}

        final_colnames = [k for k in final_col_format.keys()]

        # Initialize empty table with desired column names, descriptions and units
        empty_table = Table(names=final_colnames,
                      descriptions=final_col_descrip,
                      units=final_col_units)

        # Add formatting for each column
        for fcf_key in final_col_format.keys():
            empty_table[fcf_key].format = final_col_format[fcf_key]

        self.sources = empty_table


    def measure_sources(self, filter_name):
        """Perform aperture photometry on identified sources
        """
        if len(self.sources) == 0:
            # Report configuration values to log
            log.info("{}".format("=" * 80))
            log.info("")
            log.info("No point sources identified for photometry for")
            log.info("image name: {}".format(self.imgname))
            log.info("Generating empty point-source catalog.")
            log.info("")
            # define this attribute for use by the .write method
            self.source_cat = self.sources

            self.subset_filter_source_cat = Table(names=["ID", "MagAp2", "CI", "Flags"])
            self.subset_filter_source_cat.rename_column("MagAp2", "MagAP2_" + filter_name)
            self.subset_filter_source_cat.rename_column("CI", "CI_" + filter_name)
            self.subset_filter_source_cat.rename_column("Flags", "Flags_" + filter_name)

            return

        log.info("Performing aperture photometry on identified point-sources")
        # Open and background subtract image
        image = self.image.data.copy()

        # load in coords of sources identified in total product
        _cnames = ['X-Center', 'Y-Center'] if 'X-Center' in self.sources.colnames else ['xcentroid', 'ycentroid']
        if len(self.sources) > 1:
            positions = (self.sources[_cnames[0]], self.sources[_cnames[1]])
        else:
            positions = (Column(self.sources[_cnames[0]]), Column(self.sources[_cnames[1]]))

        pos_xy = np.vstack(positions).T

        # define list of background annulii
        bg_apers = CircularAnnulus(pos_xy,
                                   r_in=self.param_dict['skyannulus_arcsec']/self.image.imgwcs.pscale,
                                   r_out=(self.param_dict['skyannulus_arcsec'] +
                                          self.param_dict['dskyannulus_arcsec'])/self.image.imgwcs.pscale)

        # Create the list of photometric apertures to measure
        phot_apers = [CircularAperture(pos_xy, r=r) for r in self.aper_radius_list_pixels]

        # Perform aperture photometry - the input data should NOT be background subtracted
        photometry_tbl = photometry_tools.iraf_style_photometry(phot_apers,
                                                                bg_apers,
                                                                data=image,
                                                                photflam=self.image.keyword_dict['photflam'],
                                                                photplam=self.image.keyword_dict['photplam'],
                                                                error_array=self.image.bkg_rms_ra,
                                                                bg_method=self.param_dict['salgorithm'],
                                                                epadu=self.gain)

        # calculate and add RA and DEC columns to table
        ra, dec = self.transform_list_xy_to_ra_dec(photometry_tbl["X-Center"], photometry_tbl["Y-Center"], self.imgname)  # TODO: replace with all_pix2sky or somthing at a later date
        ra_col = Column(name="RA", data=ra, dtype=np.float64)
        dec_col = Column(name="DEC", data=dec, dtype=np.float64)
        photometry_tbl.add_column(ra_col, index=2)
        photometry_tbl.add_column(dec_col, index=3)
        log.info('Obtained photometry measurements for {} sources'.format(len(photometry_tbl)))

        try:
            # Calculate and add concentration index (CI) column to table
            if math.isclose(self.image.keyword_dict['photflam'], 0.0, abs_tol=constants.TOLERANCE):
                # Fill the ci_data column with the "bad data indicator" of -9999.0
                # where photometry_tbl["MagAp1"] was populated with -9999.0 in iraf_style_photometry
                ci_data = photometry_tbl["MagAp1"].data
                ci_mask = np.logical_not(ci_data > constants.FLAG + 1.0)
            else:
                ci_data = photometry_tbl["MagAp1"].data - photometry_tbl["MagAp2"].data
                ci_mask = np.logical_and(np.abs(ci_data) > 0.0, np.abs(ci_data) < 1.0e-30)
                big_bad_index = np.where(abs(ci_data) > 1.0e20)
                ci_mask[big_bad_index] = True

        except Exception as x_cept:
            log.warning("Computation of concentration index (CI) was not successful: {} - {}.".format(self.imgname, x_cept))
            log.warning("CI measurements may be missing from the output filter catalog.\n")

            # Create a column of data from an existing column of data
            # and set each value = -9999.0
            ci_data = constants.FLAG + photometry_tbl["MagAp1"] * 0.0
            ci_mask = np.logical_not(ci_data > constants.FLAG + 1.0)

        ci_col = MaskedColumn(name='CI', data=ci_data, dtype=np.float64, mask=ci_mask)
        photometry_tbl.add_column(ci_col)

        # Add zero-value "Flags" column in preparation for source flagging
        flag_col = Column(name="Flags", data=np.zeros_like(photometry_tbl['ID']), dtype=np.int64)
        photometry_tbl.add_column(flag_col)

        # build final output table
        final_col_order = ["X-Center", "Y-Center", "RA", "DEC", "ID", "MagAp1", "MagErrAp1", "FluxAp1", 
                           "FluxErrAp1", "MagAp2", "MagErrAp2", "MSkyAp2", "StdevAp2", "FluxAp2", 
                           "FluxErrAp2", "CI", "Flags"]
        output_photometry_table = photometry_tbl[final_col_order]

        # format output table columns
        final_col_format = {"X-Center": "10.3f", "Y-Center": "10.3f", "RA": "13.7f", "DEC": "13.7f", "ID": "7d",
                            "MagAp1": '7.3f', "MagErrAp1": '7.3f', "FluxAp1": '10.4f', "FluxErrAp1": '10.4f',
                            "MagAp2": '7.3f', "MagErrAp2": '7.3f', "MSkyAp2": '7.3f', "StdevAp2": '7.3f',
                            "FluxAp2": '10.4f', "FluxErrAp2": '10.4f', "CI": "7.3f", "Flags": "5d"}
        for fcf_key in final_col_format.keys():
            output_photometry_table[fcf_key].format = final_col_format[fcf_key]

        # column descriptions
        final_col_descrip = {"ID": "Catalog Object Identification Number",
                             "X-Center": "Pixel Coordinate",
                             "Y-Center": "Pixel Coordinate",
                             "RA": "Sky coordinate at epoch of observation",
                             "DEC": "Sky coordinate at epoch of observation",
                             "MagAp1": "ABMAG of source based on the inner (smaller) aperture",
                             "MagErrAp1": "Error of MagAp1",
                             "FluxAp1": "Flux of source based on the outer (smaller) aperture",
                             "FluxErrAp1": "Flux of source based on the outer (smaller) aperture",
                             "MagAp2": "ABMAG of source based on the outer (larger) aperture",
                             "MagErrAp2": "Error of MagAp2",
                             "MSkyAp2": "Sky estimate from an annulus outside Aperture 2",
                             "StdevAp2": "Standard deviation of sky estimate from annulus outside Aperture 2",
                             "FluxAp2": "Flux of source based on the outer (larger) aperture",
                             "FluxErrAp2": "Flux of source based on the outer (larger) aperture",
                             "CI": "Concentration Index",
                             "Flags": "Numeric encoding for conditions on detected sources"}
        for fcd_key in final_col_descrip.keys():
            output_photometry_table[fcd_key].description = final_col_descrip[fcd_key]

        # add units to columns
        final_col_units = {"X-Center": "pixel", "Y-Center": "pixel", "RA": "degree", "DEC": "degree",
                           "ID": "", "MagAp1": "mag(AB)", "MagErrAp1": "mag(AB)", "MagAp2": "mag(AB)",
                           "MagErrAp2": "mag(AB)", "MSkyAp2": "electron/(s pixel)", "StdevAp2": "electron/(s pixel)",
                           "FluxAp1": "electron/s", "FluxErrAp1": "electron/s","FluxAp2": "electron/s", 
                           "FluxErrAp2": "electron/s", "CI": "mag(AB)", "Flags": ""}
        for col_title in final_col_units:
            output_photometry_table[col_title].unit = final_col_units[col_title]

        # Capture specified columns in order to append to the total detection table
        self.subset_filter_source_cat = output_photometry_table["ID", "MagAp2", "CI", "Flags"]
        self.subset_filter_source_cat.rename_column("MagAp2", "MagAP2_" + filter_name)
        self.subset_filter_source_cat.rename_column("CI", "CI_" + filter_name)
        self.subset_filter_source_cat.rename_column("Flags", "Flags_" + filter_name)

        # Add the header information to the table
        self.source_cat = self.annotate_table(output_photometry_table,
                                              self.param_dict_qc,
                                              proc_type = "aperture",
                                              product=self.image.ghd_product)
        log.info("Saved photometry table with {} sources".format(len(self.source_cat)))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def write_catalog(self, reject_catalogs):
        """Write specified catalog to file on disk

        Regardless of the setting for reject_catalogs, the regions file will be written
        solely based upon the setting of diagnostic_mode.

        Parameters
        ----------
        reject_catalogs : dictionary
            Dictionary where the keys are the types of catalogs, and the values are
            bools indicating whether or not the catalogs (*.ecsv) should be written
            with their full contents or as empty catalogs.

        Returns
        -------
        Nothing!

        The total product catalog contains several columns contributed from each filter catalog.
        However, there is no guarantee that every filter catalog contains measurements for each
        source in the total product catalog which is to say the row is actually missing from the
        filter product catalog.  When the catalog tables are combined in the combined_tables method,
        the missing entries will contain masked values represented by dashes.  When these values are
        written to output ECSV files, they are written as empty ("") strings.  In order for the catalogs
        to be ingested into a database by possible downstream processing, the masked values will be
        replaced by a numeric indicator.

        ... note : A catalog can have no rows of measurements because no sources were
        found OR because the catalog was "rejected" according to the cosmic ray rejection
        criterion.
        """
        # Insure catalog has all necessary metadata
        self.source_cat = self.annotate_table(self.source_cat, self.param_dict_qc, proc_type="aperture",
                                              product=self.image.ghd_product)
        if reject_catalogs[self.catalog_type]:
            # We still want to write out empty files
            # This will delete all rows from the existing table
            self.source_cat.remove_rows(slice(0, None))

        # Fill the nans and masked values with numeric data
        self.source_cat = fill_nans_maskvalues (self.source_cat, fill_value=constants.FLAG)

        # Write out catalog to ecsv file
        # self.source_cat.meta['comments'] = \
        #     ["NOTE: The X and Y coordinates in this table are 0-indexed (i.e. the origin is (0,0))."]
        self.source_cat.write(self.sourcelist_filename, format=self.catalog_format, overwrite=True)
        log.info("Wrote catalog file '{}' containing {} sources".format(self.sourcelist_filename, len(self.source_cat)))

        # Write out region file if in diagnostic_mode.
        if self.diagnostic_mode:
            out_table = self.source_cat.copy()
            if 'xcentroid' in out_table.keys():  # for point-source source catalogs
                # Remove all other columns besides xcentroid and ycentroid
                out_table.keep_columns(['xcentroid', 'ycentroid'])
                # Add offset of 1.0 in X and Y to line up sources in region file with image displayed in ds9.
                out_table['xcentroid'].data[:] += np.float64(1.0)
                out_table['ycentroid'].data[:] += np.float64(1.0)
            elif 'X-Center' in out_table.keys():  # for aperture photometric catalogs
                # Remove all other columns besides 'X-Center and Y-Center
                out_table.keep_columns(['X-Center', 'Y-Center'])
                # Add offset of 1.0 in X and Y to line up sources in region file with image displayed in ds9.
                out_table['X-Center'].data[:] += np.float64(1.0)
                out_table['Y-Center'].data[:] += np.float64(1.0)
            else:  # Bail out if anything else is encountered.
                log.info("Error: unrecognized catalog format. Skipping region file generation.")
                return()
            reg_filename = self.sourcelist_filename.replace("." + self.catalog_suffix.split(".")[1],
                                                            self.catalog_region_suffix)
            out_table.write(reg_filename, format="ascii")
            log.info("Wrote region file '{}' containing {} sources".format(reg_filename, len(out_table)))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def transform_list_xy_to_ra_dec(self, list_of_x, list_of_y, drizzled_image):
        """Transform lists of X and Y coordinates to lists of RA and Dec coordinates
        This is a temporary solution until somthing like pix2sky or pix2world can be implemented in measure_sources.

        directly lifted from hla classic subroutine hla_sorucelist.Transform_list_xy_to_RA_Dec()

        Tested.

        Parameters
        ----------
        list_of_x : list
            list of x coordinates to convert

        list_of_y :
            list of y coordinates to convert

        drizzled_image : str
            Name of the image that corresponds to the table from DAOPhot. This image is used to re-write x and y
            coordinates in RA and Dec.

        Returns
        -------
        ra: list
            list of right ascension values

        dec : list
            list of declination values
        """
        import stwcs

        wcs1_drz = stwcs.wcsutil.HSTWCS(drizzled_image + "[1]")
        origin = 0
        # *origin* is the coordinate in the upper left corner of the
        # image.  In FITS and Fortran standards, this is 1.  In Numpy and C
        # standards this is 0.
        try:
            skyposish = wcs1_drz.all_pix2sky(list_of_x, list_of_y, origin)
        except AttributeError:
            skyposish = wcs1_drz.all_pix2world(list_of_x, list_of_y, origin)
        ra = skyposish[0]
        dec = skyposish[1]

        return ra, dec

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def combine_tables(self, subset_table):
        """Append specified measurements from the filter table to the total detection table.

        The "ID" column is used to map the filter table measurements to the total detection table

        Parameters
        ----------
        subset_table : Astropy table
            A table containing a subset of columns from a filter catalog.

        """
        # Evaluate self.sources (the total product list) even though len(self.sources) should not be possible
        if len(subset_table) == 0 or len(self.sources) == 0:
            log.error("No sources found in the current filter table nor in the total source table.")
            return

        # Keep all the rows in the original total detection table and add columns from the filter
        # table where a matching "id" key is present.  The key must match in case.
        if 'xcentroid' in self.sources.colnames:
            self.sources.rename_column('xcentroid', 'X-Center')
        if 'ycentroid' in self.sources.colnames:
            self.sources.rename_column('ycentroid', 'Y-Center')
        if 'id' in self.sources.colnames:
            self.sources.rename_column("id", "ID")
        for col2del in ['sharpness', 'roundness1', 'roundness2', 'npix', 'sky', 'peak', 'flux', 'mag']:
            if col2del in self.sources.colnames:
                self.sources.remove_column(col2del)

        # Cast the QTable back to a Table to avoid complications
        self.sources = join(Table(self.sources), subset_table, keys="ID", join_type="left")

# ----------------------------------------------------------------------------------------------------------------------


class HAPSegmentCatalog(HAPCatalogBase):
    """Generate a sourcelist for a specified image by detecting both point and extended
       sources using the image segmentation process.

       Parameters
       ----------
       image : CatalogImage object
           The white light (aka total detection) or filter drizzled image

       param_dict : dictionary
           Configuration values for catalog generation based upon input JSON files

       diagnostic_mode : bool
           Specifies whether or not to generate the regions file used for ds9 overlay

       tp_sources: dictionary
           Dictionary containing computed information for each catalog type
    """
    catalog_suffix = "_segment-cat.ecsv"
    catalog_type = 'segment'

    # Class variable which indicates to the Filter object the Total object had to determine
    # the image background by the sigma_clipped alternate algorithm
    using_sigma_clipped_bkg = False

    def __init__(self, image, param_dict, param_dict_qc, diagnostic_mode, tp_sources):
        super().__init__(image, param_dict, param_dict_qc, diagnostic_mode, tp_sources)

        # Get the instrument/detector-specific values from the self.param_dict
        self._size_source_box = self.param_dict["sourcex"]["source_box"]
        self._nlevels = self.param_dict["sourcex"]["nlevels"]
        self._contrast = self.param_dict["sourcex"]["contrast"]
        self._border = self.param_dict["sourcex"]["border"]
        self._nsigma = self.param_dict["sourcex"]["segm_nsigma"]
        self._rw2d_size = self.param_dict["sourcex"]["rw2d_size"]
        self._rw2d_nsigma = self.param_dict["sourcex"]["rw2d_nsigma"]
        self._rw2d_biggest_source = self.param_dict["sourcex"]["rw2d_biggest_source"]
        self._rw2d_biggest_pixels = self.param_dict["sourcex"]["rw2d_biggest_pixels"]
        self._rw2d_source_fraction = self.param_dict["sourcex"]["rw2d_source_fraction"]
        self._bs_deblend_limit = self.param_dict["sourcex"]["biggest_source_deblend_limit"]
        self._sf_deblend_limit = self.param_dict["sourcex"]["source_fraction_deblend_limit"]
        self._ratio_bigsource_limit = self.param_dict["sourcex"]["ratio_bigsource_limit"]
        self._kron_scaling_radius = self.param_dict["sourcex"]["kron_scaling_radius"]
        self._kron_minimum_radius = self.param_dict["sourcex"]["kron_minimum_radius"]
        self._ratio_bigsource_deblend_limit = self.param_dict["sourcex"]["ratio_bigsource_deblend_limit"]

        # Initialize attributes to be computed later
        self.segm_img = None  # Segmentation image

        # Image data convolved with the kernel
        self.convolved_img = None

        # Instance of a SourceCatalog object defined on the total image
        # but used for the filtered objects
        self.total_source_cat = None

        # Defined in measure_sources
        self.subset_filter_source_cat = None

        # Default kernel is a Gaussian 2D kernel. This may be over-ridden in identify_sources().
        self.kernel = copy.deepcopy(self.image.kernel)

        # Attribute computed when generating the segmentation image.  If the segmentation image
        # is deemed to be of poor quality, make sure to add documentation to the output catalog.
        self.is_big_island = False

        # Final sigma value used based upon Gaussian or RickerWavelet processing
        self.final_nsigma = None


    def identify_sources(self, **pars):
        """Use photutils to find sources in image based on segmentation.

        Returns
        -------

        Defines
        -------
        self.segm_img : ``photutils.segmentation.SegmentationImage``
            Two-dimensional segmentation image where found source regions are labeled with
            unique, non-zero positive integers.

        """
        # Recognize empty/blank images and treat gracefully
        if self.image.blank:
            self._define_empty_table(None)
            log.info('Total Detection Product - Segmentation source catalog is EMPTY')
            return

        # If the total product sources have not been identified, then this needs to be done!
        # Essentially, this "if" is the total product 'identify' block.  The filter product
        # inherits some variables from the total product in the "else" block.
        if not self.tp_sources:

            # Report configuration values to log
            log.info("{}".format("=" * 80))
            log.info("")
            log.info("SExtractor-like source finding settings - Photutils segmentation")
            log.info("Total Detection Product - Input Parameters")
            log.info("Image: {}".format(self.imgname))
            log.info("Gaussian FWHM (Filter for source detection in arcseconds): {}".format(self.param_dict['TWEAK_FWHMPSF']))
            log.info("size_source_box (no. of connected pixels needed for a detection): {}".format(self._size_source_box))
            log.info("Gaussian nsigma (threshold = nsigma * background_rms): {} (Round 1), {} (Round2 for Bkg2D)".format(self._nsigma, 2.0 * self._nsigma))
            log.info("nlevels (no. of multi-thresholding levels for deblending): {}".format(self._nlevels))
            log.info("contrast (frac. flux for peak to be separate object, 0=max. deblend, 1=no deblend): {}".format(self._contrast))
            log.info("RickerWavelet nsigma (threshold = nsigma * background_rms): {} (Round 1), {} (Round 2 for Bkg2D)".format(self._rw2d_nsigma, 2.0 * self._rw2d_nsigma))
            log.info("RickerWavelet kernel X- and Y-dimension: {}".format(self._rw2d_size))
            log.info("Pixel limit on biggest source (criterion for  RickerWavelet kernel): {}".format(self._rw2d_biggest_pixels))
            log.info("Percentage limit on biggest source (criterion for  RickerWavelet kernel): {}".format(100.0 * self._rw2d_biggest_source))
            log.info("Percentage limit on source fraction over the image (criterion for RickerWavelet kernel): {}".format(100.0 * self._rw2d_source_fraction))
            log.info("Percentage limit on biggest source deblending limit: {}".format(100.0 * self._bs_deblend_limit))
            log.info("Percentage limit on source fraction deblending limit: {}".format(100.0 * self._sf_deblend_limit))
            log.info("Scaling parameter of the Kron radius: {}".format(self._kron_scaling_radius))
            log.info("Kron minimum circular radius (pixels): {}".format(self._kron_minimum_radius))
            log.info("")
            log.info("{}".format("=" * 80))

            # Get the SCI image data
            imgarr = copy.deepcopy(self.image.data)
            log.debug("IMG stats: min/max: {}, {}".format(imgarr.min(), imgarr.max()))
            # Gaussian kernel
            g2d_kernel = self.image.kernel

            # Write out diagnostic data
            if self.diagnostic_mode:
                # Exclusion mask
                outname = self.imgname.replace(".fits", "_mask.fits")
                fits.PrimaryHDU(data=self.image.inv_footprint_mask.astype(np.uint16)).writeto(outname)

                # Background image
                outname = self.imgname.replace(".fits", "_bkg.fits")
                fits.PrimaryHDU(data=self.image.bkg_background_ra).writeto(outname)

                # filter kernel as well
                outname = self.imgname.replace(".fits", "_kernel.fits")
                fits.PrimaryHDU(data=g2d_kernel).writeto(outname)

            # Detect segments and evaluate the detection in terms of big sources/islands or crowded fields
            # Round 1
            ncount = 0
            log.info("")
            log.info("ROUND 1")
            log.info("Using Gaussian kernel to generate a segmentation map.")
            g_segm_img, g_is_big_crowded, g_bs, g_sf = self.detect_and_eval_segments(imgarr,
                                                                                     g2d_kernel,
                                                                                     ncount,
                                                                                     self._size_source_box,
                                                                                     self._nsigma,
                                                                                     self.image.bkg_background_ra,
                                                                                     self.image.bkg_rms_ra,
                                                                                     check_big_island_only=False,
                                                                                     rw2d_biggest_pixels=self._rw2d_biggest_pixels,
                                                                                     rw2d_biggest_source=self._rw2d_biggest_source,
                                                                                     rw2d_source_fraction=self._rw2d_source_fraction)
            segm_img_orig = copy.deepcopy(g_segm_img)

            # If the science field via the segmentation map is deemed crowded or has big sources/islands, compute the
            # RickerWavelet2DKernel and call detect_and_eval_segments() again.
            if g_is_big_crowded and g_segm_img:
                log.info("")
                log.info("The segmentation map contains big sources/islands or a large source fraction of segments.")
                log.info("Using RickerWavelet2DKernel to generate an alternate segmentation map.")

                # Convert the FWHM into a sigma value
                rw_sigma = self.image.kernel_fwhm * gaussian_fwhm_to_sigma
                rw2dk = self.ricker_matched_kernel(rw_sigma,
                                                   x_size=self._rw2d_size,
                                                   y_size=self._rw2d_size)

                # Only pass along the array to be consistent with the g2d_kernel object
                rw2d_kernel = rw2dk.array

                # Write out diagnostic data
                if self.diagnostic_mode:
                    outname = self.imgname.replace(".fits", "_rwkernel.fits")
                    fits.PrimaryHDU(data=rw2d_kernel).writeto(outname)

                # Detect segments and evaluate the detection in terms of big sources/islands or crowded fields
                # Round 1
                ncount += 1
                rw_segm_img, rw_is_big_crowded, rw_bs, rw_sf = self.detect_and_eval_segments(imgarr,
                                                                                             rw2d_kernel,
                                                                                             ncount,
                                                                                             self._size_source_box,
                                                                                             self._rw2d_nsigma,
                                                                                             self.image.bkg_background_ra,
                                                                                             self.image.bkg_rms_ra,
                                                                                             check_big_island_only=True,
                                                                                             rw2d_biggest_pixels=self._rw2d_biggest_pixels,
                                                                                             rw2d_biggest_source=self._rw2d_biggest_source,
                                                                                             rw2d_source_fraction=self._rw2d_source_fraction)

                # Check if the RickerWavelet segmentation image still seems to be problematic
                if rw_is_big_crowded and rw_segm_img:

                    # Before giving up, check the type of background computed for the detection image,
                    # and proceed based upon the type.  If a "sigma-clipped background" is in use, compute
                    # a "2D background" instead.  If a "2D background" is in use, increase the
                    # threshold for source detection.
                    log.info("")
                    log.info("RickerWavelet computed segmentation image still contains big sources/islands.")
                    log.info("Recomputing the threshold or background image for improved segmentation detection.")

                    # Make sure to be working with the unmodified image data
                    imgarr = copy.deepcopy(self.image.data)

                    # Background types: zero_background, sigma_clipped_background, twod_background
                    # Compute a twod_background
                    if (self.image.bkg_type.lower().startswith('sigma')):

                        log.info("Recomputing the background image from a sigma-clipped background to a Background2D.")

                        # In order to force the use of a background2D, some configuration values will be
                        # re-set (i.e., bkg_skew_threshold and negative_percent).
                        self.image.compute_background(self.param_dict['bkg_box_size'],
                                                      self.param_dict['bkg_filter_size'],
                                                      bkg_skew_threshold=0.0,
                                                      negative_percent=100.0)

                        if self.diagnostic_mode:
                            outname = self.imgname.replace(".fits", "_bkg1.fits")
                            fits.PrimaryHDU(data=self.image.bkg_background_ra).writeto(outname)

                        # Need to remake image kernel as it has a dependence on self.bkg_rms_ra
                        self.image.build_kernel(self.param_dict['bkg_box_size'],
                                                self.param_dict['bkg_filter_size'],
                                                self.param_dict['TWEAK_FWHMPSF'])

                        # Reset the local version of the Gaussian kernel and the RickerWavelet
                        # kernel when the background type changes
                        g2d_kernel = self.image.kernel

                        rw_sigma = self.image.kernel_fwhm * gaussian_fwhm_to_sigma
                        rw2dk = self.ricker_matched_kernel(rw_sigma,
                                                           x_size=self._rw2d_size,
                                                           y_size=self._rw2d_size)

                        rw2d_kernel = rw2dk.array
                        sigma_for_threshold = self._nsigma
                        rw2d_sigma_for_threshold = self._rw2d_nsigma

                    # Re-compute a background2D with a higher threshold by increasing the nsigma used
                    elif (self.image.bkg_type.lower().startswith('twod')):
                        log.info("Increasing the threshold by nsigma * 2.0 for improved source detection.")
                        sigma_for_threshold = self._nsigma * 2.0
                        rw2d_sigma_for_threshold = self._rw2d_nsigma * 2.0

                    # Define thresholds for empty/zero background
                    else:
                        log.info("Defining the threshold image based on nsigma for source detection.")
                        sigma_for_threshold = self._nsigma
                        rw2d_sigma_for_threshold = self._rw2d_nsigma

                    # Detect segments and evaluate the detection in terms of big sources/islands or crowded fields
                    # Round 2
                    ncount += 1
                    log.info("")
                    log.info("ROUND 2")
                    log.info("With alternate background...using Gaussian kernel to generate a segmentation map.")
                    del g_segm_img
                    g_segm_img, g_is_big_crowded, g_bs, g_sf = self.detect_and_eval_segments(imgarr,
                                                                                             g2d_kernel,
                                                                                             ncount,
                                                                                             self._size_source_box,
                                                                                             sigma_for_threshold,
                                                                                             self.image.bkg_background_ra,
                                                                                             self.image.bkg_rms_ra,
                                                                                             check_big_island_only=False,
                                                                                             rw2d_biggest_pixels=self._rw2d_biggest_pixels,
                                                                                             rw2d_biggest_source=self._rw2d_biggest_source,
                                                                                             rw2d_source_fraction=self._rw2d_source_fraction)

                    # Check again for big sources/islands or a large source fraction
                    if g_is_big_crowded:
                        log.info("")
                        log.info("The segmentation map contains big sources/islands or a large source fraction of segments.")
                        log.info("With alternate background...using RickerWavelet2DKernel to generate an alternate segmentation map.")
                        log.info("Checking of RickerWavelet2DKernel results using more generous big sources and source fraction limits.")

                        # Detect segments and evaluate the detection in terms of big sources/islands or crowded fields
                        # Note the biggest source and source fraction limits are the much larger "deblend" values.
                        # Round 2
                        ncount += 1
                        del rw_segm_img
                        rw_segm_img, rw_is_big_crowded, rw_bs, rw_sf = self.detect_and_eval_segments(imgarr,
                                                                                                     rw2d_kernel,
                                                                                                     ncount,
                                                                                                     self._size_source_box,
                                                                                                     rw2d_sigma_for_threshold,
                                                                                                     self.image.bkg_background_ra,
                                                                                                     self.image.bkg_rms_ra,
                                                                                                     check_big_island_only=True,
                                                                                                     rw2d_biggest_pixels=self._rw2d_biggest_pixels,
                                                                                                     rw2d_biggest_source=self._bs_deblend_limit,
                                                                                                     rw2d_source_fraction=self._sf_deblend_limit)

                        # Compute the ratio of big sources/islands using Gaussian vs Rickerwavelet kernel
                        # This value used as a discriminant between overlapping point sources and nebulousity fields
                        ratio_cg2rw_bigsource = 3.0
                        if rw_bs > 0.0:
                            ratio_cg2rw_bigsource = g_bs / rw_bs

                        # Last chance - The larger "deblend" limits were used in this last detection
                        # attempt based upon the the statistics of processing lots of data - looking
                        # for a balance between not being able to generate segmentation catalogs versus
                        # deblending for an unreasonable amount of time (days).
                        #
                        # Also, the ratio_cg2rw_bigsource is indicative of overlapping PSFs versus large
                        # areas of nebulousity. If this ratio is approximately > 2, then deblending can be
                        # quite efficient and successful for the overlapping PSF case.
                        #
                        # Use the Round 2 RickerWavelet segmentation image
                        log.info("Gaussian biggest source found: {}  Rickerwavelet biggest source found: {}.".format(g_bs, rw_bs))
                        log.info("Ratio of big sources found using Gaussian vs Rickerwavelet kernel: {}.".format(ratio_cg2rw_bigsource))

                        # This is the easy case.
                        if not rw_is_big_crowded:
                            self.kernel = rw2d_kernel
                            segm_img = copy.deepcopy(rw_segm_img)
                            del rw_segm_img
                            self.final_nsigma = rw2d_sigma_for_threshold
                        # The field was found to be crowded, but the biggest source second criterion is deemed OK.
                        elif (rw_is_big_crowded and (ratio_cg2rw_bigsource > self._ratio_bigsource_limit)):
                            log.info("The Round 2 of segmentation images may still contain big sources/islands.\n"
                                     "However, the ratio between the Gaussian and Rickerwavelet biggest source is\n"
                                     "indicative of overlapping PSFs vs nebulousity.")
                            log.info("Proceeding as the time to deblend should be nominal.")
                            self.kernel = rw2d_kernel
                            segm_img = copy.deepcopy(rw_segm_img)
                            del rw_segm_img
                            self.final_nsigma = rw2d_sigma_for_threshold
                        # The segmentation image is problematic and the big island/source fraction limits are exceeded,
                        # so deblending could take days, and the results would not be viable in any case.
                        else:
                            log.warning("")
                            log.warning("The Round 2 of segmentation images still contain big sources/islands or a\n"
                                     "large source fraction of segments.")
                            log.warning("The segmentation algorithm is unable to continue and no segmentation catalog will be produced.")
                            self._define_empty_table(rw_segm_img)
                            del g_segm_img
                            del rw_segm_img
                            return

                    # Use the second round Gaussian segmentation image
                    else:
                        self.kernel = g2d_kernel
                        segm_img = copy.deepcopy(g_segm_img)
                        del g_segm_img
                        self.final_nsigma = sigma_for_threshold

                # The first round RickerWavelet segmentation image is good, continue with the processing
                elif not rw_is_big_crowded and rw_segm_img:
                    self.kernel = rw2d_kernel
                    segm_img = copy.deepcopy(rw_segm_img)
                    del rw_segm_img
                    del g_segm_img
                    self.final_nsigma = self._rw2d_nsigma

                # No segments were detected in the total data product - no further processing done for this TDP,
                # but processing of another TDP should proceed.
                elif not rw_segm_img:
                    self._define_empty_table(rw_segm_img)
                    return

            # The first round Gaussian segmentation image is good, continue with the processing
            elif not g_is_big_crowded and g_segm_img:
                self.kernel = g2d_kernel
                segm_img = copy.deepcopy(g_segm_img)
                del g_segm_img
                self.final_nsigma = self._nsigma

            # No segments were detected in the total data product - no further processing done for this TDP,
            # but processing of another TDP should proceed.
            elif not g_segm_img:
                self._define_empty_table(g_segm_img)
                return

            # Report the final nsigma value used to compute the threshold above which sources are detected
            log.info("*** Final nsigma used (threshold = nsigma * background_rms): {}".format(self.final_nsigma))

            # If appropriate, deblend the segmentation image. Otherwise, use the current segmentation image
            ncount += 1
            if segm_img.big_segments is not None:
                segm_img = self.deblend_segments(segm_img,
                                                 imgarr,
                                                 ncount,
                                                 filter_kernel=self.kernel,
                                                 source_box=self._size_source_box)

            # The total product catalog consists of at least the X/Y and RA/Dec coordinates for the detected
            # sources in the total drizzled image.  All the actual measurements are done on the filtered drizzled
            # images using the coordinates determined from the total drizzled image.  Measure the coordinates now.
            log.info("Identifying sources in total detection image.")
            self.segm_img = copy.deepcopy(segm_img)
            del segm_img
            self.source_cat = SourceCatalog(imgarr, self.segm_img, background=self.image.bkg_background_ra,
                                            convolved_data = self.convolved_img, wcs=self.image.imgwcs,
                                            kron_params=[self._kron_scaling_radius, self._kron_minimum_radius])


            enforce_icrs_compatibility(self.source_cat)
            # Convert source_cat which is a SourceCatalog to an Astropy Table - need the data in tabular
            # form to filter out bad rows and correspondingly bad segments before the filter images are processed.
            total_measurements_table = Table(self.source_cat.to_table(columns=['label', 'xcentroid', 'ycentroid', 'sky_centroid_icrs']))

            # Filter the table to eliminate nans or inf based on the coordinates, then remove these segments from
            # the segmentation image too
            good_rows = []
            good_segm_rows_by_label = []
            updated_table = None
            for i, old_row in enumerate(total_measurements_table):
                if np.isfinite(old_row["xcentroid"]):
                    good_rows.append(old_row)
                    good_segm_rows_by_label.append(total_measurements_table['label'][i])

            updated_table = Table(rows=good_rows, names=total_measurements_table.colnames)

            # Need to keep an updated copy of the total image SourceCatalog object for use when
            # making measurements in the filtered images
            self.total_source_cat = self.source_cat.get_labels(good_segm_rows_by_label)

            # Clean up the existing column names, format, and descriptions
            # At this point self.source_cat has been converted into an Astropy table
            self.source_cat = self._define_total_table(updated_table)

            # self.sources needs to be passed to a filter catalog object based on code in hapsequencer.py
            # (create_catalog_products()).  This is the way the independent catalogs of total and filter products
            # process the same segmentation image.
            # BEWARE: self.sources for "segmentation" is a SegmentationImage, but for "point" it is an Astropy table
            # Keep only the good segments from the image
            self.segm_img.keep_labels(good_segm_rows_by_label, relabel=True)

            # Make a deep copy of the total "white light" segmentation image
            self.sources = self.segm_img.copy()

            log.info("Done identifying sources in total detection image for the segmentation catalog.")
            log.info("")
            log.info("{}".format("=" * 80))
            log.info("")

        # This is the filter product section -  use sources and other information identified in
        # the total detection product previously generated
        else:
            self.sources = self.tp_sources['segment']['sources']                   # segmentation image
            self.kernel = self.tp_sources['segment']['kernel']                     # image smoothing kernel
            self.total_source_table = self.tp_sources['segment']['source_cat']     # source catalog as a table
            self.total_source_cat = self.tp_sources['segment']['total_source_cat'] # SourceCatalog datatype
            self.final_nsigma = self.tp_sources['segment']['final_nsigma']         # final sigma multiplicative factor used

        # For debugging purposes only, create a "regions" files to use for ds9 overlay of the segm_img.
        # Create the image regions file here in case there is a failure.  This diagnostic portion of the
        # code should only be invoked when working on the total object catalog (self.segm_img is defined).
        if self.diagnostic_mode and self.segm_img:
            # Copy out only the X and Y coordinates to a "diagnostic_mode table" and cast as an Astropy Table
            # so a scalar can be added to the centroid coordinates
            tbl = self.source_cat["X-Centroid", "Y-Centroid"]

            # Construct the diagnostic_mode output filename and write the regions file
            indx = self.sourcelist_filename.find("ecsv")
            outname = self.sourcelist_filename[0:indx-1] + "_all.reg"

            tbl["X-Centroid"].info.format = ".12f"
            tbl["Y-Centroid"].info.format = ".12f"

            # Add one to the X and Y table values to put the data onto a one-based system,
            # particularly for display with ds9
            tbl["X-Centroid"] = tbl["X-Centroid"] + 1
            tbl["Y-Centroid"] = tbl["Y-Centroid"] + 1
            tbl.write(outname, format="ascii.commented_header")

            log.info("Wrote region file '{}' containing {} sources".format(outname, len(tbl)))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def ricker_matched_kernel(self, sigma, sigma_factor=1.5, **kw):
        """Return normalized RickerWavelet2D kernel matched to Gaussian width

        The sigma is increased by a factor sigma_factor=1.5 to match the Gaussian
        core width, and the normalization matches the RickerWavelet peak to the Gaussian
        peak.  This kernel should be applied using the 'normalize_kernel=False' parameter
        to astropy.convolution.convolve().

        Routine contributed by Rick L. White.

        Parameters
        ----------
        sigma : float
            Sigma of the RickerWavelet kernel

        sigma_factor : float
            Factor to match the Gaussian core width

        kw : int
            x_size : Size in the x-direction 
            y_size : Size in the y-direction 

        Returns
        -------
        RickerWavelet2DKernel : ˜astropy.convolution.RickerWavelet2DKernel
            Normalized RickerWavelet2DKernel
        """

        rsigma = sigma_factor*sigma
        rnorm = 0.5 * sigma_factor**4 * sigma**2
        return rnorm*RickerWavelet2DKernel(rsigma, **kw)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def detect_and_eval_segments(self, imgarr, kernel, ncount, size_source_box, nsigma_above_bkg, background_img, background_rms,
                                 check_big_island_only=False, rw2d_biggest_pixels=25000, rw2d_biggest_source=0.015,
                                 rw2d_source_fraction=0.075):

            # Compute the threshold to use for source detection
            threshold = self.compute_threshold(nsigma_above_bkg, background_rms)

            # Write out diagnostic data
            if self.diagnostic_mode:
                outname = self.imgname.replace(".fits", "_threshold" + str(ncount) + ".fits")
                fits.PrimaryHDU(data=threshold).writeto(outname)

            # Generate the segmentation map by detecting "sources" using the nominal settings.
            segm_img = self.detect_segments(imgarr,
                                            threshold,
                                            background_img,
                                            ncount,
                                            filter_kernel=kernel,
                                            source_box=size_source_box,
                                            mask=self.image.inv_footprint_mask)

            # Check if segm_image is None indicating there are no detectable sources in this
            # total detection image.  If value is None, a warning has already been issued.  Issue
            # a final message for this particular total detection product and return.
            if segm_img is None:
                log.warning("End processing for the segmentation catalog due to no sources detected with the current kernel.")
                log.warning("An empty catalog will be produced for this total detection product, {}.".format(self.imgname))
                is_big_crowded = True
                big_island = 1.0
                source_fraction = 1.0

            else:
                # Determine if the segmentation image is filled with big sources/islands (bs) or is crowded with a large
                # source fraction (sf). Depending upon these measurements, it can take a very, very long time to deblend
                # the sources.
                is_big_crowded = True
                is_big_crowded, big_island, source_fraction = self._evaluate_segmentation_image(segm_img,
                                                                                                imgarr,
                                                                                                big_island_only=check_big_island_only,
                                                                                                max_biggest_pixels=rw2d_biggest_pixels,
                                                                                                max_biggest_source=rw2d_biggest_source,
                                                                                                max_source_fraction=rw2d_source_fraction)

            return segm_img, is_big_crowded, big_island, source_fraction

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def compute_threshold(self, nsigma, bkg_rms):
        """Compute the threshold value above which sources are deemed detected.

           Parameters
           ----------
           nsigma : float
               Multiplicative factor for the background RMS

           bkg_rms : float image
               RMS of the background determined image

           Returns
           -------
           threshold: float image
               Image which defines, on a pixel-by-pixel basis, the low limit above which
               sources are detected.

           Note: The background image is NOT built into the threshold image.
        """

        log.info("Computing the threshold value used for source detection.")
        if not self.tp_masks:
            threshold = nsigma * bkg_rms
        else:
            threshold = np.zeros_like(self.tp_masks[0]['rel_weight'])
            log.info("Using WHT masks as a scale on the RMS to compute threshold detection limit.")
            for wht_mask in self.tp_masks:
                threshold_rms = bkg_rms * np.sqrt(wht_mask['mask'] / wht_mask['rel_weight'].max())
                threshold_rms_median = np.nanmedian(threshold_rms[threshold_rms > 0])
                threshold_item = nsigma * threshold_rms_median

                # Ensure the threshold image is built up, section by section, as created by
                # the wht_mask['mask']
                threshold += (threshold_item * wht_mask['mask'])
            del(threshold_rms)

        return threshold

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def detect_segments(self, imgarr, threshold, background_img, ncount, filter_kernel=None, source_box=6, mask=None):
        """Detect segments found in the input total detection (aka white light) image.

           Image regions are identified as segments in the 'total detection image' if the region
           has n-connected pixels with values greater than the 'threshold'. The resultant
           segmentation image is then evaluated to determine if the segments represent a large portion
           of the image or if any one segment is very large.

           Parameters
           ----------
           imgarr : float
               Total detection image (no background subtraction)

           threshold : float
               Image which defines, on a pixel-by-pixel basis, the low limit above which
               sources are detected.  The threshold is composed of the 'sigma * rms'.

           background_img : float
               Computed 2D background for total detection image
               
           ncount : int
               Invocation index for this method.  The index is used to create unique names
               for diagnostic output files.

           filter_kernel : astropy.convolution.Kernel2D object, optional
               Filter used to smooth the total detection image to enhance peak or
               multi-scale detection.

           source_box : int, optional
               Number of connected pixels needed to define a source

           mask : bool image, optional
               Boolean image used to define the illuminated and non-illuminated pixels.

           Returns
           -------
           segm_img : `~photutils.segmentation.SegmentationImage` or None

        """
        log.info("Detecting sources in total image product.")

        # Subtract the background image from the white light image and convolve
        #
        # Important to background subtract to get accurate positions for faint sources or
        # positions will be shifted toward the geometrical center of the segment rather 
        # than located at the flux-weighted centroid as desired.
        self.convolved_img = convolve(imgarr-background_img, filter_kernel, normalize_kernel=False)

        # The threshold should not include the background.
        segm_img = detect_sources(self.convolved_img,
                                  threshold,
                                  npixels=source_box,
                                  mask=mask)

        # If no segments were found, there are no detectable sources in the total detection image.
        # Return as there is nothing more to do.
        if segm_img is None:
            log.warning("No segments were found in Total Detection Product, {}.".format(self.imgname))
            log.warning("Processing for segmentation source catalogs for this product is ending.")
            return segm_img

        if self.diagnostic_mode:
            outname = self.imgname.replace(".fits", "_segment" + str(ncount) + ".fits")
            fits.PrimaryHDU(data=segm_img.data).writeto(outname)

        return segm_img

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def deblend_segments(self, segm_img, imgarr, ncount, filter_kernel=None, source_box=6):
        """Deblend segments found in the input total detection (aka white light) image.

           The segmentation image generated by detect_segments is deblended in an effort to
           separate overlapping sources.

           Parameters
           ----------
           segm_img : `~photutils.segmentation.SegmentationImage`
               Segmentation image created by detect_segments() based on the total detection image

           imgarr :
               Total detection image

           ncount : int
               Invocation index for this method.  The index is used to create unique names
               for diagnostic output files.

           filter_kernel : astropy.convolution.Kernel2D object, optional
               Filter used to smooth the total detection image to enhance peak or
               multi-scale detection.

           source_box : int, optional
               Number of connected pixels needed to define a source

           Updates
           -------
           segm_img : `~photutils.segmentation.SegmentationImage`
               Deblended segmentation image
        """

        log.info("Deblending segments in total image product.")
        # Note: SExtractor has "connectivity=8" which is the default for detect_sources().
        # Initialize return value in case of failure in deblending
        segm_deblended_img = segm_img
        try:
            # Deblending is a combination of multi-thresholding and watershed
            # segmentation. Sextractor uses a multi-thresholding technique.
            # npixels = minimum number of connected pixels in source
            # npixels and filter_kernel should match those used by detect_sources()
            segm_deblended_img = deblend_sources(self.convolved_img,
                                                 segm_img,
                                                 npixels=source_box,
                                                 nlevels=self._nlevels,
                                                 contrast=self._contrast,
                                                 labels=segm_img.big_segments)

            if self.diagnostic_mode:
                log.info("Deblended {} out of {} segments".format(len(segm_img.big_segments), segm_img.nlabels))
                outname = self.imgname.replace(".fits", "_segment_deblended" + str(ncount) + ".fits")
                fits.PrimaryHDU(data=segm_deblended_img.data).writeto(outname)

        except Exception as x_cept:
            log.warning("Deblending the segments in image {} was not successful: {}.".format(self.imgname,
                        x_cept))
            log.warning("Processing can continue with the non-deblended segments, but the user should\n"
                        "check the output catalog for issues.")

        # The deblending was successful, so just return the deblended SegmentationImage to calling routine.
        return segm_deblended_img
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def measure_sources(self, filter_name):
        """Use the positions of the sources identified in the white light (total detection) image to
        measure properties of these sources in the filter images.

        An instrument/detector combination may have multiple filter-level products.
        This routine is called for each filter image which is then measured to generate
        a filter-level source catalog based on object positions measured in the total
        detection product image.

        Returns
        -------

        """
        if self.sources is None or self.sources.nlabels == 0:
            # Report configuration values to log
            log.info("{}".format("=" * 80))
            log.info("")
            log.info("No segmentation sources identified for photometry for")
            log.info("image name: {}".format(self.imgname))
            log.info("Generating empty segment catalog.")
            log.info("")
            self._define_empty_table(None)

            # Capture specified filter columns in order to append to the total detection table
            self.subset_filter_source_cat = Table(names=["ID", "MagAp2", "CI", "Flags"])
            self.subset_filter_source_cat.rename_column("MagAp2", "MagAP2_" + filter_name)
            self.subset_filter_source_cat.rename_column("CI", "CI_" + filter_name)
            self.subset_filter_source_cat.rename_column("Flags", "Flags_" + filter_name)

            return

        # Get filter-level science data
        imgarr = copy.deepcopy(self.image.data)

        # Report configuration values to log
        log.info("{}".format("=" * 80))
        log.info("")
        log.info("SExtractor-like source property measurements based on Photutils segmentation")
        log.info("Filter Level Product - Input Parameters")
        log.info("image name: {}".format(self.imgname))
        log.info("Gaussian FWHM (Filter for source detection in arcseconds): {}".format(self.param_dict['TWEAK_FWHMPSF']))
        log.info("size_source_box (pixels): {}".format(self._size_source_box))
        log.info("")

        # This is the filter science data and its computed background
        imgarr_bkgsub = imgarr - self.image.bkg_background_ra

        # Compute the Poisson error of the sources...
        total_error = calc_total_error(imgarr_bkgsub, self.image.bkg_rms_ra, 1.0)

        # Columns to include from the computation of source properties to save
        # computation time from computing values which are not used
        include_filter_cols = ['area', bac_colname, 'bbox_xmax', 'bbox_xmin', 'bbox_ymax', 'bbox_ymin',
                               'covar_sigx2', 'covar_sigxy', 'covar_sigy2', 'cxx', 'cxy', 'cyy',
                               'ellipticity', 'elongation', id_colname, 'kron_radius', 'orientation', 'sky_centroid_icrs',
                               flux_colname, ferr_colname, 'xcentroid', 'ycentroid']

        # Compute source properties...
        # No convolved image needed here since the detection_cat parameter is set.
        include_filter_cols.append('fwhm')
        self.source_cat = SourceCatalog(imgarr_bkgsub, self.sources, background=self.image.bkg_background_ra,
                                        detection_cat=self.total_source_cat,
                                        error=total_error, wcs=self.image.imgwcs)

        enforce_icrs_compatibility(self.source_cat)

        filter_measurements_table = Table(self.source_cat.to_table(columns=include_filter_cols))

        # Compute aperture photometry measurements and append the columns to the measurements table
        self.do_aperture_photometry(imgarr, filter_measurements_table, self.imgname, filter_name)

        # Now clean up and prepare the filter table for output
        self.source_cat = self._define_filter_table(filter_measurements_table)

        log.info("Found and measured {} sources from segmentation map.".format(len(self.source_cat)))
        log.info("")
        log.info("{}".format("=" * 80))
        log.info("")

        # Capture specified filter columns in order to append to the total detection table
        self.subset_filter_source_cat = self.source_cat["ID", "MagAp2", "CI", "Flags"]
        self.subset_filter_source_cat.rename_column("MagAp2", "MagAP2_" + filter_name)
        self.subset_filter_source_cat.rename_column("CI", "CI_" + filter_name)
        self.subset_filter_source_cat.rename_column("Flags", "Flags_" + filter_name)

        if self.diagnostic_mode:
            # Write out a catalog which can be used as an overlay for image in ds9
            # The source coordinates are the same for the total and filter products, but the kernel is
            # total- or filter-specific and any row with nan or inf has been removed from the filter table..

            # Copy out only the X and Y coordinates to a "diagnostic_mode table" and
            # cast as an Astropy Table so a scalar can be added later
            tbl = Table(self.source_cat["X-Centroid", "Y-Centroid"])

            # Construct the diagnostic_mode output filename and write the catalog
            indx = self.sourcelist_filename.find("ecsv")
            outname = self.sourcelist_filename[0:indx-1] + "_all.reg"

            tbl["X-Centroid"].info.format = ".10f"  # optional format
            tbl["Y-Centroid"].info.format = ".10f"

            # Add one to the X and Y table values to put the data onto a one-based system,
            # particularly for display with ds9
            tbl["X-Centroid"] = tbl["X-Centroid"] + 1
            tbl["Y-Centroid"] = tbl["Y-Centroid"] + 1
            tbl.write(outname, format="ascii.commented_header")
            log.info("Wrote the diagnostic_mode version of the filter detection source catalog: {}\n".format(outname))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def do_aperture_photometry(self, input_image, filter_measurements_table, image_name, filter_name):
        """Perform aperture photometry measurements as a means to distinguish point versus extended sources.
        """

        # Convert the SkyCoord column to separate RA and Dec columns
        radec_data = SkyCoord(filter_measurements_table["sky_centroid_icrs"])
        ra_icrs = radec_data.ra.degree
        dec_icrs = radec_data.dec.degree
        rr = Column(ra_icrs, name="RA", unit="degree")
        dd = Column(dec_icrs, name="DEC", unit="degree")
        filter_measurements_table.add_columns([dd, rr])

        # Compute the MagSegment
        filter_measurements_table["MagSegment"] = photometry_tools.convert_flux_to_abmag(filter_measurements_table[flux_colname],
                                                                                         self.image.keyword_dict['photflam'],
                                                                                         self.image.keyword_dict['photplam'])

        # Determine the "good rows" as defined by the X and Y coordinates not being nans as
        # the pos_xy array cannot contain any nan values.  Note: It is possible for a filter
        # catalog to have NO sources (subarray data).  The "bad rows" will have the RA and DEC values
        # in the filter catalog replaced with the corresponding RA and Dec values from the total
        # catalog to give the user some perspective.
        good_rows_index = []
        bad_rows_index = []
        for i, old_row in enumerate(filter_measurements_table):
            if np.isfinite(old_row["xcentroid"]):
                good_rows_index.append(i)
            else:
                bad_rows_index.append(i)

        # Create placeholder columns for the output table
        self._create_table_columns(filter_measurements_table)

        # Case: there are good/measurable sources in the input table
        if good_rows_index:
            # Obtain the X and Y positions to compute the circular annulus
            positions = (filter_measurements_table["xcentroid"][good_rows_index], filter_measurements_table["ycentroid"][good_rows_index])
            pos_xy = np.vstack(positions).T

            # Define list of background annulii
            bg_apers = CircularAnnulus(pos_xy,
                                       r_in=self.param_dict['skyannulus_arcsec']/self.image.imgwcs.pscale,
                                       r_out=(self.param_dict['skyannulus_arcsec'] + self.param_dict['dskyannulus_arcsec'])/self.image.imgwcs.pscale)

            # Create list of photometric apertures to measure
            phot_apers = [CircularAperture(pos_xy, r=r) for r in self.aper_radius_list_pixels]

            # Perform aperture photometry - the input data should NOT be background subtracted
            photometry_tbl = photometry_tools.iraf_style_photometry(phot_apers,
                                                                    bg_apers,
                                                                    data=input_image,
                                                                    photflam=self.image.keyword_dict['photflam'],
                                                                    photplam=self.image.keyword_dict['photplam'],
                                                                    error_array=self.image.bkg_rms_ra,
                                                                    bg_method=self.param_dict['salgorithm'],
                                                                    epadu=self.gain)

            # Capture data computed by the photometry tools and append to the output table
            filter_measurements_table['FluxAp1'][good_rows_index] = photometry_tbl['FluxAp1']
            filter_measurements_table['FluxErrAp1'][good_rows_index] = photometry_tbl['FluxErrAp1']
            filter_measurements_table['MagAp1'][good_rows_index] = photometry_tbl['MagAp1']
            filter_measurements_table['MagErrAp1'][good_rows_index] = photometry_tbl['MagErrAp1']

            filter_measurements_table['FluxAp2'][good_rows_index] = photometry_tbl['FluxAp2']
            filter_measurements_table['FluxErrAp2'][good_rows_index] = photometry_tbl['FluxErrAp2']
            filter_measurements_table['MagAp2'][good_rows_index] = photometry_tbl['MagAp2']
            filter_measurements_table['MagErrAp2'][good_rows_index] = photometry_tbl['MagErrAp2']

            filter_measurements_table['MSkyAp2'][good_rows_index] = photometry_tbl['MSkyAp2']
            filter_measurements_table['StdevAp2'][good_rows_index] = photometry_tbl['StdevAp2']

            mag_inner_data = photometry_tbl["MagAp1"].data
            mag_outer_data = photometry_tbl["MagAp2"].data

            try:
                # Compute the Concentration Index (CI)
                if math.isclose(self.image.keyword_dict['photflam'], 0.0, abs_tol=constants.TOLERANCE):
                    # Fill the ci_data column with the "bad data indicator" of -9999.0
                    # where photometry_tbl["MagAp1"] was populated with -9999.0 in iraf_style_photometry
                    ci_data = photometry_tbl["MagAp1"].data
                    ci_mask = np.logical_not(ci_data > constants.FLAG + 1.0)
                else:
                    ci_data = mag_inner_data - mag_outer_data
                    ci_mask = np.logical_and(np.abs(ci_data) > 0.0, np.abs(ci_data) < 1.0e-30)
                    big_bad_index = np.where(abs(ci_data) > 1.0e20)
                    ci_mask[big_bad_index] = True

            except Exception as x_cept:
                log.warning("Computation of concentration index (CI) was not successful: {} - {}.".format(self.imgname, x_cept))
                log.warning("CI measurements may be missing from the output filter catalog.\n")
                # Create a column of data from an existing column of data
                # and set each value = -9999.0
                ci_data = constants.FLAG + photometry_tbl["MagAp1"] * 0.0
                ci_mask = np.logical_not(ci_data > constants.FLAG + 1.0)

            ci_col = MaskedColumn(name='CI', data=ci_data, dtype=np.float64, mask=ci_mask)

            # OK to insert *entire* column here to preserve any values which have been computed.  The
            # column already exists and contains nans.
            if isinstance(ci_col, MaskedColumn):
                filter_measurements_table['CI'][good_rows_index] = ci_col

        # Case: no good rows in the table
        # Issue a message for the case no good rows at all being found in the filter image.
        # The bad rows will be "filled in" with the same code (below) where all the rows are bad.
        else:
            log.info("There are no valid rows in the output Segmentation filter catalog for image %s (filter: %s).", image_name, filter_name)

        # Fill in any bad rows - this code fills in sporadic missing rows, as well as all rows being missing.
        # The bad rows have nan values for the xcentroid/ycentroid coordinates, as well as the RA/Dec values,
        # so recover these values from the total source catalog to make it easy for the user to map the filter
        # catalog rows back to the total detection catalog
        filter_measurements_table['xcentroid'][bad_rows_index] = self.total_source_table['X-Centroid'][bad_rows_index]
        filter_measurements_table['ycentroid'][bad_rows_index] = self.total_source_table['Y-Centroid'][bad_rows_index]
        filter_measurements_table['RA'][bad_rows_index] = self.total_source_table['RA'][bad_rows_index]
        filter_measurements_table['DEC'][bad_rows_index] = self.total_source_table['DEC'][bad_rows_index]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def _create_table_columns(self, table):
        """Create placeholder columns for the output filter table.

           The output filter table becomes the filter catalog ECSV file.

           Define the column order, data format, output column names, descriptions, and units
           for the table.

           Parameters
           ----------
           table : Astropy table

           Returns
           -------
           table : Astropy table
               A modified version of the input table which now has additional placeholder
               columns appended.
        """

        tblLen = len(table)
        ci_col = MaskedColumn(data=np.ones(tblLen)*float('nan'), name="CI")
        table.add_column(ci_col)

        flux_col = Column(data=np.ones(tblLen)*float('nan'), name="FluxAp1")
        table.add_column(flux_col)

        flux_col_err = Column(data=np.ones(tblLen)*float('nan'), name="FluxErrAp1")
        table.add_column(flux_col_err)

        mag_col = Column(data=np.ones(tblLen)*float('nan'), name="MagAp1")
        table.add_column(mag_col)

        mag_col_err = Column(data=np.ones(tblLen)*float('nan'), name="MagErrAp1")
        table.add_column(mag_col_err)

        flux_col = Column(data=np.ones(tblLen)*float('nan'), name="FluxAp2")
        table.add_column(flux_col)

        flux_col_err = Column(data=np.ones(tblLen)*float('nan'), name="FluxErrAp2")
        table.add_column(flux_col_err)

        mag_col = Column(data=np.ones(tblLen)*float('nan'), name="MagAp2")
        table.add_column(mag_col)

        mag_col_err = Column(data=np.ones(tblLen)*float('nan'), name="MagErrAp2")
        table.add_column(mag_col_err)

        msky_col = Column(data=np.ones(tblLen)*float('nan'), name="MSkyAp2")
        table.add_column(msky_col)

        stdev_col = Column(data=np.ones(tblLen)*float('nan'), name="StdevAp2")
        table.add_column(stdev_col)

        # Add zero-value "Flags" column in preparation for source flagging
        flag_col = Column(name="Flags", data=np.zeros_like(table[id_colname]))
        table.add_column(flag_col)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def _define_filter_table(self, filter_table):
        """Set the overall format for the filter output catalog.

           Define the column order, data format, output column names, descriptions, and units
           for the table.

           Parameters
           ----------
           filter_table : Astropy table
               Table which has been generated based upon SourceCatalog information and contains
               many properties calculated for each segmented source.

           Returns
           -------
           final_filter_table : Astropy table
                A modified version of the input table which has been reformatted in preparation
                for catalog generation.
        """

        # Rename columns to names used when HLA Classic catalog distributed by MAST
        final_col_names = {id_colname: "ID", "xcentroid": "X-Centroid", "ycentroid": "Y-Centroid",
                           bac_colname: "Bck", flux_colname: "FluxSegment", ferr_colname: "FluxSegmentErr",
                           "bbox_xmin": "Xmin", "bbox_ymin": "Ymin", "bbox_xmax": "Xmax", "bbox_ymax": "Ymax",
                           "cxx": "CXX", "cyy": "CYY", "cxy": "CXY",
                           "covar_sigx2": "X2", "covar_sigy2": "Y2", "covar_sigxy": "XY",
                           "orientation": "Theta",
                           "elongation": "Elongation", "ellipticity": "Ellipticity", "area": "Area",
                           "fwhm": "FWHM", "kron_radius": "KronRadius"}
        for old_col_title in final_col_names:
            filter_table.rename_column(old_col_title, final_col_names[old_col_title])

        # Define the order of the columns
        final_col_order = ["X-Centroid", "Y-Centroid", "RA", "DEC", "ID",
                           "CI", "Flags", "MagAp1", "MagErrAp1", "FluxAp1", "FluxErrAp1",
                           "MagAp2", "MagErrAp2", "FluxAp2", "FluxErrAp2", "MSkyAp2",
                           "Bck", "Area", "FWHM", "MagSegment", "FluxSegment", "FluxSegmentErr",
                           "KronRadius", "Xmin", "Ymin", "Xmax", "Ymax",
                           "X2", "Y2", "XY",
                           "CXX", "CYY", "CXY",
                           "Elongation", "Ellipticity", "Theta"]
        final_filter_table = filter_table[final_col_order]

        # Define the format
        final_col_format = {"X-Centroid": "10.3f", "Y-Centroid": "10.3f", "RA": "13.7f", "DEC": "13.7f", "ID": "7d",
                            "CI": "7.3f", "Flags": "5d",
                            "MagAp1": "7.3f", "MagErrAp1": "7.3f", "FluxAp1": "10.4f", "FluxErrAp1": "10.4f",
                            "MagAp2": "7.3f", "MagErrAp2": "7.3f", "FluxAp2": "10.4f", "FluxErrAp2": "10.4f",
                            "MSkyAp2": "7.3f", "Bck": "10.4f", "MagSegment": "7.3f", "FluxSegment": "10.4f",
                            "FluxSegmentErr": "10.4f",
                            "KronRadius": "8.4f", "Xmin": "8.0f", "Ymin": "8.0f", "Xmax": "8.0f", "Ymax": "8.0f",
                            "X2": "8.4f", "Y2": "8.4f", "XY": "8.4f",
                            "CXX": "9.5f", "CYY": "9.5f", "CXY": "9.5f",
                            "Elongation": "7.2f", "Ellipticity": "7.2f", "Theta": "8.3f", "Area": "8.3f", "FWHM": "8.3f"}
        for fcf_key in final_col_format.keys():
            final_filter_table[fcf_key].format = final_col_format[fcf_key]

        # Add description
        final_col_descrip = {"ID": "Catalog Object Identification Number",
                             "X-Centroid": "Pixel Coordinate", "Y-Centroid": "Pixel Coordinate",
                             "RA": "Sky coordinate at epoch of observation",
                             "DEC": "Sky coordinate at epoch of observation",
                             "Bck": "Background at the position of the source centroid",
                             "Area": "Total unmasked area of the source segment",
                             "FWHM": "Circularized FWHM of 2D Gaussian function having the same second-order central moments as the source",
                             "MagAp1": "ABMAG of source based on the inner (smaller) aperture",
                             "MagErrAp1": "Error of MagAp1",
                             "FluxAp1": "Flux of source based on the inner (smaller) aperture",
                             "FluxErrAp1": "Error of FluxAp1",
                             "MagAp2": "ABMAG of source based on the outer (larger) aperture",
                             "MagErrAp2": "Error of MagAp2",
                             "FluxAp2": "Flux of source based on the outer (larger) aperture",
                             "FluxErrAp2": "Error of FluxAp2",
                             "MSkyAp2": "Sky estimate from an annulus outside Aperture 2",
                             "FluxSegment": "Sum of unmasked data values in the source segment",
                             "FluxSegmentErr": "Uncertainty of FluxSegment, propagated from the input error array",
                             "KronRadius": "The unscaled first-moment Kron radius",
                             "MagSegment": "Magnitude corresponding to FluxSegment",
                             "X2": "Variance along X",
                             "Y2": "Variance along Y",
                             "XY": "Covariance of position between X and Y",
                             "CXX": "SExtractor's ellipse parameter", "CYY": "SExtractor's ellipse parameter",
                             "CXY": "SExtractor's ellipse parameter",
                             "Xmin": "Minimum X pixel within the minimal bounding box containing the source segment",
                             "Xmax": "Maximum X pixel within the minimal bounding box containing the source segment",
                             "Ymin": "Minimum Y pixel within the minimal bounding box containing the source segment",
                             "Ymax": "Maximum Y pixel within the minimal bounding box containing the source segment",
                             "Elongation": "Ratio of the lengths of the semi-major and semi-minor axes of the ellipse",
                             "Ellipticity": "Computed as 1 minus the inverse of the elongation",
                             "Theta": "Angle between the X axis and the major axis of the 2D Gaussian function that has the same second-order moments as the source.",
                             "CI": "Concentration Index",
                             "Flags": "Numeric encoding for conditions on detected sources"}
        for fcd_key in final_col_descrip.keys():
            final_filter_table[fcd_key].description = final_col_descrip[fcd_key]

        # Add units
        final_col_unit = {"X-Centroid": "pixel", "Y-Centroid": "pixel",
                          "RA": "degree", "DEC": "degree",
                          "Bck": "electron/s",
                          "Area": "pixel**2",
                          "FWHM": "pixel",
                          "MagAp1": "mag(AB)",
                          "MagErrAp1": "mag(AB)",
                          "FluxAp1": "electron/s",
                          "FluxErrAp1": "electron/s",
                          "MagAp2": "mag(AB)",
                          "MagErrAp2": "mag(AB)",
                          "FluxAp2": "electron/s",
                          "FluxErrAp2": "electron/s",
                          "MSkyAp2": "electron/(s pixel)",
                          "MagSegment": "mag(AB)",
                          "FluxSegment": "electron/s",
                          "FluxSegmentErr": "electron/s",
                          "KronRadius": "pixel",
                          "X2": "pixel**2",
                          "Y2": "pixel**2",
                          "XY": "pixel**2",
                          "CXX": "pixel**(-2)",
                          "CYY": "pixel**(-2)",
                          "CXY": "pixel**(-2)",
                          "Xmin": "pixel",
                          "Ymin": "pixel",
                          "Xmax": "pixel",
                          "Ymax": "pixel",
                          "Theta": "radian",
                          "CI": "mag(AB)",
                          "Flags": "",
                          "ID": "",
                          "Elongation": "",
                          "Ellipticity": ""}
        for fcu_key in final_col_unit.keys():
            final_filter_table[fcu_key].unit = final_col_unit[fcu_key]

        return(final_filter_table)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def _define_empty_table(self, segm_img):
        """Create basic empty table based on total_table format to signify no valid sources were found"""

        final_col_unit = {"X-Centroid": "pixel", "Y-Centroid": "pixel",
                          "RA": "degree", "DEC": "degree", "Flags": ""}
        final_col_format = {"ID": "7d",
                            "X-Centroid": "10.3f", "Y-Centroid": "10.3f",
                            "RA": "13.7f", "DEC": "13.7f",
                            "Flags": "5d"}
        final_col_descrip = {"ID": "Catalog Object Identification Number",
                             "X-Centroid": "Pixel Coordinate",
                             "Y-Centroid": "Pixel Coordinate",
                             "RA": "Sky coordinate at epoch of observation",
                             "DEC": "Sky coordinate at epoch of observation",
                             "Flags": "Numeric encoding for conditions on detected sources"}

        final_colnames = [k for k in final_col_format.keys()]
        # Initialize empty table with desired column names, descriptions and units
        empty_table = Table(names=final_colnames,
                      descriptions=final_col_descrip,
                      units=final_col_unit)

        # Add formatting for each column
        for fcf_key in final_col_format.keys():
            empty_table[fcf_key].format = final_col_format[fcf_key]

        self.source_cat = empty_table
        self.sources = copy.deepcopy(segm_img)
        if self.sources:
            self.sources.nlabels = 0  # Insure nlabels is set to 0 to indicate no valid sources


    def _define_total_table(self, updated_table):
        """Set the overall format for the total detection output catalog.

           Define the column order, data format, output column names, descriptions, and units
           for the table.

           Parameters
           ----------
           updated_table : Astropy table
               Table which has been generated based upon SourceCatalog information and contains
               many properties calculated for each segmented source.

           Returns
           ------
           table : Astropy table
                A modified version of the input table which has been reformatted in preparation
                for catalog generation.
        """

        # Extract just a few columns generated by the SourceCatalog as
        # more columns are appended to this table from the filter results.
        # Actually, the filter columns are in a table which is "database joined"
        # to the total table.  During the combine process, the new columns are renamed,
        # formatted, and described (as necessary). For now this table only has id, xcentroid,
        # ycentroid, RA, and DEC.
        table = updated_table["label", "xcentroid", "ycentroid"]

        # Convert the RA/Dec SkyCoord into separate columns
        radec_data = SkyCoord(updated_table["sky_centroid_icrs"])
        ra_icrs = radec_data.ra.degree
        dec_icrs = radec_data.dec.degree
        rr = Column(ra_icrs, name="RA", unit="degree")
        dd = Column(dec_icrs, name="DEC", unit="degree")
        table.add_columns([rr, dd])

        # Rename columns to names to those used when HLA Classic catalog distributed by MAST
        # and/or to distinguish Point and Segment catalogs
        # The columns that are appended will be renamed during the combine process
        final_col_names = {"label": "ID", "xcentroid": "X-Centroid", "ycentroid": "Y-Centroid"}
        for old_col_title in final_col_names:
            table.rename_column(old_col_title, final_col_names[old_col_title])

        # Format the current columns
        final_col_format = {"ID": "7d", "X-Centroid": "10.3f", "Y-Centroid": "10.3f", "RA": "13.7f", "DEC": "13.7f"}
        for fcf_key in final_col_format.keys():
            table[fcf_key].format = final_col_format[fcf_key]

        # Add description
        descr_str = "Sky coordinate at epoch of observation"
        final_col_descrip = {"ID": "Catalog Object Identification Number",
                             "X-Centroid": "Pixel Coordinate", "Y-Centroid": "Pixel Coordinate",
                             "RA": descr_str, "DEC": descr_str}
        for fcd_key in final_col_descrip.keys():
            table[fcd_key].description = final_col_descrip[fcd_key]

        # Add units
        final_col_unit = {"X-Centroid": "pixel", "Y-Centroid": "pixel",
                          "RA": "degree", "DEC": "degree"}
        for fcu_key in final_col_unit.keys():
            table[fcu_key].unit = final_col_unit[fcu_key]

        return(table)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def _evaluate_segmentation_image(self, segm_img, image_data, big_island_only=False,
                                     max_biggest_pixels=25000, max_biggest_source=0.015, max_source_fraction=0.075):
        """
        Determine if the largest "source" or if the total fraction of "sources" exceeds a threshold.

        Identify situations where the largest source exceeds a user-specified fraction of all image
        pixels and / or situations where the total source fraction exceeds a user-specified fraction
        of the image (aka 'big islands')

        Algorithm is essentially the standalone function written by R.L.White (STScI) for the
        Hubble Legacy Archive (HLA).

        Parameters
        ----------
        segm_img : Segmentation image
            Segmentation image created by the Photutils package by detect_sources ().

        image_data :  FITS data
            The total drizzled detection image (aka white light data).

        big_island_only : bool, optional
            Test for 'big island' situations only? (True/False)

        max_biggest_pixels : int, optional
            Maximum limit on the single largest detected "source" in pixels.

        max_biggest_source : float, optional
            Maximum limit on the single largest detected "source".

        max_source_fraction : float, optional
            Maximum limit on the fraction of pixels identified as part of a "source".

        Returns
        -------
        is_poor_quality: boolean
            True/False value indicating if the largest source or the total combination of
            all detected sources took up an abnormally high portion of the image (aka 'big islands').

        biggest_source : float
            Value of the largest source in comparison to the total number of illuminated
            pixels in the image.

        source_fraction : float
            Value of the total combination of all detected sources in comparison to the total
            number of illuminated pixels in the image.
        """

        log.info("")
        log.info("Analyzing segmentation image.")

        is_poor_quality = True
        biggest_source = 1.0
        source_fraction = 1.0

        if segm_img is None:
            log.info("Segmentation image is blank.")
            return is_poor_quality, biggest_source, source_fraction

        # If the segmentation image is not blank, start out assuming it is good.
        is_poor_quality = False

        # segm_img is a SegmentationImage, nbins must be at least 1 or segm_img == None
        nbins = segm_img.max_label
        log.info("Number of sources from segmentation map: %d", nbins)

        # Compute the biggest source identified 
        # The clip eliminates zero-divides when there are no good pixels in the image
        real_pixels = (np.isfinite(image_data) & (image_data != 0)).sum().clip(min=1)
        biggest_source_pixels = segm_img.areas.max()
        biggest_source = biggest_source_pixels/real_pixels

        # Compute which segments are larger than the kernel.
        deb_limit = self.kernel.size
        log.debug("Deblending limit set at: {} pixels".format(deb_limit))

        # Add big_segments as an attribute to SegmentationImage for use when deblending.
        # This is a Numpy array which currently contains the segment labels of
        # the segments in the image which should be deblended. Larger segments will be
        # filtered out in the "if/for block below" as the processing time to deblend
        # really enormous segments is prohibitive.
        segm_img.big_segments = None
        big_segments = np.where(segm_img.areas >= deb_limit)[0] + 1  # Segment labels are 1-based

        # The biggest_source may be > max_biggest_source indicating there are "big islands"
        # and is_poor_quality should be set to True.  The is_poor_quality is only an indicator that
        # a different kernel type or background computation could be tried for improved results.
        if biggest_source > max_biggest_source:
            log.info("Biggest source %.4f percent EXCEEDS %f percent of the image.", (100.0*biggest_source), (100.0*max_biggest_source))
            is_poor_quality = True
        else:
            log.info("Biggest source %.4f percent within the acceptable %f percent of the image.", (100.0*biggest_source), (100.0*max_biggest_source))

        # Also check if the biggest_source exceeds the limit on the number of pixels allowed
        # for the biggest source, and set is_poor_quality True when the limit is exceeded.
        if biggest_source_pixels > max_biggest_pixels:
            log.info("Biggest source of %d pixels EXCEEDS the maximim limit of %d pixels.", biggest_source_pixels, max_biggest_pixels)
            is_poor_quality = True
        else:
            log.info("Biggest source of %d pixels within the acceptable maximim limit of %d pixels.", biggest_source_pixels, max_biggest_pixels)

        # Filter the big_segments array to remove the prohibitively large segments
        # as it is a very resource consuming task to deblend these segments
        if big_segments.size > 0:
            # Sort the areas and get the indices of the sorted areas
            area_indices = np.argsort(segm_img.areas)
            areas = np.sort(segm_img.areas)
            # Extract only the largest areas which are at least 10x size of next largest segment
            # Value of 10x is set in config file as "ratio_bigsource_deblend_limit"
            # Initialize with largest source, then look for any smaller ones that meet the
            # criteria.  This guarantees that the largest source will always be ignored.
            big_areas = [area_indices[-1] + 1]
            for (i, a), (j, ind) in zip(enumerate(areas[::-1]), enumerate(area_indices[::-1])):
                # skip the first/largest segment, since it was used to initialize the list already
                if i == 0:
                    continue
                # determine ratio of area with next smallest area
                # make sure there are enough entries in areas since the first/largest segment was skipped
                if len(areas) <= 2 or (i+2) > len(areas):
                    break
                r = a / areas[-1*(i+2)]
                if r >= self._ratio_bigsource_deblend_limit:
                    # Record the segmentation ID of large area to be ignored during deblending
                    big_areas.append(ind + 1)
                else:
                    #  we are done with these areas
                    break
            # Remove largest segments from list of 'large' segments to be deblended
            big_segments = big_segments.tolist()
            for b in big_areas:
                if len(big_segments) > 0 and b in big_segments:
                    big_segments.remove(b)
            log.debug("Ignoring {} sources exceeding the deblending limit in size".format(len(big_areas)))
            # Convert back to ndarray
            big_segments = np.array(big_segments)
            if len(big_segments) > 0:
                segm_img.big_segments = big_segments
            else:
                segm_img.big_segments = None
            log.info("Total number of sources suitable for deblending: {}".format(len(big_segments)))
        else:
            log.info("There are no big segments larger than the deblending limit.")

        # Always compute the source_fraction so the value can be reported.  Setting the
        # big_island_only parameter allows control over whether the source_fraction should
        # or should not be ignored.
        source_fraction = segm_img.areas.sum()/real_pixels
        if not big_island_only:
            if source_fraction > max_source_fraction:
                log.info("Total source fraction %.4f percent EXCEEDS %f percent of the image.", (100.0*source_fraction), (100.0*max_source_fraction))
                is_poor_quality = True
            else:
                log.info("Total source fraction %.4f percent within the acceptable %f percent of the image.", (100.0*source_fraction), (100.0*max_source_fraction))
        else:
            log.info("INFORMATIONAL ONLY: Total source fraction %.4f percent.", (100.0*source_fraction))

        return is_poor_quality, biggest_source, source_fraction

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def write_catalog(self, reject_catalogs):
        """Write the specified source catalog out to disk.

        Regardless of the setting for reject_catalogs, the regions file will be written
        solely based upon the setting of diagnostic_mode.

        Parameters
        ----------
        reject_catalogs : dictionary
            Dictionary where the keys are the types of catalogs, and the values are
            bools indicating whether or not the catalogs (*.ecsv) should be written
            with their full contents or as empty catalogs.

        Returns
        -------
        Nothing

        The total product catalog contains several columns contributed from each filter catalog.
        However, there is no guarantee that every filter catalog contains measurements for each
        source in the total product catalog which is to say the row is actually missing from the
        filter product catalog.  When the catalog tables are combined in the combined_tables method,
        the missing entries will contain masked values represented by dashes.  When these values are
        written to output ECSV files, they are written as empty ("") strings.  In order for the catalogs
        to be ingested into a database by possible downstream processing, the masked values will be
        replaced by a numeric indicator.

        ... note : A catalog can have no rows of measurements because no sources were
        found OR because the catalog was "rejected" according to the cosmic ray rejection
        criterion.
        """
        self.source_cat = self.annotate_table(self.source_cat, self.param_dict_qc, proc_type="segment",
                                              product=self.image.ghd_product)

        if reject_catalogs[self.catalog_type]:
            # We still want to write out empty files
            # This will delete all rows from the existing table
            self.source_cat.remove_rows(slice(0, None))

        # Fill the nans and masked values with numeric data
        self.source_cat = fill_nans_maskvalues (self.source_cat, fill_value=constants.FLAG)

        # Write out catalog to ecsv file
        self.source_cat.write(self.sourcelist_filename, format=self.catalog_format, overwrite=True)
        log.info("Wrote catalog file '{}' containing {} sources".format(self.sourcelist_filename, len(self.source_cat)))

        # For debugging purposes only, create a "regions" files to use for ds9 overlay of the segm_img.
        # Create the image regions file here in case there is a failure.  This diagnostic portion of the
        # code should only be invoked when working on the total object catalog (self.segm_img is defined).
        if self.diagnostic_mode:
            # Copy out only the X and Y coordinates to a "diagnostic_mode table" and cast as an Astropy Table
            # so a scalar can be added to the centroid coordinates
            tbl = self.source_cat["X-Centroid", "Y-Centroid"]

            # Construct the diagnostic_mode output filename and write the regions file
            indx = self.sourcelist_filename.find("ecsv")
            outname = self.sourcelist_filename[0:indx] + "reg"

            tbl["X-Centroid"].info.format = ".10f"
            tbl["Y-Centroid"].info.format = ".10f"

            # Add one to the X and Y table values to put the data onto a one-based system,
            # particularly for display with ds9
            tbl["X-Centroid"] = tbl["X-Centroid"] + 1
            tbl["Y-Centroid"] = tbl["Y-Centroid"] + 1
            tbl.write(outname, format="ascii.commented_header")

            log.info("Wrote region file '{}' containing {} sources".format(outname, len(tbl)))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def combine_tables(self, subset_table):
        """Append specified measurements from the filter table to the total detection table.

        The "ID" column is used to map the filter table measurements to the total detection table

        Parameters
        ----------
        subset_table : Astropy table
            A table containing a subset of columns from a filter catalog.
        """
        # Evaluate self.source_cat (the total product list) even though len(self.source_cat) should not be possible
        if len(subset_table) == 0 or len(self.source_cat) == 0:
            log.error("No sources found in the current filter table nor in the total source table.")
            return

        # Keep all the rows in the original total detection table and add columns from the filter
        # table where a matching "id" key is present
        self.source_cat = join(self.source_cat, subset_table, keys="ID", join_type="left")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Utility functions supporting segmentation of total image based on WHT array
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def make_inv_mask(mask):

    invmask = ~(mask.astype(bool))
    invmask = invmask.astype(np.uint8)

    return invmask


def make_wht_masks(
    whtarr: np.ndarray, # image weight array (dtype=float32), zeros outside of footprint
    maskarr: np.ndarray, # mask array (dtype=bool), typically inverse of footprint 
    scale: float = 1.5,
    sensitivity: float = 0.95,
    kernel: tuple = (11, 11) # kernel size for maximum filter window
) -> list: # list containing a dictionary 
    """ Uses scipy's maximum_filter to create the image weight masks of floats between 0 and 1. 
    Function produces a list including a dictionary with the scale (float), wht_limit (float), 
    mask (np.ndarray of bools, dtype=int16), and  relative weight (np.ndarray, dtype=float32). 
    
    """
    # uses scipy maximum filter on image. Maximum filter selects the largest value within an ordered 
    # window of pixel values and replaces the central pixel with the largest value.
    maxwht = ndimage.filters.maximum_filter(whtarr, size=kernel)
    
    # normalized weight array
    rel_wht = maxwht / maxwht.max()

    # initialize values
    #
    # master_mask indicates pixels remaining to be included
    # initially it includes all unmasked pixels
    master_mask = (maskarr == 0) & (rel_wht > 0)
    limit = 1 / scale
    masks = []

    # minimum threshold on number of pixels in a segment
    # this is weird, trying to preserve the peculiar definition of similarity
    osum = master_mask.sum()
    pthresh = (1 - sensitivity) * osum
    
    # loop through scale values until all pixels are used
    while master_mask.any():

        # get pixels above the new limit that have not been used yet
        mask = (rel_wht > limit) & master_mask
        if mask.any():
            mask_sum = mask.sum()

            # If the segment is sufficiently big, and if the remaining pixels are not too small,
            # keep this segment.
            mmaster_sum = master_mask.sum()
            if mask_sum == mmaster_sum or (mask_sum > pthresh and (mmaster_sum - mask_sum) > pthresh):
                masks.append(
                    dict(
                        scale=limit,
                        wht_limit=limit * maxwht.max(),
                        mask=mask.astype(np.uint16),
                        rel_weight=rel_wht * mask)
                )
                master_mask &= (~mask)
            elif (mmaster_sum - mask_sum) <= pthresh:
                # few pixels remaining, drop threshold to zero so they are found in the next iteration
                limit = 0.0

        limit /= scale

    return masks

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Utility functions supporting point and segmentation catalogs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def fill_nans_maskvalues(catalog, fill_value=constants.FLAG):

    # Fill the masked values with fill_value - the value is truncated for int as datatype of column is known
    catalog = catalog.filled(fill_value)

    # Also fill any nan values with fill_value
    for col in catalog.itercols():
        np.nan_to_num(col, copy=False, nan=fill_value)

    return catalog


def enforce_icrs_compatibility(catalog):
    """This function insures that the source catalog can return ICRS sky coordinates.

    # This function guards against a non-standard coordinate system being defined in the input image headers
    # If something non-standard is found, set the value to 'ICRS' to insure that the
    # table conversion works.  Checking against the set of valid options recognized by WCSLIB and
    # documented in http://ds9.si.edu/doc/ref/region.html#RegionFileFormat is necessary since the
    # header keyword REFFRAME can be populated with anything specified by the user in the original
    # proposal.
    """
    # We need to check whether or not the catalog was generated with photutils<1.6.0 or not
    # since photutils=1.6.0 (apparently) renamed the 'catalog._wcs' to 'catalog.wcs'.
    cat_wcs = catalog.wcs if hasattr(catalog, 'wcs') else catalog._wcs
    if cat_wcs.wcs.radesys.upper() not in RADESYS_OPTIONS:
        cat_wcs.wcs.radesys = 'ICRS'
        log.warning(f"Assuming input coordinates are ICRS, instead of {cat_wcs.wcs.radesys}")
        log.warning(f"Sky coordinates of source objects may not be accurate.")
