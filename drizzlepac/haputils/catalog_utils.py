"""This script contains code to support creation of photometric sourcelists using two techniques:
aperture photometry and segmentation-map based photometry."""

import copy
import pickle  # FIX Remove
import sys
from distutils.version import LooseVersion

from astropy.io import fits as fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Column, MaskedColumn, Table, join, vstack
from astropy.convolution import RickerWavelet2DKernel
from astropy.coordinates import SkyCoord
import numpy as np
from scipy import ndimage, stats

import photutils  # needed to check version
if LooseVersion(photutils.__version__) < '1.1.0':
    OLD_PHOTUTILS = True
    from photutils.segmentation import (detect_sources, source_properties,
                                        deblend_sources)
else:
    OLD_PHOTUTILS = False
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
from . import photometry_tools
from . import deconvolve_utils as decutils
from . import processing_utils as proc_utils

try:
    from matplotlib import pyplot as plt
except Exception:
    plt = None

CATALOG_TYPES = ['aperture', 'segment']

if OLD_PHOTUTILS:
    id_colname = 'id'
    flux_colname = 'source_sum'
    ferr_colname = 'source_sum_err'
    bac_colname = 'background_at_centroid'
else:
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

        # Get the HSTWCS object from the first extension
        self.imgwcs = HSTWCS(self.imghdu, 1)

        self.keyword_dict = self._get_header_data()

        # Populated by self.compute_background()
        self.bkg_background_ra = None
        self.bkg_rms_ra = None
        self.bkg_rms_median = None
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

        if self.bkg_background_ra is None:
            self.compute_background(box_size, win_size,
                                    simple_bkg=simple_bkg,
                                    bkg_skew_threshold=bkg_skew_threshold,
                                    zero_percent=zero_percent,
                                    negative_percent=negative_percent,
                                    nsigma_clip=nsigma_clip,
                                    maxiters=maxiters)

        log.info("Attempt to determine FWHM based upon input data within a good FWHM range of {:.1f} to {:.1f}.".format(good_fwhm[0], good_fwhm[1]))
        log.info("If no good FWHM candidate is identified, a value of {:.1f} will be used instead.".format(fwhmpsf / self.imgwcs.pscale))
        k, self.kernel_fwhm = astrometric_utils.build_auto_kernel(self.data,
                                                                  self.wht_image,
                                                                  good_fwhm=good_fwhm,
                                                                  num_fwhm=30,
                                                                  threshold=self.bkg_rms_ra,
                                                                  fwhm=fwhmpsf / self.imgwcs.pscale)
        (self.kernel, self.kernel_psf) = k

    def compute_background(self, box_size, win_size,
                           bkg_estimator=SExtractorBackground, rms_estimator=StdBackgroundRMS,
                           simple_bkg=False,
                           bkg_skew_threshold=0.5,
                           zero_percent=25.0,
                           negative_percent=15.0,
                           nsigma_clip=3.0,
                           maxiters=3):
        """Use a sigma-clipped algorithm or Background2D to determine the background of the input image.

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

        # SExtractorBackground ans StdBackgroundRMS are the defaults
        bkg = None
        is_zero_background_defined = False

        # Make a local copy of the data(image) being processed in order to reset any
        # data values which equal nan (e.g., subarrays) to zero.
        imgdata = np.nan_to_num(self.data, copy=True, nan=0.0)

        # In order to compute the proper statistics on the input data, need to use the footprint
        # mask to get the actual data - illuminated portion (True), non-illuminated (False).
        footprint_mask = self.num_images_mask > 0
        self.footprint_mask = ndimage.binary_erosion(footprint_mask, iterations=10)
        self.inv_footprint_mask = np.invert(self.footprint_mask)

        # If the image contains a lot of values identically equal to zero (as in some SBC images),
        # set the two-dimensional background image to a constant of zero and the background rms to
        # the real rms of the non-zero values in the image.
        num_of_illuminated_pixels = self.footprint_mask.sum()
        num_of_zeros = np.count_nonzero(imgdata[self.footprint_mask] == 0)
        non_zero_pixels = imgdata[self.footprint_mask]

        # BACKGROUND COMPUTATION 1 (unusual case)
        # If there are too many background zeros in the image (> number_of_zeros_in_background_threshold), set the
        # background median and background rms values
        if num_of_zeros / float(num_of_illuminated_pixels) * 100.0 > zero_percent:
            self.bkg_median = 0.0
            self.bkg_rms_median = stats.tstd(non_zero_pixels, limits=[0, None], inclusive=[False, True])
            self.bkg_background_ra = np.full_like(imgdata, 0.0)
            self.bkg_rms_ra = np.full_like(imgdata, self.bkg_rms_median)
            self.bkg_type = 'zero_background'

            is_zero_background_defined = True
            log.info("Input image contains excessive zero values in the background. Median: {} RMS: {}".format(self.bkg_median, self.bkg_rms_median))

        # BACKGROUND COMPUTATION 2 (sigma_clipped_stats)
        # If the input data is not the unusual case of an "excessive zero background", compute
        # a sigma-clipped background which returns only single values for mean,
        # median, and standard deviations
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
            # based on full unmasked image
            bkg_skew = np.abs(3.0 * (bkg_mean_full - bkg_median_full) / bkg_rms_full)
            log.info("Sigma-clipped computed skewness: {0:.2f}".format(bkg_skew))


            # Refine background to better compute the median value
            imgnz = imgdata * self.footprint_mask
            imgnz = imgnz[imgnz > 0.0]  # only want non-negative values

            imgvals = imgnz[imgnz < (bkg_median_full + (bkg_rms_full * 0.1))]
            bkg_mean, bkg_median, bkg_rms = sigma_clipped_stats(imgvals,
                                                                None,
                                                                sigma=nsigma_clip,
                                                                cenfunc='median')

            log.info("Sigma-clipped Statistics - Background mean: {}  median: {}  rms: {}".format(bkg_mean, bkg_median, bkg_rms))
            log.info("")

            # Ensure the computed values are not negative
            if bkg_mean < 0.0 or bkg_median < 0.0 or bkg_rms < 0.0:
                bkg_mean = max(0, bkg_mean)
                bkg_median = max(0, bkg_median)
                bkg_rms = max(0, bkg_rms)
                log.info("UPDATED Sigma-clipped Statistics - Background mean: {}  median: {}  rms: {}".format(bkg_mean, bkg_median, bkg_rms))
                log.info("")

            # Compute a minimum rms value based upon information directly from the data
            if self.keyword_dict["detector"].upper() != "SBC":
                minimum_rms = self.keyword_dict['atodgn'] * self.keyword_dict['readnse'] \
                              * self.keyword_dict['ndrizim'] / self.keyword_dict['texpo_time']

                # Compare a minimum rms based upon input characteristics versus the one computed and use
                # the larger of the two values.
                if (bkg_rms < minimum_rms):
                    bkg_rms = minimum_rms
                    log.info("Mimimum RMS of input based upon the readnoise, gain, number of exposures, and total exposure time: {}".format(minimum_rms))
                    log.info("Sigma-clipped RMS has been updated - Background mean: {}  median: {}  rms: {}".format(bkg_mean, bkg_median, bkg_rms))
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
                                           exclude_percentile=percentile, edge_method="pad",
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
                    log.info("Percentage of negative values in the background subtracted image {0:.2f} vs low threshold of {1:.2f}.".format(100.0 * negative_ratio, negative_percent))

                    # If the background subtracted image has too many negative values which may be
                    # indicative of large negative regions, the two-dimensional computed background
                    # fit image should NOT be used.  Use the sigma-clipped data instead.
                    if negative_ratio * 100.0 > negative_percent:
                        log.info("Percentage of negative values {0:.2f} in the background subtracted image exceeds the threshold of {1:.2f}.".format(100.0 * negative_ratio, negative_percent))
                        log.info("")
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
        log.info("    Median background: {}".format(self.bkg_median))
        log.info("    Median RMS background: {}".format(self.bkg_rms_median))
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
        if keyword_dict["detector"].upper() != "SBC":
            keyword_dict["ccd_gain"] = self.imghdu[0].header["CCDGAIN"]
            keyword_dict["readnse"] = self._get_max_key_value(self.imghdu[0].header, 'READNSE')
            keyword_dict["atodgn"] = self._get_max_key_value(self.imghdu[0].header, 'ATODGN')
        keyword_dict["aperture_pa"] = self.imghdu[0].header["PA_V3"]
        keyword_dict["gain_keys"] = [self.imghdu[0].header[k[:8]] for k in self.imghdu[0].header["ATODGN*"]]

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
            else:
                keyword_dict["filter1"] = self.imghdu[0].header["FILTER"]
                keyword_dict["filter2"] = ""

        # Get the HSTWCS object from the first extension
        keyword_dict["wcs_name"] = self.imghdu[1].header["WCSNAME"]
        keyword_dict["wcs_type"] = self.imghdu[1].header["WCSTYPE"]
        keyword_dict["orientation"] = self.imghdu[1].header["ORIENTAT"]
        keyword_dict["aperture_ra"] = self.imghdu[1].header["RA_APER"]
        keyword_dict["aperture_dec"] = self.imghdu[1].header["DEC_APER"]
        keyword_dict["photflam"] = proc_utils.find_flt_keyword(self.imghdu, "PHOTFLAM")
        keyword_dict["photplam"] = proc_utils.find_flt_keyword(self.imghdu, "PHOTPLAM")

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
    crfactor = {'aperture': 300, 'segment': 150}  # CRs / hr / 4kx4k pixels

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
                                self.param_dict['dao']['TWEAK_FWHMPSF'],
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
        self.catalogs = {}
        if 'segment' in self.types:
            self.catalogs['segment'] = HAPSegmentCatalog(self.image, self.param_dict, self.param_dict_qc,
                                                         self.diagnostic_mode, tp_sources=tp_sources)

        if 'aperture' in self.types:
            self.catalogs['aperture'] = HAPPointCatalog(self.image, self.param_dict, self.param_dict_qc,
                                                        self.diagnostic_mode, tp_sources=tp_sources)

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

        ... note : If either catalog fails the following test, then both are rejected.
                        n_cat < thresh
                   where
                        thresh = crfactor * n1_exposure_time**2 / texptime
        """
        for cat_type in self.catalogs:
            crthresh_mask = None
            source_cat = self.catalogs[cat_type].sources if cat_type == 'aperture' else self.catalogs[cat_type].source_cat

            flag_cols = [colname for colname in source_cat.colnames if colname.startswith('Flag')]
            for colname in flag_cols:
                catalog_crmask = source_cat[colname] < 2
                if crthresh_mask is None:
                    crthresh_mask = catalog_crmask
                else:
                    # Combine masks for all filters for this catalog type
                    crthresh_mask = np.bitwise_or(crthresh_mask, catalog_crmask)
            source_cat.sources_num_good = len(np.where(crthresh_mask)[0])

        reject_catalogs = False

        log.info("Determining whether point and/or segment catalogs meet cosmic-ray threshold")
        log.info("  based on EXPTIME = {}sec for the n=1 filters".format(n1_exposure_time))

        for cat_type in self.catalogs:
            source_cat = self.catalogs[cat_type]
            if source_cat.sources:
                thresh = self.crfactor[cat_type] * n1_exposure_time**2 / self.image.keyword_dict['texpo_time']
                source_cat = source_cat.sources if cat_type == 'aperture' else source_cat.source_cat
                n_sources = source_cat.sources_num_good  # len(source_cat)
                all_sources = len(source_cat)
                log.info("{} catalog with {} good sources out of {} total sources :  CR threshold = {}".format(cat_type, n_sources, all_sources, thresh))
                if n_sources < thresh:
                    reject_catalogs = True
                    log.info("{} catalog FAILED CR threshold.  Rejecting both catalogs...".format(cat_type))
                    break

        return reject_catalogs

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
        reject_catalogs : bool
            Indicator as to whether or not the catalogs (*.ecsv) should be written.

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

        # Compute average gain - there will always be at least one gain value in the primary header
        gain_keys = self.image.keyword_dict['gain_keys']
        gain_values = [g for g in gain_keys if g > 0.0]
        self.gain = self.image.keyword_dict['exptime'] * np.mean(gain_values)

        # Set the gain for ACS/SBC and WFC3/IR to 1.0
        if self.image.keyword_dict["detector"].upper() in ["IR", "SBC"]:
            self.gain = 1.0

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
        self.tp_masks = make_wht_masks(self.image.wht_image, self.image.inv_footprint_mask,
                                       scale=scale,
                                       sensitivity=self.param_dict['sensitivity'],
                                       kernel=(self.param_dict['region_size'],
                                               self.param_dict['region_size']))

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
        num_sources = len(data_table)
        data_table.meta["Number of sources"] = num_sources

        proc_type = proc_type.lower()
        ci_lower = float(param_dict_qc['ci filter'][proc_type]['ci_lower_limit'])
        ci_upper = float(param_dict_qc['ci filter'][proc_type]['ci_upper_limit'])

        data_table.meta["h09"] = ["#================================================================================================="]
        data_table.meta["h10"] = ["IMPORTANT NOTES"]
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

        if proc_type is "segment":
            if self.is_big_island:
                data_table.meta["h19"] = ["WARNING: Segmentation catalog is considered to be of poor quality due to a crowded field or large segments."]

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

                        daofind = decutils.UserStarFinder(fwhm=source_fwhm,
                                                coords=user_peaks,
                                                threshold=0.0,
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
                if reg_sources is not None:
                    if sources is None:
                        sources = reg_sources
                    else:
                        sources = vstack([sources, reg_sources])

            # If there are no detectable sources in the total detection image, return as there is nothing more to do.
            if not sources:
                log.warning("No point sources were found in Total Detection Product, {}.".format(self.imgname))
                log.warning("Processing for point source catalogs for this product is ending.")
                self._define_empty_table()
                return

            log.info("Measured {} sources in {}".format(len(sources), self.image.imgname))
            log.info("   colnames: {}".format(sources.colnames))

            # calculate and add RA and DEC columns to table
            ra, dec = self.transform_list_xy_to_ra_dec(sources["xcentroid"], sources["ycentroid"], self.imgname)
            ra_col = Column(name="RA", data=ra, dtype=np.float64)
            dec_col = Column(name="DEC", data=dec, dtype=np.float64)
            sources.add_column(ra_col, index=3)
            sources.add_column(dec_col, index=4)

            for col in sources.colnames:
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
            final_col_units = {"xcentroid": "pixels", "ycentroid": "pixels", "RA": "degrees", "DEC": "degrees",
                               "id": "unitless"}
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

        final_col_units = {"xcentroid": "pixels", "ycentroid": "pixels", "RA": "degrees", "DEC": "degrees",
                           "id": "unitless", "Flags": "unitless"}

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
        try:
            positions = (self.sources['xcentroid'], self.sources['ycentroid'])
        except Exception:
            positions = (self.sources['X-Center'], self.sources['Y-Center'])

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
            ci_data = photometry_tbl["MagAp1"].data - photometry_tbl["MagAp2"].data
        except Exception:
            log.info("Wrote catalog info to file 'catalog.pickle'.")
            pickle_out = open("catalog.pickle", "wb")
            pickle.dump(photometry_tbl, pickle_out)
            pickle_out.close()

        ci_mask = np.logical_and(np.abs(ci_data) > 0.0, np.abs(ci_data) < 1.0e-30)
        big_bad_index = np.where(abs(ci_data) > 1.0e20)
        ci_mask[big_bad_index] = True
        ci_col = MaskedColumn(name="CI", data=ci_data, dtype=np.float64, mask=ci_mask)
        photometry_tbl.add_column(ci_col)

        # Add zero-value "Flags" column in preparation for source flagging
        flag_col = Column(name="Flags", data=np.zeros_like(photometry_tbl['ID']), dtype=np.int64)
        photometry_tbl.add_column(flag_col)

        # build final output table
        final_col_order = ["X-Center", "Y-Center", "RA", "DEC", "ID", "MagAp1", "MagErrAp1", "MagAp2", "MagErrAp2",
                           "MSkyAp2", "StdevAp2", "FluxAp2", "CI", "Flags"]
        output_photometry_table = photometry_tbl[final_col_order]

        # format output table columns
        final_col_format = {"X-Center": "10.3f", "Y-Center": "10.3f", "RA": "13.7f", "DEC": "13.7f", "ID": "7d",
                            "MagAp1": '7.3f', "MagErrAp1": '7.3f', "MagAp2": '7.3f',
                            "MagErrAp2": '7.3f', "MSkyAp2": '7.3f', "StdevAp2": '7.3f',
                            "FluxAp2": '10.4f', "CI": "7.3f", "Flags": "5d"}  # TODO: Standardize precision
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
                             "MagAp2": "ABMAG of source based on the outer (larger) aperture",
                             "MagErrAp2": "Error of MagAp2",
                             "MSkyAp2": "ABMAG of sky based on outer (larger) aperture",
                             "StdevAp2": "Standard deviation of sky measurement in outer (larger) aperture",
                             "FluxAp2": "Flux of source based on the outer (larger) aperture",
                             "CI": "Concentration Index",
                             "Flags": "Numeric encoding for conditions on detected sources"}
        for fcd_key in final_col_descrip.keys():
            output_photometry_table[fcd_key].description = final_col_descrip[fcd_key]

        # add units to columns
        final_col_units = {"X-Center": "pixels", "Y-Center": "pixels", "RA": "degrees", "DEC": "degrees",
                           "ID": "unitless", "MagAp1": "ABMAG", "MagErrAp1": "ABMAG", "MagAp2": "ABMAG",
                           "MagErrAp2": "ABMAG", "MSkyAp2": "ABMAG", "StdevAp2": "ABMAG",
                           "FluxAp2": "electrons/sec", "CI": "ABMAG", "Flags": "unitless"}
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
        reject_catalogs : bool
            Indicator as to whether or not the catalogs (*.ecsv) should be written.

        Returns
        -------
        Nothing!

        """
        # Insure catalog has all necessary metadata
        self.source_cat = self.annotate_table(self.source_cat, self.param_dict_qc, proc_type="aperture",
                                              product=self.image.ghd_product)
        if reject_catalogs:
            # We still want to write out empty files
            # This will delete all rows from the existing table
            self.source_cat.remove_rows(slice(0, None))

        # Write out catalog to ecsv file
        # self.source_cat.meta['comments'] = \
        #     ["NOTE: The X and Y coordinates in this table are 0-indexed (i.e. the origin is (0,0))."]
        self.source_cat.write(self.sourcelist_filename, format=self.catalog_format)
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
        self.sources = join(self.sources, subset_table, keys="ID", join_type="left")

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
        self._fwhm = self.param_dict["sourcex"]["fwhm"]
        self._size_source_box = self.param_dict["sourcex"]["source_box"]
        self._nlevels = self.param_dict["sourcex"]["nlevels"]
        self._contrast = self.param_dict["sourcex"]["contrast"]
        self._border = self.param_dict["sourcex"]["border"]
        self._nsigma = self.param_dict["sourcex"]["segm_nsigma"]
        self._rw2d_size = self.param_dict["sourcex"]["rw2d_size"]
        self._rw2d_nsigma = self.param_dict["sourcex"]["rw2d_nsigma"]
        self._rw2d_biggest_source = self.param_dict["sourcex"]["rw2d_biggest_source"]
        self._rw2d_source_fraction = self.param_dict["sourcex"]["rw2d_source_fraction"]
        self._bs_deblend_limit = self.param_dict["sourcex"]["biggest_source_deblend_limit"]
        self._sf_deblend_limit = self.param_dict["sourcex"]["source_fraction_deblend_limit"]
        self._ratio_bigsource_limit = self.param_dict["sourcex"]["ratio_bigsource_limit"]

        # Columns to include from the computation of source properties to save
        # computation time from computing values which are not used
        self.include_filter_cols = ['area', bac_colname, 'bbox_xmax', 'bbox_xmin', 'bbox_ymax', 'bbox_ymin',
                                    'covar_sigx2', 'covar_sigxy', 'covar_sigy2', 'cxx', 'cxy', 'cyy',
                                    'ellipticity', 'elongation', id_colname, 'orientation', 'sky_centroid_icrs',
                                    flux_colname, ferr_colname, 'xcentroid', 'ycentroid']

        # Initialize attributes to be computed later
        self.segm_img = None  # Segmentation image

        # Defined in measure_sources
        self.subset_filter_source_cat = None

        # Default kernel which may be the custom kernel based upon the actual image
        # data or a Gaussian 2D kernel. This may be over-ridden in identify_sources().
        self.kernel = copy.deepcopy(self.image.kernel)

        # Attribute computed when generating the segmentation image.  If the segmentation image
        # is deemed to be of poor quality, make sure to add documentation to the output catalog.
        self.is_big_island = False

    def identify_sources(self, **pars):
        """Use photutils to find sources in image based on segmentation.

        Returns
        -------
        self.sources
        self.source_cat

        Defines
        -------
        self.segm_img : `photutils.segmentation.SegmentationImage`
            Two-dimensional segmentation image where found source regions are labeled with
            unique, non-zero positive integers.
        """

        # If the total product sources have not been identified, then this needs to be done!
        if not self.tp_sources:

            # Report configuration values to log
            log.info("{}".format("=" * 80))
            log.info("")
            log.info("SExtractor-like source finding settings - Photutils segmentation")
            log.info("Total Detection Product - Input Parameters")
            log.info("Image: {}".format(self.imgname))
            log.info("FWHM: {}".format(self._fwhm))
            log.info("size_source_box (no. of connected pixels needed for a detection): {}".format(self._size_source_box))
            log.info("nsigma (threshold = nsigma * background_rms): {}".format(self._nsigma))
            log.info("nlevels (no. of multi-thresholding levels for deblending): {}".format(self._nlevels))
            log.info("contrast (frac. flux for peak to be separate object, 0=max. deblend, 1=no deblend): {}".format(self._contrast))
            log.info("RickerWavelet nsigma (threshold = nsigma * background_rms): {}".format(self._rw2d_nsigma))
            log.info("RickerWavelet kernel X- and Y-dimension: {}".format(self._rw2d_size))
            log.info("Percentage limit on biggest source (criterion for  RickerWavelet kernel): {}".format(100.0 * self._rw2d_biggest_source))
            log.info("Percentage limit on source fraction over the image (criterion for RickerWavelet kernel): {}".format(100.0 * self._rw2d_source_fraction))
            log.info("Percentage limit on biggest source deblending limit: {}".format(100.0 * self._bs_deblend_limit))
            log.info("Percentage limit on source fraction deblending limit: {}".format(100.0 * self._sf_deblend_limit))

            # Get the SCI image data
            imgarr = copy.deepcopy(self.image.data)

            # Custom or Gaussian kernel depending upon the results of CatalogImage build_kernel()
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
            log.info("Using Custom kernel or Gaussian to generate a segmentation map.")
            g_segm_img, g_is_big_crowded, g_bs, g_sf = self.detect_and_eval_segments(imgarr,
                                                                                     g2d_kernel,
                                                                                     ncount,
                                                                                     self._size_source_box,
                                                                                     self._nsigma,
                                                                                     self.image.bkg_background_ra,
                                                                                     self.image.bkg_rms_ra,
                                                                                     check_big_island_only=False,
                                                                                     rw2d_biggest_source=self._rw2d_biggest_source,
                                                                                     rw2d_source_fraction=self._rw2d_source_fraction)
            segm_img_orig = copy.deepcopy(g_segm_img)

            # If the science field via the segmentation map is deemed crowded or has big sources/islands, compute the
            # RickerWavelet2DKernel and call detect_and_eval_segments() again. Still use the custom fwhm as it
            # should be better than a generic fwhm as it is based upon the data.
            # Note: the fwhm might be a default if the custom algorithm had to fall back to a Gaussian.
            if g_is_big_crowded and g_segm_img:
                log.info("")
                log.info("The segmentation map contains big sources/islands or a large source fraction of segments.")
                log.info("Using RickerWavelet2DKernel to generate an alternate segmentation map.")
                rw2d_kernel = RickerWavelet2DKernel(self.image.kernel_fwhm,
                                                    x_size=self._rw2d_size,
                                                    y_size=self._rw2d_size)
                rw2d_kernel.normalize()

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
                                                                                             rw2d_biggest_source=self._rw2d_biggest_source,
                                                                                             rw2d_source_fraction=self._rw2d_source_fraction)

                # Compute the ratio of big sources/islands using Custom/Gaussian kernel vs Rickerwavelet kernel
                # This value can be used as a discriminant between overlapping point sources and nebulousity fields
                ratio_cg2rw_bigsource = 3.0
                if rw_bs > 0.0:
                    ratio_cg2rw_bigsource = g_bs / rw_bs

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
                                                self.param_dict['dao']['TWEAK_FWHMPSF'])

                        # Reset the local version of the Custom/Gaussian kernel and the RickerWavelet
                        # kernel when the background type changes
                        g2d_kernel = self.image.kernel

                        rw2d_kernel = RickerWavelet2DKernel(self.image.kernel_fwhm,
                                                            x_size=self._rw2d_size,
                                                            y_size=self._rw2d_size)
                        rw2d_kernel.normalize()

                        sigma_for_threshold = self._nsigma
                        rw2d_sigma_for_threshold = self._rw2d_nsigma

                    # Re-compute a background2D with a higher threshold by increasing the nsigma used
                    elif (self.image.bkg_type.lower().startswith('twod')):
                        log.info("Increasing the threshold image (bkg + nsigma * 2.0) for improved source detection.")
                        sigma_for_threshold = self._nsigma * 2.0
                        rw2d_sigma_for_threshold = self._rw2d_nsigma * 2.0

                    # Detect segments and evaluate the detection in terms of big sources/islands or crowded fields
                    # Round 2
                    ncount += 1
                    log.info("")
                    log.info("With alternate background...using Custom/Gaussian kernel to generate a segmentation map.")
                    del g_segm_img
                    g_segm_img, g_is_big_crowded, g_bs, g_sf = self.detect_and_eval_segments(imgarr,
                                                                                             g2d_kernel,
                                                                                             ncount,
                                                                                             self._size_source_box,
                                                                                             sigma_for_threshold,
                                                                                             self.image.bkg_background_ra,
                                                                                             self.image.bkg_rms_ra,
                                                                                             check_big_island_only=False,
                                                                                             rw2d_biggest_source=self._rw2d_biggest_source,
                                                                                             rw2d_source_fraction=self._rw2d_source_fraction)

                    # Check again for big sources/islands or a large source fraction
                    if g_is_big_crowded:
                        log.info("")
                        log.info("The segmentation map contains big sources/islands or a large source fraction of segments.")
                        log.info("With alternate background...using RickerWavelet2DKernel to generate an alternate segmentation map.")

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
                                                                                                     check_big_island_only=False,
                                                                                                     rw2d_biggest_source=self._bs_deblend_limit,
                                                                                                     rw2d_source_fraction=self._sf_deblend_limit)

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
                        if not rw_is_big_crowded or (rw_is_big_crowded and (ratio_cg2rw_bigsource > self._ratio_bigsource_limit)):
                            log.info("The Round 2 of segmentation images may still contain big sources/islands.\n"
                                     "However, the ratio between the Custom/Gaussian and Rickerwavelet biggest source is\n"
                                     "indicative of overlapping PSFs vs nebulousity.")
                            log.info("Proceeding as the time to deblend should be nominal.")
                            self.kernel = rw2d_kernel
                            segm_img = copy.deepcopy(rw_segm_img)
                            del rw_segm_img
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

                    # Use the second round custom/Gaussian segmentation image
                    else:
                        self.kernel = g2d_kernel
                        segm_img = copy.deepcopy(g_segm_img)
                        del g_segm_img

                # The first round RickerWavelet segmentation image is good, continue with the processing
                elif not rw_is_big_crowded and rw_segm_img:
                    self.kernel = rw2d_kernel
                    segm_img = copy.deepcopy(rw_segm_img)
                    del rw_segm_img
                    del g_segm_img

                # No segments were detected in the total data product - no further processing done for this TDP,
                # but processing of another TDP should proceed.
                elif not rw_segm_img:
                    self._define_empty_table(rw_segm_img)
                    return

            # The first round custom/Gaussian segmentation image is good, continue with the processing
            elif not g_is_big_crowded and g_segm_img:
                self.kernel = g2d_kernel
                segm_img = copy.deepcopy(g_segm_img)
                del g_segm_img

            # No segments were detected in the total data product - no further processing done for this TDP,
            # but processing of another TDP should proceed.
            elif not g_segm_img:
                self._define_empty_table(g_segm_img)
                return

            # Deblend the segmentation image
            ncount += 1
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
            if OLD_PHOTUTILS:
                self.source_cat = source_properties(imgarr, self.segm_img, background=self.image.bkg_background_ra,
                                                    filter_kernel=self.kernel, wcs=self.image.imgwcs)
            else:
                self.source_cat = SourceCatalog(imgarr, self.segm_img, background=self.image.bkg_background_ra,
                                                    kernel=self.kernel, wcs=self.image.imgwcs)

            # Convert source_cat which is a SourceCatalog to an Astropy Table - need the data in tabular
            # form to filter out bad rows and correspondingly bad segments before the filter images are processed.
            total_measurements_table = Table(self.source_cat.to_table(columns=['label', 'xcentroid', 'ycentroid', 'sky_centroid_icrs']))

            # Filter the table to eliminate nans or inf based on the coordinates, then remove these segments from
            # the segmentation image too
            good_rows = []
            bad_segm_rows_by_id = []
            updated_table = None
            for i, old_row in enumerate(total_measurements_table):
                if np.isfinite(old_row["xcentroid"]):
                    good_rows.append(old_row)
                else:
                    bad_segm_rows_by_id.append(total_measurements_table['label'][i])
            updated_table = Table(rows=good_rows, names=total_measurements_table.colnames)
            if self.diagnostic_mode and bad_segm_rows_by_id:
                log.info("Bad segments removed from segmentation image for Total detection image {}.".format(self.imgname))

            # Remove the bad segments from the image
            self.segm_img.remove_labels(bad_segm_rows_by_id, relabel=True)

            # Clean up the existing column names, format, and descriptions
            self.source_cat = self._define_total_table(updated_table)

            # self.sources needs to be passed to a filter catalog object based on code in hapsequencer.py
            # (create_catalog_products()).  This is the way the independent catalogs of total and filter products
            # process the same segmentation image.
            # BEWARE: self.sources for "segmentation" is a SegmentationImage, but for "point" it is an Astropy table
            self.sources = copy.deepcopy(self.segm_img)

            log.info("Done identifying sources in total detection image for the segmentation catalog.")
            log.info("")
            log.info("{}".format("=" * 80))
            log.info("")

        # If filter product, use sources identified in total detection product previously generated
        else:
            self.sources = self.tp_sources['segment']['sources']
            self.kernel = self.tp_sources['segment']['kernel']
            self.total_source_table = self.tp_sources['segment']['source_cat']

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

            tbl["X-Centroid"].info.format = ".10f"
            tbl["Y-Centroid"].info.format = ".10f"

            # Add one to the X and Y table values to put the data onto a one-based system,
            # particularly for display with ds9
            tbl["X-Centroid"] = tbl["X-Centroid"] + 1
            tbl["Y-Centroid"] = tbl["Y-Centroid"] + 1
            tbl.write(outname, format="ascii.commented_header")

            log.info("Wrote region file '{}' containing {} sources".format(outname, len(tbl)))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def detect_and_eval_segments(self, imgarr, kernel, ncount, size_source_box, nsigma_above_bkg, background_img, background_rms, check_big_island_only=False, rw2d_biggest_source=0.015, rw2d_source_fraction=0.075):

            # Compute the threshold to use for source detection
            threshold = self.compute_threshold(nsigma_above_bkg, background_img, background_rms)

            # Write out diagnostic data
            if self.diagnostic_mode:
                outname = self.imgname.replace(".fits", "_threshold" + str(ncount) + ".fits")
                fits.PrimaryHDU(data=threshold).writeto(outname)

            # Generate the segmentation map by detecting "sources" using the nominal settings.
            # Use all the parameters here developed for the "custom kernel".  Note: if the
            # "custom kernel" did not work out, build_auto_kernel() drops back to a Gaussian.
            # log.info('Kernel shape: {}    source_box: {}'.format(g2d.shape, self._size_source_box))
            segm_img = self.detect_segments(imgarr,
                                            threshold,
                                            ncount,
                                            filter_kernel=kernel,
                                            source_box=size_source_box,
                                            mask=self.image.inv_footprint_mask)

            # Check if custom_segm_image is None indicating there are no detectable sources in this
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
                                                                                                max_biggest_source=rw2d_biggest_source,
                                                                                                max_source_fraction=rw2d_source_fraction)

            return segm_img, is_big_crowded, big_island, source_fraction

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def compute_threshold(self, nsigma, bkg_mean, bkg_rms):
        """Compute the threshold value above which sources are deemed detected.

           Parameters
           ----------
           nsigma : float
               Multiplicative factor for the background RMS

           bkg_mean : float image
               Mean of the background determined image

           bkg_rms : float image
               RMS of the background determined image

           Returns
           -------
           threshold: float image
               Image which defines, on a pixel-by-pixel basis, the low limit above which
               sources are detected.
        """

        log.info("Computing the threshold value used for source detection.")
        if not self.tp_masks:
            threshold = bkg_mean + (nsigma * bkg_rms)
        else:
            threshold = np.zeros_like(self.tp_masks[0]['rel_weight'])
            log.info("Using WHT masks as a scale on the RMS to compute threshold detection limit.")
            for wht_mask in self.tp_masks:
                threshold_rms = bkg_rms * np.sqrt(wht_mask['scale'] * wht_mask['mask'] / wht_mask['rel_weight'].max())
                threshold_rms_median = np.nanmedian(threshold_rms[threshold_rms > 0])
                threshold_item = bkg_mean + (nsigma * threshold_rms_median)
                threshold += threshold_item
            del(threshold_item)

        return threshold

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def detect_segments(self, imgarr, threshold, ncount, filter_kernel=None, source_box=7, mask=None):
        """Detect segments found in the input total detection (aka white light) image.

           Image regions are identified as segments in the 'total detection image' if the region
           has n-connected pixels with values greater than the 'threshold'. The resultant
           segmentation image is then evaluated to determine if the segments represent a large portion
           of the image or if any one segment is very large.

           Parameters
           ----------
           imgarr :
               Total detection image (no background subtraction)

           threshold :
               Image which defines, on a pixel-by-pixel basis, the low limit above which
               sources are detected.

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
        # Note: SExtractor has "connectivity=8" which is the default for detect_sources().
        segm_img = None
        segm_img = detect_sources(imgarr,
                                  threshold,
                                  npixels=source_box,
                                  filter_kernel=filter_kernel,
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

    def deblend_segments(self, segm_img, imgarr, ncount, filter_kernel=None, source_box=5):
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
            # npixels = number of connected pixels in source
            # npixels and filter_kernel should match those used by detect_sources()
            segm_deblended_img = deblend_sources(imgarr,
                                                 segm_img,
                                                 npixels=source_box,
                                                 filter_kernel=filter_kernel,
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
        log.info("FWHM: {}".format(self._fwhm))
        log.info("size_source_box: {}".format(self._size_source_box))
        log.info("")

        # This is the filter science data and its computed background
        imgarr_bkgsub = imgarr - self.image.bkg_background_ra

        # Compute the Poisson error of the sources...
        total_error = calc_total_error(imgarr_bkgsub, self.image.bkg_rms_ra, 1.0)

        # Compute source properties...
        if OLD_PHOTUTILS:
            self.source_cat = source_properties(imgarr_bkgsub, self.sources, background=self.image.bkg_background_ra,
                                                error=total_error, filter_kernel=self.kernel, wcs=self.image.imgwcs)
        else:
            self.source_cat = SourceCatalog(imgarr_bkgsub, self.sources, background=self.image.bkg_background_ra,
                                                error=total_error, kernel=self.kernel, wcs=self.image.imgwcs)

        # Convert source_cat which is a SourceCatalog to an Astropy Table
        filter_measurements_table = Table(self.source_cat.to_table(columns=self.include_filter_cols))

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
        rr = Column(ra_icrs, name="RA", unit="degrees")
        dd = Column(dec_icrs, name="DEC", unit="degrees")
        filter_measurements_table.add_columns([dd, rr])

        # Compute the MagIso
        filter_measurements_table["MagIso"] = photometry_tools.convert_flux_to_abmag(filter_measurements_table[flux_colname],
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
                ci_data = mag_inner_data - mag_outer_data
                ci_mask = np.logical_and(np.abs(ci_data) > 0.0, np.abs(ci_data) < 1.0e-30)
                big_bad_index = np.where(abs(ci_data) > 1.0e20)
                ci_mask[big_bad_index] = True
                ci_col = MaskedColumn(name='CI', data=ci_data, dtype=np.float64, mask=ci_mask)

            except Exception as x_cept:
                log.warning("Computation of concentration index (CI) was not successful: {} - {}.".format(self.imgname, x_cept))
                log.warning("CI measurements may be missing from the output filter catalog.\n")

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

        # iso_col = Column(data=np.ones(tblLen)*float('nan')), name="MagIso")
        # table.add_column(iso_col)

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
                           bac_colname: "Bck", flux_colname: "FluxIso", ferr_colname: "FluxIsoErr",
                           "bbox_xmin": "Xmin", "bbox_ymin": "Ymin", "bbox_xmax": "Xmax", "bbox_ymax": "Ymax",
                           "cxx": "CXX", "cyy": "CYY", "cxy": "CXY",
                           "covar_sigx2": "X2", "covar_sigy2": "Y2", "covar_sigxy": "XY",
                           "orientation": "Theta",
                           "elongation": "Elongation", "ellipticity": "Ellipticity", "area": "Area"}
        for old_col_title in final_col_names:
            filter_table.rename_column(old_col_title, final_col_names[old_col_title])

        # Define the order of the columns
        final_col_order = ["X-Centroid", "Y-Centroid", "RA", "DEC", "ID",
                           "CI", "Flags", "MagAp1", "MagErrAp1", "FluxAp1", "FluxErrAp1",
                           "MagAp2", "MagErrAp2", "FluxAp2", "FluxErrAp2", "MSkyAp2",
                           "Bck", "Area", "MagIso", "FluxIso", "FluxIsoErr",
                           "Xmin", "Ymin", "Xmax", "Ymax",
                           "X2", "Y2", "XY",
                           "CXX", "CYY", "CXY",
                           "Elongation", "Ellipticity", "Theta"]
        final_filter_table = filter_table[final_col_order]

        # Define the format
        final_col_format = {"X-Centroid": "10.3f", "Y-Centroid": "10.3f", "RA": "13.7f", "DEC": "13.7f", "ID": "7d",
                            "CI": "7.3f", "Flags": "5d",
                            "MagAp1": "7.3f", "MagErrAp1": "7.3f", "FluxAp1": "10.4f", "FluxErrAp1": "10.4f",
                            "MagAp2": "7.3f", "MagErrAp2": "7.3f", "FluxAp2": "10.4f", "FluxErrAp2": "10.4f",
                            "MSkyAp2": "7.3f", "Bck": "10.4f", "MagIso": "7.3f", "FluxIso": "10.4f", "FluxIsoErr": "10.4f",
                            "Xmin": "8.0f", "Ymin": "8.0f", "Xmax": "8.0f", "Ymax": "8.0f",
                            "X2": "8.4f", "Y2": "8.4f", "XY": "8.4f",
                            "CXX": "9.5f", "CYY": "9.5f", "CXY": "9.5f",
                            "Elongation": "7.2f", "Ellipticity": "7.2f", "Theta": "8.3f", "Area": "8.3f"}
        for fcf_key in final_col_format.keys():
            final_filter_table[fcf_key].format = final_col_format[fcf_key]

        # Add description
        final_col_descrip = {"ID": "Catalog Object Identification Number",
                             "X-Centroid": "Pixel Coordinate", "Y-Centroid": "Pixel Coordinate",
                             "RA": "Sky coordinate at epoch of observation",
                             "DEC": "Sky coordinate at epoch of observation",
                             "Bck": "Background at the position of the source centroid",
                             "Area": "Total unmasked area of the source segment",
                             "MagAp1": "ABMAG of source based on the inner (smaller) aperture",
                             "MagErrAp1": "Error of MagAp1",
                             "FluxAp1": "Flux of source based on the inner (smaller) aperture",
                             "FluxErrAp1": "Error of FluxAp1",
                             "MagAp2": "ABMAG of source based on the outer (larger) aperture",
                             "MagErrAp2": "Error of MagAp2",
                             "FluxAp2": "Flux of source based on the outer (larger) aperture",
                             "FluxErrAp2": "Error of FluxAp2",
                             "MSkyAp2": "ABMAG of sky based on outer (larger) aperture",
                             "FluxIso": "Sum of unmasked data values in the source segment",
                             "FluxIsoErr": "Uncertainty of FluxIso, propagated from the input error array",
                             "MagIso": "Magnitude corresponding to FluxIso",
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
        final_col_unit = {"X-Centroid": "pixels", "Y-Centroid": "pixels",
                          "RA": "degrees", "DEC": "degrees",
                          "Bck": "electrons/s",
                          "Area": "pixels**2",
                          "MagAp1": "ABMAG",
                          "MagErrAp1": "ABMAG",
                          "FluxAp1": "electrons/s",
                          "FluxErrAp1": "electrons/s",
                          "MagAp2": "ABMAG",
                          "MagErrAp2": "ABMAG",
                          "FluxAp2": "electrons/s",
                          "FluxErrAp2": "electrons/s",
                          "MSkyAp2": "ABMAG",
                          "MagIso": "ABMAG",
                          "FluxIso": "electrons/s",
                          "FluxIsoErr": "electrons/s",
                          "X2": "pixels**2",
                          "Y2": "pixels**2",
                          "XY": "pixels**2",
                          "CXX": "pixels**2",
                          "CYY": "pixels**2",
                          "CXY": "pixels**2",
                          "Xmin": "pixels",
                          "Ymin": "pixels",
                          "Xmax": "pixels",
                          "Ymax": "pixels",
                          "Theta": "radians",
                          "CI": "ABMAG",
                          "Flags": "unitless",
                          "ID": "unitless",
                          "Elongation": "unitless",
                          "Ellipticity": "unitless"}
        for fcu_key in final_col_unit.keys():
            final_filter_table[fcu_key].unit = final_col_unit[fcu_key]

        return(final_filter_table)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def _define_empty_table(self, segm_img):
        """Create basic empty table based on total_table format to signify no valid sources were found"""

        final_col_unit = {"X-Centroid": "pixels", "Y-Centroid": "pixels",
                          "RA": "degrees", "DEC": "degrees", "Flags": "unitless"}
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
        rr = Column(ra_icrs, name="RA", unit="degrees")
        dd = Column(dec_icrs, name="DEC", unit="degrees")
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
        final_col_unit = {"X-Centroid": "pixels", "Y-Centroid": "pixels",
                          "RA": "degrees", "DEC": "degrees"}
        for fcu_key in final_col_unit.keys():
            table[fcu_key].unit = final_col_unit[fcu_key]

        return(table)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def _evaluate_segmentation_image(self, segm_img, image_data, big_island_only=False,
                                     max_biggest_source=0.015, max_source_fraction=0.075):
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

        n, binedges = np.histogram(segm_img.data, range=(1, nbins))
        real_pixels = (image_data != 0).sum()

        # Compute which segments are larger than the kernel.
        deb_limit = self.kernel.size
        log.debug("Deblending limit set at: {} pixels".format(deb_limit))
        # add as attribute to SegmentationImage for use later
        segm_img.big_segments = np.where(segm_img.areas >= deb_limit)[0] + 1  # Segment labels are 1-based

        biggest_source = n.max()/float(real_pixels)
        log.info("Biggest_source: %f", biggest_source)
        if biggest_source > max_biggest_source:
            log.info("Biggest source %.4f percent exceeds %f percent of the image", (100.0*biggest_source), (100.0*max_biggest_source))
            is_poor_quality = True

        # Always compute the source_fraction so the value can be reported.  Setting the
        # big_island_only parameter allows control over whether the source_fraction should
        # or should not be ignored.
        source_fraction = n.sum()/float(real_pixels)
        log.info("Source_fraction: %f", source_fraction)
        if not big_island_only:
            if source_fraction > max_source_fraction:
                log.info("Total source fraction %.4f percent exceeds %f percent of the image.", (100.0*source_fraction), (100.0*max_source_fraction))
                is_poor_quality = True

        return is_poor_quality, biggest_source, source_fraction

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def write_catalog(self, reject_catalogs):
        """Write the specified source catalog out to disk.

        Regardless of the setting for reject_catalogs, the regions file will be written
        solely based upon the setting of diagnostic_mode.

        Parameters
        ----------
        reject_catalogs : bool
            Indicator as to whether or not the catalogs (*.ecsv) should be written.

        Returns
        -------
        Nothing

        """
        self.source_cat = self.annotate_table(self.source_cat, self.param_dict_qc, proc_type="segment",
                                              product=self.image.ghd_product)
        if reject_catalogs:
            # We still want to write out empty files
            # This will delete all rows from the existing table
            self.source_cat.remove_rows(slice(0, None))

        # Write out catalog to ecsv file
        self.source_cat.write(self.sourcelist_filename, format=self.catalog_format)
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


def make_wht_masks(whtarr, maskarr, scale=1.5, sensitivity=0.95, kernel=(11, 11)):

    invmask = make_inv_mask(maskarr)

    maxwht = ndimage.filters.maximum_filter(whtarr, size=kernel)
    rel_wht = maxwht / maxwht.max()

    delta = 0.0
    master_mask = np.zeros(invmask.shape, dtype=np.uint16)
    limit = 1 / scale
    masks = []
    while delta < sensitivity:

        mask = rel_wht > limit
        mask = (mask.astype(np.uint16) * invmask) - master_mask

        new_delta = master_mask.sum() / mask.sum()
        if new_delta < sensitivity:
            masks.append(dict(scale=limit,
                              wht_limit=limit * maxwht.max(),
                              mask=mask,
                              rel_weight=rel_wht * mask))

        delta = new_delta
        master_mask = master_mask + mask
        limit /= scale

    return masks
