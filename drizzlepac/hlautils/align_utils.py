import os
import sys
import shutil
import datetime
import traceback
import glob

import numpy as np
from scipy import ndimage

import astropy
from astropy.io import fits
from astropy.table import Table
from astropy.nddata.bitmask import bitfield_to_boolean_mask

from photutils import Background2D, SExtractorBackground, StdBackgroundRMS

from stwcs.wcsutil import HSTWCS
from stsci.tools import logutil

from . import astrometric_utils as amutils
from . import astroquery_utils as aqutils

from .. import wcs_functions

__taskname__ = 'align_utils'

# Default background determination parameter values
BKG_BOX_SIZE = 50
BKG_FILTER_SIZE = 3
CATALOG_TYPES = ['point', 'segment']

log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout)


class AlignmentTable:
    def __init__(self, input_list, clobber=False):
        self.zero_dt = starting_dt = datetime.datetime.now()
        log.info(str(starting_dt))
        # Apply filter to input observations to insure that they meet minimum criteria for being able to be aligned
        log.info(
            "{} AlignmentTable: Filter STEP {}".format("-" * 20, "-" * 63))
        self.filtered_table = filter.analyze_data(input_list)

        if self.filtered_table['doProcess'].sum() == 0:
            log.warning("No viable images in filtered table - no processing done.\n")
            current_dt = datetime.datetime.now()
            delta_dt = (current_dt - starting_dt).total_seconds()
            log.info('Processing time of AlignmentTable Filter STEP: {} sec'.format(delta_dt))
            return

        # Get the list of all "good" files to use for the alignment
        process_list = self.filtered_table['imageName'][np.where(self.filtered_table['doProcess'])]
        self.process_list = list(process_list)  # Convert process_list from numpy list to regular python list
        log.info("SUCCESS")

        # Initialize computed attributes
        self.imglist = []

    def build_images(self):
        """Create CatalogImage objects for each input to use in alignment."""
        if self.imglist:
            return
        for img in self.process_list:
            catimg = CatalogImage(img)
            # Build image properties needed for alignment
            catimg.build_kernel()

            self.imglist.append(catimg)


    def set_catalog(self, catalog_names):
        """Define the reference catalog(s) to be used for alignment

        Parameters
        ----------
        catalog_names : list
            List of astrometric catalog names to use for alignment.
            Options would include (but not be limited to):
            "GAIADR1" and "GAIADR2".
            If a filename is provided, file would be used instead
            of deriving catalogs from standard astrometric catalogs
            like GAIADR2.
        """
        pass

    def set_method(self, method):
        """Define what alignment method(s) are to be used

        Parameters
        ----------
        method : list
            List of alignment method names to be used.
            Supported options: relative, 2dhist, threshold.
        """
        pass

    def find_sources(self):
        """Find observable sources in each input exposure."""
        self.extracted_sources = None

    def get_reference_catalog(self, catalog):
        """Return the desired reference catalog to be used for alignment"""
        pass

    def reset_group_id(self):
        for image in self.imglist:
            image.meta["group_id"] = self.group_id_dict["{}_{}".format(image.meta["filename"], image.meta["chip"])]

    def configure_fit(self):
        # Convert input images to tweakwcs-compatible FITSWCS objects and
        # attach source catalogs to them.

        self.imglist = []
        for group_id, image in enumerate(self.process_list):
            img = amutils.build_wcscat(image, group_id,
                                       self.extracted_sources[image]['catalog_table'])
            # add the name of the image to the imglist object
            for im in img:
            #    im.meta['name'] = image
                log.info('im.meta[name] = {}'.format(im.meta['name']))
            self.imglist.extend(img)

        self.group_id_dict = {}
        for image in self.imglist:
            self.group_id_dict["{}_{}".format(image.meta["filename"], image.meta["chip"])] = image.meta["group_id"]

    def perform_fit(self, method):
        """Perform fit using specified method, then determine fit quality"""
        pass

# Including this from catalog_utils
#
# TODO:  Merge differences into version in catalog_utils
#
class CatalogImage:
    def __init__(self, filename):
        if isinstance(filename, str):
            self.imghdu = fits.open(filename)
            self.imgname = filename
        else:
            self.imghdu = filename
            self.imgname = filename.filename()

        # Get header information to annotate the output catalogs
        if "total" in self.imgname:
            self.ghd_product = "tdp"
        else:
            self.ghd_product = "fdp"  # appropriate for single exposures, too

        # Fits file read
        self.data = self.imghdu[('SCI', 1)].data
        self.num_sci = amutils.countExtn(self.imghdu)

        # Get the HSTWCS object from the first extension
        self.imgwcs = HSTWCS(self.imghdu, 1)
        self.pscale = self.imgwcs.pscale

        self.keyword_dict = self._get_header_data()

        self.bkg = None
        self.imglist = None
        self.catalog_table = {}

    def close(self):
        self.imghdu.close()

    def build_kernel(self, fwhmpsf):
        if self.bkg is None:
            self.compute_background()

        self.kernel, self.kernel_fwhm = amutils.build_auto_kernel(self.data, self.wht_image,
                                                          threshold=self.bkg.background_rms, fwhm=fwhmpsf / self.pscale)

    def compute_background(self, box_size=BKG_BOX_SIZE, win_size=BKG_FILTER_SIZE,
                           bkg_estimator=SExtractorBackground, rms_estimator=StdBackgroundRMS,
                           nsigma=5., threshold_flag=None):
        """Use Background2D to determine the background of the input image.
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
        nsigma : float
            Number of sigma above background
        threshold_flag : float or None
            Value from the image which serves as the limit for determining sources.
            If None, compute a default value of (background+5*rms(background)).
            If threshold < 0.0, use absolute value as scaling factor for default value.
        Returns
        -------
        bkg : 2D ndarray
            Background image
        bkg_dao_rms : ndarry
            Background RMS image
        threshold : ndarray
            Numpy array representing the background plus RMS
        """
        # Report configuration values to log
        log.info("")
        log.info("Computation of image background - Input Parameters")
        log.info("Box size: {}".format(box_size))
        log.info("Window size: {}".format(win_size))
        log.info("NSigma: {}".format(nsigma))

        # SExtractorBackground ans StdBackgroundRMS are the defaults
        bkg = None
        bkg_dao_rms = None

        exclude_percentiles = [10, 25, 50, 75]
        for percentile in exclude_percentiles:
            log.info("")
            log.info("Percentile in use: {}".format(percentile))
            try:
                bkg = Background2D(self.data, box_size, filter_size=win_size,
                                   bkg_estimator=bkg_estimator(),
                                   bkgrms_estimator=rms_estimator(),
                                   exclude_percentile=percentile, edge_method="pad")
            except Exception:
                bkg = None
                continue

            if bkg is not None:
                # Set the bkg_rms at "nsigma" sigma above background
                bkg_rms = nsigma * bkg.background_rms
                default_threshold = bkg.background + bkg_rms
                bkg_rms_mean = bkg.background.mean() + nsigma * bkg_rms.std()
                bkg_mean = bkg.background.mean()
                bkg_dao_rms = bkg.background_rms
                if threshold_flag is None:
                    threshold = default_threshold
                elif threshold_flag < 0:
                    threshold = -1 * threshold_flag * default_threshold
                    log.info("Background threshold set to {} based on {}".format(threshold.max(), default_threshold.max()))
                    bkg_rms_mean = threshold.max()
                else:
                    bkg_rms_mean = 3. * threshold_flag
                    threshold = bkg_rms_mean

                if bkg_rms_mean < 0:
                    bkg_rms_mean = 0.
                break

        # If Background2D does not work at all, define default scalar values for
        # the background to be used in source identification
        if bkg is None:
            bkg_mean = bkg_rms_mean = max(0.01, self.data.min())
            bkg_rms = nsigma * bkg_rms_mean
            bkg_dao_rms = bkg_rms_mean
            threshold = bkg_rms_mean + bkg_rms

        # *** FIX: Need to do something for bkg if bkg is None ***

        # Report other useful quantities
        log.info("")
        log.info("Mean background: {}".format(bkg_mean))
        log.info("Mean threshold: {}".format(np.mean(threshold)))
        log.info("")
        log.info("{}".format("=" * 80))

        self.bkg = bkg
        self.bkg_dao_rms = bkg_dao_rms
        self.bkg_rms_mean = bkg_rms_mean
        self.threshold = threshold

    def find_sources(self, output=True, dqname='DQ', fwhmpsf=3.0, **alignment_pars):
        """Find sources in all chips for this exposure."""
        for chip in range(self.numSci):
            chip += 1
            # find sources in image
            if output:
                outroot = '{}_sci{}_src'.format(self.keyword_dict['rootname'], chip)

            # apply any DQ array, if available
            dqmask = None
            if self.imghdu.index_of(dqname):
                dqarr = self.imghdu[dqname, chip].data

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

            if not self.kernel:
                imgarr = self.imghdu['sci', chip].data
                num_wht = amutils.countExtn(self.imghdu, extn='WHT')
                if num_wht > 0:
                    wht_image = self.imghdu['WHT'].data.copy()
                else:
                    # Build pseudo-wht array for detection purposes
                    errarr = self.imghdu['err', chip].data
                    wht_image = errarr / errarr.max()
                    wht_image[dqmask] = 0

                self.kernel, self.kernel_fwhm = amutils.build_auto_kernel(imgarr, wht_image,
                                                          threshold=self.bkg.background_rms,
                                                          fwhm=fwhmpsf / self.pscale)

            #  TODO: replace detector_pars with dict from OO Config class
            seg_tab, segmap = amutils.extract_sources(imgarr, dqmask=dqmask, outroot=outroot,
                                                      kernel=self.kernel,
                                                      segment_threshold=self.threshold,
                                                      dao_threshold=self.bkg_rms_mean,
                                                      fwhm=self.kernel_fwhm,
                                                      **alignment_pars)

            self.catalog_table[chip] = seg_tab


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
        keyword_dict["instrument"] = self.imghdu[0].header["INSTRUME"].upper()
        keyword_dict["detector"] = self.imghdu[0].header["DETECTOR"].upper()
        keyword_dict["target_ra"] = self.imghdu[0].header["RA_TARG"]
        keyword_dict["target_dec"] = self.imghdu[0].header["DEC_TARG"]
        keyword_dict["expo_start"] = self.imghdu[0].header["EXPSTART"]
        keyword_dict["texpo_time"] = self.imghdu[0].header["TEXPTIME"]
        keyword_dict["ccd_gain"] = self.imghdu[0].header["CCDGAIN"]
        keyword_dict["aperture_pa"] = self.imghdu[0].header["PA_V3"]
        keyword_dict["rootname"] = self.imghdu[0].header["ROOTNAME"]

        # The total detection product has the FILTER keyword in
        # the primary header - read it for any instrument.
        #
        # For the filter detection product:
        # WFC3 only has FILTER, but ACS has FILTER1 and FILTER2
        # in the primary header.
        if self.ghd_product.lower() == "tdp":
            keyword_dict["filter"] = self.imghdu[0].header["FILTER"]
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
        log.info('WCSTYPE: {}'.format(keyword_dict["wcs_type"]))
        keyword_dict["orientation"] = self.imghdu[1].header["ORIENTAT"]
        keyword_dict["aperture_ra"] = self.imghdu[1].header["RA_APER"]
        keyword_dict["aperture_dec"] = self.imghdu[1].header["DEC_APER"]

        return keyword_dict
