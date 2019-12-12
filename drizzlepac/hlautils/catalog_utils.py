"""This script contains code to support creation of photometric sourcelists using two techniques: aperture photometry
segmentation-map based photometry.
"""
import sys
import pickle  # FIX Remove
import copy

import astropy.units as u
from astropy.io import fits as fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Column, MaskedColumn, Table, join
from astropy.coordinates import SkyCoord
import numpy as np
from scipy import ndimage

from photutils import CircularAperture, CircularAnnulus, DAOStarFinder
from photutils import Background2D, SExtractorBackground, StdBackgroundRMS
from photutils import detect_sources, source_properties, deblend_sources
from photutils import make_source_mask
from stsci.tools import logutil
from stwcs.wcsutil import HSTWCS

from . import astrometric_utils
from . import photometry_tools

try:
    from matplotlib import pyplot as plt
except Exception:
    plt = None

CATALOG_TYPES = ['aperture', 'segment']

__taskname__ = 'catalog_utils'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)


class CatalogImage:
    def __init__(self, filename, log_level):
        # set logging level to user-specified level
        log.setLevel(log_level)

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
            self.ghd_product = "fdp"

        # Fits file read
        self.data = self.imghdu[('SCI', 1)].data
        self.wht_image = self.imghdu['WHT'].data.copy()

        # Get the HSTWCS object from the first extension
        self.imgwcs = HSTWCS(self.imghdu, 1)

        self.keyword_dict = self._get_header_data()

        # Populated by self.compute_background()
        self.bkg = None
        self.bkg_background_ra = None
        self.bkg_rms_ra = None
        self.bkg_rms_median = None
        self.bkg_median = None

        # Populated by self.build_kernel()
        self.kernel = None
        self.kernel_fwhm = None
        self.kernel_psf = False

    def close(self):
        self.imghdu.close()

    # def build_kernel(self, box_size, win_size, fwhmpsf):
    def build_kernel(self, box_size, win_size, fwhmpsf):
        if self.bkg is None:
            self.compute_background(box_size, win_size)

        k, self.kernel_fwhm = astrometric_utils.build_auto_kernel(self.data,
                                                                  self.wht_image,
                                                                  threshold=self.bkg_rms_ra,
                                                                  fwhm=fwhmpsf / self.imgwcs.pscale)
        (self.kernel, self.kernel_psf) = k

    def compute_background(self, box_size, win_size,
                           bkg_estimator=SExtractorBackground, rms_estimator=StdBackgroundRMS):
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

        Returns
        -------
        bkg_background_ra : 2D ndarray
            Background array

        bkg_rms_ra : 2D ndarray
            RMS map array

        bkg_rms_median : float
            bkg_rms_image median value
        """
        # Report configuration values to log
        log.info("")
        log.info("Computation of image background - Input Parameters")
        log.info("File: {}".format(self.imgname))
        log.info("Box size: {}".format(box_size))
        log.info("Window size: {}".format(win_size))

        # SExtractorBackground ans StdBackgroundRMS are the defaults
        bkg = None

        exclude_percentiles = [10, 25, 50, 75]
        for percentile in exclude_percentiles:
            log.info("Percentile in use: {}".format(percentile))
            try:
                bkg = Background2D(self.data, (box_size, box_size), filter_size=(win_size, win_size),
                                   bkg_estimator=bkg_estimator(),
                                   bkgrms_estimator=rms_estimator(),
                                   exclude_percentile=percentile, edge_method="pad",
                                   mask=(self.data == 0))

            except Exception:
                bkg = None
                continue

            if bkg is not None:
                bkg_background_ra = bkg.background
                bkg_rms_ra = bkg.background_rms
                bkg_rms_median = bkg.background_rms_median
                bkg_median = bkg.background_median
                break

        # If Background2D does not work at all, define default scalar values for
        # the background to be used in source identification
        if bkg is None:
            log.info("Background2D failure detected. Using alternative background calculation instead....")
            mask = make_source_mask(self.data, nsigma=2, npixels=5, dilate_size=11)
            sigcl_mean, sigcl_median, sigcl_std = sigma_clipped_stats(self.data, sigma=3.0, mask=mask, maxiters=9)
            bkg_median = sigcl_median
            bkg_rms_median = sigcl_std
            # create background frame shaped like self.data populated with sigma-clipped median value
            bkg_background_ra = np.full_like(self.data, sigcl_median)
            # create background frame shaped like self.data populated with sigma-clipped standard deviation value
            bkg_rms_ra = np.full_like(self.data, sigcl_std)

        log.info("Computation of image background complete")
        log.info("Found: ")
        log.info("    Median background: {}".format(bkg_median))
        log.info("    Median RMS background: {}".format(bkg_rms_median))
        log.info("")

        self.bkg = bkg
        self.bkg_background_ra = bkg_background_ra
        self.bkg_rms_ra = bkg_rms_ra
        self.bkg_rms_median = bkg_rms_median
        self.bkg_median = bkg_median

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
        keyword_dict["ccd_gain"] = self.imghdu[0].header["CCDGAIN"]
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
            else:
                keyword_dict["filter1"] = self.imghdu[0].header["FILTER"]
                keyword_dict["filter2"] = ""

        # Get the HSTWCS object from the first extension
        keyword_dict["wcs_name"] = self.imghdu[1].header["WCSNAME"]
        keyword_dict["wcs_type"] = self.imghdu[1].header["WCSTYPE"]
        keyword_dict["orientation"] = self.imghdu[1].header["ORIENTAT"]
        keyword_dict["aperture_ra"] = self.imghdu[1].header["RA_APER"]
        keyword_dict["aperture_dec"] = self.imghdu[1].header["DEC_APER"]

        return keyword_dict


class HAPCatalogs:
    """Generate photometric sourcelist for specified TOTAL or FILTER product image.
    """

    def __init__(self, fitsfile, param_dict, param_dict_qc, log_level, diagnostic_mode=False, types=None,
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

        # Compute the background for this image
        self.image = CatalogImage(fitsfile, log_level)
        self.image.compute_background(self.param_dict['bkg_box_size'], self.param_dict['bkg_filter_size'])

        self.image.build_kernel(self.param_dict['bkg_box_size'], self.param_dict['bkg_filter_size'],
                                self.param_dict['dao']['TWEAK_FWHMPSF'])

        # Initialize all catalog types here...
        # This does NOT identify or measure sources to create the catalogs at this point...
        # The syntax here is EXTREMELY cludgy, but until a more compact way to do this is found,
        #  it will have to do...
        self.catalogs = {}
        if 'aperture' in self.types:
            self.catalogs['aperture'] = HAPPointCatalog(self.image, self.param_dict, self.param_dict_qc,
                                                        self.diagnostic_mode, tp_sources=tp_sources)
        if 'segment' in self.types:
            self.catalogs['segment'] = HAPSegmentCatalog(self.image, self.param_dict, self.param_dict_qc,
                                                         self.diagnostic_mode, tp_sources=tp_sources)

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

    def write(self, **pars):
        """Write catalogs for this image to output files.

        Parameters
        ----------
        types : list
            List of catalog types to be generated.  If None, build all available catalogs.
            Supported types of catalogs include: 'aperture', 'segment'.
        """
        # Make sure we at least have a default 2D background computed

        for catalog in self.catalogs.values():
            if catalog.source_cat is None:
                catalog.source_cat = catalog.sources
            catalog.write_catalog

    def combine(self, subset_dict):
        """Combine subset columns from the filter catalog with the total detection catalog.

        Parameters
        ----------
        subset_dict: dictionary
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

    def __init__(self, image, param_dict, param_dict_qc, diagnostic_mode, tp_sources):
        self.image = image
        self.imgname = image.imgname
        self.bkg = image.bkg
        self.param_dict = param_dict
        self.param_dict_qc = param_dict_qc
        self.diagnostic_mode = diagnostic_mode

        self.sourcelist_filename = self.imgname.replace(self.imgname[-9:], self.catalog_suffix)

        # Compute average gain - there will always be at least one gain value in the primary header
        gain_keys = self.image.imghdu[0].header['atodgn*']
        gain_values = [gain_keys[g] for g in gain_keys if gain_keys[g] > 0.0]
        self.gain = self.image.imghdu[0].header['exptime'] * np.mean(gain_values)

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
            log.info("self.imgname:   {}".format(self.imgname))
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

    def identify_sources(self, **pars):
        pass

    def measure_sources(self, filter_name, **pars):
        pass

    def write_catalog(self, **pars):
        pass

    def combine_tables(self, subset_dict):
        pass

    def annotate_table(self, data_table, param_dict_qc, product="tdp"):
        """Add state metadata to the top of the output source catalog.

        Parameters
        ----------
        data_table : QTable
            Table of source properties

        param_dict_qc : dictionary
            Configuration values for quality control step based upon input JSON files (used to build catalog header)

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
        data_table.meta["CCD Gain"] = self.image.keyword_dict["ccd_gain"]
        if product.lower() == "tdp" or self.image.keyword_dict["instrument"].upper() == "WFC3":
            data_table.meta["Filter 1"] = self.image.keyword_dict["filter1"]
            data_table.meta["Filter 2"] = ""
        else:
            data_table.meta["Filter 1"] = self.image.keyword_dict["filter1"]
            data_table.meta["Filter 2"] = self.image.keyword_dict["filter2"]
        num_sources = len(data_table)
        data_table.meta["Number of sources"] = num_sources

        if "X-Center" in data_table.colnames:
            proc_type = "aperture"
        else:
            proc_type = "segment"
        ci_lower = float(param_dict_qc['ci filter'][proc_type]['ci_lower_limit'])
        ci_upper = float(param_dict_qc['ci filter'][proc_type]['ci_upper_limit'])

        data_table.meta["h09"] = ["#================================================================================================="]
        data_table.meta["h10"] = ["IMPORTANT NOTES"]
        data_table.meta["h11"] = ["The X and Y coordinates in this table are 0-indexed (i.e. the origin is (0,0))."]
        data_table.meta["h12"] = ["RA and Dec values in this table are in sky coordinates (i.e. coordinates at the epoch of observation"]
        data_table.meta["h12.1"] = ["and fit to GAIADR1 (2015.0) or GAIADR2 (2015.5))."]
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
        data_table.meta["h17.3"] = ["    8 - Faint Detection Limit"]
        data_table.meta["h17.4"] = ["   16 - Hot Pixels (CI < {})".format(ci_lower)]
        data_table.meta["h17.5"] = ["   32 - False Detection Swarm Around Saturated Source"]
        data_table.meta["h17.6"] = ["   64 - False Detections Near Image Edge"]
        data_table.meta["h18"] = ["#================================================================================================="]

        return (data_table)


# --------------------------------------------------------------------------------------------------------

class HAPPointCatalog(HAPCatalogBase):
    """Generate photometric sourcelist(s) for specified image(s) using aperture photometry of point sources.
    """
    catalog_suffix = "_point-cat.ecsv"

    def __init__(self, image, param_dict, param_dict_qc, diagnostic_mode, tp_sources):
        super().__init__(image, param_dict, param_dict_qc, diagnostic_mode, tp_sources)

        self.bkg_used = None  # actual background used for source identification/measurement

        # Defined in measure_sources
        self.subset_filter_source_cat = None

    def identify_sources(self, **pars):
        """Create a master coordinate list of sources identified in the specified total detection product image
        """
        source_fwhm = self.image.kernel_fwhm
        # read in sci, wht extensions of drizzled product
        image = self.image.data.copy()

        # Input image will be background subtracted using pre-computed background, unless
        # specified explicitly by the user
        if self.param_dict["dao"]["simple_bkg"]:
            self.bkg_used = np.nanmedian(image)
            image -= self.bkg_used
        else:
            # Estimate background
            self.bkg_used = self.image.bkg_background_ra
            image -= self.bkg_used

        if not self.tp_sources:
            # Report configuration values to log
            log.info("{}".format("=" * 80))
            log.info("")
            log.info("Point-source finding settings")
            log.info("Total Detection Product - Input Parameters")
            log.info("INPUT PARAMETERS")
            log.info("{}: {}".format("self.param_dict['dao']['bkgsig_sf']", self.param_dict["dao"]["bkgsig_sf"]))
            log.info("{}: {}".format("self.param_dict['dao']['kernel_sd_aspect_ratio']",
                                     self.param_dict['dao']['kernel_sd_aspect_ratio']))
            log.info("{}: {}".format("self.param_dict['dao']['simple_bkg']", self.param_dict['dao']['simple_bkg']))
            log.info("{}: {}".format("self.param_dict['nsigma']", self.param_dict['nsigma']))
            log.info("{}: {}".format("self.image.bkg_rms_median", self.image.bkg_rms_median))
            log.info("DERIVED PARAMETERS")
            log.info("{}: {}".format("source_fwhm", source_fwhm))
            log.info("{}: {}".format("threshold", self.param_dict['nsigma'] * self.image.bkg_rms_median))
            log.info("")
            log.info("{}".format("=" * 80))

            # find ALL the sources!!!
            log.info("DAOStarFinder(fwhm={}, threshold={}*{})".format(source_fwhm, self.param_dict['nsigma'],
                                                                      self.image.bkg_rms_median))
            log.info("{}".format("=" * 80))

            daofind = DAOStarFinder(fwhm=source_fwhm,
                                    threshold=self.param_dict['nsigma'] * self.image.bkg_rms_median)

            # create mask to reject any sources located less than 10 pixels from a image/chip edge
            wht_image = self.image.data.copy()
            binary_inverted_wht = np.where(wht_image == 0, 1, 0)
            exclusion_mask = ndimage.binary_dilation(binary_inverted_wht, iterations=10)

            sources = daofind(image, mask=exclusion_mask)

            for col in sources.colnames:
                sources[col].info.format = '%.8g'  # for consistent table output

            self.sources = sources

        # if processing filter product, use sources identified by parent total drizzle product identify_sources() run
        if self.tp_sources:
            self.sources = self.tp_sources['aperture']['sources']

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def measure_sources(self, filter_name):
        """Perform aperture photometry on identified sources
        """
        log.info("Performing aperture photometry on identified point-sources")
        # Open and background subtract image
        image = self.image.data.copy()
        image -= self.bkg_used

        # load in coords of sources identified in total product
        try:
            positions = (self.sources['xcentroid'], self.sources['ycentroid'])
        except Exception:
            positions = (self.sources['X-Center'], self.sources['Y-Center'])

        pos_xy = np.vstack(positions).T

        # define list of background annulii
        bg_apers = CircularAnnulus(pos_xy,
                                   r_in=self.param_dict['skyannulus_arcsec'],
                                   r_out=self.param_dict['skyannulus_arcsec'] +
                                   self.param_dict['dskyannulus_arcsec'])

        # Create the list of photometric apertures to measure
        phot_apers = [CircularAperture(pos_xy, r=r) for r in self.aper_radius_list_pixels]
        # Perform aperture photometry
        photometry_tbl = photometry_tools.iraf_style_photometry(phot_apers,
                                                                bg_apers,
                                                                data=image,
                                                                photflam=self.image.imghdu[1].header['photflam'],
                                                                photplam=self.image.imghdu[1].header['photplam'],
                                                                error_array=self.image.bkg_rms_ra,
                                                                bg_method=self.param_dict['salgorithm'],
                                                                epadu=self.gain)

        # calculate and add RA and DEC columns to table
        ra, dec = self.transform_list_xy_to_ra_dec(photometry_tbl["X-Center"], photometry_tbl["Y-Center"], self.imgname)  # TODO: replace with all_pix2sky or somthing at a later date
        ra_col = Column(name="RA", data=ra, dtype=np.float64)
        dec_col = Column(name="DEC", data=dec, dtype=np.float64)
        photometry_tbl.add_column(ra_col, index=2)
        photometry_tbl.add_column(dec_col, index=3)

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
        final_col_format = {"X-Center": "10.3f", "Y-Center": "10.3f", "RA": "13.10f", "DEC": "13.10f", "ID": ".8g", "MagAp1": '6.3f', "MagErrAp1": '6.3f', "MagAp2": '6.3f',
                            "MagErrAp2": '6.3f', "MSkyAp2": '10.8f', "StdevAp2": '10.4f',
                            "FluxAp2": '10.8f', "CI": "7.3f", "Flags": "3d"}  # TODO: Standardize precision
        for fcf_key in final_col_format.keys():
            output_photometry_table[fcf_key].format = final_col_format[fcf_key]

        # add units to columns
        final_col_units = {"X-Center": "Pixels", "Y-Center": "Pixels", "RA": "Sky Coords", "DEC": "Sky Coords",
                           "ID": "Unitless", "MagAp1": "ABMAG", "MagErrAp1": "ABMAG", "MagAp2": "ABMAG",
                           "MagErrAp2": "ABMAG", "MSkyAp2": "ABMAG", "StdevAp2": "ABMAG",
                           "FluxAp2": "electrons/sec", "CI": "ABMAG", "Flags": "Unitless"}
        for col_title in final_col_units:
            output_photometry_table[col_title].unit = final_col_units[col_title]

        # Capture specified columns in order to append to the total detection table
        self.subset_filter_source_cat = output_photometry_table["ID", "RA", "DEC", "MagAp2", "CI", "Flags"]
        self.subset_filter_source_cat.rename_column("MagAp2", "MagAP2_" + filter_name)
        self.subset_filter_source_cat.rename_column("CI", "CI_" + filter_name)
        self.subset_filter_source_cat.rename_column("Flags", "Flags_" + filter_name)

        # Add the header information to the table
        self.source_cat = self.annotate_table(output_photometry_table,
                                              self.param_dict_qc,
                                              product=self.image.ghd_product)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    @property
    def write_catalog(self):
        """Write specified catalog to file on disk

        Parameters
        ----------
        write_region_file : Boolean
           Write ds9-compatible region file along with the catalog file? Default value = False

        Returns
        -------
        Nothing!

        """
        # Write out catalog to ecsv file
        self.source_cat = self.annotate_table(self.source_cat, self.param_dict_qc, product=self.image.ghd_product)
        # self.source_cat.meta['comments'] = \
        #     ["NOTE: The X and Y coordinates in this table are 0-indexed (i.e. the origin is (0,0))."]
        self.source_cat.write(self.sourcelist_filename, format=self.catalog_format)
        log.info("Wrote catalog file '{}' containing {} sources".format(self.sourcelist_filename, len(self.source_cat)))

        # Write out region file if input 'write_region_file' is turned on.
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
        if 'RA' in self.sources.colnames and 'DEC' in self.sources.colnames:
            subset_table.remove_columns(['RA', 'DEC'])
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

    def __init__(self, image, param_dict, param_dict_qc, diagnostic_mode, tp_sources):
        super().__init__(image, param_dict, param_dict_qc, diagnostic_mode, tp_sources)

        # Get the instrument/detector-specific values from the self.param_dict
        self._fwhm = self.param_dict["sourcex"]["fwhm"]
        self._size_source_box = self.param_dict["sourcex"]["source_box"]
        self._nlevels = self.param_dict["sourcex"]["nlevels"]
        self._contrast = self.param_dict["sourcex"]["contrast"]
        self._border = self.param_dict["sourcex"]["border"]
        self._nsigma = self.param_dict["nsigma"]

        # Initialize attributes to be computed later
        self.segm_img = None  # Segmentation image

        # Defined in measure_sources
        self.subset_filter_source_cat = None

        # FIX
        self.kernel = self.image.kernel

    def identify_sources(self, **pars):
        """Use photutils to find sources in image based on segmentation.

        Returns
        -------
        self.sources
        self.source_catalog

        Defines
        -------
        self.segm_img : `photutils.segmentation.SegmentationImage`
            Two-dimensional segmentation image where found source regions are labeled with
            unique, non-zero positive integers.
        """

        # Check the pars to see if the "n" exposure mask has been passed.  If it has, it must be
        # converted into a boolean mask where False = actual image footprint, and True = the
        # non-illuminated portion of the images.  The True indicates to detect_sources()
        # which pixels to ignore.
        mask = pars.get('mask', None)
        if hasattr(mask, 'shape'):
            mask = mask < 1
            # outname = self.imgname.replace(".fits","_mask.fits")
            # fits.PrimaryHDU(data=mask.astype(np.uint16)).writeto(outname)

        # If the total product sources have not been identified, then this needs to be done!
        if not self.tp_sources:

            # Report configuration values to log
            log.info("{}".format("=" * 80))
            log.info("")
            log.info("SExtractor-like source finding settings - Photutils segmentation")
            log.info("Total Detection Product - Input Parameters")
            log.info("FWHM: {}".format(self._fwhm))
            log.info("size_source_box (no. of connected pixels needed for a detection): {}".format(self._size_source_box))
            log.info("nsigma (sigma * background_rms): {}".format(self._nsigma))
            log.info("nlevels (no. of multi-thresholding levels for deblending): {}".format(self._nlevels))
            log.info("contrast (frac. flux for peak to be separate object, 0=max. deblend, 1=no deblend): {}".format(self._contrast))
            log.info("border (image border width where sources not detected): {}".format(self._border))
            log.info("")
            log.info("{}".format("=" * 80))

            # Get the SCI image data
            imgarr = self.image.data.copy()

            # The bkg is an object comprised of background and background_rms images, as well as
            # background_median and background_rms_median scalars.  Set a threshold above which
            # sources can be detected.
            threshold = self._nsigma * self.image.bkg_rms_ra

            # The imgarr should be background subtracted to match the threshold which has no background
            imgarr_bkgsub = imgarr - self.image.bkg_background_ra

            log.info("Detecting sources in total image product.")
            # Note: SExtractor has "connectivity=8" which is the default for detect_sources().
            self.segm_img = detect_sources(imgarr_bkgsub, threshold, npixels=self._size_source_box,
                                           filter_kernel=self.image.kernel,
                                           mask=mask)

            try:
                # Deblending is a combination of multi-thresholding and watershed
                # segmentation. Sextractor uses a multi-thresholding technique.
                # npixels = number of connected pixels in source
                # npixels and filter_kernel should match those used by detect_sources()
                segm_deblended_img = deblend_sources(imgarr_bkgsub, self.segm_img, npixels=self._size_source_box,
                                                     filter_kernel=self.image.kernel, nlevels=self._nlevels,
                                                     contrast=self._contrast)

                # The deblending was successful, so just copy the deblended sources back to the sources attribute.
                self.segm_img = copy.deepcopy(segm_deblended_img)
            except Exception as x_cept:
                log.warning("Deblending the sources in image {} was not successful: {}.".format(self.imgname,
                                                                                                         x_cept))
                log.warning("Processing can continue with the non-deblended sources, but the user should\n"
                            "check the output catalog for issues.")

            # Regardless of whether or not deblending worked, this variable can be reset to None
            segm_deblended_img = None

            # Clean up segments that are near or partially overlap the border of the image and relabel
            # the segments to be sequential
            self.segm_img.remove_border_labels(self._border, partial_overlap=True, relabel=True)

            # The total product catalog consists of at least the X/Y and RA/Dec coordinates for the detected
            # sources in the total drizzled image.  All the actual measurements are done on the filtered drizzled
            # images using the coordinates determined from the total drizzled image.
            self.source_cat = source_properties(imgarr_bkgsub, self.segm_img, background=self.image.bkg_background_ra,
                                                filter_kernel=self.image.kernel, wcs=self.image.imgwcs)

            # Convert source_cat which is a SourceCatalog to an Astropy Table - need the data in tabular
            # form to filter out bad rows and correspondingly bad segments before the filter images are processed.
            total_measurements_table = Table(self.source_cat.to_table())

            # Filter the table to eliminate nans or inf based on the coordinates, then remove these segments from
            # the segmentation image
            good_rows = []
            bad_segm_rows_by_id = []
            updated_table = None
            for i, old_row in enumerate(total_measurements_table):
                if np.isfinite(old_row["xcentroid"]):
                    good_rows.append(old_row)
                else:
                    bad_segm_rows_by_id.append(total_measurements_table['id'][i])
            updated_table = Table(rows=good_rows, names=total_measurements_table.colnames)
            if self.diagnostic_mode:
                log.info("Bad total rows: {}".format(bad_segm_rows_by_id))
            log.info("Bad segments removed from segmentation image.")

            # Remove the bad segments from the image
            self.segm_img.remove_labels(bad_segm_rows_by_id, relabel=True)

            # Clean up the existing column names, format, and descriptions
            self.source_cat = self._define_total_table(updated_table)

            # self.sources needs to be passed to a filter catalog object based on code in hapsequencer.py
            # (create_catalog_products()).  This is the way the independent catalogs of total and filter products
            # process the same segmentation image.
            # BEWARE: self.sources for "segmentation" is a SegmentationImage, but for "point" it is an Astropy table
            self.sources = copy.deepcopy(self.segm_img)

        # If filter product, use sources identified in total detection product previously generated
        else:
            self.sources = self.tp_sources['segment']['sources']
            self.kernel = self.tp_sources['segment']['kernel']

        # For debugging purposes only, create a segmentation image and a "regions" files to use for ds9 overlay
        # Create the image regions file here in case there is a failure
        if self.diagnostic_mode and self.segm_img:
            # Generate a diagnostic_mode segmentation image
            # indx = self.sourcelist_filename.find("-cat.ecsv")
            # outname = self.sourcelist_filename[0:indx] + ".fits"
            # fits.PrimaryHDU(data = self.segm_img.data).writeto(outname)

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
        # Get filter-level science data
        imgarr = self.image.data.copy()

        # Report configuration values to log
        log.info("{}".format("=" * 80))
        log.info("")
        log.info("SExtractor-like source property measurements based on Photutils segmentation")
        log.info("Filter Level Product - Input Parameters")
        log.info("FWHM: {}".format(self._fwhm))
        log.info("size_source_box: {}".format(self._size_source_box))
        log.info("")
        log.info("{}".format("=" * 80))

        # This is the filter science data and its computed background
        imgarr_bkgsub = imgarr - self.image.bkg_background_ra

        # Compute source properties...
        self.source_cat = source_properties(imgarr_bkgsub, self.sources, background=self.image.bkg_background_ra,
                                            filter_kernel=self.image.kernel, wcs=self.image.imgwcs)

        # Convert source_cat which is a SourceCatalog to an Astropy Table
        filter_measurements_table = Table(self.source_cat.to_table())

        # Compute the MagIso
        filter_measurements_table["MagIso"] = photometry_tools.convert_flux_to_abmag(filter_measurements_table["source_sum"],
                                                                                     self.image.imghdu[1].header['photflam'],
                                                                                     self.image.imghdu[1].header['photplam'])

        # Compute aperture photometry measurements and append the columns to the measurements table
        updated_table = self.do_aperture_photometry(imgarr_bkgsub, filter_measurements_table)

        # Now clean up and prepare the filter table for output
        self.source_cat = self._define_filter_table(updated_table)

        log.info("Found and measured {} sources from segmentation map.".format(len(self.source_cat)))

        # Capture specified filter columns in order to append to the total detection table
        self.subset_filter_source_cat = self.source_cat["ID", "MagAp2", "CI", "Flags"]
        self.subset_filter_source_cat.rename_column("MagAp2", "MagAp2_" + filter_name)
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
            outname = self.sourcelist_filename[0:indx] + "reg"

            tbl["X-Centroid"].info.format = ".10f"  # optional format
            tbl["Y-Centroid"].info.format = ".10f"

            # Add one to the X and Y table values to put the data onto a one-based system,
            # particularly for display with ds9
            tbl["X-Centroid"] = tbl["X-Centroid"] + 1
            tbl["Y-Centroid"] = tbl["Y-Centroid"] + 1
            tbl.write(outname, format="ascii.commented_header")
            log.info("Wrote the diagnostic_mode version of the filter detection source catalog: {}\n".format(outname))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def do_aperture_photometry(self, bkg_subtracted_image, filter_measurements_table):
        """Perform aperture photometry measurements as a means to distinguish point versus extended sources.
        """
        # Filter the table to eliminate nans or inf based on the coordinates now that
        # measurements have been done on the filter image
        good_rows = []
        bad_rows = []
        updated_table = None
        for i, old_row in enumerate(filter_measurements_table):
            if np.isfinite(old_row["xcentroid"]):
                good_rows.append(old_row)
            else:
                bad_rows.append(filter_measurements_table['id'][i])
        updated_table = Table(rows=good_rows, names=filter_measurements_table.colnames)
        # FIX What filter?
        log.info("Bad rows removed from coordinate list for filter data based on invalid positions after source property measurements.")

        positions = (updated_table["xcentroid"], updated_table["ycentroid"])
        pos_xy = np.vstack(positions).T

        # Define list of background annulii
        bg_apers = CircularAnnulus(pos_xy,
                                   r_in=self.param_dict['skyannulus_arcsec'],
                                   r_out=self.param_dict['skyannulus_arcsec'] +
                                   self.param_dict['dskyannulus_arcsec'])

        # Create list of photometric apertures to measure
        phot_apers = [CircularAperture(pos_xy, r=r) for r in self.aper_radius_list_pixels]

        # Perform aperture photometry
        photometry_tbl = photometry_tools.iraf_style_photometry(phot_apers,
                                                                bg_apers,
                                                                data=bkg_subtracted_image,
                                                                photflam=self.image.imghdu[1].header['photflam'],
                                                                photplam=self.image.imghdu[1].header['photplam'],
                                                                error_array=self.bkg.background_rms,
                                                                bg_method=self.param_dict['salgorithm'],
                                                                epadu=self.gain)

        # Capture data computed by the photometry tools and append to the output table
        try:
            flux_inner_data = photometry_tbl["FluxAp1"].data
            flux_inner_data_err = photometry_tbl["FluxErrAp1"].data
            mag_inner_data = photometry_tbl["MagAp1"].data
            mag_inner_data_err = photometry_tbl["MagErrAp1"].data

            flux_outer_data = photometry_tbl["FluxAp2"].data
            flux_outer_data_err = photometry_tbl["FluxErrAp2"].data
            mag_outer_data = photometry_tbl["MagAp2"].data
            mag_outer_data_err = photometry_tbl["MagErrAp2"].data

            mskyap2_data = photometry_tbl["MSkyAp2"].data
            stdevap2_data = photometry_tbl["StdevAp2"].data

            ci_data = mag_inner_data - mag_outer_data
            ci_mask = np.logical_and(np.abs(ci_data) > 0.0, np.abs(ci_data) < 1.0e-30)
            big_bad_index = np.where(abs(ci_data) > 1.0e20)
            ci_mask[big_bad_index] = True
            ci_col = MaskedColumn(name="CI", data=ci_data, dtype=np.float64, mask=ci_mask)
            updated_table.add_column(ci_col)

            # Append these additional photometric measurements to the filter table
            flux_col = Column(data=flux_inner_data, name="FluxAp1", dtype=np.float64)
            flux_col_err = Column(data=flux_inner_data_err, name="FluxErrAp1", dtype=np.float64)
            mag_col = Column(data=mag_inner_data, name="MagAp1", dtype=np.float64)
            mag_col_err = Column(data=mag_inner_data_err, name="MagErrAp1", dtype=np.float64)
            updated_table.add_column(flux_col)
            updated_table.add_column(flux_col_err)
            updated_table.add_column(mag_col)
            updated_table.add_column(mag_col_err)

            flux_col = Column(data=flux_outer_data, name="FluxAp2", dtype=np.float64)
            flux_col_err = Column(data=flux_outer_data_err, name="FluxErrAp2", dtype=np.float64)
            mag_col = Column(data=mag_outer_data, name="MagAp2", dtype=np.float64)
            mag_col_err = Column(data=mag_outer_data_err, name="MagErrAp2", dtype=np.float64)
            updated_table.add_column(flux_col)
            updated_table.add_column(flux_col_err)
            updated_table.add_column(mag_col)
            updated_table.add_column(mag_col_err)

            msky_col = Column(data=mskyap2_data, name="MSkyAp2", dtype=np.float64)
            stdev_col = Column(data=stdevap2_data, name="StdevAp2", dtype=np.float64)
            updated_table.add_column(msky_col)
            updated_table.add_column(stdev_col)

        except Exception as x_cept:
            log.warning("Computation of additional photometric measurements was not successful: {} - {}.".format(self.imgname, x_cept))
            log.warning("Additional measurements have not been added to the output catalog.\n")

        # Add zero-value "Flags" column in preparation for source flagging
        flag_col = Column(name="Flags", data=np.zeros_like(updated_table["id"]))
        updated_table.add_column(flag_col)

        # Protect against None for source_sum_error
        # FIX MDD
        # log.info("*************** source sum err: {}".format(updated_table["source_sum_err"]))
        # None_index = np.where(source_sum_err is None)
        # sse_mask[None_index] = True
        # ci_col = MaskedColumn(name="CI", data=source_sum_err, dtype=np.float64, mask=sse_mask)

        return updated_table

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
           ------
           final_filter_table : Astropy table
                A modified version of the input table which has been reformatted in preparation
                for catalog generation.
        """

        radec_data = SkyCoord(filter_table["sky_centroid_icrs"])
        ra_icrs = radec_data.ra.degree
        dec_icrs = radec_data.dec.degree
        rr = Column(ra_icrs, name="RA", unit=u.deg)
        dd = Column(dec_icrs, name="DEC", unit=u.deg)
        filter_table.add_columns([dd, rr])

        # Rename columns to names used when HLA Classic catalog distributed by MAST
        final_col_names = {"id": "ID", "xcentroid": "X-Centroid", "ycentroid": "Y-Centroid",
                           "background_at_centroid": "Bck", "source_sum": "FluxIso",
                           # "background_at_centroid": "Bck", "source_sum": "FluxIso", "source_sum_err": "FluxIsoErr",
                           "bbox_xmin": "Xmin", "bbox_ymin": "Ymin", "bbox_xmax": "Xmax", "bbox_ymax": "Ymax",
                           "cxx": "CXX", "cyy": "CYY", "cxy": "CXY",
                           "covar_sigx2": "X2", "covar_sigy2": "Y2", "covar_sigxy": "XY",
                           "orientation": "Theta",
                           "elongation": "Elongation", "ellipticity": "Ellipticity"}
        for old_col_title in final_col_names:
            filter_table.rename_column(old_col_title, final_col_names[old_col_title])

        # Define the order of the columns
        final_col_order = ["X-Centroid", "Y-Centroid", "RA", "DEC", "ID",
                           "CI", "Flags", "MagAp1", "MagErrAp1", "FluxAp1", "FluxErrAp1",
                           "MagAp2", "MagErrAp2", "FluxAp2", "FluxErrAp2", "MSkyAp2",
                           "Bck", "MagIso", "FluxIso",
                           # "Bck", "MagIso", "FluxIso", "FluxIsoErr",
                           "Xmin", "Ymin", "Xmax", "Ymax",
                           "X2", "Y2", "XY",
                           "CXX", "CYY", "CXY",
                           "Elongation", "Ellipticity", "Theta"]
        final_filter_table = filter_table[final_col_order]

        # Define the format
        final_col_format = {"X-Centroid": "10.3f", "Y-Centroid": "10.3f", "RA": "13.7f", "DEC": "13.7f", "ID": "6d",
                            "CI": "7.3f", "Flags": "5d",
                            "MagAp1": "8.2f", "MagErrAp1": "9.4f", "FluxAp1": "9.2f", "FluxErrAp1": "10.5f",
                            "MagAp2": "8.2f", "MagErrAp2": "9.4f", "FluxAp2": "9.2f", "FluxErrAp2": "10.5f",
                            "MSkyAp2": "8.2f", "Bck": "9.4f", "MagIso": "8.2f", "FluxIso": "9.2f",
                            # "Bck": "9.4f", "MagIso": "8.2f", "FluxIso": "9.2f", "FluxIsoErr": "10.5f",
                            "Xmin": "8.0f", "Ymin": "8.0f", "Xmax": "8.0f", "Ymax": "8.0f",
                            "X2": "8.4f", "Y2": "8.4f", "XY": "10.5f",
                            "CXX": "9.5f", "CYY": "9.5f", "CXY": "9.5f",
                            "Elongation": "7.2f", "Ellipticity": "7.2f", "Theta": "8.3f"}
        for fcf_key in final_col_format.keys():
            final_filter_table[fcf_key].format = final_col_format[fcf_key]

        # Add description
        final_col_descrip = {"ID": "Catalog Object Identification Number",
                             "X-Centroid": "Pixel Coordinate", "Y-Centroid": "Pixel Coordinate",
                             "RA": "Sky coordinate at epoch of observation and fit to GAIA",
                             "DEC": "Sky coordinate at epoch of observation and fit to GAIA",
                             "Bck": "Background at the position of the source centroid",
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
                             # "FluxIsoErr": "Uncertainty of FluxIso propagated from the input error array",
                             "MagIso": "Magnitude corresponding to FluxIso",
                             # "MagIsoErr: "Uncertainty of MagIso corresponding to FluxIsoErr"
                             "X2": "Variance along X",
                             "Y2": "Variance along Y",
                             "XY": "Covariance of position between X and Y",
                             "CXX": "SExtractor's ellipse parameter", "CYY": "SExtractor's ellipse parameter",
                             "CXY": "SExtractor's ellipse parameter",
                             "Xmin": "Minimum X pixel within the minimal bounding box containing the source segment",
                             "Xmax": "Maximum X pixel within the minimal bounding box containing the source segment",
                             "Ymin": "Minimum Y pixel within the minimal bounding box containing the source segment",
                             "Ymax": "Maximum Y pixel within the minimal bounding box containing the source segment",
                             "Elongation": "Ratio of the lengths of the semimajor and semiminor axes of the ellipse",
                             "Ellipticity": "The value '1 minus the elongation",
                             "Theta": "Angle between the X axis and the major axis of the 2D Gaussian function that has the same second-order moments as the source.",

                             "CI": "Concentration Index"}
        for fcd_key in final_col_descrip.keys():
            final_filter_table[fcd_key].description = final_col_descrip[fcd_key]

        # Add units
        final_col_unit = {"X-Centroid": u.pix, "Y-Centroid": u.pix,
                          "RA": u.deg, "DEC": u.deg,
                          "Bck": "electrons/s",
                          "MagAp1": "ABMAG",
                          "MagErrAp1": "ABMAG",
                          "FluxAp1": "electrons/s",
                          "FluxErrAp1": "electrons/s",
                          "MagAp2": "ABMAG",
                          "MagErrAp2": "ABMAG",
                          "FluxAp2": "electrons/s",
                          "FluxErrAp2": "electrons/s",
                          "MagIso": "ABMAG",
                          "FluxIso": "electrons/s",
                          "X2": "pixel**2",
                          "Y2": "pixel**2",
                          "XY": "pixel**2",
                          "CXX": "pixel**2",
                          "CYY": "pixel**2",
                          "CXY": "pixel**2",
                          "Xmin": u.pix,
                          "Ymin": u.pix,
                          "Xmax": u.pix,
                          "Ymax": u.pix,
                          "Theta": u.rad}
        for fcu_key in final_col_unit.keys():
            final_filter_table[fcu_key].unit = final_col_unit[fcu_key]

        return(final_filter_table)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

        # Extract just a few columns generated by the source_properties() as
        # more columns are appended to this table from the filter results.
        # Actually, the filter columns are in a table which is "database joined"
        # to the total table.  During the combine process, the new columns are renamed,
        # formatted, and described (as necessary). For now this table only has id, xcentroid,
        # ycentroid, RA, and DEC.
        table = updated_table["id", "xcentroid", "ycentroid"]

        # Convert the RA/Dec SkyCoord into separate columns
        radec_data = SkyCoord(updated_table["sky_centroid_icrs"])
        ra_icrs = radec_data.ra.degree
        dec_icrs = radec_data.dec.degree
        rr = Column(ra_icrs, name="RA", unit=u.deg)
        dd = Column(dec_icrs, name="DEC", unit=u.deg)
        table.add_columns([rr, dd])

        # Rename columns to names to those used when HLA Classic catalog distributed by MAST
        # and/or to distinguish Point and Segment catalogs
        # The columns that are appended will be renamed during the combine process
        final_col_names = {"id": "ID", "xcentroid": "X-Centroid", "ycentroid": "Y-Centroid"}
        for old_col_title in final_col_names:
            table.rename_column(old_col_title, final_col_names[old_col_title])

        # Format the current columns
        final_col_format = {"ID": "6d", "X-Centroid": "10.3f", "Y-Centroid": "10.3f", "RA": "13.7f", "DEC": "13.7f"}
        for fcf_key in final_col_format.keys():
            table[fcf_key].format = final_col_format[fcf_key]

        # Add description
        descr_str = "Sky coordinate at epoch of observation and " + self.image.keyword_dict["wcs_type"]
        final_col_descrip = {"ID": "Catalog Object Identification Number",
                             "X-Centroid": "Pixel Coordinate", "Y-Centroid": "Pixel Coordinate",
                             "RA": descr_str, "DEC": descr_str}
        for fcd_key in final_col_descrip.keys():
            table[fcd_key].description = final_col_descrip[fcd_key]

        # Add units
        final_col_unit = {"X-Centroid": u.pix, "Y-Centroid": u.pix,
                          "RA": u.deg, "DEC": u.deg}
        for fcu_key in final_col_unit.keys():
            table[fcu_key].unit = final_col_unit[fcu_key]

        return(table)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    @property
    def write_catalog(self):
        """Write the specified source catalog out to disk.
        """
        self.source_cat = self.annotate_table(self.source_cat, self.param_dict_qc, product=self.image.ghd_product)
        self.source_cat.write(self.sourcelist_filename, format=self.catalog_format)
        log.info("Wrote filter source catalog: {}".format(self.sourcelist_filename))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def combine_tables(self, subset_table):
        """Append specified measurements from the filter table to the total detection table.

        The "ID" column is used to map the filter table measurements to the total detection table

        Parameters
        ----------
        subset_table : Astropy table
            A table containing a subset of columns from a filter catalog.
        """

        # Keep all the rows in the original total detection table and add columns from the filter
        # table where a matching "id" key is present
        self.source_cat = join(self.source_cat, subset_table, keys="ID", join_type="left")
