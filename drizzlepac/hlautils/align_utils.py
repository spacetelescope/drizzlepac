import os
import datetime
import copy
import sys
import traceback

from collections import OrderedDict

import numpy as np
from scipy import ndimage

from astropy.io import fits
from astropy.table import Table
from astropy.nddata.bitmask import bitfield_to_boolean_mask
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u

import photutils
from photutils import Background2D

from stwcs.wcsutil import HSTWCS
from stwcs.wcsutil import headerlet

from stsci.tools import logutil
from stsci.tools import fileutil

from .. import updatehdr
from . import astrometric_utils as amutils
from . import analyze

import tweakwcs

__taskname__ = 'align_utils'

# Default background determination parameter values
BKG_BOX_SIZE = 27
BKG_FILTER_SIZE = 3
CATALOG_TYPES = ['point', 'segment']
MIN_CATALOG_THRESHOLD = 3

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)

class AlignmentTable:
    def __init__(self, input_list, clobber=False, dqname='DQ',
                 log_level=logutil.logging.NOTSET, **alignment_pars):
        """
        **alignment_pars needs to contain the following entries:
                          # kernel defining, source finding par
                          fwhmpsf=0.12,
                          # background computing pars
                          box_size=BKG_BOX_SIZE, win_size=BKG_FILTER_SIZE,
                          bkg_estimator=SExtractorBackground,
                          rms_estimator=StdBackgroundRMS,
                          nsigma=5., threshold_flag=None,
                          # object finding pars
                          source_box=7,
                          classify=True, centering_mode="starfind", nlargest=None,
                          plot=False, vmax=None, deblend=False
        """
        log.setLevel(log_level)
        # Register fit methods with the class
        self.fit_methods = {'relative': match_relative_fit,
                            '2dhist': match_2dhist_fit,
                            'default': match_default_fit}

        # Also, register fit method default parameters with the class
        self.fit_pars = {'relative': alignment_pars['match_relative_fit'],
                         '2dhist': alignment_pars['match_2dhist_fit'],
                         'default': alignment_pars['match_default_fit']}

        # merge remaining parameters for individual alignment steps into single set
        self.alignment_pars = alignment_pars['run_align'].copy()
        self.alignment_pars.update(alignment_pars['general'])
        self.alignment_pars.update(alignment_pars['generate_source_catalogs'])
        self.alignment_pars.update(alignment_pars['determine_fit_quality'])

        self.dqname = dqname

        self.zero_dt = starting_dt = datetime.datetime.now()
        log.info(str(starting_dt))
        # Apply filter to input observations to insure that they meet minimum criteria for being able to be aligned
        log.info(
            "{} AlignmentTable: Filter STEP {}".format("-" * 20, "-" * 63))
        self.filtered_table = analyze.analyze_data(input_list)

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

        fwhmpsf = self.alignment_pars.get('fwhmpsf')
        default_fwhm_set = False

        try:
            self.haplist = []
            for img in self.process_list:
                catimg = HAPImage(img)
                # Build image properties needed for alignment
                catimg.compute_background(box_size=self.alignment_pars['box_size'],
                                          win_size=self.alignment_pars['win_size'],
                                          nsigma=self.alignment_pars['nsigma'],
                                          bkg_estimator=self.alignment_pars['bkg_estimator'],
                                          rms_estimator=self.alignment_pars['rms_estimator'],
                                          threshold_flag=self.alignment_pars['threshold'])
                catimg.build_kernel(fwhmpsf)
                # Use FWHM from first good exposure as default for remainder of exposures
                if not default_fwhm_set and catimg.kernel is not None:
                    fwhmpsf = catimg.fwhmpsf
                    default_fwhm_set = True

                self.haplist.append(catimg)
        except Exception:
            log.error("Problem encountered in setting up alignment to a catalog")
            traceback.print_exc()

            self.close()

        # Initialize computed attributes
        self.imglist = []  # list of FITSWCS objects for tweakwcs
        self.reference_catalogs = {}
        self.group_id_dict = {}

        self.fit_dict = {}  # store results of all fitting in this dict for evaluation
        self.selected_fit = None  # identify fit selected by user (as 'best'?) to be applied to images WCSs

    def close(self):
        for img in self.haplist:
            img.close()

    def find_alignment_sources(self, output=True):
        """Find observable sources in each input exposure."""
        self.extracted_sources = {}
        for img in self.haplist:
            img.find_alignment_sources(output=output, dqname=self.dqname, **self.alignment_pars)
            self.extracted_sources[img.imgname] = img.catalog_table

            # Allow user to decide when and how to write out catalogs to files
            if output:
                # write out coord lists to files for diagnostic purposes. Protip: To display the sources in these files in DS9,
                # set the "Coordinate System" option to "Physical" when loading the region file.
                imgroot = os.path.basename(img.imgname).split('_')[0]
                for chip in range(1, img.num_sci + 1):
                    chip_cat = img.catalog_table[chip]
                    if chip_cat and len(chip_cat) > 0:
                        regfilename = "{}_sci{}_src.reg".format(imgroot, chip)
                        out_table = Table(chip_cat)
                        # To align with positions of sources in DS9/IRAF
                        out_table['xcentroid'] += 1
                        out_table['ycentroid'] += 1
                        out_table.write(regfilename,
                                        include_names=["xcentroid", "ycentroid"],
                                        format="ascii.fast_commented_header")
                        log.info("Wrote region file {}\n".format(regfilename))

    def reset_group_id(self, num_ref):
        for image in self.imglist:
            image.meta["group_id"] = self.group_id_dict["{}_{}".format(image.meta["filename"], image.meta["chip"])]
            image.meta['num_ref_catalog'] = num_ref

    def configure_fit(self):
        # Convert input images to tweakwcs-compatible FITSWCS objects and
        # attach source catalogs to them.

        self.imglist = []
        for group_id, image in enumerate(self.process_list):
            img = amutils.build_wcscat(image, group_id,
                                       self.extracted_sources[image])
            # add the name of the image to the imglist object
            for im in img:
            #    im.meta['name'] = image
                log.debug('im.meta[name] = {}'.format(im.meta['name']))
            self.imglist.extend(img)

        self.group_id_dict = {}
        for image in self.imglist:
            self.group_id_dict["{}_{}".format(image.meta["filename"], image.meta["chip"])] = image.meta["group_id"]

    def get_fit_methods(self):
        """Return the list of method names for all registered functions
            available for performing alignment.
        """
        return list(self.fit_methods.keys())

    def perform_fit(self, method_name, catalog_name, reference_catalog):
        """Perform fit using specified method, then determine fit quality"""
        imglist = self.fit_methods[method_name](self.imglist, reference_catalog,
                                                **self.fit_pars[method_name])

        # save fit algorithm name to dictionary key "fit method" in imglist.
        for imglist_ctr in range(0, len(imglist)):
            imglist[imglist_ctr].meta['fit method'] = method_name

        # store results for evaluation
        self.fit_dict[(catalog_name, method_name)] = copy.deepcopy(imglist)
        self.reference_catalogs[catalog_name] = reference_catalog

        return imglist

    def select_fit(self, catalog_name, method_name):
        """Select the fit that has been identified as 'best'"""
        imglist = self.selected_fit = self.fit_dict[(catalog_name, method_name)]
        if imglist[0].meta['fit_info']['status'].startswith("FAILED"):
            self.selected_fit = None

        # Protect the writing of the table within the best_fit_rms
        info_keys = OrderedDict(imglist[0].meta['fit_info']).keys()
        # Update filtered table with number of matched sources and other information
        for item in imglist:
            imgname = item.meta['name']
            index = np.where(self.filtered_table['imageName'] == imgname)[0][0]

            if not item.meta['fit_info']['status'].startswith("FAILED"):
                for tweakwcs_info_key in info_keys:
                    if not tweakwcs_info_key.startswith("matched"):
                        if tweakwcs_info_key.lower() == 'rms':
                            self.filtered_table[index]['rms_x'] = item.meta['fit_info'][tweakwcs_info_key][0]
                            self.filtered_table[index]['rms_y'] = item.meta['fit_info'][tweakwcs_info_key][1]

                self.filtered_table[index]['fit_method'] = item.meta['fit method']
                self.filtered_table[index]['catalog'] = item.meta['fit_info']['catalog']
                self.filtered_table[index]['catalogSources'] = len(self.reference_catalogs[catalog_name])
                self.filtered_table[index]['matchSources'] = item.meta['fit_info']['nmatches']
                self.filtered_table[index]['rms_ra'] = item.meta['fit_info']['RMS_RA'].value
                self.filtered_table[index]['rms_dec'] = item.meta['fit_info']['RMS_DEC'].value
                self.filtered_table[index]['fit_rms'] = item.meta['fit_info']['FIT_RMS']
                self.filtered_table[index]['total_rms'] = item.meta['fit_info']['TOTAL_RMS']
                self.filtered_table[index]['offset_x'], self.filtered_table[index]['offset_y'] = item.meta['fit_info']['shift']
                self.filtered_table[index]['scale'] = item.meta['fit_info']['scale'][0]
                self.filtered_table[index]['rotation'] = item.meta['fit_info']['rot']
            else:
                self.filtered_table[index]['fit_method'] = None


    def apply_fit(self, headerlet_filenames=None, fit_label=None):
        """Apply solution from identified fit to image WCS's

        Parameters
        ----------
        headerlet_filenames : dict, optional
            Dictionary relating exposure filenames to headerlet filenames.  If None,
            will generate headerlet filenames where _flt or _flc is replaced by
            _flt_hlet or _flc_hlet, respectively.

        fit_label : str
            Name of fit to apply to indicate how the fit was performed in
            the WCSNAME keyword.  Common options: IMG, REL, SVM.

        """
        if not self.selected_fit:
            log.error("No FIT selected for application.  Please run 'select_fit()' method.")
            raise ValueError
        # Call update_hdr_wcs()
        headerlet_dict = update_image_wcs_info(self.selected_fit,
                                               headerlet_filenames=headerlet_filenames,
                                               fit_label=fit_label)

        for table_index in range(0, len(self.filtered_table)):
            self.filtered_table[table_index]['headerletFile'] = headerlet_dict[
                self.filtered_table[table_index]['imageName']]


class HAPImage:
    """Core class defining interface for each input exposure/product

    .. note:: This class is compatible with the CatalogImage class, while including
    additional functionality and attributes required for processing beyond
    catalog generation.

    """

    def __init__(self, filename):
        if isinstance(filename, str):
            self.imghdu = fits.open(filename)
            self.imgname = filename
        else:
            self.imghdu = filename
            self.imgname = filename.filename()

        if 'rootname' in self.imghdu[0].header:
            self.rootname = self.imghdu[0].header['rootname']
        else:
            self.rootname = self.imgname.rstrip('.fits')

        # Fits file read
        self.num_sci = amutils.countExtn(self.imghdu)
        self.num_wht = amutils.countExtn(self.imghdu, extname='WHT')
        self.data = np.concatenate([self.imghdu[('SCI', i + 1)].data for i in range(self.num_sci)])
        if not self.num_wht:
            self.dqmask = self.build_dqmask()
        else:
            self.dqmask = None
        self.wht_image = self.build_wht_image()

        # Get the HSTWCS object from the first extension
        self.imgwcs = HSTWCS(self.imghdu, 1)
        self.pscale = self.imgwcs.pscale

        self._wht_image = None
        self.bkg = {}
        self.bkg_dao_rms = {}
        self.bkg_rms_mean = {}
        self.threshold = {}

        self.kernel = None
        self.kernel_fwhm = None
        self.kernel_psf = False
        self.fwhmpsf = None

        self.catalog_table = {}

    def build_wht_image(self):
        if not self.num_wht:
            # Working with a calibrated exposure, no WHT extension
            # create a substitute WHT array from ERR and DQ
            # Build pseudo-wht array for detection purposes
            errarr = np.concatenate([self.imghdu[('ERR', i + 1)].data for i in range(self.num_sci)])
            wht_image = errarr.max() / errarr
            if self.dqmask is not None:
                wht_image[self.dqmask] = 0
        else:
            wht_image = self.imghdu['WHT'].data
        return wht_image

    def close(self):
        self.imghdu.close()

    def build_kernel(self, fwhmpsf):
        """
        Parameters
        -----------
        fwhmpsf : float
            Default FWHM of PSF in units of arcseconds.
        """
        if self.bkg is None:
            self.compute_background()

        threshold_rms = np.concatenate([rms for rms in self.threshold.values()])
        bkg = np.concatenate([background for background in self.bkg.values()])
        log.info("Looking for sample PSF in {}".format(self.rootname))
        log.debug("  based on RMS of {}".format(threshold_rms.mean()))
        fwhm = fwhmpsf / self.pscale
        
        k, self.kernel_fwhm = amutils.build_auto_kernel(self.data - bkg,
                                                        self.wht_image,
                                                        threshold=threshold_rms,
                                                        fwhm=fwhm)
        if not k[1]:
            # Try one more time with slightly different background...
            self.compute_background(box_size=BKG_BOX_SIZE // 2)
            threshold_rms = np.concatenate([rms for rms in self.threshold.values()])
            bkg = np.concatenate([background for background in self.bkg.values()])
            log.info("Looking for sample PSF in {}".format(self.rootname))
            log.debug("  based on RMS of {}".format(threshold_rms.mean()))
            fwhm = fwhmpsf / self.pscale
            
            k, self.kernel_fwhm = amutils.build_auto_kernel(self.data - bkg,
                                                            self.wht_image,
                                                            threshold=threshold_rms,
                                                            fwhm=fwhm)

        self.kernel, self.kernel_psf = k
        log.info("  Found PSF with FWHM = {}".format(self.kernel_fwhm))

        self.fwhmpsf = self.kernel_fwhm * self.pscale

    def compute_background(self, box_size=BKG_BOX_SIZE, win_size=BKG_FILTER_SIZE,
                           bkg_estimator="SExtractorBackground", rms_estimator="StdBackgroundRMS",
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
        # Interpret estimator labels as photutils functions
        bkg_estimator = register_photutils_function(bkg_estimator)
        rms_estimator = register_photutils_function(rms_estimator)

        # Report configuration values to log
        log.debug("")
        log.debug("Computation of {} background - Input Parameters".format(self.rootname))
        log.debug("Box size: {}".format(box_size))
        log.debug("Window size: {}".format(win_size))
        log.debug("NSigma: {}".format(nsigma))
        log.debug("BKG Estimator: {}".format(bkg_estimator.__name__))

        # SExtractorBackground ans StdBackgroundRMS are the defaults
        exclude_percentiles = [10, 25, 50, 75]

        for chip in range(self.num_sci):
            chip += 1
            bkg = None
            bkg_dao_rms = None
            scidata = self.imghdu[('sci', chip)].data

            for percentile in exclude_percentiles:
                try:
                    bkg = Background2D(scidata, box_size, filter_size=win_size,
                                       bkg_estimator=bkg_estimator(),
                                       bkgrms_estimator=rms_estimator(),
                                       exclude_percentile=percentile, edge_method="pad")
                except Exception:
                    bkg = None
                    continue

                if bkg is not None:
                    bkg_mean = bkg.background_median
                    bkg_dao_rms = bkg.background_rms
                    # Set the bkg_rms at "nsigma" sigma above background
                    default_threshold = bkg.background + nsigma * bkg.background_rms
                    bkg_rms_mean = bkg.background_rms_median if bkg.background_rms_median > 0 else 0.

                    if threshold_flag is None:
                        threshold = default_threshold
                    elif threshold_flag < 0:
                        threshold = -1 * threshold_flag * default_threshold
                        bkg_rms_mean = -1 * threshold_flag * bkg_rms_mean
                    else:
                        bkg_rms_mean = 3. * threshold_flag
                        threshold = default_threshold

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
            log.debug("{} CHIP: {}".format(self.rootname, chip))
            log.debug("Mean background: {}".format(bkg_mean))
            log.debug("Mean threshold: {}".format(np.mean(threshold)))
            log.debug("Mean RMS      : {}".format(bkg_rms_mean))
            log.debug("")
            log.debug("{}".format("=" * 60))

            self.bkg[chip] = bkg.background
            self.bkg_dao_rms[chip] = bkg_dao_rms
            self.bkg_rms_mean[chip] = bkg_rms_mean
            self.threshold[chip] = threshold

    def build_dqmask(self, chip=None):
        # apply any DQ array, if available
        dqmask = None
        if chip:
            dqarr = self.imghdu[('DQ', chip)].data
        else:
            dqarr = np.concatenate([self.imghdu[('DQ', i + 1)].data for i in range(self.num_sci)])

        # "grow out" regions in DQ mask flagged as saturated by several
        # pixels in every direction to prevent the
        # source match algorithm from trying to match multiple sources
        # from one image to a single source in the
        # other or vice-versa.
        # Create temp DQ mask containing all pixels flagged with any value EXCEPT 256
        non_sat_mask = bitfield_to_boolean_mask(dqarr, ignore_flags=256+2048)

        # Create temp DQ mask containing saturated pixels ONLY
        sat_mask = bitfield_to_boolean_mask(dqarr, ignore_flags=~(256+2048))

        # Ignore sources where only a couple of pixels are flagged as saturated
        sat_mask = ndimage.binary_erosion(sat_mask, iterations=1)

        # Grow out saturated pixels by a few pixels in every direction
        grown_sat_mask = ndimage.binary_dilation(sat_mask, iterations=5)

        # combine the two temporary DQ masks into a single composite DQ mask.
        dqmask = np.bitwise_or(non_sat_mask, grown_sat_mask)
        return dqmask

    def find_alignment_sources(self, output=True, dqname='DQ', **alignment_pars):
        """Find sources in all chips for this exposure."""

        for chip in range(self.num_sci):
            chip += 1
            # find sources in image
            if output:
                outroot = '{}_sci{}_src'.format(self.rootname, chip)
            else:
                outroot = None

            dqmask = self.build_dqmask(chip=chip)
            sciarr = self.imghdu[("SCI", chip)].data
            #  TODO: replace detector_pars with dict from OO Config class
            extract_pars = {'classify': alignment_pars['classify'],
                            'centering_mode': alignment_pars['centering_mode'],
                            'nlargest': alignment_pars['num_sources'],
                            'deblend': alignment_pars['deblend']}

            seg_tab, segmap = amutils.extract_sources(sciarr, dqmask=dqmask,
                                                      outroot=outroot,
                                                      kernel=self.kernel,
                                                      segment_threshold=self.threshold[chip],
                                                      dao_threshold=self.bkg_rms_mean[chip],
                                                      fwhm=self.kernel_fwhm,
                                                      **extract_pars)

            self.catalog_table[chip] = seg_tab
# ----------------------------------------------------------------------------------------------------------------------


def match_relative_fit(imglist, reference_catalog, **fit_pars):
    """Perform cross-matching and final fit using relative matching algorithm

    Parameters
    ----------
    imglist : list
        List of input image `~tweakwcs.tpwcs.FITSWCS` objects with metadata and source catalogs

    reference_catalog : Table
        Astropy Table of reference sources for this field

    Returns
    --------
    imglist : list
        List of input image `~tweakwcs.tpwcs.FITSWCS` objects with metadata and source catalogs

    """
    log.info("{} (match_relative_fit) Cross matching and fitting {}".format("-" * 20, "-" * 27))
    # 0: Specify matching algorithm to use
    match = tweakwcs.TPMatch(**fit_pars)
    # match = tweakwcs.TPMatch(searchrad=250, separation=0.1,
    #                          tolerance=100, use2dhist=False)

    # Align images and correct WCS
    # NOTE: this invocation does not use an astrometric catalog. This call allows all the input images to be aligned in
    # a relative way using the first input image as the reference.
    # 1: Perform relative alignment
    tweakwcs.align_wcs(imglist, None, match=match, expand_refcat=True)

    # Set all the group_id values to be the same so the various images/chips will be aligned to the astrometric
    # reference catalog as an ensemble.
    # astrometric reference catalog as an ensemble. BEWARE: If additional iterations of solutions are to be
    # done, the group_id values need to be restored.
    for image in imglist:
        image.meta["group_id"] = 1234567
    # 2: Perform absolute alignment
    tweakwcs.align_wcs(imglist, reference_catalog, match=match)

    # 3: Interpret RMS values from tweakwcs
    interpret_fit_rms(imglist, reference_catalog)

    return imglist

# ----------------------------------------------------------------------------------------------------------


def match_default_fit(imglist, reference_catalog, **fit_pars):
    """Perform cross-matching and final fit using default tolerance matching

    Parameters
    ----------
    imglist : list
        List of input image `~tweakwcs.tpwcs.FITSWCS` objects with metadata and source catalogs

    reference_catalog : Table
        Astropy Table of reference sources for this field

    Returns
    --------
    imglist : list
        List of input image `~tweakwcs.tpwcs.FITSWCS` objects with metadata and source catalogs

    """
    log.info("{} (match_default_fit) Cross matching and fitting "
             "{}".format("-" * 20, "-" * 27))
    # Specify matching algorithm to use
    match = tweakwcs.TPMatch(**fit_pars)

    # Align images and correct WCS
    tweakwcs.align_wcs(imglist, reference_catalog, match=match, expand_refcat=False)

    # Interpret RMS values from tweakwcs
    interpret_fit_rms(imglist, reference_catalog)

    return imglist


# ----------------------------------------------------------------------------------------------------------------------


def match_2dhist_fit(imglist, reference_catalog, **fit_pars):
    """Perform cross-matching and final fit using 2dHistogram matching

    Parameters
    ----------
    imglist : list
        List of input image `~tweakwcs.tpwcs.FITSWCS` objects with metadata and source catalogs

    reference_catalog : Table
        Astropy Table of reference sources for this field

    Returns
    --------
    imglist : list
        List of input image `~tweakwcs.tpwcs.FITSWCS` objects with metadata and source catalogs

    """
    log.info("{} (match_2dhist_fit) Cross matching and fitting "
             "{}".format("-" * 20, "-" * 28))
    # Specify matching algorithm to use
    match = tweakwcs.TPMatch(**fit_pars)
    # Align images and correct WCS
    tweakwcs.align_wcs(imglist, reference_catalog, match=match, expand_refcat=False)

    # Interpret RMS values from tweakwcs
    interpret_fit_rms(imglist, reference_catalog)

    return imglist

# ----------------------------------------------------------------------------------------------------------
def interpret_fit_rms(tweakwcs_output, reference_catalog):
    """Interpret the FIT information to convert RMS to physical units

    Parameters
    ----------
    tweakwcs_output : list
        output of tweakwcs. Contains sourcelist tables, newly computed WCS info, etc. for every chip of every valid
        input image.  This list gets updated, in-place, with the new RMS values;
        specifically,

            * 'FIT_RMS': RMS of the separations between fitted image positions and reference positions
            * 'TOTAL_RMS': mean of the FIT_RMS values for all observations
            * 'NUM_FITS': number of images/group_id's with successful fits included in the TOTAL_RMS

        These entries are added to the 'fit_info' dictionary.

    reference_catalog : astropy.Table
        Table of reference source positions used for the fit

    Returns
    -------
    Nothing
    """
    # Start by collecting information by group_id
    group_ids = [info.meta['group_id'] for info in tweakwcs_output]
    # Compress the list to have only unique group_id values to avoid some unnecessary iterations
    group_ids = list(set(group_ids))
    group_dict = {'avg_RMS': None}
    obs_rms = []
    for group_id in group_ids:
        input_mag = None
        for item in tweakwcs_output:
            # When status = FAILED (fit failed) or REFERENCE (relative alignment done with first image
            # as the reference), skip to the beginning of the loop as there is no 'fit_info'.
            if item.meta['fit_info']['status'] != 'SUCCESS':
                continue
            # Make sure to store data for any particular group_id only once.
            if item.meta['group_id'] == group_id and \
               group_id not in group_dict:
                group_dict[group_id] = {'ref_idx': None, 'FIT_RMS': None,
                                        'input_mag': None, 'ref_mag': None, 'input_idx': None}

                log.debug("fit_info: {}".format(item.meta['fit_info']))

                tinfo = item.meta['fit_info']
                ref_idx = tinfo['matched_ref_idx']
                fitmask = tinfo['fitmask']
                group_dict[group_id]['ref_idx'] = ref_idx
                ref_RA = reference_catalog[ref_idx]['RA'][fitmask]
                ref_DEC = reference_catalog[ref_idx]['DEC'][fitmask]
                input_RA = tinfo['fit_RA']
                input_DEC = tinfo['fit_DEC']
                img_coords = SkyCoord(input_RA, input_DEC,
                                      unit='deg', frame='icrs')
                ref_coords = SkyCoord(ref_RA, ref_DEC, unit='deg', frame='icrs')
                dra, ddec = img_coords.spherical_offsets_to(ref_coords)
                ra_rms = np.std(dra.to(u.mas))
                dec_rms = np.std(ddec.to(u.mas))
                fit_rms = np.std(Angle(img_coords.separation(ref_coords), unit=u.mas)).value
                group_dict[group_id]['FIT_RMS'] = fit_rms
                group_dict[group_id]['RMS_RA'] = ra_rms
                group_dict[group_id]['RMS_DEC'] = dec_rms

                group_dict[group_id]['ref_mag'] = reference_catalog[ref_idx]['mag'][fitmask]

                input_mag = item.meta['catalog']['mag']
                group_dict[group_id]['input_mag'] = input_mag
                group_dict[group_id]['input_idx'] = tinfo['matched_input_idx']

                obs_rms.append(fit_rms)

            else:
                if input_mag is not None:
                    input_mag = input_mag.copy(data=np.hstack((input_mag, item.meta['catalog']['mag'])))
                    group_dict[group_id]['input_mag'] = input_mag


    # Compute RMS for entire ASN/observation set
    total_rms = np.mean(obs_rms)

    # Now, append computed results to tweakwcs_output
    for item in tweakwcs_output:
        group_id = item.meta['group_id']
        fitmask = item.meta['fit_info']['fitmask']
        if group_id in group_dict:
            fit_rms = group_dict[group_id]['FIT_RMS']
            ra_rms = group_dict[group_id]['RMS_RA']
            dec_rms = group_dict[group_id]['RMS_DEC']
            input_idx = group_dict[group_id]['input_idx']
            input_mag = group_dict[group_id]['input_mag'][input_idx][fitmask]
            ref_mag = group_dict[group_id]['ref_mag']

        else:
            fit_rms = None
            ra_rms = None
            dec_rms = None
            input_mag = None
            ref_mag = None

        item.meta['fit_info']['FIT_RMS'] = fit_rms
        item.meta['fit_info']['TOTAL_RMS'] = total_rms
        item.meta['fit_info']['NUM_FITS'] = len(group_ids)
        item.meta['fit_info']['RMS_RA'] = ra_rms
        item.meta['fit_info']['RMS_DEC'] = dec_rms
        item.meta['fit_info']['catalog'] = reference_catalog.meta['catalog']
        item.meta['fit_info']['input_mag'] = input_mag
        item.meta['fit_info']['ref_mag'] = ref_mag

# ----------------------------------------------------------------------------------------------------------------------


def update_image_wcs_info(tweakwcs_output, headerlet_filenames=None, fit_label=None):
    """Write newly computed WCS information to image headers and write headerlet files

        Parameters
        ----------
        tweakwcs_output : list
            output of tweakwcs. Contains sourcelist tables, newly computed WCS info, etc. for every chip of every valid
            every valid input image.

        headerlet_filenames : dictionary, optional
            dictionary that maps the flt/flc.fits file name to the corresponding custom headerlet filename.

        Returns
        -------
        out_headerlet_list : dictionary
            a dictionary of the headerlet files created by this subroutine, keyed by flt/flc fits filename.
        """
    out_headerlet_dict = {}
    for item in tweakwcs_output:
        image_name = item.meta['filename']
        chipnum = item.meta['chip']
        if chipnum == 1:
            chipctr = 1
            hdulist = fits.open(image_name, mode='update')
            num_sci_ext = amutils.countExtn(hdulist)

            # generate wcs name for updated image header, headerlet
            # Just in case header value 'wcs_name' is empty.
            if fit_label is None:
                if 'relative' in item.meta['fit method']:
                    fit_label = 'REL'
                else:
                    fit_label = 'IMG'

            if not hdulist['SCI', 1].header['WCSNAME'] or hdulist['SCI', 1].header['WCSNAME'] == "":
                wcs_name = "FIT_{}_{}".format(fit_label, item.meta['catalog_name'])
            else:
                wname = hdulist['sci', 1].header['wcsname']
                if "-" in wname:
                    wcs_name = '{}-FIT_{}_{}'.format(wname[:wname.index('-')],
                                                    fit_label,
                                                    item.meta['fit_info']['catalog'])
                else:
                    wcs_name = '{}-FIT_{}_{}'.format(wname, fit_label, item.meta['fit_info']['catalog'])

            # establish correct mapping to the science extensions
            sci_ext_dict = {}
            for sci_ext_ctr in range(1, num_sci_ext + 1):
                sci_ext_dict["{}".format(sci_ext_ctr)] = fileutil.findExtname(hdulist, 'sci', extver=sci_ext_ctr)

        # update header with new WCS info
        updatehdr.update_wcs(hdulist, sci_ext_dict["{}".format(item.meta['chip'])], item.wcs,
                             wcsname=wcs_name,
                             reusename=True)
        if chipctr == num_sci_ext:
            # Close updated flc.fits or flt.fits file
            hdulist.flush()
            hdulist.close()

            # Create headerlet
            out_headerlet = headerlet.create_headerlet(image_name, hdrname=wcs_name, wcsname=wcs_name,
                                                       logging=False)

            # Update headerlet
            update_headerlet_phdu(item, out_headerlet)

            # Write headerlet
            if headerlet_filenames:
                headerlet_filename = headerlet_filenames[image_name]  # Use HAP-compatible filename defined in runhlaprocessing.py
            else:
                if image_name.endswith("flc.fits"):
                    headerlet_filename = image_name.replace("flc", "flt_hlet")
                if image_name.endswith("flt.fits"):
                    headerlet_filename = image_name.replace("flt", "flt_hlet")
            out_headerlet.writeto(headerlet_filename, clobber=True)
            log.info("Wrote headerlet file {}.\n\n".format(headerlet_filename))
            out_headerlet_dict[image_name] = headerlet_filename

            # Attach headerlet as HDRLET extension
            headerlet.attach_headerlet(image_name, headerlet_filename, logging=False)

        chipctr += 1
    return (out_headerlet_dict)

# --------------------------------------------------------------------------------------------------------------
def update_headerlet_phdu(tweakwcs_item, headerlet):
    """Update the primary header data unit keywords of a headerlet object in-place

    Parameters
    ==========
    tweakwcs_item :
        Basically the output from tweakwcs which contains the cross match and fit information for every chip
        of every valid input image.

    headerlet : headerlet object
        object containing WCS information
    """

    # Get the data to be used as values for FITS keywords
    rms_ra = tweakwcs_item.meta['fit_info']['RMS_RA'].value
    rms_dec = tweakwcs_item.meta['fit_info']['RMS_DEC'].value
    fit_rms = tweakwcs_item.meta['fit_info']['FIT_RMS']
    nmatch = tweakwcs_item.meta['fit_info']['nmatches']
    catalog = tweakwcs_item.meta['fit_info']['catalog']
    fit_method = tweakwcs_item.meta['fit method']

    x_shift = (tweakwcs_item.meta['fit_info']['shift'])[0]
    y_shift = (tweakwcs_item.meta['fit_info']['shift'])[1]
    rot = tweakwcs_item.meta['fit_info']['rot']
    scale = tweakwcs_item.meta['fit_info']['scale'][0]
    skew = tweakwcs_item.meta['fit_info']['skew']

    # Update the existing FITS keywords
    primary_header = headerlet[0].header
    primary_header['RMS_RA'] = rms_ra
    primary_header['RMS_DEC'] = rms_dec
    primary_header['NMATCH'] = nmatch
    primary_header['CATALOG'] = catalog
    primary_header['FITMETH'] = fit_method

    # Create a new FITS keyword
    primary_header['FIT_RMS'] = (fit_rms, 'RMS (mas) of the 2D fit of the headerlet solution')

    # Create the set of HISTORY keywords
    primary_header['HISTORY'] = '~~~~~ FIT PARAMETERS ~~~~~'
    primary_header['HISTORY'] = '{:>15} : {:9.4f} "/pixels'.format('platescale', tweakwcs_item.wcs.pscale)
    primary_header['HISTORY'] = '{:>15} : {:9.4f} pixels'.format('x_shift', x_shift)
    primary_header['HISTORY'] = '{:>15} : {:9.4f} pixels'.format('y_shift', y_shift)
    primary_header['HISTORY'] = '{:>15} : {:9.4f} degrees'.format('rotation', rot)
    primary_header['HISTORY'] = '{:>15} : {:9.4f}'.format('scale', scale)
    primary_header['HISTORY'] = '{:>15} : {:9.4f}'.format('skew', skew)

# --------------------------------------------------------------------------------------------------------------
def register_photutils_function(name):
    """Convert photutils name as a string into a pointer to the actual photutils function"""

    if name in dir(photutils):
        func = eval("photutils.{}".format(name))
    return func
