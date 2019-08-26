import os
import sys
import datetime

import numpy as np
from scipy import ndimage

from astropy.io import fits
from astropy.table import Table
from astropy.nddata.bitmask import bitfield_to_boolean_mask
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u

from photutils import Background2D, SExtractorBackground, StdBackgroundRMS

from stwcs.wcsutil import HSTWCS
from stsci.tools import logutil

from . import astrometric_utils as amutils

from .. import tweakwcs

__taskname__ = 'align_utils'

# Default background determination parameter values
BKG_BOX_SIZE = 50
BKG_FILTER_SIZE = 3
CATALOG_TYPES = ['point', 'segment']

log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout)


class AlignmentTable:
    def __init__(self, input_list, clobber=False, dqname='DQ', **alignment_pars):
        self.alignment_pars = alignment_pars
        self.dqname = dqname

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

        self.haplist = []
        for img in self.process_list:
            catimg = HAPImage(img)
            # Build image properties needed for alignment
            catimg.compute_background(**alignment_pars)
            catimg.build_kernel(self.alignment_pars.get('fwhmpsf'))

            self.haplist.append(catimg)

        # Initialize computed attributes
        self.imglist = []  # list of FITSWCS objects for tweakwcs
        self.reference_catalogs = {}
        self.group_id_dict = {}

        self.fit_methods = {'relative': match_relative_fit,
                            '2dhist': match_2dhist_fit,
                            'default': match_default_fit}

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

    def reset_group_id(self):
        for image in self.imglist:
            image.meta["group_id"] = self.group_id_dict["{}_{}".format(image.meta["filename"], image.meta["chip"])]

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
                log.info('im.meta[name] = {}'.format(im.meta['name']))
            self.imglist.extend(img)

        self.group_id_dict = {}
        for image in self.imglist:
            self.group_id_dict["{}_{}".format(image.meta["filename"], image.meta["chip"])] = image.meta["group_id"]

    def get_fit_methods(self):
        """Return the list of method names for all registered functions
            available for performing alignment.
        """
        return self.fit_methods.keys()

    def perform_fit(self, method_name, catalog_name, reference_catalog):
        """Perform fit using specified method, then determine fit quality"""
        imglist = self.fit_methods[method_name](self.imglist, reference_catalog)

        # store results for evaluation
        self.fit_dict[(catalog_name, method_name)] = copy.deepcopy(imglist)
        self.reference_catalogs[catalog_name] = reference_catalog

        return imglist

    def select_fit(self, catalog_name, method_name):
        """Select the fit that has been identified as 'best'"""
        imglist = self.selected_fit = self.fit_list[(method_name, reference_catalog)]

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
                self.filtered_table[index]['catalogSources'] = len(reference_catalog)
                self.filtered_table[index]['matchSources'] = item.meta['fit_info']['nmatches']
                self.filtered_table[index]['rms_ra'] = item.meta['fit_info']['RMS_RA'].value
                self.filtered_table[index]['rms_dec'] = item.meta['fit_info']['RMS_DEC'].value
                self.filtered_table[index]['fit_rms'] = item.meta['fit_info']['FIT_RMS']
                self.filtered_table[index]['total_rms'] = item.meta['fit_info']['TOTAL_RMS']
                self.filtered_table[index]['offset_x'], self.filtered_table[index]['offset_y'] = item.meta['fit_info']['shift']
                self.filtered_table[index]['scale'] = item.meta['fit_info']['scale'][0]
                self.filtered_table[index]['rotation'] = item.meta['fit_info']['rot']

                # populate self.filtered_table fields "status", "compromised" and
                # "processMsg" with fit_status_dict fields "valid", "compromised"
                # and "reason".
                explicit_dict_key = "{},{}".format(item.meta['name'], item.meta['chip'])
                if fit_status_dict[explicit_dict_key]['valid'] is True:
                    self.filtered_table[index]['status'] = 0
                else:
                    self.filtered_table[index]['status'] = 1
                if fit_status_dict[explicit_dict_key]['compromised'] is False:
                    self.filtered_table['compromised'] = 0
                else:
                    self.filtered_table['compromised'] = 1

                self.filtered_table[index]['processMsg'] = fit_status_dict[explicit_dict_key]['reason']
                self.filtered_table['fit_qual'][index] = item.meta['fit quality']


    def apply_fit(self, headerlet_filenames=None):
        """Apply solution from identified fit to image WCS's

        Parameters
        ----------
        headerlet_filenames : dict, optional
            Dictionary relating exposure filenames to headerlet filenames.  If None,
            will generate headerlet filenames where _flt or _flc is replaced by
            _flt_hlet or _flc_hlet, respectively.

        """
        if not self.selected_fit:
            print("No FIT selected for application.  Please run 'select_fit()' method.")
            raise ValueError
        # Call update_hdr_wcs()
        headerlet_dict = align_utils.update_image_wcs_info(self.selected_fit, headerlet_filenames=headerlet_filenames)

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

        # Fits file read
        self.num_sci = amutils.countExtn(self.imghdu)
        self.num_wht = amutils.countExtn(self.imghdu, extname='WHT')
        self.data = np.concatenate([self.imghdu[('SCI', i + 1)].data for i in range(self.num_sci)])
        if not self.num_wht:
            self.dqmask = build_dqmask()
        else:
            self.dqmask = None

        # Get the HSTWCS object from the first extension
        self.imgwcs = HSTWCS(self.imghdu, 1)
        self.pscale = self.imgwcs.pscale

        if 'rootname' in self.imghdu[0].header:
            self.rootname = self.imghdu[0].header['rootname']
        else:
            self.rootname = self.imgname.rstrip('.fits')

        self.bkg = None
        self.kernel = None
        self.kernel_fwhm = None
        self.fwhmpsf = None

        self.catalog_table = {}

    @def wht_image():
        doc = "The wht_image property."
        def fget(self):
            return self.wht_image
        def fset(self):
            if not self.num_wht:
                # Working with a calibrated exposure, no WHT extension
                # create a substitute WHT array from ERR and DQ
                # Build pseudo-wht array for detection purposes
                errarr = np.concatenate([self.imghdu[('ERR', i + 1)].data for i in range(self.num_sci)])
                wht_image = 1.0 / errarr
                wht_image /= wht_image.max()
                wht_image *= self.imghdu[0].header['exptime']**2
                wht_image[dqmask] = 0
                self.wht_image = wht_image
            else:
                self.wht_image = self.imghdu['WHT'].data
        def fdel(self):
            del self.wht_image
        return locals()
    wht_image = property(**wht_image())

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

        self.kernel, self.kernel_fwhm = amutils.build_auto_kernel(self.data, self.wht_image,
                                                          threshold=self.bkg.background_rms, fwhm=fwhmpsf / self.pscale)
        self.fwhmpsf = self.kernel_fwhm * self.pscale

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
        non_sat_mask = bitfield_to_boolean_mask(dqarr, ignore_flags=256)

        # Create temp DQ mask containing saturated pixels ONLY
        sat_mask = bitfield_to_boolean_mask(dqarr, ignore_flags=~256)

        # Grow out saturated pixels by a few pixels in every direction
        grown_sat_mask = ndimage.binary_dilation(sat_mask, iterations=5)

        # combine the two temporary DQ masks into a single composite DQ mask.
        dqmask = np.bitwise_or(non_sat_mask, grown_sat_mask)

        self.dqmask = dqmask

    def find_alignment_sources(self, output=True, dqname='DQ', **alignment_pars):
        """Find sources in all chips for this exposure."""
        for chip in range(self.numSci):
            chip += 1
            # find sources in image
            if output:
                outroot = '{}_sci{}_src'.format(self.rootname, chip)
            else:
                outroot = None

            dqmask = build_dqmask(chip=chip)
            sciarr = self.imghdu[("SCI", chip)].data
            #  TODO: replace detector_pars with dict from OO Config class
            seg_tab, segmap = amutils.extract_sources(sciarr, dqmask=self.dqmask,
                                                      outroot=outroot,
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

# ----------------------------------------------------------------------------------------------------------------------


def generate_astrometric_catalog(imglist, **pars):
    """Generates a catalog of all sources from an existing astrometric catalog are
       in or near the FOVs of the images in the input list.

    Parameters
    ----------
    imglist : list
        List of one or more calibrated fits images that will be used for catalog generation.

    Returns
    =======
    ref_table : object
        Astropy Table object of the catalog

    """
    # generate catalog
    temp_pars = pars.copy()
    if pars['output'] is True:
        pars['output'] = 'ref_cat.ecsv'
    else:
        pars['output'] = None
    out_catalog = amutils.create_astrometric_catalog(imglist, **pars)
    pars = temp_pars.copy()
    # if the catalog has contents, write the catalog to ascii text file
    if len(out_catalog) > 0 and pars['output']:
        catalog_filename = "refcatalog.cat"
        out_catalog.write(catalog_filename, format="ascii.fast_commented_header")
        log.info("Wrote reference catalog {}.".format(catalog_filename))

    return(out_catalog)

# ----------------------------------------------------------------------------------------------------------------------


def match_relative_fit(imglist, reference_catalog):
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
    log.info("{} STEP 5b: (match_relative_fit) Cross matching and fitting {}".format("-" * 20, "-" * 27))
    # 0: Specify matching algorithm to use
    match = tweakwcs.TPMatch(searchrad=75, separation=0.1, tolerance=2, use2dhist=True)
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


def match_default_fit(imglist, reference_catalog):
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
    log.info("{} STEP 5b: (match_default_fit) Cross matching and fitting "
             "{}".format("-" * 20, "-" * 27))
    # Specify matching algorithm to use
    match = tweakwcs.TPMatch(searchrad=250, separation=0.1,
                             tolerance=100, use2dhist=False)
    # Align images and correct WCS
    tweakwcs.align_wcs(imglist, reference_catalog, match=match, expand_refcat=False)

    # Interpret RMS values from tweakwcs
    interpret_fit_rms(imglist, reference_catalog)

    return imglist


# ----------------------------------------------------------------------------------------------------------------------


def match_2dhist_fit(imglist, reference_catalog):
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
    log.info("{} STEP 5b: (match_2dhist_fit) Cross matching and fitting "
             "{}".format("-" * 20, "-" * 28))
    # Specify matching algorithm to use
    match = tweakwcs.TPMatch(searchrad=75, separation=0.1, tolerance=2.0, use2dhist=True)
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
        for item in tweakwcs_output:
            # When status = FAILED (fit failed) or REFERENCE (relative alignment done with first image
            # as the reference), skip to the beginning of the loop as there is no 'fit_info'.
            if item.meta['fit_info']['status'] != 'SUCCESS':
                continue
            # Make sure to store data for any particular group_id only once.
            if item.meta['group_id'] == group_id and \
               group_id not in group_dict:
                group_dict[group_id] = {'ref_idx': None, 'FIT_RMS': None}
                log.info("fit_info: {}".format(item.meta['fit_info']))
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

                obs_rms.append(fit_rms)
    # Compute RMS for entire ASN/observation set
    total_rms = np.mean(obs_rms)
    # total_rms = np.sqrt(np.sum(np.array(obs_rms)**2))

    # Now, append computed results to tweakwcs_output
    for item in tweakwcs_output:
        group_id = item.meta['group_id']
        if group_id in group_dict:
            fit_rms = group_dict[group_id]['FIT_RMS']
            ra_rms = group_dict[group_id]['RMS_RA']
            dec_rms = group_dict[group_id]['RMS_DEC']
        else:
            fit_rms = None
            ra_rms = None
            dec_rms = None

        item.meta['fit_info']['FIT_RMS'] = fit_rms
        item.meta['fit_info']['TOTAL_RMS'] = total_rms
        item.meta['fit_info']['NUM_FITS'] = len(group_ids)
        item.meta['fit_info']['RMS_RA'] = ra_rms
        item.meta['fit_info']['RMS_DEC'] = dec_rms
        item.meta['fit_info']['catalog'] = reference_catalog.meta['catalog']

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
                if item.meta['fit method'] == 'match_relative_fit':
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
        updatehdr.update_wcs(hdulist, sci_ext_dict["{}".format(item.meta['chip'])], item.wcs, wcsname=wcs_name,
                                 reusename=True, verbose=True)
        if chipctr == num_sci_ext:
            # Close updated flc.fits or flt.fits file
            hdulist.flush()
            hdulist.close()

            # Create headerlet
            out_headerlet = headerlet.create_headerlet(image_name, hdrname=wcs_name, wcsname=wcs_name, logging=False)

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
