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
from astropy.coordinates import SkyCoord, Angle

from photutils import Background2D, SExtractorBackground, StdBackgroundRMS

from stwcs.wcsutil import HSTWCS
from stsci.tools import logutil

from . import astrometric_utils as amutils
from . import astroquery_utils as aqutils

from .. import tweakwcs

__taskname__ = 'align_utils'

# Default background determination parameter values
BKG_BOX_SIZE = 50
BKG_FILTER_SIZE = 3
CATALOG_TYPES = ['point', 'segment']

log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout)


class AlignmentTable:
    def __init__(self, input_list, clobber=False, **alignment_pars):
        self.alignment_pars = alignment_pars

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
        self.reference_catalogs = []
        self.group_id_dict = {}

        self.fit_methods = {'relative': match_relative_fit,
                            '2dhist': match_2dhist_fit,
                            'default': match_default_fit}

    def build_images(self, **image_pars):
        """Create CatalogImage objects for each input to use in alignment."""
        if self.imglist:
            return
        fwhmpsf = self.alignment_pars.get('fwhmpsf')

        for img in self.process_list:
            catimg = CatalogImage(img)
            # Build image properties needed for alignment
            catimg.compute_background(**image_pars)
            catimg.build_kernel(fwhmpsf)

            self.imglist.append(catimg)


    def get_reference_catalog(self, catalog_name, output=True, catalog=None):
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

        output : boolean
            Specify whether or not to write out reference catalog to a file

        catalog : Table, optional
            Astrometric catalog to use for alignment, if specified.  It will be
            labelled with value given in `catalog_name` and added as an entry
            in `self.reference_catalogs`.

        Returns
        -------
        reference_catalog : Table
            Table with astrometric source positions to use for alignment.  This
            table will be saved in `self.reference_catalogs` dictionary.

        """
        if catalog_name in self.reference_catalogs:
            log.info("Using {} reference catalog from earlier this run.".format(catalog_name))
            reference_catalog = self.reference_catalogs[catalog_name]
        else:
            if catalog:
                reference_catalog = catalog
                log.info("Using custom reference catalog {};"
                         " Storing it for potential re-use later this run.".format(catalog_name))
            else:
                log.info("Generating new reference catalog for {};"
                         " Storing it for potential re-use later this run.".format(catalog_name))
                reference_catalog = generate_astrometric_catalog(self.process_list,
                                                                 catalog=catalog_name,
                                                                 output=output)
            self.reference_catalogs[catalog_name] = reference_catalog
        return reference_catalog

    def find_sources(self, output=True, dqname='DQ', fwhmpsf=3.0, **alignment_pars):
        """Find observable sources in each input exposure."""
        self.extracted_sources = {}
        for img in self.imglist:
            img.find_sources(output=output, dqname=dqname, fwhmpsf=fwhmpsf, **alignment_pars)
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

    def perform_fit(self, method_name, reference_catalog):
        """Perform fit using specified method, then determine fit quality"""
        self.fit_methods[method_name](self.imglist, reference_catalog)

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
            else:
                outroot = None

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
