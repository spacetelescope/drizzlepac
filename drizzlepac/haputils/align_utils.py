import os
import datetime
import copy
import sys
import traceback
import warnings

from distutils.version import LooseVersion
from collections import OrderedDict

import numpy as np
from scipy import ndimage

import astropy
from astropy.io import fits
from astropy.table import Table
from stsci.tools.bitmask import bitfield_to_boolean_mask
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u

from photutils import background
from photutils.background import Background2D
from photutils.utils import NoDetectionsWarning

import stwcs
from stwcs.wcsutil import HSTWCS
from stwcs.wcsutil import headerlet

import stsci.tools
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
    """ This class handles alignment operations for HST data.

    Steps managed by this class include:
        * **find_alignment_sources** : iterates over all inputs and
          identifies sources suitable for alignment.
        * **perform_fit** : cross-matches sources using user-specified method
          and performs fit to user-specified catalog.  Cross-matching options
          include: `relative`, '2dhist' or 'default' as reported by the ``get_fit_methods``
          method.  It then populates the table with the results for this fit.
        * **select_fit** : extracts the 'best' fit to be applied to the data with the
          'best' fit determined externally by the user based on the set of fit results
          stored in the `fit_dict` table.
        * **apply_fit** : Updates all input image WCSs with the result of the selected 'best' fit

    """
    def __init__(self, input_list, clobber=False, dqname='DQ',
                 log_level=logutil.logging.NOTSET, **alignment_pars):
        """
        Parameters
        ----------
        input_list : list
            List of input file names to be provided to `~analyze.analyze_data` function.

        clobber : bool, optional
            Specifies whether or not to overwrite data on disk with updated versions of the data.

        dqname : str, optional
            Allows the user to customize the name of the extension (`extname`) containing the
            data quality flags to be applied to the data during source identification.

        log_level : int, optional
            Set the logging level for this processing

        alignment_pars : dict
            Set of alignment parameters to be used for the input data.
            ``**alignment_pars`` needs to contain the following entries:

            .. code-block:: python

                {'fwhmpsf': 0.12,  # kernel defining, source finding par
                 # background computing pars
                 'box_size': BKG_BOX_SIZE,
                 'win_size': BKG_FILTER_SIZE,
                 'bkg_estimator': SExtractorBackground,
                 'rms_estimator': StdBackgroundRMS,
                 'nsigma': 5.,
                 'threshold_flag': None,
                 # object finding pars
                 'source_box': 7,
                 'classify': True,
                 'centering_mode': "starfind",
                 'nlargest': None,
                 'plot': False,
                 'vmax': None,
                 'deblend': False
                }

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
        # Make sure to use the values from "determine_fit_quality" instead of "general"
        for key in self.fit_pars:
            self.fit_pars[key]['pars'] = alignment_pars['determine_fit_quality']
            self.fit_pars[key]['pars']['minobj'] = self.alignment_pars['mosaic_fitgeom_list']

        self.dqname = dqname
        self.haplist = []
        self.process_list = None

        self.zero_dt = starting_dt = datetime.datetime.now()
        log.info(str(starting_dt))
        # Apply filter to input observations to insure that they meet minimum criteria for being able to be aligned
        log.info(
            "{} AlignmentTable: Filter STEP {}".format("-" * 20, "-" * 63))
        self.filtered_table = analyze.analyze_data(input_list)
        log.debug("Input sorted as: \n{}".format(self.filtered_table))

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
            for img in self.process_list:
                log.info("Adding {} to HAPImage list".format(img))
                catimg = HAPImage(img)
                # Build image properties needed for alignment
                catimg.compute_background(box_size=self.alignment_pars['box_size'],
                                          win_size=self.alignment_pars['win_size'],
                                          nsigma=self.alignment_pars['nsigma'],
                                          bkg_estimator=self.alignment_pars['bkg_estimator'],
                                          rms_estimator=self.alignment_pars['rms_estimator'],
                                          threshold_flag=self.alignment_pars['threshold'])
                catimg.build_kernel(fwhmpsf)
                catimg.crclean = self.alignment_pars['classify']
                log.info("CATIMG.CRCLEAN: {}".format(catimg.crclean))
                if catimg.crclean:
                    log.info("Recomputing BACKGROUND after applying single-image CR clean")
                    for chip in range(1, catimg.num_sci + 1):
                        sciarr = amutils.crclean_image(catimg.imghdu[("SCI", chip)].data, catimg.threshold[chip],
                                                       catimg.kernel, catimg.kernel_fwhm, background=catimg.bkg[chip])
                        catimg.imghdu[("SCI", chip)].data = sciarr
                    # recompute now that the CRs have been removed (set to 0) from the science arrays
                    catimg.compute_background(box_size=self.alignment_pars['box_size'],
                                              win_size=self.alignment_pars['win_size'],
                                              nsigma=self.alignment_pars['nsigma'],
                                              bkg_estimator=self.alignment_pars['bkg_estimator'],
                                              rms_estimator=self.alignment_pars['rms_estimator'],
                                              threshold_flag=self.alignment_pars['threshold'])
                    log.info("Finished computing revised BACKROUND")
                    catimg.build_kernel(fwhmpsf)
                    log.info("Finished determining revised kernel")

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


    def find_alignment_sources(self, output=True, crclean=None):
        """Find observable sources in each input exposure."""
        if crclean is None:
            crclean = [False] * len(self.haplist)

        self.extracted_sources = {}
        for img, clean in zip(self.haplist, crclean):
            self.extracted_sources[img.imgname] = {}
            if img.imghdu is not None:
                img.find_alignment_sources(output=output, dqname=self.dqname,
                                           crclean=clean,
                                           **self.alignment_pars)
                self.extracted_sources[img.imgname] = img.catalog_table

                # Allow user to decide when and how to write out catalogs to files
                if output:
                    # write out coord lists to files for diagnostic purposes.
                    # Protip: To display the sources in these files in DS9,
                    #         set the "Coordinate System" option to "Physical"
                    #         when loading the region file.
                    imgroot = "_".join(os.path.basename(img.imgname).split('_')[:-1])
                    for chip in range(1, img.num_sci + 1):
                        chip_cat = img.catalog_table[chip]
                        if chip_cat and len(chip_cat) > 0:
                            regfilename = "{}_sci{}_src.reg".format(imgroot, chip)
                            out_table = Table(chip_cat)
                            # To align with positions of sources in DS9/IRAF
                            out_table['xcentroid'] += 1
                            out_table['ycentroid'] += 1
                            out_table.write(regfilename,
                                            include_names=["xcentroid", "ycentroid", "mag"],
                                            format="ascii.fast_commented_header")
                            log.info("Wrote region file {}\n".format(regfilename))

        self.close()

    def reset_group_id(self, num_ref):
        for image in self.imglist:
            image.meta["group_id"] = self.group_id_dict["{}_{}".format(image.meta["rootname"], image.meta["chip"])]
            image.meta['num_ref_catalog'] = num_ref


    def configure_fit(self):
        # Convert input images to tweakwcs-compatible FITSWCS objects and
        # attach source catalogs to them.
        self.imglist = []
        for group_id, image in enumerate(self.process_list):
            if image in self.extracted_sources:
                img = amutils.build_wcscat(image, group_id,
                                           self.extracted_sources[image])
                # add the name of the image to the imglist object
                for im in img:
                #    im.meta['name'] = image
                    log.debug('im.meta[name] = {}'.format(im.meta['name']))
                self.imglist.extend(img)

        self.group_id_dict = {}
        for image in self.imglist:
            self.group_id_dict["{}_{}".format(image.meta["rootname"], image.meta["chip"])] = image.meta["group_id"]

    def get_fit_methods(self):
        """Return the list of method names for all registered functions
            available for performing alignment.
        """
        return list(self.fit_methods.keys())

    def perform_fit(self, method_name, catalog_name, reference_catalog,
                    fitgeom='rscale'):
        """Perform fit using specified method, then determine fit quality

        This method populates the `fit_dict` table with the results of the
        specified fit keyed by (<method_name>, <catalog_name>).

        Parameters
        -----------
        method_name : str
            Name of cross-matching/fitting to use for this fit. Options
            are reported by the ``get_fit_methods`` method.

        catalog_name : str
            Name of reference catalog to use for the fit.  These are defined
            by the user.  Examples include: 'GAIADR1', 'GAIADR2', and 'GAIAedr3'.
            This acts as the label for indexing this fit in the `fit_dict` table.

        reference_catalog : `astropy.table.Table`
            Table containing the reference sources to be used for the fit.

        fitgeom : str, optional
            Type of polynomial fit to perform to determine the correction
            to the WCS.  Options include (from more complex to simplest):
            `general`, `rscale`, `rshift`, `shift`.

        """
        # Updated fits_pars with value for fitgeom
        self.fit_pars[method_name]['fitgeom'] = fitgeom
        log.info("Setting 'fitgeom' parameter to {} for {} fit".format(fitgeom, method_name))

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
        """Select the fit that has been identified as 'best'

        Populates the `filtered_table` with a row for each input exposure
        that contains the results of the fit selected from the `fit_dict` table
        populated by ``perform_dict()``.

        Parameters
        ----------
        catalog_name : str
            Name of reference catalog used for fit to be selected

        method_name : str
            Name of cross-matching used for the selected fit

        """
        if catalog_name is None:
            self.selected_fit = None
            return

        imglist = self.selected_fit = self.fit_dict[(catalog_name, method_name)]
        if imglist[0].meta['fit_info']['status'].startswith("FAILED"):
            self.selected_fit = None
            return

        # Protect the writing of the table within the best_fit_rms
        info_keys = OrderedDict(imglist[0].meta['fit_info']).keys()
        # Update filtered table with number of matched sources and other information
        for item in imglist:
            imgname = item.meta['name']
            index = np.where(self.filtered_table['imageName'] == imgname)[0][0]
            status = item.meta['fit_info']['status']
            if not status.startswith("FAILED"):
                self.filtered_table[index]['catalog'] = item.meta['fit_info']['catalog']
                self.filtered_table[index]['fit_rms'] = item.meta['fit_info']['FIT_RMS']
                self.filtered_table[index]['total_rms'] = item.meta['fit_info']['TOTAL_RMS']
                if status.startswith('REF'):
                    continue

                for tweakwcs_info_key in info_keys:
                    if not tweakwcs_info_key.startswith("matched"):
                        if tweakwcs_info_key.lower() == 'rms':
                            self.filtered_table[index]['rms_x'] = item.meta['fit_info'][tweakwcs_info_key][0]
                            self.filtered_table[index]['rms_y'] = item.meta['fit_info'][tweakwcs_info_key][1]

                self.filtered_table[index]['fit_method'] = item.meta['fit method']
                self.filtered_table[index]['catalogSources'] = len(self.reference_catalogs[catalog_name])
                self.filtered_table[index]['matchSources'] = item.meta['fit_info']['nmatches']
                self.filtered_table[index]['rms_ra'] = item.meta['fit_info']['RMS_RA'].value
                self.filtered_table[index]['rms_dec'] = item.meta['fit_info']['RMS_DEC'].value
                self.filtered_table[index]['offset_x'], self.filtered_table[index]['offset_y'] = item.meta['fit_info']['shift']
                self.filtered_table[index]['scale'] = item.meta['fit_info']['<scale>']
                self.filtered_table[index]['rotation'] = item.meta['fit_info']['rot'][1]
            else:
                self.filtered_table = None


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
            the WCSNAME keyword.  Common options: IMG, REL, SVM, and MVM.

        """
        if not self.selected_fit:
            log.error("No FIT selected for application.  Please run 'select_fit()' method.")
            raise ValueError
        # Call update_hdr_wcs()
        headerlet_dict = update_image_wcs_info(self.selected_fit,
                                               headerlet_filenames=headerlet_filenames,
                                               fit_label=fit_label)

        for table_index in range(0, len(self.filtered_table)):
            fname = self.filtered_table[table_index]['imageName']
            row = self.filtered_table[table_index]
            row['headerletFile'] = headerlet_dict[fname] if fname in headerlet_dict else "None"


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

        # Switch to turn on/off use of single-image CR detection/removal
        self.crclean = False

    def build_wht_image(self):
        if not self.num_wht:
            # Working with a calibrated exposure, no WHT extension
            # create a substitute WHT array from ERR and DQ
            # Build pseudo-wht array for detection purposes
            errarr = np.concatenate([self.imghdu[('ERR', i + 1)].data for i in range(self.num_sci)])
            wht_image = errarr.max() / errarr
            wht_image = np.nan_to_num(wht_image, nan=0.0, posinf=0.0, neginf=0.0)
            if self.dqmask is not None:
                wht_image[self.dqmask] = 0
        else:
            wht_image = self.imghdu['WHT'].data
        return wht_image

    def close(self):
        if self.imghdu is not None:
            self.imghdu.close()
            self.imghdu = None
        self.dqmask = None
        self.wht_image = None
        self._wht_image = None
        self.bkg_rms_mean = {}
        self.bkg = {}
        self.bkg_dao_rms = {}


    def build_kernel(self, fwhmpsf):
        """
        Parameters
        -----------
        fwhmpsf : float
            Default FWHM of PSF in units of arcseconds.
        """
        if self.bkg is None or self.bkg == {}:
            self.compute_background()

        threshold_rms = np.concatenate([rms for rms in self.threshold.values()])
        bkg = np.concatenate([background for background in self.bkg.values()])
        log.info("Looking for sample PSF in {}".format(self.rootname))
        log.debug("  based on RMS of {:9.4f}".format(threshold_rms.mean()))
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
            log.debug("  based on RMS of {:9.4f}".format(threshold_rms.mean()))
            fwhm = fwhmpsf / self.pscale

            k, self.kernel_fwhm = amutils.build_auto_kernel(self.data - bkg,
                                                            self.wht_image,
                                                            threshold=threshold_rms,
                                                            fwhm=fwhm)

        self.kernel, self.kernel_psf = k
        log.info("  Found PSF with FWHM = {:9.4f}".format(self.kernel_fwhm))

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
        log.debug("NSigma: {:9.4f}".format(nsigma))
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
                bkg_rms_mean = max(0.01, self.data.min())
                bkg_dao_rms = bkg_rms_mean
                threshold = (nsigma + 1) * bkg_rms_mean

            # *** FIX: Need to do something for bkg if bkg is None ***

            # Report other useful quantities
            log.debug("{} CHIP: {}".format(self.rootname, chip))
            log.debug("Mean background: {:9.4g}".format(bkg_rms_mean))
            log.debug("Mean threshold: {:9.4g}".format(np.mean(threshold)))
            log.debug("Mean RMS      : {:9.4g}".format(bkg_rms_mean))
            log.debug("")
            log.debug("{}".format("=" * 60))

            # Make copies using deepcopy so it will work with ndarray backgrounds
            # or scalar backgrounds depending on the situation.
            self.bkg[chip] = copy.deepcopy(bkg.background)
            self.bkg_dao_rms[chip] = copy.deepcopy(bkg_dao_rms)
            self.bkg_rms_mean[chip] = copy.deepcopy(bkg_rms_mean)
            self.threshold[chip] = copy.deepcopy(threshold)

            del bkg, default_threshold, threshold, bkg_dao_rms, bkg_rms_mean

    def build_dqmask(self, chip=None):
        # apply any DQ array, if available
        dqmask = None
        if chip:
            dqarr = self.imghdu[('DQ', chip)].data.astype(np.int16)
        else:
            dqarr = np.concatenate([self.imghdu[('DQ', i + 1)].data.astype(np.int16) for i in range(self.num_sci)])

        # "grow out" regions in DQ mask flagged as saturated by several
        # pixels in every direction to prevent the
        # source match algorithm from trying to match multiple sources
        # from one image to a single source in the
        # other or vice-versa.
        # Create temp DQ mask containing all pixels flagged with any value EXCEPT 256
        non_sat_mask = bitfield_to_boolean_mask(dqarr, ignore_flags=256 + 2048)

        # Create temp DQ mask containing saturated pixels ONLY
        sat_mask = bitfield_to_boolean_mask(dqarr, ignore_flags=~(256 + 2048))

        # Ignore sources where only a couple of pixels are flagged as saturated
        sat_mask = ndimage.binary_erosion(sat_mask, iterations=1)

        # Grow out saturated pixels by a few pixels in every direction
        grown_sat_mask = ndimage.binary_dilation(sat_mask, iterations=5)

        # combine the two temporary DQ masks into a single composite DQ mask.
        dqmask = np.bitwise_or(non_sat_mask, grown_sat_mask)

        # astropy's code returned the opposite bitmask from what was originally
        # defined by stsci.tools own bitmask code.
        if LooseVersion(stsci.tools.__version__) >= '4.0.0':
            dqmask = ~dqmask

        return dqmask


    def find_alignment_sources(self, output=True, dqname='DQ', crclean=False,
                               **alignment_pars):
        """Find sources in all chips for this exposure."""
        if crclean:
            self.imghdu = fits.open(self.imgname, mode='update')

        for chip in range(self.num_sci):
            chip += 1
            # find sources in image
            if output:
                outroot = '{}_sci{}_src'.format(self.rootname, chip)
            else:
                outroot = None

            dqmask = self.build_dqmask(chip=chip)
            sciarr = self.imghdu[("SCI", chip)].data.copy()
            #  TODO: replace detector_pars with dict from OO Config class
            # Turning off 'classify' since same CRs are being removed before segmentation now
            extract_pars = {'classify': False,  # alignment_pars['classify'],
                            'centering_mode': alignment_pars['centering_mode'],
                            'nlargest': alignment_pars['MAX_SOURCES_PER_CHIP'],
                            'deblend': alignment_pars['deblend']}

            with warnings.catch_warnings():
                warnings.simplefilter('ignore', NoDetectionsWarning)
                seg_tab, segmap, crmap = amutils.extract_sources(sciarr, dqmask=dqmask,
                                                                 outroot=outroot,
                                                                 kernel=self.kernel,
                                                                 segment_threshold=self.threshold[chip],
                                                                 dao_threshold=self.bkg_rms_mean[chip],
                                                                 fwhm=self.kernel_fwhm,
                                                                 **extract_pars)
            if crclean and crmap is not None:
                i = self.imgname.replace('.fits', '')
                if log.level < logutil.logging.INFO:
                    # apply crmap to input image
                    crfile = "{}_crmap_dq{}.fits".format(i, chip)
                    fits.PrimaryHDU(data=crmap).writeto(crfile, overwrite=True)
                log.debug("Updating DQ array for {} using single-image CR identification algorithm".format(i))
                self.imghdu[(dqname, chip)].data = np.bitwise_or(self.imghdu[(dqname, chip)].data, crmap)
                del crmap

            self.catalog_table[chip] = seg_tab

        if crclean and log.level < logutil.logging.INFO:
            self.imghdu.writeto(self.imgname.replace('.fits', '_crmap.fits'), overwrite=True)

        if self.imghdu is not None:
            self.imghdu.close()
            self.imghdu = None

# ----------------------------------------------------------------------------------------------------------------------


def match_relative_fit(imglist, reference_catalog, **fit_pars):
    """Perform cross-matching and final fit using relative matching algorithm

    The fitting performed by this algorithm first aligns all input images specified
    in the `imglist` to the FIRST image in that list.  When this fit is successful,
    it insures that the images are aligned RELATIVE to each other.  This set of co-aligned
    WCSs are then fit to the specified `reference_catalog` to improve the absolute
    astrometry for all input images (when successful).

    Parameters
    ----------
    imglist : list
        List of input image `~tweakwcs.tpwcs.FITSWCS` objects with metadata and source catalogs

    reference_catalog : Table
        Astropy Table of reference sources for this field

    fit_pars : dict
        Set of parameters and values to be used for the fit.  This should include
        `fitgeom` as well as any `tweakwcs.TPMatch
        <https://tweakwcs.readthedocs.io/en/latest/matchutils.html#tweakwcs.matchutils.TPMatch>`_
        parameter which the user feels needs to be adjusted to work best with the input data.

    Returns
    --------
    imglist : list
        List of input image `~tweakwcs.tpwcs.FITSWCS` objects with metadata and source catalogs

    """
    log.info("{} (match_relative_fit) Cross matching and fitting {}".format("-" * 20, "-" * 27))
    if 'fitgeom' in fit_pars:
        fitgeom = fit_pars['fitgeom']
        del fit_pars['fitgeom']
    else:
        fitgeom = 'rscale'

    rel_fitgeom = 'rscale'

    common_pars = fit_pars['pars']
    del fit_pars['pars']

    nclip = 1 if fitgeom == 'rscale' else 0  # Only perform sigma-clipping for 'rscale'

    # 0: Specify matching algorithm to use
    match = tweakwcs.TPMatch(**fit_pars)
    # match = tweakwcs.TPMatch(searchrad=250, separation=0.1,
    #                          tolerance=100, use2dhist=False)

    # Align images and correct WCS
    # NOTE: this invocation does not use an astrometric catalog. This call
    #       allows all the input images to be aligned in
    #       a relative way using the first input image as the reference.
    # 1: Perform relative alignment
    # Setting 'minobj' to None allows 'tweakwcs' to use its own built-in
    # limits.
    match_relcat = tweakwcs.align_wcs(imglist, None,
                                      match=match,
                                      minobj=common_pars['minobj'][rel_fitgeom],
                                      expand_refcat=True,
                                      fitgeom=rel_fitgeom,
                                      nclip=1)
    # Implement a consistency check even before trying absolute alignment
    # If relative alignment in question, no use in aligning to GAIA
    if not check_consistency(imglist):
            return imglist

    log.info("Relative alignment found: ")
    for i in imglist:
        info = i.meta['fit_info']
        if 'shift' not in info:
            off = (0., 0.)
            rot = 0.
            scale = -1.
            nmatches = 0
        else:
            off = info['shift']
            rot = info['<rot>']
            scale = info['<scale>']
            nmatches = info['nmatches']
        msg = "Image {} --".format(i.meta['name'])
        msg += "\n    SHIFT:({:9.4f},{:9.4f})  NMATCHES: {} ".format(off[0], off[1], nmatches)
        msg += "\n    ROT:{:9.4f}  SCALE:{:9.4f}".format(rot, scale)
        msg += "\n Using fitgeom = '{}'".format(rel_fitgeom)
        log.info(msg)

    # This logic enables performing only relative fitting and skipping fitting to GAIA
    if reference_catalog is not None:
        # Set all the group_id values to be the same so the various images/chips will be aligned to the astrometric
        # reference catalog as an ensemble.
        # astrometric reference catalog as an ensemble. BEWARE: If additional iterations of solutions are to be
        # done, the group_id values need to be restored.
        for image in imglist:
            image.meta["group_id"] = 1234567
        # 2: Perform absolute alignment
        matched_cat = tweakwcs.align_wcs(imglist, reference_catalog,
                                         match=match,
                                         minobj=common_pars['minobj'][fitgeom],
                                         fitgeom=fitgeom,
                                         nclip=nclip)
        # Insure the expanded reference catalog has all the information needed
        # to complete processing.
        # TODO: Work out how to get the 'mag' column from input source catalog
        #       into this extended reference catalog...
        # reference_catalog = match_relcat
        # reference_catalog['mag'] = np.array([-999.9] * len(reference_catalog),
        #                                    np.float32)
        # reference_catalog.meta['catalog'] = 'relative'

        # 3: Interpret RMS values from tweakwcs
    interpret_fit_rms(imglist, reference_catalog)

    del match_relcat
    return imglist

# ----------------------------------------------------------------------------------------------------------


def match_default_fit(imglist, reference_catalog, **fit_pars):
    """Perform cross-matching and final fit using default tolerance matching

    This function performs the specified type of fit ('general', 'rscale', ...) directly
    between all input images and the reference catalog.  If successful, each input image will
    have its absolute WCS aligned to the reference catalog.

    Parameters
    ----------
    imglist : list
        List of input image `~tweakwcs.tpwcs.FITSWCS` objects with metadata and source catalogs

    reference_catalog : Table
        Astropy Table of reference sources for this field

    fit_pars : dict
        Set of parameters and values to be used for the fit.  This should include
        `fitgeom` as well as any `tweakwcs.TPMatch
        <https://tweakwcs.readthedocs.io/en/latest/matchutils.html#tweakwcs.matchutils.TPMatch>`_
        parameter which the user feels needs to be adjusted to work best with the input data.

    Returns
    --------
    imglist : list
        List of input image `~tweakwcs.tpwcs.FITSWCS` objects with metadata and source catalogs

    """
    if 'fitgeom' in fit_pars:
        fitgeom = fit_pars['fitgeom']
        del fit_pars['fitgeom']
    else:
        fitgeom = 'rscale'

    common_pars = fit_pars['pars']
    del fit_pars['pars']

    nclip = 1 if fitgeom == 'rscale' else 0  # Only perform sigma-clipping for 'rscale'

    log.info("{} (match_default_fit) Cross matching and fitting "
             "{}".format("-" * 20, "-" * 27))
    # Specify matching algorithm to use
    match = tweakwcs.TPMatch(**fit_pars)

    # Align images and correct WCS
    matched_cat = tweakwcs.align_wcs(imglist, reference_catalog,
                                     match=match,
                                     minobj=common_pars['minobj'][fitgeom],
                                     expand_refcat=False,
                                     fitgeom=fitgeom,
                                     nclip=nclip)

    # Interpret RMS values from tweakwcs
    interpret_fit_rms(imglist, reference_catalog)

    del matched_cat

    return imglist


# ----------------------------------------------------------------------------------------------------------------------


def match_2dhist_fit(imglist, reference_catalog, **fit_pars):
    """Perform cross-matching and final fit using 2dHistogram matching

    This function performs cross-matching of each separate input image to the
    sources in the reference catalog by looking for common integer offsets between all
    sources in the input list and all sources in the reference catalog.  This
    offset is then used as the starting point for the final fit to the reference
    catalog to align each input image SEPARATELY to the reference catalog.

    Parameters
    ----------
    imglist : list
        List of input image `~tweakwcs.tpwcs.FITSWCS` objects with metadata and source catalogs

    reference_catalog : Table
        Astropy Table of reference sources for this field

    fit_pars : dict
        Set of parameters and values to be used for the fit.  This should include
        `fitgeom` as well as any `tweakwcs.TPMatch
        <https://tweakwcs.readthedocs.io/en/latest/matchutils.html#tweakwcs.matchutils.TPMatch>`_
        parameter which the user feels needs to be adjusted to work best with the input data.

    Returns
    --------
    imglist : list
        List of input image `~tweakwcs.tpwcs.FITSWCS` objects with metadata and source catalogs

    """
    if 'fitgeom' in fit_pars:
        fitgeom = fit_pars['fitgeom']
        del fit_pars['fitgeom']
    else:
        fitgeom = 'rscale'

    common_pars = fit_pars['pars']
    del fit_pars['pars']

    nclip = 1 if fitgeom == 'rscale' else 0  # Only perform sigma-clipping for 'rscale'

    log.info("{} (match_2dhist_fit) Cross matching and fitting "
             "{}".format("-" * 20, "-" * 28))
    # Specify matching algorithm to use
    match = tweakwcs.TPMatch(**fit_pars)
    # Align images and correct WCS
    matched_cat = tweakwcs.align_wcs(imglist, reference_catalog,
                                     match=match,
                                     minobj=common_pars['minobj'][fitgeom],
                                     expand_refcat=False,
                                     fitgeom=fitgeom,
                                     nclip=nclip)

    # Interpret RMS values from tweakwcs
    interpret_fit_rms(imglist, reference_catalog)

    del matched_cat

    return imglist

# ----------------------------------------------------------------------------------------------------------
def check_consistency(imglist, rot_tolerance=0.1, shift_tolerance=1.0):
    """Check to see whether relative alignment solutions are realistically consistent

        Should any image end up with a relative alignment fit where rot is more than
        0.1 degree or 1 pixel off from the other images fits, this check will
        mark this fit as 'FAILED' in the '.meta["fit_info"]["status"]' field for
        all images.
    """
    is_consistent = True

    # set up arrays for relative alignment rotation and shifts
    rots = np.zeros(len(imglist)+1, np.float64)
    nmatches = np.zeros(len(imglist), np.int16)

    for i, img in enumerate(imglist):
        finfo = img.meta['fit_info']
        status = finfo['status']
        if not status.startswith("FAILED"):
            if 'rot' in finfo:
                rots[i] = finfo['proper_rot']
                nmatches[i] = finfo['nmatches']
        else:
            # 'status' == 'FAILED': Set default as negative value
            nmatches[i] = -1
            is_consistent = False

    # We should only need to check for consistency when less than 5
    # matches were used for the fit, leading to a higher potential for
    # a singular solution or mis-identified cross-matches.
    if nmatches.min() > 4:
        return is_consistent

    # compute deltas to look for outliers
    delta_rots = rots[1:] - rots[:-1]
    if delta_rots.max() >= rot_tolerance:
        for i, img in enumerate(imglist):
            msg = 'Relative consistency check failed: {}'.format(rots[i])
            img.meta['fit_info']['status'] = 'FAILED'
            img.meta['fit_info']['process_msg'] = msg

        log.info('Relative fit solution is NOT consistent!')
        fitgeom = finfo['fitgeom'] if 'fitgeom' in finfo else 'Unknown'
        log.debug('DELTAS for "{}" fit:'.format(fitgeom))
        log.debug('  max rot={:.4f}\n '.format(delta_rots.max()))
        is_consistent = False

    else:
        log.info('Relative fit solution is consistent')

    return is_consistent


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
            tinfo = item.meta['fit_info']
            # When status = FAILED (fit failed) or REFERENCE (relative alignment done with first image
            # as the reference), skip to the beginning of the loop as there is no 'fit_info'.
            if tinfo['status'] != 'SUCCESS' or (tinfo['status'] == 'SUCCESSS' and \
                'fitmask' not in tinfo):
                continue
            # Make sure to store data for any particular group_id only once.
            if item.meta['group_id'] == group_id and \
               group_id not in group_dict:
                group_dict[group_id] = {'ref_idx': None, 'FIT_RMS': None,
                                        'input_mag': None, 'ref_mag': None, 'input_idx': None}

                # Perform checks to insure that fit was successful
                # Namely, rot < 0.1 degree and scale < 1%.
                if abs(tinfo['<scale>'] - 1) > 0.01 or abs(tinfo['<rot>']) > 0.1 or \
                    tinfo['skew'] > 0.01:
                    # fit is bad
                    tinfo['status'] = 'FAILED: singularity in fit'

                log.debug("fit_info: {}".format(item.meta['fit_info']))

                ref_idx = tinfo['matched_ref_idx']
                fitmask = tinfo['fitmask']
                group_dict[group_id]['ref_idx'] = ref_idx
                input_RA = tinfo['fit_RA']
                input_DEC = tinfo['fit_DEC']
                ref_RA = reference_catalog[ref_idx]['RA'][fitmask]
                ref_DEC = reference_catalog[ref_idx]['DEC'][fitmask]

                img_coords = SkyCoord(input_RA, input_DEC,
                                      unit='deg', frame='icrs')
                ref_coords = SkyCoord(ref_RA, ref_DEC, unit='deg', frame='icrs')
                dra, ddec = img_coords.spherical_offsets_to(ref_coords)
                ra_rms = np.std(dra.to(u.mas))
                dec_rms = np.std(ddec.to(u.mas))
                fit_rms = np.std(Angle(img_coords.separation(ref_coords), unit=u.mas)).value

                # Get mag from reference catalog, or input image catalog as needed
                group_dict[group_id]['ref_mag'] = reference_catalog[ref_idx]['mag'][fitmask]


                group_dict[group_id]['FIT_RMS'] = fit_rms
                group_dict[group_id]['RMS_RA'] = ra_rms
                group_dict[group_id]['RMS_DEC'] = dec_rms

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
        fitmask = item.meta['fit_info'].get('fitmask', None)
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

        fit_label : string, optional
            Short label to use for (part of) the name of the WCS being updated in the image header.
            This label will be appended to the end of the full `WCSNAME` value generated from the
            `IDCTAB` name and type of fit.


        Returns
        -------
        out_headerlet_list : dictionary
            a dictionary of the headerlet files created by this subroutine, keyed by flt/flc fits filename.
        """
    out_headerlet_dict = {}
    for item in tweakwcs_output:
        image_name = item.meta['filename']
        chipnum = item.meta['chip']
        hdulist = fits.open(image_name, mode='update')
        # start by insuring that the current version of STWCS and Astropy are recorded
        # in the header, not the versions used to create the previous WCS
        # This logic comes from STWCS in order to be compatible with STWCS in
        # terms of where these keywords should be found in the PRIMARY header.
        upwcsver = stwcs.__version__
        pywcsver = astropy.__version__

        log.info('Updating PRIMARY header with:')
        log.info('    UPWCSVER = {}'.format(upwcsver))
        log.info('    PYWCSVER = {}'.format(pywcsver))
        if 'HISTORY' in hdulist[0].header:
            after_kw = None
            before_kw = 'HISTORY'
        elif 'ASN_MTYP' in hdulist[0].header:
            after_kw = 'ASN_MTYP'
            before_kw = None
        else:
            after_kw = hdulist[0].header.cards[-1][0]
            before_kw = None

        hdulist[0].header.set('UPWCSVER', value=upwcsver,
                        comment="Version of STWCS used to update the WCS",
                        after=after_kw, before=before_kw)
        hdulist[0].header.set('PYWCSVER', value=pywcsver,
                        comment="Version of Astropy used to update the WCS",
                        after='UPWCSVER')

        if chipnum == 1:
            chipctr = 1
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
        sci_extn = sci_ext_dict["{}".format(item.meta['chip'])]
        hdr_name = "{}_{}-hlet.fits".format(image_name.rstrip(".fits"), wcs_name)
        updatehdr.update_wcs(hdulist, sci_extn, item.wcs, wcsname=wcs_name, reusename=True)
        info = item.meta['fit_info']
        hdulist[sci_extn].header['RMS_RA'] = info['RMS_RA'].value if info['RMS_RA'] is not None else -1.0
        hdulist[sci_extn].header['RMS_DEC'] = info['RMS_DEC'].value if info['RMS_DEC'] is not None else -1.0
        hdulist[sci_extn].header['CRDER1'] = info['RMS_RA'].value if info['RMS_RA'] is not None else -1.0
        hdulist[sci_extn].header['CRDER2'] = info['RMS_DEC'].value if info['RMS_DEC'] is not None else -1.0
        hdulist[sci_extn].header['NMATCHES'] = len(info['ref_mag']) if info['ref_mag'] is not None else -1.0
        if 'HDRNAME' in hdulist[sci_extn].header:
            del hdulist[sci_extn].header['HDRNAME']
        hdulist[sci_extn].header['HDRNAME'] = hdr_name
        hdulist.flush()
        hdulist.close()

        # Create headerlet
        out_headerlet = headerlet.create_headerlet(image_name, hdrname=hdr_name, wcsname=wcs_name,
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
        out_headerlet.writeto(headerlet_filename, overwrite=True)
        log.info("Wrote headerlet file {}.\n\n".format(headerlet_filename))
        out_headerlet_dict[image_name] = headerlet_filename

        # Attach headerlet as HDRLET extension
        if headerlet.verify_hdrname_is_unique(hdulist, hdr_name):
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
    info = tweakwcs_item.meta['fit_info']
    if 'nmatches' in info:
        rms_ra = info['RMS_RA'].value if info['RMS_RA'] is not None else -1.0
        rms_dec = info['RMS_DEC'].value if info['RMS_RA'] is not None else -1.0
        fit_rms = info['FIT_RMS']
        nmatch = info['nmatches']
        catalog = info['catalog']
        fit_method = tweakwcs_item.meta['fit method']

        x_shift = (info['shift'])[0]
        y_shift = (info['shift'])[1]
        rot = info['rot'][1]  # report rotation of Y axis only
        scale = info['<scale>']
        skew = info['skew']

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

    func = getattr(background, name)  # raises AttributeError if not found
    return func
