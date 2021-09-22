""" Definition of Super and Subclasses for the mosaic output image_list

    Classes which define the total ("white light" image), filter, and exposure
    drizzle products.

    These products represent different levels of processing with the levels noted in the
    'HAPLEVEL' keyword.  The 'HAPLEVEL' values are:

      * 1 : calibrated (FLT/FLC) input images and exposure level drizzle products with improved astrometry
      * 2 : filter and total products combined using the improved astrometry, consistent pixel scale, and oriented to North.
      * 3 : (future) multi-visit mosaics aligned to common tangent plane
"""
import copy
import logging
import sys
import os
import traceback
import shutil

from astropy.io import fits
from collections import OrderedDict
from stsci.tools import logutil
from stwcs import updatewcs
import numpy as np

from .. import align
from .. import astrodrizzle
from .. import wcs_functions
from . import align_utils
from . import astrometric_utils as amutils
from . import cell_utils

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)

# Define keywords to be added to SVM products to describe the overlap of
# exposures onto the output WCS footprint as determined from the SkyFootprint masks
MASK_KWS = {"NPIXFRAC": [None, "Fraction of pixels with data"],
            "MEANEXPT": [None, "Mean exposure time per pixel with data"],
            "MEDEXPT": [None, "Median exposure time per pixel with data"],
            "MEANNEXP": [None, "Mean number of exposures per pixel with data"],
            "MEDNEXP": [None, "Median number of exposures per pixel with data"],
            }


class HAPProduct:
    """ HAPProduct is the base class for the various products generated during the
        astrometry update and mosaicing of image data.
    """
    def __init__(self, prop_id, obset_id, instrument, detector, filename, filetype, log_level):
        # set logging level to user-specified level
        log.setLevel(log_level)
        self.log_level = log_level

        # Make sure the proposal ID is a 5-character string
        self.prop_id = prop_id.zfill(5)
        self.obset_id = obset_id
        self.instrument = instrument
        self.detector = detector
        self.filetype = filetype
        self.rules_file = None

        self.basename = "hst_" + "_".join(map(str, [prop_id, obset_id, instrument, detector])) + "_"

        # exposure_name is the ipppssoo or a portion thereof
        self.exposure_name = filename[0:8]

        # TO DO: update this variable
        self.mjdutc = None

        # HAPConfig objects are created after these Product objects have been instantiated so
        # this attribute is updated in the hapsequncer.py module (run_hla_processing()).
        self.configobj_pars = None

        # Define HAPLEVEL value for this product
        self.haplevel = 1

        # Initialize attributes for use in generating the output products
        self.meta_wcs = None
        self.mask = None
        self.mask_kws = MASK_KWS.copy()

    # def print_info(self):
        # """ Generic print at this time to indicate the information used in the
        #     construction of the object for debug purposes.
        # """
        # print("Object information: {}".format(self.info))

    def generate_footprint_mask(self):
        """ Create a footprint mask for a set of exposure images

            Create a mask which is True/1/on for the illuminated portion of the image, and
            False/0/off for the remainder of the image.
        """
        footprint = cell_utils.SkyFootprint(self.meta_wcs)
        exposure_names = [element.full_filename for element in self.edp_list]
        footprint.build(exposure_names, scale=True)

        # This mask actually represents the number of chips per pixel, not True/False.
        # To have the True/False mask it should be self.mask = footprint.footprint.
        # Do not fix this until it can be verified that a change will not have repercussions.
        self.mask = footprint.total_mask

        # Compute footprint-based SVM-specific keywords for product image header
        good_pixels = self.mask > 0
        self.mask_kws['NPIXFRAC'][0] = good_pixels.sum() / self.mask.size
        self.mask_kws['MEANEXPT'][0] = np.mean(footprint.scaled_mask[good_pixels])
        self.mask_kws['MEDEXPT'][0] = np.median(footprint.scaled_mask[good_pixels])
        self.mask_kws['MEANNEXP'][0] = np.mean(self.mask[good_pixels])
        self.mask_kws['MEDNEXP'][0] = np.median(self.mask[good_pixels])

    def generate_metawcs(self):
        """ A method to build a unique WCS for each TotalProduct product which is
            generated based upon the merging of all the ExposureProducts
            which comprise the specific TotalProduct.  This
            is done on a per TotalProduct basis as the Single
            Visit Mosaics need to be represented in their native scale.

        """
        exposure_filenames = [element.full_filename for element in self.edp_list]
        log.debug("\n\nRun make_mosaic_wcs to create a common WCS.")
        log.debug("The following images will be used: ")
        log.debug("{}\n".format(exposure_filenames))

        # Set the rotation to 0.0 to force North as up
        if exposure_filenames:
            meta_wcs = wcs_functions.make_mosaic_wcs(exposure_filenames, rot=0.0)

        # Used in generation of SkyFootprints
        self.meta_wcs = meta_wcs

        return meta_wcs

    def align_to_gaia(self, catalog_list=[], output=True,
                      fit_label='SVM', align_table=None, fitgeom=''):
        """Extract the flt/flc filenames from the exposure product list, as
           well as the corresponding headerlet filenames to use legacy alignment
           routine.
        """

        # Make a deep copy of the alignment parameters for safe keeping
        alignment_pars = self.configobj_pars.get_pars('alignment')
        orig_alignment_pars = copy.deepcopy(alignment_pars)

        exposure_filenames = []
        headerlet_filenames = {}
        align_table = None
        crclean = []

        # If no catalog list has been provided, use the list defined in the configuration file
        if not catalog_list:
            mosaic_catalog_list = alignment_pars['run_align']['mosaic_catalog_list']
        else:
            mosaic_catalog_list = catalog_list
        num_cat = len(mosaic_catalog_list)

        # Fitting methods
        mosaic_method_list = alignment_pars['run_align']['mosaic_fit_list']
        if len(self.edp_list) == 1:
            # Remove 'match_relative_fit' (first entry) since there is only 1 image to align
            del mosaic_method_list[0]
        num_method = len(mosaic_method_list)

        # Fit geometry methods
        # mosaic_fitgeom_list = list(alignment_pars['run_align']['mosaic_fitgeom_list'][0].items())
        mosaic_fitgeom_list = list(alignment_pars['run_align']['mosaic_fitgeom_list'].items())

        for edp in self.edp_list:
            exposure_filenames.append(edp.full_filename)
            headerlet_filenames[edp.full_filename] = edp.headerlet_filename
            crclean.append(edp.crclean)

        try:
            # Proceed only if there is data to process
            if self.edp_list:

                # If necessary, generate the alignment table only once
                if align_table is None:
                    align_table = align_utils.AlignmentTable(exposure_filenames,
                                                             log_level=self.log_level,
                                                             **alignment_pars)
                    align_table.find_alignment_sources(output=output, crclean=crclean)

                is_good_fit = False
                # determine minimum number of sources available for alignment
                source_nums = []
                for img_cats in align_table.extracted_sources.values():
                    sources = 0
                    for chip in img_cats:
                        sources += len(img_cats[chip])
                    source_nums.append(sources)
                min_sources = min(source_nums)

                # Loop for available catlogs, the first successful fit for a
                # (catalog, fitting method, and fit geometry) is satisfactory to break out of the looping.
                for index_cat, catalog_item in enumerate(mosaic_catalog_list):

                    more_catalogs = True
                    if (index_cat + 1) == num_cat:
                        more_catalogs = False 

                    # Override the default self.refname as it really needs to be
                    # catalog-specific to be useful
                    self.refname = self.product_basename + "_" + catalog_item + "_ref_cat.ecsv"

                    log.info("Starting alignment to absolute astrometric reference frame '{}'.".format(catalog_item))
                    log.debug('Creating reference catalog {}'.format(self.refname))
                    ref_catalog = amutils.create_astrometric_catalog(align_table.process_list,
                                                                     catalog=catalog_item,
                                                                     output=self.refname,
                                                                     gaia_only=False,
                                                                     full_catalog=True)

                    # Add a weight column which is based upon proper motion measurements
                    if 'converted' in ref_catalog.meta and ref_catalog.meta['converted'] \
                            and ref_catalog['RA_error'][0] != np.nan:
                        ref_weight = np.sqrt(ref_catalog['RA_error'] ** 2 + ref_catalog['DEC_error'] ** 2)
                        ref_weight = np.nan_to_num(ref_weight, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
                    else:
                        ref_weight = np.ones_like(ref_catalog['RA'])
                    ref_catalog.add_column(ref_weight, name='weight')

                    log.debug("Abbreviated reference catalog displayed below\n{}".format(ref_catalog))
                    align_table.reference_catalogs[self.refname] = ref_catalog

                    # Will need to satisfy the minimum criterion of found sources for the proposed
                    # fit geometry, or the fit geometry must be downgraded
                    num_ref_sources = len(ref_catalog)
                    log.info("Number of sources for reference catalog '{}' is {}.".format(catalog_item, num_ref_sources))

                    # Loop for specified fit methods
                    for index_method, full_method_name in enumerate(mosaic_method_list):
                        method_name = full_method_name.split("_")[1]

                        # If fitgeom is not provided via a parameter, get it from the configuration.
                        # This is the requested (from parameter) or default (from configuration)
                        # fit geometry, but it may be downgraded depending upon the number of
                        # sources available.
                        if not fitgeom.strip():
                            mosaic_fitgeom = alignment_pars[full_method_name]['fitgeom']
                        else:
                            mosaic_fitgeom = fitgeom

                        # Get the index of the fitgeom in the list of all possible fitgeoms
                        for index, mode in enumerate(mosaic_fitgeom_list):
                            if mode[0].lower().strip() == mosaic_fitgeom:
                                mosaic_fitgeom_index = index
                                break

                        # If there are not enough references sources for the specified fitgeom,
                        # downgrade the fitgeom until a valid one is found.  Also, if the fit done
                        # with the fitgeom was unsatisfactory, downgrade if possible and try again
                        while num_ref_sources < mosaic_fitgeom_list[mosaic_fitgeom_index][1]:
                            log.warning("Not enough reference sources for alignment using catalog '{}' with fit method '{}' and fit geometry '{}'.".format(catalog_item, method_name, mosaic_fitgeom_list[mosaic_fitgeom_index][0]))
                            mosaic_fitgeom_index -= 1
                            if mosaic_fitgeom_index < 0:
                                log.warning("No further fit geometries to try. Proceeding to try another fit method.")
                                break

                        while mosaic_fitgeom_index > -1:
                            # In case of a downgrade, make sure the fit geometry "name" is based upon the index
                            mosaic_fitgeom = mosaic_fitgeom_list[mosaic_fitgeom_index][0]
                            log.info("Proceeding to use the '{}' fit geometry.".format(mosaic_fitgeom))

                            try:
                                # Always perform the configuration to ensure there is no residual cruft
                                # from any possible previous fit
                                align_table.configure_fit()

                                log.info("Trying '{}' wth method '{}' using fit geometry '{}'.".
                                         format(catalog_item, method_name, mosaic_fitgeom))
                                align_table.imglist = align_table.perform_fit(method_name, catalog_item, ref_catalog,
                                                                              fitgeom=mosaic_fitgeom)

                                # Define comparison for min cross-matches based on fitgeom used
                                alignment_pars['determine_fit_quality']['min_xmatches'] = \
                                    alignment_pars['run_align']['mosaic_fitgeom_list'][mosaic_fitgeom]

                                # turn off consistency check for SVM and MVM processing since filter-to-filter
                                # or visit-to-visit offsets could be large relative to measurement RMS.
                                alignment_pars['determine_fit_quality']['consistency_check'] = False

                                # Evaluate the quality of the fit
                                is_good_fit, _, _, _, _, _ = align.determine_fit_quality_mvm_interface(align_table.imglist,
                                                                                                       align_table.filtered_table,
                                                                                                       more_catalogs,
                                                                                                       num_cat,
                                                                                                       alignment_pars,
                                                                                                       print_fit_parameters=True,
                                                                                                       loglevel=self.log_level)

                                # Ensure the original parameters stay intact for the iterations
                                # as the perform_fit() modifies the fitgeom
                                alignment_pars = copy.deepcopy(orig_alignment_pars)

                                # The fit was good...
                                if is_good_fit:
                                    log.info("Fit done for catalog '{}' with method '{}' using fit geometry '{}' was satisfactory.".
                                             format(catalog_item, method_name, mosaic_fitgeom))

                                    # Success - Break out of "fitgeom" while loop
                                    break

                                # The fit was not good
                                else:
                                    log.info("Fit done for catalog '{}' with method '{}' using fit geometry '{}' was NOT satisfactory.".
                                             format(catalog_item, method_name, mosaic_fitgeom))

                            except Exception:
                                log.info("Problem with fit done for catalog '{}' with method '{}' using fit geometry '{}'.".
                                         format(catalog_item, method_name, mosaic_fitgeom))
                                traceback.print_exc()


                            # Try again with a different fit geometry algorithm
                            mosaic_fitgeom_index -= 1
                            log.info("Resetting mosaic index for next fit method")
                            # Only if there are too few sources to perform an alignment with the current
                            # method should the next method be attempted.  Otherwise, it is assumed that the
                            # image sources are either cosmic-rays or features that do not match across
                            # the filters (e.g., F300X vs F475W from 'id7r03') and can't be aligned.  Trying with a
                            # less strict method (rshift vs rscale) will only increase the probability of
                            # matching incorrectly and introducing an error to the relative alignment.
                            if min_sources > mosaic_fitgeom_list[mosaic_fitgeom_index][1]:
                                mosaic_fitgeom_index = -1
                                log.warning("Too many (bad?) sources to try any more fitting with this catalog.")

                        # Break out of the "fit method" for loop
                        if is_good_fit:
                            break
                        else:
                            if index_method + 1 < num_method:
                                log.warning("Proceeding to try the next fit method.")

                    # Break out of the "catalog" for loop
                    if is_good_fit:
                        break
                    else:
                        # log.warning("Not enough reference sources for absolute alignment using catalog {}.".format(catalog_item))
                        if index_cat + 1 < num_cat:
                            log.info("Proceeding to try the next catalog.")
                        else:
                            log.warning("Not enough reference sources for absolute alignment in any available catalog.")
                            break

                # Got a good fit...
                if is_good_fit:
                    align_table.select_fit(catalog_item, method_name)
                    align_table.apply_fit(headerlet_filenames=headerlet_filenames,
                                          fit_label=fit_label)
                else:
                    log.warning("No satisfactory fit found for any catalog.")
                    raise ValueError

        except Exception:
            # Report a problem with the alignment
            if fit_label.upper().strip() == 'SVM':
                log.warning("EXCEPTION encountered in align_to_gaia for the FilteredProduct.\n")
            else:
                log.warning("EXCEPTION encountered in align_to_gaia for the SkyCellProduct.\n")
            log.warning("No correction to absolute astrometric frame applied.\n")
            log.warning("Proceeding with previous best solution.\n")

            # Only write out the traceback if in "debug" mode since not being able to
            # align the data to an absolute astrometric frame is not actually a failure.
            log.debug(traceback.format_exc())
            align_table = None

            # If the align_table is None, it is necessary to clean-up reference catalogs
            # created for alignment of each filter product here.
            if self.refname and os.path.exists(self.refname):
                os.remove(self.refname)

        # Return a table which contains data regarding the alignment, as well as the
        # list of the flt/flc exposures which were part of the alignment process
        # TODO: This does not account for individual exposures which might have been
        # excluded from alignment.
        return align_table, exposure_filenames


class TotalProduct(HAPProduct):
    """ A Total Detection Product is a 'white' light mosaic comprised of
        images acquired with one instrument, one detector, all filters, and all
        exposure times.  The exposure time characteristic may change - TBD.

        The "tdp" is short hand for TotalProduct.
    """
    def __init__(self, prop_id, obset_id, instrument, detector, filename, filetype, log_level):
        super().__init__(prop_id, obset_id, instrument, detector, filename, filetype, log_level)
        self.info = '_'.join([prop_id, obset_id, instrument, detector, filename, filetype])
        self.exposure_name = filename[0:6]

        self.product_basename = self.basename + "total_" + self.exposure_name
        self.trl_logname = self.product_basename + "_trl.log"
        self.trl_filename = self.product_basename + "_trl.txt"
        self.point_cat_filename = self.product_basename + "_point-cat.ecsv"
        self.segment_cat_filename = self.product_basename + "_segment-cat.ecsv"
        self.drizzle_filename = self.product_basename + "_" + self.filetype + ".fits"
        self.ref_drizzle_filename = self.product_basename + "_ref_" + self.filetype + ".fits"

        # Generate the name for the manifest file which is for the entire visit.  It is fine
        # to create it as an attribute of a TotalProduct as it is independent of
        # the detector in use.
        # instrument_programID_obsetID_manifest.txt (e.g.,wfc3_b46_06_manifest.txt)
        self.manifest_name = '_'.join([instrument, filename[1:4], obset_id, "manifest.txt"])

        # Define HAPLEVEL value for this product
        self.haplevel = 2

        # These attributes will be populated during processing
        self.edp_list = []
        self.fdp_list = []
        self.regions_dict = {}
        self.grism_edp_list = []
        self.bkg_used = ""

        log.debug("Total detection object {}/{} created.".format(self.instrument, self.detector))

    def update_drizpars(self):
        """ Update ALL products with final name of trailer file """
        # Update this total product
        fits.setval(self.drizzle_filename, 'DRIZPARS', value=self.trl_filename)

        # Update all filter drizzle products that were actually written out
        for fdp in self.fdp_list:
            if os.path.exists(fdp.drizzle_filename):
                fits.setval(fdp.drizzle_filename, 'DRIZPARS', value=fdp.trl_filename)

        # Update all exposure drizzle products that were actually written out
        for edp in self.edp_list:
            if os.path.exists(edp.drizzle_filename):
                fits.setval(edp.drizzle_filename, 'DRIZPARS', value=edp.trl_filename)

    def find_member(self, name):
        """ Return member instance with filename 'name' """
        desired_member = None
        for member in [self] + self.edp_list + self.fdp_list:
            if name == member.drizzle_filename:
                desired_member = member
                break
        return desired_member

    def add_member(self, edp):
        """ Add an ExposureProduct object to the list - composition.
        """
        self.edp_list.append(edp)

    def add_grism_member(self, edp):
        """ Add an GrismExposureProduct object to the list - composition.

            This routine adds Grism or Prism exposures to the exposure list.
        """
        self.grism_edp_list.append(edp)

    def add_product(self, fdp):
        """ Add a FilterProduct object to the list - composition.
        """
        self.fdp_list.append(fdp)

    def wcs_drizzle_product(self, meta_wcs):
        """
            Create the drizzle-combined total image using the meta_wcs as the reference output

            .. note:: Cosmic-ray identification is NOT performed when creating the total detection image.
        """
        # Retrieve the configuration parameters for astrodrizzle
        drizzle_pars = self.configobj_pars.get_pars("astrodrizzle")
        # ...and set parameters which are computed on-the-fly
        drizzle_pars["final_refimage"] = meta_wcs
        drizzle_pars["runfile"] = self.trl_logname

        # Setting "preserve" to false so the OrIg_files directory is deleted as the purpose
        # of this directory is now obsolete.
        drizzle_pars["preserve"] = False
        drizzle_pars['rules_file'] = self.rules_file
        drizzle_pars['resetbits'] = "0"

        log.debug("The 'final_refimage' ({}) and 'runfile' ({}) configuration variables "
                  "have been updated for the drizzle step of the total drizzle product."
                  .format(meta_wcs, self.trl_logname))

        edp_filenames = [element.full_filename for element in self.edp_list]
        astrodrizzle.AstroDrizzle(input=edp_filenames,
                                  output=self.drizzle_filename,
                                  **drizzle_pars)

        # Update product with SVM-specific keywords based on the footprint
        with fits.open(self.drizzle_filename, mode='update') as hdu:
            for kw in self.mask_kws:
                hdu[("SCI", 1)].header[kw] = tuple(self.mask_kws[kw])

        # Rename Astrodrizzle log file as a trailer file
        log.debug("Total combined image {} composed of: {}".format(self.drizzle_filename, edp_filenames))
        try:
            shutil.move(self.trl_logname, self.trl_filename)
        except PermissionError:
            pass


class FilterProduct(HAPProduct):
    """ A Filter Detection Product is a mosaic comprised of images acquired
        during a single visit with one instrument, one detector, a single filter,
        and all exposure times.  The exposure time characteristic may change - TBD.

        The "fdp" is short hand for FilterProduct.
    """
    def __init__(self, prop_id, obset_id, instrument, detector, filename, filters, filetype, log_level):
        super().__init__(prop_id, obset_id, instrument, detector, filename, filetype, log_level)

        self.info = '_'.join([prop_id, obset_id, instrument, detector, filename, filters, filetype])
        if filename[0:7].lower() != "metawcs":
            self.exposure_name = filename[0:6]
            self.product_basename = self.basename + "_".join(map(str, [filters, self.exposure_name]))
        else:
            self.exposure_name = "metawcs"
            self.product_basename = self.basename + "_".join(map(str, [filetype, self.exposure_name, filters]))

        self.filters = filters

        # Trailer names .txt or .log
        self.trl_logname = self.product_basename + "_trl.log"
        self.trl_filename = self.product_basename + "_trl.txt"
        self.point_cat_filename = self.product_basename + "_point-cat.ecsv"
        self.segment_cat_filename = self.product_basename + "_segment-cat.ecsv"
        self.drizzle_filename = self.product_basename + "_" + self.filetype + ".fits"

        # This attribute is reset in align_to_gaia to distinguish different catalogs
        # which may be created
        self.refname = self.product_basename + "_ref_cat.ecsv"

        # Define HAPLEVEL value for this product
        self.haplevel = 2

        # These attributes will be populated during processing
        self.edp_list = []
        self.regions_dict = {}

        log.info("Filter object {}/{}/{} created.".format(self.instrument, self.detector, self.filters))

    def find_member(self, name):
        """ Return member instance with filename 'name' """
        desired_member = None
        for member in [self] + self.edp_list:
            if name == member.drizzle_filename:
                desired_member = member
                break
        return desired_member

    def add_member(self, edp):

        """ Add an ExposureProduct object to the list - composition.
        """
        self.edp_list.append(edp)

    def wcs_drizzle_product(self, meta_wcs):
        """
            Create the drizzle-combined filter image using the meta_wcs as the reference output
        """
        # This insures that keywords related to the footprint are generated for this
        # specific object to use in updating the output drizzle product.
        self.meta_wcs = meta_wcs
        if self.mask is None:
            self.generate_footprint_mask()

        # Retrieve the configuration parameters for astrodrizzle
        drizzle_pars = self.configobj_pars.get_pars("astrodrizzle")
        # ...and set parameters which are computed on-the-fly
        drizzle_pars["final_refimage"] = meta_wcs
        drizzle_pars["runfile"] = self.trl_logname
        # Setting "preserve" to false so the OrIg_files directory is deleted as the purpose
        # of this directory is now obsolete.
        drizzle_pars["preserve"] = False
        drizzle_pars['rules_file'] = self.rules_file

        log.debug("The 'final_refimage' ({}) and 'runfile' ({}) configuration variables "
                  "have been updated for the drizzle step of the filter drizzle product."
                  .format(meta_wcs, self.trl_logname))

        edp_filenames = [element.full_filename for element in self.edp_list]

        if len(edp_filenames) == 1:
            drizzle_pars['resetbits'] = "0"  # Use any pixels already flagged as CRs

        astrodrizzle.AstroDrizzle(input=edp_filenames,
                                  output=self.drizzle_filename,
                                  **drizzle_pars)

        # Update product with SVM-specific keywords based on the footprint
        with fits.open(self.drizzle_filename, mode='update') as hdu:
            for kw in self.mask_kws:
                hdu[("SCI", 1)].header[kw] = tuple(self.mask_kws[kw])

        # Rename Astrodrizzle log file as a trailer file
        log.debug("Filter combined image {} composed of: {}".format(self.drizzle_filename, edp_filenames))
        try:
            shutil.move(self.trl_logname, self.trl_filename)
        except PermissionError:
            # TODO:  trailer filename should be saved for moving later...
            pass


class ExposureProduct(HAPProduct):
    """ An Exposure Product is an individual exposure/image (flt/flc).

        The "edp" is short hand for ExposureProduct.
    """
    def __init__(self, prop_id, obset_id, instrument, detector, filename, filters, filetype, log_level):
        super().__init__(prop_id, obset_id, instrument, detector, filename, filetype, log_level)

        self.info = '_'.join([prop_id, obset_id, instrument, detector, filename, filters, filetype])
        self.filters = filters
        self.full_filename = self.copy_exposure(filename)

        # Open the input FITS file to mine some header information.
        hdu_list = fits.open(filename)
        self.mjdutc = hdu_list[0].header['EXPSTART']
        self.exptime = hdu_list[0].header['EXPTIME']
        hdu_list.close()

        self.product_basename = self.basename + "_".join(map(str, [filters, self.exposure_name]))
        self.drizzle_filename = self.product_basename + "_" + self.filetype + ".fits"
        self.headerlet_filename = self.product_basename + "_hlet.fits"
        self.trl_logname = self.product_basename + "_trl.log"
        self.trl_filename = self.product_basename + "_trl.txt"

        self.regions_dict = {}

        # Define HAPLEVEL value for this product
        self.haplevel = 1

        # This attribute is set in poller_utils.py
        self.is_singleton = False

        # Flag whether this exposure is being processed for the 'first' time or not
        self.new_process = True

        # Flag whether to use single-image CR identification with this exposure
        self.crclean = False

        log.info("Exposure object {} created.".format(self.full_filename[0:9]))

    def find_member(self, name):
        """ Return member instance with filename 'name' """
        if name == self.drizzle_filename:
            return self
        else:
            return None

    def __getattribute__(self, name):
        if name in ["generate_footprint_mask", "generate_metawcs", "meta_wcs", "mask_kws", "mask"]:
            raise AttributeError(name)
        else:
            return super(ExposureProduct, self).__getattribute__(name)

    def __dir__(self):
        class_set = (set(dir(self.__class__)) | set(self.__dict__.keys()))
        unwanted_set = set(["generate_footprint_mask", "generate_metawcs", "meta_wcs", "mask_kws", "mask"])
        return sorted(class_set - unwanted_set)

    def wcs_drizzle_product(self, meta_wcs):
        """
            Create the drizzle-combined exposure image using the meta_wcs as the reference output
        """
        # Retrieve the configuration parameters for astrodrizzle
        drizzle_pars = self.configobj_pars.get_pars("astrodrizzle")
        # ...and set parameters which are computed on-the-fly
        drizzle_pars["final_refimage"] = meta_wcs
        drizzle_pars["runfile"] = self.trl_logname
        # Setting "preserve" to false so the OrIg_files directory is deleted as the purpose
        # of this directory is now obsolete.
        drizzle_pars["preserve"] = False
        drizzle_pars['rules_file'] = self.rules_file
        drizzle_pars['resetbits'] = "0"

        log.debug("The 'final_refimage' ({}) and 'runfile' ({}) configuration variables "
                  "have been updated for the drizzle step of the exposure drizzle product."
                  .format(meta_wcs, self.trl_logname))

        astrodrizzle.AstroDrizzle(input=self.full_filename,
                                  output=self.drizzle_filename,
                                  **drizzle_pars)

        # Rename Astrodrizzle log file as a trailer file
        log.debug("Exposure image {}".format(self.drizzle_filename))
        try:
            shutil.move(self.trl_logname, self.trl_filename)
        except PermissionError:
            pass

    def copy_exposure(self, filename):
        """
            Create a copy of the original input to be renamed and used for single-visit processing.

            New exposure filename needs to follow the convention:
            hst_<propid>_<obsetid>_<instr>_<detector>_<filter>_<ipppssoo>_fl[ct].fits

            Parameters
            ----------
            filename : str
                Original pipeline filename for input exposure

            Returns
            -------
            edp_filename : str
                New SVM-compatible HAP filename for input exposure

        """
        suffix = filename.split("_")[1]
        edp_filename = self.basename + \
            "_".join(map(str, [self.filters, filename[:8], suffix]))

        log.info("Copying {} to SVM input: \n    {}".format(filename, edp_filename))
        try:
            shutil.copy(filename, edp_filename)
        except PermissionError:
            pass

        # Add HAP keywords as required by pipeline processing
        with fits.open(edp_filename, mode='update') as edp_hdu:
            edp_hdu[0].header['HAPLEVEL'] = (0, 'Classification level of this product')
            edp_hdu[0].header['IPPPSSOO'] = edp_hdu[0].header['ROOTNAME']
            edp_hdu[0].header['FILENAME'] = edp_filename

        return edp_filename


class GrismExposureProduct(HAPProduct):
    """ A Grism Exposure Product is an individual Grism/Prism exposure/image (flt/flc).

        The "grism_edp" is short hand for GrismExposureProduct.
    """
    def __init__(self, prop_id, obset_id, instrument, detector, filename, filters, filetype, log_level):
        super().__init__(prop_id, obset_id, instrument, detector, filename, filetype, log_level)

        self.info = '_'.join([prop_id, obset_id, instrument, detector, filename, filters, filetype])
        self.filters = filters
        self.full_filename = self.copy_exposure(filename)

        # Open the input FITS file to mine some header information.
        # and make sure the WCS is up-to-date
        hdu_list = fits.open(filename)
        self.mjdutc = hdu_list[0].header['EXPSTART']
        self.exptime = hdu_list[0].header['EXPTIME']

        drizcorr = hdu_list[0].header['DRIZCORR']
        if drizcorr == "OMIT":
            updatewcs.updatewcs(self.full_filename, use_db=True)
        hdu_list.close()

        self.product_basename = self.basename + "_".join(map(str, [filters, self.exposure_name]))
        self.drizzle_filename = self.product_basename + "_" + self.filetype + ".fits"
        self.headerlet_filename = self.product_basename + "_hlet.fits"
        self.trl_logname = self.product_basename + "_trl.log"
        self.trl_filename = self.product_basename + "_trl.txt"

        self.regions_dict = {}

        # Define HAPLEVEL value for this product
        self.haplevel = 1

        log.info("Grism Exposure object {} created.".format(self.full_filename))

    def copy_exposure(self, filename):
        """
            Create a copy of the original input to be renamed and used for single-visit processing.

            New exposure filename needs to follow the convention:
            hst_<propid>_<obsetid>_<instr>_<detector>_<filter>_<ipppssoo>_fl[ct].fits

            Parameters
            ----------
            filename : str
                Original pipeline filename for input exposure

            Returns
            -------
            edp_filename : str
                New SVM-compatible HAP filename for input exposure

        """
        suffix = filename.split("_")[1]
        edp_filename = self.basename + \
            "_".join(map(str, [self.filters, filename[:8], suffix]))

        log.info("Copying {} to SVM input: \n    {}".format(filename, edp_filename))
        try:
            shutil.copy(filename, edp_filename)
        except PermissionError:
            pass

        # Add HAP keywords as required by pipeline processing
        with fits.open(edp_filename, mode='update') as edp_hdu:
            edp_hdu[0].header['HAPLEVEL'] = (0, 'Classification level of this product')
            edp_hdu[0].header['IPPPSSOO'] = edp_hdu[0].header['ROOTNAME']
            edp_hdu[0].header['FILENAME'] = edp_filename

        return edp_filename


class SkyCellExposure(HAPProduct):
    """ An SkyCell Exposure Product is an individual exposure/image (flt/flc).

        The "sce" is short hand for ExposureProduct.
    """
    def __init__(self, prop_id, obset_id, instrument, detector, filename, layer, filetype, log_level):
        super().__init__(prop_id, obset_id, instrument, detector, filename, filetype, log_level)

        filter_str = layer[0]

        # parse layer information into filename layer_str
        # layer: [filter_str, pscale_str, exptime_str, epoch_str]
        # e.g.: [f160w, coarse, all, all]
        #
        if filename.startswith("hst"):
            exposure_name = filename.split("_")[-2]
        else:
            exposure_name = self.exposure_name
        if layer[1] == 'coarse':
            layer_vals = [layer[1], layer[2], exposure_name]
        else:
            layer_vals = ['all', exposure_name]
        layer_str = '-'.join(layer_vals)

        cell_id = "p{}{}".format(prop_id, obset_id)
        self.basename = "hst_skycell-" + "_".join(map(str, [cell_id, instrument, detector])) + "_"

        self.info = self.basename + '_'.join([filter_str, layer_str, filename, filetype])
        self.filters = filter_str

        # Open the input FITS file to mine some header information.
        hdu_list = fits.open(filename)
        self.mjdutc = hdu_list[0].header['EXPSTART']
        self.exptime = hdu_list[0].header['EXPTIME']
        hdu_list.close()

        self.product_basename = self.basename + "_".join(map(str, [filter_str, layer_str]))
        self.drizzle_filename = self.product_basename + "_" + self.filetype + ".fits"
        self.headerlet_filename = self.product_basename + "_hlet.fits"
        self.trl_logname = self.product_basename + "_trl.log"
        self.trl_filename = self.product_basename + "_trl.txt"

        self.full_filename = self.copy_exposure(filename)

        self.regions_dict = {}

        # Define HAPLEVEL value for this product
        self.haplevel = 1

        # This attribute is set in poller_utils.py
        self.is_singleton = False

        # Flag whether this exposure is being processed for the 'first' time or not
        self.new_process = True

        # Flag whether to use single-image CR identification with this exposure
        self.crclean = False

        log.info("Create SkyCellExposure object:\n    {}".format(self.full_filename))

    def find_member(self, name):
        """ Return member instance with filename 'name' """
        if name == self.drizzle_filename:
            return self
        else:
            return None

    def __getattribute__(self, name):
        if name in ["generate_footprint_mask", "generate_metawcs", "meta_wcs", "mask_kws", "mask"]:
            raise AttributeError(name)
        else:
            return super(SkyCellExposure, self).__getattribute__(name)

    def __dir__(self):
        class_set = (set(dir(self.__class__)) | set(self.__dict__.keys()))
        unwanted_set = set(["generate_footprint_mask", "generate_metawcs", "meta_wcs", "mask_kws", "mask"])
        return sorted(class_set - unwanted_set)

    def wcs_drizzle_product(self, meta_wcs):
        """
            Create the drizzle-combined exposure image using the meta_wcs as the reference output
        """
        # Retrieve the configuration parameters for astrodrizzle
        drizzle_pars = self.configobj_pars.get_pars("astrodrizzle")
        # ...and set parameters which are computed on-the-fly
        drizzle_pars["final_refimage"] = meta_wcs
        drizzle_pars["runfile"] = self.trl_logname
        # Setting "preserve" to false so the OrIg_files directory is deleted as the purpose
        # of this directory is now obsolete.
        drizzle_pars["preserve"] = False
        drizzle_pars['rules_file'] = self.rules_file
        drizzle_pars['resetbits'] = "0"
        # Set 'context' to False in order to turn off creation of CTX extensions in MVM products
        # to save on disk space and memory use during processing.
        drizzle_pars['context'] = False

        log.debug("The 'final_refimage' ({}) and 'runfile' ({}) configuration variables "
                  "have been updated for the drizzle step of the exposure drizzle product."
                  .format(meta_wcs, self.trl_logname))

        astrodrizzle.AstroDrizzle(input=self.full_filename,
                                  output=self.drizzle_filename,
                                  **drizzle_pars)

        # Rename Astrodrizzle log file as a trailer file
        log.debug("Exposure image {}".format(self.drizzle_filename))

        try:
            shutil.move(self.trl_logname, self.trl_filename)
        except PermissionError:
            pass

    def copy_exposure(self, filename):
        """
            Create a copy of the original input to be renamed and used for multi-visit processing.

            New exposure filename needs to follow the MVM naming convention:
            hst_skycell-p<PPPP>x<XX>y<YY>_<instr>_<detector>_<filter>-<layer>_<ipppssoo>_fl[ct].fits

            Parameters
            ----------
            filename : str
                Original pipeline filename for input exposure

            Returns
            -------
            sce_filename : str
                New MVM-compatible HAP filename for input exposure

        """
        if filename.startswith("hst"):
            sce_filename = "{}_{}".format(self.product_basename, filename.split("_")[-1])
        else:
            suffix = filename.split("_")[1]
            sce_filename = '_'.join([self.product_basename, suffix])
        log.info("Copying {} to MVM input: \n    {}".format(filename, sce_filename))
        try:
            shutil.copy(filename, sce_filename)
        except PermissionError:
            pass

        # Add HAPLEVEL keyword as required by pipeline processing
        fits.setval(sce_filename, 'HAPLEVEL', value=0, comment='Classification level of this product')

        return sce_filename


class SkyCellProduct(HAPProduct):
    """ A SkyCell Product is a mosaic comprised of images acquired
        during a multiple visits with one instrument, one detector, a single filter,
        and all exposure times.  The exposure time characteristic may change - TBD.

        The "scp" is short hand for SkyCellProduct.
    """
    def __init__(self, prop_id, obset_id, instrument, detector, skycell_name, layer, filetype, log_level):
        super().__init__(prop_id, obset_id, instrument, detector, skycell_name, filetype, log_level)
        # May need to exclude 'filter' component from layer_str
        filter_str = layer[0]
        self.filters = filter_str

        # parse layer information into filename layer_str
        # layer: [filter_str, pscale_str, exptime_str, epoch_str]
        # e.g.: [f160w, coarse, all, all]
        #
        if layer[1] == 'coarse':
            layer_str = '-'.join([layer[1], layer[2]])
        else:
            layer_str = 'all'

        layer_scale = layer[1]

        self.info = '_'.join(['hst', skycell_name, instrument, detector, filter_str, layer_str])
        self.exposure_name = skycell_name
        self.product_basename = self.info

        # Trailer names .txt or .log
        self.trl_logname = self.product_basename + "_trl.log"
        self.trl_filename = self.product_basename + "_trl.txt"
        # Initialize the output catalog name attributes (for the time when they are created)
        self.point_cat_filename = None
        self.segment_cat_filename = None

        self.drizzle_filename = '_'.join([self.product_basename, self.filetype]) + ".fits"
        self.refname = self.product_basename + "_ref_cat.ecsv"

        # Generate the name for the manifest file which is for the entire multi-visit.  It is fine
        # to use only one of the SkyCellProducts to generate the manifest name as the name
        # is only dependent on the sky cell.
        # Example: hst_skycell-p<PPPP>x<XX>y<YY>_manifest.txt (e.g., hst_skycell-p0797x12y05_manifest.txt)
        self.manifest_name = '_'.join(['hst', skycell_name, 'manifest.txt'])

        # Define HAPLEVEL value for this product
        self.haplevel = 4

        # These attributes will be populated during processing
        self.edp_list = []
        self.new_to_layer = 0
        self.regions_dict = {}
        self.skycell = cell_utils.SkyCell.from_name(skycell_name, scale=layer_scale)
        self.configobj_pars = None
        self.meta_wcs = None

        self.all_mvm_exposures = []
        self.meta_bounded_wcs = None

        log.debug("SkyCell object {}/{}/{} created.".format(self.instrument, self.detector, self.filters))

    def find_member(self, name):
        """ Return member instance with filename 'name' """
        desired_member = None
        for member in [self] + self.edp_list:
            if name == member.drizzle_filename:
                desired_member = member
                break
        return desired_member

    def add_member(self, edp):
        """ Add an ExposureProduct object to the list - composition.
        """
        self.edp_list.append(edp)
        self.new_to_layer += edp.new_process

    def add_all_mvm_exposures_list(self, exp_list):
        """ Add a list containing all the MVM FLT or FLC filenames, even the
            filenames for exposures which have been previously processed.
        """
        self.all_mvm_exposures = exp_list

    def generate_metawcs(self, custom_limits=None):
        """Generate a meta wcs

        Parameters
        ----------
        custom_limits : list, optional.
            a 4-element list containing the mosaic bounding rectangle X min and max and Y min and max values
            This input argument is only used for creation of custom mosaics. These coordinates are in the
            frame of reference of the projection cell, at fine (platescale = 0.04 arcsec/pixel) resolution.
        """
        if custom_limits:  # for creation of custom multi-skycell mosaics, base meta_wcs on projection cell WCS
            wcs = copy.deepcopy(self.skycell.projection_cell.wcs)
            ratio = self.skycell.projection_cell.wcs.pscale / self.skycell.scale

            # Store unscaled and scaled versions of the custom mosaic bounding box X and Y limits
            xmin_unscaled = int(np.rint(custom_limits[0]))
            xmax_unscaled = int(np.rint(custom_limits[1]))
            ymin_unscaled = int(np.rint(custom_limits[2]))
            ymax_unscaled = int(np.rint(custom_limits[3]))
            xmin_scaled = int(np.rint(custom_limits[0] * ratio))
            xmax_scaled = int(np.rint(custom_limits[1] * ratio))
            ymin_scaled = int(np.rint(custom_limits[2] * ratio))
            ymax_scaled = int(np.rint(custom_limits[3] * ratio))

            # Adjust WCS values as needed
            crpix1 = int(np.rint((wcs.wcs.crpix[0] - xmin_unscaled) * ratio))
            crpix2 = int(np.rint((wcs.wcs.crpix[1] - ymin_unscaled) * ratio))
            wcs.wcs.crpix = [crpix1, crpix2]
            wcs.wcs.crval = wcs.wcs.crval
            wcs.wcs.cd = wcs.wcs.cd / ratio
            wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
            wcs.pscale = wcs.pscale / ratio
            wcs._naxis[0] = int(np.rint(wcs._naxis[0] * ratio))
            wcs._naxis[1] = int(np.rint(wcs._naxis[1] * ratio))
            self.bounding_box = [slice(ymin_scaled, ymax_scaled), slice(xmin_scaled, xmax_scaled)]
            wcs.pixel_shape = [xmax_scaled - xmin_scaled + 1, ymax_scaled - ymin_scaled + 1]
        else:  # For regular MVM processing (single skycell), base meta_wcs on skycell wcs
            wcs = copy.deepcopy(self.skycell.wcs)

        # This is the exposure-independent WCS.
        self.meta_wcs = wcs

        # Create footprint on the sky for all input exposures using the skycell wcs
        # This footprint includes all the exposures in the visit, NEW exposures, as well
        # as exposures which have been previously processed (all are listed in the original
        # poller file).
        mvm_footprint = cell_utils.SkyFootprint(wcs)
        mvm_footprint.build(self.all_mvm_exposures)

        # This is the exposure-dependent WCS.
        self.meta_bounded_wcs = mvm_footprint.bounded_wcs
        return self.meta_bounded_wcs

    def wcs_drizzle_product(self, meta_wcs):
        """
            Create the drizzle-combined sky-cell layer image using the meta_wcs as the reference output
        """
        # This insures that keywords related to the footprint are generated for this
        # specific object to use in updating the output drizzle product.
        self.meta_wcs = meta_wcs
        if self.mask is None:
            self.generate_footprint_mask()

        # Retrieve the configuration parameters for astrodrizzle
        drizzle_pars = self.configobj_pars.get_pars("astrodrizzle")
        # ...and set parameters which are computed on-the-fly
        drizzle_pars["final_refimage"] = meta_wcs
        drizzle_pars["runfile"] = self.trl_logname
        # Setting "preserve" to false so the OrIg_files directory is deleted as the purpose
        # of this directory is now obsolete.
        drizzle_pars["preserve"] = False
        # Turn off generation of CTX extension (not useful for these products)
        drizzle_pars["context"] = False
        log.debug("The 'final_refimage' ({}) and 'runfile' ({}) configuration variables "
                  "have been updated for the drizzle step of the filter drizzle product."
                  .format(meta_wcs, self.trl_logname))

        edp_filenames = [element.full_filename for element in self.edp_list]

        astrodrizzle.AstroDrizzle(input=edp_filenames,
                                  output=self.drizzle_filename,
                                  **drizzle_pars)

        # Update product with SVM-specific keywords based on the footprint
        with fits.open(self.drizzle_filename, mode='update') as hdu:
            for kw in self.mask_kws:
                hdu[("SCI", 1)].header[kw] = tuple(self.mask_kws[kw])

        # Rename Astrodrizzle log file as a trailer file
        log.debug("Sky-cell layer image {} composed of: {}".format(self.drizzle_filename, edp_filenames))
        try:
            shutil.move(self.trl_logname, self.trl_filename)
        except PermissionError:
            pass
