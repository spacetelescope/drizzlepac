""" Definition of Super and Subclasses for the mosaic output image_list

    Classes which define the total ("white light" image), filter, and exposure
    drizzle products.

    These products represent different levels of processing with the levels noted in the
    'HAPLEVEL' keyword.  The 'HAPLEVEL' values are:

      * 1 : calibrated (FLT/FLC) input images and exposure level drizzle products with improved astrometry
      * 2 : filter and total products combined using the improved astrometry,
            consistent pixel scale, and oriented to North.
      * 3 : (future) multi-visit mosaics aligned to common tangent plane
"""
import logging
import sys
import os
import traceback
import shutil

from stsci.tools import logutil
from astropy.io import fits
import numpy as np

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
            generated based upon the merging of all the ExposureProducts which comprise the
            specific TotalProduct.  This is done on a per TotalProduct basis as the Single
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

        log.debug("Total detection object {}/{} created.".format(self.instrument, self.detector))

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
        self.refname = self.product_basename + "_ref_cat.ecsv"

        # Define HAPLEVEL value for this product
        self.haplevel = 2

        # These attributes will be populated during processing
        self.edp_list = []
        self.regions_dict = {}

        log.debug("Filter object {}/{}/{} created.".format(self.instrument, self.detector, self.filters))

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

    def align_to_gaia(self, catalog_name='GAIADR2', headerlet_filenames=None, output=True,
                        fit_label='EVM', align_table=None, fitgeom='rscale'):
        """Extract the flt/flc filenames from the exposure product list, as
           well as the corresponding headerlet filenames to use legacy alignment
           routine.
        """
        log.info('Starting alignment to absolute astrometric reference frame {}'.format(catalog_name))
        alignment_pars = self.configobj_pars.get_pars('alignment')

        # Only perform the relative alignment
        method_name = 'relative'

        exposure_filenames = []
        headerlet_filenames = {}
        align_table = None
        crclean = []

        try:
            if self.edp_list:
                for edp in self.edp_list:
                    exposure_filenames.append(edp.full_filename)
                    headerlet_filenames[edp.full_filename] = edp.headerlet_filename
                    crclean.append(edp.crclean)

                if align_table is None:
                    align_table = align_utils.AlignmentTable(exposure_filenames,
                                                             log_level=self.log_level,
                                                             **alignment_pars)

                    align_table.find_alignment_sources(output=output, crclean=crclean)
                    align_table.configure_fit()
                log.debug('Creating reference catalog {}'.format(self.refname))
                ref_catalog = amutils.create_astrometric_catalog(align_table.process_list,
                                                                 catalog=catalog_name,
                                                                 output=self.refname,
                                                                 gaia_only=False)

                log.debug("Abbreviated reference catalog displayed below\n{}".format(ref_catalog))
                align_table.reference_catalogs[self.refname] = ref_catalog
                if len(ref_catalog) > align_utils.MIN_CATALOG_THRESHOLD:
                    fit_again = False
                    try:
                        align_table.perform_fit(method_name, catalog_name, ref_catalog,
                                                fitgeom=fitgeom)
                    except Exception:
                        log.info('Problem with initial fit...')
                        fit_again = True

                    # If there was a problem with the first fit...
                    if fit_again:
                        log.info("Trying 'DEFAULT' fit to {} instead.".format(catalog_name))
                        # Try one more time with slightly different cross-matching algorithm
                        method_name = 'default'
                        align_table.configure_fit()
                        align_table.perform_fit(method_name, catalog_name, ref_catalog,
                                                fitgeom=fitgeom)

                    align_table.select_fit(catalog_name, method_name)
                    align_table.apply_fit(headerlet_filenames=headerlet_filenames,
                                         fit_label=fit_label)
                else:
                    log.warning("Not enough reference sources for absolute alignment...")
                    raise ValueError

        except Exception:
            # Report a problem with the alignment
            log.warning("EXCEPTION encountered in align_to_gaia for the FilteredProduct.\n")
            log.warning("No correction to absolute astrometric frame applied.\n")
            log.warning("Proceeding with previous best solution.\n")

            # Only write out the traceback if in "debug" mode since not being able to
            # align the data to an absolute astrometric frame is not actually a failure.
            log.info(traceback.format_exc())
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


class SkyCellProduct(HAPProduct):
    """ A SkyCell Product is a mosaic comprised of images acquired
        during a multiple visits with one instrument, one detector, a single filter,
        and all exposure times.  The exposure time characteristic may change - TBD.

        The "scp" is short hand for SkyCellProduct.
    """
    def __init__(self, prop_id, obset_id, instrument, detector, skycell_name, layer, filetype, log_level):
        super().__init__(prop_id, obset_id, instrument, detector, skycell_name, filetype, log_level)
        layer_str = '-'.join(layer)

        self.info = '_'.join([skycell_name, instrument, detector, layer_str])

        self.exposure_name = skycell_name
        filters = layer[0]
        self.product_basename = self.info

        self.filters = filters

        # Trailer names .txt or .log
        self.trl_logname = self.product_basename + "_trl.log"
        self.trl_filename = self.product_basename + "_trl.txt"
        self.point_cat_filename = None
        self.segment_cat_filename = None
        self.drizzle_filename = '_'.join(['hst', self.product_basename, self.filetype]) + ".fits"
        self.refname = self.product_basename + "_ref_cat.ecsv"
        # Generate the name for the manifest file which is for the entire visit.  It is fine
        # to create it as an attribute of a TotalProduct as it is independent of
        # the detector in use.
        # instrument_programID_obsetID_manifest.txt (e.g.,wfc3_b46_06_manifest.txt)
        self.manifest_name = '_'.join(['hst', self.product_basename, "manifest.txt"])

        # Define HAPLEVEL value for this product
        self.haplevel = 3

        # These attributes will be populated during processing
        self.edp_list = []
        self.new_to_layer = 0
        self.regions_dict = {}
        self.skycell = cell_utils.SkyCell(name=skycell_name)
        self.configobj_pars = None

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

    def generate_metawcs(self):
        self.meta_wcs = self.skycell.wcs
        return self.meta_wcs

    def align_to_gaia(self, catalog_name='GAIADR2', headerlet_filenames=None, output=True,
                        fit_label='MVM', align_table=None, fitgeom='rscale'):
        """Extract the flt/flc filenames from the exposure product list, as
           well as the corresponding headerlet filenames to use legacy alignment
           routine.
        """
        log.info('Starting alignment to absolute astrometric reference frame {}'.format(catalog_name))
        alignment_pars = self.configobj_pars.get_pars('alignment')

        # Only perform the relative alignment
        method_name = 'relative'

        exposure_filenames = []
        headerlet_filenames = {}
        align_table = None
        try:
            if self.edp_list:
                for edp in self.edp_list:
                    exposure_filenames.append(edp.full_filename)
                    headerlet_filenames[edp.full_filename] = edp.headerlet_filename

                if align_table is None:
                    align_table = align_utils.AlignmentTable(exposure_filenames,
                                                             log_level=self.log_level,
                                                             **alignment_pars)
                    align_table.find_alignment_sources(output=output)
                    align_table.configure_fit()
                log.debug('Creating reference catalog {}'.format(self.refname))
                ref_catalog = amutils.create_astrometric_catalog(align_table.process_list,
                                                                 catalog=catalog_name,
                                                                 output=self.refname,
                                                                 gaia_only=False)

                log.debug("Abbreviated reference catalog displayed below\n{}".format(ref_catalog))
                align_table.reference_catalogs[self.refname] = ref_catalog
                if len(ref_catalog) > align_utils.MIN_CATALOG_THRESHOLD:
                    align_table.perform_fit(method_name, catalog_name, ref_catalog,
                                           fitgeom=fitgeom)
                    align_table.select_fit(catalog_name, method_name)
                    align_table.apply_fit(headerlet_filenames=headerlet_filenames,
                                         fit_label=fit_label)
                else:
                    log.warning("Not enough reference sources for absolute alignment...")
                    raise ValueError

        except Exception:
            # Report a problem with the alignment
            log.warning("EXCEPTION encountered in align_to_gaia for the FilteredProduct.\n")
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
