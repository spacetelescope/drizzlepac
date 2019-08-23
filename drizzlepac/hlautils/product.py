""" Definition of Super and Subclasses for the mosaic output image_list

    Classes which define the total ("white light" image), filter, and exposure
    drizzle products.
"""
import sys
import traceback
import shutil

from drizzlepac import wcs_functions
from drizzlepac import alignimages
from drizzlepac import astrodrizzle

class HAPProduct:
    """ HAPProduct is the base class for the various products generated during the
        astrometry update and mosaicing of image data.
    """
    def __init__(self, prop_id, obset_id, instrument, detector, filename, filetype):
        # Make sure the proposal ID is a 5-character string
        self.prop_id = prop_id.zfill(5)
        self.obset_id = obset_id
        self.instrument = instrument
        self.detector = detector
        self.filetype = filetype

        self.basename = "hst_" + "_".join(map(str, [prop_id, obset_id, instrument, detector])) + "_"

        # exposure_name is the ipppssoo or a portion thereof
        self.exposure_name = filename[0:8]

        # TO DO: update this variable
        self.mjdutc = None

        # class variables to be set later
        self.pars = None

    # def print_info(self):
        # """ Generic print at this time to indicate the information used in the
        #     construction of the object for debug purposes.
        # """
        # print("Object information: {}".format(self.info))

    def update_config(self, param_dict):
        """ Method reserved for potential update of configuration information
            for the particular data product in question.
        """
        pass

class TotalProduct(HAPProduct):
    """ A Total Detection Product is a 'white' light mosaic comprised of
        images acquired with one instrument, one detector, all filters, and all
        exposure times.  The exposure time characteristic may change - TBD.

        The "tdp" is short hand for TotalProduct.
    """
    def __init__(self, prop_id, obset_id, instrument, detector, filename, filetype):
        super().__init__(prop_id, obset_id, instrument, detector, filename, filetype)

        self.info = '_'.join([prop_id, obset_id, instrument, detector, filename, filetype])
        self.exposure_name = filename[0:6]

        self.product_basename = self.basename + "total_" + self.exposure_name
        self.trl_filename = self.product_basename + "_trl.log"
        self.point_cat_filename = self.product_basename + "_point-cat.ecsv"
        self.segment_cat_filename = self.product_basename + "_segment-cat.ecsv"
        self.drizzle_filename = self.product_basename + "_" + self.filetype + ".fits"
        self.ref_drizzle_filename = self.product_basename + "_ref_" + self.filetype + ".fits"

        # Generate the name for the manifest file which is for the entire visit.  It is fine
        # to create it as an attribute of a TotalProduct as it is independent of
        # the detector in use.
        # instrument_programID_obsetID_manifest.txt (e.g.,wfc3_b46_06_manifest.txt)
        self.manifest_name = '_'.join([instrument, filename[1:4], obset_id, "manifest.txt"])

        # These attributes will be populated during processing
        self.edp_list = []
        self.fdp_list = []
        self.regions_dict = {}
        self.meta_wcs = None

    def add_member(self, edp):
        """ Add an ExposureProduct object to the list - composition.
        """
        self.edp_list.append(edp)

    def add_product(self, fdp):
        """ Add a FilterProduct object to the list - composition.
        """
        self.fdp_list.append(fdp)

    def build_metawcs(self):
        """ A reserved method to build a unique WCS for each TotalProduct product which is
            generated based upon the merging of all the ExposureProductss which comprise the
            specific TotalProduct.

        image_list = [element.full_filename for element in edp_list]
        meta_wcs = wcs_functions.make_mosaic_wcs(image_list)
        outnx = meta_wcs.pixel_shape[1]
        outny = meta_wcs.pixel_shape[0]
        """
        pass

    def wcs_drizzle_product(self, meta_wcs, configobj):
        """
        Create the drizzle-combined total image using the meta_wcs as the reference output
        """
        edp_filenames = [element.full_filename for element in self.edp_list]
        astrodrizzle.AstroDrizzle(input=edp_filenames,
                                  output=self.drizzle_filename,
                                  final_refimage=meta_wcs,
                                  final_outnx=meta_wcs.pixel_shape[1],
                                  final_outny=meta_wcs.pixel_shape[0],
                                  runfile=self.trl_filename,
                                  configobj=configobj)

        # Rename Astrodrizzle log file as a trailer file
        shutil.move(self.trl_filename, self.trl_filename.replace('.log', '.txt'))

    def image_drizzle_product(self, ref_image, total_shape, configobj):
        """
        Create the drizzle-combined total image using the reference image as the reference output
        """
        edp_filenames = [element.full_filename for element in self.edp_list]
        astrodrizzle.AstroDrizzle(input=edp_filenames,
                                  output=self.drizzle_filename,
                                  final_refimage=ref_image,
                                  final_outnx=total_shape[1],
                                  final_outny=total_shape[0],
                                  runfile=self.trl_filename,
                                  configobj=configobj)

        # Rename Astrodrizzle log file as a trailer file
        shutil.move(self.trl_filename, self.trl_filename.replace('.log', '.txt'))

        # log.info("Total combined image consists of ... {} {}".format(self.drizzle_filename, edp_filenames))


class FilterProduct(HAPProduct):
    """ A Filter Detection Product is a mosaic comprised of images acquired
        during a single visit with one instrument, one detector, a single filter,
        and all exposure times.  The exposure time characteristic may change - TBD.

        The "fdp" is short hand for FilterProduct.
    """
    def __init__(self, prop_id, obset_id, instrument, detector, filename, filters, filetype):
        super().__init__(prop_id, obset_id, instrument, detector, filename, filetype)

        self.info = '_'.join([prop_id, obset_id, instrument, detector, filename, filters, filetype])
        self.exposure_name = filename[0:6]
        self.filters = filters

        self.product_basename = self.basename + "_".join(map(str, [filters, self.exposure_name]))
        # Trailer names .txt or .log
        self.trl_filename = self.product_basename + "_trl.log"
        self.point_cat_filename = self.product_basename + "_point-cat.ecsv"
        self.segment_cat_filename = self.product_basename + "_segment-cat.ecsv"
        self.drizzle_filename = self.product_basename + "_" + self.filetype + ".fits"

        # These attributes will be populated during processing
        self.edp_list = []
        self.regions_dict = {}

    def add_member(self, edp):
        """ Add an ExposureProduct object to the list - composition.
        """
        self.edp_list.append(edp)

    def align_to_gaia(self):
        """Extract the flt/flc filenames from the exposure product list, as
           well as the corresponding headerlet filenames to use legacy alignment
           routine.
        """
        exposure_filenames = []
        headerlet_filenames = {}
        align_table = None
        try:
            if self.edp_list:
                for edp in self.edp_list:
                    exposure_filenames.append(edp.full_filename)
                    headerlet_filenames[edp.full_filename] = edp.headerlet_filename

                align_table = alignimages.perform_align(exposure_filenames,
                                                        debug=False,
                                                        runfile="alignimages.log",
                                                        update_hdr_wcs=True,
                                                        headerlet_filenames=headerlet_filenames)

        except Exception:
            # TODO: Fix up the logging
            # Report a problem with the alignment
            # log.info("EXCEPTION encountered in alignimages.\n")
            exc_type, exc_value, exc_tb = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)
            # log.info("No correction to absolute astrometric frame applied.\n")
            align_table = None

        # Return a table which contains data regarding the alignment, as well as the
        # list of the flt/flc exposures which were part of the alignment process
        # TODO: This does not account for individual exposures which might have been
        # excluded from alignment.
        return align_table, exposure_filenames

    def wcs_drizzle_product(self, meta_wcs, configobj):
        """
        Create the drizzle-combined filter image using the meta_wcs as the reference output
        """
        edp_filenames = [element.full_filename for element in self.edp_list]
        astrodrizzle.AstroDrizzle(input=edp_filenames,
                                  output=self.drizzle_filename,
                                  final_refimage=meta_wcs,
                                  final_outnx=meta_wcs.pixel_shape[1],
                                  final_outny=meta_wcs.pixel_shape[0],
                                  runfile=self.trl_filename,
                                  configobj=configobj)

        # Rename Astrodrizzle log file as a trailer file
        shutil.move(self.trl_filename, self.trl_filename.replace('.log', '.txt'))

    def image_drizzle_product(self, ref_image, total_shape, configobj):
        """
        Create the drizzle-combined filter image using the reference image as the reference output
        """
        edp_filenames = [element.full_filename for element in self.edp_list]
        astrodrizzle.AstroDrizzle(input=edp_filenames,
                                  output=self.drizzle_filename,
                                  final_refimage=ref_image,
                                  final_outnx=total_shape[1],
                                  final_outny=total_shape[0],
                                  runfile=self.trl_filename,
                                  configobj=configobj)

        # Rename Astrodrizzle log file as a trailer file
        shutil.move(self.trl_filename, self.trl_filename.replace('.log', '.txt'))

        # log.info("Filter combined image... {} {}".format(self.drizzle_filename, edp_filenames))

class ExposureProduct(HAPProduct):
    """ An Exposure Product is an individual exposure/image (flt/flc).

        The "edp" is short hand for ExposureProduct.
    """
    def __init__(self, prop_id, obset_id, instrument, detector, filename, filters, filetype):
        super().__init__(prop_id, obset_id, instrument, detector, filename, filetype)

        self.info = '_'.join([prop_id, obset_id, instrument, detector, filename, filters, filetype])
        self.full_filename = filename
        self.filters = filters

        # TO DO: update this variable
        # self.exptime = exptime

        self.product_basename = self.basename + "_".join(map(str, [filters, self.exposure_name]))
        self.drizzle_filename = self.product_basename + "_" + self.filetype + ".fits"
        self.headerlet_filename = self.product_basename + "_hlet.fits"
        self.trl_filename = self.product_basename + "_trl.log"

        self.regions_dict = {}

    def wcs_drizzle_product(self, meta_wcs, configobj):
        """
            Create the drizzle-combined exposure image using the meta_wcs as the reference output
        """
        # AstroDrizzle will soon take a meta_wcs object which contains outnx, outny
        astrodrizzle.AstroDrizzle(input=self.full_filename,
                                  output=self.drizzle_filename,
                                  final_refimage=meta_wcs,
                                  final_outnx=meta_wcs.pixel_shape[1],
                                  final_outny=meta_wcs.pixel_shape[0],
                                  runfile=self.trl_filename,
                                  configobj=configobj)

        # Rename Astrodrizzle log file as a trailer file
        shutil.move(self.trl_filename, self.trl_filename.replace('.log', '.txt'))

    def image_drizzle_product(self, ref_image, total_shape, configobj):
        """
            Create the drizzle-combined exposure image using the reference image as the reference output
        """
        astrodrizzle.AstroDrizzle(input=self.full_filename,
                                  output=self.drizzle_filename,
                                  final_refimage=ref_image,
                                  final_outnx=total_shape[1],
                                  final_outny=total_shape[0],
                                  runfile=self.trl_filename,
                                  configobj=configobj)

        # Rename Astrodrizzle log file as a trailer file
        shutil.move(self.trl_filename, self.trl_filename.replace('.log', '.txt'))

        # log.info("Filter combined image... {} {}".format(self.drizzle_filename, self.full_filename))
