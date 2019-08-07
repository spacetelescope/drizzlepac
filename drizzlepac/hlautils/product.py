# ia1s70jrq_flt.fits,11150,A1S,70,149.232269,F110W,IR,/ifs/archive/ops/hst/public/ia1s/ia1s70jrq/ia1s70jrq_flt.fits
# filename,prop_id,prog_id,obset_id,exptime,filters,detector,path
# Make sure the proposal_id is a 5-character string
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
    def __init__(self, prop_id, obset_id, instrument, detector, filename):
        self.prop_id = prop_id.zfill(5)
        self.obset_id = obset_id
        self.instrument = instrument
        self.detector = detector

        self.basename = "hst_" + "_".join(map(str, [prop_id, obset_id, instrument, detector]))

        # exposure_name is the ipppssoo or a portion thereof
        self.exposure_name = filename[0:8]
        self.mjdutc = None

        # class variables to be set later
        self.param_dict = {}

    def update_config(self, param_dict):
        pass

class TotalProduct(HAPProduct):
    """ A Total Detection Product is a 'white' light mosaic comprised of
        images acquired with one instrument, one detector, all filters, and all
        exposure times

        The "tdp" is short hand for TotalProduct.
    """
    def __init__(self, prop_id, obset_id, instrument, detector, filename):
        super().__init__(prop_id, obset_id, instrument, detector, filename)

        self.exposure_name = filename[0:6]

        self.product_basename = self.basename + "_total_" + self.exposure_name
        self.trl_filename = self.product_basename + "_trl.log"
        self.point_cat_filename = self.product_basename + "_point-cat.ecsv"
        self.segment_cat_filename = self.product_basename + "_segment-cat.ecsv"
        self.drizzle_filename = ""

        # How exactly is this used as it is set the number of times a TDP is created for
        # an obset
        # Manifest is created for the obset and instrument, regarding of the detector
        self.manifest_name = '_'.join([instrument, filename[1:4], obset_id, "manifest.txt"])

        # These attributes will be set later
        self.edp_list = []
        self.fdp_list = []
        self.regions_dict = {}
        self.meta_wcs = None

    def add_member(self, edp):
        self.edp_list.append(edp)

    def add_product(self, fdp):
        self.fdp_list.append(fdp)

    def create_drizzle_filename(self):
        self._create_drizzle_filename_helper(self.edp_list[0])
        return self.drizzle_filename

    def _create_drizzle_filename_helper(self, edp):
        self.drizzle_filename = self.product_basename + "_" + edp.filetype + ".fits"

    # There is a unique WCS for each TotalProduct product which is generated based upon
    # the merging of all the ExposureProductss which comprise the specific TotalProduct.
    def build_metawcs(self):
        """
        image_list = [element.full_filename for element in edp_list]
        meta_wcs = wcs_functions.make_mosaic_wcs(image_list)
        """
        pass

    def drizzle_product(self, meta_wcs, configobj):
        """
        Create the drizzle-combined total image using the meta_wcs as the reference output
        """
        edp_filenames = [element.full_filename for element in self.edp_list]
        astrodrizzle.AstroDrizzle(input=edp_filenames,
                                  output=self.drizzle_filename,
                                  final_refimage=meta_wcs,
                                  runfile=self.trl_filename,
                                  configobj=configobj)

        # Rename Astrodrizzle log file as a trailer file
        shutil.move(self.trl_filename, self.trl_filename.replace('.log', '.txt'))

        # log.info("Total combined image... {} {}".format(self.drizzle_filename, edp_filenames))


class FilterProduct(HAPProduct):
    """ A Filter Detection Product is a mosaic comprised of images acquired
        during a single visit with one instrument, one detector, a single filter,
        and all exposure times.

        The "fdp" is short hand for FilterProduct.
    """
    def __init__(self, prop_id, obset_id, instrument, detector, filename, filters):
        super().__init__(prop_id, obset_id, instrument, detector, filename)
        self.exposure_name = filename[0:6]
        self.filters = filters

        self.product_basename = self.basename + "_" + "_".join(map(str, [filters, self.exposure_name]))
        # Trailer names .txt or .log
        self.trl_filename = self.product_basename + "_trl.log"
        self.point_cat_filename = self.product_basename + "_point-cat.ecsv"
        self.segment_cat_filename = self.product_basename + "_segment-cat.ecsv"
        self.drizzle_filename = ""

        # These attributes will be set later
        self.edp_list = []
        self.regions_dict = {}

    def add_member(self, edp):
        self.edp_list.append(edp)

    def create_drizzle_filename(self):
        self._create_drizzle_filename_helper(self.edp_list[0])
        return self.drizzle_filename

    def _create_drizzle_filename_helper(self, edp):
        self.drizzle_filename = self.product_basename + "_" + edp.filetype + ".fits"

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
                    # print("exposure_filenames: {}".format(exposure_filenames))
                    # print("headerlet_filenames: {}".format(headerlet_filenames))

                align_table = alignimages.perform_align(exposure_filenames,
                                                        debug=False,
                                                        runfile="alignimages.log",
                                                        update_hdr_wcs=True,
                                                        headerlet_filenames=headerlet_filenames)

                # FIX - Update exposure_filenames here?
        except Exception:
            # *** FIX Not sure about the logging here
            # Report a problem with the alignment
            # log.info("EXCEPTION encountered in alignimages.\n")
            exc_type, exc_value, exc_tb = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)
            # log.info("No correction to absolute astrometric frame applied.\n")
            align_table = None

        # Return a table which contains data regarding the alignment, as well as the
        # list of the flt/flc exposures which were part of the alignment process
        # FIX - This does not account for individual exposures which might have been
        # excluded from alignment.
        return align_table, exposure_filenames

    def drizzle_product(self, meta_wcs, configobj):
        """
        Create the drizzle-combined filter image using the meta_wcs as the reference output
        """
        edp_filenames = [element.full_filename for element in self.edp_list]
        astrodrizzle.AstroDrizzle(input=edp_filenames,
                                  output=self.drizzle_filename,
                                  final_refimage=meta_wcs,
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
        super().__init__(prop_id, obset_id, instrument, detector, filename)

        self.full_filename = filename
        # self.exptime = exptime
        self.filters = filters
        self.filetype = filetype

        self.product_basename = self.basename + "_" + "_".join(map(str, [filters, self.exposure_name]))
        self.drizzle_filename = self.product_basename + "_" + self.filetype + ".fits"
        self.headerlet_filename = self.product_basename + "_hlet.fits"
        self.trl_filename = self.product_basename + "_trl.log"

        self.regions_dict = {}

    def drizzle_product(self, meta_wcs, configobj):
        """
            Create the drizzle-combined exposure image using the meta_wcs as the reference output
        """
        # AstroDrizzle will soon take a meta_wcs object which contains outnx, outny
        astrodrizzle.AstroDrizzle(input=self.full_filename,
                                  output=self.drizzle_filename,
                                  final_refimage=meta_wcs,
                                  runfile=self.trl_filename,
                                  configobj=configobj)

        # Rename Astrodrizzle log file as a trailer file
        shutil.move(self.trl_filename, self.trl_filename.replace('.log', '.txt'))

        # log.info("Filter combined image... {} {}".format(self.drizzle_filename, self.full_filename))
