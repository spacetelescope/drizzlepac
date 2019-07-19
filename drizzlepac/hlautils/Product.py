# ia1s70jrq_flt.fits,11150,A1S,70,149.232269,F110W,IR,/ifs/archive/ops/hst/public/ia1s/ia1s70jrq/ia1s70jrq_flt.fits
# filename,prop_id,prog_id,obset_id,exptime,filters,detector,path
# Make sure the proposal_id is a 5-character string
import sys
import traceback

from drizzlepac import wcs_functions
from drizzlepac import alignimages

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
    """
    def __init__(self, prop_id, obset_id, instrument, detector, filename):
        super().__init__(prop_id, obset_id, instrument, detector, filename)

        self.exposure_name = filename[0:6]

        self.product_basename = self.basename + "_total_" + self.exposure_name
        self.trl_filename = self.product_basename + "_trl.txt"
        self.point_cat_filename = self.product_basename + "_point-cat.ecsv"
        self.segment_cat_filename = self.product_basename + "_segment-cat.ecsv"
        self.drizzle_filename = ""
        # How exactly is this used as it is set the number of times a TDP is created for
        # an obset
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

    def create_product_name(self, edp):
        self.drizzle_filename = self.product_basename + "_" + edp.filetype + ".fits"

    # There is a unique WCS for each TotalProduct product which is generated based upon
    # the merging of all the ExposureProductss which comprise the specific TotalProduct.
    def build_metawcs(self):
        """
        image_list = [element.full_filename for element in edp_list]
        meta_wcs = wcs_functions.make_mosaic_wcs(image_list)
        """
        pass

    # make outnx, outny, and wcs attribute for TotalProduct


class FilterProduct(HAPProduct):
    """ A Filter Detection Product is a mosaic comprised of images acquired
        during a single visit with one instrument, one detector, a single filter,
        and all exposure times.
    """
    def __init__(self, prop_id, obset_id, instrument, detector, filename, filters):
        super().__init__(prop_id, obset_id, instrument, detector, filename)

        self.exposure_name = filename[0:6]
        self.filters = filters

        self.product_basename = self.basename + "_".join(map(str, [filters, self.exposure_name]))
        self.trl_filename = self.product_basename + "_trl.txt"
        self.point_cat_filename = self.product_basename + "_point-cat.ecsv"
        self.segment_cat_filename = self.product_basename + "_segment-cat.ecsv"
        self.drizzle_filename = ""

        # These attributes will be set later
        self.edp_list = []
        self.regions_dict = {}

    def add_member(self, edp):
        self.edp_list.append(edp)

    def create_drizzle_name(self, edp):
        self.drizzle_filename = self.product_basename + "_" + edp.filetype + ".fits"
        pass

    def align_to_GAIA(self):
        """Extract the flt/flc filenames from the exposure product list, as
           well as the corresponding headerlet filenames to use legacy alignment
           routine.
        """
        print("Here is filterProduct")
        exposure_filenames = []
        headerlet_filenames = {}
        align_table = None
        try:
            if self.edp_list:
                for edp in self.edp_list:
                    exposure_filenames.append(edp.full_filename)
                    headerlet_filenames[edp.full_filename] = edp.headerlet_filename
                    print("exposure_filenames: {}".format(exposure_filenames))
                    print("headerlet_filenames: {}".format(headerlet_filenames))
                align_table = alignimages.perform_align(exposure_filenames,
                                                        debug=False,
                                                        runfile="alignimages.log",
                                                        update_hdr_wcs=True,
                                                        headerlet_filenames=headerlet_filenames)
        except Exception:
            # *** FIX Not sure about the logging here
            # Report a problem with the alignment
            # log.info("EXCEPTION encountered in alignimages.\n")
            exc_type, exc_value, exc_tb = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)
            # log.info("No correction to absolute astrometric frame applied.\n")
            align_table = None

        return align_table

    def drizzle_fdp(self):
        # for ExposureProduct in edp_list:
        pass

class ExposureProduct(HAPProduct):
    """ An Exposure Product is an individual exposure/image (flt/flc).
    """
    def __init__(self, prop_id, obset_id, instrument, detector, filename, filters, filetype):
        super().__init__(prop_id, obset_id, instrument, detector, filename)

        # self.exptime = exptime
        self.filters = filters
        self.filetype = filetype
        self.full_filename = filename

        self.product_basename = self.basename + "_".join(map(str, [filters, self.exposure_name]))
        self.drizzle_filename = self.product_basename + "_" + self.filetype + ".fits"
        self.headerlet_filename = self.product_basename + "_hlet.fits"
        self.trl_filename = self.product_basename + "_trl.txt"

        self.regions_dict = {}
