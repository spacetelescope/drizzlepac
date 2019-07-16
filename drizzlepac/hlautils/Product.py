# ia1s70jrq_flt.fits,11150,A1S,70,149.232269,F110W,IR,/ifs/archive/ops/hst/public/ia1s/ia1s70jrq/ia1s70jrq_flt.fits
# filename,prop_id,prog_id,obset_id,exptime,filters,detector,path

# Define the mapping between the first character of the filename and the associated instrument
# INSTRUMENT_DICT = {"i": "WFC3", "j": "ACS", "o": "STIS", "u": "WFPC2", "x": "FOC", "w": "WFPC"}

class HAPProduct:
    def __init__(self, prop_id, obset_id, instrument, detector, filename):
        self.prop_id = prop_id.zfill(5)
        self.obset_id = obset_id
        self.instrument = instrument
        self.detector = detector
        # self.instrument = INSTRUMENT_DICT[filename[0])

        self.basename = "hst_" + "_".join(map(str, [prop_id, obset_id, instrument, detector])) + "_"

        # exposure_name is the ipppssoo or a portion thereof
        self.exposure_name = filename[0:8]
        self.mjdutc = None

        # class variables to be set later
        self.param_dict = {}

    def update_config(self, param_dict):
        pass

class TDP(HAPProduct):
    """ A Total Detection Product (TDP) is a 'white' light mosaic comprised of
        images acquired with one instrument, one detector, all filters, and all
        exposure times
    """
    def __init__(self, prop_id, obset_id, instrument, detector, filename):
        super().__init__(prop_id, obset_id, instrument, detector, filename)

        self.exposure_name = filename[0:6]

        self.product_basename = self.basename + "_total_" + self.exposure_name
        self.trl_filename = self.product_basename + "trl.txt"
        self.point_cat_filename = self.product_basename + "_point-cat.ecsv"
        self.segment_cat_filename = self.product_basename + "_segment-cat.ecsv"
        self.drizzle_filename = ""

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

    # There is a unique WCS for each TDP product which is generated based upon
    # the merging of all the EDPs which comprise the specific TDP.
    # Is the list the flt/flcs or the associated headerlets for generating the wcs?
    def build_metawcs(self):
        # image_list = [element.full_filename for element in edp_list]
        # meta_wcs = wcs_functions.make_mosaic_wcs(image_list)
        pass

    # make outnx, outny, and wcs attribute for TDP


class FDP(HAPProduct):
    """ A Filter Detection Product (FDP) is a mosaic comprised of images acquired
        during a single visit with one instrument, one detector, a single filter,
        and all exposure times.
    """
    def __init__(self, prop_id, obset_id, instrument, detector, filename, filters):
        super().__init__(prop_id, obset_id, instrument, detector, filename)

        self.exposure_name = filename[0:6]
        self.filters = filters

        self.product_basename = self.basename + "_".join(map(str, [filters, self.exposure_name]))
        self.trl_filename = self.product_basename + "trl.txt"
        self.point_cat_filename = self.product_basename + "_point-cat.ecsv"
        self.segment_cat_filename = self.product_basename + "_segment-cat.ecsv"
        self.drizzle_filename = ""

        # These attributes will be set later
        self.edp_list = []
        self.regions_dict = {}

    def add_member(self, edp):
        self.edp_list.append(edp)

    def create_drizzle_name(self, edp):
        # self.drizzle_filename = self.product_basename + "_" + edp.filetype + ".fits"
        pass

    def align_to_GAIA(self):
        pass

    def drizzle_fdp(self):
        # for EDP in edp_list:
        pass

class EDP(HAPProduct):
    def __init__(self, prop_id, obset_id, instrument, detector, filename, filters, filetype):
        super().__init__(prop_id, obset_id, instrument, detector, filename)

        # ipppssoot_xxx.fits
        self.full_filename = filename
        # self.exptime = exptime
        self.filters = filters
        self.filetype = filetype

        self.product_basename = self.basename + "_".join(map(str, [filters, self.exposure_name]))
        self.drizzle_filename = self.product_basename + "_" + self.filetype + ".fits"
        self.headerlet_filename = self.product_basename + "_hlet.fits"
        self.trl_filename = self.product_basename + "trl.txt"

        self.regions_dict = {}
