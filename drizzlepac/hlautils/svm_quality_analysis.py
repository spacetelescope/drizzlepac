#!/usr/bin/env python

import pdb  # TODO: remove once everything is working
import pickle  # TODO: remove once everything is working
import sys


from drizzlepac.hlautils import astrometric_utils
from drizzlepac.hlautils import diagnostic_utils
from stsci.tools import logutil


__taskname__ = 'svm_quality_analysis'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)


def run_find_gaia_sources(hap_obj, log_level=logutil.logging.NOTSET):
    """Creates a catalog of all GAIA sources in the footprint of a specified HAP final product image, and
    stores the GAIA object catalog as a hap diagnostic json file.

    Parameters
    ----------
    hap_obj : drizzlepac.hlautils.Product.TotalProduct, drizzlepac.hlautils.Product.FilterProduct, or
        drizzlepac.hlautils.Product.ExposureProduct, depending on input.
        hap product object to process

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the .log file.
        Default value is 'NOTSET'.

    Returns
    -------
    Nothing.
    """
    log.setLevel(log_level)

    # Gather list of input flc/flt images
    img_list = []
    log.info("GAIA catalog will be created using the following input images:")
    if hasattr(hap_obj,"edp_list"):  # for total and filter product objects
        for edp_item in hap_obj.edp_list:
            parse_info = edp_item.info.split("_")
            imgname = "{}_{}".format(parse_info[4], parse_info[5])
            log.info(imgname)
            img_list.append(imgname)
    else:  # For single-exposure product objects
        parse_info = hap_obj.info.split("_")
        imgname = "{}_{}".format(parse_info[4], parse_info[5])
        log.info(imgname)
        img_list.append(imgname)

    # generate catalog of GAIA sources
    ref_table = astrometric_utils.create_astrometric_catalog(img_list)
    ref_table.remove_columns(['objID', 'GaiaID'])
    if len(ref_table) == 0:
        log.warning("No GAIA sources were found!")
    elif len(ref_table) == 1:
        log.info("1 GAIA source was found.")
    else:
        log.info("{} GAIA sources were found.".format(len(ref_table)))

    # write catalog to HapDiagnostic-formatted .json file.
    diag_obj = diagnostic_utils.HapDiagnostic(log_level=log_level)
    diag_obj.instantiate_from_hap_obj(hap_obj, data_source="run_find_gaia_sources", description="A table of GAIA sources in image footprint")
    diag_obj.add_data_item(ref_table, "GAIA sources")
    diag_obj.write_json_file(hap_obj.drizzle_filename+"_gaia_sources.json", clobber=True)



# ======================================================================================================================


if __name__ == "__main__":
    # Testing
    pfile = "total_obj_list_full.pickle"
    filehandler = open(pfile, 'rb')
    total_obj_list = pickle.load(filehandler)
    run_find_gaia_sources(total_obj_list[0].edp_list[0], log_level=logutil.logging.DEBUG)