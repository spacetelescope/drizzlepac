"""Code that evaluates the quality of the SVM products generated by the drizzlepac package.

The JSON files generated here can be converted directly into a Pandas DataFrame
using the syntax:

>>> import json
>>> import pandas as pd
>>> with open("<rootname>_astrometry_resids.json") as jfile:
>>>     resids = json.load(jfile)
>>> pdtab = pd.DataFrame(resids)

These DataFrames can then be concatenated using:

>>> allpd = pdtab.concat([pdtab2, pdtab3])

where 'pdtab2' and 'pdtab3' are DataFrames generated from other datasets.  For
more information on how to merge DataFrames, see

https://pandas.pydata.org/pandas-docs/stable/user_guide/merging.html

Visualization of these Pandas DataFrames with Bokeh can follow the example
from:

https://programminghistorian.org/en/lessons/visualizing-with-bokeh

"""

# Standard library imports
import collections
import json
import os
import pdb
import sys

# external imports
import numpy as np
from scipy.spatial import KDTree

# Local application imports
from drizzlepac.hlautils import astrometric_utils as au
import drizzlepac.hlautils.diagnostic_utils as du
from stsci.tools import logutil
from stwcs.wcsutil import HSTWCS

__taskname__ = 'svm_quality_analysis'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)
# ----------------------------------------------------------------------------------------------------------------------

def characterize_gaia_distribution(hap_obj, log_level=logutil.logging.NOTSET):
    """Statistically describe distribution of GAIA sources in footprint.

    Computes and writes the file to a json file:

    - Number of GAIA sources
    - X centroid location
    - Y centroid location
    - X offset of centroid from image center
    - Y offset of centroid from image center
    - X standard deviation
    - Y standard deviation
    - minimum closest neighbor distance
    - maximum closest neighbor distance
    - mean closest neighbor distance
    - standard deviation of closest neighbor distances

    Parameters
    ----------
    hap_obj : drizzlepac.hlautils.Product.FilterProduct
        hap product object to process

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the .log file.
        Default value is 'NOTSET'.

    Returns
    -------
    Nothing
    """
    log.setLevel(log_level)

    # get table of GAIA sources in footprint
    gaia_table = generate_gaia_catalog(hap_obj, columns_to_remove=['mag', 'objID', 'GaiaID'])

    # if log_level is either 'DEBUG' or 'NOTSET', write out GAIA sources to DS9 region file
    if log_level <= 10:
        reg_file = "gaia_sources.reg"
        gaia_table.write(reg_file, format='ascii.csv')
        log.debug("Wrote GAIA source RA and Dec positions to DS9 region file '{}'".format(reg_file))

    # convert RA, Dec to image X, Y
    outwcs = HSTWCS(hap_obj.drizzle_filename + "[1]")
    x, y = outwcs.all_world2pix(gaia_table['RA'], gaia_table['DEC'], 1)

    # compute stats for the distribution
    centroid = [np.mean(x), np.mean(y)]
    centroid_offset = []
    for idx in range(0, 2):
        centroid_offset.append(outwcs.wcs.crpix[idx] - centroid[idx])
    std_dev = [np.std(x), np.std(y)]

    # Find straight-line distance to the closest neighbor for each GAIA source
    xys = np.array([x, y])
    xys = xys.reshape(len(x), 2)
    tree = KDTree(xys)
    neighborhood = tree.query(xys, 2)
    min_seps = np.empty([0])
    for sep_pair in neighborhood[0]:
        min_seps = np.append(min_seps, sep_pair[1])

    # add statistics to out_dict
    out_dict = collections.OrderedDict()
    out_dict["units"] = "pixels"
    out_dict["Number of GAIA sources"] = len(gaia_table)
    axis_list = ["X", "Y"]
    title_list = ["centroid", "offset of centroid from image center", "standard deviation"]
    for item_value, item_title in zip([centroid, centroid_offset, std_dev], title_list):
        for axis_item in enumerate(axis_list):
            log.info("{} {} ({}): {}".format(axis_item[1], item_title, out_dict["units"], item_value[axis_item[0]]))
            out_dict["{} {}".format(axis_item[1], item_title)] = item_value[axis_item[0]]
    min_sep_stats = [min_seps.min(), min_seps.max(), min_seps.mean(), min_seps.std()]
    min_sep_title_list = ["minimum closest neighbor distance",
                          "maximum closest neighbor distance",
                          "mean closest neighbor distance",
                          "standard deviation of closest neighbor distances"]
    for item_value, item_title in zip(min_sep_stats, min_sep_title_list):
        log.info("{} ({}): {}".format(item_title, out_dict["units"], item_value))
        out_dict[item_title] = item_value

    # write catalog to HapDiagnostic-formatted .json file.
    diag_obj = du.HapDiagnostic(log_level=log_level)
    diag_obj.instantiate_from_hap_obj(hap_obj,
                                      data_source="{}.characterize_gaia_distribution".format(__taskname__),
                                      description="A statistical characterization of the distribution of GAIA sources in image footprint")
    diag_obj.add_data_item(out_dict, "distribution characterization statistics")
    diag_obj.write_json_file(hap_obj.drizzle_filename[:-9] + "_svm_gaia_distribution_characterization.json", clobber=True)


# ----------------------------------------------------------------------------------------------------------------------

def compare_num_sources(catalog_list, drizzle_list, log_level=logutil.logging.NOTSET):
    """Determine the number of viable sources actually listed in SVM output catalogs.

    Parameters
    ----------
    catalog_list: list of strings
        Set of files on which to actually perform comparison.  Catalogs, Point and
        Segment, are generated for all of the Total data products in a single visit.
        The catalogs are detector-dependent.

    drizzle_list: list of strings
        Drizzle files for tht Total products which were mined to generate the output catalogs.

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the .log file.
        Default value is 'NOTSET'.

    Returns
    --------
    results : string
        Name of JSON file containing all the extracted results from the comparisons
        being performed.
    """
    log.setLevel(log_level)

    pnt_suffix = '_point-cat.ecsv'
    seg_suffix = '_segment-cat.ecsv'

    # Generate a separate JSON file for each detector
    # Drizzle filename example: hst_11665_06_wfc3_ir_total_ib4606_drz.fits
    # The filename is all lower-case by design.
    for drizzle_file in drizzle_list:
        tokens = drizzle_file.split('_')
        detector = tokens[4]
        ipppss = tokens[6]

        sources_dict = {'detector': detector, 'point': 0, 'segment': 0}

        # Construct the output JSON filename
        json_filename = ipppss + '_' + detector + '_svm_num_sources.json'

        # Construct catalog names for catalogs that should have been produced
        prefix = '_'.join(tokens[0:-1])
        cat_names = [prefix + pnt_suffix, prefix + seg_suffix]

        # If the catalog were actually produced, get the number of sources.
        # A catalog may not be produced because it was not requested, or there
        # was an error.  However, for the purposes of this program, it is OK
        # that no catalog was produced.
        for catalog in cat_names:
            does_exist = any(catalog in x for x in catalog_list)

            # if the catalog exists, open it and find the number of sources string
            num_sources = -1
            cat_type = ""
            if does_exist:
                file = open(catalog, 'r')
                for line in file:
                    sline = line.strip()

                    # All the comments are grouped at the start of the file. When
                    # the first non-comment line is found, there is no need to look further.
                    if not sline.startswith('#'):
                        log.info("Number of sources not reported in Catalog: {}.".format(catalog))
                        break

                    # When the matching comment line is found, get the value.
                    if sline.find('Number of sources') != -1:
                        num_sources = sline.split(' ')[-1][0:-1]
                        log.info("Catalog: {} Number of sources: {}.".format(catalog, num_sources))
                        break

                cat_type = 'point' if catalog.find("point") != -1 else 'segment'
                sources_dict[cat_type] = int(num_sources)

        # Set up the diagnostic object and write out the results
        diagnostic_obj = du.HapDiagnostic()
        diagnostic_obj.instantiate_from_fitsfile(drizzle_file,
                                                 data_source="{}.compare_num_sources".format(__taskname__),
                                                 description="Number of sources in Point and Segment catalogs")
        diagnostic_obj.add_data_item(sources_dict, 'number_of_sources')
        diagnostic_obj.write_json_file(json_filename)
        log.info("Generated quality statistics (number of sources) as {}.".format(json_filename))

        # Clean up
        del diagnostic_obj

# ----------------------------------------------------------------------------------------------------------------------

def find_gaia_sources(hap_obj, log_level=logutil.logging.NOTSET):
    """Creates a catalog of all GAIA sources in the footprint of a specified HAP final product image, and
    stores the GAIA object catalog as a hap diagnostic json file. The catalog contains RA, Dec and magnitude
    of each identified source. The catalog is sorted in decending order by brightness.

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
    gaia_table = generate_gaia_catalog(hap_obj, columns_to_remove=['objID', 'GaiaID'])
    # write catalog to HapDiagnostic-formatted .json file.
    diag_obj = du.HapDiagnostic(log_level=log_level)
    diag_obj.instantiate_from_hap_obj(hap_obj,
                                      data_source="{}.find_gaia_sources".format(__taskname__),
                                      description="A table of GAIA sources in image footprint")
    diag_obj.add_data_item(gaia_table, "GAIA sources")  # write catalog of identified GAIA sources
    diag_obj.add_data_item(len(gaia_table), "Number of GAIA sources")  # write the number of identified GAIA sources
    diag_obj.write_json_file(hap_obj.drizzle_filename[:-9]+"_svm_gaia_sources.json", clobber=True)

    # Clean up
    del diag_obj
    del gaia_table

# ----------------------------------------------------------------------------------------------------------------------

def generate_gaia_catalog(hap_obj, columns_to_remove = None):
    """Uses astrometric_utils.create_astrometric_catalog() to create a catalog of all GAIA sources in the
    image footprint. This catalog contains right ascension, declination, and magnitude values, and is sorted
    in descending order by brightness.

    Parameters
    ----------
    hap_obj : drizzlepac.hlautils.Product.TotalProduct, drizzlepac.hlautils.Product.FilterProduct, or
        drizzlepac.hlautils.Product.ExposureProduct, depending on input.
        hap product object to process

    Returns
    -------
    gaia_table : astropy table
        table containing right ascension, declination, and magnitude of all GAIA sources identified in the
        image footprint, sorted in descending order by brightness.
    """
    # Gather list of input flc/flt images
    img_list = []
    log.debug("GAIA catalog will be created using the following input images:")
    if hasattr(hap_obj, "edp_list"):  # for total and filter product objects
        for edp_item in hap_obj.edp_list:
            parse_info = edp_item.info.split("_")
            imgname = "{}_{}".format(parse_info[4], parse_info[5])
            log.debug(imgname)
            img_list.append(imgname)
    else:  # For single-exposure product objects
        parse_info = hap_obj.info.split("_")
        imgname = "{}_{}".format(parse_info[4], parse_info[5])
        log.debug(imgname)
        img_list.append(imgname)

    # generate catalog of GAIA sources
    gaia_table = au.create_astrometric_catalog(img_list, gaia_only =True)

    if columns_to_remove:
        gaia_table.remove_columns(columns_to_remove)
    # outwcs = HSTWCS(hap_obj.drizzle_filename,ext=1)
    # outwcs = au.build_self_reference(img_list[0] , clean_wcs=True)    # #
    # x, y = outwcs.all_world2pix(gaia_table['RA'], gaia_table['DEC'], 1)
    # foo,bar = au.within_footprint(hap_obj.drizzle_filename,outwcs,x,y)
    # # pdb.set_trace()
    if len(gaia_table) == 0:
        log.warning("No GAIA sources were found!")
    elif len(gaia_table) == 1:
        log.info("1 GAIA source was found.")
    else:
        log.info("{} GAIA sources were found.".format(len(gaia_table)))
    return gaia_table

# ============================================================================================================
if __name__ == "__main__":
    # Testing
    import pickle

    pfile = sys.argv[1]
    filehandler = open(pfile, 'rb')
    total_obj_list = pickle.load(filehandler)

    log_level = logutil.logging.DEBUG

    # Test compare_num_sources
    if False:
        total_catalog_list = []
        total_drizzle_list = []
        for total_obj in total_obj_list:
            total_drizzle_list.append(total_obj.drizzle_filename)
            total_catalog_list.append(total_obj.point_cat_filename)
            total_catalog_list.append(total_obj.segment_cat_filename)
        compare_num_sources(total_catalog_list, total_drizzle_list, log_level=log_level)

    # test find_gaia_sources
    if False:
        for total_obj in total_obj_list:
            find_gaia_sources(total_obj, log_level=log_level)
            for filter_obj in total_obj.fdp_list:
                find_gaia_sources(filter_obj, log_level=log_level)
                for exp_obj in filter_obj.edp_list:
                    find_gaia_sources(exp_obj, log_level=log_level)

    # test characterize_gaia_distribution
    if True:
        for total_obj in total_obj_list:
            for filter_obj in total_obj.fdp_list:
                characterize_gaia_distribution(filter_obj, log_level=log_level)