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

# Non-standard library imports
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
import numpy as np

# Local application imports
from drizzlepac.hlautils import astrometric_utils
import drizzlepac.hlautils.diagnostic_utils as du
import drizzlepac.devutils.comparison_tools.compare_sourcelists as csl
from stsci.tools import logutil

__taskname__ = 'svm_quality_analysis'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)


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

# ------------------------------------------------------------------------------------------------------------

def compare_ra_dec_crossmatches(hap_obj, log_level=logutil.logging.NOTSET):
    """Compare the equatorial coordinates of cross-matches sources between the Point and Segment catalogs.
    The results .json file contains the following information:

        - image header information
        - cross-match details (input catalog lengths, number of cross-matched sources, coordinate system)
        - catalog containing RA and dec values of cross-matched point catalog sources
        - catalog containing RA and dec values of cross-matched segment catalog sources
        - Statistics describing the on-sky seperation of the cross-matched point and segment catalogs
        (non-clipped and sigma-clipped mean, median and standard deviation values)

    Parameters
    ----------
    hap_obj : drizzlepac.hlautils.Product.FilterProduct
        hap filter product object to process

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the .log file.
        Default value is 'NOTSET'.

    Returns
    --------
    nothing.
    """
    log.setLevel(log_level)
    slNames = [hap_obj.point_cat_filename,hap_obj.segment_cat_filename]
    imgNames = [hap_obj.drizzle_filename, hap_obj.drizzle_filename]
    good_flag_sum = 255 # all bits good

    diag_obj = du.HapDiagnostic(log_level=log_level)
    diag_obj.instantiate_from_hap_obj(hap_obj,
                                      data_source="{}.compare_ra_dec_crossmatches".format(__taskname__),
                                      description="matched point and segment catalog RA and Dec values")
    json_results_dict = collections.OrderedDict()
    # add reference and comparision catalog filenames as header elements
    json_results_dict["point catalog filename"] = slNames[0]
    json_results_dict["segment catalog filename"] = slNames[1]

    # 1: Read in sourcelists files into astropy table or 2-d array so that individual columns from each sourcelist can be easily accessed later in the code.
    point_data, seg_data = csl.slFiles2dataTables(slNames)
    log.info("Valid point data columns:   {}".format(list(point_data.keys())))
    log.info("Valid segment data columns: {}".format(list(seg_data.keys())))
    log.info("\n")
    log.info("Data columns to be compared:")
    columns_to_compare = list(set(point_data.keys()).intersection(set(seg_data.keys())))
    for listItem in sorted(columns_to_compare):
        log.info(listItem)
    log.info("\n")
    # 2: Run starmatch_hist to get list of matched sources common to both input sourcelists
    slLengths = [len(point_data['RA']), len(seg_data['RA'])]
    json_results_dict['point catalog length'] = slLengths[0]
    json_results_dict['segment catalog length'] = slLengths[1]
    matching_lines_ref, matching_lines_img = csl.getMatchedLists(slNames, imgNames, slLengths, log_level=log_level)
    json_results_dict['number of cross-matches'] = len(matching_lines_ref)
    if len(matching_lines_ref) == 0 or len(matching_lines_img) == 0:
        log.critical("*** Comparisons cannot be computed. No matching sources were found. ***")
        return ("ERROR")
    # 2: Create masks to remove missing values or values not considered "good" according to user-specified good bit values
    # 2a: create mask that identifies lines any value from any column is missing
    missing_mask = csl.mask_missing_values(point_data, seg_data, matching_lines_ref, matching_lines_img, columns_to_compare)
    # 2b: create mask based on flag values
    matched_values = csl.extractMatchedLines("FLAGS", point_data, seg_data, matching_lines_ref, matching_lines_img)
    bitmask = csl.make_flag_mask(matched_values, good_flag_sum, missing_mask)

    matched_values_ra = csl.extractMatchedLines("RA", point_data, seg_data, matching_lines_ref, matching_lines_img,
                                                bitmask=bitmask)
    matched_values_dec = csl.extractMatchedLines("DEC", point_data, seg_data, matching_lines_ref, matching_lines_img,
                                                 bitmask=bitmask)

    if len(matched_values_ra) > 0 and len(matched_values_ra) == len(matched_values_dec):
        # get coordinate system type from fits headers

        point_frame = fits.getval(imgNames[0], "radesys", ext=('sci', 1)).lower()
        seg_frame = fits.getval(imgNames[1], "radesys", ext=('sci', 1)).lower()
        # Add 'ref_frame' and 'comp_frame" values to header so that will SkyCoord() execute OK
        json_results_dict["point frame"] = point_frame
        json_results_dict["segment frame"] = seg_frame

        # convert reference and comparision RA/Dec values into SkyCoord objects
        matched_values_point = SkyCoord(matched_values_ra[0, :], matched_values_dec[0, :], frame=point_frame,
                                        unit="deg")
        matched_values_seg = SkyCoord(matched_values_ra[1, :], matched_values_dec[1, :], frame=seg_frame,
                                      unit="deg")
        # convert to ICRS coord system
        if point_frame != "icrs":
            matched_values_point = matched_values_point.icrs
        if seg_frame != "icrs":
            matched_values_seg = matched_values_seg.icrs

        # compute on-sky separations in arcseconds
        sep = matched_values_seg.separation(matched_values_point).arcsec

        # Compute and store statistics  on separations
        sep_stat_dict=collections.OrderedDict()
        sep_stat_dict["units"] = "arcseconds"
        sep_stat_dict["Non-clipped min"] = np.min(sep)
        sep_stat_dict["Non-clipped max"] = np.max(sep)
        sep_stat_dict["Non-clipped mean"] = np.mean(sep)
        sep_stat_dict["Non-clipped median"] = np.median(sep)
        sep_stat_dict["Non-clipped standard deviation"] = np.std(sep)
        sigma = 3
        maxiters = 3
        clippedStats = sigma_clipped_stats(sep, sigma=sigma, maxiters=maxiters)
        sep_stat_dict["{}x{} sigma-clipped mean".format(maxiters, sigma)] = clippedStats[0]
        sep_stat_dict["{}x{} sigma-clipped median".format(maxiters, sigma)] = clippedStats[1]
        sep_stat_dict["{}x{} sigma-clipped standard deviation".format(maxiters, sigma)] = clippedStats[2]

        # Create output catalogs for json file
        out_cat_point = Table([matched_values_ra[0], matched_values_dec[0]], names=("Right ascension", "Declination"))
        out_cat_seg = Table([matched_values_ra[1], matched_values_dec[1]], names=("Right ascension", "Declination"))
        for table_item in [out_cat_point,out_cat_seg]:
            for col_name in ["Right ascension", "Declination"]:
                table_item[col_name].unit = "degrees"  # Add correct units

        # add various data items to diag_obj
        diag_obj.add_data_item(json_results_dict, "Cross-match details")
        diag_obj.add_data_item(out_cat_point, "Cross-matched point catalog")
        diag_obj.add_data_item(out_cat_seg, "Cross-matched segment catalog")
        diag_obj.add_data_item(sep_stat_dict, "Segment - point on-sky separation statistics")

        # write everything out to the json file
        json_filename = hap_obj.drizzle_filename[:-9]+"_point_segment_crossmatch.json"
        diag_obj.write_json_file(json_filename, clobber=True)


# ------------------------------------------------------------------------------------------------------------

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
    ref_table = astrometric_utils.create_astrometric_catalog(img_list)
    ref_table.remove_columns(['objID', 'GaiaID'])
    if len(ref_table) == 0:
        log.warning("No GAIA sources were found!")
    elif len(ref_table) == 1:
        log.info("1 GAIA source was found.")
    else:
        log.info("{} GAIA sources were found.".format(len(ref_table)))

    # write catalog to HapDiagnostic-formatted .json file.
    diag_obj = du.HapDiagnostic(log_level=log_level)
    diag_obj.instantiate_from_hap_obj(hap_obj,
                                      data_source="{}.find_gaia_sources".format(__taskname__),
                                      description="A table of GAIA sources in image footprint")
    diag_obj.add_data_item(ref_table, "GAIA sources")  # write catalog of identified GAIA sources
    diag_obj.add_data_item(len(ref_table), "Number of GAIA sources")  # write the number of identified GAIA sources
    diag_obj.write_json_file(hap_obj.drizzle_filename+"_gaia_sources.json", clobber=True)

    # Clean up
    del diag_obj
    del ref_table

# ============================================================================================================
if __name__ == "__main__":
    # Testing
    import pickle

    pfile = sys.argv[1]
    filehandler = open(pfile, 'rb')
    total_obj_list = pickle.load(filehandler)

    log_level = logutil.logging.INFO

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

    # test compare_ra_dec_crossmatches
    if True:
        for total_obj in total_obj_list:
            for filter_obj in total_obj.fdp_list:
                compare_ra_dec_crossmatches(filter_obj, log_level=log_level)
