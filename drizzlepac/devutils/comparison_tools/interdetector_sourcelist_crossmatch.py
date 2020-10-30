#!/usr/bin/env python


"""perform interdetector sourcelist crossmatch. This can be either HAP vs. HAP or HAP vs. HLA or
HLA vs. HLA."""

import argparse
import os
import pdb
import sys

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
import numpy as np

from drizzlepac.devutils.comparison_tools import compare_sourcelists
from drizzlepac.haputils import hla_flag_filter
import drizzlepac.haputils.comparison_utils as cu
import drizzlepac.haputils.svm_quality_analysis as svmqa
from stsci.stimage import xyxymatch
from stsci.tools import logutil

__taskname__ = 'interdetector_sourcelist_crossmatch'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)

# =======================================================================================================================
def run(sl_names, img_names, diagnostic_mode=False, log_level=logutil.logging.INFO):
    """main running subroutine

    Parameters
    ----------
    sl_names : list
        A list containing the reference sourcelist filename and the comparison sourcelist filename, in that
        order.

    img_names : list
        A list containing the reference image filename and the comparison image filename, in that order.

    diagnostic_mode : Bool, optional
        If this option is set to Boolean 'True', region files will be created to test the quality of the
        coordinate transformation. Default value is Boolean 'False'.

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the
        .log file. Default value is 'INFO'.

    Returns
    -------
    Nothing.

    """

    log.setLevel(log_level)
    # 0: correct ra, dec, x, y in HLA sourcelists
    if sl_names[1].endswith("phot.txt"):
        if not os.path.exists(sl_names[1].replace("phot.txt", "phot_corrected.txt")):
            corrected_hla_slname = svmqa.correct_hla_classic_ra_dec(sl_names[1], img_names[0],
                                                                    sl_names[1].split("_")[-1][:-8], log_level=log_level)
        sl_names[1] = corrected_hla_slname

    # 1: get sourcelist data from files
    ref_data, comp_data = cu.slFiles2dataTables(sl_names)

    # 2: stack up RA and DEC columns
    ref_ra_dec_values = np.stack((ref_data['RA'], ref_data['DEC']), axis=1)
    comp_ra_dec_values = np.stack((comp_data['RA'], comp_data['DEC']), axis=1)

    # 3: transform comp frame RA, DEC into Ref frame X, Y values and ref RA, DEC into comp X, Y values
    ref_xy_in_comp_frame = hla_flag_filter.rdtoxy(ref_ra_dec_values, img_names[1], "[1]", origin=0)
    comp_xy_in_ref_frame = hla_flag_filter.rdtoxy(comp_ra_dec_values, img_names[0], "[1]", origin=0)

    # 4: (diagnostic only) write out region files to check coordinate transformation
    if diagnostic_mode:
        for data_table, reg_filename in zip([ref_data, comp_data, comp_xy_in_ref_frame], ["ref_orig.reg", "comp_orig.reg", "comp_xform.reg"]):
            write_region_file(data_table, reg_filename)

    new_comp_xy = np.stack((comp_xy_in_ref_frame[:, 0], comp_xy_in_ref_frame[:, 1]), axis=1)
    ref_xy = np.stack((ref_data['X'], ref_data['Y']), axis=1)
    comp_xy = np.stack((comp_data['X'], comp_data['Y']), axis=1)

    matches = xyxymatch(new_comp_xy, ref_xy, tolerance=5.0, separation=1.0)
    # Report number and percentage of the total number of detected ref and comp sources that were matched
    log.info("Sourcelist Matching Results")
    log.info(
        "Reference sourcelist:  {} of {} total sources matched ({} %)".format(len(matches),
                                                                              len(ref_xy),
                                                                              100.0 * (float(len(matches)) / float(len(ref_xy)))))
    log.info(
        "Comparison sourcelist: {} of {} total sources matched ({} %)".format(len(matches),
                                                                              len(new_comp_xy),
                                                                              100.0 * (float(len(matches)) / float(len(new_comp_xy)))))
    # extract indices of the matching ref and comp lines and then use them to compile lists of matched X, Y,
    # RA and DEC for calculation of differences
    matched_lines_comp = []
    matched_lines_ref = []
    for item in matches:
        matched_lines_comp.append(item[2])
        matched_lines_ref.append(item[5])

    matching_values_ref_x = ref_data['X'][[matched_lines_ref]]
    matching_values_ref_y = ref_data['Y'][[matched_lines_ref]]
    matching_values_ref_ra = ref_data['RA'][[matched_lines_ref]]
    matching_values_ref_dec = ref_data['DEC'][[matched_lines_ref]]
    matching_values_comp_x = new_comp_xy[:, 0][[matched_lines_comp]]
    matching_values_comp_y = new_comp_xy[:, 1][[matched_lines_comp]]
    matching_values_comp_ra = comp_ra_dec_values[:, 0][[matched_lines_comp]]
    matching_values_comp_dec = comp_ra_dec_values[:, 1][[matched_lines_comp]]


    # get coordinate system type from fits headers
    ref_frame = fits.getval(img_names[0], "radesys", ext=('sci', 1)).lower()
    comp_frame = fits.getval(img_names[1], "radesys", ext=('sci', 1)).lower()

    # force RA and Dec values to be the correct type for SkyCoord() call
    if str(type(matching_values_ref_ra)) == "<class 'astropy.table.column.Column'>":
        matching_values_ref_ra = matching_values_ref_ra.tolist()
    if str(type(matching_values_ref_dec)) == "<class 'astropy.table.column.Column'>":
        matching_values_ref_dec = matching_values_ref_dec.tolist()
    if str(type(matching_values_comp_ra)) == "<class 'astropy.table.column.Column'>":
        matching_values_comp_ra = matching_values_comp_ra.tolist()
    if str(type(matching_values_comp_dec)) == "<class 'astropy.table.column.Column'>":
        matching_values_comp_dec = matching_values_comp_dec.tolist()

    # convert reference and comparison RA/Dec values into SkyCoord objects
    matching_values_ref_rd = SkyCoord(matching_values_ref_ra, matching_values_ref_dec, frame=ref_frame, unit="deg")
    matching_values_comp_rd = SkyCoord(matching_values_comp_ra, matching_values_comp_dec, frame=comp_frame, unit="deg")
    # convert to ICRS coord system
    if ref_frame != "icrs":
        matching_values_ref_rd = matching_values_ref_rd.icrs
    if comp_frame != "icrs":
        matching_values_comp_rd = matching_values_comp_rd.icrs

    # compute mean-subtracted differences
    diff_x = (matching_values_comp_x - matching_values_ref_x)
    diff_x -= sigma_clipped_stats(diff_x, sigma=3, maxiters=3)[0]
    diff_y = (matching_values_comp_y - matching_values_ref_y)
    diff_y -= sigma_clipped_stats(diff_y, sigma=3, maxiters=3)[0]
    diff_xy = np.sqrt(diff_x**2 + diff_y**2)
    diff_rd = matching_values_comp_rd.separation(matching_values_ref_rd).arcsec


    diff_list = [diff_x, diff_y, diff_xy, diff_rd]
    title_list = ["X-axis differences", "Y-axis differences", "Seperation", "On-sky separation"]
    units_list = ["HAP WFC3/UVIS pixels", "HAP WFC3/UVIS pixels", "HAP WFC3/UVIS pixels", "Arcseconds"]
    for diff_ra, title, units in zip(diff_list, title_list, units_list):
        log.info("Comparison - reference {} statistics ({})".format(title, units))

        compute_stats(diff_ra, title)

    generate_sorted_region_file(diff_xy, ref_xy_in_comp_frame[matched_lines_ref], comp_xy[matched_lines_comp], ref_data['FLAGS'][matched_lines_ref], comp_data['FLAGS'][matched_lines_comp])
# =======================================================================================================================
def compute_stats(diff_ra, title):
    """Compute linear statistics on specified differences

    Parameters
    ----------
    diff_ra : numpy or astropy array
       1xn array containing difference values that will be used to compute stats.

    title : string
       title of the differences

    Returns
    -------
    Nothing!
    """
    # 'sigma' and 'iters' input values used for various np.sigma_clipped_stats() runs
    sigma = 3
    n_iters = 3
    clipped_stats = sigma_clipped_stats(diff_ra, sigma=sigma, maxiters=n_iters)
    sigma_percentages = []
    for sig_val in [1.0, 2.0, 3.0]:
        sigma_percentages.append((float(np.shape(np.where((diff_ra >= (clipped_stats[0] - sig_val * clipped_stats[2])) & (diff_ra <= (clipped_stats[0] + sig_val * clipped_stats[2]))))[1])/float(np.shape(diff_ra)[0])) * 100.0)
    log.info("            Non-Clipped Statistics")
    log.info("Non-clipped minimum........................... {}".format(np.min(diff_ra)))
    log.info("Non-clipped maximum........................... {}".format(np.max(diff_ra)))
    log.info("Non-clipped mean.............................. {}".format(np.mean(diff_ra)))
    log.info("Non-clipped median............................ {}".format(np.median(diff_ra)))
    log.info("Non-clipped standard deviation................ {}".format(np.std(diff_ra)))
    log.info(
        "Non-clipped mean in units of SD............... {}\n".format(np.divide(np.mean(diff_ra), np.std(diff_ra))))
    log.info(
        "       Sigma-clipped Statistics; \u03C3 = {}, Number of clipping steps = {}".format(sigma,
                                                                                             n_iters))
    log.info("Sigma-clipped mean............................ {}".format(clipped_stats[0]))
    log.info("Sigma-clipped median.......................... {}".format(clipped_stats[1]))
    log.info("Sigma-clipped standard deviation.............. {}".format(clipped_stats[2]))
    for sig_val, pct_val in zip([1.0, 2.0, 3.0], sigma_percentages):
        log.info(
            "% all diff values within {}\u03C3 of clipped mean... {}%".format(int(sig_val), pct_val))
    log.info("\n\n")
    
# =======================================================================================================================
def generate_sorted_region_file(diff_ra, ref_xy, comp_xy, ref_flags, comp_flags):
    # #subtact off 3x3 sigma-clipped mean to eliminate any large-scale systemic offsets
    # sigma = 3
    # n_iters = 3
    # clipped_stats = sigma_clipped_stats(diff_ra, sigma=sigma, maxiters=n_iters)
    # diff_ra_meansub = diff_ra - clipped_stats[0]
    #
    # # get indicies of above array sorted by absolute value
    # # actual array still maintains sign (i.e. positive or negitve value)
    # sorted_idx = np.argsort(abs(diff_ra_meansub))[::-1]
    sorted_idx = np.argsort(abs(diff_ra))[::-1]
    region_filename = "testout.reg"
    f = open(region_filename, "w")
    ctr = 1
    cutoff = 50
    for idx in sorted_idx:
        print("{} {}      {} {}     {} {}      {} {}".format(ctr,diff_ra[idx], ref_xy[idx][0], ref_xy[idx][1], comp_xy[idx][0], comp_xy[idx][1], ref_flags[idx], comp_flags[idx]))
        f.write("line {} {} {} {} #line=0 1  text={{{} {}}}\n".format(ref_xy[idx][0], ref_xy[idx][1], comp_xy[idx][0], comp_xy[idx][1], ctr, diff_ra[idx]))

        if ctr == cutoff:
            break
        ctr +=1
    # ctr = 1
    # for diff_value, ref_line, comp_line in zip(diff_ra[sorted_idx], ref_xy[sorted_idx], comp_xy[sorted_idx]):
    #     print("{} {}      {} {}     {} {}".format(ctr,diff_value, ref_line[0], ref_line[1], comp_line[0], comp_line[1]))
    #     f.write("line {} {} {} {} #line=0 1\n".format(ref_line[0], ref_line[1], comp_line[0], comp_line[1]))
    #     if ctr == cutoff:
    #         break
    #     ctr +=1
    f.close()
    log.info("Wrote ds9 region file {}".format(region_filename))

# =======================================================================================================================

def write_region_file(coords, region_filename):
    """Writes specified x, y coordinates to specified region file.

    Parameters
    ----------
    coords : ndarray
        list of X, Y coordinates to write

    region_filename : str
        Name of the output region file

    Returns
    -------
    Nothing!
    """
    f = open(region_filename, "w")
    for line in coords:
        f.write("{} {}\n".format(line[0], line[1]))
    f.close()
    log.info("Wrote ds9 region file {}".format(region_filename))


# =======================================================================================================================
if __name__ == "__main__":
    # process command-line inputs with argparse
    parser = argparse.ArgumentParser(description='convert the X, Y coords in the comparison sourcelist into '
                                                 'the frame of reference of the reference sourcelist and run '
                                                 'compare_sourcelists.py')
    parser.add_argument('sl_list', nargs=2,
                        help='A space-separated pair of sourcelists to compare. The first sourcelist is '
                             'assumed to be the reference sourcelist that the second is being compared to.')
    parser.add_argument('-i', '--img_list', nargs=2, required=True,
                        help='A space-seperated list of containing the reference and comparison images '
                             'that correspond to reference and comparison sourcelists')
    parser.add_argument('-d', '--diagnostic_mode', required=False, action='store_true',
                        help='If this option is turned on, region files will be created to test the quality '
                             'of the coordinate transformation')
    parser.add_argument('-l', '--log_level', required=False, default='info',
                        choices=['critical', 'error', 'warning', 'info', 'debug'],
                        help='The desired level of verboseness in the log statements displayed on the screen '
                             'and written to the .log file. The level of verboseness from left to right, and '
                             'includes all log statements with a log_level left of the specified level. '
                             'Specifying "critical" will only record/display "critical" log statements, and '
                             'specifying "error" will record/display both "error" and "critical" log '
                             'statements, and so on.')
    user_args = parser.parse_args()

    # set up logging
    log_dict = {"critical": logutil.logging.CRITICAL,
                "error": logutil.logging.ERROR,
                "warning": logutil.logging.WARNING,
                "info": logutil.logging.INFO,
                "debug": logutil.logging.DEBUG}
    log_level = log_dict[user_args.log_level]
    log.setLevel(log_level)


    run(user_args.sl_list, user_args.img_list, user_args.diagnostic_mode, log_level=log_level)
