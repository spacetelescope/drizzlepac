#!/usr/bin/env python


"""perform interdetector sourcelist crossmatch. This can be either HAP vs. HAP or HAP vs. HLA or
HLA vs. HLA."""

import argparse
import os
import pdb
import sys

from astropy.table import Table
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

    # 2: stack up RA and DEC columns from comp SL
    comp_ra_dec_values = np.stack((comp_data['RA'], comp_data['DEC']), axis=1)

    # 3: transform comp frame RA, DEC into Ref frame X, Y values
    comp_xy_in_ref_frame = hla_flag_filter.rdtoxy(comp_ra_dec_values, img_names[0], "[1]", origin=0)

    # 4: (diagnostic only) write out region files to check coordinate transformation
    if diagnostic_mode:
        for data_table, reg_filename in zip([ref_data, comp_data, comp_xy_in_ref_frame], ["ref_orig.reg", "comp_orig.reg", "comp_xform.reg"]):
            f = open(reg_filename, "w")
            for line in data_table:
                f.write("{} {}\n".format(line[0], line[1]))
            f.close()
            log.info("Wrote ds9 region file {}".format(reg_filename))



    new_comp_xy = np.stack((comp_xy_in_ref_frame[:, 0], comp_xy_in_ref_frame[:, 1]), axis=1)
    ref_xy = np.stack((ref_data['X'], ref_data['Y']), axis=1)

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
    matched_lines_comp = []
    matched_lines_ref = []
    for item in matches:
        matched_lines_comp.append(item[2])
        matched_lines_ref.append(item[5])
    pdb.set_trace()
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