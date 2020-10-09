#!/usr/bin/env python


"""perform interdetector sourcelist crossmatch. This can be either HAP vs. HAP or HAP vs. HLA or
HLA vs. HLA."""

import argparse
import os
import pdb
import sys

import numpy as np

from drizzlepac.devutils.comparison_tools import compare_sourcelists
from drizzlepac.haputils import hla_flag_filter
import drizzlepac.haputils.comparison_utils as cu
import drizzlepac.haputils.svm_quality_analysis as svmqa
from stsci.tools import logutil

__taskname__ = 'interdetector_sourcelist_crossmatch'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)

# =======================================================================================================================
def run(sl_names, img_names, log_level=logutil.logging.INFO):

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

    for data_table, reg_filename in zip([ref_data, comp_data, comp_xy_in_ref_frame], ["ref_orig.reg", "comp_orig.reg", "comp_xform.reg"]):
        f = open(reg_filename, "w")
        for line in data_table:
            f.write("{} {}\n".format(line[0], line[1]))
        f.close()
        log.info("Wrote ds9 region file {}".format(reg_filename))
# =======================================================================================================================
if __name__ == "__main__":
    # process command-line inputs with argparse
    parser = argparse.ArgumentParser(description='convert the X, Y coords in the comparision sourcelist into '
                                                 'the frame of reference of the reference sourcelist and run '
                                                 'compare_sourcelists.py')
    parser.add_argument('sourcelistNames', nargs=2,
                        help='A space-separated pair of sourcelists to compare. The first sourcelist is '
                             'assumed to be the reference sourcelist that the second is being compared to.')
    parser.add_argument('-i', '--imageName', nargs=2, required=True,
                        help='A space-seperated list of containing the reference and comparison images '
                             'that correspond to reference and comparison sourcelists')
    user_args = parser.parse_args()

    run(user_args.sourcelistNames, user_args.imageNames, log_level=logutil.logging.INFO)