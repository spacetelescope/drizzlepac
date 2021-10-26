#!/usr/bin/env python

"""Quantify how well MVM products are aligned to GAIA sources found in the image footprint"""

# Standard library imports
import argparse
import glob
import pdb
import sys

# Related third party imports
from astropy.io import fits
import numpy as np
from photutils.detection import DAOStarFinder

# Local application imports
from drizzlepac.haputils import astrometric_utils as amutils
import stwcs
# ============================================================================================================

def perform(mosaic_imgname, flcflt_list):
    """ Statistically quantify quality of GAIA MVM alignment

    Parameters
    ----------
    mosaic_imgname : str
        Name of the MVM-processed mosaic image to process

    flcflt_list : list
        lList of calibrated flc.fits and/or flt.fits images to process

    Returns
    -------
    Nothing!
    """
    # 0: read in flc/flt fits files from user-specifed fits file
    with open(flcflt_list, mode='r') as imgfile:
        imglist = imgfile.readlines()
    for x in range(0, len(imglist)): imglist[x] = imglist[x].strip()

    # 1: generate WCS obj. for custom mosaic image
    mosaic_wcs = stwcs.wcsutil.HSTWCS(mosaic_imgname, ext=1)

    # 2a: generate table of all gaia sources in frame
    gaia_tab = amutils.create_astrometric_catalog(imglist, existing_wcs=mosaic_wcs, catalog='GAIAedr3', output='test_gaia_edr3.ecsv', use_footprint=True)

    # 2b: compute new WCS from stack of input flc.fits images to remove gaia sources outside image footprint
    uncal_imglist = []
    for imgname in imglist:
        uncal_imglist.append(fits.getval(imgname, "FILENAME"))

    # 2c: use WCS computed in 2B to remove remove gaia sources outside image footprint from table (see svm_quality_analysis.py, line #1114)
    # 3: Use step 1 WCS to compute detector X, Y coords from RA, DEC positions in trimmed gaia table (see svm_quality_analysis.py, line #1114)
    # 4: feed x, y coords into photutils.detection.daostarfinder() as initial guesses to get actual centroid positions of gaia sources
    # 5: convert daostarfinder output x, y centroid positions to RA, DEC using step 1 WCS info
    # 6: compute and report statistics based on X, Y and RA, DEC position residuals. Some of what's needed
    # 7 here can be pulled from svm_quality_analysis.characterize_gaia_distribution() and also from compare_sourcelists() or comparision_utils.

    pdb.set_trace()
# ============================================================================================================

if __name__ == "__main__":
    # Parse command-line input args
    parser = argparse.ArgumentParser(description='Statistically quantify quality of GAIA MVM alignment')
    parser.add_argument('mosaic_imgname', help='Name of the MVM-processed mosaic image to process')
    parser.add_argument('flcflt_list', help='list of calibrated flc.fits and/or flt.fits images to process')
    input_args = parser.parse_args()

    # Perform analysis
    perform(input_args.mosaic_imgname, input_args.flcflt_list)


