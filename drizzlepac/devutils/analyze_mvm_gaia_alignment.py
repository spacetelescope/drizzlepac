#!/usr/bin/env python

"""Quantify how well MVM products are aligned to GAIA sources found in the image footprint"""

# Standard library imports
import argparse
import glob
import pdb
import sys

# Related third party imports
from astropy.io import fits
from astropy.table import Column
import numpy as np
from photutils.detection import DAOStarFinder

# Local application imports
from drizzlepac import wcs_functions
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
    # 0: read in flc/flt fits files from user-specified fits file
    with open(flcflt_list, mode='r') as imgfile:
        imglist = imgfile.readlines()
    for x in range(0, len(imglist)): imglist[x] = imglist[x].strip()

    # 1: generate WCS obj. for custom mosaic image
    mosaic_wcs = stwcs.wcsutil.HSTWCS(mosaic_imgname, ext=1)

    # 2a: generate table of all gaia sources in frame
    gaia_table = amutils.create_astrometric_catalog(imglist, existing_wcs=mosaic_wcs, catalog='GAIAedr3', use_footprint=True)
    gaia_table.remove_columns(["mag", "objID"])  # remove "mag" and "objID" columns leaving just "RA" and "DEC"

    # 2b: Remove gaia sources outside footprint of input flc/flt images, add X and Y coord columns
    mosaic_hdu = fits.open(mosaic_imgname)
    drc_wht_array = np.zeros_like(mosaic_hdu["WHT"].data)
    drc_list = glob.glob("hst*{}*drc.fits".format(fits.getval(imglist[0], "FILENAME")[:6]))
    for wht_ctr, item in zip(range(0, len(drc_list)), drc_list):
        print("{}/{}: Adding weight image from {} to combined weight image".format(wht_ctr + 1, len(drc_list), item))
        drc_hdu = fits.open(item)
        drc_wht_array += drc_hdu["WHT"].data
        drc_hdu.close()
    x, y = mosaic_wcs.all_world2pix(gaia_table['RA'], gaia_table['DEC'], 1)
    x_col = Column(name="X", data=x, dtype=np.float64)
    y_col = Column(name="Y", data=y, dtype=np.float64)
    gaia_table.add_column(x_col, index=3)
    gaia_table.add_column(y_col, index=4)
    drc_wht_array = np.where(drc_wht_array == 0, np.nan, drc_wht_array)
    mask = amutils.within_footprint(drc_wht_array, mosaic_wcs, x, y)
    gaia_table.write("gaia_edr3_untrimmed.reg", format="ascii.ecsv", overwrite=True)  # TODO: DIAGNOSTIC LINE REMOVE PRIOR TO DEPLOYMENT
    gaia_table = gaia_table[mask]
    gaia_table.write("gaia_edr3_trimmed.reg", format="ascii.ecsv", overwrite=True)  # TODO: DIAGNOSTIC LINE REMOVE PRIOR TO DEPLOYMENT

    # 3: feed x, y coords into photutils.detection.daostarfinder() as initial guesses to get actual centroid positions of gaia sources
    # 4: convert daostarfinder output x, y centroid positions to RA, DEC using step 1 WCS info
    # 5: compute and report statistics based on X, Y and RA, DEC position residuals. Some of what's needed
    # 6 here can be pulled from svm_quality_analysis.characterize_gaia_distribution() and also from compare_sourcelists() or comparision_utils.

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


