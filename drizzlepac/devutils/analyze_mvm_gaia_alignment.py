#!/usr/bin/env python

"""Quantify how well MVM products are aligned to GAIA sources found in the image footprint"""
### NOTE: BORROW HEAVILY FROM SVM_QUALILITY_ANALYSIS
import argparse
import pdb
import sys

import numpy as np
from photutils.detection import DAOStarFinder


# ============================================================================================================
def perform(imagename, gaia_list_name):
    # 1: generate WCS obj. for custom mosaic image
    # 2a: generate table of all gaia sources in frame
    # 2b: compute new WCS from stack of input flc.fits images to remove gaia sources outside image footprint
    # 2c: use WCS computed in 2B to remove remove gaia sources outside image footprint from table (see svm_quality_analysis.py, line #1114)
    # 3: Use step 1 WCS to compute detector X, Y coords from RA, DEC positions in trimmed gaia table (see svm_quality_analysis.py, line #1114)
    # 4: feed x, y coords into photutils.detection.daostarfinder() as initial guesses to get actual centroid positions of gaia sources
    # 5: convert daostarfinder output x, y centroid positions to RA, DEC using step 1 WCS info
    # 6: compute and report statistics based on X, Y and RA, DEC position residuals. Some of what's needed
    # here can be pulled from svm_quality_analysis.characterize_gaia_distribution() and also from compare_sourcelists() or comparision_utils.