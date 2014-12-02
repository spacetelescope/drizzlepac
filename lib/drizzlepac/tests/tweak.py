#!/usr/bin/env python

import drizzlepac
from drizzlepac import tweakreg

tweakreg.TweakReg(['output_blot_grid.fits','reference_blot_grid.fits'],
threshold=30.0, updatehdr=False, updatewcs=False, writecat=False,
clean=True, interactive=False, use2dhist=True, see2dplot=False,
fitgeometry='rscale', residplot='both')
