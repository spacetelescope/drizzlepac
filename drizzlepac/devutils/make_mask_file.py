#!/usr/bin/env python
from astropy.io import fits
import numpy as np
import pdb

in_imgname = "/Users/dulude/Documents/HLAtransition/runhlaprocessing_testing/acs_10595_06_flag_testing/hst_10595_06_acs_wfc_f435w_wht.fits"
out_imgname = in_imgname.replace("wht","msk")

in_hdulist = fits.open(in_imgname)
print(in_hdulist[0].data.dtype.name)
out_hdulist = in_hdulist.copy()
out_hdulist[0].data = np.where(in_hdulist[0].data == 0, 1, 0)
# out_hdulist[0].data.dtype = np.int16

out_hdulist.writeto(out_imgname)
# create frame shaped like in_image filled with 1s