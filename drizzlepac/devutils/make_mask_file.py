#!/usr/bin/env python
from astropy.io import fits
import numpy as np
import pdb

def make_mask_file_old(in_imgname):
    from scipy import ndimage
    maskfilename = in_imgname.replace("_drz.fits","_msk.fits")
    in_hdulist = fits.open(in_imgname)
    out_hdulist = fits.PrimaryHDU(in_hdulist[2].data)



    out_hdulist.data = np.where(out_hdulist.data == 0, 1, 0)
    ndimage.binary_fill_holes(out_hdulist.data)


    out_hdulist.writeto(maskfilename)
    # create frame shaped like in_image filled with 1s

def make_mask_file(drz_image,input_list):
    from drizzlepac.hlautils import cell_utils as cutils
    from stwcs.wcsutil import HSTWCS


    s = cutils.SkyFootprint(HSTWCS(drz_image, ext = 1))
    s.build(input_list)
    hdu = fits.PrimaryHDU(s.total_mask)
    hdul = fits.HDUList([hdu])
    maskfilename = drz_image.replace("_drc.fits","_msk.fits")
    hdul.writeto(maskfilename)
    hdul.close()
    print("Wrote {}".format(maskfilename))

if __name__ == "__main__":
    in_imgname = "hst_10265_01_acs_wfc_f606w_j92c01_drc.fits"
    in_imgname = "hst_10265_01_acs_wfc_total_drz.fits"
    input_list = ["j92c01b4q_flc.fits", "j92c01b5q_flc.fits", "j92c01b7q_flc.fits", "j92c01b9q_flc.fits"]
    make_mask_file(in_imgname,input_list)
    # make_mask_file_old("hst_10265_01_acs_wfc_total_drz.fits")
