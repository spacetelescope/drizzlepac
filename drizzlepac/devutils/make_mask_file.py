#!/usr/bin/env python
from astropy.io import fits as pyfits
import numpy as np
import pdb

def make_mask_file_old(in_imgname):
    from scipy import ndimage
    out_imgname = in_imgname.replace("wht","msk")

    in_wht_hdulist = fits.open(in_imgname)
    in_sat_hdulist = fits.open(in_imgname.replace("wht","sat"))

    out_hdulist = in_wht_hdulist.copy()
    sat_data = in_sat_hdulist[0].data
    ndimage.binary_dilation(sat_data,iterations=11)
    out_hdulist[0].data = np.where(in_wht_hdulist[0].data+sat_data == 0, 1, 0)


    out_hdulist.writeto(out_imgname)

def make_mask_file(drz_image,input_list):
    from drizzlepac.hlautils import cell_utils as cutils
    from stwcs.wcsutil import HSTWCS

    s = cutils.SkyFootprint(HSTWCS(drz_image, ext = 1))
    s.build(input_list)
    hdu = pyfits.PrimaryHDU(s.total_mask)
    hdu.data = numpy.where(hdu.data  == 0, 1, 0)

    hdul = pyfits.HDUList([hdu])

    maskfilename = drz_image.replace(drz_image[-9:],"_msk.fits")
    hdul.writeto(maskfilename)
    hdul.close()
    print("Wrote {}".format(maskfilename))

def make_mask_file_actual(imfile):
    import scipy.ndimage
    import os
    mask = pyfits.open(imfile)[1].data != 0
    dilate = scipy.ndimage.morphology.binary_dilation
    erode = scipy.ndimage.morphology.binary_erosion
    kernel1 = np.ones((25, 25), dtype=int)
    kernel2 = np.ones((31, 31), dtype=int)
    # add padding around the edge so pixels close to image boundary are correct
    padding = 13
    bigmask = np.pad(mask, padding, 'constant')
    # strip the padding back off after creating mask
    mask = (erode(dilate(bigmask, kernel1), kernel2) == 0)[padding:-padding, padding:-padding]

    flagfile = imfile.replace(imfile[-9:],"_msk2.fits")
    if os.path.exists(flagfile):
        os.remove(flagfile)
    pyfits.writeto(flagfile, mask.astype(np.int16))



if __name__ == "__main__":
    in_imgname = "hst_10265_01_acs_wfc_f606w_j92c01_drc.fits"
    input_list = ["j92c01b4q_flc.fits", "j92c01b5q_flc.fits", "j92c01b7q_flc.fits", "j92c01b9q_flc.fits"]
    # make_mask_file(in_imgname,input_list)
    # make_mask_file_old("hst_10265_01_acs_wfc_f606w_wht.fits")
    make_mask_file_actual("hst_10265_01_acs_wfc_f606w_drz.fits")