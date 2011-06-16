import pyfits
import numpy as np

import pywcs

from pytools import fileutil
from stwcs import wcsutil, updatewcs
from stwcs.wcsutil import wcscorr

import linearfit

wcs_keys = ['CRVAL1','CRVAL2','CD1_1','CD1_2','CD2_1','CD2_2','CRPIX1','CRPIX2','ORIENTAT']

def update_from_shiftfile(shiftfile,force=False):
    """ Update headers of images specified in shiftfile with shifts from shiftfile.
    """
    f = open(shiftfile)
    shift_lines = f.readlines()
    f.close()

    # interpret header of shift file
    for line in shift_lines:
        if 'refimage' in line:
            refimage = line.split(':')[-1]
            refimage = refimage[:refimage.find('[wcs]')]
            break

    # Now read in numerical values from shiftfile
    type_list = {'names':('fnames','xsh','ysh','rot','scale'),
                 'formats':('S24','f4','f4','f4','f4')}
    sdict = np.loadtxt(shiftfile,dtype=type_list,unpack=True)
    # loadtxt now returns a list of ndarray columns... need to fix this!
    for img in sdict:
        updatewcs_with_shift(img['fnames'],refimage,
                rot=img['rot'],scale=img['scale'],
                xsh=img['xsh'],ysh=img['ysh'],force=force)

def updatewcs_with_shift(image,reference,wcsname=None,rot=0.0,scale=1.0,xsh=0.0,ysh=0.0,
                            verbose=False,force=False):

    """ Update the SCI headers in 'image' based on the fit provided as determined
        in the WCS specified by 'reference'.  The fit should be a 2-D matrix as
        generated for use with 'make_vector_plot()'.

        Algorithm
        ----------
        Apply the following steps to the WCS of each of the input image's chips:
            1- compute RA/Dec with full distortion correction for
                reference point as (Rc_i,Dc_i)
            2- find the Xc,Yc for each Rc_i,Dc_i and get the difference from the
                CRPIX position for the reference WCS as (dXc_i,dYc_i)
            3- apply fit (rot&scale) to (dXc_i,dYc_i) then apply shift, then add
                CRPIX back to get new (Xcs_i,Ycs_i) position
            4- compute (Rcs_i,Dcs_i) as the sky coordinates for (Xcs_i,Ycs_i)
            5- compute delta of (Rcs_i-Rc_i, Dcs_i-Dcs_i) as (dRcs_i,dDcs_i)
            6- apply the fit to the chip's undistorted CD matrix, the apply linear
                distortion terms back in to create a new CD matrix
            7- add (dRcs_i,dDcs_i) to CRVAL of the reference chip's WCS
            8- update header with new WCS values
    """
    # if input reference is a ref_wcs file from tweakshifts, use it
    if isinstance(reference, wcsutil.HSTWCS) or isinstance(reference, pywcs.pywcs.WCS):
        wref = reference
    else:
        refimg = pyfits.open(reference)
        wref = None
        for extn in refimg:
            if extn.header.has_key('extname') and extn.header['extname'] == 'WCS':
                wref = pywcs.WCS(refimg['wcs'].header)
                break
        refimg.close()
        # else, we have presumably been provided a full undistorted image
        # as a reference, so use it with HSTWCS instead
        if wref is None:
            wref = wcsutil.HSTWCS(reference)

    # Now that we are sure we have a good reference WCS to use, continue with the update
    print '\n....Updating header for ',image,'...\n'

    # reset header WCS keywords to original (OPUS generated) values
    numextn = fileutil.countExtn(image)
    archive_wcsname = ""
    if numextn > 0:
        # Create initial WCSCORR extension
        wcscorr.init_wcscorr(image,force=force)

        extlist = []
        for extn in xrange(1,numextn+1):
            extlist.append(('SCI',extn))
    else:
        extlist = [0]
        archive_wcsname = "DRZ_"+fileutil.getDate()

    fimg = pyfits.open(image,mode='update')
    # Process MEF images...
    for ext in extlist:
        if verbose:
            print 'Processing %s[',ext,']'
        chip_wcs = wcsutil.HSTWCS(image,ext=ext)

        update_refchip_with_shift(chip_wcs,wref,rot=rot,scale=scale,xsh=xsh,ysh=ysh)
        if wcsname in [None,' ','','INDEF']:
            wcsname = 'TWEAK'
        # Update FITS file with newly updated WCS for this chip
        if numextn > 0:
            extnum = fileutil.findExtname(fimg,ext[0],ext[1])
        else:
            extnum = ext
        update_wcs(fimg,extnum,chip_wcs,wcsname=wcsname,verbose=verbose)
        
    if numextn > 0:
        # Update WCSCORR table with new WCS information
        wcscorr.update_wcscorr(fimg,wcs_id=wcsname)    
    fimg.close()

def updatewcs_with_fit(image,reference,wcsname=None,rot=0.0,scale=1.0,xsh=0.0,ysh=0.0,
                            verbose=False,force=False):
    """ Update the SCI headers in 'image' based on the fit provided as determined
        in the WCS specified by 'reference'.  The fit should be a 2-D matrix as
        generated for use with 'make_vector_plot()'.

        Algorithm
        ----------
        Apply the following steps to the input image's reference chip WCS:
            1- compute RA/Dec with full distortion correction for
                reference point as (Rc_i,Dc_i)
            2- find the Xc,Yc for each Rc_i,Dc_i and get the difference from the
                CRPIX position for the reference WCS as (dXc_i,dYc_i)
            3- apply fit (rot&scale) to (dXc_i,dYc_i) then apply shift, then add
                CRPIX back to get new (Xcs_i,Ycs_i) position
            4- compute (Rcs_i,Dcs_i) as the sky coordinates for (Xcs_i,Ycs_i)
            5- compute delta of (Rcs_i-Rc_i, Dcs_i-Dcs_i) as (dRcs_i,dDcs_i)
            6- apply the fit to the chip's undistorted CD matrix, the apply linear
                distortion terms back in to create a new CD matrix
            7- add (dRcs_i,dDcs_i) to CRVAL of the reference chip's WCS
            8- update header with new WCS values
            9- run updatewcs to update all chips to be consistent with
                input images's updated reference chip's WCS

    """
    # if input reference is a ref_wcs file from tweakshifts, use it
    if isinstance(reference, wcsutil.HSTWCS) or isinstance(reference, pywcs.pywcs.WCS):
        wref = reference
    else:
        refimg = pyfits.open(reference)
        wref = None
        for extn in refimg:
            if extn.header.has_key('extname') and extn.header['extname'] == 'WCS':
                wref = pywcs.WCS(refimg['wcs'].header)
                break
        refimg.close()
        # else, we have presumably been provided a full undistorted image
        # as a reference, so use it with HSTWCS instead
        if wref is None:
            wref = wcsutil.HSTWCS(reference)

    # Now that we are sure we have a good reference WCS to use, continue with the update
    print '\n....Updating header for ',image,'...\n'


    # reset header WCS keywords to original (OPUS generated) values
    numextn = fileutil.countExtn(image)
    archive_wcsname = ""
    if numextn > 0:
        # Create initial WCSCORR extension
        wcscorr.init_wcscorr(image,force=force)

        extlist = []
        for extn in xrange(1,numextn+1):
            extlist.append(('SCI',extn))
        wcsutil.altwcs.restoreWCS(image,extlist,wcskey='O',clobber=True)

    else:
        next_key = wcsutil.altwcs.next_wcskey(pyfits.getheader(image))
        archive_wcsname = "DRZ_"+fileutil.getDate()
        
    # archive and update PA_V3
    fimg= pyfits.open(image,mode='update')

    #fimg[0].header.update('HPA_V3',fimg[0].header['PA_V3'])
    # Archive OPUS generated value of PA_V3, if it has not already been done
    if not fimg[0].header.has_key('PA_V3O'): 
        fimg[0].header['PA_V3O'] = fimg[0].header['PA_V3']
    # Update value of PA_V3 based on rotation from fit/shiftfile
    pav3 = (fimg[0].header['PA_V3'] + rot)%360
    fimg[0].header.update('PA_V3', pav3)
    fimg.flush()

    # for each chip in image, apply algorithm
    if numextn > 0:
        nchip,extn = updatewcs.getNrefchip(fimg)
        extver = ('sci',fimg[extn].header['extver'])
        update_key = 'O'
    else:
        extn = 0
        extver = (0)

    if verbose:
        print 'Processing ',extver
    chip_wcs = wcsutil.HSTWCS(image,extver)
    
    if numextn == 0:
        chipwcs_hdr = chip_wcs.wcs2header(idc2hdr=False)
        for key in chipwcs_hdr:
            fimg[extn].header.update(key[:7]+next_key,chipwcs_hdr[key])
        fimg[extn].header.update('WCSNAME'+next_key,archive_wcsname)

        update_key = wcsutil.altwcs.next_wcskey(fimg[0].header)
    fimg.close()

    update_refchip_with_shift(chip_wcs,wref,rot=rot,scale=scale,xsh=xsh,ysh=ysh)
    if wcsname in [None,' ','','INDEF']:
        wcsname = 'TWEAK'
    
    # step 8
    # Update the 'O' WCS (OPUS generated values) in the header
    idchdr = True
    if chip_wcs.idcscale is None:
        idchdr = False

    chipwcs_hdr = chip_wcs.wcs2header(idc2hdr=idchdr)
        
    fimg = pyfits.open(image,mode='update')
    for key in chipwcs_hdr:
        fimg[extn].header.update(key[:7]+update_key,chipwcs_hdr[key])
        fimg[extn].header.update(key,chipwcs_hdr[key])
    if update_key != 'O':
        fimg[extn].header.update('WCSNAME'+update_key,wcsname)
    fimg.close()

    # step 9
    # Apply shifted reference WCS to remainder of chips (if any)
    # and archive the new primary WCS as a new keyed WCS
        
    if numextn > 0:
        updatewcs.updatewcs(image,checkfiles=False,wcsname=wcsname)

        # Restore the 'O' WCS (OPUS generated values) from the WCSCORR table
        wcscorr.restore_file_from_wcscorr(image,id='OPUS',wcskey='O')

        # Record updated values in WCSCORR extension now
        wcscorr.archive_wcs_file(image)


def apply_db_fit(data,fit,xsh=0.0,ysh=0.0):
    xy1x = data[0]
    xy1y = data[1]
    if xsh != 0.0:
        xy1x += xsh
    if ysh != 0.0:
        xy1y += ysh
    numpts = xy1x.shape[0]
    if fit is not None:
        xy1 = np.zeros((xy1x.shape[0],2),np.float64)
        xy1[:,0] = xy1x
        xy1[:,1] = xy1y
        xy1 = np.dot(xy1,fit)
        xy1x = xy1[:,0]
        xy1y = xy1[:,1]
    return xy1x,xy1y

def update_refchip_with_shift(chip_wcs, wcslin, rot=0.0,scale=1.0,xsh=0.0,ysh=0.0):
    # compute the matrix for the scale and rotation correction
    fit = np.linalg.inv(fileutil.buildRotMatrix(-1*rot)*scale)

    # step 1
    xpix = [chip_wcs.wcs.crpix[0],chip_wcs.wcs.crpix[0]+1,chip_wcs.wcs.crpix[0]]
    ypix = [chip_wcs.wcs.crpix[1],chip_wcs.wcs.crpix[1],chip_wcs.wcs.crpix[1]+1]

    # This full transformation includes all parts of model, excluding DGEO/NPOL
    #Rc_i,Dc_i = chip_wcs.all_pix2sky(xpix,ypix,1)
    dpx,dpy = chip_wcs.det2im(xpix,ypix,1)  #Apply the detector to image correction
    dpx = xpix + (dpx[0] - xpix[0])
    dpy = ypix + (dpy[0] - ypix[0])
    spx,spy = chip_wcs.sip_pix2foc(dpx,dpy,1)
    fx = dpx + (spx - dpx)
    fy = dpy + (spy - dpy)
    Rc_i,Dc_i = chip_wcs.wcs_pix2sky(fx,fy,1)

    # step 2
    Xc_i,Yc_i = wcslin.wcs_sky2pix([Rc_i],[Dc_i],1)
    Xc_i -= wcslin.wcs.crpix[0]
    Yc_i -= wcslin.wcs.crpix[1]
    # step 3
    Xcs_i,Ycs_i = apply_db_fit([Xc_i,Yc_i],fit,xsh=-1*xsh,ysh=-1*ysh)
    Xcs_i += wcslin.wcs.crpix[0]
    Ycs_i += wcslin.wcs.crpix[1]

    chip_fit = fit
    # step 4
    Rcs_i,Dcs_i = wcslin.wcs_pix2sky(Xcs_i,Ycs_i,1)
    # step 5
    # new crval should be first member
    new_crval1 = Rcs_i[0]
    new_crval2 = Dcs_i[0]
    chip_wcs.wcs.crval = np.array([new_crval1,new_crval2])
    chip_wcs.wcs.set()
    # step 6
    # compute new sky positions (with full model) based on new CRVAL
    """
    Rc_iu,Dc_iu = chip_wcs.all_pix2sky(xpix,ypix,1) # may need to remove 'dgeo' from this step
    """
    dpx,dpy = chip_wcs.det2im(xpix,ypix,1)  #Apply the detector to image correction
    dpx = xpix + (dpx[0] - xpix[0])
    dpy = ypix + (dpy[0] - ypix[0])
    spx,spy = chip_wcs.sip_pix2foc(dpx,dpy,1)
    fx = dpx + (spx - dpx)
    fy = dpy + (spy - dpy)
    Rc_iu,Dc_iu = chip_wcs.wcs_pix2sky(fx,fy,1)

    Xc_iu,Yc_iu = wcslin.wcs_sky2pix([Rc_iu],[Dc_iu],1)
    # step 7
    # Perform rscale (linear orthogonal) fit between previously updated positions
    # and newly updated positions
    XYc_iu = np.transpose([Xc_iu,Yc_iu])
    XYcs_i = np.transpose([Xcs_i,Ycs_i])
    rfit = linearfit.fit_all(XYcs_i,XYc_iu,mode='rscale',center=[new_crval1,new_crval2],verbose=False)
    rmat = fileutil.buildRotMatrix(rfit['rot'])*rfit['scale'][0]
    
    # Step 8
    # apply final fit to CD matrix
    chip_wcs.wcs.cd = np.dot(chip_wcs.wcs.cd,rmat)
    

###
### Header keyword prefix related archive functions
###
def update_wcs(image,extnum,new_wcs,wcsname="",verbose=False):
    """ Updates the WCS of the specified extension number with the new WCS
        after archiving the original WCS.

        The value of 'new_wcs' needs to be the full HSTWCS object.
    """
    # Start by insuring that the correct value of 'orientat' has been computed
    new_wcs.setOrient()

    fimg_open=False
    if not isinstance(image,pyfits.HDUList):
        fimg = pyfits.open(image,mode='update')
        fimg_open = True
    else:
        fimg = image
        
    idchdr = True
    if new_wcs.idcscale is None:
        idchdr = False
    # Open the file for updating the WCS
    try:
        print 'Updating header for ',fimg.filename(),'[',extnum,']'
        hdr = fimg[extnum].header

        if verbose:
            print 'Updating header with new values...'
        # Insure that if a copy of the WCS has not been created yet, it will be now
        wcs_hdr = new_wcs.wcs2header(idc2hdr=idchdr)

        for key in wcs_hdr:
            hdr.update(key,wcs_hdr[key])
        hdr.update('ORIENTAT',new_wcs.orientat)
        
        if wcsname not in ['',' ',None,'INDEF','N/A']:
            # Save the newly updated WCS as an alternate WCS as well
            next_key = wcsutil.altwcs.next_wcskey(hdr)
            wcsutil.altwcs.archiveWCS(fimg,[extnum],wcskey=next_key,wcsname=wcsname)
    finally:
        if fimg_open:
            # finish up by closing the file now
            fimg.close()

def archive_wcs(hdr,suffix='O',wcsver="",force=False):
    """ Stores a copy of relevant WCS keywords as keywords with a different
    suffix (for example, CRVAL1 -> CRVAL1O)
    It will also make the copy from 'CRVAL1'+wcsver instead if 'wcsver' is specified.
    """
    # Check to see whether 'O' WCS already exists
    if suffix == 'O' and hdr.has_key(wcs_keys[0]+suffix):
        # This set of keywords already exists, so do not over-write
        return

    # Look to see whether this header was updated using prefix or (FITS compliant)
    # suffix WCS keywords
    if hdr.has_key(suffix+wcs_keys[0]):
        wcsmode='prefix'
        check_key = suffix+wcs_keys[0]
    else:
        wcsmode='suffix'
        check_key = wcs_keys[0]+suffix

    if not hdr.has_key(check_key) or (hdr.has_key(check_key) and force):
        for key in wcs_keys:
            if mode == 'prefix': check_key = (suffix+key)[:8]
            else: check_key = (key+suffix)[:8]
            # Always write out WCS keywords using FITS Paper I suffix mode
            hdr.update((key+suffix)[:8],hdr[check_key])

def archive_prefix_wcs(hdr,prefix='H',wcsprefix="",force=False):
    """ Stores a copy of relevant WCS keywords as keywords with a different
    prefix (for example, CRVAL1 -> HCRVAL1)
    It will also make the copy from OCRVAL1 instead if 'wcsprefix' is specified.
    """
    if not hdr.has_key(prefix+wcs_keys[0]) or (hdr.has_key(prefix+wcs_keys[0]) and force):
        for key in wcs_keys:
            hdr.update((prefix+key)[:8],hdr[(wcsprefix+key)[:8]])

def restore_prefix_wcs(hdr,prefix='H',wcsprefix=""):
    """ Restores the WCS archived with 'prefix' to the WCS keywords
    starting with 'wcsprefix'.  By default, it will restore the "H" keywords
    as the standard WCS keywords (no prefix).
    """
    if hdr.has_key(prefix+wcs_keys[0]):
        for key in wcs_keys:
            hdr.update((wcsprefix+key)[:8],hdr[(prefix+key)[:8]])

def restore_wcs_file(filename,prefix='H',wcsprefix=""):
    fimg = pyfits.open(filename,mode='update')
    if fimg[0].header.has_key(prefix+'PA_V3'):
        fimg[0].header['PA_V3'] = fimg[0].header[prefix+'PA_V3']
    nextns = count_chips(filename)
    for extn in range(1,nextns+1):
        hdr = fimg['sci',extn].header
        restore_wcs(hdr,prefix=prefix,wcsprefix=wcsprefix)
        restore_wcs(hdr,prefix='O',wcsprefix='')
    fimg.close()

