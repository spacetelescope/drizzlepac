import pyfits
import numpy as np

import pywcs

from pytools import fileutil
from stwcs import wcsutil, updatewcs
from stwcs.wcsutil import wcscorr

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
    if wcsname is None or wcsname == 'TWEAK':
        wcsname = 'TWEAK_'+fileutil.getDate()
    elif wcsname == '':
        wcsname = ' '
    
    # step 8
    # Update the 'O' WCS (OPUS generated values) in the header
    idchdr = True
    if chip_wcs.idcscale is None:
        idchdr = False

    chipwcs_hdr = chip_wcs.wcs2header(idc2hdr=idchdr)
    fimg = pyfits.open(image,mode='update')
    for key in chipwcs_hdr:
        fimg[extn].header.update(key[:7]+update_key,chipwcs_hdr[key])
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
    numpts = xy1x.shape[0]
    if fit is not None:
        xy1 = np.zeros((xy1x.shape[0],2),np.float64)
        xy1[:,0] = xy1x
        xy1[:,1] = xy1y
        xy1 = np.dot(xy1,fit)
        xy1x = xy1[:,0] + xsh
        xy1y = xy1[:,1] + ysh
    return xy1x,xy1y

def update_refchip_with_shift(chip_wcs, wcslin, rot=0.0,scale=1.0,xsh=0.0,ysh=0.0):
    # compute the matrix for the scale and rotation correction
    fit = scale*fileutil.buildRotMatrix(rot)

    # step 1
    xpix = [chip_wcs.wcs.crpix[0],chip_wcs.wcs.crpix[0]+1,chip_wcs.wcs.crpix[0]]
    ypix = [chip_wcs.wcs.crpix[1],chip_wcs.wcs.crpix[1],chip_wcs.wcs.crpix[1]+1]

    # This full transformation includes all parts of model, including DGEO
    Rc_i,Dc_i = chip_wcs.all_pix2sky(xpix,ypix,1)

    # step 2
    Xc_i,Yc_i = wcslin.wcs_sky2pix([Rc_i],[Dc_i],1)
    Xc_i -= wcslin.wcs.crpix[0]#+xsh
    Yc_i -= wcslin.wcs.crpix[1]#+ysh
    # step 3
    Xcs_i,Ycs_i = apply_db_fit([Xc_i,Yc_i],fit,xsh=-1*xsh,ysh=-1*ysh)
    Xcs_i += wcslin.wcs.crpix[0]#+xsh
    Ycs_i += wcslin.wcs.crpix[1]#+ysh

    chip_fit = fit
    # step 4
    Rcs_i,Dcs_i = wcslin.wcs_pix2sky(Xcs_i,Ycs_i,1)
    # step 5
    # new crval should be first member
    chip_wcs.wcs.crval = np.array([Rcs_i[0],Dcs_i[0]])
    new_crval1 = Rcs_i[0]
    new_crval2 = Dcs_i[0]
    # step 6
    # see about computing the CD matrix directly from the 3 points around
    # the shifted/rotated CRPIX position in the output frame as projected
    # back onto the sky
    am1 = Rcs_i[1]-new_crval1
    bm1 = Rcs_i[2]-new_crval1
    cm1 = Dcs_i[1]-new_crval2
    dm1 = Dcs_i[2]-new_crval2
    chip_wcs.wcs.cd = np.array([[am1*np.cos(new_crval2*np.pi/180),
                                bm1*np.cos(new_crval2*np.pi/180)],[cm1,dm1]],
                                dtype=np.float64)

###
### Header keyword prefix related archive functions
###
def update_wcs(image,extver,new_wcs,extname='SCI',wcsprefix="",verbose=False):
    """ Updates the WCS of the specified (extname,extver) with the new WCS
        after archiving the original WCS.

        The value of 'new_wcs' can either be the full HSTWCS object or just
        the pywcs object within the HSTWCS object.
    """
    # If an HSTWCS object is passed in for 'new_wcs', we only need the
    # pywcs object within
    if isinstance(new_wcs,wcsutil.HSTWCS):
        new_wcs = new_wcs.wcs

    # Open the file for updating the WCS
    try:
        print 'Updating header for ',image,'[',extname,',',extver,']'
        fimg = pyfits.open(image,mode='update')
        hdr = fimg[extname,extver].header

        if verbose:
            print 'Updating header with new values...'
        # Insure that if a copy of the WCS has not been created yet, it will be now
        #archive_wcs(hdr,wcsprefix=wcsprefix)

        # update the values in the WCS
        hdr[(wcsprefix+'CRVAL1')[:8]] = new_wcs.crval[0]
        hdr[(wcsprefix+'CRVAL2')[:8]] = new_wcs.crval[1]
        hdr[(wcsprefix+'CD1_1')[:8]] = new_wcs.cd[0][0]
        hdr[(wcsprefix+'CD1_2')[:8]] = new_wcs.cd[0][1]
        hdr[(wcsprefix+'CD2_1')[:8]] = new_wcs.cd[1][0]
        hdr[(wcsprefix+'CD2_2')[:8]] = new_wcs.cd[1][1]
        # Recompute and update ORIENTAT keyword
        orientat = fileutil.RADTODEG(np.arctan2(new_wcs.cd[0][1],new_wcs.cd[1][1]))

        # Reset DGEOEXT in reference chip header to get updatewcs to reset
        # the DGEO extensions based on the updated WCS keywords
        # but only if we are updating the archived version of the keywords
        if wcsprefix is not '' and fimg[0].header['NPOLFILE'] not in ['',' ','N/A']:
            hdr['NPOLEXT'] = ''

    finally:
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

