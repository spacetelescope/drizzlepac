import re

import pyfits
import numpy as np

import pywcs

from stsci.tools import fileutil, logutil
from stwcs import wcsutil, updatewcs
from stwcs.wcsutil import wcscorr

from . import util

import linearfit

log = logutil.create_logger(__name__)

wcs_keys = ['CRVAL1','CRVAL2','CD1_1','CD1_2','CD2_1','CD2_2','CRPIX1','CRPIX2','ORIENTAT']

def update_from_shiftfile(shiftfile,wcsname=None,force=False):
    """
    Update headers of all images specified in shiftfile with shifts
    from shiftfile.

    Parameters
    ----------
    shiftfile : str
        Filename of shiftfile.

    wcsname : str
        Label to give to new WCS solution being created by this fit. If
        a value of None is given, it will automatically use 'TWEAK' as the
        label. [Default =None]

    force : bool
        Update header even though WCS already exists with this solution or
        wcsname? [Default=False]

    """
    f = open(fileutil.osfn(shiftfile))
    shift_lines = f.readlines()
    f.close()

    # interpret header of shift file
    for line in shift_lines:
        if 'refimage' in line:
            refimage = line.split(':')[-1]
            refimage = refimage[:refimage.find('[wcs]')]
            break

    # Now read in numerical values from shiftfile
    type_list = {'names':('fnames','xsh','ysh','rot','scale','xrms','yrms'),
                 'formats':('S24','f4','f4','f4','f4','f4','f4')}
    try:
        sdict = np.loadtxt(shiftfile,dtype=type_list,unpack=False)
    except IndexError:
        tlist = {'names':('fnames','xsh','ysh','rot','scale'),
                     'formats':('S24','f4','f4','f4','f4')}
        s = np.loadtxt(shiftfile,dtype=tlist,unpack=False)
        sdict = np.zeros([s['fnames'].shape[0],],dtype=type_list)
        for sname in s.dtype.names:
            sdict[sname] = s[sname]

    for img in sdict:
        updatewcs_with_shift(img['fnames'], refimage, wcsname=wcsname,
                rot=img['rot'], scale=img['scale'],
                xsh=img['xsh'], ysh=img['ysh'],
                xrms=img['xrms'], yrms=img['yrms'],
                force=force)

def updatewcs_with_shift(image,reference,wcsname=None,
                        rot=0.0,scale=1.0,xsh=0.0,ysh=0.0,fit=None,
                        xrms=None, yrms = None,
                        verbose=False,force=False,sciext='SCI'):

    """
    Update the SCI headers in 'image' based on the fit provided as determined
    in the WCS specified by 'reference'.  The fit should be a 2-D matrix as
    generated for use with 'make_vector_plot()'.

    Notes
    -----
    The algorithm used to apply the provided fit solution to the image
    involves applying the following steps to the WCS of each of the
    input image's chips:

    1. compute RA/Dec with full distortion correction for
            reference point as (Rc_i,Dc_i)

    2. find the Xc,Yc for each Rc_i,Dc_i and get the difference from the
            CRPIX position for the reference WCS as (dXc_i,dYc_i)

    3. apply fit (rot&scale) to (dXc_i,dYc_i) then apply shift, then add
            CRPIX back to get new (Xcs_i,Ycs_i) position

    4. compute (Rcs_i,Dcs_i) as the sky coordinates for (Xcs_i,Ycs_i)

    5. compute delta of (Rcs_i-Rc_i, Dcs_i-Dcs_i) as (dRcs_i,dDcs_i)

    6. apply the fit to the chip's undistorted CD matrix, the apply linear
            distortion terms back in to create a new CD matrix

    7. add (dRcs_i,dDcs_i) to CRVAL of the reference chip's WCS

    8. update header with new WCS values

    Parameters
    ----------
    image : str or PyFITS.HDUList object
        Filename, or PyFITS object, of image with WCS to be updated.
        All extensions with EXTNAME matches the value of the 'sciext'
        parameter value (by default, all 'SCI' extensions) will be updated.

    reference : str
        Filename of image/headerlet (FITS file) which contains the WCS
        used to define the tangent plane in which all the fit parameters
        (shift, rot, scale) were measured.

    wcsname : str
        Label to give to new WCS solution being created by this fit. If
        a value of None is given, it will automatically use 'TWEAK' as the
        label. If a WCS has a name with this specific value, the code will
        automatically append a version ID using the format '_n', such as
        'TWEAK_1', 'TWEAK_2',or 'TWEAK_update_1'.
        [Default =None]

    rot : float
        Amount of rotation measured in fit to be applied.
        [Default=0.0]

    scale : float
        Amount of scale change measured in fit to be applied.
        [Default=1.0]

    xsh : float
        Offset in X pixels from defined tangent plane to be applied to image.
        [Default=0.0]

    ysh : float
        Offset in Y pixels from defined tangent plane to be applied to image.
        [Default=0.0]

    fit : arr
        Linear coefficients for fit
        [Default = None]

    xrms : float
        RMS of fit in RA (in decimal degrees) that will be recorded as
        CRDER1 in WCS and header
        [Default = None]

    yrms : float
        RMS of fit in Dec (in decimal degrees) that will be recorded as
        CRDER2 in WCS and header
        [Default = None]

    verbose : bool
        Print extra messages during processing? [Default=False]

    force : bool
        Update header even though WCS already exists with this solution or
        wcsname? [Default=False]

    sciext : string
        Value of FITS EXTNAME keyword for extensions with WCS headers to
        be updated with the fit values. [Default='SCI']

    """
    # if input reference is a ref_wcs file from tweakshifts, use it
    if isinstance(reference, wcsutil.HSTWCS) or isinstance(reference, pywcs.WCS):
        wref = reference
    else:
        refimg = pyfits.open(reference)
        wref = None
        for extn in refimg:
            if 'extname' in extn.header and extn.header['extname'] == 'WCS':
                wref = pywcs.WCS(refimg['wcs'].header)
                break
        refimg.close()
        # else, we have presumably been provided a full undistorted image
        # as a reference, so use it with HSTWCS instead
        if wref is None:
            wref = wcsutil.HSTWCS(reference)

    if isinstance(image,pyfits.HDUList):
        open_image = False
        filename = image.filename()
        if image.fileinfo(0)['filemode'] is 'update':
            image_update = True
        else:
            image_update = False
    else:
        open_image = True
        filename = image
        image_update = None

    # Now that we are sure we have a good reference WCS to use, continue with the update
    logstr = '....Updating header for %s...'%filename
    if verbose:
        print '\n'+logstr+'\n'
    else:
        log.info(logstr)

    # reset header WCS keywords to original (OPUS generated) values
    numextn = fileutil.countExtn(image,extname=sciext)

    if numextn > 0:
        if image_update is True:
            # Create initial WCSCORR extension
            wcscorr.init_wcscorr(image,force=force)

        extlist = []
        for extn in xrange(1,numextn+1):
            extlist.append((sciext,extn))
    else:
        extlist = [0]

    # insure that input PRIMARY WCS has been archived before overwriting
    # with new solution
    if open_image:
        fimg = pyfits.open(image,mode='update')
        image_update = True
    else:
        fimg = image

    if image_update is True:
        wcsutil.altwcs.archiveWCS(fimg,extlist,reusekey=True)

    # Process MEF images...
    for ext in extlist:
        logstr = 'Processing %s[%s]'%(fimg.filename(),str(ext))
        if verbose:
            print '\n'+logstr+'\n'
        else:
            log.info(logstr)
        chip_wcs = wcsutil.HSTWCS(fimg,ext=ext)

        update_refchip_with_shift(chip_wcs, wref,
                    rot=rot, scale=scale, xsh=xsh, ysh=ysh,
                    fit=fit, xrms=xrms, yrms=yrms)
        if wcsname in [None,' ','','INDEF']:
            wcsname = 'TWEAK'
        # Update FITS file with newly updated WCS for this chip
        if numextn > 0:
            extnum = fileutil.findExtname(fimg,ext[0],ext[1])
        else:
            extnum = ext

        update_wcs(fimg,extnum,chip_wcs,wcsname=wcsname,verbose=verbose)

#    if numextn > 0:
#        # Update WCSCORR table with new WCS information
#        wcscorr.update_wcscorr(fimg,wcs_id=wcsname)
    if open_image:
        fimg.close()

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

def update_refchip_with_shift(chip_wcs, wcslin,
                            rot=0.0, scale=1.0, xsh=0.0, ysh=0.0,
                            fit=None, xrms=None, yrms=None):
    # compute the matrix for the scale and rotation correction
    if fit is None:
        fitmat = fileutil.buildRotMatrix(-1*rot)*scale
    else:
        fitmat = fit

    rmode='rscale'
    fit = np.linalg.inv(fitmat)

    # step 1
    xpix = [chip_wcs.wcs.crpix[0],chip_wcs.wcs.crpix[0]+1,chip_wcs.wcs.crpix[0]]
    ypix = [chip_wcs.wcs.crpix[1],chip_wcs.wcs.crpix[1],chip_wcs.wcs.crpix[1]+1]

    # This full transformation includes all parts of model, excluding DGEO/NPOL
    Rc_i,Dc_i = chip_wcs.wcs_pix2sky(xpix,ypix,1)

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
    Rc_iu,Dc_iu = chip_wcs.wcs_pix2sky(xpix,ypix,1)
    Xc_iu,Yc_iu = wcslin.wcs_sky2pix([Rc_iu],[Dc_iu],1)
    # step 7
    # Perform rscale (linear orthogonal) fit between previously updated positions
    # and newly updated positions
    XYc_iu = np.transpose([Xc_iu,Yc_iu])
    XYcs_i = np.transpose([Xcs_i,Ycs_i])
    #rfit = linearfit.fit_all(XYcs_i,XYc_iu,mode=rmode,center=[new_crval1,new_crval2],verbose=False)
    #rmat = rfit['fit_matrix']
    rfit = linearfit.fit_all(XYcs_i,XYc_iu,mode='rscale',center=[new_crval1,new_crval2],verbose=False)
    rmat = fileutil.buildRotMatrix(rfit['rot'])*rfit['scale'][0]

    # Step 8
    # apply final fit to CD matrix
    chip_wcs.wcs.cd = np.dot(chip_wcs.wcs.cd,rmat)

    # Step 9
    # record any reported error for this fit
    if xrms is not None:
        chip_wcs.wcs.crder = np.array([xrms,yrms])
###
### Header keyword prefix related archive functions
###
def update_wcs(image,extnum,new_wcs,wcsname="",verbose=False):
    """
    Updates the WCS of the specified extension number with the new WCS
    after archiving the original WCS.

    The value of 'new_wcs' needs to be the full
    HSTWCS object.

    Parameters
    ----------
    image : str
        Filename of image with WCS that needs to be updated

    extnum : int
        Extension number for extension with WCS to be updated/replaced

    new_wcs : object
        Full HSTWCS object which will replace/update the existing WCS

    wcsname : str
        Label to give newly updated WCS

    verbose : bool, int
        Print extra messages during processing? [Default: False]

    """
    # Start by insuring that the correct value of 'orientat' has been computed
    new_wcs.setOrient()

    fimg_open=False
    if not isinstance(image,pyfits.HDUList):
        fimg = pyfits.open(image,mode='update')
        fimg_open = True
        fimg_update = True
    else:
        fimg = image
        if fimg.fileinfo(0)['filemode'] is 'update':
            fimg_update = True
        else:
            fimg_update = False

    # Determine final (unique) WCSNAME value, either based on the default or
    # user-provided name
    if wcsname in ['',' ',None,'INDEF','N/A']:
        wcsname = 'TWEAK'
    wcsname = create_unique_wcsname(fimg, extnum, wcsname)

    idchdr = True
    if new_wcs.idcscale is None:
        idchdr = False
    # Open the file for updating the WCS
    try:
        logstr = 'Updating header for %s[%s]'%(fimg.filename(),str(extnum))
        if verbose:
            print(logstr)
        else:
            log.info(logstr)

        hdr = fimg[extnum].header

        if verbose:
            log.info('    with WCS of')
            new_wcs.printwcs()
            print "WCSNAME  : ",wcsname

        # Insure that if a copy of the WCS has not been created yet, it will be now
        wcs_hdr = new_wcs.wcs2header(idc2hdr=idchdr)

        for key in wcs_hdr:
            hdr.update(key,wcs_hdr[key])
        hdr.update('ORIENTAT',new_wcs.orientat)
        hdr.update('WCSNAME',wcsname)
        util.updateNEXTENDKw(fimg)

        # Only if this image was opened in update mode should this
        # newly updated WCS be archived, as it will never be written out
        # to a file otherwise.
        if fimg_update:
            # Save the newly updated WCS as an alternate WCS as well
            wkey = wcsutil.altwcs.next_wcskey(fimg,ext=extnum)
            # wcskey needs to be specified so that archiveWCS will create a
            # duplicate WCS with the same WCSNAME as the Primary WCS
            wcsutil.altwcs.archiveWCS(fimg,[extnum],wcsname=wcsname,wcskey=wkey)

    finally:
        if fimg_open:
            # finish up by closing the file now
            fimg.close()

def create_unique_wcsname(fimg, extnum ,wcsname):
    """
    This function evaluates whether the specified wcsname value has
    already been used in this image.  If so, it automatically modifies
    the name with a simple version ID using wcsname_NNN format.

    Parameters
    ----------
    fimg : obj
        PyFITS object of image with WCS information to be updated

    extnum : int
        Index of extension with WCS information to be updated

    wcsname : str
        Value of WCSNAME specified by user for labelling the new WCS

    Returns
    -------
    uniqname : str
        Unique WCSNAME value

    """
    wnames = wcsutil.altwcs.wcsnames(fimg, ext=extnum).values()
    if wcsname not in wnames:
        uniqname = wcsname
    else:
        # setup pattern to match
        rpatt = re.compile(wcsname+'_\d')
        index = 0
        for wname in wnames:
            rmatch = rpatt.match(wname)
            if rmatch:
                # get index
                n = int(wname[wname.rfind('_')+1:])
                if n > index: index = 1
        index += 1 # for use with new name
        uniqname = "%s_%d"%(wcsname,index)
    return uniqname
