"""

:Authors: Warren Hack, Mihai Cara

:License: :doc:`LICENSE`

"""
from __future__ import absolute_import, division, print_function
import re
import math
import warnings

from astropy.io import fits
import numpy as np

from astropy import wcs as pywcs

from stwcs import wcsutil, updatewcs
from stwcs.wcsutil import wcscorr

__version__ = '0.3.0'
__version_date__ = '26-Oct-2018'


wcs_keys = ['CRVAL1','CRVAL2','CD1_1','CD1_2','CD2_1','CD2_2',
            'CRPIX1','CRPIX2','ORIENTAT']
blank_list = [None, '', ' ',"None","INDEF"]

if hasattr(np, 'float128'):
    ndfloat128 = np.float128
elif hasattr(np, 'float96'):
    ndfloat128 = np.float96
else:
    ndfloat128 = np.float64


def updatewcs_with_shift(image,reference,wcsname=None, reusename=False,
                         fitgeom='rscale',
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

    reusename : bool
        User can specify whether or not to over-write WCS with same name.
        [Default: False]

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
        refimg = fits.open(reference, memmap=False)
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

    if isinstance(image, fits.HDUList):
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

    # Now that we are sure we have a good reference WCS to use,
    # continue with the update
    logstr = "....Updating header for {:s}...".format(filename)
    if verbose:
        print("\n{:s}\n".format(logstr))
    else:
        print(logstr)

    # reset header WCS keywords to original (OPUS generated) values
    extlist = get_ext_list(image, extname='SCI')
    if extlist:
        if image_update:
            # Create initial WCSCORR extension
            wcscorr.init_wcscorr(image,force=force)
    else:
        extlist = [0]

    # insure that input PRIMARY WCS has been archived before overwriting
    # with new solution
    if open_image:
        fimg = fits.open(image, mode='update', memmap=False)
        image_update = True
    else:
        fimg = image

    if image_update:
        wcsutil.altwcs.archiveWCS(fimg,extlist,reusekey=True)

    # Process MEF images...
    for ext in extlist:
        logstr = "Processing {:s}[{:s}]".format(fimg.filename(),
                                                ext2str(ext))
        if verbose:
            print("\n{:s}\n".format(logstr))
        else:
            print(logstr)
        chip_wcs = wcsutil.HSTWCS(fimg,ext=ext)

        update_refchip_with_shift(chip_wcs, wref, fitgeom=fitgeom,
                    rot=rot, scale=scale, xsh=xsh, ysh=ysh,
                    fit=fit, xrms=xrms, yrms=yrms)

        # Update FITS file with newly updated WCS for this chip
        extnum = fimg.index(fimg[ext])
        update_wcs(fimg, extnum, chip_wcs, wcsname=wcsname,
                   reusename=reusename, verbose=verbose)

    if open_image:
        fimg.close()


def linearize(wcsim, wcsima, wcsref, imcrpix, f, shift, hx=1.0, hy=1.0):
    """ linearization using 5-point formula for first order derivative

    """
    x0 = imcrpix[0]
    y0 = imcrpix[1]
    p = np.asarray([[x0, y0],
                    [x0 - hx, y0],
                    [x0 - hx * 0.5, y0],
                    [x0 + hx * 0.5, y0],
                    [x0 + hx, y0],
                    [x0, y0 - hy],
                    [x0, y0 - hy * 0.5],
                    [x0, y0 + hy * 0.5],
                    [x0, y0 + hy]],
                   dtype=np.float64)
    # convert image coordinates to reference image coordinates:
    p = wcsref.wcs_world2pix(wcsim.wcs_pix2world(p, 1), 1).astype(ndfloat128)
    # apply linear fit transformation:
    p = np.dot(f, (p - shift).T).T
    # convert back to image coordinate system:
    p = wcsima.wcs_world2pix(
        wcsref.wcs_pix2world(p.astype(np.float64), 1), 1).astype(ndfloat128)

    # derivative with regard to x:
    u1 = ((p[1] - p[4]) + 8 * (p[3] - p[2])) / (6*hx)
    # derivative with regard to y:
    u2 = ((p[5] - p[8]) + 8 * (p[7] - p[6])) / (6*hy)

    return (np.asarray([u1, u2]).T, p[0])


def _inv2x2(x):
    assert(x.shape == (2,2))
    inv = x.astype(ndfloat128)
    det = inv[0,0]*inv[1,1] - inv[0,1]*inv[1,0]
    if np.abs(det) < np.finfo(np.float64).tiny:
        raise ArithmeticError('Singular matrix.')
    a = inv[0, 0]
    d = inv[1, 1]
    inv[1, 0] *= -1.0
    inv[0, 1] *= -1.0
    inv[0, 0] = d
    inv[1, 1] = a
    inv /= det
    inv = inv.astype(np.float64)
    if not np.all(np.isfinite(inv)):
        raise ArithmeticError('Singular matrix.')
    return inv


def update_refchip_with_shift(chip_wcs, wcslin, fitgeom='rscale',
                              rot=0.0, scale=1.0, xsh=0.0, ysh=0.0,
                              fit=None, xrms=None, yrms=None):
    """ Compute the matrix for the scale and rotation correction

    Parameters
    ----------
    chip_wcs: wcs object
        HST of the input image
    wcslin: wcs object
        Reference WCS from which the offsets/rotations are determined
    fitgeom: str
        NOT USED
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
        """
    # compute the matrix for the scale and rotation correction
    if fit is None:
        fit = buildFitMatrix(rot, scale)

    shift = np.asarray([xsh, ysh]) - np.dot(wcslin.wcs.crpix, fit) + wcslin.wcs.crpix

    fit = _inv2x2(fit).T if fit.shape == (2,2) else np.linalg.inv(fit).T

    cwcs = chip_wcs.deepcopy()
    cd_eye = np.eye(chip_wcs.wcs.cd.shape[0], dtype=ndfloat128)
    zero_shift = np.zeros(2, dtype=ndfloat128)

    # estimate precision necessary for iterative processes:
    maxiter = 100
    crpix2corners = np.dstack([i.flatten() for i in np.meshgrid(
        [1,chip_wcs._naxis1],
        [1,chip_wcs._naxis2])])[0] - chip_wcs.wcs.crpix
    maxUerr = 1.0e-5 / np.amax(np.linalg.norm(crpix2corners, axis=1))

    # estimate step for numerical differentiation. We need a step
    # large enough to avoid rounding errors and small enough to get a
    # better precision for numerical differentiation.
    # TODO: The logic below should be revised at a later time so that it
    # better takes into account the two competing requirements.
    hx = max(1.0, min(20.0, (chip_wcs.wcs.crpix[0] - 1.0)/100.0,
                      (chip_wcs._naxis1 - chip_wcs.wcs.crpix[0])/100.0))
    hy = max(1.0, min(20.0, (chip_wcs.wcs.crpix[1] - 1.0)/100.0,
                      (chip_wcs._naxis2 - chip_wcs.wcs.crpix[1])/100.0))

    # compute new CRVAL for the image WCS:
    crpixinref = wcslin.wcs_world2pix(
        chip_wcs.wcs_pix2world([chip_wcs.wcs.crpix],1),1)
    crpixinref = np.dot(fit, (crpixinref - shift).T).T
    chip_wcs.wcs.crval = wcslin.wcs_pix2world(crpixinref, 1)[0]
    chip_wcs.wcs.set()

    # initial approximation for CD matrix of the image WCS:
    (U, u) = linearize(cwcs, chip_wcs, wcslin, chip_wcs.wcs.crpix,
                       fit, shift, hx=hx, hy=hy)
    err0 = np.amax(np.abs(U-cd_eye)).astype(np.float64)
    chip_wcs.wcs.cd = np.dot(chip_wcs.wcs.cd.astype(ndfloat128), U).astype(np.float64)
    chip_wcs.wcs.set()

    # NOTE: initial solution is the exact mathematical solution (modulo numeric
    # differentiation). However, e.g., due to rounding errors, approximate
    # numerical differentiation, the solution may be improved by performing
    # several iterations. The next step will try to perform
    # fixed-point iterations to "improve" the solution
    # but this is not really required.

    # Perform fixed-point iterations to improve the approximation
    # for CD matrix of the image WCS (actually for the U matrix).
    for i in range(maxiter):
        (U, u) = linearize(chip_wcs, chip_wcs, wcslin, chip_wcs.wcs.crpix,
                           cd_eye, zero_shift, hx=hx, hy=hy)
        err = np.amax(np.abs(U-cd_eye)).astype(np.float64)
        if err > err0:
            break
        chip_wcs.wcs.cd = np.dot(chip_wcs.wcs.cd, U).astype(np.float64)
        chip_wcs.wcs.set()
        if err < maxUerr:
            break
        err0 = err

    if xrms is not None:
        chip_wcs.wcs.crder = np.array([xrms,yrms])


###
### Header keyword prefix related archive functions
###
def update_wcs(image,extnum,new_wcs,wcsname="",reusename=False,verbose=False):
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

    reusename : bool
        User can choose whether to over-write WCS with same name or not.
        [Default: False]

    verbose : bool, int
        Print extra messages during processing? [Default: False]

    """
    # Start by insuring that the correct value of 'orientat' has been computed
    new_wcs.setOrient()

    fimg_open=False
    if not isinstance(image, fits.HDUList):
        fimg = fits.open(image, mode='update', memmap=False)
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
    if is_blank(wcsname):
        wcsname = 'TWEAK'
    if not reusename:
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
            print(logstr)

        hdr = fimg[extnum].header

        if verbose:
            print('    with WCS of')
            new_wcs.printwcs()
            print("WCSNAME  : ",wcsname)

        # Insure that if a copy of the WCS has not been created yet, it will be now
        wcs_hdr = new_wcs.wcs2header(idc2hdr=idchdr, relax=True)

        for key in wcs_hdr:
            hdr[key] = wcs_hdr[key]
        hdr['ORIENTAT'] = new_wcs.orientat
        hdr['WCSNAME'] = wcsname
        updateNEXTENDKw(fimg)

        # Only if this image was opened in update mode should this
        # newly updated WCS be archived, as it will never be written out
        # to a file otherwise.
        if fimg_update:
            if not reusename:
                # Save the newly updated WCS as an alternate WCS as well
                wkey = wcsutil.altwcs.next_wcskey(fimg,ext=extnum)
            else:
                wkey = wcsutil.altwcs.getKeyFromName(hdr,wcsname)

            # wcskey needs to be specified so that archiveWCS will create a
            # duplicate WCS with the same WCSNAME as the Primary WCS
            wcsutil.altwcs.archiveWCS(fimg,[extnum],wcsname=wcsname,
                wcskey=wkey, reusekey=reusename)
    finally:
        if fimg_open:
            # finish up by closing the file now
            fimg.close()


def create_unique_wcsname(fimg, extnum, wcsname):
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
    wnames = list(wcsutil.altwcs.wcsnames(fimg, ext=extnum).values())
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
    
    
##########################################
#
#  Functions from drizzlepac.util
#
##########################################

def is_blank(val):
    """ Determines whether or not a value is considered 'blank'.
    """
    return val in blank_list

def updateNEXTENDKw(fobj):
    """ Update NEXTEND keyword in PRIMARY header (if present) to accurately
        reflect the number of extensions in the MEF file.
    """
    if 'nextend' in fobj[0].header:
        fobj[0].header['nextend'] = len(fobj)-1

def buildFitMatrix(rot, scale=1):
    if hasattr(rot, '__iter__'):
        rx = rot[0]
        ry = rot[1]
    else:
        rx = float(rot)
        ry = rx
    if hasattr(scale, '__iter__'):
        sx = scale[0]
        sy = scale[1]
    else:
        sx = float(scale)
        sy = sx
    m = np.array(
        [
            [ sx*np.cos(np.deg2rad(rx)), -sx*np.sin(np.deg2rad(rx)) ],
            [ sy*np.sin(np.deg2rad(ry)),  sy*np.cos(np.deg2rad(ry)) ]
        ]
    )
    return m


##########################################
#
#  Functions from skypac.utils
#
##########################################

def ext2str(ext, compact=False, default_extver=1):
    """
    Return a string representation of an extension specification.

    Parameters
    ----------
    ext : tuple, int, str
        Extension specification can be a tuple of the form (str,int), e.g.,
        ('sci',1), an integer (extension number), or a string (extension
        name).

    compact : bool (Default = False)
        If `compact`\ =\ `True` the returned string will have extension
        name quoted and separated by a comma from the extension number,
        e.g., "'sci',1".
        If `compact`\ =\ `False` the returned string will have extension
        version immediately follow the extension name, e.g., 'sci1'.

    default_extver : int (Default = 1)
        Specifies the extension version to be used when the `ext` parameter
        is a string (extension name).

    Returns
    -------
    strext : str
        String representation of extension specification `ext`.

    Raises
    ------
    TypeError
        Unexpected extension type.

    Examples
    --------
    >>> ext2str('sci',compact=False,default_extver=6)
    "'sci',6"
    >>> ext2str(('sci',2))
    "'sci',2"
    >>> ext2str(4)
    '4'
    >>> ext2str('dq')
    "'dq',1"
    >>> ext2str('dq',default_extver=2)
    "'dq',2"
    >>> ext2str('sci',compact=True,default_extver=2)
    'sci2'

    """
    if isinstance(ext, tuple) and len(ext) == 2 and \
        isinstance(ext[0], str) and isinstance(ext[1], int):
        if compact:
            return "{:s}{:d}".format(ext[0], ext[1])
        else:
            return "\'{:s}\',{:d}".format(ext[0], ext[1])

    elif isinstance(ext, int):
        return "{:d}".format(ext)

    elif isinstance(ext,str):
        if compact:
            extver = '' if default_extver is None else '{:d}'.format(default_extver)
            return "{:s}{:s}".format(ext, extver)
        else:
            extver = '' if default_extver is None else ',{:d}'.format(default_extver)
            return "\'{:s}\'{:s}".format(ext, extver)

    else:
        raise TypeError("Unexpected extension type.")

def get_ext_list(img, extname='SCI'):
    """
    Return a list of all extension versions of `extname` extensions.
    `img` can be either a file name or a `astropy.io.fits.HDUList` object.

    This function is similar to :py:func:`get_extver_list`, the main
    difference being that it returns a list of fully qualified extensions: \
    either tuples of the form (\ `extname`\ ,\ `extver`\ ) or integer extension
    numbers (when `extname`\ =\ `None`).

    See Also
    --------
    get_extver_list, count_extensions

    Examples
    --------
    >>> get_ext_list('j9irw1rqq_flt.fits',extname='SCI')
    [('SCI', 1), ('SCI', 2)]
    >>> get_ext_list('j9irw1rqq_flt.fits',extname=None)
    [1, 2, 3, 4, 5, 6, 8, 9, 10, 11]

    """
    if isinstance(img, str):
        hdulist = fits.open(img)
    else:
        hdulist = img # assume it is already a FITS HDUList object
        
    # when extver is None - return the range of all 'image'-like FITS extensions
    if extname is None:
        extn = []
        for i in range(len(hdulist)):
            hdr = hdulist[i].header
            if not ('NAXIS' in hdr and hdr['NAXIS'] == 2):
                continue
            if 'XTENSION' in hdr and \
               hdr['XTENSION'].upper().strip() == 'IMAGE':
                extn.append(i)
            elif 'SIMPLE' in hdr:
                extn.append(i)
        if doRelease:
            img.release()
        return extn
    
    extname = extname.upper()

    extver = []
    for e in hdulist:
        #if not isinstance(e, fits.ImageHDU): continue
        #hkeys = map(str.upper, e.header.keys())
        if 'EXTNAME' in e.header and e.header['EXTNAME'].upper() == extname:
            extver.append(e.header['EXTVER'] if 'EXTVER' in e.header else 1)

    if extname is None:
        return extver

    extlist = [ (extname,extv) for extv in extver ]

    return extlist


##########################################
#
#  Functions from fitsblender.blendheaders
#
##########################################

def remove_distortion_keywords(hdr):
    """
    Remove WCS distortion related keywords from the input header
    """
    distortion_kws = ['TDDALPHA','TDDBETA',
                        'D2IMEXT','D2IMERR',
                        'DGEOEXT','NPOLEXT']

    # Remove any reference to TDD correction from
    #    distortion-corrected products
    # We also need to remove the D2IM* keywords so that HSTWCS/PyWCS
    # does not try to look for non-existent extensions
    for kw in distortion_kws:
        if kw in hdr:
            try:
                del hdr[kw]
            except KeyError:
                pass

    # Remove '-SIP' from CTYPE for output product
    if 'ctype1' in hdr and hdr['ctype1'].find('SIP') > -1:
            hdr['ctype1'] = hdr['ctype1'][:-4]
            hdr['ctype2'] = hdr['ctype2'][:-4]

    # Remove SIP coefficients from DRZ product
    for k in list(hdr.items()):
        if (k[0][:2] in ['A_','B_']) or (k[0][:3] in ['IDC','SCD'] and k[0] != 'IDCTAB') or \
        (k[0][:6] in ['SCTYPE','SCRVAL','SNAXIS','SCRPIX']):
            try:
                del hdr[k[0]]
            except KeyError:
                pass

    # Remove paper IV related keywords related to the
    #   DGEO correction here
    for k in list(hdr.items()):
        if (k[0][:2] == 'DP'):
            try:
                del hdr[k[0]+'*']
                del hdr[k[0]+'.*']
                del hdr[k[0]+'.*.*']
            except:
                print("ERROR (bleandheaders.remove_distortion_keywords) trying to delete \'{:s}\' in the header.".format(k[0]+'*'))
                pass
        if (k[0][:2] == 'CP'):
            try:
                del hdr[k[0]]
            except KeyError:
                pass
