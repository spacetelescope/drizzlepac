"""

:Authors: Warren Hack, Mihai Cara

:License: :doc:`LICENSE`

"""
import re
import sys

from astropy.io import fits
import numpy as np

from astropy import wcs as pywcs

from stsci.tools import fileutil, logutil
from stwcs import wcsutil
from stwcs.wcsutil import wcscorr, altwcs
from stsci.skypac.utils import get_ext_list, ext2str

from . import util
from . import linearfit

__version__ = '0.4.0'
__version_date__ = '13-July-2020'

log = logutil.create_logger(__name__, level=logutil.logging.NOTSET)

wcs_keys = ['CRVAL1', 'CRVAL2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2',
            'CRPIX1', 'CRPIX2', 'ORIENTAT']

_SHIFT_COLNAMES = ['xsh', 'ysh', 'rot', 'scale', 'xrms', 'yrms']


def update_from_shiftfile(shiftfile, wcsname=None, force=False):
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
    with open(fileutil.osfn(shiftfile)) as f:
        lines = f.readlines()

    refimage = None
    shift_info = {}

    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue

        if refimage is not None and ('refimage' in line or
                                     'reference' in line):
            refimage = (line.split(':')[-1]).strip()
            idx = refimage.find('[wcs]')
            if idx >= 0:
                refimage = refimage[:idx].lstrip()
            continue

        cols = list(map(str.strip, line.split()))

        if len(cols) not in [5, 7]:
            raise ValueError("Unsupported shift file format: invalid number "
                             "of columns.")

        shift_info[cols[0]] = {
            k: float(v) for k, v in zip(_SHIFT_COLNAMES, cols[1:])
        }

    for filename, pars in shift_info:
        updatewcs_with_shift(filename, refimage, wcsname=wcsname,
                             force=force, **pars)


def updatewcs_with_shift(image, reference, wcsname='TWEAK', reusename=False,
                         fitgeom='rscale', rot=0.0, scale=1.0,
                         xsh=0.0, ysh=0.0, fit=None, xrms=None, yrms=None,
                         verbose=False, force=False, sciext='SCI'):
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

    wcsname : str, None, optional
        Label to give to new WCS solution being created by this fit. If
        a value of None is given, it will automatically use 'TWEAK' as the
        label.
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
        for extn in refimg:
            if 'extname' in extn.header and extn.header['extname'] == 'WCS':
                wref = pywcs.WCS(refimg['wcs'].header)
                break
        else:
            # else, we have presumably been provided a full undistorted image
            # as a reference, so use it with HSTWCS instead
            wref = wcsutil.HSTWCS(reference)
        refimg.close()

    if isinstance(image, fits.HDUList):
        open_image = False
        filename = image.filename()
        if filename is None:
            filename = 'In-memory HDUList'
        image_update = image.fileinfo(0)['filemode'] == 'update'
    else:
        open_image = True
        filename = image
        image_update = None

    # Now that we are sure we have a good reference WCS to use,
    # continue with the update
    logstr = f"....Updating header for '{filename:s}' ..."
    if verbose:
        print(f"\n{logstr:s}\n")
    else:
        log.info(logstr)

    # reset header WCS keywords to original (OPUS generated) values
    extlist = get_ext_list(image, extname='SCI')
    if extlist:
        if image_update:
            # Create initial WCSCORR extension
            wcscorr.init_wcscorr(image, force=force)
    else:
        extlist = [0]

    # insure that input PRIMARY WCS has been archived before overwriting
    # with new solution
    if open_image:
        fimg = fits.open(image, mode='update', memmap=False)
    else:
        fimg = image

    # Process MEF images...
    for ext in extlist:
        logstr = f"Processing {filename:s}[{ext2str(ext):s}]"
        if verbose:
            print(f"\n{logstr:s}\n")
        else:
            log.info(logstr)

        chip_wcs = wcsutil.HSTWCS(fimg, ext=ext)

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
    """ linearization using 5-point formula for first order derivative """
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
    p = wcsref.wcs_world2pix(wcsim.wcs_pix2world(p, 1), 1).astype(np.longdouble)
    # apply linear fit transformation:
    p = np.dot(f, (p - shift).T).T
    # convert back to image coordinate system:
    p = wcsima.wcs_world2pix(
        wcsref.wcs_pix2world(p.astype(np.float64), 1), 1).astype(np.longdouble)

    # derivative with regard to x:
    u1 = ((p[1] - p[4]) + 8 * (p[3] - p[2])) / (6 * hx)
    # derivative with regard to y:
    u2 = ((p[5] - p[8]) + 8 * (p[7] - p[6])) / (6 * hy)

    return (np.asarray([u1, u2]).T, p[0])


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
        fit = linearfit.buildFitMatrix(rot, scale)

    shift = np.asarray([xsh, ysh]) - np.dot(wcslin.wcs.crpix, fit) + wcslin.wcs.crpix

    fit = np.linalg.inv(fit).T

    cwcs = chip_wcs.deepcopy()
    cd_eye = np.eye(chip_wcs.wcs.cd.shape[0], dtype=np.longdouble)
    zero_shift = np.zeros(2, dtype=np.longdouble)

    naxis1, naxis2 = chip_wcs.pixel_shape

    # estimate precision necessary for iterative processes:
    maxiter = 100
    crpix2corners = np.dstack([i.flatten() for i in np.meshgrid(
        [1, naxis1], [1, naxis2])])[0] - chip_wcs.wcs.crpix
    maxUerr = 1.0e-5 / np.amax(np.linalg.norm(crpix2corners, axis=1))

    # estimate step for numerical differentiation. We need a step
    # large enough to avoid rounding errors and small enough to get a
    # better precision for numerical differentiation.
    # TODO: The logic below should be revised at a later time so that it
    # better takes into account the two competing requirements.
    hx = max(1.0, min(20.0, (chip_wcs.wcs.crpix[0] - 1.0) / 100.0,
                      (naxis1 - chip_wcs.wcs.crpix[0]) / 100.0))
    hy = max(1.0, min(20.0, (chip_wcs.wcs.crpix[1] - 1.0) / 100.0,
                      (naxis2 - chip_wcs.wcs.crpix[1]) / 100.0))

    # compute new CRVAL for the image WCS:
    crpixinref = wcslin.wcs_world2pix(
        chip_wcs.wcs_pix2world([chip_wcs.wcs.crpix], 1), 1)
    crpixinref = np.dot(fit, (crpixinref - shift).T).T
    chip_wcs.wcs.crval = wcslin.wcs_pix2world(crpixinref, 1)[0]
    chip_wcs.wcs.set()

    # initial approximation for CD matrix of the image WCS:
    (U, u) = linearize(cwcs, chip_wcs, wcslin, chip_wcs.wcs.crpix,
                       fit, shift, hx=hx, hy=hy)
    err0 = np.amax(np.abs(U - cd_eye)).astype(np.float64)
    chip_wcs.wcs.cd = np.dot(chip_wcs.wcs.cd.astype(np.longdouble), U).astype(np.float64)
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
        err = np.amax(np.abs(U - cd_eye)).astype(np.float64)
        if err > err0:
            break
        chip_wcs.wcs.cd = np.dot(chip_wcs.wcs.cd, U).astype(np.float64)
        chip_wcs.wcs.set()
        if err < maxUerr:
            break
        err0 = err

    if xrms is not None:
        chip_wcs.wcs.crder = np.array([xrms, yrms])


###
### Header keyword prefix related archive functions
###
def update_wcs(image, extnum, new_wcs, wcsname="", reusename=False, verbose=False):
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

    if isinstance(image, fits.HDUList):
        close_file = False
        fname = image.filename()
    else:
        fname = image
        image = fits.open(image, mode='update', memmap=False)
        close_file = True

    hdr = image[extnum].header

    # Name of the updated primary WCS
    if util.is_blank(wcsname):
        wcsname = 'TWEAK'

    # Auto-rename old primary WCS when archiving it if an alternate WCS with
    # the same name already exists:
    if 'WCSNAME' in hdr:
        pri_wcsname = hdr['WCSNAME']
        pri_wcsname_u = pri_wcsname.upper()
    else:
        pri_wcsname = None
        pri_wcsname_u = None

    if pri_wcsname_u == wcsname.upper():
        if not reusename:
            raise ValueError(
                f"WCSNAME '{wcsname}' already present in '{fname}'. A unique "
                "value for the 'wcsname' parameter needs to be specified."
            )
    elif close_file or image.fileinfo(0) is None or image.fileinfo(0)['filemode'] == 'update':
        altwcs.archive_wcs(image, [extnum], wcsname=pri_wcsname, mode=altwcs.ArchiveMode.AUTO_RENAME)

    # Update Primary WCS:
    try:
        logstr = f'Updating header for {image.filename()}[{extnum}]'
        if verbose:
            print(logstr)
            log.info('    with WCS of')
            new_wcs.printwcs()
            print("WCSNAME  : ", wcsname)
        else:
            log.info(logstr)

        wcs_hdr = new_wcs.wcs2header(idc2hdr=new_wcs.idcscale is not None, relax=True)
        wcs_hdr.set('WCSNAME', wcsname, before=0)
        wcs_hdr.set('WCSTYPE', interpret_wcsname_type(wcsname), after=0)
        wcs_hdr.set('ORIENTAT', new_wcs.orientat, after=len(wcs_hdr))
        hdr.update(wcs_hdr)
        util.updateNEXTENDKw(image)

    finally:
        if close_file:
            image.close()

def interpret_wcsname_type(wcsname):
    """Interpret WCSNAME as a standardized human-understandable description """
    wcstype = ''
    fit_terms = {'REL': 'relatively aligned to {}',
                 'IMG': 'aligned image-by-image to {}',
                 'EVM': 'aligned by visit-exposures to {}',
                 'SVM': 'aligned by visit to {}'}
    post_fit = 'a posteriori solution '
    default_fit = 'a priori solution based on {}'
    base_terms = {'IDC': 'undistorted ',
                  'OPU': 'pipeline default '}
    no_fit = 'not aligned'

    if wcsname is None:
        return no_fit

    wcsname = wcsname.upper()  # make this comparison case-insensitive
    wcsname_list = wcsname.split('-')
    wcsname_term = wcsname_list[0][:3]
    if wcsname_term not in base_terms:
        wcstype = 'user-defined'  # Set to user-defined default
        return wcstype

    # Interpret base terms
    wcstype += base_terms[wcsname_term]

    # Interpret fit term (if any)
    fit_term = wcsname_list[1] if len(wcsname_list) > 1 else None

    if len(wcsname_list) == 1:
        wcstype += no_fit
    else:
        if 'FIT' not in fit_term:
            wcstype += default_fit.format(fit_term)
        else:
            wcstype += post_fit
            postfit_type = fit_term.split('_')
            wcstype += fit_terms[postfit_type[1]].format(postfit_type[2])
    return wcstype

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
    wnames = list(altwcs.wcsnames(fimg, ext=extnum).values())
    if wcsname not in wnames:
        uniqname = wcsname
    else:
        # setup pattern to match
        rpatt = re.compile(wcsname + '_\d')
        index = 0
        for wname in wnames:
            rmatch = rpatt.match(wname)
            if rmatch:
                # get index
                n = int(wname[wname.rfind('_') + 1:])
                if n > index: index = 1
        index += 1  # for use with new name
        uniqname = "%s_%d" % (wcsname, index)
    return uniqname
