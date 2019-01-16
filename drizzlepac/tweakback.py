"""
`tweakback` - propagate the "tweaked" solutions back to the original
input files.

Version 0.4.0 - replaced previous algorithm that used fitting of WCS
footprints to reconstruct the transformation that was applied to the
old drizzled image (to align it with another image) to obtain the new
drizzled image WCS with an algorithm that is based on linearization of
the exact compound operator that transforms current image coordinates to
the "aligned" (to the new drizzled WCS) image coordinates.

:Authors: Warren Hack, Mihai Cara

:License: :doc:`LICENSE`

"""
import os

import numpy as np
from astropy.io import fits

from stwcs import wcsutil
from stsci.tools import parseinput, logutil
from stsci.skypac.utils import get_ext_list, ext2str

from . import updatehdr
from . import linearfit
from . import util


__taskname__ = 'tweakback' # unless someone comes up with anything better

# This is specifically NOT intended to match the package-wide version information.
__version__ = '0.4.0'
__version_date__ = '14-Oct-2014'


log = logutil.create_logger(__name__, level=logutil.logging.NOTSET)

if hasattr(np, 'float128'):
    ndfloat128 = np.float128
elif hasattr(np, 'float96'):
    ndfloat128 = np.float96
else:
    ndfloat128 = np.float64


#### Primary function
def tweakback(drzfile, input=None,  origwcs = None,
                newname = None, wcsname = None,
                extname='SCI', force=False, verbose=False):
    """
    Apply WCS solution recorded in drizzled file to distorted input images
    (``_flt.fits`` files) used to create the drizzled file.  This task relies on
    the original WCS and updated WCS to be recorded in the drizzled image's
    header as the last 2 alternate WCSs.

    Parameters
    ----------
    drzfile : str (Default = '')
        filename of undistorted image which contains the new WCS
        and WCS prior to being updated

    newname : str (Default = None)
        Value of ``WCSNAME`` to be used to label the updated solution in the
        output (eq., ``_flt.fits``) files.  If left blank or None, it will
        default to using the current ``WCSNAME`` value from the input drzfile.

    input : str (Default = '')
        filenames of distorted images to be updated using new WCS
        from 'drzfile'.  These can be provided either as an ``@-file``,
        a comma-separated list of filenames or using wildcards.

        .. note:: A blank value will indicate that the task should derive the
           filenames from the 'drzfile' itself, if possible. The filenames will be
           derived from the ``D*DATA`` keywords written out by
           ``AstroDrizzle``. If they can not be found, the task will quit.

    origwcs : str (Default = None)
        Value of ``WCSNAME`` keyword prior to the drzfile image being updated
        by ``TweakReg``.  If left blank or None, it will default to using the
        second to last ``WCSNAME*`` keyword value found in the header.

    wcsname : str (Default = None)
        Value of WCSNAME for updated solution written out by ``TweakReg`` as
        specified by the `wcsname` parameter from ``TweakReg``.  If this is
        left blank or `None`, it will default to the current ``WCSNAME``
        value from the input drzfile.

    extname : str (Default = 'SCI')
        Name of extension in `input` files to be updated with new WCS

    force : bool  (Default = False)
        This parameters specified whether or not to force an update of the WCS
        even though WCS already exists with this solution or `wcsname`?

    verbose : bool (Default = False)
        This parameter specifies whether or not to print out additional
        messages during processing.


    Notes
    -----
    The algorithm used by this function is based on linearization of
    the exact compound operator that converts input image coordinates
    to the coordinates (in the input image) that would result in
    alignment with the new drizzled image WCS.

    If no input distorted files are specified as input, this task will attempt
    to generate the list of filenames from the drizzled input file's own
    header.

    EXAMPLES
    --------
    An image named ``acswfc_mos2_drz.fits`` was created from 4 images using
    astrodrizzle. This drizzled image was then aligned to another image using
    tweakreg and the header was updated using the ``WCSNAME`` = ``TWEAK_DRZ``.
    The new WCS can then be used to update each of the 4 images that were
    combined to make up this drizzled image using:

    >>> from drizzlepac import tweakback
    >>> tweakback.tweakback('acswfc_mos2_drz.fits')

    If the same WCS should be applied to a specific set of images, those images
    can be updated using:

    >>> tweakback.tweakback('acswfc_mos2_drz.fits',
    ...                     input='img_mos2a_flt.fits,img_mos2e_flt.fits')

    See Also
    --------
    stwcs.wcsutil.altwcs: Alternate WCS implementation

    """
    print("TweakBack Version {:s}({:s}) started at: {:s}\n"
          .format(__version__,__version_date__,util._ptime()[0]))

    # Interpret input list/string into list of filename(s)
    fltfiles = parseinput.parseinput(input)[0]

    if fltfiles is None or len(fltfiles) == 0:
        # try to extract the filenames from the drizzled file's header
        fltfiles = extract_input_filenames(drzfile)
        if fltfiles is None:
            print('*'*60)
            print('*')
            print('* ERROR:')
            print('*    No input filenames found! ')
            print('*    Please specify "fltfiles" or insure that input drizzled')
            print('*    image contains D*DATA keywords. ')
            print('*')
            print('*'*60)
            raise ValueError

    if not isinstance(fltfiles,list):
        fltfiles = [fltfiles]

    sciext = determine_extnum(drzfile, extname='SCI')
    scihdr = fits.getheader(drzfile, ext=sciext, memmap=False)

    ### Step 1: Read in updated and original WCS solutions
    # determine keys for all alternate WCS solutions in drizzled image header
    wkeys = wcsutil.altwcs.wcskeys(drzfile, ext=sciext)
    wnames = wcsutil.altwcs.wcsnames(drzfile, ext=sciext)
    if not util.is_blank(newname):
        final_name = newname
    else:
        final_name = wnames[wkeys[-1]]

    # Read in HSTWCS objects for final,updated WCS and previous WCS from
    # from drizzled image header
    # The final solution also serves as reference WCS when using updatehdr
    if not util.is_blank(wcsname):
        for k in wnames:
            if wnames[k] == wcsname:
                wcskey = k
                break
    else:
        wcskey = wkeys[-1]
    final_wcs = wcsutil.HSTWCS(drzfile, ext=sciext, wcskey=wkeys[-1])

    if not util.is_blank(origwcs):
        for k in wnames:
            if wnames[k] == origwcs:
                orig_wcskey = k
                orig_wcsname = origwcs
                break
    else:
        orig_wcsname,orig_wcskey = determine_orig_wcsname(scihdr,wnames,wkeys)

    orig_wcs = wcsutil.HSTWCS(drzfile,ext=sciext,wcskey=orig_wcskey)

    # read in RMS values reported for new solution
    crderr1kw = 'CRDER1'+wkeys[-1]
    crderr2kw = 'CRDER2'+wkeys[-1]

    if crderr1kw in scihdr:
        crderr1 = fits.getval(drzfile, crderr1kw, ext=sciext, memmap=False)
    else:
        crderr1 = 0.0

    if crderr2kw in scihdr:
        crderr2 = fits.getval(drzfile, crderr2kw, ext=sciext, memmap=False)
    else:
        crderr2 = 0.0
    del scihdr

    ### Step 2: Apply solution to input file headers
    for fname in fltfiles:
        logstr = "....Updating header for {:s}...".format(fname)
        if verbose:
            print("\n{:s}\n".format(logstr))
        else:
            log.info(logstr)

        # reset header WCS keywords to original (OPUS generated) values
        imhdulist = fits.open(fname, mode='update', memmap=False)
        extlist = get_ext_list(imhdulist, extname='SCI')
        if not extlist:
            extlist = [0]

        # insure that input PRIMARY WCS has been archived before overwriting
        # with new solution
        wcsutil.altwcs.archiveWCS(imhdulist, extlist, reusekey=True)

        # Process MEF images...
        for ext in extlist:
            logstr = "Processing {:s}[{:s}]".format(imhdulist.filename(),
                                                    ext2str(ext))
            if verbose:
                print("\n{:s}\n".format(logstr))
            else:
                log.info(logstr)
            chip_wcs = wcsutil.HSTWCS(imhdulist, ext=ext)

            update_chip_wcs(chip_wcs, orig_wcs, final_wcs,
                            xrms=crderr1, yrms = crderr2)

            # Update FITS file with newly updated WCS for this chip
            extnum = imhdulist.index(imhdulist[ext])
            updatehdr.update_wcs(imhdulist, extnum, chip_wcs,
                                 wcsname=final_name, reusename=False,
                                 verbose=verbose)

        imhdulist.close()


def linearize(wcsim, wcsima, wcs_olddrz, wcs_newdrz, imcrpix, hx=1.0, hy=1.0):
    # linearization using 5-point formula for first order derivative
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
    # convert image coordinates to old drizzled image coordinates:
    p = wcs_olddrz.wcs_world2pix(wcsim.wcs_pix2world(p, 1), 1)
    # convert to sky coordinates using the new drizzled image's WCS:
    p = wcs_newdrz.wcs_pix2world(p, 1)
    # convert back to image coordinate system using partially (CRVAL only)
    # aligned image's WCS:
    p = wcsima.wcs_world2pix(p, 1).astype(ndfloat128)

    # derivative with regard to x:
    u1 = ((p[1] - p[4]) + 8 * (p[3] - p[2])) / (6*hx)
    # derivative with regard to y:
    u2 = ((p[5] - p[8]) + 8 * (p[7] - p[6])) / (6*hy)

    return (np.asarray([u1, u2]).T, p[0])


def update_chip_wcs(chip_wcs, drz_old_wcs, drz_new_wcs,
                    xrms=None, yrms=None):
    cd_eye = np.eye(chip_wcs.wcs.cd.shape[0])

    # estimate precision necessary for iterative processes:
    naxis1, naxis2 = chip_wcs.pixel_shape
    maxiter = 100
    crpix2corners = np.dstack([i.flatten() for i in np.meshgrid(
        [1, naxis1], [1, naxis2])])[0] - chip_wcs.wcs.crpix
    maxUerr = 1.0e-5 / np.amax(np.linalg.norm(crpix2corners, axis=1))

    # estimate step for numerical differentiation. We need a step
    # large enough to avoid rounding errors and small enough to get a
    # better precision for numerical differentiation.
    # TODO: The logic below should be revised at a later time so that it
    # better takes into account the two competing requirements.
    hx = max(1.0, min(20.0, (chip_wcs.wcs.crpix[0] - 1.0)/100.0,
                      (naxis1 - chip_wcs.wcs.crpix[0])/100.0))
    hy = max(1.0, min(20.0, (chip_wcs.wcs.crpix[1] - 1.0)/100.0,
                      (naxis2 - chip_wcs.wcs.crpix[1])/100.0))

    # compute new CRVAL for the image WCS:
    chip_wcs_orig = chip_wcs.deepcopy()
    crpix_in_old_drz = drz_old_wcs.wcs_world2pix([chip_wcs.wcs.crval], 1)
    chip_wcs.wcs.crval = drz_new_wcs.wcs_pix2world(crpix_in_old_drz, 1)[0]
    chip_wcs.wcs.set()

    # initial approximation for CD matrix of the image WCS:
    (U, u) = linearize(chip_wcs_orig, chip_wcs, drz_old_wcs, drz_new_wcs,
                       chip_wcs_orig.wcs.crpix, hx=hx, hy=hy)
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
        (U, u) = linearize(chip_wcs_orig, chip_wcs, drz_old_wcs, drz_new_wcs,
                           chip_wcs_orig.wcs.crpix, hx=hx, hy=hy)
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


#### TEAL Interfaces to run this task


def run(configobj):
    # Interpret user-input from TEAL GUI and call function
    tweakback(configobj['drzfile'], newname = configobj['newname'],
            input=configobj['input'], origwcs = configobj['origwcs'],
            wcsname = configobj['wcsname'],
            extname=configobj['extname'],verbose=configobj['verbose'],
            force=configobj['force'])


#### Utility functions
#
def extract_input_filenames(drzfile):
    """
    Generate a list of filenames from a drizzled image's header
    """
    data_kws = fits.getval(drzfile, 'd*data', ext=0, memmap=False)
    if len(data_kws) == 0:
        return None
    fnames = []
    for kw in data_kws.cards:
        f = kw.value.split('[')[0]
        if f not in fnames:
            fnames.append(f)

    return fnames

def determine_extnum(drzfile, extname='SCI'):
    # Determine what kind of drizzled file input has been provided: MEF or single
    hdulist = fits.open(drzfile, memmap=False)
    numext = len(hdulist)
    sciext = 0
    for e,i in zip(hdulist,list(range(numext))):
        if 'extname' in e.header and e.header['extname'] == extname:
            sciext = i
            break
    hdulist.close()

    return sciext

def determine_orig_wcsname(header, wnames, wkeys):
    """
    Determine the name of the original, unmodified WCS solution
    """
    orig_wcsname = None
    orig_key = None
    if orig_wcsname is None:
        for k,w in wnames.items():
            if w[:4] == 'IDC_':
                orig_wcsname = w
                orig_key = k
                break
    if orig_wcsname is None:
        # No IDC_ wcsname found... revert to second to last if available
        if len(wnames) > 1:
            orig_key = wkeys[-2]
            orig_wcsname = wnames[orig_key]
    return orig_wcsname,orig_key


def help(file=None):
    """
    Print out syntax help for running astrodrizzle

    Parameters
    ----------
    file : str (Default = None)
        If given, write out help to the filename specified by this parameter
        Any previously existing file with this name will be deleted before
        writing out the help.

    """
    helpstr = getHelpAsString(docstring=True, show_ver = True)
    if file is None:
        print(helpstr)
    else:
        if os.path.exists(file): os.remove(file)
        f = open(file, mode = 'w')
        f.write(helpstr)
        f.close()


def getHelpAsString(docstring = False, show_ver = True):
    """
    return useful help from a file in the script directory called
    __taskname__.help

    """
    install_dir = os.path.dirname(__file__)
    taskname = util.base_taskname(__taskname__, '')
    htmlfile = os.path.join(install_dir, 'htmlhelp', taskname + '.html')
    helpfile = os.path.join(install_dir, taskname + '.help')

    if docstring or (not docstring and not os.path.exists(htmlfile)):
        if show_ver:
            helpString = os.linesep + \
                ' '.join([__taskname__, 'Version', __version__,
                ' updated on ', __version_date__]) + 2*os.linesep
        else:
            helpString = ''
        if os.path.exists(helpfile):
            helpString += teal.getHelpFileAsString(taskname, __file__)
        else:
            if __doc__ is not None:
                helpString += tweakback.__doc__ + os.linesep
    else:
        helpString = 'file://' + htmlfile

    return helpString
