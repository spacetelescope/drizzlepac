import os

import numpy as np
import pyfits

from stwcs import wcsutil
from stsci.tools import parseinput

from . import updatehdr
from . import linearfit
from . import util


__taskname__ = 'tweakback' # unless someone comes up with anything better

# This is specifically NOT intended to match the package-wide version information.
__version__ = '0.3.0'
__vdate__ = '1-Feb-2012'

#### Primary function
def tweakback(drzfile, input=None,  origwcs = None,
                newname = None, wcsname = None,
                extname='SCI', force=False, verbose=False):
    """
    Apply WCS solution recorded in drizzled file to distorted input images
    (_flt.fits files) used to create the drizzled file.  This task relies on
    the original WCS and updated WCS to be recorded in the drizzled image's
    header as the last 2 alternate WCSs.

    Parameters
    ----------
    drzfile : str (Default: '')
        filename of undistorted image which contains the new WCS
        and WCS prior to being updated

    newname : str (Default: None)
        Value of WCSNAME to be used to label the updated solution in the
        output (eq., _flt.fits) files.  If left blank or None, it will
        default to using the current WCSNAME value from the input drzfile.

    input : str (Default: '')
        filenames of distorted images to be updated using new WCS
        from 'drzfile'.  These can be provided either as an '@'-file,
        a comma-separated list of filenames or using wildcards.

        .. note:: A blank value will indicate that the task should derive the
        filenames from the 'drzfile' itself, if possible. The filenames will be
        derived from the D*DATA keywords written out by astrodrizzle
        (or MultiDrizzle or drizzle).
        If they can not be found, the task will quit.

    origwcs : str (Default: None)
        Value of WCSNAME keyword prior to the drzfile image being updated
        by tweakreg.  If left blank or None, it will default to using the
        second to last WCSNAME* keyword value found in the header.

    wcsname : str (Default: None)
        Value of WCSNAME for updated solution written out by 'tweakreg' as
        specified by the 'wcsname' parameter from 'tweakreg'.  If this is
        left blank or None, it will default to the current WCSNAME value from the
        input drzfile.

    extname : str (Default: 'SCI')
        Name of extension in 'input' files to be updated with new WCS

    force: bool  (Default: False)
        This parameters specified whether or not to force an update of the WCS
        even though WCS already exists with this solution or wcsname?

    verbose: bool (Default: False)
        This parameter specifies whether or not to print out additional
        messages during processing.


    Notes
    -----
    The algorithm used by this function follows these steps::

    0. Verify or determine list of distorted image's that need to be updated
        with final solution from drizzled image
    1. Read in HSTWCS objects for last 2 alternate WCS solutions
    2. Generate footprints using .calFootprint() for each WCS
    3. Create pixel positions for corners of each WCS's footprint
        by running the .wcs_sky2pix() method for the last (updated) WCS
    4. Perform linear 'rscale' fit between the 2 sets of X,Y coords
    5. Update each input image WCS with fit using 'updatehdr_with_shift()'

    If no input distorted files are specified as input, this task will attempt
    to generate the list of filenames from the drizzled input file's own
    header.

    EXAMPLES
    --------
    An image named 'acswfc_mos2_drz.fits' was created from 4 images using
    astrodrizzle. This drizzled image was then aligned to another image using
    tweakreg and the header was updated using the WCSNAME = TWEAK_DRZ.
    The new WCS can then be used to update each of the 4 images that were
    combined to make up this drizzled image using::

    >>> from drizzlepac import tweakback
    >>> tweakback.tweakback('acswfc_mos2_drz.fits')

    If the same WCS should be applied to a specific set of images, those images
    can be updated using::

    >>> tweakback.tweakback('acswfc_mos2_drz.fits',
            input='img_mos2a_flt.fits,img_mos2e_flt.fits')

    See Also
    --------
    stwcs.wcsutil.altwcs: Alternate WCS implementation

    """
    print 'TweakBack Version %s(%s) started at: %s \n'%(
                    __version__,__vdate__,util._ptime()[0])
    
    # Interpret input list/string into list of filename(s)
    fltfiles = parseinput.parseinput(input)[0]

    if fltfiles is None or len(fltfiles) == 0:
        # try to extract the filenames from the drizzled file's header
        fltfiles = extract_input_filenames(drzfile)
        if fltfiles is None:
            print '*'*60
            print '*'
            print '* ERROR:'
            print '*    No input filenames found! '
            print '*    Please specify "fltfiles" or insure that input drizzled'
            print '*    image contains D*DATA keywords. '
            print '*'
            print '*'*60
            raise ValueError

    if not isinstance(fltfiles,list):
        fltfiles = [fltfiles]

    sciext = determine_extnum(drzfile, extname='SCI')
    scihdr = pyfits.getheader(drzfile,ext=sciext)

    ### Step 1: Read in updated and original WCS solutions
    # determine keys for all alternate WCS solutions in drizzled image header
    wkeys = wcsutil.altwcs.wcskeys(drzfile,ext=sciext)
    wnames = wcsutil.altwcs.wcsnames(drzfile,ext=sciext)
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
    final_wcs = wcsutil.HSTWCS(drzfile,ext=sciext,wcskey=wkeys[-1])

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
        crderr1 = pyfits.getval(drzfile,crderr1kw,ext=sciext)
    else:
        crderr1 = 0.0

    if crderr2kw in scihdr:
        crderr2 = pyfits.getval(drzfile,crderr2kw,ext=sciext)
    else:
        crderr2 = 0.0
    del scihdr
    ### Step 2: Generate footprints for each WCS
    final_fp = final_wcs.calcFootprint()
    orig_fp = orig_wcs.calcFootprint()

    ### Step 3: Create pixel positions in final WCS for each footprint
    final_xy_fp = final_wcs.wcs_sky2pix(final_fp,1)
    orig_xy_fp = final_wcs.wcs_sky2pix(orig_fp,1)

    ### Step 4: Perform fit between footprint X,Y positions
    wfit = linearfit.iter_fit_all(orig_xy_fp,final_xy_fp,range(4),range(4),
                                mode='rscale',nclip=0, verbose=verbose,
                                center=final_wcs.wcs.crpix)

    ### Step 5: Apply solution to input file headers
    for fname in fltfiles:
        updatehdr.updatewcs_with_shift(fname,final_wcs,wcsname=final_name,
                        rot=wfit['rot'],scale=wfit['scale'][0],
                        xsh=wfit['offset'][0],ysh=wfit['offset'][1],
                        fit=wfit['fit_matrix'],
                        xrms=crderr1, yrms = crderr2,
                        verbose=verbose,force=force,sciext=extname)

#### TEAL Interfaces to run this task
def getHelpAsString(docstring=False):
    """
    return useful help from a file in the script directory called __taskname__.help
    """
    install_dir = os.path.dirname(__file__)
    htmlfile = os.path.join(install_dir,'htmlhelp',__taskname__+'.html')
    helpfile = os.path.join(install_dir,__taskname__+'.help')
    if docstring or (not docstring and not os.path.exists(htmlfile)):
        helpString = __taskname__+' Version '+__version__+' updated on '+__vdate__+'\n\n'
        if os.path.exists(helpfile):
            helpString += teal.getHelpFileAsString(__taskname__,__file__)
        else:
            helpString += tweakback.__doc__
    else:
        helpString = 'file://'+htmlfile

    return helpString

def run(configobj):
    # Interpret user-input from TEAL GUI and call function
    tweakback(configobj['drzfile'], newname = configobj['newname'],
            input=configobj['input'], origwcs = configobj['origwcs'],
            wcsname = configobj['wcsname'],
            extname=configobj['extname'],verbose=configobj['verbose'],
            force=configobj['force'])

def help(file=None):
    """
    Print out syntax help for running tweakback

    Parameters
    ----------
    file : str (Default = None)
        If given, write out help to the filename specified by this parameter
        Any previously existing file with this name will be deleted before
        writing out the help.

    """
    helpstr = getHelpAsString(docstring=True)
    if file is None:
        print helpstr
    else:
        if os.path.exists(file): os.remove(file)
        f = open(file,mode='w')
        f.write(helpstr)
        f.close()
#
#### Utility functions
#
def extract_input_filenames(drzfile):
    """
    Generate a list of filenames from a drizzled image's header
    """
    data_kws = pyfits.getval(drzfile,'d*data',ext=0)
    if len(data_kws) == 0:
        return None
    fnames = []
    for kw in data_kws.ascard:
        f = kw.value.split('[')[0]
        if f not in fnames:
            fnames.append(f)

    return fnames

def determine_extnum(drzfile, extname='SCI'):
    # Determine what kind of drizzled file input has been provided: MEF or single
    hdulist = pyfits.open(drzfile)
    numext = len(hdulist)
    sciext = 0
    for e,i in zip(hdulist,range(numext)):
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
