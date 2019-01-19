""" pixtosky - A module to perform coordinate transformation from pixel to sky coordinates.

    :Authors: Warren Hack

    :License: :doc:`LICENSE`

    PARAMETERS
    ----------
    input : str
        full filename with path of input image, an extension name ['sci',1] should be
        provided if input is a multi-extension FITS file

    Optional Parameters
    -------------------
    x : float or list or array, optional
        X position from input image for a single or multiple sources
    y : float or list or array, optional
        Y position from input image for a single or multiple sources
    coords : str, deprecated
        [DEPRECATED] full filename with path of file with x,y coordinates
        Filename given here will be *ignored* if a file has been specified
        in `coordfile` parameter.
    coordfile : str, optional
        full filename with path of file with x,y coordinates
    colnames : str, optional
        comma separated list of column names or list of column name strings
        from 'coordfile' files containing x,y coordinates, respectively.
        This parameter will default to first two columns if None are specified.
        Column names for ASCII files will use 'c1','c2',... convention.
        Valid syntax: ['c1','c3'] or 'c1,c3'
    separator : str, optional
        non-blank separator used as the column delimiter in the coordfile file
    hms : bool, optional
        Produce output in HH:MM:SS.S format instead of decimal degrees? (default: False)
    precision : int, optional
        Number of floating-point digits in output values
    output : str, optional
        Name of output file with results, if desired
    verbose : bool
        Print out full list of transformation results (default: False)

    RETURNS
    -------
    ra : float or array
        Right Ascension of pixel. If more than 1 input value, then it will be a
        numpy array.
    dec : float or array
        Declination of pixel. If more than 1 input value, then it will be a
        numpy array.

    NOTES
    -----
    This task performs a full distortion-correction coordinate transformation
    based on all WCS keywords and any recognized distortion keywords from the
    input image header. The transformation recognizes the conventions for
    describing distortion implemented as part of the SIP and Paper IV
    conventions used with ``AstroDrizzle``. Input images can be updated to use
    these conventions through the use of the ``updatewcs`` module the ``STWCS``
    package.


    See Also
    --------
    `stwcs`

    EXAMPLES
    --------
    1. The following command will transform the position 256,256 into a
       position on the sky for the image 'input_flt.fits[sci,1]' using::

       >>> from drizzlepac import pixtosky
       >>> r,d = pixtosky.xy2rd("input_file_flt.fits[sci,1]", 256,256)


    2. The set of X,Y positions from 'input_flt.fits[sci,1]' stored as
       the 3rd and 4th columns from the ASCII file 'xy_sci1.dat'
       will be transformed and written out to 'radec_sci1.dat' using::

       >>> from drizzlepac import pixtosky
       >>> r,d = pixtosky.xy2rd("input_flt.fits[sci,1]", coordfile='xy_sci1.dat',
       ...                      colnames=['c3','c4'], output="radec_sci1.dat")

"""
import os,copy
import warnings
import numpy as np

from stsci.tools import fileutil, teal
from . import util
from . import wcs_functions
import stwcs
from stwcs import distortion, wcsutil

# This is specifically NOT intended to match the package-wide version information.
__version__ = '0.1'
__version_date__ = '20-Jan-2011'

__taskname__ = 'pixtosky'

blank_list = [None, '', ' ']

def xy2rd(input,x=None,y=None,coords=None, coordfile=None,colnames=None,separator=None,
            hms=True, precision=6,output=None,verbose=True):
    """ Primary interface to perform coordinate transformations from
        pixel to sky coordinates using STWCS and full distortion models
        read from the input image header.
    """
    single_coord = False
    # Only use value provided in `coords` if nothing has been specified for coordfile
    if coords is not None and coordfile is None:
        coordfile = coords
        warnings.simplefilter('always',DeprecationWarning)
        warnings.warn("Please update calling code to pass in `coordfile` instead of `coords`.",
            category=DeprecationWarning)
        warnings.simplefilter('default',DeprecationWarning)

    if coordfile is not None:
        if colnames in blank_list:
            colnames = ['c1','c2']
        # Determine columns which contain pixel positions
        cols = util.parse_colnames(colnames,coordfile)
        # read in columns from input coordinates file
        xyvals = np.loadtxt(coordfile,usecols=cols,delimiter=separator)
        if xyvals.ndim == 1:  # only 1 entry in coordfile
            xlist = [xyvals[0].copy()]
            ylist = [xyvals[1].copy()]
        else:
            xlist = xyvals[:,0].copy()
            ylist = xyvals[:,1].copy()
        del xyvals
    else:
        if isinstance(x, np.ndarray):
            xlist = x.tolist()
            ylist = y.tolist()
        elif not isinstance(x,list):
            xlist = [x]
            ylist = [y]
            single_coord = True
        else:
            xlist = x
            ylist = y

    # start by reading in WCS+distortion info for input image
    inwcs = wcsutil.HSTWCS(input)
    if inwcs.wcs.is_unity():
        print("####\nNo valid WCS found in {}.\n  Results may be invalid.\n####\n".format(input))

    # Now, convert pixel coordinates into sky coordinates
    dra,ddec = inwcs.all_pix2world(xlist,ylist,1)

    # convert to HH:MM:SS.S format, if specified
    if hms:
        ra,dec = wcs_functions.ddtohms(dra,ddec,precision=precision)
        rastr = ra
        decstr = dec
    else:
        # add formatting based on precision here...
        rastr = []
        decstr = []
        fmt = "%."+repr(precision)+"f"
        for r,d in zip(dra,ddec):
            rastr.append(fmt%r)
            decstr.append(fmt%d)

        ra = dra
        dec = ddec

    if verbose or (not verbose and util.is_blank(output)):
        print('# Coordinate transformations for ',input)
        print('# X      Y         RA             Dec\n')
        for x,y,r,d in zip(xlist,ylist,rastr,decstr):
            print("%.4f  %.4f    %s  %s"%(x,y,r,d))

    # Create output file, if specified
    if output:
        f = open(output,mode='w')
        f.write("# Coordinates converted from %s\n"%input)
        for r,d in zip(rastr,decstr):
            f.write('%s    %s\n'%(r,d))
        f.close()
        print('Wrote out results to: ',output)

    if single_coord:
        ra = ra[0]
        dec = dec[0]

    return ra,dec

#--------------------------
# TEAL Interface functions
#--------------------------
def run(configObj):

    if 'coords' in configObj:
        coords = util.check_blank(configObj['coords'])
    else:
        coords = None
    coordfile = util.check_blank(configObj['coordfile'])
    colnames = util.check_blank(configObj['colnames'])
    sep = util.check_blank(configObj['separator'])
    outfile = util.check_blank(configObj['output'])

    xy2rd(configObj['input'],
            x = configObj['x'], y = configObj['y'], coords=coords,
            coordfile = coordfile, colnames = colnames,
            separator= sep, hms = configObj['hms'], precision= configObj['precision'],
            output= outfile, verbose = configObj['verbose'])


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
                helpString += __doc__ + os.linesep
    else:
        helpString = 'file://' + htmlfile

    return helpString


__doc__ = getHelpAsString(docstring = True, show_ver = False)
