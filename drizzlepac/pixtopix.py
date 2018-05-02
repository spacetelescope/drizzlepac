""" pixtosky - A module to perform coordinate transformation from pixel coordinates
                in one image to pixel coordinates in another frame

    :Authors: Warren Hack

    :License: :doc:`LICENSE`

    PARAMETERS
    ----------
    inimage : str
        full filename with path of input image, an extension name ['sci',1] should be
        provided if input is a multi-extension FITS file
    outimage : str, optional
        full filename with path of output image, an extension name ['sci',1] should be
        provided if output is a multi-extension FITS file. If no image gets
        specified, the input image will be used to generate a default output
        WCS using stwcs.distortion.util.output_wcs().
    direction : str
        Direction of transform (forward or backward). The 'forward' transform
        takes the pixel positions (assumed to be from the 'input' image) and determines
        their position in the 'output' image. The 'backward' transform converts
        the pixel positions (assumed to be from the 'output' image) into pixel
        positions in the 'input' image.


    Optional Parameters
    -------------------
    x : float, optional
        X position from image
    y : float, optional
        Y position from image
    coords : str, deprecated
        [DEPRECATED] full filename with path of file with x,y coordinates
        Filename given here will be *ignored* if a file has been specified
        in `coordfile` parameter.
    coordfile : str, optional
        full filename with path of file with starting x,y coordinates
    colnames : str, optional
        comma separated list of column names from 'coordfile' files
        containing x,y coordinates, respectively. Will default to
        first two columns if None are specified. Column names for ASCII
        files will use 'c1','c2',... convention.
    separator : str, optional
        non-blank separator used as the column delimiter in the coordfile file
    precision : int, optional
        Number of floating-point digits in output values
    output : str, optional
        Name of output file with results, if desired
    verbose : bool
        Print out full list of transformation results (default: False)

    RETURNS
    -------
    outx : float
        X position of transformed pixel. If more than 1 input value, then it
        will be a numpy array.
    outy : float
        Y position of transformed pixel. If more than 1 input value, then it
        will be a numpy array.

    NOTES
    -----
    This module performs a full distortion-corrected coordinate transformation
    based on all WCS keywords and any recognized distortion keywords from the
    input image header.

    Usage
    -----
    It can be called from within Python using the syntax::

        >>> from drizzlepac import pixtopix
        >>> outx,outy = pixtopix.tran("input_flt.fits[sci,1]",
                        "output_drz.fits[sci,1],"forward",100,100)

    EXAMPLES
    --------

    1. The following command will transform the position 256,256 from
        'input_flt.fits[sci,1]' into a position on the output image
        'output_drz.fits[sci,1]' using::

            >>> from drizzlepac import pixtopix
            >>> outx,outy = pixtopix.tran("input_file_flt.fits[sci,1]",
                        "output_drz.fits[sci,1],"forward", 256,256)


    2. The set of X,Y positions from 'output_drz.fits[sci,1]' stored as
        the 3rd and 4th columns from the ASCII file 'xy_sci1.dat'
        will be transformed into pixel positions from 'input_flt.fits[sci,1]'
        and written out to 'xy_flt1.dat' using::

            >>> from drizzlepac import pixtopix
            >>> x,y = pixtopix.tran("input_flt.fits[sci,1]", "output_drz.fits[sci,1]",
                    "backward", coordfile='xy_sci1.dat', colnames=['c3','c4'],
                    output="xy_flt1.dat")

"""
from __future__ import absolute_import, division, print_function # confidence medium

import os,copy
import warnings
import numpy as np

from stsci.tools import fileutil, teal
from . import wcs_functions
from . import util
from stwcs import wcsutil, distortion

# This is specifically NOT intended to match the package-wide version information.
__version__ = '0.1'
__version_date__ = '1-Mar-2011'

__taskname__ = 'pixtopix'

def tran(inimage,outimage,direction='forward',x=None,y=None,
        coords=None, coordfile=None,colnames=None,separator=None,
        precision=6, output=None,verbose=True):
    """ Primary interface to perform coordinate transformations in pixel
        coordinates between 2 images using STWCS and full distortion models
        read from each image's header.
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
        if colnames in util.blank_list:
            colnames = ['c1','c2']
        # Determine columns which contain pixel positions
        cols = util.parse_colnames(colnames,coordfile)
        # read in columns from input coordinates file
        xyvals = np.loadtxt(coordfile,usecols=cols,delimiter=separator)
        if xyvals.ndim == 1: # only 1 entry in coordfile
            xlist = [xyvals[0].copy()]
            ylist = [xyvals[1].copy()]
        else:
            xlist = xyvals[:,0].copy()
            ylist = xyvals[:,1].copy()
        del xyvals
    else:
        if isinstance(x,np.ndarray):
            xlist = x.tolist()
            ylist = y.tolist()
        elif not isinstance(x,list):
            xlist = [x]
            ylist = [y]
            single_coord = True
        else:
            xlist = x
            ylist = y

    # start by reading in WCS+distortion info for each image
    im1wcs = wcsutil.HSTWCS(inimage)
    if im1wcs.wcs.is_unity():
        print("####\nNo valid input WCS found in {}.\n  Results may be invalid.\n####\n".format(inimage))

    if util.is_blank(outimage):
        fname,fextn = fileutil.parseFilename(inimage)
        numsci = fileutil.countExtn(fname)
        chips = []
        for e in range(1,numsci+1):
            chips.append(wcsutil.HSTWCS(fname,ext=('sci',e)))
        if len(chips) == 0:
            chips = [im1wcs]
        im2wcs = distortion.utils.output_wcs(chips)
    else:
        im2wcs = wcsutil.HSTWCS(outimage)

    if im2wcs.wcs.is_unity():
        print("####\nNo valid output WCS found in {}.\n  Results may be invalid.\n####\n".format(outimage))

    # Setup the transformation
    p2p = wcs_functions.WCSMap(im1wcs,im2wcs)

    if direction[0].lower() == 'f':
        outx,outy = p2p.forward(xlist,ylist)
    else:
        outx,outy = p2p.backward(xlist,ylist)

    if isinstance(outx,np.ndarray):
        outx = outx.tolist()
        outy = outy.tolist()

    # add formatting based on precision here...
    xstr = []
    ystr = []
    fmt = "%."+repr(precision)+"f"
    for ox,oy in zip(outx,outy):
        xstr.append(fmt%ox)
        ystr.append(fmt%oy)

    if verbose or (not verbose and util.is_blank(output)):
        print('# Coordinate transformations for ',inimage)
        print('# X(in)      Y(in)             X(out)         Y(out)\n')
        for xs,ys,a,b in zip(xlist,ylist,xstr,ystr):
            print("%.4f  %.4f    %s  %s"%(xs,ys,a,b))

    # Create output file, if specified
    if output:
        f = open(output,mode='w')
        f.write("# Coordinates converted from %s\n"%inimage)
        for xs,ys in zip(xstr,ystr):
            f.write('%s    %s\n'%(xs,ys))
        f.close()
        print('Wrote out results to: ',output)

    if single_coord:
        outx = outx[0]
        outy = outy[0]
    return outx,outy


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

    outimage = util.check_blank(configObj['outimage'])
    tran(configObj['inimage'], outimage,direction=configObj['direction'],
            x = configObj['x'], y = configObj['y'], coords=coords,
            coordfile = coordfile, colnames = colnames,
            separator= sep, precision= configObj['precision'],
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
