""" pixtosky - A module to perform coordinate transformation from pixel coordinates
                in one image to pixel coordinates in another frame

    License:
        http://www.stsci.edu/resources/software_hardware/pyraf/LICENSE

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
    coords : str, optional
        full filename with path of file with starting x,y coordinates
    colnames : str, optional
        comma separated list of column names from 'coords' files
        containing x,y coordinates, respectively. Will default to
        first two columns if None are specified. Column names for ASCII
        files will use 'c1','c2',... convention.
    separator : str, optional
        non-blank separator used as the column delimiter in the coords file
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
                    "backward", coords='xy_sci1.dat', colnames=['c3','c4'],
                    output="xy_flt1.dat")

"""
from __future__ import division # confidence medium

import os,copy
import numpy as np

import pyfits
from stsci.tools import fileutil, teal
import wcs_functions
import util
from stwcs import wcsutil, distortion

# This is specifically NOT intended to match the package-wide version information.
__version__ = '0.1'
__vdate__ = '1-Mar-2011'

__taskname__ = 'pixtopix'

def tran(inimage,outimage,direction='forward',x=None,y=None,
        coords=None,colnames=None,separator=None,
        precision=6, output=None,verbose=True):
    """ Primary interface to perform coordinate transformations in pixel
        coordinates between 2 images using STWCS and full distortion models
        read from each image's header.
    """
    if coords is not None:
        if colnames in util.blank_list:
            colnames = ['c1','c2']
        # Determine columns which contain pixel positions
        cols = util.parse_colnames(colnames,coords)
        # read in columns from input coordinates file
        xyvals = np.loadtxt(coords,usecols=cols,delimiter=separator)
        xlist = xyvals[:,0].copy()
        ylist = xyvals[:,1].copy()
        del xyvals
    else:
        if not isinstance(x,list):
            xlist = [x]
            ylist = [y]
        else:
            xlist = x
            ylist = y

    # start by reading in WCS+distortion info for each image
    im1wcs = wcsutil.HSTWCS(inimage)

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

    # Setup the transformation
    p2p = wcs_functions.WCSMap(im1wcs,im2wcs)

    if direction[0].lower() == 'f':
        outx,outy = p2p.forward(xlist,ylist)
    else:
        outx,outy = p2p.backward(xlist,ylist)

    # add formatting based on precision here...
    xstr = []
    ystr = []
    fmt = "%."+repr(precision)+"f"
    for x,y in zip(outx,outy):
        xstr.append(fmt%x)
        ystr.append(fmt%y)

    if verbose or (not verbose and util.is_blank(output)):
        print '# Coordinate transformations for ',inimage
        print '# X(in)      Y(in)             X(out)         Y(out)\n'
        for x,y,a,b in zip(xlist,ylist,xstr,ystr):
            print "%.4f  %.4f    %s  %s"%(x,y,a,b)

    # Create output file, if specified
    if output:
        f = open(output,mode='w')
        f.write("# Coordinates converted from %s\n"%inimage)
        for x,y in zip(xstr,ystr):
            f.write('%s    %s\n'%(x,y))
        f.close()
        print 'Wrote out results to: ',output

    return outx,outy


#--------------------------
# TEAL Interface functions
#--------------------------
def run(configObj):

    coords = util.check_blank(configObj['coords'])
    colnames = util.check_blank(configObj['colnames'])
    sep = util.check_blank(configObj['separator'])
    outfile = util.check_blank(configObj['output'])

    outimage = util.check_blank(configObj['outimage'])
    tran(configObj['inimage'], outimage,direction=configObj['direction'],
            x = configObj['x'], y = configObj['y'],
            coords = coords, colnames = colnames,
            separator= sep, precision= configObj['precision'],
            output= outfile, verbose = configObj['verbose'])

def getHelpAsString():
    helpString = ''
    if teal:
        helpString += teal.getHelpFileAsString(__taskname__,__file__)

    if helpString.strip() == '':
        helpString += __doc__+'\n'

    return helpString
