""" skytopix - A module to perform coordinate transformation from sky to pixel coordinates.

    License:
        http://www.stsci.edu/resources/software_hardware/pyraf/LICENSE

    PARAMETERS
    ----------
    input : str
        full filename with path of input image, an extension name ['sci',1] should be
        provided if input is a multi-extension FITS file

    Optional Parameters
    -------------------
    ra : string, optional
        RA from input image
    dec : string, optional
        Dec from input image
    coordfile : str, optional
        full filename with path of file with sky coordinates
    colnames : str, optional
        comma separated list of column names from 'coordfile' files
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
    x : float
        X position of pixel. If more than 1 input value, then it will be a
        numpy array.
    y : float
        Y position of pixel. If more than 1 input value, then it will be a
        numpy array.

    NOTES
    -----
    This module performs a full distortion-corrected coordinate transformation
    based on all WCS keywords and any recognized distortion keywords from the
    input image header.

    Usage
    -----
    It can be called from within Python using the syntax::

        >>> from drizzlepac import skytopix
        >>> x,y = skytopix.rd2xy("input_flt.fits[sci,1]","00:22:36.79","-72:4:9.0")

    EXAMPLES
    --------

    1. The following command will transform the position 00:22:36.79 -72:4:9.0 into a
        position on the image 'input_flt.fits[sci,1]' using::

            >>> from drizzlepac import skytopix
            >>> x,y = skytopix.rd2xy("input_file_flt.fits[sci,1]", "00:22:36.79","-72:4:9.0")


    2. The set of sky positions from 'input_flt.fits[sci,1]' stored as
        the 3rd and 4th columns from the ASCII file 'radec_sci1.dat'
        will be transformed and written out to 'xy_sci1.dat' using::

            >>> from drizzlepac import skytopix
            >>> x,y = skytopix.rd2xy("input_flt.fits[sci,1]", coordfile='radec_sci1.dat',
                colnames=['c3','c4'], output="xy_sci1.dat")

"""
from __future__ import division # confidence medium

import os,copy
import numpy as np

import pyfits
from stsci.tools import fileutil, teal
import util,wcs_functions,tweakutils
import stwcs
from stwcs import distortion,wcsutil

# This is specifically NOT intended to match the package-wide version information.
__version__ = '0.1'
__vdate__ = '25-Feb-2011'

__taskname__ = 'skytopix'

blank_list = [None, '', ' ']

def rd2xy(input,ra=None,dec=None,coordfile=None,colnames=None,
            precision=6,output=None,verbose=True):
    """ Primary interface to perform coordinate transformations from
        pixel to sky coordinates using STWCS and full distortion models
        read from the input image header.
    """
    if coordfile is not None:
        if colnames in blank_list:
            colnames = ['c1','c2']
        else:
            colnames = colnames.split(',')
        # convert input file coordinates to lists of decimal degrees values
        xlist,ylist = tweakutils.readcols(coordfile,cols=colnames)
    else:
        # convert input value into decimal degrees value
        xval,yval = tweakutils.parse_skypos(ra,dec)
        xlist = [xval]
        ylist = [yval]

    # start by reading in WCS+distortion info for input image
    inwcs = wcsutil.HSTWCS(input)
    # Now, convert pixel coordinates into sky coordinates
    try:
        outx,outy = inwcs.all_sky2pix(xlist,ylist,1)
    except RuntimeError:
        outx,outy = inwcs.wcs_sky2pix(xlist,ylist,1)

    # add formatting based on precision here...
    xstr = []
    ystr = []
    fmt = "%."+repr(precision)+"f"
    for x,y in zip(outx,outy):
        xstr.append(fmt%x)
        ystr.append(fmt%y)

    if verbose or (not verbose and util.is_blank(output)):
        print '# Coordinate transformations for ',input
        print '# X      Y         RA             Dec\n'
        for x,y,r,d in zip(xstr,ystr,xlist,ylist):
            print "%s  %s    %s  %s"%(x,y,r,d)

    # Create output file, if specified
    if output:
        f = open(output,mode='w')
        f.write("# Coordinates converted from %s\n"%input)
        for x,y in zip(xstr,ystr):
            f.write('%s    %s\n'%(x,y))
        f.close()
        print 'Wrote out results to: ',output

    return outx,outy

    """ Convert colnames input into list of column numbers
    """
    cols = []
    if not isinstance(colnames,list):
        colnames = colnames.split(',')

    # parse column names from coords file and match to input values
    if coordfile is not None and fileutil.isFits(coordfile)[0]:
        # Open FITS file with table
        ftab = pyfits.open(coordfile)
        # determine which extension has the table
        for extn in ftab:
            if isinstance(extn,pyfits.BinTableHDU):
                # parse column names from table and match to inputs
                cnames = extn.columns.names
                if colnames is not None:
                    for c in colnames:
                        for name,i in zip(cnames,xrange(len(cnames))):
                            if c == name.lower(): cols.append(i)
                    if len(cols) < len(colnames):
                        errmsg = "Not all input columns found in table..."
                        ftab.close()
                        raise ValueError, errmsg
                else:
                    cols = cnames[:2]
                break
        ftab.close()
    else:
        for c in colnames:
            if isinstance(c, str):
                if c[0].lower() == 'c': cols.append(int(c[1:])-1)
                else:
                    cols.append(int(c))
            else:
                if isinstance(c, int):
                    cols.append(c)
                else:
                    errmsg = "Unsupported column names..."
                    raise ValueError, errmsg
    return cols


#--------------------------
# TEAL Interface functions
#--------------------------
def run(configObj):

    coordfile = util.check_blank(configObj['coordfile'])
    colnames = util.check_blank(configObj['colnames'])
    outfile = util.check_blank(configObj['output'])

    rd2xy(configObj['input'],
            ra = configObj['ra'], dec = configObj['dec'],
            coordfile = coordfile, colnames = colnames,
            precision= configObj['precision'],
            output= outfile, verbose = configObj['verbose'])

def getHelpAsString():
    helpString = ''
    if teal:
        helpString += teal.getHelpFileAsString(__taskname__,__file__)

    if helpString.strip() == '':
        helpString += __doc__+'\n'

    return helpString
