""" Pixreplace -- Replace pixels which have one value with another value

    :License: :doc:`LICENSE`


    PARAMETERS
    -----------
    input : str, @-file, list of filenames
        Filename(s) of image to be processed.

    pixvalue : float
        Pixel value from `input` file to be replaced.
        [Default: np.nan]

    newvalue : float
        New pixel value to use to replace `pixvalue`.
        [Default: 0.0]

    ext : int, list of ints, None
        Extensions from `input` file to process with new pixel values.
        If None (default), all image extensions (and only image extensions)
        will be processed.

    Usage
    -----
    It can be called from within Python using the syntax::

        >>> from drizzlepac import pixreplace
        >>> pixreplace.replace('adriz_nanSCI_drz.fits')
        or
        >>> epar pixreplace

    EXAMPLES
    ---------
    1.  Replace all pixels in all extensions which have a value of NaN
        in 'adriz_nanSCI_drz.fits' with a constant value of 0.0.
    >>> from drizzlepac import pixreplace
    >>> pixreplace.replace('adriz_nanSCI_drz.fits')

    2.  Replace pixels in first extension only which have a value of NaN
        in both 'j8c061vnq_drc.fits' and 'j8c061nyq_drc.fits'
        with a constant value of 0.0.
    >>> pixreplace.replace('j8c061vnq_drc.fits,j8c061nyq_drc.fits', ext=1)

    3.  Replace pixels in first and fourth extensions which have a value of 0.0
        in 'adriz_nanSCI_drz.fits' with a value of -999.
    >>> pixreplace.replace('adriz_nanSCI_drz.fits',pixvalue=0.0, newvalue=-999, ext=[1,4])


"""

import os

import numpy as np
from astropy.io import fits

from stsci.tools import parseinput

from . import util

# This is specifically NOT intended to match the package-wide version information.
__version__ = '0.1'
__version_date__ = '27-May-2015'

__taskname__ = 'pixreplace'

def replace(input, **pars):
    """ Replace pixels in `input` that have a value of `pixvalue`
        with a value given by `newvalue`.


    """
    pixvalue = pars.get('pixvalue', np.nan)
    if pixvalue is None: pixvalue = np.nan # insure that None == np.nan

    newvalue = pars.get('newvalue', 0.0)
    ext = pars.get('ext',None)
    if ext in ['',' ','None',None]:
        ext = None

    files = parseinput.parseinput(input)[0]

    for f in files:
        fimg = fits.open(f, mode='update', memmap=False)

        if ext is None:
            # replace pixels in ALL extensions
            extn = [i for i in fimg]
        else:
            if type(ext) == type([]):
                extn = [fimg[e] for e in ext]
            else:
                extn = [fimg[ext]]

        for e in extn:
            if e.data is not None and e.is_image: # ignore empty Primary HDUs
                print("Converting {}[{},{}] value of {} to {}".format(
                        f,e.name,e.ver,pixvalue,newvalue))
                if np.isnan(pixvalue):
                    e.data[np.isnan(e.data)] = newvalue
                else:
                    e.data[np.where(e.data == pixvalue)] = newvalue

        fimg.close()


##############################
#
# TEAL Interfaces
#
##############################
def run(configobj):

    ext = util.check_blank(configobj['ext'])
    if ext is not None:
        ext = ext.split(',')
        ext = [int(e) for e in ext]

    replace(configobj['input'],
            pixvalue=configobj['pixvalue'],
            newvalue=configobj['newvalue'],
            ext=ext)

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
