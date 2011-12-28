""" AstroDither - Python tasks to 'dither' combine images

This package provides the tools to align, manage and combine images using 
algorithms originally developed as part of IRAF's 'dither' package. Tasks
supported by this package include: 
  * astrodrizzle
  * tweakreg
  * resetbits
  * updatenpol
  * pixtosky, skytopix and pixtopix 

**Output**: The primary output from this task is the distortion-corrected,
cosmic-ray cleaned, and combined image as a FITS file.

"""
from __future__ import division # confidence high

import os,string,glob

import astrodrizzle
import buildmask
import adrizzle
import ablot
import imageObject
import outputimage
import processInput,mdzhandler
import sky
import createMedian
import drizCR
import staticMask
import util
import wcs_functions
import resetbits

# These modules provide the user-interfaces to coordinate transformation tasks
import pixtosky
import skytopix
import pixtopix

# The following modules are for 'tweakreg' and are included here to make
# it easier to get to this code interactively
try:
    import tweakreg, catalogs, imgclasses, tweakutils
    import imagefindpars
except:
    print 'The libraries needed for "tweakreg" were not available!'
    print 'None of the code related to that task can be used at this time.'

# Add updatenpol to the list of tasks imported automatically here
import updatenpol
import buildwcs

# These lines allow TEAL to print out the names of TEAL-enabled tasks
# upon importing this package.
from stsci.tools import teal
teal.print_tasknames(__name__, os.path.dirname(__file__))

# Begin Version Information -------------------------------------------
if False :
    __version__ = ''
    __svn_version__ = 'Unable to determine SVN revision'
    __full_svn_info__ = ''
    __setup_datetime__ = None

    try:
        __version__ = __import__('pkg_resources').\
                            get_distribution('astrodither').version
    except:
        pass

else :
    __version__ = '4.2.2dev'

__vdate__ = '25-Oct-2011'
# Revision based version info
# End Version Information ---------------------------------------------
try:
    from astrodither.svninfo import (__svn_version__, __full_svn_info__,
                                     __setup_datetime__)
except ImportError:
    pass
