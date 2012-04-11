""" DrizzlePac - Python tasks to 'dither' combine images

This package provides the tools to align, manage and combine images using
algorithms originally developed as part of IRAF's 'dither' package. Tasks
supported by this package include:
  * astrodrizzle
  * tweakreg
  * resetbits
  * updatenpol
  * pixtosky, skytopix and pixtopix
  * tweakback

**Output**: The primary output from this task is the distortion-corrected,
cosmic-ray cleaned, and combined image as a FITS file.

"""

from __future__ import division  # confidence high

import os

from . import ablot
from . import adrizzle
from . import astrodrizzle
from . import buildmask
from . import createMedian
from . import drizCR
from . import imageObject
from . import mdzhandler
from . import outputimage
from . import processInput
from . import resetbits
from . import sky
from . import staticMask
from . import util
from . import wcs_functions


# These modules provide the user-interfaces to coordinate transformation tasks
from . import pixtosky
from . import skytopix
from . import pixtopix

# The following modules are for 'tweakreg' and are included here to make
# it easier to get to this code interactively
try:
    from . import tweakreg, catalogs, imgclasses, tweakutils, imagefindpars
except:
    print 'The libraries needed for "tweakreg" were not available!'
    print 'None of the code related to that task can be used at this time.'

# Add updatenpol to the list of tasks imported automatically here
from . import updatenpol
from . import buildwcs

# This module supports applying WCS from _drz to _flt files
from . import tweakback

# These lines allow TEAL to print out the names of TEAL-enabled tasks
# upon importing this package.
from stsci.tools import teal

teal.print_tasknames(__name__, os.path.dirname(__file__),
                     hidden=['adrizzle','ablot','buildwcs'])

# Begin Version Information -------------------------------------------
if False :
    __version__ = ''
    __svn_version__ = 'Unable to determine SVN revision'
    __full_svn_info__ = ''
    __setup_datetime__ = None

    try:
        __version__ = __import__('pkg_resources').\
                            get_distribution('drizzlepac').version
    except:
        pass

else :
    __version__ = '4.3.0dev'


__vdate__ = '01-Mar-2012'


# Revision based version info
# End Version Information ---------------------------------------------
try:
    from drizzlepac.svninfo import (__svn_version__, __full_svn_info__,
                                     __setup_datetime__)
except ImportError:
    pass

def help():
    msg = \
""" The DrizzlePac package contains a suite of tasks that allow users to align HST images, combine them, and perform coordinate transformations on source positions.

drizzlepac:
       astrodrizzle - primary task for combining images, removing cosmic rays, and removing distortion
           tweakreg - task to compute offsets in WCS between images and a reference image or reference frame
      imagefindpars - sub-task containing parameters to find point sources used by tweakreg to build source catalogs for each tweakreg input image
          tweakback - apply an updated WCS solution created by tweakreg for a drizzled image to the constituent distorted (flt.fits) images
           pixtopix - task to convert pixel positions from an input image to pixel positions in an output WCS or image
           pixtosky - task to convert pixel positions from an input image to sky coordinates with full distortion correction as appropriate
           skytopix - task to convert sky positions to pixel positions in an image
          resetbits - sub-task to reset specified flt.fits data quality (DQ) values to 0
         updatenpol - task to add the names of the new ACS distortion reference files NPOLFILE and D2IMFILE then update headers to include residual distortion corrections as image extensions

fitsblender:
       blendheaders - task to merge the keywords from all images used to create a drizzle product into a single header with a table extension using rules defined for each instrument

stwcs:
    apply_headerlet - apply a headerlet to a file
  archive_headerlet - save a WCS solution as a headerlet extension and write it out as a headerlet FITS file
   attach_headerlet - attach a headerlet as an extension to a file
   delete_headerlet - delete a headerlet extension from a file
  extract_headerlet - write out a headerlet extension as a separate FITS file
  headerlet_summary - print a summary of all headerlet extensions in a file
  restore_headerlet - replace current WCS solution with the WCS solution from a headerlet extension
    write_headerlet - save a WCS solution as a separate headerlet FITS file
          updatewcs - recompute the WCS keywords and import the distortion model from the reference files
"""
    print msg
