"""
Functions to build mask files for PyDrizzle.
    - buildMaskImage(rootname,bitvalue,extname='DQ',extver=None):
        This function takes the DQ array(or array from extension given)
        and builds a bit mask as an input weighting mask for 'drizzle'.

    - buildShadowMaskImage(rootname,detnum,replace=no):
        This function builds a weighting image for use with drizzle from
        the WFPC2 shadow mask functions derived from 'wmosaic'.

:License: :doc:`/LICENSE`

"""

#
# Revision History:
#
#   14-Nov-2001 - Revised to produce proper simple FITS image starting
#                   with 'SIMPLE=T' keyword instead of 'XTENSION'. WJH
#
#   25-Jan-2002 - Converted to work with PyFITS V0.6.1 and numarray. WJH
#
#   12-June-2002 - Added check to remove old versions of or existing mask
#                  files before building a new one.
#   23-July-2002 - Saved space in mask file by using only Int8 instead of Int16
#
#   7-Mar-2003   - Updated to use memory mapping in PyFITS where possible
#
#   4-Jun-2003   - WFPC2 mask file now incorporates .c1h file if present
#
#  13-Jan-2005   - Simplified filename generation into one function and
#                   added output filename as input parameter to 'buildMask'
#                   functions. WJH
#
import os

from stsci.tools import fileutil
from stsci.tools.bitmask import bitfield_to_boolean_mask

from astropy.io import fits
import numpy as np

from . import processInput, util

# Define global variable for number of buffer pixels (rows/columns) used to define shadow mask overlap
SHADOW_BUFFER = 5  # pixels

__taskname__ = 'buildmask'
#
# Interactive interface
#


def run(configObj=None, input_dict={}, loadOnly=False):
    """ Build DQ masks from all input images, then apply static mask(s).
    """
    # gets configObj defaults using EPAR/TEAL.
    #
    # Also insure that the input_dict (user-specified values) are folded in
    # with a fully populated configObj instance.
    configObj = util.getDefaultConfigObj(__taskname__,configObj,input_dict,loadOnly=loadOnly)
    if configObj is None:
        return
    # Define list of imageObject instances and output WCSObject instance
    # based on input paramters
    imgObjList,outwcs = processInput.setCommonInput(configObj)

    # Build DQ masks for all input images.
    buildMask(imgObjList,configObj)


def buildDQMasks(imageObjectList,configObj):
    """ Build DQ masks for all input images.
    """
    # Insure that input imageObject is a list
    if not isinstance(imageObjectList, list):
        imageObjectList = [imageObjectList]

    for img in imageObjectList:
        img.buildMask(configObj['single'], configObj['bits'])


def buildMask(dqarr, bitvalue):
    """ Builds a bit-mask from an input DQ array and a bitvalue flag """
    return bitfield_to_boolean_mask(dqarr.astype(np.uint16), ignore_flags=bitvalue, good_mask_value=1,
                                    dtype=np.uint8)


def buildMaskImage(rootname, bitvalue, output, extname='DQ', extver=1):
    """ Builds mask image from rootname's DQ array
        If there is no valid 'DQ' array in image, then return
        an empty string.
    """

    # If no bitvalue is set or rootname given, assume no mask is desired
    # However, this name would be useful as the output mask from
    # other processing, such as MultiDrizzle, so return it anyway.
    #if bitvalue == None or rootname == None:
    #    return None

    # build output name
    maskname = output

    # If an old version of the maskfile was present, remove it and rebuild it.
    if fileutil.findFile(maskname):
        fileutil.removeFile(maskname)

    # Open input file with DQ array
    fdq = fileutil.openImage(rootname, mode='readonly', memmap=False)
    try:
        _extn = fileutil.findExtname(fdq, extname, extver=extver)
        if _extn is not None:
            # Read in DQ array
            dqarr = fdq[_extn].data
        else:
            dqarr = None

        # For the case where there is no DQ array,
        # create a mask image of all ones.
        if dqarr is None:
            # We need to get the dimensions of the output DQ array
            # Since the DQ array is non-existent, look for the SCI extension
            _sci_extn = fileutil.findExtname(fdq,'SCI',extver=extver)
            if _sci_extn is not None:
                _shape = fdq[_sci_extn].data.shape
                dqarr = np.zeros(_shape,dtype=np.uint16)
            else:
                raise Exception
        # Build mask array from DQ array
        maskarr = buildMask(dqarr,bitvalue)
        #Write out the mask file as simple FITS file
        fmask = fits.open(maskname, mode='append', memmap=False)
        maskhdu = fits.PrimaryHDU(data = maskarr)
        fmask.append(maskhdu)

        #Close files
        fmask.close()
        del fmask
        fdq.close()
        del fdq

    except:
        fdq.close()
        del fdq
        # Safeguard against leaving behind an incomplete file
        if fileutil.findFile(maskname):
            os.remove(maskname)
        _errstr = "\nWarning: Problem creating MASK file for "+rootname+".\n"
        #raise IOError, _errstr
        print(_errstr)
        return None

    # Return the name of the mask image written out
    return maskname


# Utility functions for generating X and Y arrays for
# creating the WFPC2 shadow masks.
# Coefficients are from WMOSAIC function 'get_edge2.x'
"""
_coeffs = {'1':[[52.20921,0.009072887,0.009072887],[42.62779,0.009122855,-1.106709E-5]],
          '2':[[21.77283,0.01842164,-1.398300E-5],[47.68184,0.00265608,-1.468158E-5]],
          '3':[[44.18944,0.0138938,-1.412296E-5],[30.23626,0.008156041,-1.491324E-5]],
          '4':[[44.56632,0.003509023,-1.278723E-5],[40.91462,0.01273679,-1.063462E-5]]}
"""


def _func_Shadow_WF1y(x,y):
    return y + 0.5 - (42.62779 + 0.009122855*y - 1.106709E-5*y**2)


def _func_Shadow_WF1x(x,y):
    return x + 0.5 - (52.20921 + 0.009072887*x - 9.941337e-6*x**2)


def _func_Shadow_WF2y(x,y):
    return y + 0.5 - (47.68184 + 0.00265608*y - 1.468158E-5*y**2)


def _func_Shadow_WF2x(x,y):
    return x + 0.5 - (21.77283 + 0.01842164*x - 1.398300E-5*x**2)


def _func_Shadow_WF3y(x,y):
    return y + 0.5 - (30.23626 + 0.008156041*y - 1.491324E-5*y**2 )


def _func_Shadow_WF3x(x,y):
    return x + 0.5 - (44.18944 + 0.0138938*x - 1.412296E-5*x**2)


def _func_Shadow_WF4y(x,y):
    return y + 0.5 - (40.91462 + 0.01273679*y - 1.063462E-5*y**2)


def _func_Shadow_WF4x(x,y):
    return x + 0.5 - (44.56632 + 0.003509023*x - 1.278723E-5*x**2)


# Function for creating the weighting image for WFPC2 drizzling
# from the 'wmosaic' shadow mask polynomials.


def buildShadowMaskImage(dqfile, detnum, extnum, maskname, bitvalue=None, binned=1):
    """ Builds mask image from WFPC2 shadow calibrations.
        detnum - string value for 'DETECTOR' detector

        This function uses the region defined by the calibration functions above to
        identify the functional shadow mask pixels from the mask created using the DQ array with all values
        except DQ=2.  This shadow mask then applied to the DQ-based mask and then expanded by a few pixels,
        as defined by ``SHADOW_BUFFER``, before being used to expand shadow mask region in the user-defined
        DQ-based mask.  This logic insures that the pixels defined as GOOD by the user in the imaging
        portion of the chip outside the shadowed region gets used while reducing the number of pixels
        in the shaowed region to increase the amount of unmasked pixels overlapping from chip to chip.

        Parameters
        ----------
        dqfile : str
            Filename of file containing the exposures DQ array

        detnum : int
            Detector number from exposure header

        extnum : int
            Extension number for SCI and DQ arrays

        maskname : str
            Filename to use for writing out computed mask

        bitvalue : int, optional
            Bits user defines as GOOD from DQ array for defining final output mask.  If None,
            DQ array is not used, and only functional shadow mask gets used.

        binned : int, optional
            Specifies whether or not exposure was taken in binned mode.

        Returns
        --------
        maskname : str
            Filename of mask file written out by function.  If there was an error in
            accessing the exposures DQ array, this will be None.

    """
    # insure detnum is a string
    if not isinstance(detnum, str):
        detnum = repr(detnum)

    _funcroot = '_func_Shadow_WF'

    # build template shadow mask's filename

    # If an old version of the maskfile was present, remove it and rebuild it.
    if fileutil.findFile(maskname):
        fileutil.removeFile(maskname)

    _flt_file = dqfile.endswith('flt.fits')  # check to see if DQ file is just an extension from an FLT file
    # Check for existance of input .c1h file for use in making inmask file
    # Check to see if file exists...
    # If not, create the file.
    # This takes a long time to run, so it should be done
    # only when absolutely necessary...
    _funcx = _funcroot+detnum+'x'
    _funcy = _funcroot+detnum+'y'

    _xarr = np.clip(np.fromfunction(eval(_funcx), (800, 800)), 0.0, 1.0).astype(np.uint8)
    _yarr = np.clip(np.fromfunction(eval(_funcy), (800, 800)), 0.0, 1.0).astype(np.uint8)
    maskarr = _xarr * _yarr

    if binned !=1:
        bmaskarr = maskarr[::2, ::2]
        bmaskarr *= maskarr[1::2, ::2]
        bmaskarr *= maskarr[::2, 1::2]
        bmaskarr *= maskarr[1::2, 1::2]
        maskarr = bmaskarr.copy()
        del bmaskarr

    if bitvalue is not None:
        #
        # Build full mask based on .c1h combined with shadow mask
        #
        fdq = fileutil.openImage(dqfile, mode='readonly', memmap=False)
        if _flt_file:
            _extnum = ('DQ', int(extnum))
        else:
            _extnum = int(extnum)
        try:
            # Read in DQ array from .c1h and from shadow mask files
            dqarr = fdq[_extnum].data

            # Build mask array from DQ array
            dqmaskarr = buildMask(dqarr, bitvalue)

        except Exception:
            # Safeguard against leaving behind an incomplete file
            if fileutil.findFile(maskname):
                os.remove(maskname)
            _errstr = "\nWarning: Problem creating DQMASK file for " + dqfile + ".\n"
            print(_errstr)
            return None
        finally:
            fdq.close()
            del fdq
    else:
        # simply use the functional shadow mask
        dqmaskarr = maskarr

    # Write out the mask file as simple FITS file
    fdqmask = fits.open(maskname, mode='append', memmap=False)
    maskhdu = fits.PrimaryHDU(data=dqmaskarr)
    fdqmask.append(maskhdu)

    # Close files
    fdqmask.close()
    del fdqmask

    # Return the name of the mask image written out
    return maskname
