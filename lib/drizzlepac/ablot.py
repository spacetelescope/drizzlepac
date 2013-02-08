from __future__ import division  # confidence medium

import os
import sys
import util
import numpy as np
from stsci.tools import fileutil, teal, logutil
import outputimage, wcs_functions, processInput, util
import stwcs
from stwcs import distortion

from .version import *

try:
    import cdriz
except ImportError:
    cdriz = None
    print '\n Coordinate transformation and image resampling library NOT found!'
    print '\n Please check the installation of this package to insure C code was built successfully.'
    raise ImportError


__taskname__ = 'drizzlepac.ablot'
_blot_step_num_ = 5

log = logutil.create_logger(__name__)


#
#### User level interface run from TEAL
#

def blot(data, outdata, configObj=None, wcsmap=wcs_functions.WCSMap,
         editpars=False, **input_dict):
    if input_dict is None:
        input_dict = {}
    input_dict['data'] = data
    input_dict['outdata'] = outdata

    # If called from interactive user-interface, configObj will not be
    # defined yet, so get defaults using EPAR/TEAL.
    #
    # Also insure that the input_dict (user-specified values) are folded in
    # with a fully populated configObj instance.
    configObj = util.getDefaultConfigObj(__taskname__, configObj,
                                         input_dict, loadOnly=(not editpars))
    if configObj is None:
        return

    if not editpars:
        run(configObj, wcsmap=wcsmap)


def run(configObj,wcsmap=None):
    """
    Run the blot task based on parameters provided interactively by the user.
    """

    # Insure all output filenames specified have .fits extensions
    if configObj['outdata'][-5:] != '.fits': configObj['outdata'] += '.fits'

    scale_pars = configObj['Data Scaling Parameters']
    user_wcs_pars = configObj['User WCS Parameters']

    # PyFITS can be used here as it will always operate on
    # output from PyDrizzle (which will always be a FITS file)
    # Open the input (drizzled?) image
    _fname,_sciextn = fileutil.parseFilename(configObj['data'])
    _inimg = fileutil.openImage(_fname)
    _expin = fileutil.getKeyword(configObj['data'],scale_pars['expkey'],handle=_inimg)

    # Return the PyFITS HDU corresponding to the named extension
    _scihdu = fileutil.getExtn(_inimg,_sciextn)
    _insci = _scihdu.data.copy()

    _inexptime = 1.0
    if scale_pars['in_units'] == 'counts':
        if scale_pars['expkey'] in _inimg['PRIMARY'].header:
            _inexptime = _inimg['PRIMARY'].header[scale_pars['expkey']]
        elif 'DRIZEXPT' in _inimg['PRIMARY'].header:
            # Try keyword written out by new 'drizzle' if no valid 'expkey' was given
            _inexptime = _inimg['PRIMARY'].header['DRIZEXPT']
        else:
            raise ValueError('No valid exposure time keyword could be found '
                             'for input %s' % configObj['data'])
    # always convert input to 'cps' for blot() algorithm
    if _inexptime != 0.0 or _inexptime != 1.0:
        np.divide(_insci, _inexptime, _insci)

    _inimg.close()
    del _inimg, _scihdu

    # read in WCS from source (drizzled) image
    source_wcs = stwcs.wcsutil.HSTWCS(configObj['data'])

    # define blot_wcs
    blot_wcs = None
    if os.path.exists(configObj['outdata']):
        # read in WCS from pre-existing output image
        blot_wcs = stwcs.wcsutil.HSTWCS(configObj['outdata'])
    else:
        # define blot WCS based on input images or specified reference WCS values
        if not util.is_blank(user_wcs_pars['refimage']):
            blot_wcs = stwcs.wcsutil.HSTWCS(user_wcs_pars['refimage'])
        else:
            if not util.is_blank(user_wcs_pars['outscale']):
                blot_wcs = wcs_functions.build_hstwcs(
                    user_wcs_pars['raref'], user_wcs_pars['decref'],
                    user_wcs_pars['xrefpix'], user_wcs_pars['yrefpix'],
                    user_wcs_pars['outnx'], user_wcs_pars['outny'],
                    user_wcs_pars['outscale'], user_wcs_pars['orient'] )
                configObj['coeffs'] = None

        # If blot_wcs is still not defined at this point, we have a problem...
        if blot_wcs is None:
            blot_wcs = stwcs.distortion.utils.output_wcs([source_wcs],undistort=False)

    # perform blotting operation now
    _outsci = do_blot(_insci, source_wcs, blot_wcs, _expin, coeffs=configObj['coeffs'],
                    interp=configObj['interpol'], sinscl=configObj['sinscl'],
            stepsize=configObj['stepsize'], wcsmap=wcsmap)

    # create output with proper units and exptime-scaling
    if scale_pars['out_units'] == 'counts':
        if scale_pars['expout'] == 'input':
            _outscale = _expin
        else:
            _outscale = float(scale_pars['expout'])
        np.multiply(_outsci, _outscale, _outsci)

    # Add sky back in to the blotted image, as specified by the user
    if configObj['addsky']:
        skyval = _scihdu.header['MDRIZSKY']
    else:
        skyval = configObj['skyval']
    _outsci += skyval

    # Write output Numpy objects to a PyFITS file
    # Blotting only occurs from a drizzled SCI extension
    # to a blotted SCI extension...
    outputimage.writeSingleFITS(_outsci,blot_wcs, configObj['outdata'],configObj['data'],blot=True)


def help():
    print getHelpAsString()


def getHelpAsString():
    """
    Return useful help from a file in the script directory called module.help
    """
    helpString = 'ABLOT Version '+__version__+' Revision date: '+__vdate__
    try:
        helpString += teal.getHelpFileAsString(__taskname__,__file__)
    except IndexError:
        pass
    return helpString


#
#### Top-level interface from inside MultiDrizzle
#
def runBlot(imageObjectList, output_wcs, configObj={},
            wcsmap=wcs_functions.WCSMap, procSteps=None):
    if procSteps is not None:
        procSteps.addStep('Blot')

    blot_name = util.getSectionName(configObj, _blot_step_num_)

    # This can be called directly from MultiDrizle, so only execute if
    # switch has been turned on (no guarantee MD will check before calling).
    if configObj[blot_name]['blot']:
        paramDict = buildBlotParamDict(configObj)

        log.info('USER INPUT PARAMETERS for Blot Step:')
        util.printParams(paramDict, log=log)

        run_blot(imageObjectList, output_wcs.single_wcs, paramDict,
                 wcsmap=wcsmap)
    else:
        log.info('Blot step not performed.')

    if procSteps is not None:
        procSteps.endStep('Blot')


# Run 'drizzle' here...
#
def buildBlotParamDict(configObj):
    blot_name = util.getSectionName(configObj,_blot_step_num_)

    paramDict = {'blot_interp':configObj[blot_name]['blot_interp'],
                'blot_sinscl':configObj[blot_name]['blot_sinscl'],
                'blot_addsky':configObj[blot_name]['blot_addsky'],
                'blot_skyval':configObj[blot_name]['blot_skyval'],
                'coeffs':configObj['coeffs']}
    return paramDict

def _setDefaults(configObj={}):
    """ set up the default parameters to run drizzle
        build,single,units,wt_scl,pixfrac,kernel,fillval,
        rot,scale,xsh,ysh,blotnx,blotny,outnx,outny,data
    """

    paramDict={"build":True,
              "single":True,
              "in_units":"cps",
              "wt_scl":1.,
              "pixfrac":1.,
              "kernel":"square",
              "fillval":999.,
              "rot":0.,
              "scale":1.,
              "xsh":0.,
              "ysh":0.,
              "blotnx":2048,
              "blotny":2048,
              "outnx":4096,
              "outny":4096,
              "data":None,
              "driz_separate":True,
              "driz_combine":False}

    if(len(configObj) !=0):
        for key in configObj.keys():
            paramDict[key]=configObj[key]

    return paramDict

def run_blot(imageObjectList,output_wcs,paramDict,wcsmap=wcs_functions.WCSMap):
    """ Perform the blot operation on the list of images.
    """
    # Insure that input imageObject is a list
    if not isinstance(imageObjectList, list):
        imageObjectList = [imageObjectList]
    #
    # Setup the versions info dictionary for output to PRIMARY header
    # The keys will be used as the name reported in the header, as-is
    #
    _versions = {'AstroDrizzle':__version__,
                 'PyFITS':util.__pyfits_version__,
                 'Numpy':util.__numpy_version__}

    _hdrlist = []

    for img in imageObjectList:

        for chip in img.returnAllChips(extname=img.scienceExt):

            print '    Blot: creating blotted image: ',chip.outputNames['data']

            #### Check to see what names need to be included here for use in _hdrlist
            chip.outputNames['driz_version'] = _versions['AstroDrizzle']
            outputvals = chip.outputNames.copy()
            outputvals.update(img.outputValues)
            outputvals['blotnx'] = chip.wcs.naxis1
            outputvals['blotny'] = chip.wcs.naxis2
            _hdrlist.append(outputvals)

            plist = outputvals.copy()
            plist.update(paramDict)

            # PyFITS can be used here as it will always operate on
            # output from PyDrizzle (which will always be a FITS file)
            # Open the input science file
            medianPar = 'outMedian'
            outMedianObj = img.getOutputName(medianPar)
            if img.inmemory:
                outMedian = img.outputNames[medianPar]
                _fname,_sciextn = fileutil.parseFilename(outMedian)
                _inimg = outMedianObj
            else:
                outMedian = outMedianObj
                _fname,_sciextn = fileutil.parseFilename(outMedian)
                _inimg = fileutil.openImage(_fname)

            # Return the PyFITS HDU corresponding to the named extension
            _scihdu = fileutil.getExtn(_inimg,_sciextn)
            _insci = _scihdu.data.copy()
            _inimg.close()
            del _inimg, _scihdu

            _outsci = do_blot(_insci, output_wcs,
                   chip.wcs, chip._exptime, coeffs=paramDict['coeffs'],
                   interp=paramDict['blot_interp'], sinscl=paramDict['blot_sinscl'],
                   wcsmap=wcsmap)
            # Apply sky subtraction and unit conversion to blotted array to
            # match un-modified input array
            if paramDict['blot_addsky']:
                skyval = chip.computedSky
            else:
                skyval = paramDict['blot_skyval']
            _outsci /= chip._conversionFactor
            _outsci += skyval
            log.info('Applying sky value of %0.6f to blotted image %s'%
                        (skyval,chip.outputNames['data']))

            # Write output Numpy objects to a PyFITS file
            # Blotting only occurs from a drizzled SCI extension
            # to a blotted SCI extension...

            _outimg = outputimage.OutputImage(_hdrlist, paramDict, build=False, wcs=chip.wcs, blot=True)
            _outimg.outweight = None
            _outimg.outcontext = None
            outimgs = _outimg.writeFITS(plist['data'],_outsci,None,
                                versions=_versions,blend=False,
                                virtual=img.inmemory)

            img.saveVirtualOutputs(outimgs)
            #_buildOutputFits(_outsci,None,plist['outblot'])
            _hdrlist = []

            del _outsci

        del _outimg

def do_blot(source, source_wcs, blot_wcs, exptime, coeffs = True, interp='poly5', sinscl=1.0,
            stepsize=10, wcsmap=None):
    """ Core functionality of performing the 'blot' operation to create a single
        blotted image from a single source image.
        All distortion information is assumed to be included in the WCS specification
        of the 'output' blotted image given in 'blot_wcs'.

        Parameters
        ----------
        source
            Input numpy array of undistorted source image in units of 'cps'.
        source_wcs
            HSTWCS object representing source image WCS.
        blot_wcs
            HSTWCS object representing the blotted image WCS.
        exptime

    """
    _outsci = np.zeros((blot_wcs.naxis2,blot_wcs.naxis1),dtype=np.float32)

    # Now pass numpy objects to callable version of Blot...
    build=False
    misval = 0.0
    kscale = 1.0

    xmin = 1
    xmax = source_wcs.naxis1
    ymin = 1
    ymax = source_wcs.naxis2

    # compute the undistorted 'natural' plate scale for this chip
    if coeffs:
        wcslin = distortion.utils.undistortWCS(blot_wcs)
    else:
        wcslin = blot_wcs
        blot_wcs.sip = None
        blot_wcs.cpdis1 = None
        blot_wcs.cpdis2 = None
        blot_wcs.det2im = None

    if wcsmap is None and cdriz is not None:
        """
        Use default C mapping function.
        """
        print 'Using default C-based coordinate transformation...'
        mapping = cdriz.DefaultWCSMapping(blot_wcs,source_wcs,int(blot_wcs.naxis1),int(blot_wcs.naxis2),stepsize)
        pix_ratio = source_wcs.pscale/wcslin.pscale
    else:
        #
        ##Using the Python class for the WCS-based transformation
        #
        # Use user provided mapping function
        print 'Using coordinate transformation defined by user...'
        if wcsmap is None:
            wcsmap = wcs_functions.WCSMap
        wmap = wcsmap(blot_wcs,source_wcs)
        mapping = wmap.forward
        pix_ratio = source_wcs.pscale/wcslin.pscale

    t = cdriz.tblot(
        source, _outsci,xmin,xmax,ymin,ymax,
        pix_ratio, kscale, 1.0, 1.0,
        'center',interp, exptime,
        misval, sinscl, 1, mapping)
    del mapping

    return _outsci
