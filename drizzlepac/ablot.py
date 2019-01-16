"""
In this step the median image now gets blotted back to create median-cleaned
images which can be compared directly with each input image to identify
cosmic-rays.

:Authors: Warren Hack

:License: :doc:`LICENSE`

"""
import os
import sys
import numpy as np
from stsci.tools import fileutil, teal, logutil
from . import outputimage
from . import wcs_functions
from . import processInput
from . import util
import stwcs
from stwcs import distortion

try:
    from . import cdriz
except ImportError:
    cdriz = None
    print('\n Coordinate transformation and image resampling library NOT found!')
    print('\n Please check the installation of this package to insure C code was built successfully.')
    raise ImportError

from .version import *

__all__ = ['blot', 'runBlot', 'help', 'getHelpAsString']

__taskname__ = 'drizzlepac.ablot'
_blot_step_num_ = 5

log = logutil.create_logger(__name__, level=logutil.logging.NOTSET)


#
#### User level interface run from TEAL
#

def blot(data, reference, outdata, configObj=None, wcsmap=wcs_functions.WCSMap,
         editpars=False, **input_dict):
    if input_dict is None:
        input_dict = {}
    input_dict['data'] = data
    input_dict['reference'] = reference
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
    _inimg = fileutil.openImage(_fname, memmap=False)
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
    del _inimg

    # read in WCS from source (drizzled) image
    source_wcs = stwcs.wcsutil.HSTWCS(configObj['data'])
    if source_wcs.wcs.is_unity():
        print("WARNING: No valid WCS found for input drizzled image: {}!".format(configObj['data']))

    # define blot_wcs
    blot_wcs = None
    _refname,_refextn = fileutil.parseFilename(configObj['reference'])
    if os.path.exists(_refname):
        # read in WCS from pre-existing output image
        blot_wcs = stwcs.wcsutil.HSTWCS(configObj['reference'])
        if blot_wcs.wcs.is_unity():
            print("WARNING: No valid WCS found for output image: {} !".format(configObj['reference']))

    # define blot WCS based on input images or specified reference WCS values
    if user_wcs_pars['user_wcs']:
        blot_wcs = wcs_functions.build_hstwcs(
            user_wcs_pars['raref'], user_wcs_pars['decref'],
            user_wcs_pars['xrefpix'], user_wcs_pars['yrefpix'],
            user_wcs_pars['outnx'], user_wcs_pars['outny'],
            user_wcs_pars['outscale'], user_wcs_pars['orient'] )
        configObj['coeffs'] = None

    # If blot_wcs is still not defined at this point, we have a problem...
    if blot_wcs is None:
        blot_wcs = stwcs.distortion.utils.output_wcs([source_wcs],undistort=False)

    out_wcs = blot_wcs.copy()
    # perform blotting operation now
    _outsci = do_blot(_insci, source_wcs, out_wcs, _expin, coeffs=configObj['coeffs'],
                    interp=configObj['interpol'], sinscl=configObj['sinscl'],
            stepsize=configObj['stepsize'], wcsmap=wcsmap)
    # create output with proper units and exptime-scaling
    if scale_pars['out_units'] == 'counts':
        if scale_pars['expout'] == 'input':
            _outscale = fileutil.getKeyword(configObj['reference'],scale_pars['expkey'])
            #_outscale = _expin
        else:
            _outscale = float(scale_pars['expout'])
        print("Output blotted images scaled by exptime of {}".format(_outscale))
        np.multiply(_outsci, _outscale, _outsci)

    # Add sky back in to the blotted image, as specified by the user
    if configObj['addsky']:
        skyval = _scihdu.header['MDRIZSKY']
    else:
        skyval = configObj['skyval']
    print("Added {} counts back in to blotted image as sky.".format(skyval))
    _outsci += skyval

    del _scihdu

    # Write output Numpy objects to a PyFITS file
    # Blotting only occurs from a drizzled SCI extension
    # to a blotted SCI extension...
    outputimage.writeSingleFITS(_outsci,blot_wcs, configObj['outdata'],configObj['reference'])


#
#### Top-level interface from inside AstroDrizzle
#
def runBlot(imageObjectList, output_wcs, configObj={},
            wcsmap=wcs_functions.WCSMap, procSteps=None):
    """
    runBlot(imageObjectList, output_wcs, configObj={},
            wcsmap=wcs_functions.WCSMap, procSteps=None)
    """
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
    """
    run_blot(imageObjectList, output_wcs, paramDict, wcsmap=wcs_functions.WCSMap)

    Perform the blot operation on the list of images.
    """
    # Insure that input imageObject is a list
    if not isinstance(imageObjectList, list):
        imageObjectList = [imageObjectList]
    #
    # Setup the versions info dictionary for output to PRIMARY header
    # The keys will be used as the name reported in the header, as-is
    #
    _versions = {'AstroDrizzle':__version__,
                 'PyFITS':util.__fits_version__,
                 'Numpy':util.__numpy_version__}

    _hdrlist = []

    for img in imageObjectList:

        for chip in img.returnAllChips(extname=img.scienceExt):

            print('    Blot: creating blotted image: ',chip.outputNames['data'])

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
                _inimg = fileutil.openImage(_fname, memmap=False)

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
            if skyval is not None:
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


def do_blot(source, source_wcs, blot_wcs, exptime, coeffs = True,
            interp='poly5', sinscl=1.0, stepsize=10, wcsmap=None):
    """ Core functionality of performing the 'blot' operation to create a single
        blotted image from a single source image.
        All distortion information is assumed to be included in the WCS specification
        of the 'output' blotted image given in 'blot_wcs'.

        This is the simplest interface that can be called for stand-alone
        use of the blotting function.

        Parameters
        ----------
        source
            Input numpy array of undistorted source image in units of 'cps'.
        source_wcs
            HSTWCS object representing source image distortion-corrected WCS.
        blot_wcs
            (py)wcs.WCS object representing the blotted image WCS.
        exptime
            exptime to use for scaling output blot image. A value of 1 will
            result in output blot image in units of 'cps'.
        coeffs
            Flag to specify whether or not to use distortion coefficients
            associated with blot_wcs. If False, do not apply any distortion
            model.
        interp
            Form of interpolation to use when blotting pixels. Valid options::

                "nearest","linear","poly3", "poly5"(default), "spline3", "sinc"
        sinscl
            Scale for sinc interpolation kernel (in output, blotted pixels)
        stepsize
            Number of pixels for WCS interpolation
        wcsmap
            Custom mapping class to use to provide transformation from
            drizzled to blotted WCS.  Default will be to use
            `drizzlepac.wcs_functions.WCSMap`.

    """
    _outsci = np.zeros(blot_wcs.array_shape, dtype=np.float32)

    # Now pass numpy objects to callable version of Blot...
    build=False
    misval = 0.0
    kscale = 1.0

    xmin = 1
    ymin = 1
    xmax, ymax = source_wcs.pixel_shape

    # compute the undistorted 'natural' plate scale for this chip
    if coeffs:
        wcslin = distortion.utils.make_orthogonal_cd(blot_wcs)
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
        print('Using default C-based coordinate transformation...')
        mapping = cdriz.DefaultWCSMapping(
            blot_wcs, source_wcs,
            blot_wcs.pixel_shape[0], blot_wcs.pixel_shape[1],
            stepsize
        )
        pix_ratio = source_wcs.pscale/wcslin.pscale
    else:
        #
        ##Using the Python class for the WCS-based transformation
        #
        # Use user provided mapping function
        print('Using coordinate transformation defined by user...')
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
    taskname = util.base_taskname(__taskname__, __package__)
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


blot.__doc__ = getHelpAsString(docstring = True, show_ver = False)
