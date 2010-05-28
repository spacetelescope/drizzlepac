from __future__ import division # confidence medium

import sys,types,os
import util
from util import _ptime
import numpy as np
from pytools import fileutil
import outputimage,wcs_functions,processInput,util
from stwcs import distortion

try:
    import cdriz as arrdriz
except ImportError:
    arrdriz = None
    
__taskname__ = 'betadrizzle.blot'
_blot_step_num_ = 5

#
#### User level interface run from TEAL
#
def run(configObj,wcsmap=wcs_functions.WCSMap):

    # Define list of imageObject instances and output WCSObject instance
    # based on input paramters

    # Explicitly turn off making copies so as to not over-write any analysis
    # already performed on the data.
    configObj['workinplace'] = False

    imgObjList,outwcs = processInput.setCommonInput(configObj)

    runblot(imgObjList,outwcs,configObj,wcsmap=wcsmap)

def help():
    print getHelpAsString()
    
def getHelpAsString():
    """ 
    return useful help from a file in the script directory called module.help
    """
    #get the local library directory where the code is stored
    localDir=os.path.split(__file__)
    helpfile=__taskname__.split(".")
    helpfile=localDir[0]+"/"+helpfile[1]+".help"
    
    if os.access(helpfile,os.R_OK):
        fh=open(helpfile,'r')
        ss=fh.readlines()
        fh.close()
        helpString=""
        for line in ss:
            helpString+=line
    else:    
        helpString=__doc__

    return helpString

# 
#### Interactive interface for running drizzle tasks separately
#
def blot(configObj=None,wcsmap=wcs_functions.WCSMap,editpars=False,**input_dict):

    # If called from interactive user-interface, configObj will not be 
    # defined yet, so get defaults using EPAR/TEAL.
    #
    # Also insure that the input_dict (user-specified values) are folded in
    # with a fully populated configObj instance.
    configObj = util.getDefaultConfigObj(__taskname__,configObj,input_dict,loadOnly=(not editpars))
    if configObj is None:
        return
    
    if not editpars:
        run(configObj,wcsmap=wcsmap)

#
#### Top-level interface from inside MultiDrizzle
#
def runBlot(imageObjectList, output_wcs, configObj={},wcsmap=wcs_functions.WCSMap,procSteps=None):
    if procSteps is not None:
        procSteps.addStep('Blot')

    blot_name = util.getSectionName(configObj,_blot_step_num_)

    # This can be called directly from MultiDrizle, so only execute if
    # switch has been turned on (no guarantee MD will check before calling).
    if configObj[blot_name]['blot']:
        paramDict = buildBlotParamDict(configObj)
        run_blot(imageObjectList, output_wcs.single_wcs, paramDict, wcsmap=wcsmap)
    else:
        print 'Blot step not performed.'

    if procSteps is not None:
        procSteps.endStep('Blot')

        
# Run 'drizzle' here...
#
def buildBlotParamDict(configObj):
    blot_name = util.getSectionName(configObj,_blot_step_num_)

    paramDict = {'blot_interp':configObj[blot_name]['blot_interp'],
                'blot_sinscl':configObj[blot_name]['blot_sinscl']}
    return paramDict

def _setDefaults(configObj={}):
    """set up the default parameters to run drizzle
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
    _versions = {'PyDrizzle':util.__version__,'PyFITS':util.__pyfits_version__,'Numpy':util.__numpy_version__}

    _hdrlist = []

    print 'MultiDrizzle: blot task started at ',_ptime()
    
    for img in imageObjectList:
        
        for chip in img.returnAllChips(extname=img.scienceExt):

            print '    Blot: creating blotted image: ',chip.outputNames['data']
            
            _insci = np.zeros((img.outputValues['outny'],img.outputValues['outnx']),dtype=np.float32)
            _outsci = np.zeros((chip.wcs.naxis2,chip.wcs.naxis1),dtype=np.float32)

            #### Check to see what names need to be included here for use in _hdrlist
            chip.outputNames['driz_version'] = _versions
            outputvals = chip.outputNames.copy()
            outputvals.update(img.outputValues)
            outputvals['blotnx'] = chip.wcs.naxis1
            outputvals['blotny'] = chip.wcs.naxis2
            _hdrlist.append(outputvals)

            plist = outputvals.copy()
            plist.update(paramDict)
            
            _data = img.outputNames['outMedian']
            
            # The following type of logic belongs in the user callable (modular)
            # interface to the blot routine, as the user may not want to simply
            # blot the median image, but rather another image altogether (like 
            # the single_drizzle product). 
            #
            # Determine which product was created and should be blotted back
            #if plist['outsingle'] != plist['outdata']:
            #    _data = plist['outsingle']
            #else:
            #    _data = plist['outdata']

            # PyFITS can be used here as it will always operate on
            # output from PyDrizzle (which will always be a FITS file)
            # Open the input science file
            _fname,_sciextn = fileutil.parseFilename(_data)
            _inimg = fileutil.openImage(_fname)

            # Return the PyFITS HDU corresponding to the named extension
            _scihdu = fileutil.getExtn(_inimg,_sciextn)
            _insci = _scihdu.data.copy()

            # DGEO arrays are assumed to be included in the WCS specification
            # Read in the distortion correction arrays, if specified
            #_pxg,_pyg = plist['exposure'].getDGEOArrays()

            # Now pass numpy objects to callable version of Blot...
            #runBlot(plist)
            build=False
            misval = 0.0
            kscale = 1.0

            xmin = 1
            xmax = img.outputValues['outnx']
            ymin = 1
            ymax = img.outputValues['outny']
            # compute the undistorted 'natural' plate scale for this chip
            wcslin = distortion.utils.undistortWCS(chip.wcs)

            if wcsmap is None and arrdriz is not None:
                """
                # Use default C mapping function
                #
                # Convert shifts to input units
                #
                xsh = plist['xsh'] * img.outputValues['scale']
                ysh = plist['ysh'] * img.outputValues['scale']
                # ARRDRIZ.TBLOT needs to be updated to support 'poly5' interpolation,
                # and exptime scaling of output image.
                #
                if (_insci.dtype > np.float32):
                    #WARNING: Input array recast as a float32 array
                    _insci = _insci.astype(np.float32)
                mapping = arrdriz.DefaultMapping(
                    _outsci.shape[1], _outsci.shape[0],
                    _insci.shape[1], _insci.shape[0],
                    xsh, ysh, 'output', 'input', plist['rot'],
                    plist['scale'], 0.0, 0.0, 1.0, 1.0, 0.0,
                    'output', _pxg, _pyg, 'center', plist['coeffs'],
                    None, plist['alpha'], plist['beta'])
                pix_ratio = plist['scale']
                """
                print 'Using default C-based coordinate transformation...'
                wcs_functions.applyShift_to_WCS(img,chip.wcs,output_wcs)
                mapping = arrdriz.DefaultWCSMapping(chip.wcs,output_wcs,int(chip.size1),int(chip.size2),2.0)
                pix_ratio = output_wcs.pscale/wcslin.pscale
            else:
                #
                ##Using the Python class for the WCS-based transformation 
                #
                # Use user provided mapping function
                print 'Using coordinate transformation defined by user...'
                wmap = wcsmap(chip.wcs,output_wcs)
                wmap.applyShift(img)
                mapping = wmap.forward                
                pix_ratio = output_wcs.pscale/wcslin.pscale
                
            t = arrdriz.tblot(
                _insci, _outsci,xmin,xmax,ymin,ymax,
                pix_ratio, kscale, 1.0, 1.0,
                'center',paramDict['blot_interp'], chip._exptime,
                misval, paramDict['blot_sinscl'], 1, mapping)
            del mapping

            # Write output Numpy objects to a PyFITS file
            # Blotting only occurs from a drizzled SCI extension
            # to a blotted SCI extension...
            #_header = fileutil.getHeader(plist['data'])
            #_wcs = wcsutil.WCSObject(plist['data'],header=_header)
            #_wcs = chip.wcs

            _outimg = outputimage.OutputImage(_hdrlist, paramDict, build=False, wcs=chip.wcs, blot=True)
            _outimg.outweight = None
            _outimg.outcontext = None
            _outimg.writeFITS(plist['data'],_outsci,None,versions=_versions)

            #_buildOutputFits(_outsci,None,plist['outblot'])
            _insci *= 0.
            _outsci *= 0.
            _inimg.close()
            del _inimg
            _hdrlist = []

            #del _pxg,_pyg

            del _insci,_outsci
        del _outimg
    
    print 'MultiDrizzle: blot task completed at ',_ptime()
   
