import sys,types,os
import util
from util import _ptime
import numpy as np
from pytools import fileutil
import outputimage,wcs_functions,processInput,util
try:
    import arrdriz
except ImportError:
    arrdriz = None

__taskname__ = "BigBlackBox.drizzle"
_single_step_num_ = 3
_final_step_num_ = 7

#
####  User level interface to run drizzle tasks from TEAL
#
def run(configObj=None,wcsmap=wcs_functions.WCSMap):
    
    # Define list of imageObject instances and output WCSObject instance
    # based on input paramters
    imgObjList,outwcs = processInput.setCommonInput(configObj)

    # Parse out which mode is to be run: single drizzle or final drizzle
    # Call only the mode of interest
    single_step = util.getSectionName(configObj,_single_step_num_)
    if configObj[single_step]['driz_separate']:
        drizSeparate(imgObjList,outwcs,configObj,wcsmap=wcsmap)
    else:
        drizFinal(imgObjList,outwcs,configObj,wcsmap=wcsmap)
def getHelpAsString():
    return "Drizzle Help"


# 
#### Interactive interface for running drizzle tasks separately
#
def drizzle(input=None,output=None,configObj=None,wcsmap=wcs_functions.WCSMap,loadOnly=False,**input_dict):
    # Now, merge required input parameters into input_dict
    if input is not None:
        input_dict['input'] = input
    input_dict['output'] = output
    # If called from interactive user-interface, configObj will not be 
    # defined yet, so get defaults using EPAR/TEAL.
    #
    # Also insure that the input_dict (user-specified values) are folded in
    # with a fully populated configObj instance.
    configObj = util.getDefaultConfigObj(__taskname__,configObj,input_dict,loadOnly=loadOnly)
    
    run(configObj,wcsmap=wcsmap)

#
#### Top-level interface from inside MultiDrizzle
#
def drizSeparate(imageObjectList,output_wcs,configObj,wcsmap=wcs_functions.WCSMap):
    # ConfigObj needs to be parsed specifically for driz_separate set of parameters
    single_step = util.getSectionName(configObj,_single_step_num_)
    # This can be called directly from MultiDrizle, so only execute if
    # switch has been turned on (no guarantee MD will check before calling).
    if configObj[single_step]['driz_separate']:
        paramDict = buildDrizParamDict(configObj)
        paramDict['crbit'] = None
        run_driz(imageObjectList, output_wcs.single_wcs, paramDict, single=True, wcsmap=wcsmap)
    
def drizFinal(imageObjectList, output_wcs, configObj,wcsmap=wcs_functions.WCSMap):
    # ConfigObj needs to be parsed specifically for driz_final set of parameters
    final_step = util.getSectionName(configObj,_final_step_num_)
    # This can be called directly from MultiDrizle, so only execute if
    # switch has been turned on (no guarantee MD will check before calling).
    if configObj[final_step]['driz_combine']:
        paramDict = buildDrizParamDict(configObj,single=False)
        paramDict['crbit'] = configObj['crbit']
        run_driz(imageObjectList, output_wcs.final_wcs, paramDict, single=False, wcsmap=wcsmap)

# Run 'drizzle' here...
#
def getWeightMask(maskname,imgObject,chip,bits):
    dqarr = mergeDQarray(maskname,imgObject,chip,bits)
    _inwht = dqarr.astype(np.float32)
    return _inwht

def mergeDQarray(maskname,imageObject,chip,bits):
    """ Merge static or CR mask with mask created from DQ array on-the-fly here.
    """
    dqarr = imageObject.buildMask(chip,bits)
    if maskname is not None and os.path.exists(maskname):
        mask = fileutil.openImage(maskname)
        maskarr = mask[0].data
        dqarr = np.bitwise_and(dqarr,maskarr)
        mask.close()
    return dqarr

def updateInputDQArray(dqname,dq_extn, crmaskname,cr_bits_value):
    if not os.path.exists(crmaskname):
        print 'WARNING: No CR mask file found! Input DQ array not updated.'
        return 
    if cr_bits_value == None:
        print 'WARNING: Input DQ array not updated!'
        return
    crmask = fileutil.openImage(crmaskname)
    if os.path.exists(dqname):
        infile = fileutil.openImage(dqname,mode='update')
        __bitarray = np.logical_not(crmask[0].data).astype(np.int16) * cr_bits_value
        np.bitwise_or(infile[dq_extn].data,__bitarray,infile[dq_extn].data)
        infile.close()
        crmask.close()

def buildDrizParamDict(configObj,single=True):
    chip_pars = ['units','wt_scl','pixfrac','kernel','fillval','bits']
    # Initialize paramDict with global parameter(s)
    paramDict = {'build':configObj['build']}
    # build appro
    if single:
        driz_prefix = 'driz_sep_'
        stepnum = 3
    else:
        driz_prefix = 'final_'
        stepnum = 7
    section_name = util.getSectionName(configObj,stepnum)
    # Copy values from configObj for the appropriate step to paramDict
    for par in chip_pars:
        if par == 'units':
            if single:
                # Hard-code single-drizzle to always returns 'cps'
                paramDict[par] = 'cps'
            else:
                paramDict[par] = configObj[section_name][driz_prefix+par]
        else:
            paramDict[par] = configObj[section_name][driz_prefix+par]
            
    return paramDict
def _setDefaults(configObj={}):
    """set up the default parameters to run drizzle
        build,single,units,wt_scl,pixfrac,kernel,fillval,
        rot,scale,xsh,ysh,blotnx,blotny,outnx,outny,data
        
        Used exclusively for unit-testing, if any are defined.
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

def run_driz(imageObjectList,output_wcs,paramDict,single,wcsmap=None):
    """Perform drizzle operation on input to create output.
     The input parameters originally was a list
     of dictionaries, one for each input, that matches the
     primary parameters for an IRAF drizzle task.

     This method would then loop over all the entries in the
     list and run 'drizzle' for each entry. 
    
    Parameters required for input in paramDict:
        build,single,units,wt_scl,pixfrac,kernel,fillval,
        rot,scale,xsh,ysh,blotnx,blotny,outnx,outny,data
    """    
    print 'MultiDrizzle: drizzle task started at ',_ptime()
    # Insure that input imageObject is a list
    if not isinstance(imageObjectList, list):
        imageObjectList = [imageObjectList]
        
    # Create a list which points to all the chips being combined 
    # by extracting all the chips from each of the input imageObjects
    #chiplist = []
    #for img in imageObjectList:
    #    for chip in range(1,img._numchips+1):
    #        chiplist.append(img._image[img.scienceExt,chip])    

    #
    # Setup the versions info dictionary for output to PRIMARY header
    # The keys will be used as the name reported in the header, as-is
    #
    _versions = {'PyDrizzle':util.__version__,'PyFITS':util.__pyfits_version__,'Numpy':util.__numpy_version__}

    # Interpret input parameters for use in drizzling
    build = paramDict['build']
    crbit = paramDict['crbit']
    bits = paramDict['bits']
    
    # Check for existance of output file.
    if single == False and build == True and fileutil.findFile(imageObjectList[0].outputNames['outFinal']):
        print 'Removing previous output product...'
        os.remove(imageObjectList[0].outputNames['outFinal'])

    # Set parameters for each input and run drizzle on it here.
    #
    # Perform drizzling...
    #

    numctx = 0
    for img in imageObjectList:
        numctx += img._numchips
    _numctx = {'all':numctx}
    
    #            if single:
    # Determine how many chips make up each single image
    for img in imageObjectList:
        for chip in img.returnAllChips(extname=img.scienceExt):
            plsingle = chip.outputNames['outSingle']
            if _numctx.has_key(plsingle): _numctx[plsingle] += 1
            else: _numctx[plsingle] = 1
    #
    # A image buffer needs to be setup for converting the input
    # arrays (sci and wht) from FITS format to native format
    # with respect to byteorder and byteswapping.
    # This buffer should be reused for each input.
    #
    _outsci = np.zeros((output_wcs.naxis2,output_wcs.naxis1),dtype=np.float32)
    _outwht = np.zeros((output_wcs.naxis2,output_wcs.naxis1),dtype=np.float32)

    # Compute how many planes will be needed for the context image.
    _nplanes = int((_numctx['all']-1) / 32) + 1
    # For single drizzling or when context is turned off,
    # minimize to 1 plane only...
    if single or imageObjectList[0][1].outputNames['outContext'] == '' or imageObjectList[0][1].outputNames['outContext'] == None:
        _nplanes = 1

    # Always initialize context images to a 3-D array
    # and only pass the appropriate plane to drizzle as needed
    _outctx = np.zeros((_nplanes,output_wcs.naxis2,output_wcs.naxis1),dtype=np.int32)

    # Keep track of how many chips have been processed
    # For single case, this will determine when to close
    # one product and open the next.
    _numchips = 0
    _nimg = 0
    _hdrlist = []

    for img in imageObjectList:
        for chip in img.returnAllChips(extname=img.scienceExt):
            # Open the SCI image
            _expname = chip.outputNames['data']
            _handle = fileutil.openImage(_expname,mode='readonly',memmap=0)
            #_extn = chip.header['extname']+str(chip.header['extver'])
            #_sciext = fileutil.getExtn(_handle,extn=_extn)
            _sciext = _handle[chip.header['extname'],chip.header['extver']]

            ####
            #
            # Put the units keyword handling in the imageObject class
            #
            ####
            # Determine output value of BUNITS
            # and make sure it is not specified as 'ergs/cm...'
            _bunit = chip._bunit
            
            _bindx = _bunit.find('/')

            if paramDict['units'] == 'cps':
                # If BUNIT value does not specify count rate already...
                if _bindx < 1:
                    # ... append '/SEC' to value
                    _bunit += '/S'
                else:
                    # reset _bunit here to None so it does not
                    #    overwrite what is already in header
                    _bunit = None
            else:
                if _bindx > 0:
                    # remove '/S'
                    _bunit = _bunit[:_bindx]
                else:
                    # reset _bunit here to None so it does not
                    #    overwrite what is already in header
                    _bunit = None

            # Compute what plane of the context image this input would
            # correspond to:
            # _numchips increments as each chip is drizzled
            _planeid = int(_numchips /32)

            # Select which mask needs to be read in for drizzling
            ####
            #
            # Actually need to generate mask file here 'on-demand'
            # and combine it with the static_mask for single_drizzle case...
            #        
            ####
            if single:
                _inwht = getWeightMask(chip.outputNames['staticMask'],img,chip._chip,bits)
            else:
                _inwht = getWeightMask(chip.outputNames['crmaskImage'],img,chip._chip,bits)
                updateInputDQArray(chip.dqname,chip.dq_extn,
                                    chip.outputNames['crmaskImage'],crbit)

            if paramDict['wt_scl'] != None:
                if isinstance(paramDict['wt_scl'],types.StringType):
                    if  paramDict['wt_scl'].isdigit() == False :
                        # String passed in as value, check for 'exptime' or 'expsq'
                        _wtscl_float = None
                        try:
                            _wtscl_float = float(paramDict['wt_scl'])
                        except ValueError:
                            _wtscl_float = None
                        if _wtscl_float != None:
                            _wtscl = _wtscl_float
                        elif paramDict['wt_scl'] == 'expsq':
                            _wtscl = chip._exptime*chip._exptime
                        else:
                            # Default to the case of 'exptime', if
                            #   not explicitly specified as 'expsq'
                            _wtscl = chip._exptime
                    else:
                        # int value passed in as a string, convert to float
                        _wtscl = float(paramDict['wt_scl'])
                else:
                    # We have a non-string value passed in...
                    _wtscl = float(paramDict['wt_scl'])
            else:
                # Default case: wt_scl = exptime
                _wtscl = chip._exptime

            # Set additional parameters needed by 'drizzle'
            _in_units = chip.in_units
            if _in_units == 'cps':
                _expin = 1.0
            else:
                _expin = chip._exptime
            _shift_fr = 'output'
            _shift_un = 'output'
            _uniqid = _numchips + 1
            ystart = 0
            nmiss = 0
            nskip = 0

            _con = True
            _imgctx = _numctx['all']
            if single:
                _imgctx = _numctx[chip.outputNames['outSingle']]

            if _nplanes == 1:
                _con = False
                # We need to reset what gets passed to TDRIZ
                # when only 1 context image plane gets generated
                # to prevent overflow problems with trying to access
                # planes that weren't created for large numbers of inputs.
                _planeid = 0
                _uniqid = ((_uniqid-1) % 32) + 1

            #
            # This call to 'arrdriz.tdriz' uses the new C syntax
            #
            _dny = _sciext.data.shape[0]
            # Call 'drizzle' to perform image combination
            if (_sciext.data.dtype > np.float32):
                #WARNING: Input array recast as a float32 array
                _sciext.data = _sciext.data.astype(np.float32)

            _pxg = np.zeros([2,2],dtype=np.float32)
            _pyg = np.zeros([2,2],dtype=np.float32)


            if wcsmap is None and arrdriz is not None:
                # Use default C mapping function
                _inwcs = np.zeros([8],dtype=np.float64)
                _inwcs = wcs_functions.convertWCS(output_wcs.wcs,_inwcs)
                print 'Default mapping sciext: ',_sciext.data.shape
                print 'Default mapping outsci: ',_outsci.shape

                mapping = arrdriz.DefaultMapping(
                    _sciext.data.shape[1], _sciext.data.shape[0],
                    _outsci.shape[1], _outsci.shape[0],
                    plist['xsh'], plist['ysh'], 'output', 'output',
                    plist['rot'], plist['scale'], 0.0, 0.0, 1.0, 1.0,
                    0.0, 'output', _pxg, _pyg, 'center', plist['coeffs'], _inwcs,
                    plist['alpha'], plist['beta'])

                print 'Default Mapping results: ',mapping(np.array([1,4096]),np.array([1,2048]))
            else:
                # Use user provided mapping function
                wmap = wcsmap(chip.wcs,output_wcs)
                mapping = wmap.forward
            
            _vers,nmiss,nskip = arrdriz.tdriz(_sciext.data,_inwht, _outsci, _outwht,
                _outctx[_planeid], _uniqid, ystart, 1, 1, _dny,
                1.0, 1.0, 1.0, 'center', paramDict['pixfrac'],
                paramDict['kernel'], _in_units, _expin,_wtscl,
                str(paramDict['fillval']), nmiss, nskip, 1, mapping)
            

            # Set up information for generating output FITS image
            #### Check to see what names need to be included here for use in _hdrlist
            chip.outputNames['driz_version'] = _vers
            outputvals = chip.outputNames.copy()
            # Update entries for names/values based on final output
            outputvals.update(img.outputValues)
            outputvals.update(img.outputNames)
            _hdrlist.append(outputvals)

            if nmiss > 0:
                print '! Warning, ',nmiss,' points were outside the output image.'
            if nskip > 0:
                print '! Note, ',nskip,' input lines were skipped completely.'
            # Close image handle
            _handle.close()
            del _handle,_sciext
            del _inwht

            # Remember the name of the first image that goes into
            # this particular product
            # This will insure that the header reports the proper
            # values for the start of the exposure time used to make
            # this product; in particular, TIME-OBS and DATE-OBS.
            if _numchips == 0:
                _template = chip.outputNames['data']

            # Increment number of chips processed for single output
            _numchips += 1
            if _numchips == _imgctx:
                ###########################
                #
                #   IMPLEMENTATION REQUIREMENT:
                #
                # Need to implement scaling of the output image
                # from 'cps' to 'counts' in the case where 'units'
                # was set to 'counts'... 21-Mar-2005
                #
                ###########################
                # Start by determining what exposure time needs to be used
                # to rescale the product.
                if single:
                    _expscale = chip._exptime
                else:
                    _expscale = img.outputValues['texptime']

                #If output units were set to 'counts', rescale the array in-place
                if paramDict['units'] == 'counts':
                    np.multiply(_outsci, _expscale, _outsci)

                #
                # Write output arrays to FITS file(s) and reset chip counter
                #                
                _outimg = outputimage.OutputImage(_hdrlist, paramDict, build=build, wcs=output_wcs, single=single)
                _outimg.set_bunit(_bunit)
                _outimg.set_units(paramDict['units'])

                _outimg.writeFITS(_template,_outsci,_outwht,ctxarr=_outctx,versions=_versions)
                del _outimg
                #
                # Reset chip counter for next output image...
                #
                _numchips = 0
                _nimg = 0
                np.multiply(_outsci,0.,_outsci)
                np.multiply(_outwht,0.,_outwht)
                np.multiply(_outctx,0,_outctx)

                _hdrlist = []
            else:
                _nimg += 1

        del _outsci,_outwht,_outctx, _hdrlist
        # end of loop over each chip


    print 'PyDrizzle drizzling completed at ',_ptime()
