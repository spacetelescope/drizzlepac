import sys,types
import util
from util import _ptime
import numpy as np
from pytools import fileutil
import outputimage,imageObject,wcs_functions
try:
    import arrdriz
except ImportError:
    arrdriz = None
#
#### Top-level interface from inside MultiDrizzle
#
def drizSeparate(imageObjectList,output_wcs,configObj={},wcsmap=wcs_functions.WCSMap):
    if configObj['driz_separate']:
        run_driz(imageObjectList,output_wcs,configObj,single=True,wcsmap=wcsmap)
    
def drizFinal(imageObjectList, output_wcs, configObj={},wcsmap=wcs_functions.WCSMap):
    if configObj['driz_combine']:
        run_driz(imageObjectList, output_wcs, configObj,single=False,wcsmap=wcsmap)
    
# Run 'drizzle' here...
#

def _setDefaults(configObj={}):
    """set up the default parameters to run drizzle
        build,single,units,wt_scl,pixfrac,kernel,fillval,
        rot,scale,xsh,ysh,blotnx,blotny,outnx,outny,data
    """

    paramDict={"build":True,
              "single":True,
              "units":"cps",
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
              "data":None }

    if(len(configObj) !=0):
        for key in configObj.keys():
            paramDict[key]=configObk[key]

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

    # Check for existance of output file.
    if single == False and build == True and fileutil.findFile(output_wcs._filename):
        print 'Removing previous output product...'
        os.remove(output_wcs._filename)

    # Set parameters for each input and run drizzle on it here.
    #
    # Perform drizzling...
    #
    # Only work on a copy of the product WCS, so that when
    # this gets updated for the output image, it does not
    # modify the original WCS computed by PyDrizzle
    #_wcs = observation.product.geometry.wcs.copy()
    _wcs = output_wcs.wcs

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
    _outsci = np.zeros((_wcs.naxis2,_wcs.naxis1),dtype=np.float32)
    _outwht = np.zeros((_wcs.naxis2,_wcs.naxis1),dtype=np.float32)

    # Compute how many planes will be needed for the context image.
    _nplanes = int((_numctx['all']-1) / 32) + 1
    # For single drizzling or when context is turned off,
    # minimize to 1 plane only...
    if single or imageObjectList[0][1].outputNames['outContext'] == '' or imageObjectList[0][1].outputNames['outContext'] == None:
        _nplanes = 1

    # Always initialize context images to a 3-D array
    # and only pass the appropriate plane to drizzle as needed
    _outctx = np.zeros((_nplanes,_wcs.naxis2,_wcs.naxis1),dtype=np.int32)

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
                _mask = chip.outputNames['singleDrizMask']
            else:
                _mask = chip.outputNames['drizMask']

            # Check to see whether there is a mask_array at all to use...
            if isinstance(_mask,types.StringType):
                if _mask != None and _mask != '':
                    _wht_handle = fileutil.openImage(_mask,mode='readonly',memmap=0)
                    _inwht = _wht_handle[0].data.astype(np.float32)
                    _wht_handle.close()
                    del _wht_handle
                else:
                    print 'No weight or mask file specified!  Assuming all pixels are good.'
                    _inwht = np.ones((_sciext.data.shape[0],_sciext.data.shape[1]),dtype=np.float32)
            elif _mask != None:
                _inwht = _mask.astype(np.float32)
            else:
                print 'No weight or mask file specified!  Assuming all pixels are good.'
                _inwht = np.ones((_sciext.data.shape[0],_sciext.data.shape[1]),dtype=np.float32)

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
            _in_units = paramDict['in_units']
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
                _inwcs = wcs_functions.convertWCS(_wcs.wcs,_inwcs)
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
                wmap = wcsmap(chip.wcs,output_wcs.wcs)
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
                    _expscale = output_wcs._exptime

                #If output units were set to 'counts', rescale the array in-place
                if paramDict['units'] == 'counts':
                    np.multiply(_outsci, _expscale, _outsci)

                #
                # Write output arrays to FITS file(s) and reset chip counter
                #                
                _outimg = outputimage.OutputImage(_hdrlist, paramDict, build=build, wcs=_wcs, single=single)
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

    
    for img in imageObjectList:
        
        for chip in img.returnAllChips(extname=img.scienceExt):

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
            scale = 1.0

            xmin = 1
            xmax = img.outputValues['outnx']
            ymin = 1
            ymax = img.outputValues['outny']
            if wcsmap is None and arrdriz is not None:
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
            else:
                # Use user provided mapping function
                wmap = wcsmap(chip.wcs,output_wcs.wcs)
                mapping = wmap.forward

            t = arrdriz.tblot(
                _insci, _outsci,xmin,xmax,ymin,ymax,
                scale, kscale, 1.0, 1.0,
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
    
    
