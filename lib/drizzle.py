import sys
import util
from util import _ptime

# Default mapping function based on PyWCS 
class WCSMap:
    def __init__(self,input,output):
        # Verify that we have valid WCS input objects
        self.checkWCS(input,'Input')
        self.checkWCS(output,'Output')

        self.input = input
        self.output = output

    def checkWCS(self,obj,name):
        try:
            assert isinstance(obj, pywcs.WCS)
        except AssertionError:
            print name +' object needs to be an instance or subclass of a PyWCS object.'
            raise
    def forward(self,pixx,pixy):
        return self.output.wcs_sky2pix_fits(*self.input.all_pix2sky_fits(pixx,pixy))
    def backward(self,pixx,pixy):
        return self.input.wcs_sky2pix_fits(*self.output.all_pix2sky_fits(pixx,pixy))

#
#### Top-level interface from inside MultiDrizzle
#
def drizSeparate(imageObjectList,output_wcs,configObj={}):
    namelist = []
    hstwcs_list = []
    for img in imageObjectList:
        namelist += img.getKeywordList(outputNames)
        hstwcs_list += img.getKeywordList(wcs)
    parlist = configObj
    run_driz(namelist,hstwcs_list,output_wcs,parlist,wcsmap=WCSMap)
    
def drizFinal(imageObjectList, output_wcs, configObj={}):
    drizSeparate(imageObjectList, output_wcs, configObj={})
    
# Run 'drizzle' here...
#

def run_driz(namelist,hstwcs_list,output_wcs,parlist,wcsmap=WCSMap):
    #(save=no,build=yes,blot=no,single=no,clean=no,interp='linear',sinscl=1.0, debug=no):
    """Perform drizzle operation on input to create output.
     The input parameters originally was a list
     of dictionaries, one for each input, that matches the
     primary parameters for an IRAF drizzle task.

     This method would then loop over all the entries in the
     list and run 'drizzle' for each entry. 
    
    Parameters required for input in parlist:
        build,single,interp,sinscl,debug,
    """    
    print 'MultiDrizzle: drizzle task started at ',_ptime()

    #
    # Setup the versions info dictionary for output to PRIMARY header
    # The keys will be used as the name reported in the header, as-is
    #
    _versions = {'PyDrizzle':__version__,'PyFITS':pyfits.__version__,'Numpy':np.__version__}

    # Interpret input parameters for use in drizzling
    build = parlist['build']
    single= parlist['single']
    interp = parlist['interp']
    sinscl = parlist['sinscl']
    debug = parlist['debug']

    # Store the value of build set by the user for use, if desired,
    # in the 'clean()' method.
    #
    #### This will be replaced with an attribute of the output_wcs object
    output=namelist[0]['outFinal']

    # Check for existance of output file.
    if single == no and build == yes and fileutil.findFile(output):
        print 'Removing previous output product...'
        os.remove(output)

    # Set parameters for each input and run drizzle on it here.
    #
    # Perform drizzling...
    #
    # Only work on a copy of the product WCS, so that when
    # this gets updated for the output image, it does not
    # modify the original WCS computed by PyDrizzle
    #_wcs = observation.product.geometry.wcs.copy()
    _wcs = output_wcs

    _numctx = {'all':len(hstwcs_list)}
    #            if single:
    # Determine how many chips make up each single image
    for plist in namelist:
        plsingle = plist['outSingle']
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
    if single or namelist[0]['outContext'] == '' or namelist[0]['outContext'] == None:
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

    chips = len(hstwcs_list)
    for chip in range(chips):
        # Open the SCI image
        _expname = namelist[chip]['inData']
        _handle = fileutil.openImage(_expname,mode='readonly',memmap=0)
        _extn = 'sci,'+str(namelist[chip]['extver'])
        _sciext = fileutil.getExtn(_handle,extn=_extn)

        ####
        #
        # Put the units keyword handling in the imageObject class
        #
        ####
        # Determine output value of BUNITS
        # and make sure it is not specified as 'ergs/cm...'
        _bunit = sci_chip._bunit
        
        _bindx = _bunit.find('/')

        if parlist['units'] == 'cps':
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
            _mask = namelist[chip]['singleDrizMask']
        else:
            _mask = namelist[chip]['drizMask']

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

        if parlist['wt_scl'] != None:
            if isinstance(parlist['wt_scl'],types.StringType):
                if  parlist['wt_scl'].isdigit() == False :
                    # String passed in as value, check for 'exptime' or 'expsq'
                    _wtscl_float = None
                    try:
                        _wtscl_float = float(parlist['wt_scl'])
                    except ValueError:
                        _wtscl_float = None
                    if _wtscl_float != None:
                        _wtscl = _wtscl_float
                    elif parlist['wt_scl'] == 'expsq':
                        _wtscl = hstwcs_list[chip].exptime*hstwcs_list[chip].exptime
                    else:
                        # Default to the case of 'exptime', if
                        #   not explicitly specified as 'expsq'
                        _wtscl = hstwcs_list[chip].exptime
                else:
                    # int value passed in as a string, convert to float
                    _wtscl = float(parlist['wt_scl'])
            else:
                # We have a non-string value passed in...
                _wtscl = float(parlist['wt_scl'])
        else:
            # Default case: wt_scl = exptime
            _wtscl = hstwcs_list[chip].exptime

        #print 'WT_SCL: ',plist['wt_scl'],' _wtscl: ',_wtscl
        # Set additional parameters needed by 'drizzle'
        _in_units = parlist['in_units']
        if _in_units == 'cps':
            _expin = 1.0
        else:
            _expin = hstwcs_list[chip].exptime
        _shift_fr = 'output'
        _shift_un = 'output'
        _uniqid = _numchips + 1
        ystart = 0
        nmiss = 0
        nskip = 0

        _con = yes
        _imgctx = _numctx['all']
        if single:
            _imgctx = _numctx[namelist[chip]['outSingle']]
        #if single or (plist['outcontext'] == '' and single == yes):
        if _nplanes == 1:
            _con = no
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

        wmap = wcsmap(hstwcs_list[chip],output_wcs)
       
        _vers,nmiss,nskip = arrdriz.tdriz(_sciext.data,_inwht, _outsci, _outwht,
            _outctx[_planeid], _uniqid, ystart, 1, 1, _dny,
            1.0, 1.0, 1.0, 'center', parlist['pixfrac'],
            parlist['kernel'], _in_units, _expin,_wtscl,
            str(parlist['fillval']), nmiss, nskip, 1, wmap.forward)
        del wmap

        # Set up information for generating output FITS image
        namelist[chip]['driz_version'] = _vers
        _hdrlist.append(namelist[chip])

        if nmiss > 0:
            print '! Warning, ',nmiss,' points were outside the output image.'
        if nskip > 0:
            print '! Note, ',nskip,' input lines were skipped completely.'
        # Close image handle
        _handle.close()
        del _handle,_extn,_sciext
        del _inwht

        # Remember the name of the first image that goes into
        # this particular product
        # This will insure that the header reports the proper
        # values for the start of the exposure time used to make
        # this product; in particular, TIME-OBS and DATE-OBS.
        if _numchips == 0:
            _template = namelist[chip]['inData']

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
                _expscale = hstwcs_list[chip].exptime
            else:
                _expscale = output_wcs.exptime

            #If output units were set to 'counts', rescale the array in-place
            if parlist['units'] == 'counts':
                np.multiply(_outsci, _expscale, _outsci)

            #
            # Write output arrays to FITS file(s) and reset chip counter
            #
            _outimg = outputimage.OutputImage(_hdrlist, build=build, wcs=_wcs, single=single)
            _outimg.set_bunit(_bunit)
            _outimg.set_units(parlist['units'])

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
        #del _outsci,_outwht,_inwcs,_outctx, _hdrlist
    # end of loop over each chip


    print 'PyDrizzle drizzling completed at ',_ptime()

