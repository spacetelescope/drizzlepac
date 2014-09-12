"""
Interfaces to main drizzle functions.

:Authors: Warren Hack

:License: `<http://www.stsci.edu/resources/software_hardware/pyraf/LICENSE>`_

"""
from __future__ import division # confidence medium

import copy, math, os, sys, time
import util
import numpy as np
from astropy.io import fits
from stsci.tools import fileutil, logutil, mputil, teal
import outputimage, wcs_functions, processInput, util
import stwcs
from stwcs import distortion

from .version import *

try:
    import cdriz
except ImportError:
    cdriz = None
    print '\n Coordinate transformation and image resampling library, cdriz, NOT found!'
    print '\n Please check the installation of this package to insure C code was built successfully.'
    raise ImportError

if util.can_parallel:
    import multiprocessing

__taskname__ = "drizzlepac.adrizzle"
_single_step_num_ = 3
_final_step_num_ = 7

log = logutil.create_logger(__name__)

time_pre_all = []
time_driz_all = []
time_post_all = []
time_write_all = []

#
#### Interactive interface for running drizzle tasks separately
#

def drizzle(input, outdata, wcsmap=None, editpars=False, configObj=None, **input_dict):

    # Pass along values of input and outdata as members of input_dict
    if input_dict is None:
        input_dict = {}
    input_dict['input'] = input
    input_dict['outdata'] = outdata

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
####  User level interface to run drizzle tasks from TEAL
#
def run(configObj, wcsmap=None):
    """ Interface for running `wdrizzle` from TEAL or Python command-line.

    This code performs all file ``I/O`` to set up the use of the drizzle code for
    a single exposure to replicate the functionality of the original `wdrizzle`.
    """

    # Insure all output filenames specified have .fits extensions
    if configObj['outdata'][-5:] != '.fits': configObj['outdata'] += '.fits'
    if not util.is_blank(configObj['outweight']) and configObj['outweight'][-5:] != '.fits': configObj['outweight'] += '.fits'
    if not util.is_blank(configObj['outcontext']) and configObj['outcontext'][-5:] != '.fits': configObj['outcontext'] += '.fits'

    # Keep track of any files we need to open
    in_sci_handle = None
    in_wht_handle = None
    out_sci_handle = None
    out_wht_handle = None
    out_con_handle = None

    _wcskey = configObj['wcskey']
    if util.is_blank(_wcskey):
        _wcskey = ' '

    scale_pars = configObj['Data Scaling Parameters']
    user_wcs_pars = configObj['User WCS Parameters']

    # Open the SCI (and WHT?) image
    # read file to get science array
    insci = get_data(configObj['input'])
    expin = fileutil.getKeyword(configObj['input'],scale_pars['expkey'])
    in_sci_phdr = fits.getheader(fileutil.parseFilename(configObj['input'])[0])

    # we need to read in the input WCS
    input_wcs = stwcs.wcsutil.HSTWCS(configObj['input'],wcskey=_wcskey)

    if not util.is_blank(configObj['inweight']):
        inwht = get_data(configObj['inweight']).astype(np.float32)
    else:
        # Generate a default weight map of all good pixels
        inwht = np.ones(insci.shape,dtype=insci.dtype)

    output_exists = False
    outname = fileutil.osfn(fileutil.parseFilename(configObj['outdata'])[0])
    if os.path.exists(outname):
        output_exists = True
    # Output was specified as a filename, so open it in 'update' mode
    outsci = get_data(configObj['outdata'])

    if output_exists:
        # we also need to read in the output WCS from pre-existing output
        output_wcs = stwcs.wcsutil.HSTWCS(configObj['outdata'])

        out_sci_hdr = fits.getheader(outname)
        outexptime = out_sci_hdr['DRIZEXPT']
        if 'ndrizim' in out_sci_hdr:
            uniqid = out_sci_hdr['ndrizim']+1
        else:
            uniqid = 1

    else:  # otherwise, define the output WCS either from user pars or refimage
        if util.is_blank(configObj['User WCS Parameters']['refimage']):
            # Define a WCS based on user provided WCS values
            # NOTE:
            #   All parameters must be specified, not just one or a few
            if not util.is_blank(user_wcs_pars['outscale']):
                output_wcs = wcs_functions.build_hstwcs(
                    user_wcs_pars['raref'], user_wcs_pars['decref'],
                    user_wcs_pars['xrefpix'], user_wcs_pars['yrefpix'],
                    user_wcs_pars['outnx'], user_wcs_pars['outny'],
                    user_wcs_pars['outscale'], user_wcs_pars['orient'] )
            else:
                # Define default WCS based on input image
                applydist = True
                if input_wcs.sip is None:
                    applydist = False
                output_wcs = stwcs.distortion.utils.output_wcs([input_wcs],undistort=applydist)
        else:
            # Define the output WCS based on a user specified reference image WCS
            output_wcs = stwcs.wcsutil.HSTWCS(configObj['User WCS Parameters']['refimage'])

        # Initialize values used for combining results
        outexptime = 0.0
        uniqid = 1

    # Set up the output data array and insure that the units for that array is 'cps'
    if outsci is None:
        # Define a default blank array based on definition of output_wcs
        outsci = np.zeros((output_wcs._naxis2,output_wcs._naxis1),dtype=np.float32)
    else:
        # Convert array to units of 'cps', if needed
        if outexptime != 0.0:
            np.divide(outsci, outexptime, outsci)
        outsci = outsci.astype(np.float32)

    # Now update output exposure time for additional input file
    outexptime += expin

    outwht = None
    if not util.is_blank(configObj['outweight']):
        outwht = get_data(configObj['outweight'])

    if outwht is None:
        outwht = np.zeros((output_wcs._naxis2,output_wcs._naxis1),dtype=np.float32)
    else:
        outwht = outwht.astype(np.float32)

    outcon = None
    keep_con = False

    if not util.is_blank(configObj['outcontext']):
        outcon = get_data(configObj['outcontext'])
        keep_con = True
        if outcon is None:
            outcon = np.zeros((1,output_wcs._naxis2,output_wcs._naxis1),dtype=np.int32)
        else:
            outcon = outcon.astype(np.int32)
    # Interpret wt_scl parameter
    if configObj['wt_scl'] == 'exptime':
        wt_scl = expin
    elif configObj['wt_scl'] == 'expsq':
        wt_scl = expin*expin
    else:
        wt_scl = float(configObj['wt_scl'])

    # Interpret coeffs parameter to determine whether to apply coeffs or not
    undistort = True
    if not configObj['coeffs'] or input_wcs.sip is None:
        undistort = False
    # turn off use of coefficients if undistort is False (coeffs == False)
    if not undistort:
        input_wcs.sip = None
        input_wcs.cpdis1 = None
        input_wcs.cpdis2 = None
        input_wcs.det2im = None

    wcslin = distortion.utils.output_wcs([input_wcs],undistort=undistort)

    # Perform actual drizzling now...
    _vers = do_driz(insci, input_wcs, inwht,
            output_wcs, outsci, outwht, outcon,
            expin, scale_pars['in_units'],
            wt_scl, wcslin_pscale=wcslin.pscale ,uniqid=uniqid,
            pixfrac=configObj['pixfrac'], kernel=configObj['kernel'],
            fillval=scale_pars['fillval'], stepsize=configObj['stepsize'],
            wcsmap=None, rqstd_num_cores=configObj.get('num_cores'))

    out_sci_handle,outextn = create_output(configObj['outdata'])
    if not output_exists:
        # Also, define default header based on input image Primary header
        out_sci_handle[outextn].header = in_sci_phdr.copy()

    # Update header of output image with exptime used to scale the output data
    # if out_units is not counts, this will simply be a value of 1.0
    # the keyword 'exptime' will always contain the total exposure time
    # of all input image regardless of the output units
    out_sci_handle[outextn].header['EXPTIME'] = outexptime

    # create CTYPE strings
    ctype1 = input_wcs.wcs.ctype[0]
    ctype2 = input_wcs.wcs.ctype[1]
    if ctype1.find('-SIP'): ctype1 = ctype1.replace('-SIP','')
    if ctype2.find('-SIP'): ctype2 = ctype2.replace('-SIP','')

    # Update header with WCS keywords
    out_sci_handle[outextn].header['ORIENTAT'] = output_wcs.orientat
    out_sci_handle[outextn].header['CD1_1'] = output_wcs.wcs.cd[0][0]
    out_sci_handle[outextn].header['CD1_2'] = output_wcs.wcs.cd[0][1]
    out_sci_handle[outextn].header['CD2_1'] = output_wcs.wcs.cd[1][0]
    out_sci_handle[outextn].header['CD2_2'] = output_wcs.wcs.cd[1][1]
    out_sci_handle[outextn].header['CRVAL1'] = output_wcs.wcs.crval[0]
    out_sci_handle[outextn].header['CRVAL2'] = output_wcs.wcs.crval[1]
    out_sci_handle[outextn].header['CRPIX1'] = output_wcs.wcs.crpix[0]
    out_sci_handle[outextn].header['CRPIX2'] = output_wcs.wcs.crpix[1]
    out_sci_handle[outextn].header['CTYPE1'] = ctype1
    out_sci_handle[outextn].header['CTYPE2'] = ctype2
    out_sci_handle[outextn].header['VAFACTOR'] = 1.0


    if scale_pars['out_units'] == 'counts':
        np.multiply(outsci, outexptime, outsci)
        out_sci_handle[outextn].header['DRIZEXPT'] = outexptime

    else:
        out_sci_handle[outextn].header['DRIZEXPT'] = 1.0

    # Update header keyword NDRIZIM to keep track of how many images have
    # been combined in this product so far
    out_sci_handle[outextn].header['NDRIZIM'] = uniqid

    #define keywords to be written out to product header
    drizdict = outputimage.DRIZ_KEYWORDS.copy()

    # Update drizdict with current values
    drizdict['VER']['value'] = _vers[:44]
    drizdict['DATA']['value'] = configObj['input'][:64]
    drizdict['DEXP']['value'] = expin
    drizdict['OUDA']['value'] = configObj['outdata'][:64]
    drizdict['OUWE']['value'] = configObj['outweight'][:64]
    drizdict['OUCO']['value'] = configObj['outcontext'][:64]
    drizdict['MASK']['value'] = configObj['inweight'][:64]
    drizdict['WTSC']['value'] = wt_scl
    drizdict['KERN']['value'] = configObj['kernel']
    drizdict['PIXF']['value'] = configObj['pixfrac']
    drizdict['OUUN']['value'] = scale_pars['out_units']
    drizdict['FVAL']['value'] = scale_pars['fillval']
    drizdict['WKEY']['value'] = configObj['wcskey']
    outputimage.writeDrizKeywords(out_sci_handle[outextn].header,uniqid,drizdict)

    # add output array to output file
    out_sci_handle[outextn].data = outsci
    out_sci_handle.close()

    if not util.is_blank(configObj['outweight']):
        out_wht_handle,outwhtext = create_output(configObj['outweight'])
        out_wht_handle[outwhtext].header = out_sci_handle[outextn].header.copy()
        out_wht_handle[outwhtext].data = outwht
        out_wht_handle.close()

    if keep_con:
        out_con_handle,outconext = create_output(configObj['outcontext'])
        out_con_handle[outconext].data = outcon
        out_con_handle.close()


#
# drizzlepac based interfaces: relying on imageObject instances and drizzlepac internals
#
#
#### Top-level interface from inside MultiDrizzle
#
def drizSeparate(imageObjectList,output_wcs,configObj,wcsmap=None,procSteps=None):
    if procSteps is not None:
        procSteps.addStep('Separate Drizzle')

    # ConfigObj needs to be parsed specifically for driz_separate set of parameters
    single_step = util.getSectionName(configObj,_single_step_num_)
    # This can be called directly from MultiDrizle, so only execute if
    # switch has been turned on (no guarantee MD will check before calling).
    if configObj[single_step]['driz_separate']:
        paramDict = buildDrizParamDict(configObj)
        paramDict['crbit'] = None
        paramDict['proc_unit'] = 'electrons'
        paramDict['wht_type'] = None
        # Force 'build' to always be False, so that this step always generates
        # simple FITS files as output for compatibility with 'createMedian'
        paramDict['build'] = False
        # Record whether or not intermediate files should be deleted when finished
        paramDict['clean'] = configObj['STATE OF INPUT FILES']['clean']
        paramDict['num_cores'] = configObj.get('num_cores')

        log.info('USER INPUT PARAMETERS for Separate Drizzle Step:')
        util.printParams(paramDict, log=log)

        # override configObj[build] value with the value of the build parameter
        # this is necessary in order for MultiDrizzle to always have build=False
        # for single-drizzle step when called from the top-level.
        run_driz(imageObjectList, output_wcs.single_wcs, paramDict, single=True,
                 build=False, wcsmap=wcsmap)
    else:
        log.info('Single drizzle step not performed.')

    if procSteps is not None:
        procSteps.endStep('Separate Drizzle')


def drizFinal(imageObjectList, output_wcs, configObj,build=None,wcsmap=None,procSteps=None):

    if procSteps is not None:
        procSteps.addStep('Final Drizzle')
    # ConfigObj needs to be parsed specifically for driz_final set of parameters
    final_step = util.getSectionName(configObj,_final_step_num_)
    # This can be called directly from MultiDrizle, so only execute if
    # switch has been turned on (no guarantee MD will check before calling).
    if configObj[final_step]['driz_combine']:
        paramDict = buildDrizParamDict(configObj,single=False)
        paramDict['crbit'] = configObj['crbit']
        paramDict['proc_unit'] = configObj['proc_unit']
        paramDict['wht_type'] = configObj[final_step]['final_wht_type']
        paramDict['num_cores'] = configObj.get('num_cores')

        # override configObj[build] value with the value of the build parameter
        # this is necessary in order for MultiDrizzle to always have build=False
        # for single-drizzle step when called from the top-level.
        if build is None:
            build = paramDict['build']
        # Record whether or not intermediate files should be deleted when finished
        paramDict['clean'] = configObj['STATE OF INPUT FILES']['clean']

        log.info('USER INPUT PARAMETERS for Final Drizzle Step:')
        util.printParams(paramDict, log=log)

        run_driz(imageObjectList, output_wcs.final_wcs, paramDict, single=False,
                 build=build, wcsmap=wcsmap)
    else:
        log.info('Final drizzle step not performed.')

    if procSteps is not None:
        procSteps.endStep('Final Drizzle')

# Run 'drizzle' here...
#
def mergeDQarray(maskname,dqarr):
    """ Merge static or CR mask with mask created from DQ array on-the-fly here.
    """
    maskarr = None
    if maskname is not None:
        if isinstance(maskname, str):
            # working with file on disk (default case)
            if os.path.exists(maskname):
                mask = fileutil.openImage(maskname)
                maskarr = mask[0].data.astype(np.bool)
                mask.close()
        else:
            if isinstance(maskname, fits.HDUList):
                # working with a virtual input file
                maskarr = maskname[0].data.astype(np.bool)
            else:
                maskarr = maskname.data.astype(np.bool)

        if maskarr is not None:
            # merge array with dqarr now
            np.bitwise_and(dqarr,maskarr,dqarr)

def updateInputDQArray(dqfile,dq_extn,chip, crmaskname,cr_bits_value):
    if not isinstance(crmaskname, fits.HDUList) and not os.path.exists(crmaskname):
        log.warning('No CR mask file found! Input DQ array not updated.')
        return
    if cr_bits_value == None:
        log.warning('Input DQ array not updated!')
        return
    if isinstance(crmaskname, fits.HDUList):
        # in_memory case
        crmask = crmaskname
    else:
        crmask = fileutil.openImage(crmaskname)

    if os.path.exists(dqfile):
        fullext=dqfile+"["+dq_extn+str(chip)+"]"
        infile = fileutil.openImage(fullext,mode='update')
        __bitarray = np.logical_not(crmask[0].data).astype(np.int16) * cr_bits_value
        np.bitwise_or(infile[dq_extn,chip].data,__bitarray,infile[dq_extn,chip].data)
        infile.close()
        crmask.close()

def buildDrizParamDict(configObj,single=True):
    chip_pars = ['units','wt_scl','pixfrac','kernel','fillval','bits']
    cfunc_pars = {'pixfrac':float}

    # Initialize paramDict with global parameter(s)
    paramDict = {'build':configObj['build'],'stepsize':configObj['stepsize'],
                'coeffs':configObj['coeffs'],'wcskey':configObj['wcskey']}

    # build appro
    if single:
        driz_prefix = 'driz_sep_'
        stepnum = 3
    else:
        driz_prefix = 'final_'
        stepnum = 7
    section_name = util.getSectionName(configObj,stepnum)
    # Copy values from configObj for the appropriate step to paramDict
    for p in configObj[section_name].keys()+[driz_prefix+'units']:
        if p.startswith(driz_prefix):
            par = p[len(driz_prefix):]
            if par == 'units':
                if single:
                    # Hard-code single-drizzle to always returns 'cps'
                    paramDict[par] = 'cps'
                else:
                    paramDict[par] = configObj[section_name][driz_prefix+par]
            else:
                val = configObj[section_name][driz_prefix+par]
                if par in cfunc_pars:
                    val = cfunc_pars[par](val)
                paramDict[par] = val
    return paramDict

def _setDefaults(configObj={}):
    """set up the default parameters to run drizzle
       build,single,units,wt_scl,pixfrac,kernel,fillval,
       rot,scale,xsh,ysh,blotnx,blotny,outnx,outny,data

       Used exclusively for unit-testing, if any are defined.

    """

    paramDict={"build":True,
              "single":True,
              "stepsize":10,
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

def run_driz(imageObjectList,output_wcs,paramDict,single,build,wcsmap=None):
    """ Perform drizzle operation on input to create output.
    The input parameters originally was a list
    of dictionaries, one for each input, that matches the
    primary parameters for an ``IRAF`` `drizzle` task.

    This method would then loop over all the entries in the
    list and run `drizzle` for each entry.

    Parameters required for input in paramDict:
        build,single,units,wt_scl,pixfrac,kernel,fillval,
        rot,scale,xsh,ysh,blotnx,blotny,outnx,outny,data
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

    # Set sub-sampling rate for drizzling
    #stepsize = 2.0
    log.info('  **Using sub-sampling value of %s for kernel %s' %
             (paramDict['stepsize'], paramDict['kernel']))

    outwcs = copy.deepcopy(output_wcs)

    # Check for existance of output file.
    if single == False and build == True and fileutil.findFile(
                                imageObjectList[0].outputNames['outFinal']):
        log.info('Removing previous output product...')
        os.remove(imageObjectList[0].outputNames['outFinal'])

    # print out parameters being used for drizzling
    log.info("Running Drizzle to create output frame with WCS of: ")
    output_wcs.printwcs()

    # Will we be running in parallel?
    pool_size = util.get_pool_size(paramDict.get('num_cores'), len(imageObjectList))
    will_parallel = single and pool_size > 1
    if will_parallel:
        log.info('Executing %d parallel workers' % pool_size)
    else:
        if single: # note only for single driz; final uses different scheme
            log.info('Executing serially')

    # Set parameters for each input and run drizzle on it here.
    #
    # Perform drizzling...

    numctx = 0
    for img in imageObjectList:
        numctx += img._nmembers
    _numctx = {'all':numctx}

    #            if single:
    # Determine how many chips make up each single image
    for img in imageObjectList:
        for chip in img.returnAllChips(extname=img.scienceExt):
            plsingle = chip.outputNames['outSingle']
            if plsingle in _numctx: _numctx[plsingle] += 1
            else: _numctx[plsingle] = 1

    # Compute how many planes will be needed for the context image.
    _nplanes = int((_numctx['all']-1) / 32) + 1
    # For single drizzling or when context is turned off,
    # minimize to 1 plane only...
    if single or imageObjectList[0][1].outputNames['outContext'] in [None,'',' ']:
        _nplanes = 1

    #
    # An image buffer needs to be setup for converting the input
    # arrays (sci and wht) from FITS format to native format
    # with respect to byteorder and byteswapping.
    # This buffer should be reused for each input if possible.
    #
    _outsci = _outwht = _outctx = _hdrlist = None
    if (not single) or \
       ( (single) and (not will_parallel) and (not imageObjectList[0].inmemory) ):
        # Note there are four cases/combinations for single drizzle alone here:
        # (not-inmem, serial), (not-inmem, parallel), (inmem, serial), (inmem, parallel)
        _outsci=np.zeros((output_wcs._naxis2,output_wcs._naxis1),dtype=np.float32)
        _outwht=np.zeros((output_wcs._naxis2,output_wcs._naxis1),dtype=np.float32)
        # initialize context to 3-D array but only pass appropriate plane to drizzle as needed
        _outctx=np.zeros((_nplanes,output_wcs._naxis2,output_wcs._naxis1),dtype=np.int32)
        _hdrlist = []

    # Keep track of how many chips have been processed
    # For single case, this will determine when to close
    # one product and open the next.
    _chipIdx = 0

    # Remember the name of the 1st image that goes into this particular product
    # Insure that the header reports the proper values for the start of the
    # exposure time used to make this; in particular, TIME-OBS and DATE-OBS.
    template = None

    #
    # Work on each image
    #
    subprocs = []
    for img in imageObjectList:

        chiplist = img.returnAllChips(extname=img.scienceExt)

        # How many inputs should go into this product?
        num_in_prod = _numctx['all']
        if single:
            num_in_prod = _numctx[chiplist[0].outputNames['outSingle']]

        # The name of the 1st image
        fnames = []
        for chip in chiplist:
            fnames.append(chip.outputNames['data'])

        if _chipIdx == 0:
            template = fnames
        else:
            template.extend(fnames)

        # Work each image, possibly in parallel
        if will_parallel:
            # use multiprocessing.Manager only if in parallel and in memory
            if img.inmemory:
                manager = multiprocessing.Manager()
                dproxy = manager.dict(img.virtualOutputs) # copy & wrap it in proxy
                img.virtualOutputs = dproxy

            # parallelize run_driz_img (currently for separate drizzle only)
            p = multiprocessing.Process(target=run_driz_img,
                name='adrizzle.run_driz_img()', # for err msgs
                args=(img,chiplist,output_wcs,outwcs,template,paramDict,
                      single,num_in_prod,build,_versions,_numctx,_nplanes,
                      _chipIdx,None,None,None,None,wcsmap))
            subprocs.append(p)
        else:
            # serial run_driz_img run (either separate drizzle or final drizzle)
            run_driz_img(img,chiplist,output_wcs,outwcs,template,paramDict,
                         single,num_in_prod,build,_versions,_numctx,_nplanes,
                         _chipIdx,_outsci,_outwht,_outctx,_hdrlist,wcsmap)

        # Increment/reset master chip counter
        _chipIdx += len(chiplist)
        if _chipIdx == num_in_prod:
            _chipIdx = 0

    # do the join if we spawned tasks
    if will_parallel:
        mputil.launch_and_wait(subprocs, pool_size) # blocks till all done

    del _outsci,_outwht,_outctx,_hdrlist
    # have looped over each img/chip


#
# Still to check:
#    - why have both output_wcs and outwcs?

def run_driz_img(img,chiplist,output_wcs,outwcs,template,paramDict,single,
                 num_in_prod,build,_versions,_numctx,_nplanes,chipIdxCopy,
                 _outsci,_outwht,_outctx,_hdrlist,wcsmap):
    """ Perform the drizzle operation on a single image.
    This is separated out from :py:func:`run_driz` so as to keep together
    the entirety of the code which is inside the loop over
    images.  See the :py:func:`run_driz` code for more documentation.
    """

    # Check for unintialized inputs
    here = _outsci==None and _outwht==None and _outctx==None
    if _outsci is None:
        _outsci=np.zeros((output_wcs._naxis2,output_wcs._naxis1),dtype=np.float32)
    if _outwht is None:
        _outwht=np.zeros((output_wcs._naxis2,output_wcs._naxis1),dtype=np.float32)
    if _outctx is None:
        _outctx = np.zeros((_nplanes,output_wcs._naxis2,output_wcs._naxis1),dtype=np.int32)
    if _hdrlist is None:
        _hdrlist = []

    # Work on each chip - note that they share access to the arrays above
    for chip in chiplist:
        # See if we will be writing out data
        doWrite = chipIdxCopy == num_in_prod-1

#       debuglog('#chips='+str(chipIdxCopy)+', num_in_prod='+\
#                 str(num_in_prod)+', single='+str(single)+', write='+\
#                 str(doWrite)+', here='+str(here))

        # run_driz_chip
        run_driz_chip(img,chip,output_wcs,outwcs,template,paramDict,
                      single,doWrite,build,_versions,_numctx,_nplanes,
                      chipIdxCopy,_outsci,_outwht,_outctx,_hdrlist,wcsmap)

        # Increment chip counter (also done outside of this function)
        chipIdxCopy += 1

    #
    # Reset for next output image...
    #
    if here:
        del _outsci,_outwht,_outctx,_hdrlist
    elif single:
        np.multiply(_outsci,0.,_outsci)
        np.multiply(_outwht,0.,_outwht)
        np.multiply(_outctx,0,_outctx)
        # this was "_hdrlist=[]", but we need to preserve the var ptr itself
        while len(_hdrlist)>0: _hdrlist.pop()
    # else, these were intended to live and be used beyond this function call

    # img.saveVirtualOutputs() has already been done in run_driz_chip (but
    # only if single and doWrite)


def run_driz_chip(img,chip,output_wcs,outwcs,template,paramDict,single,
                  doWrite,build,_versions,_numctx,_nplanes,_numchips,
                  _outsci,_outwht,_outctx,_hdrlist,wcsmap):
    """ Perform the drizzle operation on a single chip.
    This is separated out from `run_driz_img` so as to keep together
    the entirety of the code which is inside the loop over
    chips.  See the `run_driz` code for more documentation.
    """
    global time_pre_all, time_driz_all, time_post_all, time_write_all

    epoch = time.time()

    # Look for sky-subtracted product
    if os.path.exists(chip.outputNames['outSky']):
        chipextn = '['+chip.header['extname']+','+str(chip.header['extver'])+']'
        _expname = chip.outputNames['outSky']+chipextn
    else:
        # If sky-subtracted product does not exist, use regular input
        _expname = chip.outputNames['data']
    log.info('-Drizzle input: %s' % _expname)

    # Open the SCI image
    _handle = fileutil.openImage(_expname,mode='readonly',memmap=0)
    _sciext = _handle[chip.header['extname'],chip.header['extver']]

    # Apply sky subtraction and unit conversion to input array
    log.info("Applying sky value of %0.6f to %s"%(chip.computedSky,_expname))
    _insci = _sciext.data - chip.computedSky
    _insci *= chip._effGain

    # Set additional parameters needed by 'drizzle'
    _in_units = chip.in_units.lower()
    if _in_units == 'cps':
        _expin = 1.0
    else:
        _expin = chip._exptime

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

    _uniqid = _numchips + 1
    if _nplanes == 1:
        # We need to reset what gets passed to TDRIZ
        # when only 1 context image plane gets generated
        # to prevent overflow problems with trying to access
        # planes that weren't created for large numbers of inputs.
        _uniqid = ((_uniqid-1) % 32) + 1

    # Select which mask needs to be read in for drizzling
    ####
    #
    # Actually need to generate mask file here 'on-demand'
    # and combine it with the static_mask for single_drizzle case...
    #
    ####
    # Build basic DQMask from DQ array and bits value
    dqarr = img.buildMask(chip._chip,bits=paramDict['bits'])

    # get correct mask filenames/objects
    staticMaskName = chip.outputNames['staticMask']
    crMaskName = chip.outputNames['crmaskImage']

    if img.inmemory:
        if staticMaskName in img.virtualOutputs:
            staticMaskName = img.virtualOutputs[staticMaskName]
        if crMaskName in img.virtualOutputs:
            crMaskName = img.virtualOutputs[crMaskName]

    # Merge appropriate additional mask(s) with DQ mask
    if single:
        mergeDQarray(staticMaskName,dqarr)
        if dqarr.sum() == 0:
            log.warning('All pixels masked out when applying static mask!')
    else:
        mergeDQarray(staticMaskName,dqarr)

        if dqarr.sum() == 0:
            log.warning('All pixels masked out when applying static mask!')
        else:
            # Only apply cosmic-ray mask when some good pixels remain after
            # applying the static mask
            mergeDQarray(crMaskName,dqarr)

            if dqarr.sum() == 0:
                log.warning('WARNING: All pixels masked out when applying '
                            'cosmic ray mask to %s' % _expname)
        updateInputDQArray(chip.dqfile,chip.dq_extn,chip._chip,
                           crMaskName, paramDict['crbit'])

    img.set_wtscl(chip._chip,paramDict['wt_scl'])

    pix_ratio = outwcs.pscale/chip.wcslin_pscale

    # Convert mask to a datatype expected by 'tdriz'
    # Also, base weight mask on ERR or IVM file as requested by user
    wht_type = paramDict['wht_type']

    if wht_type == 'ERR':
        _inwht = img.buildERRmask(chip._chip,dqarr,pix_ratio)
    elif wht_type == 'IVM':
        _inwht = img.buildIVMmask(chip._chip,dqarr,pix_ratio)
    elif wht_type == 'EXP':
        _inwht = img.buildEXPmask(chip._chip,dqarr)
    else:  # wht_type == None, used for single drizzle images
        _inwht = chip._exptime * dqarr.astype(np.float32)

    if not(paramDict['clean']):
        # Write out mask file if 'clean' has been turned off
        if single:
            step_mask = 'singleDrizMask'
        else:
            step_mask = 'finalMask'

        _outmaskname = chip.outputNames[step_mask]
        if os.path.exists(_outmaskname): os.remove(_outmaskname)
        pimg = fits.PrimaryHDU(data=_inwht)
        img.saveVirtualOutputs({step_mask:pimg})
        # Only write out mask files if in_memory=False
        if not img.inmemory:
            pimg.writeto(_outmaskname)
            del pimg
            log.info('Writing out mask file: %s' % _outmaskname)

    time_pre = time.time() - epoch; epoch = time.time()
    # New interface to performing the drizzle operation on a single chip/image
    _vers = do_driz(_insci, chip.wcs, _inwht, outwcs, _outsci, _outwht, _outctx,
                _expin, _in_units, chip._wtscl,
                wcslin_pscale=chip.wcslin_pscale, uniqid=_uniqid,
                pixfrac=paramDict['pixfrac'], kernel=paramDict['kernel'],
                fillval=paramDict['fillval'], stepsize=paramDict['stepsize'],
                wcsmap=wcsmap, rqstd_num_cores=paramDict.get('num_cores'),
                for_final=(not single))
    time_driz = time.time() - epoch; epoch = time.time()

    # Set up information for generating output FITS image
    #### Check to see what names need to be included here for use in _hdrlist
    chip.outputNames['driz_version'] = _vers
    chip.outputNames['driz_wcskey'] = paramDict['wcskey']
    outputvals = chip.outputNames.copy()

    # Update entries for names/values based on final output
    outputvals.update(img.outputValues)
    for kw in img.outputNames:
        if kw[:3] == 'out':
            outputvals[kw] = img.outputNames[kw]
    outputvals['exptime'] = chip._exptime
    outputvals['expstart'] = chip._expstart
    outputvals['expend'] = chip._expend

    outputvals['wt_scl_val'] = chip._wtscl

    _hdrlist.append(outputvals)
    time_post = time.time() - epoch; epoch = time.time()

    if doWrite:
        ###########################
        #
        #   IMPLEMENTATION REQUIREMENT:
        #
        # Need to implement scaling of the output image
        # from 'cps' to 'counts' in the case where 'units'
        # was set to 'counts'... 21-Mar-2005
        #
        ###########################

        # Convert output data from electrons/sec to counts/sec as specified
        native_units = img.native_units
        if paramDict['proc_unit'].lower() == 'native' and native_units.lower()[:6] == 'counts':
            np.divide(_outsci, chip._gain, _outsci)
            _bunit = native_units.lower()
            if paramDict['units'] == 'counts':
                indx = _bunit.find('/')
                if indx > 0: _bunit = _bunit[:indx]

        # record IDCSCALE for output to product header
        paramDict['idcscale'] = chip.wcs.idcscale
        #If output units were set to 'counts', rescale the array in-place
        if paramDict['units'] == 'counts':
            #determine what exposure time needs to be used
            # to rescale the product.
            if single:
                _expscale = chip._exptime
            else:
                _expscale = img.outputValues['texptime']
            np.multiply(_outsci, _expscale, _outsci)
        #
        # Write output arrays to FITS file(s)
        #
        if not single:
            img.inmemory = False

        _outimg = outputimage.OutputImage(_hdrlist, paramDict, build=build,
                                          wcs=output_wcs, single=single)
        _outimg.set_bunit(_bunit)
        _outimg.set_units(paramDict['units'])
        # performance-wise, this next line takes 99% of the file-writing time
        outimgs = _outimg.writeFITS(template,_outsci,_outwht,ctxarr=_outctx,
                                    versions=_versions,virtual=img.inmemory)
        del _outimg

        # update imageObject with product in memory
        if single:
            img.saveVirtualOutputs(outimgs)

    # this is after the doWrite
    time_write = time.time() - epoch; epoch = time.time()
#   if False and not single: # turn off all this perf reporting for now
    if not single: # !!!
        time_pre_all.append(time_pre)
        time_driz_all.append(time_driz)
        time_post_all.append(time_post)
        time_write_all.append(time_write)

#       log.info('chip time pre-drizzling:  %6.3f' % time_pre)
#       log.info('chip time drizzling:      %6.3f' % time_driz)
#       log.info('chip time post-drizzling: %6.3f' % time_post)
#       log.info('chip time writing output: %6.3f' % time_write)

        if doWrite:
            tot_pre = sum(time_pre_all)
            tot_driz = sum(time_driz_all)
            tot_post = sum(time_post_all)
            tot_write = sum(time_write_all)
            tot = tot_pre+tot_driz+tot_post+tot_write
            log.info('time:')
            log.info('chip total pre-drizzling:  %7.3fs (%4.1f%%)' % (tot_pre,   (100.*tot_pre/tot)))
            log.info('chip total drizzling:      %7.3fs (%4.1f%%)' % (tot_driz,  (100.*tot_driz/tot)))
            log.info('chip total post-drizzling: %7.3fs (%4.1f%%)' % (tot_post,  (100.*tot_post/tot)))
            log.info('chip total writing output: %7.3fs (%4.1f%%)' % (tot_write, (100.*tot_write/tot)))


def do_driz(insci, input_wcs, inwht,
            output_wcs, outsci, outwht, outcon,
            expin, in_units, wt_scl,
            wcslin_pscale=1.0, uniqid=1, pixfrac=1.0, kernel='square',
            fillval="INDEF", stepsize=10, wcsmap=None,
            rqstd_num_cores=None, for_final=False):
    """
    Core routine for performing 'drizzle' operation on a single input image
    All input values will be Python objects such as ndarrays, instead
    of filenames.
    File handling (input and output) will be performed by calling routine.

    """
    # Insure that the fillval parameter gets properly interpreted for use with tdriz
    if util.is_blank(fillval):
        fillval = 'INDEF'
    else:
        fillval = str(fillval)

    if in_units == 'cps':
        expscale = 1.0
    else:
        expscale = expin

    # Compute what plane of the context image this input would
    # correspond to:
    _planeid = int((uniqid-1) /32)
    # Compute how many planes will be needed for the context image.
    _nplanes = _planeid + 1

    if outcon is not None and (outcon.ndim < 3 or (outcon.ndim == 3 and outcon.shape[0] < _nplanes)):
        # convert context image to 3-D array and pass along correct plane for drizzling
        if outcon.ndim == 3:
            nplanes = outcon.shape[0]+1
        else:
            nplanes = 1
        # We need to expand the context image here to accomodate the addition of
        # this new image
        newcon = np.zeros((nplanes,output_wcs._naxis2,output_wcs._naxis1),dtype=np.int32)
        # now copy original outcon arrays into new array
        if outcon.ndim == 3:
            for n in range(outcon.shape[0]):
                newcon[n] = outcon[n].copy()
        else:
            newcon[0] = outcon.copy()
    else:
        if outcon is None:
            outcon = np.zeros((1,output_wcs._naxis2,output_wcs._naxis1),dtype=np.int32)
            _planeid = 0
        newcon = outcon

    # At this point, newcon will always be a 3-D array, so only pass in
    # correct plane to drizzle code
    outctx = newcon[_planeid]

    pix_ratio = output_wcs.pscale/wcslin_pscale

    if wcsmap is None and cdriz is not None:
        log.info('Using WCSLIB-based coordinate transformation...')
        log.info('stepsize = %s' % stepsize)
        mapping = cdriz.DefaultWCSMapping(
            input_wcs, output_wcs, int(input_wcs._naxis1), int(input_wcs._naxis2),
            stepsize)
    else:
        #
        ##Using the Python class for the WCS-based transformation
        #
        # Use user provided mapping function
        log.info('Using coordinate transformation defined by user...')
        if wcsmap is None:
            wcsmap = wcs_functions.WCSMap
        wmap = wcsmap(input_wcs, output_wcs)
        mapping = cdriz.generate_mapping(
            wmap.forward, int(input_wcs._naxis1), int(input_wcs._naxis2), stepsize)

    _shift_fr = 'output'
    _shift_un = 'output'
    ystart = 0
    nmiss = 0
    nskip = 0
    #
    # This call to 'cdriz.tdriz' uses the new C syntax
    #
    _dny = insci.shape[0]
    # Call 'drizzle' to perform image combination
    if (insci.dtype > np.float32):
        #WARNING: Input array recast as a float32 array
        insci = insci.astype(np.float32)

    pool_size = util.get_pool_size(rqstd_num_cores, None)
    will_parallel = for_final and pool_size > 1
    if will_parallel:
        # debug the whole parallelization setup?
        mode_str = 'in parallel'
        run_sequentially_to_debug = 'ASTRODRIZ_SQTL_DBG' in os.environ # False
        if run_sequentially_to_debug:
            mode_str = 'SEQUENTIALLY in MOCK-parallel (DEBUG)'

        # set tile geometry
        NTILES_X, NTILES_Y = mputil.best_tile_layout(pool_size)
        log.info('Drizzling %d tiles %s, using a %dx%d grid.' % \
                 ((NTILES_X*NTILES_Y), mode_str, NTILES_X, NTILES_Y))
        tilexsz = int(math.ceil(outsci.shape[1]/NTILES_X))
        tileysz = int(math.ceil(outsci.shape[0]/NTILES_Y))

        # Start Manager server to keep track of I/O numpy arrays.  In the future
        # we might consider using multiprocessing.sharedctypes instead of a
        # manager server, as they may be faster in IPC.
        manager = multiprocessing.Manager()

        # We also personally keep track of each lproxy we create - need them at end
        all_lproxies = {}
        subprocs = []

        total_t_launch = time.time()
        for i in range(NTILES_Y):
            for j in range(NTILES_X):

                log.debug('drizzling tile: %d' % ((i*NTILES_X) + j + 1))
                tilesci = outsci[tileysz * i:tileysz * (i+1), tilexsz * j:tilexsz * (j+1)].copy()
                tilewht = outwht[tileysz * i:tileysz * (i+1), tilexsz * j:tilexsz * (j+1)].copy()
                tilectx = outctx[tileysz * i:tileysz * (i+1), tilexsz * j:tilexsz * (j+1)].copy()

                # lproxy: carrier list object used by Manager to sync data between processes
                lproxy = manager.list()
                lproxy.append(None) # tilesci
                lproxy.append(None) # tilewht
                lproxy.append(None) # tilectx
                lproxy.append('') # _vers
                lproxy.append(0)  # nmiss
                lproxy.append(0)  # nskip
                lproxy.append(0)  # iotime
                all_lproxies[(i,j)] = lproxy

                p = multiprocessing.Process(target=parallel_tdriz,
                    name='adrizzle.parallel_tdriz()', # for err msgs
                    args=(insci, inwht, tilesci, tilewht, tilectx, lproxy,
                    uniqid, ystart,
                    # These two arguments are the "logical" offset of the output
                    # image, i.e. the origin of the tile within the larger output
                    # frame.  As they have to do with WCS, they are 1-based.
                    tilexsz * j + 1, tileysz * i + 1, _dny,
                    pix_ratio, 1.0, 1.0, 'center', pixfrac,
                    kernel, in_units, expscale, wt_scl,
                    fillval, 0, 0, 1, mapping)) # nmiss & nskip are output only here
                subprocs.append(p)
                if not run_sequentially_to_debug:
                    p.start()

        t = total_t_wait = time.time()
        total_t_launch = t - total_t_launch

        # The Join: wait until all parallel computations are finished
        if run_sequentially_to_debug:
            # The following loop defeats the entire purpose of parallelizing the
            # call to tdriz(), and would indeed be slower than if will_parallel were
            # False (due just to overhead), but should be useful if one wishes to
            # debug all the multiprocessing mechanisms in place, since this incurs
            # no potential concurrency conflicts.  Run _only_ if debugging.
            for p in subprocs:
                p.start()
                p.join()
        else:
            # All tasks are now running in parallel; wait for them to finish
            for p in subprocs:
                p.join()

        total_t_wait = time.time() - total_t_wait

        # Collect all our results back out of the individual lproxy objects.
        # Apply all the tile array data back onto outsci, outwht, and outctx.
        total_io_inside = 0
        total_t_collect = time.time()
        for key in all_lproxies:
            i,j = key
            # tiles
            tilesci = all_lproxies[key][0]
            tilewht = all_lproxies[key][1]
            tilectx = all_lproxies[key][2]
            # update full image
            outsci[tileysz * i:tileysz * (i+1), tilexsz * j:tilexsz * (j+1)] = tilesci
            outwht[tileysz * i:tileysz * (i+1), tilexsz * j:tilexsz * (j+1)] = tilewht
            outctx[tileysz * i:tileysz * (i+1), tilexsz * j:tilexsz * (j+1)] = tilectx
            # _vers (note: each tile's ver string overwrites previous tiles'. eez ok.)
            _vers = all_lproxies[key][3]
            # accumulations to nmiss, nskip
            nmiss += all_lproxies[key][4]
            nskip += all_lproxies[key][5]
            total_io_inside += all_lproxies[key][6]

        total_t_collect = time.time() - total_t_collect
#       log.info('Time spent launching: %.2f secs' % total_t_launch)
#       log.info('Time spent calculating: %.2f secs' % (total_t_wait-total_io_inside) )
#       log.info('Time spent inside the sub-processes moving data: %.2f secs' % total_io_inside)
#       log.info('Time spent collecting results: %.2f secs' % total_t_collect)
        log.info('Time spent altogether on this chip: %.2f secs' % (total_t_launch+total_t_wait+total_t_collect) )
    else:
        if for_final:
            # note only if final; single uses different scheme & notes elsewhere
            log.info('Executing serially')

        # run serially
        _vers,nmiss,nskip = cdriz.tdriz(insci, inwht, outsci, outwht,
            outctx, uniqid, ystart, 1, 1, _dny,
            pix_ratio, 1.0, 1.0, 'center', pixfrac,
            kernel, in_units, expscale, wt_scl,
            fillval, nmiss, nskip, 1, mapping)

    if nmiss > 0:
        if will_parallel: # this is expected because of the tiling
            log.debug('%s points were outside the output image.' % nmiss)
        else:
            log.warning('! %s points were outside the output image.' % nmiss)
    if nskip > 0:
        log.debug('! Note, %s input lines were skipped completely.' % nskip)

    return _vers


def parallel_tdriz(insci, inwht, tilesci, tilewht, tilectx, lproxy,
                   uniqid, ystart, xmin, ymin, dny,
                   pix_ratio, xscale, yscale,
                   align_str, pixfrac, kernel, in_units,
                   expscale, wt_scl, fillval,
                   nmiss, nskip, vflag, mapping):

    """ Wrap direct call to tdriz.cdriz() for parallel processing of it.
        See docs in do_driz and cdriz for more info on tdriz() itself.
        Try to keep as little logic in here as possible.
    """
    # call tdriz
    _vers,nmiss,nskip = cdriz.tdriz(insci, inwht,
        tilesci, tilewht, tilectx,
        uniqid, ystart, xmin, ymin, dny,
        pix_ratio, xscale, yscale,
        align_str, pixfrac, kernel, in_units,
        expscale, wt_scl, fillval,
        nmiss, nskip, vflag, mapping)

    t = time.time()

    # re-assign lproxy values in sub-process; gets them communicated to the parent
    # every other arg is an input to tdriz so we need not communicate it
    lproxy[0] = tilesci
    lproxy[1] = tilewht
    lproxy[2] = tilectx
    lproxy[3] = _vers
    lproxy[4] = nmiss
    lproxy[5] = nskip
    iotime = time.time() - t
    lproxy[6] = iotime


def get_data(filename):
    fileroot,extn = fileutil.parseFilename(filename)
    extname = fileutil.parseExtn(extn)
    if extname[0] == '': extname = "PRIMARY"
    if os.path.exists(fileroot):
        handle = fileutil.openImage(filename)
        data = handle[extname].data
        handle.close()
    else:
        data = None
    return data

def create_output(filename):
    fileroot,extn = fileutil.parseFilename(filename)
    extname = fileutil.parseExtn(extn)
    if extname[0] == '': extname = "PRIMARY"

    if not os.path.exists(fileroot):
        # We need to create the new file
        pimg = fits.HDUList()
        phdu = fits.PrimaryHDU()
        phdu.header['NDRIZIM'] = 1
        pimg.append(phdu)
        if extn is not None:
            # Create a MEF file with the specified extname
            ehdu = fits.ImageHDU(data=arr)
            ehdu.header['EXTNAME'] = extname[0]
            ehdu.header['EXTVER'] = extname[1]
            pimg.append(ehdu)
        log.info('Creating new output file: %s' % fileroot)
        pimg.writeto(fileroot)
        del pimg
    else:
        log.info('Updating existing output file: %s' % fileroot)

    handle = fits.open(fileroot, mode='update')

    return handle,extname


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
    return useful help from a file in the script directory
    called ``__taskname__.help``

    """
    install_dir = os.path.dirname(__file__)
    taskname = util.base_taskname(__taskname__, __package__)
    htmlfile = os.path.join(install_dir, 'htmlhelp', taskname + '.html')
    helpfile = os.path.join(install_dir, taskname + '.help')

    if docstring or (not docstring and not os.path.exists(htmlfile)):
        if show_ver:
            helpString = os.linesep + \
                ' '.join([__taskname__, 'Version', __version__,
                ' updated on ', __vdate__]) + 2*os.linesep
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


drizzle.__doc__ = getHelpAsString(docstring = True, show_ver = False)
