"""
Mask blemishes in dithered data by comparison of an image with a model
image and the derivative of the model image.

:Authors: Warren Hack

:License: :doc:`LICENSE`

"""
from __future__ import absolute_import, division, print_function # confidence medium

import numpy as np
import stsci.convolve as NC
from astropy.io import fits
import os
from . import quickDeriv
from . import util
from stsci.tools import fileutil, logutil, mputil, teal


if util.can_parallel:
    import multiprocessing

from .version import *

__taskname__= "drizzlepac.drizCR"  # looks in drizzlepac for sky.cfg
_step_num_ = 6  # this relates directly to the syntax in the cfg file


log = logutil.create_logger(__name__, level=logutil.logging.NOTSET)


#this is the user access function
def drizCR(input=None, configObj=None, editpars=False, **inputDict):
    """
        Look for cosmic rays.
    """

    log.debug(inputDict)
    inputDict["input"] = input
    configObj = util.getDefaultConfigObj(__taskname__, configObj, inputDict,
                                         loadOnly=(not editpars))
    if configObj is None:
        return

    if not editpars:
        run(configObj)


#this is the function that will be called from TEAL
def run(configObj):
    # outwcs is not neaded here
    imgObjList,outwcs = processInput.setCommonInput(configObj,
                                                    createOutwcs=False)
    rundrizCR(imgObjList, configObj)


#the final function that calls the workhorse
def rundrizCR(imgObjList,configObj,procSteps=None):

    if procSteps is not None:
        procSteps.addStep('Driz_CR')

    step_name = util.getSectionName(configObj,_step_num_)
    if not configObj[step_name]['driz_cr']:
        log.info('Cosmic-ray identification (driz_cr) step not performed.')
        return
    paramDict = configObj[step_name]
    paramDict['crbit'] = configObj['crbit']
    paramDict['inmemory'] = imgObjList[0].inmemory

    log.info("USER INPUT PARAMETERS for Driz_CR Step:")
    util.printParams(paramDict, log=log)

    # if we have the cpus and s/w, ok, but still allow user to set pool size
    pool_size = util.get_pool_size(configObj.get('num_cores'), len(imgObjList))
    if imgObjList[0].inmemory:
        pool_size = 1 # reason why is output in drizzle step

    subprocs = []
    if pool_size > 1:
        log.info('Executing %d parallel workers' % pool_size)
        for image in imgObjList:
            manager = multiprocessing.Manager()
            mgr = manager.dict({})
            #mgr = manager.dict(image.virtualOutputs)

            p = multiprocessing.Process(target=_drizCr,
                name='drizCR._drizCr()', # for err msgs
                args=(image, mgr, paramDict.dict()))
            subprocs.append(p)
            image.virtualOutputs.update(mgr)
        mputil.launch_and_wait(subprocs, pool_size) # blocks till all done
    else:
        log.info('Executing serially')
        for image in imgObjList:
            _drizCr(image,image.virtualOutputs,paramDict)

    if procSteps is not None:
        procSteps.endStep('Driz_CR')


#the workhorse function
def _drizCr(sciImage, virtual_outputs, paramDict):
    """mask blemishes in dithered data by comparison of an image
    with a model image and the derivative of the model image.

    sciImage is an imageObject which contains the science data
    blotImage is inferred from the sciImage object here which knows the name of its blotted image :)
    chip should be the science chip that corresponds to the blotted image that was sent
    paramDict contains the user parameters derived from the full configObj instance
    dgMask is inferred from the sciImage object, the name of the mask file to combine with the generated Cosmic ray mask

    here are the options you can override in configObj

    gain     = 7               # Detector gain, e-/ADU
    grow     = 1               # Radius around CR pixel to mask [default=1 for 3x3 for non-NICMOS]
    ctegrow  = 0               # Length of CTE correction to be applied
    rn       = 5               # Read noise in electrons
    snr      = "4.0 3.0"       # Signal-to-noise ratio
    scale    = "0.5 0.4"       # scaling factor applied to the derivative
    backg    = 0              # Background value
    expkey   = "exptime"        # exposure time keyword

    blot images are saved out to simple fits files with 1 chip in them
    so for example in ACS, there will be 1 image file with 2 chips that is
    the original image and 2 blotted image files, each with 1 chip

    so I'm imagining calling this function twice, once for each chip,
    but both times with the same original science image file, output files
    and some input (output from previous steps) are referenced in the imageobject
    itself

    """

    grow = paramDict["driz_cr_grow"]
    ctegrow = paramDict["driz_cr_ctegrow"]

#    try:
#        assert(chip is not None), 'Please specify a chip to process for blotting'
#        assert(sciImage is not None), 'Please specify a science image object for blotting'

#    except AssertionError:
#        print "Problem with value of chip or sciImage to drizCR"
#        print sciImage
#        raise # raise orig error
    crcorr_list =[]
    crMaskDict = {}

    for chip in range(1, sciImage._numchips + 1, 1):
        exten = sciImage.scienceExt + ',' +str(chip)
        scienceChip = sciImage[exten]

        if scienceChip.group_member:
            blotImagePar = 'blotImage'
            blotImageName = scienceChip.outputNames[blotImagePar]
            if sciImage.inmemory:
                __blotImage = sciImage.virtualOutputs[blotImageName]
            else:
                try:
                    os.access(blotImageName,os.F_OK)
                except IOError:
                    print("Could not find the Blotted image on disk:",blotImageName)
                    raise # raise orig error

                try:
                    __blotImage = fits.open(blotImageName, mode="readonly", memmap=False)
                except IOError:
                    print("Problem opening blot images")
                    raise

            #blotImageName=scienceChip.outputNames["blotImage"] # input file
            crMaskImage=scienceChip.outputNames["crmaskImage"] # output file
            ctedir=scienceChip.cte_dir

            #check that sciImage and blotImage are the same size?

            #grab the actual image from disk
            __inputImage=sciImage.getData(exten)

            # Apply any unit conversions to input image here for comparison
            # with blotted image in units of electrons
            __inputImage *= scienceChip._conversionFactor

            #make the derivative blot image
            __blotData=__blotImage[0].data*scienceChip._conversionFactor #simple fits
            __blotDeriv = quickDeriv.qderiv(__blotData)
            if not sciImage.inmemory:
                __blotImage.close()

            #this grabs the original dq mask from the science image
            # This mask needs to take into account any crbits values
            # specified by the user to be ignored. A call to the
            # buildMask() method may work better here...
            #__dq = sciImage.maskExt + ',' + str(chip)
            #__dqMask=sciImage.getData(__dq)
            __dqMask = sciImage.buildMask(chip,paramDict['crbit']) # both args are ints

            #parse out the SNR information
            __SNRList=(paramDict["driz_cr_snr"]).split()
            __snr1=float(__SNRList[0])
            __snr2=float(__SNRList[1])

            #parse out the scaling information
            __scaleList = (paramDict["driz_cr_scale"]).split()
            __mult1 = float(__scaleList[0])
            __mult2 = float(__scaleList[1])

            __gain=scienceChip._effGain
            __rn=scienceChip._rdnoise
            __backg = scienceChip.subtractedSky*scienceChip._conversionFactor

            # Define output cosmic ray mask to populate
            __crMask = np.zeros(__inputImage.shape,dtype=np.uint8)

            # Set scaling factor (used by MultiDrizzle) to 1 since scaling has
            # already been accounted for in blotted image
            __expmult = 1.

        ##################   COMPUTATION PART I    ###################
            # Create a temporary array mask
            __t1 = np.absolute(__inputImage - __blotData)
            __ta = np.sqrt(__gain * np.absolute(__blotData * __expmult + __backg * __expmult) + __rn * __rn)
            __tb = ( __mult1 * __blotDeriv + __snr1 * __ta / __gain )
            del __ta
            __t2 = __tb / __expmult
            del __tb
            __tmp1 = np.logical_not(np.greater(__t1, __t2))
            del __t1
            del __t2

            # Create a convolution kernel that is 3 x 3 of 1's
            __kernel = np.ones((3,3),dtype=np.uint8)
            # Create an output tmp file the same size as the input temp mask array
            __tmp2 = np.zeros(__tmp1.shape,dtype=np.int16)
            # Convolve the mask with the kernel
            NC.convolve2d(__tmp1,__kernel,output=__tmp2,fft=0,mode='nearest',cval=0)
            del __kernel
            del __tmp1

        ##################   COMPUTATION PART II    ###################
            # Create the CR Mask
            __xt1 = np.absolute(__inputImage - __blotData)
            __xta = np.sqrt(__gain * np.absolute(__blotData * __expmult + __backg * __expmult) + __rn * __rn)
            __xtb = ( __mult2 *__blotDeriv + __snr2 * __xta / __gain )
            del __xta
            __xt2 = __xtb / __expmult
            del __xtb
            # It is necessary to use a bitwise 'and' to create the mask with numarray objects.
            __crMask = np.logical_not(np.greater(__xt1, __xt2) & np.less(__tmp2,9) )

            del __xt1
            del __xt2
            del __tmp2



        ##################   COMPUTATION PART III    ###################
        #flag additional cte 'radial' and 'tail' pixels surrounding CR pixels as CRs

            # In both the 'radial' and 'length' kernels below, 0->good and 1->bad, so that upon
            # convolving the kernels with __crMask, the convolution output will have low->bad and high->good
            # from which 2 new arrays are created having 0->bad and 1->good. These 2 new arrays are then 'anded'
            # to create a new __crMask.

            # recast __crMask to int for manipulations below; will recast to Bool at end
            __crMask_orig_bool= __crMask.copy()
            __crMask= __crMask_orig_bool.astype( np.int8 )

            # make radial convolution kernel and convolve it with original __crMask
            cr_grow_kernel = np.ones((grow, grow))     # kernel for radial masking of CR pixel
            cr_grow_kernel_conv = __crMask.copy()   # for output of convolution
            NC.convolve2d( __crMask, cr_grow_kernel, output = cr_grow_kernel_conv)

            # make tail convolution kernel and convolve it with original __crMask
            cr_ctegrow_kernel = np.zeros((2*ctegrow+1,2*ctegrow+1))  # kernel for tail masking of CR pixel
            cr_ctegrow_kernel_conv = __crMask.copy()  # for output convolution

            # which pixels are masked by tail kernel depends on sign of ctedir (i.e.,readout direction):
            if ( ctedir == 1 ):  # HRC: amp C or D ; WFC: chip = sci,1 ; WFPC2
                cr_ctegrow_kernel[ 0:ctegrow, ctegrow ]=1    #  'positive' direction
            if ( ctedir == -1 ): # HRC: amp A or B ; WFC: chip = sci,2
                cr_ctegrow_kernel[ ctegrow+1:2*ctegrow+1, ctegrow ]=1    #'negative' direction
            if ( ctedir == 0 ):  # NICMOS: no cte tail correction
                pass

            # do the convolution
            NC.convolve2d( __crMask, cr_ctegrow_kernel, output = cr_ctegrow_kernel_conv)

            # select high pixels from both convolution outputs; then 'and' them to create new __crMask
            where_cr_grow_kernel_conv    = np.where( cr_grow_kernel_conv < grow*grow,0,1 )        # radial
            where_cr_ctegrow_kernel_conv = np.where( cr_ctegrow_kernel_conv < ctegrow, 0, 1 )     # length

            __crMask = np.logical_and( where_cr_ctegrow_kernel_conv, where_cr_grow_kernel_conv) # combine masks
            __crMask = __crMask.astype(np.uint8) # cast back to Bool

            del __crMask_orig_bool
            del cr_grow_kernel
            del cr_grow_kernel_conv
            del cr_ctegrow_kernel
            del cr_ctegrow_kernel_conv
            del where_cr_grow_kernel_conv
            del where_cr_ctegrow_kernel_conv

            # Apply CR mask to the DQ array in place
            np.bitwise_and(__dqMask,__crMask,__dqMask)

            ####### Create the corr file
            __corrFile = np.where(__dqMask, __inputImage, __blotData)
            __corrFile /= scienceChip._conversionFactor
            __corrDQMask = np.where(np.equal(__dqMask,0),
                                    paramDict['crbit'],0).astype(np.uint16)

            if paramDict['driz_cr_corr']:
                crcorr_list.append({'sciext':fileutil.parseExtn(exten),
                                'corrFile':__corrFile.copy(),
                                'dqext':fileutil.parseExtn(scienceChip.dq_extn),
                                'dqMask':__corrDQMask.copy()})


            ######## Save the cosmic ray mask file to disk
            _cr_file = np.zeros(__inputImage.shape,np.uint8)
            _cr_file = np.where(__crMask,1,0).astype(np.uint8)

            if not paramDict['inmemory']:
                outfile = crMaskImage
                # Always write out crmaskimage, as it is required input for
                # the final drizzle step. The final drizzle step combines this
                # image with the DQ information on-the-fly.
                #
                # Remove the existing mask file if it exists
                if(os.access(crMaskImage, os.F_OK)):
                    os.remove(crMaskImage)
                    print("Removed old cosmic ray mask file:",crMaskImage)
                print('Creating output : ',outfile)
            else:
                print('Creating in-memory(virtual) FITS file...')
                outfile = None

            _pf = util.createFile(_cr_file, outfile=outfile, header = None)

            if paramDict['inmemory']:
                crMaskDict[crMaskImage] = _pf

    if paramDict['driz_cr_corr']:
        #util.createFile(__corrFile,outfile=crCorImage,header=None)
        createCorrFile(sciImage.outputNames["crcorImage"],
                        crcorr_list, sciImage._filename)
    del crcorr_list
    if paramDict['inmemory']:
        sciImage.saveVirtualOutputs(crMaskDict)
        virtual_outputs = sciImage.virtualOutputs

#### Create _cor file based on format of original input image
def createCorrFile(outfile, arrlist, template):
    """
    Create a _cor file with the same format as the original input image

    The DQ array will be replaced with the mask array used to create the _cor
    file.
    """
    # Remove the existing cor file if it exists
    if(os.access(outfile, os.F_OK)):
        os.remove(outfile)
        print("Removing old corr file:",outfile)

    ftemplate = fits.open(template, memmap=False)
    for arr in arrlist:
        ftemplate[arr['sciext']].data = arr['corrFile']
        if arr['dqext'][0] != arr['sciext'][0]:
            ftemplate[arr['dqext']].data = arr['dqMask']
    ftemplate.writeto(outfile)
    print('Created CR corrected file: ',outfile)

def setDefaults(configObj={}):
    """ Return a dictionary of the default parameters
        which also been updated with the user overrides.
    """
    gain     = 7               # Detector gain, e-/ADU
    grow     = 1               # Radius around CR pixel to mask [default=1 for 3x3 for non-NICMOS]
    ctegrow  = 0               # Length of CTE correction to be applied
    rn       = 5               # Read noise in electrons
    snr      = "4.0 3.0"       # Signal-to-noise ratio
    scale    = "0.5 0.4"       # scaling factor applied to the derivative
    backg    = 0              # Background value
    expkey   = "exptime"        # exposure time keyword

    paramDict={"gain":gain,
                "grow": grow,
                "ctegrow":ctegrow,
                "rn":rn,
                "snr":snr,
                "scale":scale,
                "backg":backg,
                "expkey":expkey}

    if (len(configObj) != 0):
        for key in configObj:
            paramDict[key]=configObj[key]


    return paramDict


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


drizCR.__doc__ = getHelpAsString(docstring = True, show_ver = False)
