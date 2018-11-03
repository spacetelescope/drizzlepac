"""
Mask blemishes in dithered data by comparison of an image with a model
image and the derivative of the model image.

:Authors: Warren Hack

:License: :doc:`LICENSE`

"""
from __future__ import absolute_import, division, print_function # confidence medium

import os
import re

import numpy as np
from scipy import signal
from astropy.io import fits
from stsci.tools import fileutil, logutil, mputil, teal

from . import quickDeriv
from . import util

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
    crcorr_list =[]
    cr_mask_dict = {}

    for chip in range(1, sciImage._numchips + 1, 1):
        exten = sciImage.scienceExt + ',' +str(chip)
        sci_chip = sciImage[exten]

        if not sci_chip.group_member:
            continue

        blot_image_name = sci_chip.outputNames['blotImage']

        if sciImage.inmemory:
            blot_image = sciImage.virtualOutputs[blot_image_name]
        else:
            if not os.path.isfile(blot_image_name):
                raise IOError("Blotted image not found: {:s}"
                              .format(blot_image_name))

            try:
                blot_image = fits.open(blot_image_name, memmap=False)
            except IOError:
                print("Problem opening blot images")
                raise

        input_image = sciImage.getData(exten)

        # Apply any unit conversions to input image here for comparison
        # with blotted image in units of electrons
        input_image *= sci_chip._conversionFactor

        #make the derivative blot image
        blot_data = blot_image[0].data * sci_chip._conversionFactor
        blot_deriv = quickDeriv.qderiv(blot_data)
        if not sciImage.inmemory:
            blot_image.close()

        # Boolean mask needs to take into account any crbits values
        # specified by the user to be ignored when converting DQ array.
        dq_mask = sciImage.buildMask(chip, paramDict['crbit'])

        #parse out the SNR information
        snr1, snr2 = map(
            float, filter(None, re.split("[,;\s]+", paramDict["driz_cr_snr"]))
        )

        # parse out the scaling information
        mult1, mult2 = map(
            float, filter(None, re.split("[,;\s]+", paramDict["driz_cr_scale"]))
        )

        gain = sci_chip._effGain
        rn = sci_chip._rdnoise
        backg = sci_chip.subtractedSky * sci_chip._conversionFactor

        # Set scaling factor (used by MultiDrizzle) to 1 since scaling has
        # already been accounted for in blotted image
        # expmult = 1.

        ##################   COMPUTATION PART I    ###################
        # Create a temporary array mask
        t1 = np.absolute(input_image - blot_data)
        # ta = np.sqrt(gain * np.abs((blot_data + backg) * expmult) + rn**2)
        ta = np.sqrt(gain * np.abs(blot_data + backg) + rn**2)
        tb = mult1 * blot_deriv + snr1 * ta / gain
        t2 = tb # / expmult
        tmp1 = t1 <= t2

        # Create a convolution kernel that is 3 x 3 of 1's
        kernel = np.ones((3, 3), dtype=np.uint16)
        # Convolve the mask with the kernel
        tmp2 = signal.convolve2d(tmp1, kernel, boundary='symm', mode='same')

        ##################   COMPUTATION PART II    ###################
        # Create the CR Mask
        tb = mult2 * blot_deriv + snr2 * ta / gain
        t2 = tb # / expmult
        cr_mask = (t1 <= t2) | (tmp2 >= 9)

        ##################   COMPUTATION PART III    ##################
        # flag additional cte 'radial' and 'tail' pixels surrounding CR pixels
        # as CRs

        # In both the 'radial' and 'length' kernels below, 0->good and 1->bad,
        # so that upon convolving the kernels with cr_mask, the convolution
        # output will have low->bad and high->good from which 2 new arrays are
        # created having 0->bad and 1->good. These 2 new arrays are then
        # 'anded' to create a new cr_mask.

        # make radial convolution kernel and convolve it with original cr_mask
        cr_grow_kernel = np.ones((grow, grow), dtype=np.uint16)
        cr_grow_kernel_conv = signal.convolve2d(
            cr_mask, cr_grow_kernel, boundary='symm', mode='same'
        )

        # make tail convolution kernel and convolve it with original cr_mask
        cr_ctegrow_kernel = np.zeros((2 * ctegrow + 1, 2 * ctegrow + 1))

        # which pixels are masked by tail kernel depends on sign of
        # sci_chip.cte_dir (i.e.,readout direction):
        if sci_chip.cte_dir == 1:
            # 'positive' direction:  HRC: amp C or D; WFC: chip = sci,1; WFPC2
            cr_ctegrow_kernel[0:ctegrow, ctegrow] = 1
        elif sci_chip.cte_dir == -1:
            # 'negative' direction:  HRC: amp A or B; WFC: chip = sci,2
            cr_ctegrow_kernel[ctegrow+1:2*ctegrow+1, ctegrow] = 1

        # do the convolution
        cr_ctegrow_kernel_conv = signal.convolve2d(
            cr_mask, cr_ctegrow_kernel, boundary='symm', mode='same'
        )

        # select high pixels from both convolution outputs;
        # then 'and' them to create new cr_mask
        cr_grow_mask = cr_grow_kernel_conv >= grow**2 # radial
        cr_ctegrow_mask = cr_ctegrow_kernel_conv >= ctegrow # length
        cr_mask = cr_grow_mask & cr_ctegrow_mask

        # Apply CR mask to the DQ array in place
        dq_mask &= cr_mask

        ####### Create the corr file
        corrFile = np.where(dq_mask, input_image, blot_data)
        corrFile /= sci_chip._conversionFactor
        corrDQMask = np.where(dq_mask, 0, paramDict['crbit']).astype(np.uint16)

        if paramDict['driz_cr_corr']:
            crcorr_list.append({
                'sciext': fileutil.parseExtn(exten),
                'corrFile': corrFile.copy(),
                'dqext': fileutil.parseExtn(sci_chip.dq_extn),
                'dqMask': corrDQMask
            })

        # Save the cosmic ray mask file to disk
        cr_mask_image = sci_chip.outputNames["crmaskImage"]
        if paramDict['inmemory']:
            print('Creating in-memory(virtual) FITS file...')
            _pf = util.createFile(cr_mask.astype(np.uint8),
                                  outfile=None, header=None)
            cr_mask_dict[cr_mask_image] = _pf
            sciImage.saveVirtualOutputs(cr_mask_dict)

        else:
            # Always write out crmaskimage, as it is required input for
            # the final drizzle step. The final drizzle step combines this
            # image with the DQ information on-the-fly.
            #
            # Remove the existing mask file if it exists
            if os.path.isfile(cr_mask_image):
                os.remove(cr_mask_image)
                print("Removed old cosmic ray mask file: '{:s}'"
                      .format(cr_mask_image))
            print("Creating output: {:s}".format(cr_mask_image))
            util.createFile(cr_mask.astype(np.uint8),
                            outfile=cr_mask_image, header=None)

    if paramDict['driz_cr_corr']:
        createCorrFile(sciImage.outputNames["crcorImage"], crcorr_list,
                       sciImage._filename)


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
