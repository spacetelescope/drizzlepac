"""
Mask blemishes in dithered data by comparison of an image with a model
image and the derivative of the model image.

:Authors: Warren Hack

:License: :doc:`/LICENSE`

"""
import os
import re

import numpy as np
from scipy import signal
from astropy.io import fits
from stsci.tools import fileutil, logutil, mputil


from . import quickDeriv
from . import util
from . import processInput
from . import version as __version__

if util.can_parallel:
    import multiprocessing

__all__ = ['rundrizCR', 'drizCR']


__taskname__ = "drizCR"  # looks in drizzlepac for sky.cfg
STEP_NUM = 6  # this relates directly to the syntax in the cfg file
PROCSTEPS_NAME = "Driz_CR"


log = logutil.create_logger(__name__, level=logutil.logging.NOTSET)


def drizCR(input=None, configObj=None, editpars=False, **inputDict):
    """ 
    The blotted median images that are now transformed back into the original
    reference frame, get compared to the original input images to detect any
    spurious pixels (which may include pixels impacted by cosmic rays) in
    each input. Those spurious pixels then get flagged as 'bad' in the output
    cosmic ray mask files, which get used as input for the final combination
    so that they do not show up in the final product.
    The identified bad pixels get flagged by updating the input mask files.
    Optionally, copies of the original images with the bad pixels removed
    can be created through the use of the ``driz_cr_corr`` parameter.

    Parameters
    ----------

    input : str or list of str (Default = None)
        A Python list of blotted median image filenames, or just a single
        filename.

    configObj : configObject (Default = None)
        An instance of ``configObject`` which overrides default parameter settings.

    editpars : bool (Default = False)
        A parameter that allows user to edit input parameters by hand in the GUI.
        True to use the GUI to edit parameters.

    inputDict : dict, optional
        An optional list of parameters specified by the user, which can also
        be used to override the defaults.


    Notes
    -----
    The following are additional parameters that can be set in the configObj
    for this step:

    driz_cr : bool (Default = False)
        Perform cosmic-ray detection? If set to ``True``, cosmic-rays will be
        detected and used to create cosmic-ray masks based on the algorithms
        from ``deriv`` and ``driz_cr``.

    driz_cr_corr : bool (Default = False)
        Create a cosmic-ray cleaned input image? I set to ``True``, a cosmic-ray
        cleaned ``_cor`` image will be generated directly from the input image,
        and a corresponding ``_crmask`` file will be written to document detected
        pixels affected by cosmic-rays.

    driz_cr_snr : list of floats (Default = '3.5 3.0')
        The values for this parameter specify the signal-to-noise ratios for
        the ``driz_cr`` task to be used in detecting cosmic rays.
        See the help file for ``driz_cr`` for further discussion of this parameter.

    driz_cr_grow : int (Default = 1)
        The radius, in pixels, around each detected cosmic-ray,
        in which more stringent detection criteria for additional cosmic
        rays will be used.

    driz_cr_ctegrow : int (Default = 0)
        Length, in pixels, of the CTE tail that should be masked in
        the drizzled output.

    driz_cr_scale : str (Default = '1.2 0.7')
        Scaling factor applied to the derivative in ``driz_cr`` when detecting
        cosmic-rays. See the help file for ``driz_cr`` for further discussion
        of this parameter.


    These tasks are designed to work together seemlessly when run in the
    full ``AstroDrizzle`` interface. More advanced users may wish to create
    specialized scripts for their own datasets, making use of only a subset
    of the predefined ``AstroDrizzle`` tasks, or add additional processing,
    which may be usefull for their particular data. In these cases, individual
    access to the tasks is important.

    Something to keep in mind is that the full ``AstroDrizzle`` interface will
    make backup copies of your original files and place them in
    the ``OrIg/`` directory of your current working directory. If you are
    working with the stand alone interfaces, it is assumed that the user
    has already taken care of backing up their original datafiles as the
    input file with be directly altered.


    Examples
    --------
    Basic example of how to call ``drizCR`` yourself from a Python command line
    using the default parameters for the task.

    >>> from drizzlepac import drizCR
    >>> drizCR.drizCR('*flt.fits')

    """
    log.debug(inputDict)
    inputDict["input"] = input
    configObj = util.getDefaultConfigObj(__taskname__, configObj, inputDict,
                                         loadOnly=(not editpars))
    if configObj is None:
        return

    if not editpars:
        run(configObj)

#--------------------------------
# TEAL Interface functions
# (these functions are deprecated)
#---------------------------------
def run(configObj):
    # outwcs is not neaded here
    imgObjList, outwcs = processInput.setCommonInput(configObj,
                                                     createOutwcs=False)
    rundrizCR(imgObjList, configObj)


def rundrizCR(imgObjList, configObj, procSteps=None):
    """Execute the ``driz_cr`` cosmic-ray rejection step.

    Parameters
    ----------
    imgObjList : list
        Iterable of image objects created by ``processInput`` that contain the
        science data, masks, and bookkeeping information for each exposure.
    configObj : ConfigObj-like
        Configuration structure holding AstroDrizzle parameter sections. The
        ``driz_cr`` section supplies thresholds and flags for the rejection
        run.
    procSteps : drizzlepac.util.ProcessingSteps, optional
        Optional progress tracker used by the pipeline harness; when provided
        the step is registered and marked complete automatically.

    Returns
    -------
    None

    """

    if procSteps is not None:
        procSteps.addStep(PROCSTEPS_NAME)

    step_name = util.getSectionName(configObj, STEP_NUM)
    if not configObj[step_name]['driz_cr']:
        log.info('Cosmic-ray identification (driz_cr) step not performed.')
        if procSteps is not None:
            procSteps.endStep(PROCSTEPS_NAME, reason="off", delay_msg=True)
        return

    paramDict = configObj[step_name]
    paramDict['crbit'] = configObj['crbit']
    paramDict['inmemory'] = imgObjList[0].inmemory

    log.info(f"USER INPUT PARAMETERS for {PROCSTEPS_NAME} Step:")
    util.printParams(paramDict, log=log)

    # if we have the cpus and s/w, ok, but still allow user to set pool size
    pool_size = util.get_pool_size(configObj.get('num_cores'), len(imgObjList))
    if imgObjList[0].inmemory:
        pool_size = 1  # reason why is output in drizzle step

    subprocs = []
    if pool_size > 1:
        log.info('Executing {:d} parallel workers'.format(pool_size))
        mp_ctx = multiprocessing.get_context('fork')
        for image in imgObjList:
            manager = mp_ctx.Manager()
            mgr = manager.dict({})

            p = mp_ctx.Process(
                target=_driz_cr,
                name='drizCR._driz_cr()',  # for err msgs
                args=(image, mgr, paramDict.dict())
            )
            subprocs.append(p)
            image.virtualOutputs.update(mgr)
        mputil.launch_and_wait(subprocs, pool_size)  # blocks till all done

    else:
        log.info('Executing serially')
        for image in imgObjList:
            _driz_cr(image, image.virtualOutputs, paramDict)

    if procSteps is not None:
        procSteps.endStep(PROCSTEPS_NAME)


def _driz_cr(sciImage, virtual_outputs, paramDict):
    """mask blemishes in dithered data by comparison of an image
    with a model image and the derivative of the model image.

    - ``sciImage`` is an imageObject which contains the science data
    - ``blotImage`` is inferred from the ``sciImage`` object here which knows
        the name of its blotted image
    - ``chip`` should be the science chip that corresponds to the blotted
        image that was sent
    - ``paramDict`` contains the user parameters derived from the full
        ``configObj`` instance
    - ``dqMask`` is inferred from the ``sciImage`` object, the name of the mask
        file to combine with the generated Cosmic ray mask

    Here are the options you can override in ``configObj``

    ``gain``     = 7          # Detector gain, e-/ADU
    ``grow``     = 1          # Radius around CR pixel to mask
                              # [default=1 for 3x3 for non-NICMOS]
    ``ctegrow``  = 0          # Length of CTE correction to be applied
    ``rn``       = 5          # Read noise in electrons
    ``snr``      = "4.0 3.0"  # Signal-to-noise ratio
    ``scale``    = "0.5 0.4"  # scaling factor applied to the derivative
    ``backg``    = 0          # Background value
    ``expkey``   = "exptime"  # exposure time keyword

    Blot images are saved out to simple fits files with 1 chip in them
    so for example in ACS, there will be 1 image file with 2 chips that is
    the original image and 2 blotted image files, each with 1 chip

    So I'm imagining calling this function twice, once for each chip,
    but both times with the same original science image file, output files
    and some input (output from previous steps) are referenced in the
    imageobject itself

    """
    grow = paramDict["driz_cr_grow"]
    ctegrow = paramDict["driz_cr_ctegrow"]
    crcorr_list = []
    cr_mask_dict = {}

    for chip in range(1, sciImage._numchips + 1, 1):
        exten = sciImage.scienceExt + ',' + str(chip)
        sci_chip = sciImage[exten]

        if not sci_chip.group_member:
            continue

        blot_image_name = sci_chip.outputNames['blotImage']

        if sciImage.inmemory:
            blot_data = sciImage.virtualOutputs[blot_image_name][0].data
        else:
            if not os.path.isfile(blot_image_name):
                raise IOError("Blotted image not found: {:s}"
                              .format(blot_image_name))

            try:
                blot_data = fits.getdata(blot_image_name, ext=0)
            except IOError:
                print("Problem opening blot images")
                raise
        # Scale blot image, as needed, to match original input data units.
        blot_data *= sci_chip._conversionFactor

        input_image = sciImage.getData(exten)

        # Apply any unit conversions to input image here for comparison
        # with blotted image in units of electrons
        input_image *= sci_chip._conversionFactor

        # make the derivative blot image
        blot_deriv = quickDeriv.qderiv(blot_data)

        # Boolean mask needs to take into account any crbits values
        # specified by the user to be ignored when converting DQ array.
        dq_mask = sciImage.buildMask(chip, paramDict['crbit'])

        # parse out the SNR information
        snr1, snr2 = map(
            float, filter(None, re.split("[,;\s]+", paramDict["driz_cr_snr"]))
        )

        # parse out the scaling information
        mult1, mult2 = map(
            float, filter(
                None, re.split("[,;\s]+", paramDict["driz_cr_scale"])
            )
        )

        gain = sci_chip._effGain
        rn = sci_chip._rdnoise
        backg = sci_chip.subtractedSky * sci_chip._conversionFactor

        # Set scaling factor (used by MultiDrizzle) to 1 since scaling has
        # already been accounted for in blotted image
        # expmult = 1.

        # #################   COMPUTATION PART I    ###################
        # Create a temporary array mask
        t1 = np.absolute(input_image - blot_data)
        # ta = np.sqrt(gain * np.abs((blot_data + backg) * expmult) + rn**2)
        ta = np.sqrt(gain * np.abs(blot_data + backg) + rn**2)
        t2 = (mult1 * blot_deriv + snr1 * ta / gain)  # / expmult
        tmp1 = t1 <= t2

        # Create a convolution kernel that is 3 x 3 of 1's
        kernel = np.ones((3, 3), dtype=np.uint16)
        # Convolve the mask with the kernel
        tmp2 = signal.convolve2d(tmp1, kernel, boundary='symm', mode='same')

        # #################   COMPUTATION PART II    ###################
        # Create the CR Mask
        t2 = (mult2 * blot_deriv + snr2 * ta / gain)  # / expmult
        cr_mask = (t1 <= t2) | (tmp2 >= 9)

        # #################   COMPUTATION PART III    ##################
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
        cr_grow_mask = cr_grow_kernel_conv >= grow**2  # radial
        cr_ctegrow_mask = cr_ctegrow_kernel_conv >= ctegrow  # length
        cr_mask = cr_grow_mask & cr_ctegrow_mask

        # Apply CR mask to the DQ array in place
        dq_mask &= cr_mask

        # Create the corr file
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


def createCorrFile(outfile, arrlist, template):
    """
    Create a _cor file with the same format as the original input image.

    The DQ array will be replaced with the mask array used to create the _cor
    file.
    """
    # Remove the existing cor file if it exists
    if os.path.isfile(outfile):
        os.remove(outfile)
        print("Removing old corr file: '{:s}'".format(outfile))

    with fits.open(template, memmap=False) as ftemplate:
        for arr in arrlist:
            ftemplate[arr['sciext']].data = arr['corrFile']
            if arr['dqext'][0] != arr['sciext'][0]:
                ftemplate[arr['dqext']].data = arr['dqMask']
        ftemplate.writeto(outfile)
        print("Created CR corrected file: '{:s}'".format(outfile))


def setDefaults(configObj={}):
    """ Return a dictionary of the default parameters
        which also been updated with the user overrides.
    """
    paramDict = {
        'gain': 7,           # Detector gain, e-/ADU
        'grow': 1,           # Radius around CR pixel to mask [default=1 for
                             #     3x3 for non-NICMOS]
        'ctegrow': 0,        # Length of CTE correction to be applied
        'rn': 5,             # Read noise in electrons
        'snr': '4.0 3.0',    # Signal-to-noise ratio
        'scale': '0.5 0.4',  # scaling factor applied to the derivative
        'backg': 0,          # Background value
        'expkey': 'exptime'  # exposure time keyword
    }

    if len(configObj) > 0:
        for key in configObj:
            paramDict[key] = configObj[key]

    return paramDict