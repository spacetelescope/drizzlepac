"""
Create a median image from the singly drizzled images.

:Authors: Warren Hack, Mihai Cara

:License: :doc:`LICENSE`

"""
import os
import sys
import math
import numpy as np
from astropy.io import fits

from stsci.imagestats import ImageStats
from stsci.image import numcombine
from stsci.tools import iterfile, teal, logutil

from . import imageObject
from . import util
from .minmed import min_med
from . import processInput
from .adrizzle import _single_step_num_

from .version import *

# look in drizzlepac for createMedian.cfg:
__taskname__ = "drizzlepac.createMedian"
_step_num_ = 4  # this relates directly to the syntax in the cfg file

BUFSIZE = 1024*1024   # 1MB cache size

log = logutil.create_logger(__name__, level=logutil.logging.NOTSET)


# this is the user access function
def median(input=None, configObj=None, editpars=False, **inputDict):
    """
    Create a median image from the seperately drizzled images.
    """
    if input is not None:
        inputDict["input"] = input
    else:
        raise ValueError("Please supply an input image")

    configObj = util.getDefaultConfigObj(__taskname__, configObj, inputDict,
                                         loadOnly=(not editpars))
    if configObj is None:
        return

    if not editpars:
        run(configObj)


# this is the function that will be called from TEAL
def run(configObj):
    imgObjList, outwcs = processInput.setCommonInput(
        configObj,
        createOutwcs=False
    )  # outwcs is not needed here
    createMedian(imgObjList, configObj)


# ###################################################
# ## Top-level interface from inside AstroDrizzle  ##
# ###################################################
def createMedian(imgObjList, configObj, procSteps=None):
    """ Top-level interface to createMedian step called from top-level
    AstroDrizzle.

    This function parses the input parameters then calls the `_median()`
    function to median-combine the input images into a single image.

    """
    if imgObjList is None:
        msg = "Please provide a list of imageObjects to the median step"
        print(msg, file=sys.stderr)
        raise ValueError(msg)

    if procSteps is not None:
        procSteps.addStep('Create Median')

    step_name = util.getSectionName(configObj, _step_num_)
    if not configObj[step_name]['median']:
        log.info('Median combination step not performed.')
        return

    paramDict = configObj[step_name]
    paramDict['proc_unit'] = configObj['proc_unit']

    # include whether or not compression was performed
    driz_sep_name = util.getSectionName(configObj, _single_step_num_)
    driz_sep_paramDict = configObj[driz_sep_name]
    paramDict['compress'] = driz_sep_paramDict['driz_sep_compress']

    log.info('USER INPUT PARAMETERS for Create Median Step:')
    util.printParams(paramDict, log=log)

    _median(imgObjList, paramDict)

    if procSteps is not None:
        procSteps.endStep('Create Median')


# this is the internal function, the user called function is below
def _median(imageObjectList, paramDict):
    """Create a median image from the list of image Objects
       that has been given.
    """
    newmasks = paramDict['median_newmasks']
    comb_type = paramDict['combine_type'].lower()
    nlow = paramDict['combine_nlow']
    nhigh = paramDict['combine_nhigh']
    grow = paramDict['combine_grow'] if 'minmed' in comb_type else 0
    maskpt = paramDict['combine_maskpt']
    proc_units = paramDict['proc_unit']
    compress = paramDict['compress']
    bufsizeMB = paramDict['combine_bufsize']

    sigma = paramDict["combine_nsigma"]
    sigmaSplit = sigma.split()
    nsigma1 = float(sigmaSplit[0])
    nsigma2 = float(sigmaSplit[1])

    if paramDict['combine_lthresh'] is None:
        lthresh = None
    else:
        lthresh = float(paramDict['combine_lthresh'])

    if paramDict['combine_hthresh'] is None:
        hthresh = None
    else:
        hthresh = float(paramDict['combine_hthresh'])

    # the name of the output median file isdefined in the output wcs object and
    # stuck in the image.outputValues["outMedian"] dict of every imageObject
    medianfile = imageObjectList[0].outputNames["outMedian"]

    # Build combined array from single drizzled images.

    # Start by removing any previous products...
    if os.access(medianfile, os.F_OK):
        os.remove(medianfile)

    # Define lists for instrument specific parameters, these should be in
    # the image objects need to be passed to the minmed routine
    readnoiseList = []
    exposureTimeList = []
    backgroundValueList = []  # list of  MDRIZSKY *platescale values
    singleDrizList = []  # these are the input images
    singleWeightList = []  # pointers to the data arrays
    wht_mean = []  # Compute the mean value of each wht image

    single_hdr = None
    virtual = None

    # for each image object
    for image in imageObjectList:
        if virtual is None:
            virtual = image.inmemory

        det_gain = image.getGain(1)
        img_exptime = image._image['sci', 1]._exptime
        native_units = image.native_units
        native_units_lc = native_units.lower()

        if proc_units.lower() == 'native':
            if native_units_lc not in ['counts', 'electrons', 'counts/s',
                                       'electrons/s']:
                raise ValueError("Unexpected native units: '{}'"
                                 .format(native_units))

            if lthresh is not None:
                if native_units_lc.startswith('counts'):
                    lthresh *= det_gain
                if native_units_lc.endswith('/s'):
                    lthresh *= img_exptime

            if hthresh is not None:
                if native_units_lc.startswith('counts'):
                    hthresh *= det_gain
                if native_units_lc.endswith('/s'):
                    hthresh *= img_exptime

        singleDriz = image.getOutputName("outSingle")
        singleDriz_name = image.outputNames['outSingle']
        singleWeight = image.getOutputName("outSWeight")
        singleWeight_name = image.outputNames['outSWeight']

        # If compression was used, reference ext=1 as CompImageHDU only writes
        # out MEF files, not simple FITS.
        if compress:
            wcs_ext = '[1]'
            wcs_extnum = 1
        else:
            wcs_ext = '[0]'
            wcs_extnum = 0

        if not virtual:
            if isinstance(singleDriz, str):
                iter_singleDriz = singleDriz + wcs_ext
                iter_singleWeight = singleWeight + wcs_ext
            else:
                iter_singleDriz = singleDriz[wcs_extnum]
                iter_singleWeight = singleWeight[wcs_extnum]
        else:
            iter_singleDriz = singleDriz_name + wcs_ext
            iter_singleWeight = singleWeight_name + wcs_ext

        # read in WCS from first single drizzle image to use as WCS for
        # median image
        if single_hdr is None:
            if virtual:
                single_hdr = singleDriz[wcs_extnum].header
            else:
                single_hdr = fits.getheader(singleDriz_name, ext=wcs_extnum,
                                            memmap=False)

        single_image = iterfile.IterFitsFile(iter_singleDriz)
        if virtual:
            single_image.handle = singleDriz
            single_image.inmemory = True

        singleDrizList.append(single_image)  # add to an array for bookkeeping

        # If it exists, extract the corresponding weight images
        if (not virtual and os.access(singleWeight, os.F_OK)) or (
                virtual and singleWeight):
            weight_file = iterfile.IterFitsFile(iter_singleWeight)
            if virtual:
                weight_file.handle = singleWeight
                weight_file.inmemory = True

            singleWeightList.append(weight_file)
            try:
                tmp_mean_value = ImageStats(weight_file.data, lower=1e-8,
                                            fields="mean", nclip=0).mean
            except ValueError:
                tmp_mean_value = 0.0
            wht_mean.append(tmp_mean_value * maskpt)

            # Extract instrument specific parameters and place in lists

            # If an image has zero exposure time we will
            # redefine that value as '1'.  Although this will cause inaccurate
            # scaling of the data to occur in the 'minmed' combination
            # algorith, this is a necessary evil since it avoids divide by
            # zero exceptions.  It is more important that the divide by zero
            # exceptions not cause AstroDrizzle to crash in the pipeline than
            # it is to raise an exception for this obviously bad data even
            # though this is not the type of data you would wish to process
            # with AstroDrizzle.
            #
            # Get the exposure time from the InputImage object
            #
            # MRD 19-May-2011
            # Changed exposureTimeList to take exposure time from img_exptime
            # variable instead of hte image._exptime attribute, since
            # image._exptime was just giving 1.
            #
            exposureTimeList.append(img_exptime)

            # Use only "commanded" chips to extract subtractedSky and rdnoise:
            rdnoise = 0.0
            nchips = 0
            bsky = None  # minimum sky across **used** chips

            for chip in image.returnAllChips(extname=image.scienceExt):
                # compute sky value as sky/pixel using the single_drz
                # pixel scale:
                if bsky is None or bsky > chip.subtractedSky:
                    bsky = chip.subtractedSky * chip._conversionFactor

                # Extract the readnoise value for the chip
                rdnoise += chip._rdnoise**2
                nchips += 1

            if bsky is None:
                bsky = 0.0

            if nchips > 0:
                rdnoise = math.sqrt(rdnoise/nchips)

            backgroundValueList.append(bsky)
            readnoiseList.append(rdnoise)

            print("reference sky value for image '{}' is {}"
                  .format(image._filename, backgroundValueList[-1]))
        #
        # END Loop over input image list
        #

    # create an array for the median output image, use the size of the first
    # image in the list. Store other useful image characteristics:
    single_driz_data = singleDrizList[0].data
    data_item_size = single_driz_data.itemsize
    single_data_dtype = single_driz_data.dtype
    imrows, imcols = single_driz_data.shape

    medianImageArray = np.zeros_like(single_driz_data)

    del single_driz_data

    if comb_type == "minmed" and not newmasks:
        # Issue a warning if minmed is being run with newmasks turned off.
        print('\nWARNING: Creating median image without the application of '
              'bad pixel masks!\n')

    # The overlap value needs to be set to 2*grow in order to
    # avoid edge effects when scrolling down the image, and to
    # insure that the last section returned from the iterator
    # has enough rows to span the kernel used in the boxcar method
    # within minmed.
    overlap = 2 * grow
    buffsize = BUFSIZE if bufsizeMB is None else (BUFSIZE * bufsizeMB)
    section_nrows = min(imrows, int(buffsize / (imcols * data_item_size)))

    if section_nrows == 0:
        buffsize = imcols * data_item_size
        print("WARNING: Buffer size is too small to hold a single row.\n"
              "         Buffer size size will be increased to minimal "
              "required: {}MB".format(float(buffsize) / 1048576.0))
        section_nrows = 1

    if section_nrows < overlap + 1:
        new_grow = int((section_nrows - 1) / 2)
        if section_nrows == imrows:
            print("'grow' parameter is too large for actual image size. "
                  "Reducing 'grow' to {}".format(new_grow))
        else:
            print("'grow' parameter is too large for requested buffer size. "
                  "Reducing 'grow' to {}".format(new_grow))
        grow = new_grow
        overlap = 2 * grow

    nbr = section_nrows - overlap
    nsec = (imrows - overlap) // nbr
    if (imrows - overlap) % nbr > 0:
        nsec += 1

    for k in range(nsec):
        e1 = k * nbr
        e2 = e1 + section_nrows
        u1 = grow
        u2 = u1 + nbr

        if k == 0:  # first section
            u1 = 0

        if k == nsec - 1:  # last section
            e2 = min(e2, imrows)
            e1 = min(e1, e2 - overlap - 1)
            u2 = e2 - e1

        imdrizSectionsList = np.empty(
            (len(singleDrizList), e2 - e1, imcols),
            dtype=single_data_dtype
        )
        for i, w in enumerate(singleDrizList):
            imdrizSectionsList[i, :, :] = w[e1:e2]

        if singleWeightList:
            weightSectionsList = np.empty(
                (len(singleWeightList), e2 - e1, imcols),
                dtype=single_data_dtype
            )
            for i, w in enumerate(singleWeightList):
                weightSectionsList[i, :, :] = w[e1:e2]
        else:
            weightSectionsList = None

        weight_mask_list = None

        if newmasks and weightSectionsList is not None:
            # Build new masks from single drizzled images.
            # Generate new pixel mask file for median step.
            # This mask will be created from the single-drizzled
            # weight image for this image.

            # The mean of the weight array will be computed and all
            # pixels with values less than 0.7 of the mean will be flagged
            # as bad in this mask. This mask will then be used when
            # creating the median image.
            # 0 means good, 1 means bad here...
            weight_mask_list = np.less(
                weightSectionsList,
                np.asarray(wht_mean)[:, None, None]
            ).astype(np.uint8)

        if 'minmed' in comb_type:  # Do MINMED
            # set up use of 'imedian'/'imean' in minmed algorithm
            fillval = comb_type.startswith('i')

            # Create the combined array object using the minmed algorithm
            result = min_med(
                imdrizSectionsList,
                weightSectionsList,
                readnoiseList,
                exposureTimeList,
                backgroundValueList,
                weight_masks=weight_mask_list,
                combine_grow=grow,
                combine_nsigma1=nsigma1,
                combine_nsigma2=nsigma2,
                fillval=fillval
            )

        else:  # DO NUMCOMBINE
            # Create the combined array object using the numcombine task
            result = numcombine.num_combine(
                imdrizSectionsList,
                masks=weight_mask_list,
                combination_type=comb_type,
                nlow=nlow,
                nhigh=nhigh,
                upper=hthresh,
                lower=lthresh
            )

        # Write out the processed image sections to the final output array:
        medianImageArray[e1+u1:e1+u2, :] = result[u1:u2, :]

    # Write out the combined image
    # use the header from the first single drizzled image in the list
    pf = _writeImage(medianImageArray, inputHeader=single_hdr)

    if virtual:
        mediandict = {}
        mediandict[medianfile] = pf
        for img in imageObjectList:
            img.saveVirtualOutputs(mediandict)
    else:
        try:
            print("Saving output median image to: '{}'".format(medianfile))
            pf.writeto(medianfile)
        except IOError:
            msg = "Problem writing file '{}'".format(medianfile)
            print(msg)
            raise IOError(msg)

    # Always close any files opened to produce median image; namely,
    # single drizzle images and singly-drizzled weight images
    #
    for img in singleDrizList:
        if not virtual:
            img.close()

    # Close all singly drizzled weight images used to create median image.
    for img in singleWeightList:
        if not virtual:
            img.close()


def _writeImage(dataArray=None, inputHeader=None):
    """ Writes out the result of the combination step.
        The header of the first 'outsingle' file in the
        association parlist is used as the header of the
        new image.

        Parameters
        ----------
        dataArray : arr
            Array of data to be written to a fits.PrimaryHDU object

        inputHeader : obj
            fits.header.Header object to use as basis for the PrimaryHDU header

    """
    prihdu = fits.PrimaryHDU(data=dataArray, header=inputHeader)
    pf = fits.HDUList()
    pf.append(prihdu)
    return pf


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
    helpstr = getHelpAsString(docstring=True, show_ver=True)
    if file is None:
        print(helpstr)
    else:
        if os.path.exists(file):
            os.remove(file)
        f = open(file, mode='w')
        f.write(helpstr)
        f.close()


def getHelpAsString(docstring=False, show_ver=True):
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


median.__doc__ = getHelpAsString(docstring=True, show_ver=False)
