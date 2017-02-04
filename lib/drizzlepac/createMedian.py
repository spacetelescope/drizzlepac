"""
Create a median image from the singly drizzled images.

:Authors: Warren Hack

:License: `<http://www.stsci.edu/resources/software_hardware/pyraf/LICENSE>`_

"""

# Import external packages
from __future__ import absolute_import, division, print_function # confidence medium

import sys
import numpy as np
from astropy.io import fits
import os, math

from stsci.imagestats import ImageStats
from stsci.image import numcombine
from stsci.tools import iterfile, nimageiter, teal, logutil

from . import imageObject
from . import util
from .minmed import minmed
from . import processInput
from .adrizzle import _single_step_num_


from .version import *

__taskname__= "drizzlepac.createMedian" #looks in drizzlepac for createMedian.cfg
_step_num_ = 4  #this relates directly to the syntax in the cfg file

log = logutil.create_logger(__name__)


#this is the user access function
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


#this is the function that will be called from TEAL
def run(configObj):

    imgObjList,outwcs = processInput.setCommonInput(configObj,createOutwcs=False) #outwcs is not needed here

    createMedian(imgObjList,configObj)

#
#### Top-level interface from inside MultiDrizzle
#
def createMedian(imgObjList,configObj,procSteps=None):
    """ Top-level interface to createMedian step called from top-level MultiDrizzle.

    This function parses the input parameters then calls the `_median()` function
    to median-combine the input images into a single image.

    """
    if imgObjList is None:
        msg = "Please provide a list of imageObjects to the median step"
        print(msg, file=sys.stderr)
        raise ValueError(msg)

    if procSteps is not None:
        procSteps.addStep('Create Median')

    step_name = util.getSectionName(configObj,_step_num_)
    if not configObj[step_name]['median']:
        log.info('Median combination step not performed.')
        return

    paramDict=configObj[step_name]
    paramDict['proc_unit'] = configObj['proc_unit']

    # include whether or not compression was performed
    driz_sep_name = util.getSectionName(configObj,_single_step_num_)
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
    comb_type = paramDict['combine_type']
    nlow = paramDict['combine_nlow']
    nhigh = paramDict['combine_nhigh']
    grow = paramDict['combine_grow']
    maskpt = paramDict['combine_maskpt']
    proc_units = paramDict['proc_unit']
    compress = paramDict['compress']
    bufsizeMb = paramDict['combine_bufsize']

    sigma=paramDict["combine_nsigma"]
    sigmaSplit=sigma.split()
    nsigma1 = float(sigmaSplit[0])
    nsigma2 = float(sigmaSplit[1])

    #print "Checking parameters:"
    #print comb_type,nlow,nhigh,grow,maskpt,nsigma1,nsigma2
    if paramDict['combine_lthresh'] is None:
        lthresh = None
    else:
        lthresh = float(paramDict['combine_lthresh'])

    if paramDict['combine_hthresh'] is None:
        hthresh = None
    else:
        hthresh = float(paramDict['combine_hthresh'])

    #the name of the output median file isdefined in the output wcs object
    #and stuck in the image.outputValues["outMedian"] dict of every imageObject
    medianfile=imageObjectList[0].outputNames["outMedian"]


    """ Builds combined array from single drizzled images."""
    # Start by removing any previous products...
    if(os.access(medianfile,os.F_OK)):
        os.remove(medianfile)


    # Define lists for instrument specific parameters, these should be in the image objects
    # need to be passed to the minmed routine
    readnoiseList = []
    exposureTimeList = []
    backgroundValueList = [] #list of  MDRIZSKY *platescale values
    singleDrizList=[] #these are the input images
    singleWeightList=[] #pointers to the data arrays
    #skylist=[] #the list of platescale values for the images
    _wht_mean = [] # Compute the mean value of each wht image

    _single_hdr = None
    virtual = None

    #for each image object
    for image in imageObjectList:
        if virtual is None:
            virtual = image.inmemory

        det_gain = image.getGain(1)
        img_exptime = image._image['sci',1]._exptime
        native_units = image.native_units
        if lthresh is not None:
            if proc_units.lower() == 'native':
                if native_units.lower() == "counts":
                    lthresh = lthresh * det_gain
                    if native_units.lower() == "counts/s":
                        lthresh = lthresh * img_exptime
        if hthresh is not None:
            if proc_units.lower() == 'native':
                if native_units.lower() == "counts":
                    hthresh = hthresh * det_gain
                    if native_units.lower() == "counts/s":
                        hthresh = hthresh * img_exptime

        singleDriz = image.getOutputName("outSingle")
        singleDriz_name = image.outputNames['outSingle']
        singleWeight = image.getOutputName("outSWeight")
        singleWeight_name = image.outputNames['outSWeight']
        #singleDriz=image.outputNames["outSingle"] #all chips are drizzled to a single output image
        #singleWeight=image.outputNames["outSWeight"]

        # If compression was used, reference ext=1 as CompImageHDU only writes
        # out MEF files, not simple FITS.
        if compress:
            wcs_ext = '[1]'
            wcs_extnum = 1
        else:
            wcs_ext = '[0]'
            wcs_extnum = 0
        if not virtual:
            if isinstance(singleDriz,str):
                iter_singleDriz = singleDriz + wcs_ext
                iter_singleWeight = singleWeight + wcs_ext
            else:
                iter_singleDriz = singleDriz[wcs_extnum]
                iter_singleWeight = singleWeight[wcs_extnum]
        else:
            iter_singleDriz = singleDriz_name + wcs_ext
            iter_singleWeight = singleWeight_name + wcs_ext

        # read in WCS from first single drizzle image to use as WCS for median image
        if _single_hdr is None:
            if virtual:
                _single_hdr = singleDriz[wcs_extnum].header
            else:
                _single_hdr = fits.getheader(singleDriz_name,ext=wcs_extnum)

        _singleImage=iterfile.IterFitsFile(iter_singleDriz)
        if virtual:
            _singleImage.handle = singleDriz
            _singleImage.inmemory = True

        singleDrizList.append(_singleImage) #add to an array for bookkeeping

        # If it exists, extract the corresponding weight images
        if (not virtual and os.access(singleWeight,os.F_OK)) or (
                virtual and singleWeight):
            _weight_file=iterfile.IterFitsFile(iter_singleWeight)
            if virtual:
                _weight_file.handle = singleWeight
                _weight_file.inmemory = True

            singleWeightList.append(_weight_file)
            try:
                tmp_mean_value = ImageStats(_weight_file.data, lower=1e-8,
                    fields="mean", nclip=0).mean
            except ValueError:
                tmp_mean_value = 0.0
            _wht_mean.append(tmp_mean_value * maskpt)

            # Extract instrument specific parameters and place in lists

            # If an image has zero exposure time we will
            # redefine that value as '1'.  Although this will cause inaccurate scaling
            # of the data to occur in the 'minmed' combination algorith, this is a
            # necessary evil since it avoids divide by zero exceptions.  It is more
            # important that the divide by zero exceptions not cause Multidrizzle to
            # crash in the pipeline than it is to raise an exception for this obviously
            # bad data even though this is not the type of data you would wish to process
            # with Multidrizzle.
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
            nchips  = 0
            bsky    = None # minimum sky across **used** chips

            for chip in image.returnAllChips(extname=image.scienceExt):
                # compute sky value as sky/pixel using the single_drz pixel scale
                if bsky is None or bsky > chip.subtractedSky:
                    bsky = chip.subtractedSky * chip._conversionFactor

                # Extract the readnoise value for the chip
                rdnoise += (chip._rdnoise)**2
                nchips  += 1

            if bsky is None:
                bsky = 0.0

            if nchips > 0:
                rdnoise = math.sqrt(rdnoise/nchips)

            backgroundValueList.append(bsky)
            readnoiseList.append(rdnoise)

            ## compute sky value as sky/pixel using the single_drz pixel scale
            #bsky = image._image[image.scienceExt,1].subtractedSky# * (image.outputValues['scale']**2)
            #backgroundValueList.append(bsky)

            ## Extract the readnoise value for the chip
            #sci_chip = image._image[image.scienceExt,1]
            #readnoiseList.append(sci_chip._rdnoise) #verify this is calculated correctly in the image object

            print("reference sky value for image ",image._filename," is ", backgroundValueList[-1])
        #
        # END Loop over input image list
        #

    # create an array for the median output image, use the size of the first image in the list
    medianImageArray = np.zeros(singleDrizList[0].shape,dtype=singleDrizList[0].type())

    if ( comb_type.lower() == "minmed") and not newmasks:
        # Issue a warning if minmed is being run with newmasks turned off.
            print('\nWARNING: Creating median image without the application of bad pixel masks!\n')

    # create the master list to be used by the image iterator
    masterList = []
    masterList.extend(singleDrizList)
    masterList.extend(singleWeightList)

    print('\n')

    # Specify the location of the drz image sections
    startDrz = 0
    endDrz = len(singleDrizList)+startDrz

    # Specify the location of the wht image sections
    startWht = len(singleDrizList)+startDrz
    endWht = startWht + len(singleWeightList)
    _weight_mask_list = None

    # Fire up the image iterator
    #
    # The overlap value needs to be set to 2*grow in order to
    # avoid edge effects when scrolling down the image, and to
    # insure that the last section returned from the iterator
    # has enough row to span the kernel used in the boxcar method
    # within minmed.
    _overlap = 2*int(grow)

    #Start by computing the buffer size for the iterator
    _imgarr = masterList[0].data
    _bufsize = nimageiter.BUFSIZE
    if bufsizeMb is not None:
        _bufsize *= bufsizeMb
    _imgrows = _imgarr.shape[0]
    _nrows = nimageiter.computeBuffRows(_imgarr)
#        _overlaprows = _nrows - (_overlap+1)
#        _niter = int(_imgrows/_nrows)
#        _niter = 1 + int( (_imgrows - _overlaprows)/_nrows)
    niter = nimageiter.computeNumberBuff(_imgrows,_nrows,_overlap)
    #computeNumberBuff actually returns (niter,buffrows)
    _niter= niter[0]
    _nrows = niter[1]
    _lastrows = _imgrows - (_niter*(_nrows-_overlap))

    # check to see if this buffer size will leave enough rows for
    # the section returned on the last iteration
    if _lastrows < _overlap+1:
        _delta_rows = (_overlap+1 - _lastrows)//_niter
        if _delta_rows < 1 and _delta_rows >= 0: _delta_rows = 1
        _bufsize += (_imgarr.shape[1]*_imgarr.itemsize) * _delta_rows

    if not virtual:
        masterList[0].close()
    del _imgarr

    for imageSectionsList,prange in nimageiter.FileIter(masterList,overlap=_overlap,bufsize=_bufsize):

        if newmasks:
            """ Build new masks from single drizzled images. """
            _weight_mask_list = []
            listIndex = 0
            for _weight_arr in imageSectionsList[startWht:endWht]:
                # Initialize an output mask array to ones
                # This array will be reused for every output weight image
                _weight_mask = np.zeros(_weight_arr.shape,dtype=np.uint8)

                """ Generate new pixel mask file for median step.
                This mask will be created from the single-drizzled
                weight image for this image.

                The mean of the weight array will be computed and all
                pixels with values less than 0.7 of the mean will be flagged
                as bad in this mask.  This mask will then be used when
                creating the median image.
                """
                # Compute image statistics
                _mean = _wht_mean[listIndex]

                # 0 means good, 1 means bad here...
                np.putmask(_weight_mask, np.less(_weight_arr,_mean), 1)
                #_weight_mask.info()
                _weight_mask_list.append(_weight_mask)
                listIndex += 1

        # Do MINMED
        if ( "minmed" in comb_type.lower()):
            if comb_type.lower()[0] == 'i':
                # set up use of 'imedian'/'imean' in minmed algorithm
                fillval = True
            else:
                fillval = False

            if (_weight_mask_list in [None,[]]):
                _weight_mask_list = None

            # Create the combined array object using the minmed algorithm
            result = minmed(imageSectionsList[startDrz:endDrz],  # list of input data to be combined.
                                imageSectionsList[startWht:endWht],# list of input data weight images to be combined.
                                readnoiseList,                         # list of readnoise values to use for the input images.
                                exposureTimeList,                      # list of exposure times to use for the input images.
                                backgroundValueList,                   # list of image background values to use for the input images
                                weightMaskList = _weight_mask_list,  # list of imput data weight masks to use for pixel rejection.
                                combine_grow = grow,                   # Radius (pixels) for neighbor rejection
                                combine_nsigma1 = nsigma1,             # Significance for accepting minimum instead of median
                                combine_nsigma2 = nsigma2,              # Significance for accepting minimum instead of median
                                fillval=fillval                      # turn on use of imedian/imean
                            )
#              medianOutput[prange[0]:prange[1],:] = result.out_file1
#             minOutput[prange[0]:prange[1],:] = result.out_file2

        # DO NUMCOMBINE
        else:
            # Create the combined array object using the numcombine task
            result = numcombine.numCombine(imageSectionsList[startDrz:endDrz],
                                    numarrayMaskList=_weight_mask_list,
                                    combinationType=comb_type.lower(),
                                    nlow=nlow,
                                    nhigh=nhigh,
                                    upper=hthresh,
                                    lower=lthresh
                                )

        # We need to account for any specified overlap when writing out
        # the processed image sections to the final output array.
        if prange[1] != _imgrows:
            medianImageArray[prange[0]:prange[1]-_overlap,:] = result.combArrObj[:-_overlap,:]
        else:
            medianImageArray[prange[0]:prange[1],:] = result.combArrObj


    del result
    del _weight_mask_list
    _weight_mask_list = None

    # Write out the combined image
    # use the header from the first single drizzled image in the list
    #header=fits.getheader(imageObjectList[0].outputNames["outSingle"])
    _pf = _writeImage(medianImageArray, inputHeader=_single_hdr)

    if virtual:
        mediandict = {}
        mediandict[medianfile] = _pf
        for img in imageObjectList:
            img.saveVirtualOutputs(mediandict)
    else:
        try:
            print("Saving output median image to: ",medianfile)
            _pf.writeto(medianfile)
        except IOError:
            msg = "Problem writing file: "+medianfile
            print(msg)
            raise IOError(msg)

    del _pf

    # Always close any files opened to produce median image; namely,
    # single drizzle images and singly-drizzled weight images
    #

    for img in singleDrizList:
        if not virtual:
            img.close()
    singeDrizList = []

    # Close all singly drizzled weight images used to create median image.
    for img in singleWeightList:
        if not virtual:
            img.close()
    singleWeightList = []

    # If new median masks was turned on, close those files
    if _weight_mask_list:
        for arr in _weight_mask_list:
            del arr
        _weight_mask_list = None

    del masterList
    del medianImageArray

def _writeImage( dataArray=None, inputHeader=None):
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

    #_fname =inputFilename
    #_file = fits.open(_fname, mode='readonly')
    #_prihdu = fits.PrimaryHDU(header=_file[0].header,data=dataArray)

    _prihdu = fits.PrimaryHDU(data=dataArray, header=inputHeader)
    """
    if inputHeader is None:
        #use a general primary HDU
        _prihdu = fits.PrimaryHDU(data=dataArray)

    else:
        _prihdu = inputHeader
        _prihdu.data=dataArray
    """
    _pf = fits.HDUList()
    _pf.append(_prihdu)

    return _pf


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


median.__doc__ = getHelpAsString(docstring = True, show_ver = False)
