"""
This module provides functions and classes that manage the creation
of the global static masks.

For ``staticMask``, the user interface function is :py:func:`createMask`.

:Authors: Ivo Busko, Christopher Hanley, Warren Hack, Megan Sosey

:License: :doc:`/LICENSE`

"""
import os
import sys

import numpy as np
from stsci.tools import fileutil, logutil
from astropy.io import fits
from stsci.imagestats import ImageStats
from . import util
from . import processInput

__taskname__ = "staticMask"
STEP_NUM = 1
PROCSTEPS_NAME = "Static Mask"

log = logutil.create_logger(__name__, level=logutil.logging.NOTSET)


# this is called by the user
def createMask(input=None, static_sig=4.0, group=None, editpars=False, configObj=None, **inputDict):
    """ The user can input a list of images if they like to create static masks
        as well as optional values for static_sig and inputDict.

        The configObj.cfg file will set the defaults and then override them
        with the user options.

        Create a static mask for all input images. The mask contains pixels that fall
        more than ``static_sig`` RMS below the mode for a given chip or extension.
        Those severely negative, or low pixels, might result from oversubtraction
        of bad pixels in the dark image, or high sky levels during calibration.
        For example, each ACS WFC image contains a separate image for each of 2 CCDs,
        and seperate masks will be generated for each chip accordingly.

        The final static mask for each chip contains all of the bad pixels that meet
        this criteria from all of the input images along with any bad pixels that
        satisfy the final_bits value specified by the user, and found in the images
        DQ mask.

        Users should consider the details of their science image and decide whether
        or not creating this mask is appropriate for their resulting science.
        For example, if your field is very crowded, or contains mostly nebulous
        or extended objects, then the statistcs could be heavily skewed and the mask
        could end up containing sources.

        The generated static masks are saved to disk for use in later steps with
        the following naming convention:

            [Instrument][Detector]_[xsize]x[ysize]_[detector number]_staticMask.fits

        so an ACS image would produce a static mask with the name:

            ACSWFC_2048x4096_1_staticMask.fits

        and this would be the only file saved to disk, storing the logic and of all
        the badpixel masks created for each acs image in the set.

        For more information on the science applications of the static mask task,
        see the `DrizzlePac Handbook <http://drizzlepac.stsci.edu>`_


        Parameters
        ----------
        input : str, None (Default = None)
            A list of images or associations you would like to use to compute
            the mask.

        static : bool (Default = True)
            Create a static bad-pixel mask from the data?  This mask flags all pixels
            that deviate by more than a value of ``static_sig`` sigma below the image
            median, since these pixels are typically the result of bad pixel
            oversubtraction in the dark image during calibration.

        static_sig : float (Default = 4.0)
            The number of sigma below the RMS to use as the clipping limit for
            creating the static mask.

        editpars : bool (Default = False)
            Set to `True` if you would like to edit the parameters using the GUI
            interface.


        Examples
        --------
        These tasks are designed to work together seemlessly when run in the full
        ``AstroDrizzle`` interface. More advanced users may wish to create specialized
        scripts for their own datasets, making use of only a subset of the predefined
        ``AstroDrizzle`` tasks, or add additional processing, which may be usefull for
        their particular data. In these cases, individual access to the tasks is
        important.

        Something to keep in mind is that the full ``AstroDrizzle`` interface will
        make backup copies of your original files and place them in
        the ``OrIg/`` directory of your current working directory. If you are working
        with the stand alone interfaces, it is assumed that the user has already
        taken care of backing up their original datafiles as the input file with
        be directly altered.

        Basic example of how to call static yourself from a python command line,
        using the default parameters for the task.

        >>> from drizzlepac import staticMask
        >>> staticMask.createMask('*flt.fits')
    """

    if input is not None:
        inputDict["static_sig"]=static_sig
        inputDict["group"]=group
        inputDict["updatewcs"]=False
        inputDict["input"]=input
    else:
        print >> sys.stderr, "Please supply an input image\n"
        raise ValueError

    #this accounts for a user-called init where config is not defined yet
    configObj = util.getDefaultConfigObj(__taskname__,configObj,inputDict,loadOnly=(not editpars))
    if configObj is None:
        return

    if not editpars:
        run(configObj)

#--------------------------------
# TEAL Interface functions
# (these functions are deprecated)
#---------------------------------
def run(configObj):

    #now we really just need the imageObject list created for the dataset
    filelist,output,ivmlist,oldasndict=processInput.processFilenames(configObj['input'],None)

    imageObjList=processInput.createImageObjectList(filelist,instrpars={},group=configObj['group'])
    createStaticMask(imageObjList,configObj)


# this is the workhorse function called by MultiDrizzle
def createStaticMask(imageObjectList=[],configObj=None,procSteps=None):
    if procSteps is not None:
        procSteps.addStep(PROCSTEPS_NAME)

    step_name = util.getSectionName(configObj,STEP_NUM)

    if not configObj[step_name]['static']:
        log.info(f"{PROCSTEPS_NAME} step not performed.")
        procSteps.endStep(PROCSTEPS_NAME)
        return

    if not isinstance(imageObjectList, list) or len(imageObjectList) == 0:
        procSteps.skipStep(PROCSTEPS_NAME, reason="aborted")
        msg = "Invalid image object list given to static mask"
        print(msg, file=sys.stderr)
        raise ValueError(msg)

    log.info(f"USER INPUT PARAMETERS for {PROCSTEPS_NAME} Step:")
    util.printParams(configObj[step_name], log=log)

    #create a static mask object
    myMask = staticMask(configObj)

    for image in imageObjectList:
        myMask.addMember(image) # create tmp filename here...

    #save the masks to disk for later access
    myMask.saveToFile(imageObjectList)
    myMask.close()

    if procSteps is not None:
        procSteps.endStep(PROCSTEPS_NAME)

def constructFilename(signature):
    """Construct an output filename for the given signature::

         signature=[instr+detector,(nx,ny),detnum]

    The signature is in the image object.
    """
    suffix = buildSignatureKey(signature)
    filename = os.path.join('.', suffix)
    return filename

def buildSignatureKey(signature):
    """
    Build static file filename suffix used by mkstemp()
    """
    return signature[0]+"_"+str(signature[1][0])+"x"+str(signature[1][1])+"_"+str(signature[2])+"_staticMask.fits"

class staticMask:
    """
    This class manages the creation of the global static mask which
    masks pixels that are unwanted in the SCI array.
    A static mask  object gets created for each global
    mask needed, one for each chip from each instrument/detector.
    Each static mask array has type Int16, and resides in memory.

    :Notes:
        Class that manages the creation of a global static
        mask which is used to mask pixels that are some
        sigma BELOW the mode computed for the image.

    """
    def __init__(self, configObj=None):

        # the signature is created in the imageObject class

        self.masklist={}
        self.masknames = {}
        self.step_name=util.getSectionName(configObj, STEP_NUM)
        if configObj is not None:
            self.static_sig = configObj[self.step_name]['static_sig']
        else:
            self.static_sig = 4. # define a reasonable number
            log.warning('Using default of 4. for static mask sigma.')

    def addMember(self, imagePtr=None):
        """
        Combines the input image with the static mask that
        has the same signature.

        Parameters
        ----------
        imagePtr : object
            An imageObject reference

        Notes
        -----
        The signature parameter consists of the tuple::

            (instrument/detector, (nx,ny), chip_id)

        The signature is defined in the image object for each chip

        """
        numchips=imagePtr._numchips
        log.info("Computing static mask:\n")

        chips = imagePtr.group
        if chips is None:
            chips = imagePtr.getExtensions()

        # for chip in range(1,numchips+1,1):
        for chip in chips:
            chipid=imagePtr.scienceExt + ','+ str(chip)
            chipimage=imagePtr.getData(chipid)
            signature=imagePtr[chipid].signature

            # If this is a new signature, create a new Static Mask file which is empty
            # only create a new mask if one doesn't already exist
            if ((signature not in self.masklist) or (len(self.masklist) == 0)):
                self.masklist[signature] = self._buildMaskArray(signature)
                maskname =  constructFilename(signature)
                self.masknames[signature] = maskname
            else:
                chip_sig = buildSignatureKey(signature)
                for s in self.masknames:
                    if chip_sig in self.masknames[s]:
                        maskname  = self.masknames[s]
                        break
            imagePtr[chipid].outputNames['staticMask'] = maskname
            stats = ImageStats(
                chipimage,
                nclip=3,
                fields="mode",
                lower=np.nanmin(chipimage),
                upper=np.nanmax(chipimage),
            )
            mode = stats.mode
            rms  = stats.stddev
            nbins = len(stats.histogram)
            del stats

            log.info('  mode = %9f;   rms = %7f;   static_sig = %0.2f' %
                     (mode, rms, self.static_sig))

            if nbins >= 2: # only combine data from new image if enough data to mask
                sky_rms_diff = mode - (self.static_sig*rms)
                np.bitwise_and(self.masklist[signature],
                               np.logical_not(np.less(chipimage, sky_rms_diff)),
                               self.masklist[signature])
            del chipimage

    def _buildMaskArray(self,signature):
        """ Creates empty  numpy array for static mask array signature. """
        return np.ones(signature[1],dtype=np.int16)

    def getMaskArray(self, signature):
        """ Returns the appropriate StaticMask array for the image. """
        if signature in self.masklist:
            mask =  self.masklist[signature]
        else:
            mask = None
        return mask

    def getFilename(self,signature):
        """Returns the name of the output mask file that
        should reside on disk for the given signature. """

        filename=constructFilename(signature)

        if(fileutil.checkFileExists(filename)):
            return filename
        else:
            print("\nmMask file for ", str(signature), " does not exist on disk", file=sys.stderr)
            return None

    def getMaskname(self,chipid):
        """Construct an output filename for the given signature::

             signature=[instr+detector,(nx,ny),detnum]

        The signature is in the image object and the
        name of the static mask file is saved as sci_chip.outputNames["staticMask"].
        """
        return self._image[chipid].outputNames["staticMask"]

    def close(self):
        """ Deletes all static mask objects. """

        for key in self.masklist.keys():
            self.masklist[key] = None
        self.masklist = {}

    def deleteMask(self,signature):
        """ Delete just the mask that matches the signature given."""
        if signature in self.masklist:
            self.masklist[signature] = None
        else:
            log.warning("No matching mask")

    def saveToFile(self,imageObjectList):
        """ Saves the static mask to a file
            it uses the signatures associated with each
            mask to contruct the filename for the output mask image.
        """
        virtual = imageObjectList[0].inmemory

        for key in self.masklist.keys():
            # check to see if the file already exists on disk
            filename = self.masknames[key]
            # create a new fits image with the mask array and a standard header
            # open a new header and data unit
            newHDU = fits.PrimaryHDU()
            newHDU.data = self.masklist[key]

            if virtual:
                for img in imageObjectList:
                    img.saveVirtualOutputs({filename:newHDU})

            else:
                try:
                    newHDU.writeto(filename, overwrite=True)
                    log.info("Saving static mask to disk: %s" % filename)

                except IOError:
                    log.error("Problem saving static mask file: %s to "
                              "disk!\n" % filename)
                    raise IOError


createMask.__doc__ = util._def_help_functions(
    locals(), module_file=__file__, task_name=__taskname__, module_doc=__doc__
)
