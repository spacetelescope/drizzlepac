"""
This module provides functions and classes that manage the creation
of the global static masks.

For `staticMask`, the user interface function is :py:func:`createMask`.

:Authors: Ivo Busko, Christopher Hanley, Warren Hack, Megan Sosey

:License: :doc:`LICENSE`

"""
import os
import sys
from distutils.version import LooseVersion

import numpy as np
from stsci.tools import fileutil, teal, logutil
import astropy
from astropy.io import fits
from stsci.imagestats import ImageStats
from . import util
from . import processInput

__taskname__ = "drizzlepac.staticMask"
_step_num_ = 1


log = logutil.create_logger(__name__, level=logutil.logging.NOTSET)


#this is called by the user
def createMask(input=None, static_sig=4.0, group=None, editpars=False, configObj=None, **inputDict):
    """ The user can input a list of images if they like to create static masks
        as well as optional values for static_sig and inputDict.

        The configObj.cfg file will set the defaults and then override them
        with the user options.
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

#this is called by the TEAL interface
def run(configObj):

    #now we really just need the imageObject list created for the dataset
    filelist,output,ivmlist,oldasndict=processInput.processFilenames(configObj['input'],None)

    imageObjList=processInput.createImageObjectList(filelist,instrpars={},group=configObj['group'])
    createStaticMask(imageObjList,configObj)


#this is the workhorse function called by MultiDrizzle
def createStaticMask(imageObjectList=[],configObj=None,procSteps=None):
    if procSteps is not None:
        procSteps.addStep('Static Mask')

    step_name = util.getSectionName(configObj,_step_num_)

    if not configObj[step_name]['static']:
        log.info('Static Mask step not performed.')
        procSteps.endStep('Static Mask')
        return

    if (not isinstance(imageObjectList,list) or (len(imageObjectList) ==0)):
        msg = "Invalid image object list given to static mask"
        print(msg, file=sys.stderr)
        raise ValueError(msg)

    log.info('USER INPUT PARAMETERS for Static Mask Step:')
    util.printParams(configObj[step_name], log=log)

    #create a static mask object
    myMask = staticMask(configObj)

    for image in imageObjectList:
        myMask.addMember(image) # create tmp filename here...


    #save the masks to disk for later access
    myMask.saveToFile(imageObjectList)
    myMask.close()

    if procSteps is not None:
        procSteps.endStep('Static Mask')

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
        self.step_name=util.getSectionName(configObj,_step_num_)
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

        #for chip in range(1,numchips+1,1):
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

            stats = ImageStats(chipimage,nclip=3,fields='mode')
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
            #check to see if the file already exists on disk
            filename = self.masknames[key]
            #create a new fits image with the mask array and a standard header
            #open a new header and data unit
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


createMask.__doc__ = getHelpAsString(docstring = True, show_ver = False)
