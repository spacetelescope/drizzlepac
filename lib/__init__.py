""" betadrizzle - test implementation of MultiDrizzle: The Next Generation

MultiDrizzle automates the process of aligning images in an output frame, identifying cosmic-rays, removing distortion, and then combining the images while removing the identified cosmic-rays.  

This process involves a number of steps; namely:
  1.  Processing the input images and input parameters
  2.  Creating a static mask
  3.  Performing sky subtraction
  4.  Drizzling onto separate output images
  5.  Creating the median image
  6.  Blotting the median image
  7.  Identifying and flagging cosmic-rays
  8.  Final combination

A full description of this process can be found in the MultiDrizzle Handbook available online at:

http://stsdas.stsci.edu/multidrizzle

**Output**: The primary output from this task is the distortion-corrected, cosmic-ray cleaned combined image as a FITS file.

This task requires numerous user-settable parameters to control the primary aspects of each of the processing steps.  

"""
from __future__ import division # confidence high

import os,string

import buildmask
import drizzle
import blot
import imageObject
import outputimage
import processInput,mdzhandler
import sky
import createMedian
import drizCR
import staticMask
import util
import wcs_functions
import resetbits
import wcscorr
import pixtosky

# The following modules are for 'tweakreg' and are included here to make
# it easier to get to this code interactively
import tweakreg, catalogs, imgclasses, tweakutils, wcscorr
import imagefindpars, sextractorpars

# Add updatenpol to the list of tasks imported automatically here
import updatenpol

# Add TEAL interface to 'updatewcs' here
import wcsupdate

__taskname__ = "betadrizzle"

# These lines allow TEAL to print out the names of TEAL-enabled tasks 
# upon importing this package.
from pytools import teal
teal.print_tasknames(__name__, os.path.dirname(__file__))

# Begin Version Information -------------------------------------------
# Revision based version info
try:
    import svn_version
    __svn_version__ = svn_version.__svn_version__
except:
    __svn_version__ = 'Unable to determine SVN revision'

__version__ = '4.0.9dev11906'
__vdate__ = "10-Feb-2011"
# End Version Information ---------------------------------------------

# Pointer to the included Python class for WCS-based coordinate transformations
PYTHON_WCSMAP = wcs_functions.WCSMap

#
#### Interactive user interface (functional form)
#
def MultiDrizzle(files, editpars=False, configObj=None, wcsmap=None, **input_dict):
    """
    """
    # support input of filenames from command-line without a parameter name
    # then copy this into input_dict for merging with TEAL ConfigObj parameters
    if files not in ['',' ','INDEF', None]:
        if input_dict is None:
            input_dict = {}
        input_dict['input'] = files
    
    # If called from interactive user-interface, configObj will not be 
    # defined yet, so get defaults using EPAR/TEAL.
    #
    # Also insure that the input_dict (user-specified values) are folded in
    # with a fully populated configObj instance.
    configObj = util.getDefaultConfigObj(__taskname__,configObj,input_dict,loadOnly=(not editpars))
    if configObj is None:
        return
    # If 'editpars' was set to True, util.getDefaultConfigObj() will have already
    # called 'run()'.
    if editpars == False:
        run(configObj,wcsmap=wcsmap)
        
#
#### Interfaces used by TEAL
#

def getHelpAsString():
    """ 
    return useful help from a file in the script directory called __taskname__.help
    """
    helpString = __taskname__+' Version '+__version__+'\n\n'
    helpString += teal.getHelpFileAsString(__taskname__,__file__)
    
    return helpString

MultiDrizzle.__doc__ = getHelpAsString()

def run(configObj=None,wcsmap=None):
    """    
    Initial example by Nadia ran MD with configObj EPAR using:
    It can be run in one of two ways:

        from pytools import teal

        1. Passing a config object to teal

        from runmdz import mdriz
        mdobj = mdriz('multidrizzle/pars/mdriz.cfg')
        teal.teal(mdobj)


        2. Passing a task  name:

        teal.teal('multidrizzle')

        The example config files are in multidrizzle/pars

    """
    #
    # turn on logging, redirecting stdout/stderr messages to a log file
    # while also printing them out to stdout as well
    # also, initialize timing of processing steps
    # 
    util.init_logging(logfile=configObj['runfile'])
    procSteps = util.ProcSteps()
    print '[betadrizzle] MultiDrizzle Version '+__version__+' started at: ',util._ptime()[0],'\n'
    try:
        try:
            # Define list of imageObject instances and output WCSObject instance
            # based on input paramters
            procSteps.addStep('Initialization')
            imgObjList = None
            imgObjList,outwcs = processInput.setCommonInput(configObj)
            procSteps.endStep('Initialization')
            if not imgObjList:
                return
            
            print "\nUSER INPUT PARAMETERS common to all Processing Steps:"
            util.printParams(configObj)
                      
            # Call rest of MD steps...
            #create static masks for each image
            staticMask.createStaticMask(imgObjList,configObj,procSteps=procSteps)
            
            #subtract the sky
            sky.subtractSky(imgObjList,configObj,procSteps=procSteps)
            
            #drizzle to separate images
            drizzle.drizSeparate(imgObjList,outwcs,configObj,wcsmap=wcsmap,procSteps=procSteps)
            
            #create the median images from the driz sep images
            createMedian.createMedian(imgObjList,configObj,procSteps=procSteps)
            
            #blot the images back to the original reference frame
            blot.runBlot(imgObjList, outwcs, configObj,wcsmap=wcsmap,procSteps=procSteps)
            
            #look for cosmic rays
            drizCR.rundrizCR(imgObjList,configObj,saveFile=configObj["clean"],procSteps=procSteps)
            
            #Make your final drizzled image
            drizzle.drizFinal(imgObjList, outwcs, configObj,wcsmap=wcsmap,procSteps=procSteps)
            
            print '\n[betadrizzle]MultiDrizzle Version '+__version__+' is all finished at ',util._ptime()[0],' !\n'

        except:  
            raise
    finally:
        procSteps.reportTimes()

        if imgObjList:
            for image in imgObjList:
                if configObj['clean']:
                    image.clean()
                image.close()
                
            del imgObjList
            del outwcs

        # Turn off logging now
        util.end_logging()
        
