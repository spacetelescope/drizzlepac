""" betadrizzle - test implementation of MultiDrizzle: The Next Generation

"""
import os

import buildmask
import drizzle
import blot
import imageObject
import MultiDrizzleObject
import outputimage
import processInput,mdzhandler
import sky
import createMedian
import drizCR
import staticMask
import util
import wcs_functions

__taskname__ = "betadrizzle"

# Begin Version Information -------------------------------------------
# Revision based version info
try:
    import svn_version
    __svn_version__ = svn_version.__svn_version__
except:
    __svn_version__ = 'Unable to determine SVN revision'

__version__ = '4.0.1dev8577'
# End Version Information ---------------------------------------------

# Pointer to the included Python class for WCS-based coordinate transformations
PYTHON_WCSMAP = wcs_functions.WCSMap

#
#### Interactive user interface (functional form)
#
def MultiDrizzle(editpars=False, configObj=None, wcsmap=None, **input_dict):

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
    # Does NOT work with TEAL/teal.teal()
    helpString = __doc__+'\n'
    helpString += 'Version '+__version__+'\n'

    """ 
    return useful help from a file in the script directory called module.help
    """
    #get the local library directory where the code is stored
    localDir=os.path.split(__file__)
    helpfile=__taskname__.split(".")
    
    helpfile=localDir[0]+"/"+helpfile[0]+".help"
    
    if os.access(helpfile,os.R_OK):
        fh=open(helpfile,'r')
        ss=fh.readlines()
        fh.close()
        #helpString=""
        for line in ss:
            helpString+=line
    else:    
        helpString=__doc__

    return helpString

    
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
    print '[betadrizzle] mdriz started at: ',util._ptime(),'\n'
    try:
        # Define list of imageObject instances and output WCSObject instance
        # based on input paramters
        imgObjList,outwcs = processInput.setCommonInput(configObj)

        # Call rest of MD steps...
        print 'Finished interpreting configObj...\n'
        #create static masks for each image
        staticMask._staticMask(imgObjList,configObj)
        
        #subtract the sky
        sky.subtractSky(imgObjList,configObj)
        
        #drizzle to separate images
        drizzle.drizSeparate(imgObjList,outwcs,configObj,wcsmap=wcsmap)
        
        #create the median images from the driz sep images
        createMedian._median(imgObjList,configObj,configObj["clean"])
        
        #blot the images back to the original reference frame
        blot.runBlot(imgObjList, outwcs, configObj,wcsmap=wcsmap)
        
        #look for cosmic rays
        drizCR.rundrizCR(imgObjList,configObj,saveFile=configObj["clean"])
        
        #Make your final drizzled image
        drizzle.drizFinal(imgObjList, outwcs, configObj,wcsmap=wcsmap)
        
        print '\n[betadrizzle] mdriz is all finished at ',util._ptime(),' !\n'

    finally:
        for image in imgObjList:
            if configObj['clean']:
                image.clean()
            image.close()
            
        del imgObjList
        del outwcs
        