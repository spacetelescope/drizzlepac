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
__version__ = '4.0.0dev'
# End Version Information ---------------------------------------------
# Revision based version info
try:
    import svn_version
    __svn_version__ = svn_version.__svn_version__
except:
    __svn_version__ = 'Unable to determine SVN revision'
    
#
#### Interactive user interface (functional form)
#
def MultiDrizzle(editpars=False, configObj=None, wcsmap=wcs_functions.WCSMap, **input_dict):

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
        help_str = __doc__+'\n'
        help_str += 'Version '+__version__+'\n'
        return help_str
    
def run(configObj=None,wcsmap=wcs_functions.WCSMap):
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
    print '[betadrizzle] mdriz is NOW running... \n'

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
    
    print '\n[betadrizzle] mdriz is all finished!\n'
    
    for image in imgObjList:
        image.close()
        
    del imgObjList
    del outwcs
    
    
    
    
    
