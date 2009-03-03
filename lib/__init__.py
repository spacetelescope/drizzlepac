""" BigBlackBox - test implementation of MultiDrizzle: The Next Generation

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

__taskname__ = "BigBlackBox"

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
def MultiDrizzle(input = '*flt.fits',output = None, shiftfile = None, updatewcs = True, editpars=True,
                configObj=None, wcsmap=wcs_functions.WCSMap, **input_dict):

    # Only create an updated input_dict if there is NO configObj provided to
    # avoid having the positional parameters values override those from the
    # configObj input.
    # Now, merge required input parameters into input_dict
    input_dict['input'] = input
    input_dict['output'] = output
    input_dict['shiftfile'] = shiftfile
    input_dict['updatewcs'] = updatewcs
    # If called from interactive user-interface, configObj will not be 
    # defined yet, so get defaults using EPAR/TEAL.
    #
    # Also insure that the input_dict (user-specified values) are folded in
    # with a fully populated configObj instance.
    configObj = util.getDefaultConfigObj(__taskname__,configObj,input_dict,loadOnly=(not editpars))
    
    if editpars == False:
        run(configObj,wcsmap=wcsmap)

#
#### Interfaces used by TEAL
#

def getHelpAsString():
    # Does NOT work with TEAL/cfgepar.epar()
        help_str = __doc__+'\n'
        help_str += 'Version '+__version__+'\n'
        return help_str
    
def run(configObj=None,wcsmap=wcs_functions.WCSMap):
    """    
    Initial example by Nadia ran MD with configObj EPAR using:
    It can be run in one of two ways:

        from pytools import cfgepar

        1. Passing a config object to epar

        from runmdz import mdriz
        mdobj = mdriz('multidrizzle/pars/mdriz.cfg')
        cfgepar.epar(mdobj)


        2. Passing a task  name:

        cfgepar.epar('multidrizzle')

        The example config files are in multidrizzle/pars

    """
    print '[BigBlackBox] mdriz is NOW running... \n'

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
    drizzle.drizSeparate(imgObjList,outwcs,configObj,wcsmap=wcs_functions.WCSMap)
    
    #create the median images from the driz sep images
    createMedian._median(imgObjList,configObj,configObj["clean"])
    
    #blot the images back to the original reference frame
    blot.runBlot(imgObjList, outwcs, configObj,wcsmap=wcs_functions.WCSMap)
    
    #look for cosmic rays
    drizCR.rundrizCR(imgObjList,configObj,saveFile=configObj["clean"])
    
    #Make your final drizzled image
    drizzle.drizFinal(imgObjList, outwcs, configObj,wcsmap=wcs_functions.WCSMap)
    
    print '\n[BigBlackBox] mdriz is all finished!\n'
    
    for image in imgObjList:
        image.close()
        
    del imgObjList
    del outwcs
    
    
    
    
    
