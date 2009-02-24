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
import staticMask
import util
import wcs_functions
from pytools import cfgpars

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
def MultiDrizzle(input = '*flt.fits',output = None, shiftfile = None, updatewcs = True, 
                configObj=None, **input_dict):

    # Only create an updated input_dict if there is NO configObj provided to
    # avoid having the positional parameters values override those from the
    # configObj input.
    # Now, merge required input parameters into input_dict
    input_dict['input'] = input
    input_dict['output'] = output
    input_dict['shiftfile'] = shiftfile
    input_dict['updatewcs'] = updatewcs
    
    run(configObj,input_dict=input_dict)

#
#### Interfaces used by TEAL
#
class mdriz(cfgpars.ConfigObjPars):
    """ This needs to be called using the following syntax:

        mdobj = BigBlackBox.mdriz()
        cfgepar.epar(mdobj)

    """
    def __init__(self, cfgFileName):
        if cfgFileName is None:
            cfgFileName = __cfg_file__
        cfgpars.ConfigObjPars.__init__(self, cfgFileName)

    def run(self, *args, **kw):
        # Place your code to invoke Multidrizzle here
        print "running MultiDrizzle from TEAL..."        
        MultiDrizzle(configObj=self)
    def getHelpAsString(self):
        getHelpAsString()
            
def getHelpAsString():
    # Does NOT work with TEAL/cfgepar.epar()
        help_str = __doc__+'\n'
        help_str += 'Version '+__version__+'\n'
        return help_str
    
def run(configObj=None,input_dict={},wcsmap=wcs_functions.WCSMap,loadOnly=False):
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
    # If called from interactive user-interface, configObj will not be 
    # defined yet, so get defaults using EPAR/TEAL.
    #
    # Also insure that the input_dict (user-specified values) are folded in
    # with a fully populated configObj instance.
    configObj = util.getDefaultConfigObj(__taskname__,configObj,input_dict,loadOnly=loadOnly)
    print '[BigBlackBox] mdriz is NOW running... '

    # Define list of imageObject instances and output WCSObject instance
    # based on input paramters
    imgObjList,outwcs = processInput.setCommonInput(configObj)
    
    # Call rest of MD steps...
    print 'Finished interpreting configObj...'
