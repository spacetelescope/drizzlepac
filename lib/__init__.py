import buildmask
import drizzle
import imageObject
import MultiDrizzle
import outputimage
import processInput
import sky
import staticMask
import util
import wcs_functions

__cfg_file__ = "pars/mdriz.cfg"

# Begin Version Information -------------------------------------------
__version__ = '3.3.0dev'
# End Version Information ---------------------------------------------
# Revision based version info
try:
    import svn_version
    __svn_version__ = svn_version.__svn_version__
except:
    __svn_version__ = 'Unable to determine SVN revision'

# Interfaces used by TEAL
def run(configObj):
    MultiDrizzle(configObj=configObj)

# Interactive user interface (functional form)
def Multidrizzle(input = '*flt.fits',output = None, shiftfile = None, updatewcs = True, 
                configObj=None, **input_dict):

    # Only create an updated input_dict if there is NO configObj provided to
    # avoid having the positional parameters values override those from the
    # configObj input.
    if configObj is None:        
        # Now, merge required input parameters into input_dict
        input_dict['input'] = input
        input_dict['output'] = output
        input_dict['shiftfile'] = shiftfile
        input_dict['updatewcs'] = updatewcs
    
    imgObjList,outwcs = processInput.processCommonInput(input_dict,configObj,cfg_file=__cfg_file__)
    
    # Call rest of MD steps...