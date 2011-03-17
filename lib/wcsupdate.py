import os
import stwcs
from stwcs import updatewcs
from pytools import parseinput,fileutil
import pyfits

import convertwcs
import util

allowed_corr_dict = {'vacorr':'VACorr','tddcorr':'TDDCorr','npolcorr':'NPOLCorr','d2imcorr':'DET2IMCorr'}

__taskname__ = 'wcsupdate'
__version__ = stwcs.__version__
#
#### Interactive user interface (functional form)
#
def update(files, editpars=False, configObj=None, **input_dict):
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
        run(configObj)
        
#
#### Interfaces used by TEAL
#
def getHelpAsString():
    """ 
    return useful help from a file in the script directory called __taskname__.help
    """
    helpString = __taskname__+' Version '+__version__+'\n\n'

    # Start by using the docstring for the underlying task for this docstring
    helpString += updatewcs.updatewcs.__doc__

    return helpString

update.__doc__ = getHelpAsString()

def run(configObj=None):

    # Interpret primary parameters from configObj instance
    extname = configObj['extname']
    input = configObj['input']

    # create dictionary of remaining parameters, deleting extraneous ones
    # such as those above
    cdict = configObj.dict()
    # remove any rules defined for the TEAL interface
    if cdict.has_key("_RULES_"): del cdict['_RULES_']
    del cdict['_task_name_']
    del cdict['input']
    del cdict['extname']
    
    # parse input 
    input,altfiles = parseinput.parseinput(configObj['input'])

    # Insure that all input files have a correctly archived 
    #    set of OPUS WCS keywords
    # Legacy files from OTFR, like all WFPC2 data from OTFR, will only
    #   have the OPUS WCS keywords archived using a prefix of 'O'
    # These keywords need to be converted to the Paper I alternate WCS
    #   standard using a wcskey (suffix) of 'O'
    # If an alternate WCS with wcskey='O' already exists, this will copy
    #   the values from the old prefix-'O' WCS keywords to insure the correct
    #   OPUS keyword values get archived for use with updatewcs.
    #
    for file in input:
        # Check to insure that there is a valid reference file to be used
        idctab = pyfits.getval(file,'idctab')
        if not os.path.exists(fileutil.osfn(idctab)):
            print 'No valid distortion reference file ',idctab,' found in ',file,'!'
            raise ValueError

    # Re-define 'cdict' to only have switches for steps supported by that instrument
    # the set of supported steps are defined by the dictionary 
    #    updatewcs.apply_corrections.allowed_corrections
    #
    for file in input:
        # get instrument name from input file
        instr = pyfits.getval(file,'INSTRUME')
        # make copy of input parameters dict for this file
        fdict = cdict.copy()
        # Remove any parameter that is not part of this instrument's allowed corrections
        for step in allowed_corr_dict:
            if allowed_corr_dict[step] not in updatewcs.apply_corrections.allowed_corrections[instr]:
                fdict[step]
        # Call 'updatewcs' on correctly archived file
        updatewcs.updatewcs(file,**fdict)
    