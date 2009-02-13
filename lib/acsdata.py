#!/usr/bin/env python

"""
These are functions related directly to ACS images
that grab instrument specific informationrigh
"""

from pytools import fileutil
import numpy as np


def getACSInfo(imageSet=None):
    """
    This takes a pyfits primary header instance and pulls out
    the basic set of keywords we would like to have which 
    are specific to ACS and then returns a dictionary of key value pairs

    The keyword list is not passed in because they
    can be different from instrument to instrument
    """
    keywords={} #store keyword name specifics
    finalDict={} #store the final keyword:values to be returned, 
                #where they keywords be the same across instruments
    
    #the flatfile keyword name           
    keywords["flatname"]="PFLTFILE"
    keywords["cr_bits_value"] = 4096

    # Effective gain to be used in the driz_cr step.  Since the
    # ACS images have already been converted to electrons,
    # the effective gain is 1.    
    keywords["effGain"]=1.

    try:    
        keywords["darkcurrent"]=primaryHeader["MEANDARK"] #in header units        
    except:                                                         
        str =  "#############################################\n"        
        str += "#                                           #\n"        
        str += "# Error:                                    #\n"        
        str += "#   Cannot find the value for 'MEANDARK'    #\n"        
        str += "#   in the image header.  ACS input images  #\n"        
        str += "#   are expected to have this header        #\n"        
        str += "#   keyword.                                #\n"        
        str += "#                                           #\n"        
        str += "# Error occured in getASCInfo               #\n"        
        str += "#                                           #\n"        
        str += "#############################################\n"        
        raise ValueError, str                                           
    
    
    keywords["GAIN"]=float(imageSet._image["PRIMARY"].header["ATODGNA"])
    keywords["READNOISE"]=float(imageSet._image["PRIMARY"].header["READNSEA"])
    
    keywords["PLATESCALE"]=0.04
    
    
    finalDict["flatfield"]=imageSet._image["PRIMARY"].header[keywords["flatname"]]
    #there could be different gains for each chip
    for chip in range(1,imageSet._numchips+1,1):
        finalDict["gain"]=[imageSet._image[
    
    return keywords

def getCTEdirection(extension):
    """given the science extension, get the direction of the CTE
    """
    if ( extension == 'sci,1') : # get cte direction, which depends on which chip but is independent of amp 
        cte_dir = -1    
    if ( extension == 'sci,2') : 
        cte_dir = 1   
    
    return cte_dir
    

def getHRCInfo(finalDict={}):

        instrument = 'ACS/HRC'        
        full_shape = (1024,1024)
        platescale = platescale

        if ( amp == 'A' or amp == 'B' ) : # cte direction depends on amp (but is independent of chip)
            cte_dir = 1   
        if ( amp == 'C' or amp == 'D' ) :
            cte_dir = -1   

def getSBCInfo (finalDict={}):

    full_shape = (1024,1024)
    platescale = platescale
    instrument = 'ACS/SBC'

    # no cte correction for SBC so set cte_dir=0.
    print('\nWARNING: No cte correction will be made for this SBC data.\n')
    cte_dir = 0       

    # Set the default readnoise or gain values based upon the amount of user input given.

    # Case 1: User supplied no gain or readnoise information
    if usingDefaultReadnoise and usingDefaultGain:
        # Set the default gain and readnoise values
        _setSBCchippars()
    # Case 2: The user has supplied a value for gain
    elif usingDefaultReadnoise and not usingDefaultGain:
        # Set the default readnoise value
        _setDefaultSBCReadnoise()
    # Case 3: The user has supplied a value for readnoise 
    elif not usingDefaultReadnoise and usingDefaultGain:
        # Set the default gain value
        _setDefaultSBCGain()
    else:
        # In this case, the user has specified both a gain and readnoise values.  Just use them as is.
        pass

def _setSBCchippars(finalDict={}):
    _setDefaultSBCGain()
    _setDefaultSBCReadnoise()

def _setDefaultSBCGain(finalDict{}):
    _gain = 1

def _setDefaultSBCReadnoise(finalDict{}):
    _rdnoise = 0



