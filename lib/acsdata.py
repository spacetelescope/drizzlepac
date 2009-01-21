#!/usr/bin/env python

from pytools import fileutil
import numpy as np

"""
These are functions related directly to ACS images
"""

def getACSInfo(primaryHeader):
    """
    This takes a pyfits primary header instance and pulls out
    the basic set of keywords we would like to have which 
    are specific to ACS and then returns a dictionary of key value pairs

    The keyword list is not passed in because they
    can be different from instrument to instrument
    """
    keywords={}
    
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

    #the number of science chips in the image
    #ACS has SCI,ERR,DQ for each science chip
    #so,  nextend / 3 gives you the number of science chips
    #can't think of a better way to get this right now
    keywords["NUMCHIPS"]=int(primaryHeader["NEXTEND"]) / 3
    keywords["GAIN"]=float(primaryHeader["ATODGNA"])
    keywords["READNOISE"]=float(primaryHeader["READNSEA"])
    
    keywords["PLATESCALE"]=0.04
    
    return keywords

def getCTEdirection(extension):
    """given the science extension, get the direction of the CTE
    """
    if ( extension == 'sci,1') : # get cte direction, which depends on which chip but is independent of amp 
        cte_dir = -1    
    if ( extension == 'sci,2') : 
        cte_dir = 1   
    
    return cte_dir
    
"""
def getflat(self):

    Purpose
    =======
    Method for retrieving a detector's flat field.

    This method will return an array the same shape as the
    image.

    :units: electrons


    # The keyword for ACS flat fields in the primary header of the flt
    # file is pfltfile.  This flat file is already in the required 
    # units of electrons.

    filename = self.header['PFLTFILE']

    try:
        handle = fileutil.openImage(filename,mode='readonly',memmap=0)
        hdu = fileutil.getExtn(handle,extn=self.extn)
        data = hdu.data[self.ltv2:self.size2,self.ltv1:self.size1]
    except:
        try:
            handle = fileutil.openImage(filename[5:],mode='readonly',memmap=0)
            hdu = fileutil.getExtn(handle,extn=self.extn)
            data = hdu.data[self.ltv2:self.size2,self.ltv1:self.size1]
        except:
            data = np.ones(self.image_shape,dtype=self.image_dtype)
            str = "Cannot find file "+filename+".  Treating flatfield constant value of '1'.\n"
            print str
    flat = data
    return flat

class HRCInputImage (ACSInputImage):

    def __init__(self, input, dqname, platescale, memmap=0,proc_unit="native"):
        ACSInputImage.__init__(self, input, dqname, platescale,memmap=0,proc_unit=proc_unit)
        self.instrument = 'ACS/HRC'        
        self.full_shape = (1024,1024)
        self.platescale = platescale

        if ( self.amp == 'A' or self.amp == 'B' ) : # cte direction depends on amp (but is independent of chip)
            self.cte_dir = 1   
        if ( self.amp == 'C' or self.amp == 'D' ) :
            self.cte_dir = -1   

class SBCInputImage (ACSInputImage):

    def __init__(self, input, dqname, platescale, memmap=0,proc_unit="native"):
        ACSInputImage.__init__(self,input,dqname,platescale,memmap=0,proc_unit=proc_unit)
        self.full_shape = (1024,1024)
        self.platescale = platescale
        self.instrument = 'ACS/SBC'

        # no cte correction for SBC so set cte_dir=0.
        print('\nWARNING: No cte correction will be made for this SBC data.\n')
        self.cte_dir = 0       

        # Set the default readnoise or gain values based upon the amount of user input given.

    # Case 1: User supplied no gain or readnoise information
    if usingDefaultReadnoise and usingDefaultGain:
        # Set the default gain and readnoise values
        self._setSBCchippars()
    # Case 2: The user has supplied a value for gain
    elif usingDefaultReadnoise and not usingDefaultGain:
        # Set the default readnoise value
        self._setDefaultSBCReadnoise()
    # Case 3: The user has supplied a value for readnoise 
    elif not usingDefaultReadnoise and usingDefaultGain:
        # Set the default gain value
        self._setDefaultSBCGain()
    else:
        # In this case, the user has specified both a gain and readnoise values.  Just use them as is.
        pass

def _setSBCchippars(self):
    self._setDefaultSBCGain()
    self._setDefaultSBCReadnoise()

def _setDefaultSBCGain(self):
    self._gain = 1

def _setDefaultSBCReadnoise(self):
    self._rdnoise = 0
"""
