#!/usr/bin/env python

from pytools import fileutil
import numpy as np
from input_image import InputImage


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
    
    #the flatfile keyword name           
    keywords["flatname"]="PFLTFILE"
    cr_bits_value = 4096

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
        
    
    
class ACSdata(filename):

    SEPARATOR = '_'

    def __init__(self, input,dqname,platescale,memmap=0,proc_unit="native"):
        InputImage.__init__(self,input,dqname,platescale,memmap=0,proc_unit=proc_unit)
        # define the cosmic ray bits value to use in the dq array
        self.cr_bits_value = 4096
        self.platescale = platescale
        self._effGain = 1

    def doUnitConversions(self):
        # Effective gain to be used in the driz_cr step.  Since the
        # ACS images have already been converted to electrons,
        # the effective gain is 1.
        self._effGain = 1
        
    def _isSubArray(self):
        _subarray = False
        _ltv1 = float(fileutil.getKeyword(parlist['data'],'LTV1'))
        _ltv2 = float(fileutil.getKeyword(parlist['data'],'LTV2'))
        if (_ltv1 != 0.) or (_ltv2 != 0.):
            _subarray = True
        _naxis1 = float(fileutil.getKeyword(parlist['data'],'NAXIS1'))
        _naxis2 = float(fileutil.getKeyword(parlist['data'],'NAXIS2'))
        if (_naxis1 < self.full_shape[0]) or (_naxis2 < self.full_shape[0]):
            _subarray = True
        return _subarray

    def setInstrumentParameters(self, instrpars, pri_header):
        """ This method overrides the superclass to set default values into
            the parameter dictionary, in case empty entries are provided.
        """
        if self._isNotValid (instrpars['gain'], instrpars['gnkeyword']):
            instrpars['gnkeyword'] = 'ATODGNA,ATODGNB,ATODGNC,ATODGND'
        if self._isNotValid (instrpars['rdnoise'], instrpars['rnkeyword']):
            instrpars['rnkeyword'] = 'READNSEA,READNSEB,READNSEC,READNSED'
        if self._isNotValid (instrpars['exptime'], instrpars['expkeyword']):
            instrpars['expkeyword'] = 'EXPTIME'
        if instrpars['crbit'] == None:
            instrpars['crbit'] = self.cr_bits_value
         
        self._gain      = self.getInstrParameter(instrpars['gain'], pri_header,
                                                 instrpars['gnkeyword'])
        self._rdnoise   = self.getInstrParameter(instrpars['rdnoise'], pri_header,
                                                 instrpars['rnkeyword'])
        self._exptime   = self.getInstrParameter(instrpars['exptime'], pri_header,
                                                 instrpars['expkeyword'])
        self._crbit     = instrpars['crbit']

        if self._gain == None or self._rdnoise == None or self._exptime == None:
            print 'ERROR: invalid instrument task parameter'
            raise ValueError

        # Convert the science data to electrons if specified by the user.  Each
        # instrument class will need to define its own version of doUnitConversions
        if self.proc_unit == "electrons":
            self.doUnitConversions()

    def getflat(self):
        """

        Purpose
        =======
        Method for retrieving a detector's flat field.
        
        This method will return an array the same shape as the
        image.
        
        :units: electrons

        """

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

class WFCInputImage (ACSInputImage):

    def __init__(self, input, dqname, platescale, memmap=0,proc_unit="native"):
        ACSInputImage.__init__(self,input,dqname,platescale,memmap=0,proc_unit=proc_unit)
        self.instrument = 'ACS/WFC'
        self.full_shape = (4096,2048)
        self.platescale = platescale

        if ( self.extn == 'sci,1') : # get cte direction, which depends on which chip but is independent of amp 
            self.cte_dir = -1    
        if ( self.extn == 'sci,2') : 
            self.cte_dir = 1   

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

    def setInstrumentParameters(self, instrpars, pri_header):
        """ Sets the instrument parameters.
        """
        if self._isNotValid (instrpars['gain'], instrpars['gnkeyword']):
            instrpars['gnkeyword'] = None
        if self._isNotValid (instrpars['rdnoise'], instrpars['rnkeyword']):
            instrpars['rnkeyword'] = None
        if self._isNotValid (instrpars['exptime'], instrpars['expkeyword']):
            instrpars['expkeyword'] = 'EXPTIME'
        if instrpars['crbit'] == None:
            instrpars['crbit'] = self.cr_bits_value
      
        self._exptime   = self.getInstrParameter(instrpars['exptime'], pri_header,
                                                 instrpars['expkeyword'])
        self._crbit     = instrpars['crbit']

        if self._exptime == None:
            print 'ERROR: invalid instrument task parameter'
            raise ValueError

        # We need to treat Read Noise and Gain as a special case since it is 
        # not populated in the SBC primary header for the MAMA
        if (instrpars['rnkeyword'] != None):
            self._rdnoise   = self.getInstrParameter(instrpars['rdnoise'], pri_header,
                                                     instrpars['rnkeyword'])                                                 
        else:
            self._rdnoise = None
        if (instrpars['gnkeyword'] != None):
            self._gain = self.getInstrParameter(instrpars['gain'], pri_header,
                                                     instrpars['gnkeyword'])
        else:
            self._gain = None
 

        if self._exptime == None:
            print 'ERROR: invalid instrument task parameter'
            raise ValueError

        # We need to determine if the user has used the default readnoise/gain value
        # since if not, they will need to supply a gain/readnoise value as well                
        usingDefaultGain = False
        usingDefaultReadnoise = False
        if (instrpars['gnkeyword'] == None):
            usingDefaultGain = True
        if (instrpars['rnkeyword'] == None):
            usingDefaultReadnoise = True

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
