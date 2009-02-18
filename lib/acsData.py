#
#   Authors: Christopher Hanley, Warren Hack, Ivo Busko, David Grumm
#   Program: acs_input.py
#   Purpose: Class used to model ACS specific instrument data.

from pytools import fileutil
import numpy as np
from imageObject import imageObject



class ACSInputImage(imageObject):

    SEPARATOR = '_'

    def __init__(self,filename=None,proc_unit="native"):
        imageObject.__init__(self,filename)
        # define the cosmic ray bits value to use in the dq array
        self.cr_bits_value = 4096
        self._instrument=self._image["PRIMARY"].header["INSTRUME"]
    
        for chip in range(1,self._numchips+1,1):
            self._image[self.scienceExt,chip].darkcurrent=self.getdarkcurrent(chip)

    def _assignSignature(self, chip):
        """assign a unique signature for the image based 
           on the  instrument, detector, chip, and size
           this will be used to uniquely identify the appropriate
           static mask for the image
           
           this also records the filename for the static mask to the outputNames dictionary
           
        """
        ny=self._image[self.scienceExt,chip]._naxis1
        nx=self._image[self.scienceExt,chip]._naxis2
        detnum = self._image[self.scienceExt,chip].detnum
        
        self._image[self.scienceExt,chip].signature=(self._instrument+'/'+self._detector,(nx,ny),detnum) #signature is a tuple

        
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

    def setInstrumentParameters(self):
        """ This method overrides the superclass to set default values into
            the parameter dictionary, in case empty entries are provided.
            
            this should probably be moved to the sub instrument classes
            for each detector type?
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
        
        filename = self._image["PRIMARY"].header['PFLTFILE']
        
        try:
            handle = fileutil.openImage(filename,mode='readonly',memmap=0)
            hdu = fileutil.getExtn(handle,extn=self.extn)
            data = hdu.data[self.ltv2:self.size2,self.ltv1:self.size1]
            handle.close()
        except:
            try:
                #see if jref$ was appended to the filename
                handle = fileutil.openImage(filename[5:],mode='readonly',memmap=0)
                hdu = fileutil.getExtn(handle,extn=self.extn)
                data = hdu.data[self.ltv2:self.size2,self.ltv1:self.size1]
                handle.close()
            except:
                data = np.ones(self.image_shape,dtype=self.image_dtype)
                str = "Cannot find file "+filename+".  Treating flatfield constant value of '1'.\n"
                print str
        flat = data
        return flat


    def getdarkcurrent(self,chip):
        """
        
        Purpose
        =======
        Return the dark current for the ACS detector.  This value
        will be contained within an instrument specific keyword.
        The value in the image header will be converted to units
        of electrons.
        
        :units: electrons
        
        """
              
        darkcurrent=0.        
        try:
            darkcurrent = self._image[self.scienceExt,chip].header['MEANDARK']
        except:
            str =  "#############################################\n"
            str += "#                                           #\n"
            str += "# Error:                                    #\n"
            str += "#   Cannot find the value for 'MEANDARK'    #\n"
            str += "#   in the image header.  ACS input images  #\n"
            str += "#   are expected to have this header        #\n"
            str += "#   keyword.                                #\n"
            str += "#                                           #\n"
            str += "# Error occured in the ACSInputImage class  #\n"
            str += "#                                           #\n"
            str += "#############################################\n"
            raise ValueError, str
        
        
        return darkcurrent


class WFCInputImage(ACSInputImage):

    def __init__(self,filename=None):
        ACSInputImage.__init__(self,filename)
        self._detector=self._image['PRIMARY'].header["DETECTOR"]
        self.full_shape = (4096,2048)

        # get cte direction, which depends on which chip but is independent of amp 
        for chip in range(1,self._numchips+1,1):
            self._assignSignature(chip) #this is used in the static mask
            
            if ( chip == 1) :
                self._image[self.scienceExt,chip].cte_dir = -1    
            if ( chip == 2) : 
                self._image[self.scienceExt,chip].cte_dir = 1   

            
class HRCInputImage (ACSInputImage):

    def __init__(self, filename=None):
        ACSInputImage.__init__(self, input, filename)
        self.detector=self._image['PRIMARY'].header["DETECTOR"]
        self.full_shape = (1024,1024)

        for chip in range(1,self._numchips+1,1):
            self._assignSignature(chip) #this is used in the static mask
            
            amp=self._image[self.scienceExt,chip].header["CCDAMP"]
            if ( amp == 'A' or amp == 'B' ) : # cte direction depends on amp (but is independent of chip)
                self._image[self.scienceExt,chip].cte_dir = 1   
            if ( amp == 'C' or amp == 'D' ) :
                self._image[self.scienceExt,chip].cte_dir = -1   

class SBCInputImage (ACSInputImage):

    def __init__(self, filename=None):
        ACSInputImage.__init__(self,filename)
        self.full_shape = (1024,1024)
        self._detector=self._image['PRIMARY'].header["DETECTOR"]

        for chip in range(1,self._numchips+1,1):
            self._assignSignature(chip) #this is used in the static mask
            # no cte correction for SBC so set cte_dir=0.
            print('\nWARNING: No cte correction will be made for this SBC data.\n')
            self._image[self.scienceExt,chip].cte_dir = 0       


    def _setSBCchippars(self):
        self._setDefaultSBCGain()
        self._setDefaultSBCReadnoise()
     
    def _setDefaultSBCGain(self):
        self._gain = 1

    def _setDefaultSBCReadnoise(self):
        self._rdnoise = 0

    def _setInstrumentParameters(self):
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
