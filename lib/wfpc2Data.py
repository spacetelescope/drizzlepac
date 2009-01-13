#
#   Authors: Warren Hack, Ivo Busko, Christopher Hanley
#   Program: wfpc2_input.py
#   Purpose: Class used to model WFPC2 specific instrument data.

from pytools import fileutil
from input_image import InputImage
import numpy as np

class WFPC2InputImage (InputImage):

    SEPARATOR = '_'

    def __init__(self, input,dqname,platescale,memmap=0,proc_unit="native"):
        InputImage.__init__(self,input,dqname,platescale,memmap=0,proc_unit=proc_unit)
        # define the cosmic ray bits value to use in the dq array
        self.cr_bits_value = 4096
        
        self.platescale = platescale 

        self.cte_dir = -1    # independent of amp, chip   
        
        self._effGain = 1

        # Attribute defining the pixel dimensions of WFPC2 chips.
        self.full_shape = (800,800)
        
        # Reference Plate Scale used for updates to MDRIZSKY
        self.refplatescale = 0.0996 # arcsec / pixel

    def getEffGain(self):
        """
        
        Purpose
        =======
        Method used to return the effective gain of a instrument's
        detector.
        
        This method returns a single floating point value.

        """

        return self._effGain

    def setInstrumentParameters(self, instrpars, pri_header):
        """ This method overrides the superclass to set default values into
            the parameter dictionary, in case empty entries are provided.
        """
                
        if self._isNotValid (instrpars['gain'], instrpars['gnkeyword']):
            instrpars['gnkeyword'] = 'ATODGAIN'
    
        if self._isNotValid (instrpars['rdnoise'], instrpars['rnkeyword']):
            instrpars['rnkeyword'] = None
            
        if self._isNotValid (instrpars['exptime'], instrpars['expkeyword']):
            instrpars['expkeyword'] = 'EXPTIME'

        if instrpars['crbit'] == None:
            instrpars['crbit'] = self.cr_bits_value

        self._headergain      = self.getInstrParameter(instrpars['gain'], pri_header,
                                                 instrpars['gnkeyword'])    
        self._exptime   = self.getInstrParameter(instrpars['exptime'], pri_header,
                                                 instrpars['expkeyword'])
        self._crbit     = instrpars['crbit']
        
        # We need to treat Read Noise as a special case since it is 
        # not populated in the WFPC2 primary header
        if (instrpars['rnkeyword'] != None):
            self._rdnoise   = self.getInstrParameter(instrpars['rdnoise'], pri_header,
                                                     instrpars['rnkeyword'])                                                 
        else:
            self._rdnoise = None

        if self._headergain == None or self._exptime == None:
            print 'ERROR: invalid instrument task parameter'
            raise ValueError

        # We need to determine if the user has used the default readnoise/gain value
        # since if not, they will need to supply a gain/readnoise value as well        
        
        usingDefaultGain = False
        usingDefaultReadnoise = False
        if (instrpars['gnkeyword'] == 'ATODGAIN'):
            usingDefaultGain = True
        if (instrpars['rnkeyword'] == None):
            usingDefaultReadnoise = True

            
        # If the user has specified either the readnoise or the gain, we need to make sure
        # that they have actually specified both values.  In the default case, the readnoise
        # of the system depends on what the gain
        
        if usingDefaultReadnoise and usingDefaultGain:
            self._setchippars()
        elif usingDefaultReadnoise and not usingDefaultGain:
            raise ValueError, "ERROR: You need to supply readnoise information\n when not using the default gain for WFPC2."
        elif not usingDefaultReadnoise and usingDefaultGain:
            raise ValueError, "ERROR: You need to supply gain information when\n not using the default readnoise for WFPC2." 
        else:
            # In this case, the user has specified both a gain and readnoise values.  Just use them as is.
            self._gain = self._headergain
            print "Using user defined values for gain and readnoise"

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
        

        """

        # The keyword for WFPC2 flat fields in the primary header of the flt
        # file is FLATFILE.  This flat file is *not* already in the required 
        # units of electrons.
        
        filename = self.header['FLATFILE']
        
        try:
            handle = fileutil.openImage(filename,mode='readonly',writefits=False,memmap=0)
            hdu = fileutil.getExtn(handle,extn=self.grp)
            data = hdu.data[self.ltv2:self.size2,self.ltv1:self.size1]
        except:
            try:
                handle = fileutil.openImage(filename[5:],mode='readonly',writefits=False,memmap=0)
                hdu = fileutil.getExtn(handle,extn=self.grp)
                data = hdu.data[self.ltv2:self.size2,self.ltv1:self.size1]
            except:
                data = np.ones(self.image_shape,dtype=self.image_dtype)
                str = "Cannot find file "+filename+".  Treating flatfield constant value of '1'.\n"
                print str
        # For the WFPC2 flat we need to invert
        # for use in Multidrizzle
        flat = (1.0/data)
        return flat

    def doUnitConversions(self): 
        # Image information 
        _handle = fileutil.openImage(self.name,mode='update',memmap=0) 
        _sciext = fileutil.getExtn(_handle,extn=self.extn)         

        # Multiply the values of the sci extension pixels by the gain. 
        print "Converting %s from COUNTS to ELECTRONS"%(self.name) 
        # If the exptime is 0 the science image will be zeroed out. 
        np.multiply(_sciext.data,self.getGain(),_sciext.data)

        # Set the BUNIT keyword to 'electrons'
        _handle[1].header.update('BUNIT','ELECTRONS')

        # Update the PHOTFLAM value
        photflam = _handle[1].header['PHOTFLAM']
        _handle[1].header.update('PHOTFLAM',(photflam/self.getGain()))

        # Close the files and clean-up
        _handle.close() 

    def getdarkcurrent(self):
        """
        
        Purpose
        =======
        Return the dark current for the WFPC2 detector.  This value
        will be contained within an instrument specific keyword.
        The value in the image header will be converted to units
        of electrons.
        
        :units: counts/electrons
        
        """
        darkrate = 0.005 # electrons / s
        if self.proc_unit == 'native':
            darkrate = darkrate / self.getGain() #count/s
        
        try:
            darkcurrent = self.header['DARKTIME'] * darkrate
            
        except:
            str =  "#############################################\n"
            str += "#                                           #\n"
            str += "# Error:                                    #\n"
            str += "#   Cannot find the value for 'DARKTIME'    #\n"
            str += "#   in the image header.  WFPC2 input       #\n"
            str += "#   images are expected to have this header #\n"
            str += "#   keyword.                                #\n"
            str += "#                                           #\n"
            str += "# Error occured in the WFPC2InputImage class#\n"
            str += "#                                           #\n"
            str += "#############################################\n"
            raise ValueError, str
        
        
        return darkcurrent

    def getReadNoise(self):
        """
        
        Purpose
        =======
        Method for returning the readnoise of a detector (in counts).
        
        :units: counts/electrons
        
        """
        
        rn = self._rdnoise
        if self.proc_unit == 'native':
            rn = self._rdnoise / self.getGain()
        return rn

    def _setchippars(self):
        pass

    def _getCalibratedGain(self):
        return self._gain

    def getComputedSky(self):
        return (self._computedsky * (self.refplatescale / self.platescale)**2 )
        
    def setSubtractedSky(self,newValue):
        self._subtractedsky = (newValue / (self.refplatescale /  self.platescale)**2)
            
    def subtractSky(self):
        try:
            try:
                _handle = fileutil.openImage(self.name,mode='update',memmap=0)
                _sciext = fileutil.getExtn(_handle,extn=self.extn)
                print "%s (computed sky,subtracted sky)  : (%f,%f)"%(self.name,self.getComputedSky(),self.getSubtractedSky()*(self.refplatescale / self.platescale)**2)
                np.subtract(_sciext.data,self.getSubtractedSky(),_sciext.data)
            except:
                raise IOError, "Unable to open %s for sky subtraction"%self.name
        finally:
            _handle.close()
            del _sciext,_handle

    def updateMDRIZSKY(self,filename=None):
    
        if (filename == None):
            filename = self.name
            
        try:
            _handle = fileutil.openImage(filename,mode='update',memmap=0)
        except:
            raise IOError, "Unable to open %s for sky level computation"%filename

         # Compute the sky level subtracted from all the WFPC2 detectors based upon the reference plate scale.
        skyvalue = (self.getSubtractedSky()  * (self.refplatescale/self.platescale)**2)
        if self.proc_unit == 'electrons':
            skyvalue = skyvalue / self.getGain()

        print "Updating MDRIZSKY keyword in primary header with value %f"%skyvalue
        _handle[0].header.update('MDRIZSKY',skyvalue, comment="Sky value subtracted by Multidrizzle")
        _handle.close()


class WF2InputImage (WFPC2InputImage):

    def __init__(self, input, dqname, platescale, memmap=0,proc_unit="native"):
        WFPC2InputImage.__init__(self,input,dqname,platescale, memmap=0,proc_unit=proc_unit)
        self.instrument = 'WFPC2/WF2'
        self.platescale = platescale #0.0996 #arcsec / pixel
        
    def _setchippars(self):
        if self._headergain == 7:
            self._gain    = 7.12
            self._rdnoise = 5.51  # electrons
        elif self._headergain == 15:
            self._gain    = 14.50
            self._rdnoise = 7.84 # electrons
        else:
            raise ValueError, "! Header gain value is not valid for WFPC2"

class WF3InputImage (WFPC2InputImage):

    def __init__(self, input, dqname, platescale, memmap=0,proc_unit="native"):
        WFPC2InputImage.__init__(self, input, dqname, platescale, memmap=0,proc_unit=proc_unit)
        self.instrument = 'WFPC2/WF3'
        self.platescale = platescale #0.0996 #arcsec / pixel

    def _setchippars(self):
        if self._headergain == 7:
            self._gain    = 6.90
            self._rdnoise = 5.22  #electrons
        elif self._headergain == 15:
            self._gain    = 13.95
            self._rdnoise = 6.99 #electrons
        else:
            raise ValueError, "! Header gain value is not valid for WFPC2"

class WF4InputImage (WFPC2InputImage):

    def __init__(self, input, dqname, platescale, memmap=0,proc_unit="native"):
        WFPC2InputImage.__init__(self, input, dqname, platescale, memmap=0,proc_unit=proc_unit)
        self.instrument = 'WFPC2/WF4'
        self.platescale = platescale #0.0996 #arcsec / pixel

    def _setchippars(self):
        if self._headergain == 7:
            self._gain    = 7.10
            self._rdnoise = 5.19   #electrons
        elif self._headergain == 15:
            self._gain    = 13.95
            self._rdnoise = 8.32 # electrons
        else:
            raise ValueError, "! Header gain value is not valid for WFPC2"

class PCInputImage (WFPC2InputImage):

    def __init__(self, input, dqname, platescale, memmap=0,proc_unit="native"):
        WFPC2InputImage.__init__(self,input,dqname,platescale,memmap=0,proc_unit=proc_unit)
        self.instrument = 'WFPC2/PC'
        self.platescale = platescale #0.0455 #arcsec / pixel

    def _setchippars(self):
        if self._headergain == 7:
            self._gain    = 7.12
            self._rdnoise = 5.24   # electrons
        elif self._headergain == 15:
            self._gain    = 13.99
            self._rdnoise = 7.02 # electrons
        else:
            raise ValueError, "! Header gain value is not valid for WFPC2"
