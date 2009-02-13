#
#   Authors: Christopher Hanley, David Grumm
#   Program: nicmos_input.py
#   Purpose: Class used to model NICMOS specific instrument data.

from pytools import fileutil
from nictools import readTDD
import numpy as np

from ir_input import IRInputImage
from input_image import InputImage


class NICMOSInputImage(imageObject):

    SEPARATOR = '_'

    def __init__(self, input,dqname,platescale,memmap=0,proc_unit="native"):
        IRInputImage.__init__(self,input,dqname,platescale,memmap=0,proc_unit=proc_unit)
        
        # define the cosmic ray bits value to use in the dq array
        self.cr_bits_value = 4096

        # Detector parameters
        self.platescale = platescale
        self.full_shape = (256,256)
         
        # no cte correction for NICMOS so set cte_dir=0.
        self.cte_dir = 0   

        self._effGain = 1

    def updateMDRIZSKY(self,filename=None): 
        if (filename == None): 
            filename = self.name     
        try: 
            _handle = fileutil.openImage(filename,mode='update',memmap=0) 
        except IOError:
            raise IOError, "Unable to open %s for sky level computation"%filename 
        # Get the exposure time for the image.  If the exposure time of the image 
        # is 0, set the MDRIZSKY value to 0.  Otherwise update the MDRIZSKY value 
        # in units of counts per second. 
        if (self.getExpTime() == 0.0): 
            str =  "*******************************************\n" 
            str += "*                                         *\n" 
            str += "* ERROR: Image EXPTIME = 0.               *\n" 
            str += "* MDRIZSKY header value cannot be         *\n" 
            str += "* converted to units of 'counts/s'        *\n" 
            str += "* MDRIZSKY will be set to a value of '0'  *\n" 
            str += "*                                         *\n" 
            str =  "*******************************************\n" 
            _handle[0].header['MDRIZSKY'] = 0 
            print str 
        else:
            # Assume the MDRIZSKY keyword is in the primary header.  Try to update 
            # the header value
            if (_handle[0].header['UNITCORR'].strip() == 'PERFORM'): 
                skyvalue = self.getSubtractedSky()/self.getExpTime() 
            else: 
                skyvalue = self.getSubtractedSky() 
            # We need to convert back to native units if computations were done in electrons
            if self.proc_unit != "native":
                skyvalue = skyvalue/self.getGain()
            print "Updating MDRIZSKY keyword to primary header with value %f"%(skyvalue) 
            _handle[0].header.update('MDRIZSKY',skyvalue)  
        _handle.close() 

    def doUnitConversions(self): 
        # Image information        
        _handle = fileutil.openImage(self.name,mode='update',memmap=0) 
        _sciext = fileutil.getExtn(_handle,extn=self.extn)         

        # Determine if Multidrizzle is in units of counts/second or counts 
        # 
        # Counts per second case 
        if (_handle[0].header['UNITCORR'].strip() == 'PERFORM'):         
            # Multiply the values of the sci extension pixels by the gain. 
            print "Converting %s from COUNTS/S to ELECTRONS"%(self.name) 
            # If the exptime is 0 the science image will be zeroed out. 
            conversionFactor = (self.getExpTime() * self.getGain())

        # Counts case 
        else:
            # Multiply the values of the sci extension pixels by the gain. 
            print "Converting %s from COUNTS to ELECTRONS"%(self.name) 
            # If the exptime is 0 the science image will be zeroed out. 
            conversionFactor = (self.getGain())  

        np.multiply(_sciext.data,conversionFactor,_sciext.data)
        
        # Set the BUNIT keyword to 'electrons'
        _handle[0].header.update('BUNIT','ELECTRONS')

        # Update the PHOTFLAM value
        photflam = _handle[0].header['PHOTFLAM']
        _handle[0].header.update('PHOTFLAM',(photflam/self.getGain()))
        
        # Close the files and clean-up
        _handle.close() 

    def setInstrumentParameters(self, instrpars, pri_header):
        """ This method overrides the superclass to set default values into
            the parameter dictionary, in case empty entries are provided.
        """
        if self._isNotValid (instrpars['gain'], instrpars['gnkeyword']):
            instrpars['gnkeyword'] = 'ADCGAIN'
        if self._isNotValid (instrpars['rdnoise'], instrpars['rnkeyword']):
            instrpars['rnkeyword'] = None
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

        if self._gain == None or self._exptime == None:
            print 'ERROR: invalid instrument task parameter'
            raise ValueError

        # We need to treat Read Noise as a special case since it is 
        # not populated in the NICMOS primary header
        if (instrpars['rnkeyword'] != None):
            self._rdnoise   = self.getInstrParameter(instrpars['rdnoise'], pri_header,
                                                     instrpars['rnkeyword'])                                                 
        else:
            self._rdnoise = None


        # We need to determine if the user has used the default readnoise/gain value
        # since if not, they will need to supply a gain/readnoise value as well        
        
        usingDefaultReadnoise = False
        if (instrpars['rnkeyword'] == None):
            usingDefaultReadnoise = True
            
        # Set the default readnoise values based upon the amount of user input given.
        
        # User supplied no readnoise information
        if usingDefaultReadnoise:
            # Set the default gain and readnoise values
            self._setchippars()
        
        # Set the darkrate for the chips
        self._setDarkRate()
        
        # Convert the science data to electrons if specified by the user.  Each
        # instrument class will need to define its own version of doUnitConversions
        if self.proc_unit == "electrons":
            self.doUnitConversions()

    def _setchippars(self):
        self._setDefaultReadnoise()
                
    def getflat(self):
        """

        Purpose
        =======
        Method for retrieving a detector's flat field.
        
        This method will return an array the same shape as the
        image.

        :units: cps

        """

        # The keyword for NICMOS flat fields in the primary header of the flt
        # file is FLATFILE.  This flat file is not already in the required 
        # units of electrons.
        
        filename = self.header['FLATFILE']
        
        try:
            handle = fileutil.openImage(filename,mode='readonly',memmap=0)
            hdu = fileutil.getExtn(handle,extn=self.grp)
            data = hdu.data[self.ltv2:self.size2,self.ltv1:self.size1]
        except:
            try:
                handle = fileutil.openImage(filename[5:],mode='readonly',memmap=0)
                hdu = fileutil.getExtn(handle,extn=self.grp)
                data = hdu.data[self.ltv2:self.size2,self.ltv1:self.size1]
            except:
                data = np.ones(self.image_shape,dtype=self.image_dtype)
                str = "Cannot find file "+filename+".  Treating flatfield constant value of '1'.\n"
                print str

        flat = (1.0/data) # The flat field is normalized to unity.

        return flat
        

    def getdarkcurrent(self):
        """
        
        Purpose
        =======
        Return the dark current for the NICMOS detectors.
        
        :units: cps
        
        """
                
        try:
            darkcurrent = self.header['exptime'] * self.darkrate
            
        except:
            str =  "#############################################\n"
            str += "#                                           #\n"
            str += "# Error:                                    #\n"
            str += "#   Cannot find the value for 'EXPTIME'     #\n"
            str += "#   in the image header.  NICMOS input      #\n"
            str += "#   images are expected to have this header #\n"
            str += "#   keyword.                                #\n"
            str += "#                                           #\n"
            str += "#Error occured in the NICMOSInputImage class#\n"
            str += "#                                           #\n"
            str += "#############################################\n"
            raise ValueError, str
        
        
        return darkcurrent
        
    def getdarkimg(self):
        """
        
        Purpose
        =======
        Return an array representing the dark image for the detector.
        
        :units: cps
        
        """

        # Read the temperature dependeant dark file.  The name for the file is taken from
        # the TEMPFILE keyword in the primary header.
        tddobj = readTDD.fromcalfile(self.name)

        if tddobj == None:
            return np.ones(self.image_shape,dtype=self.image_dtype)*self.getdarkcurrent()
        else:
            # Create Dark Object from AMPGLOW and Lineark Dark components
            darkobj = tddobj.getampglow() + tddobj.getlindark()
                        
            # Return the darkimage taking into account an subarray information available
            return darkobj[self.ltv2:self.size2,self.ltv1:self.size1]
        
    
class NIC1InputImage(NICMOSInputImage):

    def __init__(self, input, dqname, platescale, memmap=0,proc_unit="native"):
        NICMOSInputImage.__init__(self,input,dqname,platescale,memmap=0,proc_unit=proc_unit)
        self.instrument = 'NICMOS/1'
        
    def _setDarkRate(self):
        self.darkrate = 0.08 #electrons/s
        if self.proc_unit == 'native':
            self.darkrate = self.darkrate / self.getGain() # DN/s

    def _setDefaultReadnoise(self):
        self._rdnoise = 26.0 # electrons
        if self.proc_unit == 'native':
            self._rdnoise = self._rdnoise / self.getGain() # ADU

class NIC2InputImage(NICMOSInputImage):
    def __init__(self, input, dqname, platescale, memmap=0,proc_unit="native"):
        NICMOSInputImage.__init__(self,input,dqname,platescale,memmap=0,proc_unit=proc_unit)
        self.instrument = 'NICMOS/2'
        
    def _setDarkRate(self):
        self.darkrate = 0.08 #electrons/s
        if self.proc_unit == 'native':
            self.darkrate = self.darkrate / self.getGain() # DN/s

    def _setDefaultReadnoise(self):
        self._rdnoise = 26.0 #electrons
        if self.proc_unit == 'native':
            self._rdnoise = self._rdnoise/self.getGain() #ADU

class NIC3InputImage(NICMOSInputImage):
    def __init__(self, input, dqname, platescale, memmap=0,proc_unit="native"):
        NICMOSInputImage.__init__(self,input,dqname,platescale,memmap=0,proc_unit=proc_unit)
        self.instrument = 'NICMOS/3'
        
    def _setDarkRate(self):
        self.darkrate = 0.15 #electrons/s
        if self.proc_unit == 'native':
            self.darkrate = self.darkrate/self.getGain() #DN/s

    def _setDefaultReadnoise(self):
        self._rdnoise = 29.0 # electrons
        if self.proc_unit == 'native':
            self._rdnoise = self._rdnoise/self.getGain() #ADU
