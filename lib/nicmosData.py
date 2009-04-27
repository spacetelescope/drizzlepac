#
#   Authors: Christopher Hanley, David Grumm, Megan Sosey
#   Program: nicmos_input.py
#   Purpose: Class used to model NICMOS specific instrument data.

from pytools import fileutil
from nictools import readTDD
import numpy as np
from imageObject import imageObject
from staticMask import constructFilename

class NICMOSInputImage(imageObject):

    SEPARATOR = '_'

    def __init__(self, filename=None):
        imageObject.__init__(self,filename)
        
        # define the cosmic ray bits value to use in the dq array
        self.cr_bits_value = 4096

        # Detector parameters, nic only has 1 detector in each file
        self.full_shape = (256,256)
        self._instrument=self._image['PRIMARY'].header["INSTRUME"]
         
        for chip in range(1,self._numchips+1,1):
            self._assignSignature(chip) #this is used in the static mask, static mask name also defined here, must be done after outputNames
            self._image[self.scienceExt,chip].cte_dir = 0   #no correction for nicmos
       
        self._effGain = 1. #get the specific gain from the detector subclass
            

    def _assignSignature(self, chip):
        """assign a unique signature for the image based 
           on the  instrument, detector, chip, and size
           this will be used to uniquely identify the appropriate
           static mask for the image
           
           this also records the filename for the static mask to the outputNames dictionary
           
        """
        instr=self._instrument
        ny=self._image[self.scienceExt,chip]._naxis1
        nx=self._image[self.scienceExt,chip]._naxis2
        detnum = self._image[self.scienceExt,chip].detnum
        
        sig=(instr+self._detector,(nx,ny),detnum) #signature is a tuple
        self._image[self.scienceExt,chip].signature=sig
        filename=constructFilename(sig)
        self._image[self.scienceExt,chip].outputNames["staticMask"]=filename #this is the name of the static mask file
        

    def doUnitConversions(self):
        """convert the data to electrons
        
        This converts all science data extensions and saves
        the results back to disk. We need to make sure
        the data inside the chips already in memory is altered as well
        
        """
        

         # Image information 
        _handle = fileutil.openImage(self._filename,mode='update',memmap=0) 

        for det in range(1,self._numchips,1):

            chip=self._image[self.scienceExt,det]
            
            if chip._gain != None:

                # Multiply the values of the sci extension pixels by the gain. 
                print "Converting %s from COUNTS to ELECTRONS"%(self._filename) 

                # If the exptime is 0 the science image will be zeroed out. 
                np.multiply(_handle[self.scienceExt,det].data,chip._gain,_handle[self.scienceExt,det].data)
                chip.data=_handle[det].data

                # Set the BUNIT keyword to 'electrons'
                _handle[det].header.update('BUNIT','ELECTRONS')

                # Update the PHOTFLAM value
                photflam = _handle[det].header['PHOTFLAM']
                _handle[det].header.update('PHOTFLAM',(photflam/self._gain()))
                
                chip._effGain = 1.
            
            else:
                print "Invalid gain value for data, no conversion done"
                return ValueError

        # Close the files and clean-up
        _handle.close() 

        self._effGain = 1.


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
        
        filename = self._image["PRIMARY"].header['FLATFILE']
        
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

        flat = (1.0/data) # The reference flat field is inverted

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
            return np.ones(self.full_shape,dtype=self.image_dtype)*self.getdarkcurrent()
        else:
            # Create Dark Object from AMPGLOW and Lineark Dark components
            darkobj = tddobj.getampglow() + tddobj.getlindark()
                        
            # Return the darkimage taking into account an subarray information available
            return darkobj[self.ltv2:self.size2,self.ltv1:self.size1]
 
    def isCountRate(self):
        """
        isCountRate: Method or IRInputObject used to indicate if the
        science data is in units of counts or count rate.  This method
        assumes that the keyword 'BUNIT' is in the header of the input
        FITS file.
        """
        
        if self.header.has_key('BUNIT'):       
            if self.header['BUNIT'].find("/") != -1:
                return True
        else:
            return False
        
    
class NIC1InputImage(NICMOSInputImage):

    def __init__(self, filename=None):
        NICMOSInputImage.__init__(self,filename)
        self._effGain = 1. #get the gain from the detector subclass        
        self._detector=self._image["PRIMARY"].header["CAMERA"]

        
    def _setDarkRate(self):
        self.darkrate = 0.08 #electrons/s
        if self.proc_unit == 'native':
            self.darkrate = self.darkrate / self._gain # DN/s

    def _setDefaultReadnoise(self):
        """ this could be updated to calculate the readnoise from the NOISFILE            
        """
        self._rdnoise = 26.0 # electrons
        if self.proc_unit == 'native':
            self._rdnoise = self._rdnoise / self._gain # ADU

    def setInstrumentParameters(self, instrpars):
        """ This method overrides the superclass to set default values into
            the parameter dictionary, in case empty entries are provided.
        """
        pri_header = self._image[0].header
        self.proc_unit = instrpars['proc_unit']

        if self._isNotValid (instrpars['gain'], instrpars['gnkeyword']):
            instrpars['gnkeyword'] = 'ADCGAIN'
        if self._isNotValid (instrpars['rdnoise'], instrpars['rnkeyword']):
            instrpars['rnkeyword'] = None
        if self._isNotValid (instrpars['exptime'], instrpars['expkeyword']):
            instrpars['expkeyword'] = 'EXPTIME'
        #if instrpars['crbit'] == None:
        #    instrpars['crbit'] = self.cr_bits_value
   
        for chip in self.returnAllChips(extname=self.scienceExt):
            #self._gain      = self.getInstrParameter(instrpars['gain'], pri_header,
            #                                         instrpars['gnkeyword'])
            self._gain=5.4 #measured
            
            self._rdnoise   = self.getInstrParameter(instrpars['rdnoise'], pri_header,
                                                     instrpars['rnkeyword'])
            self._exptime   = self.getInstrParameter(instrpars['exptime'], pri_header,
                                                     instrpars['expkeyword'])
            #self._crbit     = instrpars['crbit']

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

        self._effgain= self._gain


class NIC2InputImage(NICMOSInputImage):
    def __init__(self,filename=None):
        NICMOSInputImage.__init__(self,filename)
        self._effgain=1. #measured
        self._detector=self._image["PRIMARY"].header["CAMERA"]

    def _setDarkRate(self):
        self.darkrate = 0.08 #electrons/s
        if self.proc_unit == 'native':
            self.darkrate = self.darkrate / self._gain # DN/s

    def _setDefaultReadnoise(self):
        self._rdnoise = 26.0 #electrons
        if self.proc_unit == 'native':
            self._rdnoise = self._rdnoise/self._gain #ADU

    def setInstrumentParameters(self, instrpars):
        """ This method overrides the superclass to set default values into
            the parameter dictionary, in case empty entries are provided.
        """
        pri_header = self._image[0].header
        self.proc_unit = instrpars['proc_unit']

        if self._isNotValid (instrpars['gain'], instrpars['gnkeyword']):
            instrpars['gnkeyword'] = 'ADCGAIN'
        if self._isNotValid (instrpars['rdnoise'], instrpars['rnkeyword']):
            instrpars['rnkeyword'] = None
        if self._isNotValid (instrpars['exptime'], instrpars['expkeyword']):
            instrpars['expkeyword'] = 'EXPTIME'
        #if instrpars['crbit'] == None:
        #    instrpars['crbit'] = self.cr_bits_value
   
        for chip in self.returnAllChips(extname=self.scienceExt):
            #self._gain      = self.getInstrParameter(instrpars['gain'], pri_header,
            #                                         instrpars['gnkeyword'])
            self._gain=5.4 #measured
            
            self._rdnoise   = self.getInstrParameter(instrpars['rdnoise'], pri_header,
                                                     instrpars['rnkeyword'])
            self._exptime   = self.getInstrParameter(instrpars['exptime'], pri_header,
                                                     instrpars['expkeyword'])
            #self._crbit     = instrpars['crbit']

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

        self._effgain = self._gain


class NIC3InputImage(NICMOSInputImage):
    def __init__(self,filename=None):
        NICMOSInputImage.__init__(self,filename)
        self._detector=self._image["PRIMARY"].header["CAMERA"]

    def _setDarkRate(self):
        self.darkrate = 0.15 #electrons/s
        if self.proc_unit == 'native':
            self.darkrate = self.darkrate/self._gain #DN/s

    def _setDefaultReadnoise(self):
        self._rdnoise = 29.0 # electrons
        if self.proc_unit == 'native':
            self._rdnoise = self._rdnoise/self._gain #ADU

    def setInstrumentParameters(self, instrpars):
        """ This method overrides the superclass to set default values into
            the parameter dictionary, in case empty entries are provided.
        """
        pri_header = self._image[0].header
        self.proc_unit = instrpars['proc_unit']

        if self._isNotValid (instrpars['gain'], instrpars['gnkeyword']):
            instrpars['gnkeyword'] = 'ADCGAIN'
        if self._isNotValid (instrpars['rdnoise'], instrpars['rnkeyword']):
            instrpars['rnkeyword'] = None
        if self._isNotValid (instrpars['exptime'], instrpars['expkeyword']):
            instrpars['expkeyword'] = 'EXPTIME'
        #if instrpars['crbit'] == None:
        #    instrpars['crbit'] = self.cr_bits_value
   
        for chip in self.returnAllChips(extname=self.scienceExt):
            #self._gain      = self.getInstrParameter(instrpars['gain'], pri_header,
            #                                         instrpars['gnkeyword'])
            self._gain=6.5 #measured
            
            self._rdnoise   = self.getInstrParameter(instrpars['rdnoise'], pri_header,
                                                     instrpars['rnkeyword'])
            self._exptime   = self.getInstrParameter(instrpars['exptime'], pri_header,
                                                     instrpars['expkeyword'])
            #self._crbit     = instrpars['crbit']

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

        self._effgain = self._gain
