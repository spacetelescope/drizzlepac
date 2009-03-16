#
#   Authors: Warren Hack, Ivo Busko, Christopher Hanley
#   Program: wfpc2_input.py
#   Purpose: Class used to model WFPC2 specific instrument data.

from pytools import fileutil
import numpy as np
from imageObject import imageObject
from staticMask import constructFilename

#### Calibrated gain and readnoise values for each chip
WFPC2_GAINS = { 1:{7:[7.12,5.24],15:[13.99,7.02]},
                2:{7:[7.12,5.51],15:[14.50,7.84]},
                3:{7:[6.90,5.22],15:[13.95,6.99]},
                4:{7:[7.10,5.19],15:[13.95,8.32]}}
class WFPC2InputImage (imageObject):

    SEPARATOR = '_'

    def __init__(self, filename, group=None):
        imageObject.__init__(self,filename, group=group)
        # define the cosmic ray bits value to use in the dq array
        self.cr_bits_value = 4096
        self._instrument=self._image["PRIMARY"].header["INSTRUME"]        

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

    def setInstrumentParameters(self, instrpars):
        """ This method overrides the superclass to set default values into
            the parameter dictionary, in case empty entries are provided.
        """
        pri_header = self._image[0].header
        self.proc_unit = instrpars['proc_unit']
                
        if self._isNotValid (instrpars['gain'], instrpars['gnkeyword']):
            instrpars['gnkeyword'] = 'ATODGAIN'
    
        if self._isNotValid (instrpars['rdnoise'], instrpars['rnkeyword']):
            instrpars['rnkeyword'] = None
            
        if self._isNotValid (instrpars['exptime'], instrpars['expkeyword']):
            instrpars['expkeyword'] = 'EXPTIME'

        for chip in self.returnAllChips(extname=self.scienceExt): 
            chip._headergain    = self.getInstrParameter(instrpars['gain'], pri_header,
                                                     instrpars['gnkeyword'])    
            chip._exptime       = self.getInstrParameter(instrpars['exptime'], pri_header,
                                                     instrpars['expkeyword'])
        
            # We need to treat Read Noise as a special case since it is 
            # not populated in the WFPC2 primary header
            if (instrpars['rnkeyword'] != None):
                chip._rdnoise   = self.getInstrParameter(instrpars['rdnoise'], pri_header,
                                                         instrpars['rnkeyword'])                                                 
            else:
                chip._rdnoise = None

            if chip._headergain == None or chip._exptime == None:
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
            for chip in self.returnAllChips(extname=self.scienceExt): 
                chip._gain = chip._headergain
            print "Using user defined values for gain and readnoise"

        # Convert the science data to electrons if specified by the user.  Each
        # instrument class will need to define its own version of doUnitConversions
        if self.proc_unit == "electrons":
            self.doUnitConversions()

    def getflat(self,exten):
        """

        Purpose
        =======
        Method for retrieving a detector's flat field.
        
        This method will return an array the same shape as the
        image.
        

        """

        extnum = self.interpretExten(exten)
        chip = self._image[exten]
        
        # The keyword for WFPC2 flat fields in the primary header of the flt
        # file is FLATFILE.  This flat file is *not* already in the required 
        # units of electrons.
        
        filename = self._image["PRIMARY"].header['FLATFILE']
        
        try:
            handle = fileutil.openImage(filename,mode='readonly',writefits=False,memmap=0)
            hdu = fileutil.getExtn(handle,extn=extnum)
            data = hdu.data[chip.ltv2:chip.size2,chip.ltv1:chip.size1]
            handle.close()
        except:
            try:
                handle = fileutil.openImage(filename[5:],mode='readonly',writefits=False,memmap=0)
                hdu = fileutil.getExtn(handle,extn=self.grp)
                data = hdu.data[chip.ltv2:chip.size2,chip.ltv1:chip.size1]
                handle.close()
            except:
                data = np.ones((chip._naxis2,chip._naxis1),dtype=chip.image_dtype)
                str = "Cannot find file "+filename+".  Treating flatfield constant value of '1'.\n"
                print str
        # For the WFPC2 flat we need to invert
        # for use in Multidrizzle
        flat = (1.0/data)
        return flat

    def doUnitConversions(self):
        """ Apply unit conversions to all the chips, ignoring the group parameter.
            This insures that all the chips get the same conversions when this 
            gets done, even if only 1 chip was specified to be processed.
        """
        for chip in range(1,numchips+1,1):
            myext=self.scienceExt+","+str(chip)
            
            image = self._image[myext]
            #add the data back into the chip, leave it there til the end of this function          
            image.data = self.getData(myext)
            
            # Multiply the values of the sci extension pixels by the gain. 
            print "Converting %s from COUNTS to ELECTRONS"%(self._filename) 
            # If the exptime is 0 the science image will be zeroed out. 
            np.multiply(image.data,image._gain,image.data)

            # Set the BUNIT keyword to 'electrons'
            image.header.update('BUNIT','ELECTRONS')
            # Update the PHOTFLAM value
            photflam = image.header['PHOTFLAM']
            image.header.update('PHOTFLAM',(photflam/image._gain))
            
            # Write out converted data array to original FLT image
            self.updateData(myext,image.data)
            
        # Delete the converted arrays from memory now that they have been
        # written out 
        self.close()

    def getdarkcurrent(self,exten):
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
            darkrate = darkrate / self.getGain(exten) #count/s
        
        try:
            chip = self._image[exten]
            darkcurrent = chip.header['DARKTIME'] * darkrate
            
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

    def getReadNoise(self,exten):
        """
        
        Purpose
        =======
        Method for returning the readnoise of a detector (in counts).
        
        :units: counts/electrons
        
        """
        
        rn = self._image[exten]._rdnoise
        if self.proc_unit == 'native':
            rn = self._rdnoise / self.getGain(exten)
        return rn

    def _setchippars(self):
        for chip in self.returnAllChips(extname=self.scienceExt): 
            try:
                chip._gain,chip._rdnoise = WFPC2_GAINS[chip.detnum][chip._headergain]
            except KeyError:
                raise ValueError, "! Header gain value is not valid for WFPC2"

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


