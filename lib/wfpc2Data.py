#
#   Authors: Warren Hack, Ivo Busko, Christopher Hanley
#   Program: wfpc2_input.py
#   Purpose: Class used to model WFPC2 specific instrument data.
import pyfits
import numpy as np

from pytools import fileutil

from imageObject import imageObject
from staticMask import constructFilename
import buildmask

# Translation table for any image that does not use the DQ extension of the MEF
# for the DQ array.
DQ_EXTNS = {'c0h':'sdq','c0f':'sci'}

#### Calibrated gain and readnoise values for each chip
WFPC2_GAINS = { 1:{7:[7.12,5.24],15:[13.99,7.02]},
                2:{7:[7.12,5.51],15:[14.50,7.84]},
                3:{7:[6.90,5.22],15:[13.95,6.99]},
                4:{7:[7.10,5.19],15:[13.95,8.32]}}
WFPC2_DETECTOR_NAMES = {1:"PC",2:"WF2",3:"WF3",4:"WF4"}

class WFPC2InputImage (imageObject):

    SEPARATOR = '_'

    def __init__(self, filename, group=None):
        imageObject.__init__(self,filename, group=group)
        # define the cosmic ray bits value to use in the dq array
        self.cr_bits_value = 4096
        self._instrument=self._image["PRIMARY"].header["INSTRUME"]        
        self._detector=self._image["PRIMARY"].header["DETECTOR"]  
        self._effGain = 1

        # Attribute defining the pixel dimensions of WFPC2 chips.
        self.full_shape = (800,800)
        
        # Reference Plate Scale used for updates to MDRIZSKY, we should get this from the wcs class
        #self.refplatescale = 0.0996 # arcsec / pixel

        for chip in range(1,self._numchips+1,1):
            self._assignSignature(chip) #this is used in the static mask
            self._image[self.scienceExt,chip].cte_dir = -1 # independent of amp, chip   
            
    def find_DQ_extension(self):
        ''' Return the suffix for the data quality extension and the name of the file
            which that DQ extension should be read from.
        '''
        dqfile = None
        # Look for additional file with DQ array, primarily for WFPC2 data
        indx = self._filename.find('.fits')
        suffix = self._filename[indx-4:indx]
        dqfile = self._filename.replace(suffix[:3],'_c1')
        dq_suffix = DQ_EXTNS[suffix[1:]]

        return dqfile,dq_suffix

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
        if (instrpars['rnkeyword'] == None or instrpars['rnkeyword'] == 'None'):
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

        extnum = self._interpretExten(exten)
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
        
    def buildMask(self,chip,bits=0,write=False):
        """ Build masks as specified in the user parameters found in the 
            configObj object.
        """
        
        sci_chip = self._image[self.scienceExt,chip]
        ### For WFPC2 Data, build mask files using:
        maskname = sci_chip.dqrootname+'_dqmask.fits'
        dqmask_name = buildmask.buildShadowMaskImage(sci_chip.dqfile,sci_chip.detnum,sci_chip.extnum,maskname,bitvalue=bits,binned=sci_chip.binned)
        sci_chip.dqmaskname = dqmask_name
        dqmask = pyfits.getdata(dqmask_name,0)
        return dqmask

    def _assignSignature(self, chip):
        """assign a unique signature for the image based 
           on the  instrument, detector, chip, and size
           this will be used to uniquely identify the appropriate
           static mask for the image
           
           this also records the filename for the static mask to the outputNames dictionary
           
        """
        sci_chip = self._image[self.scienceExt,chip] 
        ny=sci_chip._naxis1
        nx=sci_chip._naxis2
        detnum = sci_chip.detnum
        instr=self._instrument
        
        sig=(instr+WFPC2_DETECTOR_NAMES[detnum],(nx,ny),chip) #signature is a tuple
        sci_chip.signature=sig #signature is a tuple
        filename=constructFilename(sig)
        sci_chip.outputNames["staticMask"]=filename #this is the name of the static mask file

    def _setchippars(self):
        for chip in self.returnAllChips(extname=self.scienceExt): 
            try:
                chip._gain,chip._rdnoise = WFPC2_GAINS[chip.detnum][chip._headergain]
            except KeyError:
                raise ValueError, "! Header gain value is not valid for WFPC2"

       
