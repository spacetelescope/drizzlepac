#
#   Authors: Christopher Hanley
#   Program: wfc3_input.py
#   Purpose: Class used to model WFC3 specific instrument data.


from pytools import fileutil
from nictools import readTDD
from imageObject import imageObject
from staticMask import constructFilename

class WFC3InputImage(imageObject):

    SEPARATOR = '_'

    def __init__(self,filename=None,group=None):
        imageObject.__init__(self,filename,group=group)

        # define the cosmic ray bits value to use in the dq array
        self.cr_bits_value = 4096
        self._instrument=self._image["PRIMARY"].header["INSTRUME"]

    def doUnitConversions(self):
        # Effective gain to be used in the driz_cr step.  Since the
        # WFC3 images have already been converted to electrons,
        # the effective gain is 1.
        for chip in self.returnAllChips(extname=self.scienceExt): 
            chip._effGain=1.
            
                
    def _isSubArray(self):
        _subarray = False
        _ltv1 = float(self._image["PRIMARY"].header["LTV1"])
        _ltv2 = float(self._image["PRIMARY"].header["LTV2"])

        if (_ltv1 != 0.) or (_ltv2 != 0.):
            _subarray = True
 
        _naxis1 = float(self._image["PRIMARY"].header["NAXIS1"])
        _naxis2 = float(self._image["PRIMARY"].header["NAXIS2"])       
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
 
        for chip in self.returnAllChips(extname=self.scienceExt): 

            chip._gain      = self.getInstrParameter(instrpars['gain'], pri_header,
                                                     instrpars['gnkeyword'])
            chip._rdnoise   = self.getInstrParameter(instrpars['rdnoise'], pri_header,
                                                     instrpars['rnkeyword'])
            chip._exptime   = self.getInstrParameter(instrpars['exptime'], pri_header,
                                                     instrpars['expkeyword'])
            chip._crbit     = instrpars['crbit']

            if self._gain == None or self._rdnoise == None or self._exptime == None:
                print 'ERROR: invalid instrument task parameter'
                raise ValueError

            chip._assignSignature(chip) #this is used in the static mask                     

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

        # The keyword for WFC3 UVIS flat fields in the primary header of the flt
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

        return data

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
        
        sig=(instr+self._detector,(nx,ny),int(chip)) #signature is a tuple
        sci_chip.signature=sig #signature is a tuple
        filename=constructFilename(sig)
        sci_chip.outputNames["staticMask"]=filename #this is the name of the static mask file


class WFC3UVISInputImage(imageObject):

    def __init__(self,filename=None,group=None):
        WFC3InputImage.__init__(self,filename,group=group)

        # define the cosmic ray bits value to use in the dq array
        self.full_shape = (4096,2051)
        
        # get cte direction, which depends on which chip but is independent of amp 
        if(chip.extnum  == 1):
            chip.cte_dir = -1
        if(self.extnum  == 2):
            chip.cte_dir = 1
 
    def getdarkcurrent(self):
        """
        
        Purpose
        =======
        Return the dark current for the WFC3 UVIS detector.  This value
        will be contained within an instrument specific keyword.
        The value is in units of electrons.
        
        :units: electrons
        
        """
        
        darkcurrent = 0.
        
        try:
            darkcurrent = self._image["PRIMARY"].header['MEANDARK']
        except:
            str =  "#############################################\n"
            str += "#                                           #\n"
            str += "# Error:                                    #\n"
            str += "#   Cannot find the value for 'MEANDARK'    #\n"
            str += "#   in the image header.  WFC3 input images #\n"
            str += "#   are expected to have this header        #\n"
            str += "#   keyword.                                #\n"
            str += "#                                           #\n"
            str += "# Error occured in WFC3UVISInputImage class #\n"
            str += "#                                           #\n"
            str += "#############################################\n"
            raise ValueError, str
        
        
        return darkcurrent
 



class WFC3IRInputImage(imageObject):

    def __init__(self,filename=None,group=None):
        WFC3InputImage.__init__(self,filename,group=group)

        # define the cosmic ray bits value to use in the dq array
        self.full_shape = (1024,1024)
              
        # Effective gain to be used in the driz_cr step.  Since the
        # WFC3 images have already been converted to electrons the 
        # effective gain is 1.
        self._effGain = 1
 
        # no cte correction for WFC3/IR so set cte_dir=0.
        self.cte_dir = 0   

    def doUnitConversions(self):
        """WF3 IR data come out in electrons, and I imagine  the 
         photometry keywords will be calculated as such, so no image
         manipulation needs be done between native and electrons """
        for chip in self.returnAllChips(extname=self.scienceExt): 
            chip._effGain = 1
         

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
 
        for chip in self.returnAllChips(extname=self.scienceExt): 

            chip._gain      = self.getInstrParameter(instrpars['gain'], pri_header,
                                                     instrpars['gnkeyword'])
            chip._rdnoise   = self.getInstrParameter(instrpars['rdnoise'], pri_header,
                                                     instrpars['rnkeyword'])
            chip._exptime   = self.getInstrParameter(instrpars['exptime'], pri_header,
                                                     instrpars['expkeyword'])
            chip._crbit     = instrpars['crbit']

            if self._gain == None or self._rdnoise == None or self._exptime == None:
                print 'ERROR: invalid instrument task parameter'
                raise ValueError

           
            chip.cte_dir = 0 #no cte
           
 
            chip._assignSignature(chip) #this is used in the static mask                     

        # Convert the science data to electrons if specified by the user.  Each
        # instrument class will need to define its own version of doUnitConversions
        if self.proc_unit == "electrons":
            self.doUnitConversions()

    def getdarkimg(self):
        """
        
        Purpose
        =======
        Return an array representing the dark image for the detector.
        
        :units: cps
        
        """
        
        # First attempt to get the dark image specified by the "DARKFILE"
        # keyword in the primary keyword of the science data.
        try:
            filename = self.header["DARKFILE"]
            handle = fileutil.openImage(filename,mode='readonly',memmap=0)
            hdu = fileutil.getExtn(handle,extn="sci")
            darkobj = hdu.data[self.ltv2:self.size2,self.ltv1:self.size1]
            
        # If the darkfile cannot be located, create the dark image from
        # what we know about the detector dark current and assume a
        # constant dark current for the whole image.
        except:
            darkobj = np.ones(self.image_shape,dtype=self.image_dtype)*self.getdarkcurrent()
 
 
        return darkobj


    def getdarkcurrent(self):
        """
        
        Purpose
        =======
        Return the dark current for the WFC3/IR detector.  This value
        will be contained within an instrument specific keyword.
        
        :units: electrons
        
        """
        
        darkcurrent = 0
        
        try:
            darkcurrent = self._image["PRIMARY"].header['MEANDARK']
        except:
            str =  "#############################################\n"
            str += "#                                           #\n"
            str += "# Error:                                    #\n"
            str += "#   Cannot find the value for 'MEANDARK'    #\n"
            str += "#   in the image header.  WFC3 input images #\n"
            str += "#   are expected to have this header        #\n"
            str += "#   keyword.                                #\n"
            str += "#                                           #\n"
            str += "# Error occured in WFC3IRInputImage class   #\n"
            str += "#                                           #\n"
            str += "#############################################\n"
            raise ValueError, str
        
        
        return darkcurrent
        
