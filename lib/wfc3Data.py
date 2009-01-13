#
#   Authors: Christopher Hanley
#   Program: wfc3_input.py
#   Purpose: Class used to model WFC3 specific instrument data.

from pytools import fileutil
import numpy as n
from input_image import InputImage
from ir_input import IRInputImage

class WFC3UVISInputImage(InputImage):

    SEPARATOR = '_'

    def __init__(self, input,dqname,platescale,memmap=0,proc_unit="native"):
        InputImage.__init__(self,input,dqname,platescale,memmap=0,proc_unit=proc_unit)

        # define the cosmic ray bits value to use in the dq array
        self.cr_bits_value = 4096
        self.platescale = platescale
        
        # Effective gain to be used in the driz_cr step.  Since the
        # WFC3 UVIS images have already been converted to electrons,
        # the effective gain is 1.
        self._effGain = 1

        self.instrument = 'WFC3/UVIS'
        self.full_shape = (4096,2048)
        self.platescale = platescale

        # get direction, which depends on which chip but is independent of amp 
        if ( self.extn == 'sci,1') : 
            self.cte_dir = -1    
        if ( self.extn == 'sci,2') : 
            self.cte_dir = 1   
        
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

        # The keyword for WFC3 UVIS flat fields in the primary header of the flt
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


    def getdarkcurrent(self):
        """
        
        Purpose
        =======
        Return the dark current for the WFC3 UVIS detector.  This value
        will be contained within an instrument specific keyword.
        The value is in units of electrons.
        
        :units: electrons
        
        """
        
        darkcurrent = 0
        
        try:
            darkcurrent = self.header['MEANDARK']
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

class WFC3IRInputImage(IRInputImage):

    def __init__(self, input, dqname, platescale, memmap=0,proc_unit="native"):
        IRInputImage.__init__(self,input,dqname,platescale,memmap=0,proc_unit=proc_unit)
        
        # define the cosmic ray bits value to use in the dq array
        self.cr_bits_value = 4096
        
        # Effective gain to be used in the driz_cr step.  Since the
        # NICMOS images have already been converted to electrons the 
        # effective gain is 1.
        self._effGain = 1
 
        # no cte correction for WFC3/IR so set cte_dir=0.
        self.cte_dir = 0   

        self.instrument = 'WFC3/IR'
        self.full_shape = (1000,1000)
        self.platescale = platescale

    def doUnitConversions(self):
        # Image information 
        _handle = fileutil.openImage(self.name,mode='update',memmap=0) 
        _sciext = fileutil.getExtn(_handle,extn=self.extn)         

        # Set the BUNIT keyword to 'electrons'
        _handle[1].header.update('BUNIT','ELECTRONS')

        # Counts case 
        np.multiply(_sciext.data,self.getExpTime(),_sciext.data)
        
        # Close the files and clean-up
        _handle.close() 

    def updateMDRIZSKY(self,filename=None): 
        if (filename == None): 
            filename = self.name     
        try: 
            _handle = fileutil.openImage(filename,mode='update',memmap=0) 
        except IOError:
            raise IOError, "Unable to open %s for sky level computation"%filename 
        # Get the exposure time for the image.  If the exposure time of the image 
        # is 0, set the MDRIZSKY value to 0.  Otherwise update the MDRIZSKY value 
        # in units of electrons per second. 
        if (self.getExpTime() == 0.0): 
            str =  "*******************************************\n" 
            str += "*                                         *\n" 
            str += "* ERROR: Image EXPTIME = 0.               *\n" 
            str += "* MDRIZSKY header value cannot be         *\n" 
            str += "* converted to units of 'electrons/s'     *\n" 
            str += "* MDRIZSKY will be set to a value of '0'  *\n" 
            str += "*                                         *\n" 
            str =  "*******************************************\n" 
            _handle[0].header['MDRIZSKY'] = 0 
            print str 
        else:
            # Assume the MDRIZSKY keyword is in the primary header.  Try to update 
            # the header value
            if self.proc_unit == 'electrons': 
                skyvalue = self.getSubtractedSky()/self.getExpTime() 
            else: 
                skyvalue = self.getSubtractedSky() 
            print "Updating MDRIZSKY keyword to primary header with value %f"%(skyvalue) 
            _handle[0].header.update('MDRIZSKY',skyvalue)  
        _handle.close() 
        
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

        # The keyword for WFC3 IR flat fields in the primary header of the flt
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
            darkcurrent = self.header['MEANDARK']
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
        
