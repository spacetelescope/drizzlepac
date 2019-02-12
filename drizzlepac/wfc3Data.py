"""
`wfc3Data` module provides classes used to import WFC3 specific instrument data.

:Authors: Megan Sosey, Christopher Hanley

:License: :doc:`LICENSE`

"""
from stsci.tools import fileutil
from nictools import readTDD
from .imageObject import imageObject
import numpy as np

class WFC3InputImage(imageObject):

    SEPARATOR = '_'

    def __init__(self, filename=None, group=None):
        super().__init__(filename, group=group)

        # define the cosmic ray bits value to use in the dq array
        self.cr_bits_value = 4096
        self._instrument=self._image["PRIMARY"].header["INSTRUME"]
        self.flatkey = 'PFLTFILE'

    def _assignSignature(self, chip):
        """assign a unique signature for the image based
           on the  instrument, detector, chip, and size
           this will be used to uniquely identify the appropriate
           static mask for the image

           this also records the filename for the static mask to the outputNames dictionary

        """
        sci_chip = self._image[self.scienceExt,chip]
        ny = sci_chip._naxis1
        nx = sci_chip._naxis2
        detnum = sci_chip.detnum
        instr = self._instrument
        sig = (instr + self._detector, (nx, ny), int(chip)) #signature is a tuple
        sci_chip.signature=sig #signature is a tuple


class WFC3UVISInputImage(WFC3InputImage):
    def __init__(self,filename=None,group=None):
        super().__init__(filename, group=group)

        # define the cosmic ray bits value to use in the dq array
        self.full_shape = (4096,2051)
        self._detector = self._image["PRIMARY"].header["DETECTOR"]

        # get cte direction, which depends on which chip but is independent of amp
        for chip in range(1,self._numchips+1,1):
            self._assignSignature(chip) #this is used in the static mask

            if ( chip == 1) :
                self._image[self.scienceExt,chip].cte_dir = -1
            if ( chip == 2) :
                self._image[self.scienceExt,chip].cte_dir = 1
            self._image[self.scienceExt,chip].darkcurrent = self.getdarkcurrent(chip)



    def doUnitConversions(self):
        # Effective gain to be used in the driz_cr step.  Since the
        # WFC3 images have already been converted to electrons,
        # the effective gain is 1.
        for chip in self.returnAllChips(extname=self.scienceExt):
            chip._effGain=1.

    def setInstrumentParameters(self, instrpars):
        """ This method overrides the superclass to set default values into
            the parameter dictionary, in case empty entries are provided.
        """
        pri_header = self._image[0].header

        if len(instrpars) == 0:
            instrpars['proc_unit']='native'
            instrpars['gain']=''
            instrpars['rdnoise']=''
            instrpars['exptime']=''
            instrpars['gnkeyword']=''
            instrpars['rnkeyword']=''
            instrpars['expkeyword']=''

        if self._isNotValid (instrpars['gain'], instrpars['gnkeyword']):
            instrpars['gnkeyword'] = 'ATODGNA,ATODGNB,ATODGNC,ATODGND'
        if self._isNotValid (instrpars['rdnoise'], instrpars['rnkeyword']):
            instrpars['rnkeyword'] = 'READNSEA,READNSEB,READNSEC,READNSED'
        if self._isNotValid (instrpars['exptime'], instrpars['expkeyword']):
            instrpars['expkeyword'] = 'EXPTIME'

        self.proc_unit = instrpars['proc_unit']

        for chip in self.returnAllChips(extname=self.scienceExt):

            chip._gain      = self.getInstrParameter(instrpars['gain'], pri_header,
                                                     instrpars['gnkeyword'])
            chip._rdnoise   = self.getInstrParameter(instrpars['rdnoise'], pri_header,
                                                     instrpars['rnkeyword'])
            chip._exptime   = self.getInstrParameter(instrpars['exptime'], pri_header,
                                                     instrpars['expkeyword'])
            chip._effGain=chip._gain

            if chip._gain is None or chip._rdnoise is None or chip._exptime is None:
                print('ERROR: invalid instrument task parameter')
                raise ValueError

        # Convert the science data to electrons.
        self.doUnitConversions()

    def getdarkcurrent(self,chip):
        """
        Return the dark current for the WFC3 UVIS detector.  This value
        will be contained within an instrument specific keyword.

        Returns
        -------
        darkcurrent: float
            The dark current value with **units of electrons**.
        """
        darkcurrent = 0.

        try:
            darkcurrent = self._image[self.scienceExt, chip].header['MEANDARK']

        except:
            msg =  "#############################################\n"
            msg += "#                                           #\n"
            msg += "# Error:                                    #\n"
            msg += "#   Cannot find the value for 'MEANDARK'    #\n"
            msg += "#   in the image header.  WFC3 input images #\n"
            msg += "#   are expected to have this header        #\n"
            msg += "#   keyword.                                #\n"
            msg += "#                                           #\n"
            msg += "# Error occured in WFC3UVISInputImage class #\n"
            msg += "#                                           #\n"
            msg += "#############################################\n"
            raise ValueError(msg)

        return darkcurrent


class WFC3IRInputImage(WFC3InputImage):

    def __init__(self, filename=None, group=None):
        super().__init__(filename, group=group)
        self.timeExt = 'TIME'

        # define the cosmic ray bits value to use in the dq array
        self.full_shape = (1024, 1024)
        self._detector=self._image["PRIMARY"].header["DETECTOR"]
        self.native_units = 'ELECTRONS/S'

        # Effective gain to be used in the driz_cr step.  Since the
        # WFC3 images have already been converted to electrons the
        # effective gain is 1.
        self._effGain = 1.

        # no cte correction for WFC3/IR so set cte_dir=0.
        self.cte_dir = 0
        self._image[self.scienceExt,1].cte_dir = 0
        self._image[self.scienceExt,1].darkcurrent = self.getdarkcurrent()

    def doUnitConversions(self):
        """WF3 IR data come out in electrons, and I imagine  the
         photometry keywords will be calculated as such, so no image
         manipulation needs be done between native and electrons """
         # Image information
        _handle = fileutil.openImage(self._filename, mode='readonly', memmap=False)

        for chip in self.returnAllChips(extname=self.scienceExt):
            conversionFactor = 1.0
            if '/S' in chip._bunit:
                conversionFactor = chip._exptime
            else:
                print("Input %s[%s,%d] already in units of ELECTRONS"
                      %(self._filename,self.scienceExt,chip._chip))

            chip._effGain = 1.0# chip._gain #1.
            chip._conversionFactor = conversionFactor #1.

        _handle.close()
        self._effGain= 1.0 #conversionFactor #1.0

    def setInstrumentParameters(self, instrpars):
        """ This method overrides the superclass to set default values into
            the parameter dictionary, in case empty entries are provided.
        """
        pri_header = self._image[0].header

        if len(instrpars) == 0:
            instrpars['proc_unit']='native'
            instrpars['gain']=''
            instrpars['rdnoise']=''
            instrpars['exptime']=''
            instrpars['gnkeyword']=''
            instrpars['rnkeyword']=''
            instrpars['expkeyword']=''

        if self._isNotValid (instrpars['gain'], instrpars['gnkeyword']):
            instrpars['gnkeyword'] = 'ATODGNA,ATODGNB,ATODGNC,ATODGND'
        if self._isNotValid (instrpars['rdnoise'], instrpars['rnkeyword']):
            instrpars['rnkeyword'] = 'READNSEA,READNSEB,READNSEC,READNSED'
        if self._isNotValid (instrpars['exptime'], instrpars['expkeyword']):
            instrpars['expkeyword'] = 'EXPTIME'

        self.proc_unit = instrpars['proc_unit']

        for chip in self.returnAllChips(extname=self.scienceExt):
            chip._gain      = self.getInstrParameter(instrpars['gain'], pri_header,
                                                     instrpars['gnkeyword'])
            chip._rdnoise   = self.getInstrParameter(instrpars['rdnoise'], pri_header,
                                                     instrpars['rnkeyword'])
            chip._exptime   = self.getInstrParameter(instrpars['exptime'], pri_header,
                                                     instrpars['expkeyword'])
            chip._effGain= 1

            if chip._gain is None or chip._rdnoise is None or chip._exptime is None:
                print('ERROR: invalid instrument task parameter')
                raise ValueError

            self._assignSignature(chip.extnum) #this is used in the static mask

        #Convert from ELECTRONS/S to ELECTRONS
        self.doUnitConversions()

    def getexptimeimg(self,chip):
        """
        Return an array representing the exposure time per pixel for the detector.

        Returns
        -------
        dark: array
            Exposure time array in the same shape as the input image

        """
        return self._image[self.timeExt,chip].data

    def getdarkimg(self,chip):
        """
        Return an array representing the dark image for the detector.

        Returns
        -------
        dark: array
            Dark image array in the same shape as the input image with **units of cps**

        """
        sci_chip = self._image[self.scienceExt,chip]

        # First attempt to get the dark image specified by the "DARKFILE"
        # keyword in the primary keyword of the science data.
        try:
            filename = self.header["DARKFILE"]
            handle = fileutil.openImage(filename, mode='readonly', memmap=False)
            hdu = fileutil.getExtn(handle,extn="sci,1")
            darkobj = hdu.data[sci_chip.ltv2:sci_chip.size2,sci_chip.ltv1:sci_chip.size1]

        # If the darkfile cannot be located, create the dark image from
        # what we know about the detector dark current and assume a
        # constant dark current for the whole image.
        except:
            darkobj = (np.ones(sci_chip.image_shape,
                               dtype=sci_chip.image_dtype) *
                       self.getdarkcurrent())
        return darkobj

    def getskyimg(self,chip):
        """
        Notes
        =====
        Return an array representing the sky image for the detector.  The value
        of the sky is what would actually be subtracted from the exposure by
        the skysub step.

        :units: electrons

        """
        sci_chip = self._image[self.scienceExt,chip]
        skyimg = np.ones(sci_chip.image_shape,dtype=sci_chip.image_dtype)*sci_chip.subtractedSky
        if sci_chip._conversionFactor != 1.0: # If units are not already ELECTRONS
            skyimg *= self.getexptimeimg(chip)
        return skyimg

    def getdarkcurrent(self):
        """
        Return the dark current for the WFC3/IR detector.  This value
        will be contained within an instrument specific keyword.

        Returns
        -------
        darkcurrent: float
            The dark current value in **units of electrons**.
        """
        darkcurrent = 0
        try:
            darkcurrent = self._image[self.scienceExt,1].header['MEANDARK']
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
            raise ValueError(str)

        return darkcurrent
