"""
Class used to model NICMOS specific instrument data.

:Authors: Christopher Hanley, David Grumm, Megan Sosey

:License: :doc:`LICENSE`

"""
from stsci.tools import fileutil
from nictools import readTDD
import numpy as np
from .imageObject import imageObject

class NICMOSInputImage(imageObject):

    SEPARATOR = '_'

    def __init__(self, filename=None):
        super().__init__(filename)
        self.timeExt = 'TIME'

        # define the cosmic ray bits value to use in the dq array
        self.cr_bits_value = 4096

        # Detector parameters, nic only has 1 detector in each file
        self.full_shape = (256,256)
        self._instrument=self._image['PRIMARY'].header["INSTRUME"]
        self.native_units = 'COUNTS/S'

        self.flatkey = 'FLATFILE'

        for chip in range(1,self._numchips+1,1):
            self._image[self.scienceExt,chip].cte_dir = 0   #no correction for nicmos

        self._effGain = 1. #get the specific gain from the detector subclass

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

        sig=(instr+str(self._detector),(nx,ny),int(detnum)) #signature is a tuple
        sci_chip.signature=sig #signature is a tuple


    def doUnitConversions(self):
        """Convert the data to electrons

        This converts all science data extensions and saves
        the results back to disk. We need to make sure
        the data inside the chips already in memory is altered as well.

        """

         # Image information
        _handle = fileutil.openImage(self._filename, mode='readonly', memmap=False)

        for det in range(1,self._numchips+1,1):

            chip=self._image[self.scienceExt,det]

            if chip._gain is not None:

                #conversionFactor = (self.getExpTime() * self.getGain())
                conversionFactor = chip._gain
                if self.isCountRate():
                    conversionFactor *= chip._exptime
                    counts_str = 'COUNTS/S'
                else:
                    counts_str = 'COUNTS'

                # Multiply the values of the sci extension pixels by the gain.
                print("Converting %s[%s,%d] from %s to ELECTRONS"%(self._filename,self.scienceExt,det,counts_str))
                """
                # If the exptime is 0 the science image will be zeroed out.
                np.multiply(_handle[self.scienceExt,det].data,conversionFactor,_handle[self.scienceExt,det].data)
                #chip.data=_handle[self.scienceExt,det].data.copy()

                # Set the BUNIT keyword to 'electrons'
                chip.header.update('BUNIT','ELECTRONS')
                _handle[0].header.update('BUNIT','ELECTRONS')

                # Update the PHOTFLAM value
                photflam = _handle[0].header['PHOTFLAM']
                _handle[0].header.update('PHOTFLAM',(photflam/chip._gain))

                chip._effGain = 1.0
                """
                chip._effGain = chip._gain
                chip._conversionFactor = conversionFactor
            else:
                msg = "Invalid gain value for data, no conversion done"
                print(msg)
                raise ValueError(msg)

        # Close the files and clean-up
        _handle.close()

        self._effGain = conversionFactor #1.0


    def _setchippars(self):
        self._setDefaultReadnoise()

    def getexptimeimg(self,chip):
        """
        Return an array representing the exposure time per pixel for the detector.

        Returns
        -------
        dark: array
            Exposure time array in the same shape as the input image

        """
        return self._image[self.timeExt,chip].data

    def getflat(self, chip):
        """
        Method for retrieving a detector's flat field.

        Returns
        -------
        flat : array
            The flat field array in the same shape as the input image with **units of cps**.
        """
        # The reference flat field is inverted:
        flat = 1.0 / super().getflat(chip)
        return flat

    def getdarkcurrent(self):
        """
        Return the dark current for the NICMOS detectors.

        Returns
        -------
        darkcurrent : float
            Dark current value with **units of cps**.

        """

        try:
            darkcurrent = self._image[0].header['exptime'] * \
                            self._image[self.scienceExt,1]._darkrate
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
            raise ValueError(str)

        return darkcurrent

    def getdarkimg(self,chip):
        """
        Return an array representing the dark image for the detector.

        Returns
        -------
        dark : array
            The dark array in the same shape as the image with **units of cps**.

        """

        # Read the temperature dependeant dark file.  The name for the file is taken from
        # the TEMPFILE keyword in the primary header.
        tddobj = readTDD.fromcalfile(self.name)

        if tddobj is None:
            return np.ones(self.full_shape, dtype=self.image_dtype) * self.getdarkcurrent()
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
        has_bunit = False
        if 'BUNIT' in self._image['sci',1].header :
            has_bunit = True

        countrate = False
        if (self._image[0].header['UNITCORR'].strip() == 'PERFORM') or \
            (has_bunit and self._image['sci',1].header['bunit'].find('/') != -1) :
            countrate = True

        return countrate


class NIC1InputImage(NICMOSInputImage):
    def __init__(self, filename=None):
        super().__init__(filename)
        self._effGain = 1. #get the gain from the detector subclass
        self._detector = self._image["PRIMARY"].header["CAMERA"]
        self.proc_unit = "native"

    def _getDarkRate(self):
        _darkrate = 0.08 #electrons/s
        if self.proc_unit == 'native':
            _darkrate = _darkrate / self._effGain # DN/s

        return _darkrate

    def _getDefaultReadnoise(self):
        """ This could be updated to calculate the readnoise from the NOISFILE.
        """
        _rdnoise = 26.0 # electrons
        if self.proc_unit == 'native':
            _rdnoise = _rdnoise / self._effGain # ADU

        return _rdnoise

    def setInstrumentParameters(self, instrpars):
        """ This method overrides the superclass to set default values into
            the parameter dictionary, in case empty entries are provided.
        """
        pri_header = self._image[0].header
        self.proc_unit = instrpars['proc_unit']

        if self._isNotValid (instrpars['gain'], instrpars['gnkeyword']):
            instrpars['gnkeyword'] = 'ADCGAIN' #gain has been hardcoded below

        if self._isNotValid (instrpars['rdnoise'], instrpars['rnkeyword']):
            instrpars['rnkeyword'] = None

        if self._isNotValid (instrpars['exptime'], instrpars['expkeyword']):
            instrpars['expkeyword'] = 'EXPTIME'

        for chip in self.returnAllChips(extname=self.scienceExt):
            chip._gain= 5.4 #measured gain
            chip._rdnoise   = self.getInstrParameter(instrpars['rdnoise'], pri_header,
                                                     instrpars['rnkeyword'])
            chip._exptime   = self.getInstrParameter(instrpars['exptime'], pri_header,
                                                     instrpars['expkeyword'])

            if chip._gain is None or self._exptime is None:
                print('ERROR: invalid instrument task parameter')
                raise ValueError

            # We need to treat Read Noise as a special case since it is
            # not populated in the NICMOS primary header
            if chip._rdnoise is None:
                chip._rdnoise = self._getDefaultReadnoise()

            chip._darkrate=self._getDarkRate()
            chip.darkcurrent = self.getdarkcurrent()

            chip._effGain = chip._gain
            self._assignSignature(chip._chip) #this is used in the static mask, static mask name also defined here, must be done after outputNames

        # Convert the science data to electrons if specified by the user.
        self.doUnitConversions()


class NIC2InputImage(NICMOSInputImage):
    def __init__(self,filename=None):
        super().__init__(filename)
        self._effGain=1. #measured
        self._detector=self._image["PRIMARY"].header["CAMERA"]
        self.proc_unit = "native"

    def _getDarkRate(self):
        _darkrate = 0.08 #electrons/s
        if self.proc_unit == 'native':
            _darkrate = _darkrate / self._effGain # DN/s

        return _darkrate

    def _getDefaultReadnoise(self):
        _rdnoise = 26.0 #electrons
        if self.proc_unit == 'native':
            _rdnoise = _rdnoise/self._effGain #ADU

        return _rdnoise

    def setInstrumentParameters(self, instrpars):
        """ This method overrides the superclass to set default values into
            the parameter dictionary, in case empty entries are provided.
        """
        pri_header = self._image[0].header
        self.proc_unit = instrpars['proc_unit']

        if self._isNotValid (instrpars['gain'], instrpars['gnkeyword']):
            instrpars['gnkeyword'] = 'ADCGAIN' #gain has been hardcoded below

        if self._isNotValid (instrpars['rdnoise'], instrpars['rnkeyword']):
            instrpars['rnkeyword'] = None

        if self._isNotValid (instrpars['exptime'], instrpars['expkeyword']):
            instrpars['expkeyword'] = 'EXPTIME'

        for chip in self.returnAllChips(extname=self.scienceExt):
            chip._gain= 5.4 #measured gain
            chip._rdnoise   = self.getInstrParameter(
                instrpars['rdnoise'], pri_header, instrpars['rnkeyword']
            )
            chip._exptime   = self.getInstrParameter(
                instrpars['exptime'], pri_header, instrpars['expkeyword']
            )

            if chip._gain is None or self._exptime is None:
                print('ERROR: invalid instrument task parameter')
                raise ValueError

            # We need to treat Read Noise as a special case since it is
            # not populated in the NICMOS primary header
            if chip._rdnoise is None:
                chip._rdnoise = self._getDefaultReadnoise()

            chip._darkrate=self._getDarkRate()
            chip.darkcurrent = self.getdarkcurrent()

            chip._effGain = chip._gain
            # this is used in the static mask, static mask name also defined
            # here, must be done after outputNames
            self._assignSignature(chip._chip)

        # Convert the science data to electrons if specified by the user.
        self.doUnitConversions()

    def createHoleMask(self):
        """Add in a mask for the coronographic hole to the general static
        pixel mask. """
        pass


class NIC3InputImage(NICMOSInputImage):
    def __init__(self, filename=None):
        super().__init__(filename)
        self._detector=self._image["PRIMARY"].header["CAMERA"] #returns 1,2,3
        self._effGain = 1.
        self.proc_unit = "native"

    def _getDarkRate(self):
        _darkrate = 0.15 #electrons/s
        if self.proc_unit == 'native':
            _darkrate = _darkrate/self._effGain #DN/s

        return _darkrate

    def _getDefaultReadnoise(self):
        _rdnoise = 29.0 # electrons
        if self.proc_unit == 'native':
            _rdnoise = _rdnoise/self._effGain #ADU

        return _rdnoise

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

        for chip in self.returnAllChips(extname=self.scienceExt):
            chip._gain= 6.5 #measured gain
            chip._rdnoise   = self.getInstrParameter(
                instrpars['rdnoise'], pri_header, instrpars['rnkeyword']
            )
            chip._exptime   = self.getInstrParameter(
                instrpars['exptime'], pri_header, instrpars['expkeyword']
            )

            if chip._gain is None or self._exptime is None:
                print('ERROR: invalid instrument task parameter')
                raise ValueError

            # We need to treat Read Noise as a special case since it is
            # not populated in the NICMOS primary header
            if chip._rdnoise is None:
                chip._rdnoise = self._getDefaultReadnoise()

            chip._darkrate=self._getDarkRate()
            chip.darkcurrent = self.getdarkcurrent()

            chip._effGain = chip._gain
            self._assignSignature(chip._chip) #this is used in the static mask, static mask name also defined here, must be done after outputNames

        # Convert the science data to electrons if specified by the user.
        self.doUnitConversions()
