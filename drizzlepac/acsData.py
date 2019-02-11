"""
Class used to model ACS specific instrument data.

:Authors: Christopher Hanley, Warren Hack, Ivo Busko, David Grumm

:License: :doc:`LICENSE`

"""
from stsci.tools import fileutil
import numpy as np
from .imageObject import imageObject


class ACSInputImage(imageObject):

    SEPARATOR = '_'

    def __init__(self,filename=None,group=None):
        super().__init__(filename, group=group)
        # define the cosmic ray bits value to use in the dq array
        self.cr_bits_value = 4096
        self._instrument=self._image["PRIMARY"].header["INSTRUME"]
        self._effGain=1.
        self.flatkey = 'PFLTFILE'

        for chip in range(1,self._numchips+1,1):
            if self._image[self.scienceExt,chip].group_member:
                self._image[self.scienceExt,chip].darkcurrent=self.getdarkcurrent(chip)

    def doUnitConversions(self):
        # Effective gain to be used in the driz_cr step.  Since the
        # ACS images have already been converted to electrons,
        # the effective gain is 1.
        for chip in self.returnAllChips(extname=self.scienceExt):
            chip._effGain = 1.0 #chip._effGain is was drizCr uses
            chip._conversionFactor = 1.0

        self._effGain=1.0

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

        sig = (instr + self._detector, (nx, ny), int(chip))  # signature is a tuple
        sci_chip.signature=sig  # signature is a tuple

    def getdarkcurrent(self,extver):
        """
        Return the dark current for the ACS detector.  This value
        will be contained within an instrument specific keyword.
        The value in the image header will be converted to units
        of electrons.

        Returns
        -------
        darkcurrent: float
            Dark current value for the ACS detector in **units of electrons**.

        """
        darkcurrent=0.
        try:
            darkcurrent = self._image[self.scienceExt,extver].header['MEANDARK']
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
            raise ValueError(str)

        return darkcurrent


class WFCInputImage(ACSInputImage):
    def __init__(self, filename=None, group=None):
        super().__init__(filename, group=group)
        self.full_shape = (4096, 2048)
        self._detector=self._image["PRIMARY"].header["DETECTOR"]

        # get cte direction, which depends on which chip but is independent of amp
        for chip in range(1, self._numchips + 1, 1):
            self._assignSignature(chip) #this is used in the static mask
            if chip == 1:
                self._image[self.scienceExt,chip].cte_dir = -1
            elif chip == 2:
                self._image[self.scienceExt,chip].cte_dir = 1

    def setInstrumentParameters(self,instrpars):
        """ This method overrides the superclass to set default values into
            the parameter dictionary, in case empty entries are provided.

            This method gets called from processInput.
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

        self.proc_unit = instrpars['proc_unit']

        if self._isNotValid (instrpars['gain'], instrpars['gnkeyword']):
            instrpars['gnkeyword'] = 'ATODGNA,ATODGNB,ATODGNC,ATODGND'
        if self._isNotValid (instrpars['rdnoise'], instrpars['rnkeyword']):
            instrpars['rnkeyword'] = 'READNSEA,READNSEB,READNSEC,READNSED'
        if self._isNotValid (instrpars['exptime'], instrpars['expkeyword']):
            instrpars['expkeyword'] = 'EXPTIME'

        for chip in self.returnAllChips(extname=self.scienceExt):
            chip._gain      = self.getInstrParameter(instrpars['gain'], pri_header,
                                                     instrpars['gnkeyword'])
            chip._rdnoise   = self.getInstrParameter(instrpars['rdnoise'], pri_header,
                                                     instrpars['rnkeyword'])
            chip._exptime   = self.getInstrParameter(instrpars['exptime'], pri_header,
                                                     instrpars['expkeyword'])
            chip._effGain = 1.

            if (chip._gain is None or chip._rdnoise is None or
                chip._exptime is None):
                print('ERROR: invalid instrument task parameter')
                raise ValueError

        # Convert the science data to electrons if specified by the user.
        self.doUnitConversions()


class HRCInputImage (ACSInputImage):
    def __init__(self, filename=None, group=None):
        super().__init__(filename, group=group)
        self._detector = self._image['PRIMARY'].header["DETECTOR"]
        self.full_shape = (1024, 1024)
        amp = self._image['PRIMARY'].header["CCDAMP"]
        self._detector = self._image['PRIMARY'].header["DETECTOR"]

        for chip in range(1,self._numchips+1,1):
            self._assignSignature(chip) #this is used in the static mask

            # cte direction depends on amp (but is independent of chip):
            if amp in ['A', 'B']:
                self._image[self.scienceExt, chip].cte_dir = 1
            elif amp in ['C', 'D']:
                self._image[self.scienceExt, chip].cte_dir = -1

    def setInstrumentParameters(self,instrpars):
        """ This method overrides the superclass to set default values into
            the parameter dictionary, in case empty entries are provided.

            This method gets called from processInput.
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

        self.proc_unit = instrpars['proc_unit']

        if self._isNotValid (instrpars['gain'], instrpars['gnkeyword']):
            instrpars['gnkeyword'] = 'ATODGNA,ATODGNB,ATODGNC,ATODGND'
        if self._isNotValid (instrpars['rdnoise'], instrpars['rnkeyword']):
            instrpars['rnkeyword'] = 'READNSEA,READNSEB,READNSEC,READNSED'
        if self._isNotValid (instrpars['exptime'], instrpars['expkeyword']):
            instrpars['expkeyword'] = 'EXPTIME'

        for chip in self.returnAllChips(extname=self.scienceExt):
            chip._gain      = self.getInstrParameter(instrpars['gain'], pri_header,
                                                     instrpars['gnkeyword'])
            chip._rdnoise   = self.getInstrParameter(instrpars['rdnoise'], pri_header,
                                                     instrpars['rnkeyword'])
            chip._exptime   = self.getInstrParameter(instrpars['exptime'], pri_header,
                                                     instrpars['expkeyword'])
            chip._effGain = chip._gain

            if (chip._gain is None or chip._rdnoise is None or
                chip._exptime is None):
                print('ERROR: invalid instrument task parameter')
                raise ValueError

        # Convert the science data to electrons if specified by the user.
        self.doUnitConversions()


class SBCInputImage (ACSInputImage):
    def __init__(self, filename=None, group=None):
        super().__init__(filename, group=group)
        self.full_shape = (1024, 1024)
        self._detector = self._image['PRIMARY'].header["DETECTOR"]

        # no cte correction for SBC so set cte_dir=0.
        print('WARNING: No cte correction will be made for this SBC data.')
        for chip in range(1,self._numchips+1,1):
            self._assignSignature(chip) #this is used in the static mask
            self._image[self.scienceExt,chip].cte_dir = 0

    def _setSBCchippars(self):
        self._setDefaultSBCGain()
        self._setDefaultSBCReadnoise()

    def _setDefaultSBCGain(self):
        self._gain = 1

    def _setDefaultSBCReadnoise(self):
        self._rdnoise = 0

    def setInstrumentParameters(self,instrpars):
        """ Sets the instrument parameters.
        """
        pri_header = self._image[0].header

        if self._isNotValid (instrpars['gain'], instrpars['gnkeyword']):
            instrpars['gnkeyword'] = None
        if self._isNotValid (instrpars['rdnoise'], instrpars['rnkeyword']):
            instrpars['rnkeyword'] = None
        if self._isNotValid (instrpars['exptime'], instrpars['expkeyword']):
            instrpars['expkeyword'] = 'EXPTIME'

        # We need to treat Read Noise and Gain as a special case since it is
        # not populated in the SBC primary header for the MAMA
        for chip in self.returnAllChips(extname=self.scienceExt):
            chip._gain      = 1.0 #self.getInstrParameter("", pri_header,
                                  #                  instrpars['gnkeyword'])
            chip._rdnoise   = 0.0 #self.getInstrParameter("", pri_header,
                                  #                   instrpars['rnkeyword'])
            chip._exptime   = self.getInstrParameter(instrpars['exptime'], pri_header,
                                                     instrpars['expkeyword'])
            if chip._exptime is None:
                print('ERROR: invalid instrument task parameter')
                raise ValueError

        # We need to determine if the user has used the default readnoise/gain value
        # since if not, they will need to supply a gain/readnoise value as well
        usingDefaultGain = instrpars['gnkeyword'] is None
        usingDefaultReadnoise = instrpars['rnkeyword'] is None

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
