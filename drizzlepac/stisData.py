"""
``stisData`` module provides classes used to import STIS specific instrument
data.

:Authors: Megan Sosey, Christopher Hanley, Mihai Cara

:License: :doc:`LICENSE`

"""
from stsci.tools import fileutil
import numpy as np
from .imageObject import imageObject


class STISInputImage (imageObject):

    SEPARATOR = '_'

    def __init__(self,filename=None, group=None):
        super().__init__(filename, group=group)

        # define the cosmic ray bits value to use in the dq array
        self.cr_bits_value = 8192
        self._effGain = 1.
        self._instrument = self._image["PRIMARY"].header["INSTRUME"]  # this just shows instrument, not detector, detector asigned by subclass
        self.native_units = 'COUNTS'

    def getflat(self, chip):
        """
        Method for retrieving a detector's flat field. For STIS there are three.
        This method will return an array the same shape as the image.

        """
        sci_chip = self._image[self.scienceExt,chip]
        exten = self.errExt+','+str(chip)

        # The keyword for STIS flat fields in the primary header of the flt

        lflatfile = fileutil.osfn(self._image["PRIMARY"].header['LFLTFILE'])
        pflatfile = fileutil.osfn(self._image["PRIMARY"].header['PFLTFILE'])

        # Try to open the file in the location specified by LFLTFILE.
        try:
            handle = fileutil.openImage(lflatfile, mode='readonly', memmap=False)
            hdu = fileutil.getExtn(handle,extn=exten)
            lfltdata = hdu.data
            if lfltdata.shape != self.full_shape:
                lfltdata = expand_image(lfltdata, self.full_shape)
        except IOError:
            lfltdata = np.ones(self.full_shape, dtype=sci_chip.data.dtype)
            print("Cannot find file '{:s}'. Treating flatfield constant value "
                  "of '1'.\n".format(lflatfile))

        # Try to open the file in the location specified by PFLTFILE.
        try:
            handle = fileutil.openImage(pflatfile, mode='readonly', memmap=False)
            hdu = fileutil.getExtn(handle,extn=exten)
            pfltdata = hdu.data
        except IOError:
            pfltdata = np.ones(self.full_shape, dtype=sci_chip.data.dtype)
            print("Cannot find file '{:s}'. Treating flatfield constant value "
                  "of '1'.\n".format(pflatfile))

        flat = lfltdata * pfltdata

        return flat

    def doUnitConversions(self):
        """Convert the data to electrons.

        This converts all science data extensions and saves
        the results back to disk. We need to make sure
        the data inside the chips already in memory is altered as well.

        """
         # Image information
        _handle = fileutil.openImage(self._filename, mode='readonly', memmap=False)

        for det in range(1,self._numchips+1,1):
            chip=self._image[self.scienceExt,det]
            if chip._gain is not None:

                conversionFactor = chip._gain
                chip._effGain = chip._gain #1.
                chip._conversionFactor = conversionFactor #1.

            else:
                msg = "Invalid gain value for data, no conversion done"
                print(msg)
                raise ValueError(msg)

        # Close the files and clean-up
        _handle.close()

        self._effGain = conversionFactor # 1.0

    def _assignSignature(self, chip):
        """Assign a unique signature for the image based
           on the  instrument, detector, chip, and size
           this will be used to uniquely identify the appropriate
           static mask for the image.

           This also records the filename for the static mask to the outputNames dictionary.

        """
        sci_chip = self._image[self.scienceExt,chip]
        ny=sci_chip._naxis1
        nx=sci_chip._naxis2
        detnum = sci_chip.detnum
        instr=self._instrument

        sig=(instr+self._detector,(nx,ny),int(detnum)) #signature is a tuple
        sci_chip.signature=sig #signature is a tuple


class CCDInputImage(STISInputImage):
    def __init__(self, filename=None, group=None):
        super().__init__(filename, group=group)

        self.full_shape = (1024, 1024)
        self._detector=self._image["PRIMARY"].header["DETECTOR"]

        for chip in range(1,self._numchips+1,1):
            self._image[self.scienceExt,chip].cte_dir = 1
            self._image[self.scienceExt,chip].darkcurrent = self.getdarkcurrent()

        self.cte_dir =  1

    def getdarkcurrent(self):
        """
        Returns the dark current for the STIS CCD chip.

        Returns
        -------
        darkcurrent : float
            Dark current value in **units of electrons** (or counts, if proc_unit=='native').
        """
        darkcurrent = 0.009 #electrons/sec
        if self.proc_unit == 'native':
            return darkcurrent / self._gain()
        return darkcurrent

    def getReadNoise(self):
        """
        Method for returning the readnoise of a detector (in DN).

        :units: DN

        This should work on a chip, since different chips to be consistant with other
        detector classes where different chips have different gains.

        """
        if self.proc_unit == 'native':
            return self._rdnoise / self._gain()
        return self._rdnoise

    def setInstrumentParameters(self, instrpars):
        """ This method overrides the superclass to set default values into
            the parameter dictionary, in case empty entries are provided.
        """
        pri_header = self._image[0].header

        if self._isNotValid (instrpars['gain'], instrpars['gnkeyword']):
            instrpars['gnkeyword'] = 'ATODGAIN'
        if self._isNotValid (instrpars['rdnoise'], instrpars['rnkeyword']):
            instrpars['rnkeyword'] = 'READNSE'
        if self._isNotValid (instrpars['exptime'], instrpars['expkeyword']):
            instrpars['expkeyword'] = 'EXPTIME'

        for chip in self.returnAllChips(extname=self.scienceExt):

            chip._gain      = self.getInstrParameter(instrpars['gain'], pri_header,
                                                     instrpars['gnkeyword'])
            chip._rdnoise   = self.getInstrParameter(instrpars['rdnoise'], pri_header,
                                                     instrpars['rnkeyword'])
            chip._exptime   = self.getInstrParameter(instrpars['exptime'], chip.header,
                                                     instrpars['expkeyword'])

            if chip._gain is None or chip._rdnoise is None or chip._exptime is None:
                print('ERROR: invalid instrument task parameter')
                raise ValueError

            chip._effGain = chip._gain
            self._assignSignature(chip._chip) #this is used in the static mask

        self.doUnitConversions()


class NUVInputImage(STISInputImage):
    def __init__(self, filename, group=None):
        self.effGain = 1.0

        super().__init__(filename, group=None)

        self._detector=self._image["PRIMARY"].header["DETECTOR"]

        # no cte correction for STIS/NUV-MAMA so set cte_dir=0.
        print('WARNING: No cte correction will be made for this STIS/NUV-MAMA data.')

        for chip in range(1,self._numchips+1,1):
            self._image[self.scienceExt,chip].cte_dir = 0
            self._image[self.scienceExt,chip].darkcurrent = self.getdarkcurrent()

    def setInstrumentParameters(self, instrpars):
        """ This method overrides the superclass to set default values into
            the parameter dictionary, in case empty entries are provided.
        """
        pri_header = self._image[0].header

        if self._isNotValid (instrpars['gain'], instrpars['gnkeyword']):
            instrpars['gnkeyword'] = None
        if self._isNotValid (instrpars['rdnoise'], instrpars['rnkeyword']):
            instrpars['rnkeyword'] = None
        if self._isNotValid (instrpars['exptime'], instrpars['expkeyword']):
            instrpars['expkeyword'] = 'EXPTIME'

        # We need to determine if the user has used the default readnoise/gain value
        # since if not, they will need to supply a gain/readnoise value as well
        usingDefaultGain = instrpars['gnkeyword'] is None
        usingDefaultReadnoise = instrpars['rnkeyword'] is None

        for chip in self.returnAllChips(extname=self.scienceExt):
            #pri_header=chip.header
            chip.cte_dir=0
            # We need to treat Read Noise and Gain as a special case since it is
            # not populated in the STIS primary header for the MAMAs
            if instrpars['rnkeyword'] is not None:
                chip._rdnoise   = self.getInstrParameter(
                    instrpars['rdnoise'], pri_header, instrpars['rnkeyword']
                )
            else:
                chip._rdnoise = None

            if instrpars['gnkeyword'] is not None:
                chip._gain = self.getInstrParameter(
                    instrpars['gain'], pri_header, instrpars['gnkeyword']
                )
            else:
                chip._gain = None

            # Set the default readnoise or gain values based upon the amount of user input given.

            if usingDefaultReadnoise:
                chip._rdnoise= self._setMAMADefaultReadnoise()

            if usingDefaultGain:
                chip._gain = self._setMAMADefaultGain()

            self._assignSignature(chip._chip) #this is used in the static mask
            chip._exptime   = self.getInstrParameter(instrpars['exptime'], chip.header,
                                                     instrpars['expkeyword'])

            if chip._exptime is None:
                print('ERROR: invalid instrument task parameter')
                raise ValueError
        # Convert the science data to electrons if specified by the user.
        self.doUnitConversions()

    def _setMAMAchippars(self):
        self._setMAMADefaultGain()
        self._setMAMADefaultReadnoise()

    def _setMAMADefaultGain(self):
        self._gain = 1
        self.effGain = 1
        return self._gain

    def _setMAMADefaultReadnoise(self):
        self._rdnoise = 0
        return self._rdnoise

    def getdarkcurrent(self):
        """
        Returns the dark current for the STIS NUV detector.

        Returns
        -------
        darkcurrent : float
            Dark current value in **units of electrons** (or counts, if proc_unit=='native').
        """
        darkcurrent = 0.0013 #electrons/sec
        if self.proc_unit == 'native':
            return darkcurrent / self._gain()
        return darkcurrent

    def doUnitConversions(self):
        """Convert the data to electrons.

        This converts all science data extensions and saves
        the results back to disk. We need to make sure
        the data inside the chips already in memory is altered as well.

        """
        for det in range(1,self._numchips+1,1):

            chip=self._image[self.scienceExt,det]

            conversionFactor = self.effGain
            chip._gain = self.effGain #1.
            chip.effGain = self.effGain
            chip._conversionFactor = conversionFactor #1.


class FUVInputImage(STISInputImage):
    def __init__(self,filename=None,group=None):
        self.effGain=1.0
        super().__init__(filename, group=group)
        self._detector=self._image["PRIMARY"].header["DETECTOR"]

        # no cte correction for STIS/FUV-MAMA so set cte_dir=0.
        print('WARNING: No cte correction will be made for this STIS/FUV-MAMA data.')
        for chip in range(1,self._numchips+1,1):
            self._image[self.scienceExt,chip].cte_dir = 0
            self._image[self.scienceExt,chip].darkcurrent = self.getdarkcurrent()

    def setInstrumentParameters(self, instrpars):
        """ This method overrides the superclass to set default values into
            the parameter dictionary, in case empty entries are provided.
        """
        pri_header = self._image[0].header
        usingDefaultGain = False
        usingDefaultReadnoise = False

        if self._isNotValid (instrpars['gain'], instrpars['gnkeyword']):
            instrpars['gnkeyword'] = None
        if self._isNotValid (instrpars['rdnoise'], instrpars['rnkeyword']):
            instrpars['rnkeyword'] = None
        if self._isNotValid (instrpars['exptime'], instrpars['expkeyword']):
            instrpars['expkeyword'] = 'EXPTIME'

        for chip in self.returnAllChips(extname=self.scienceExt):
            #pri_header=chip.header #stis stores stuff in the science data header

            chip.cte_dir=0

            chip._exptime = self.getInstrParameter(
                instrpars['exptime'], chip.header, instrpars['expkeyword']
            )

            if chip._exptime is None:
                print('ERROR: invalid instrument task parameter')
                raise ValueError

            if instrpars['rnkeyword'] is not None:
                chip._rdnoise   = self.getInstrParameter(
                    instrpars['rdnoise'], pri_header, instrpars['rnkeyword']
                )
            else:
                chip._rdnoise = None
                usingDefaultReadnoise = True

            if instrpars['gnkeyword'] is not None:
                chip._gain = self.getInstrParameter(
                    instrpars['gain'], pri_header, instrpars['gnkeyword']
                )
            else:
                chip._gain = None
                usingDefaultGain = True

            if chip._exptime is None:
                print('ERROR: invalid instrument task parameter')
                raise ValueError

            # We need to determine if the user has used the default readnoise/gain value
            # since if not, they will need to supply a gain/readnoise value as well

            if usingDefaultReadnoise:
                chip._rdnoise= self._setMAMADefaultReadnoise()

            if usingDefaultGain:
                chip._gain = self._setMAMADefaultGain()

            self._assignSignature(chip._chip) #this is used in the static mask
            chip._effGain=chip._gain

        # Convert the science data to electrons if specified by the user.
        self.doUnitConversions()

    def getdarkcurrent(self):
        """
        Returns the dark current for the STIS FUV detector.

        Returns
        -------
        darkcurrent : float
            Dark current value in **units of electrons** (or counts, if proc_unit=='native').
        """
        darkcurrent = 0.07 #electrons/sec
        if self.proc_unit == 'native':
            return darkcurrent / self._gain()
        return darkcurrent

    def _setMAMADefaultGain(self):
        return 1

    def _setMAMADefaultReadnoise(self):
        return 0

    def doUnitConversions(self):
        """Convert the data to electrons.

        This converts all science data extensions and saves
        the results back to disk. We need to make sure
        the data inside the chips already in memory is altered as well.

        """
        for det in range(1,self._numchips+1,1):
            chip=self._image[self.scienceExt,det]
            conversionFactor = self.effGain
            chip._gain = self.effGain #1.
            chip.effGain = self.effGain
            chip._conversionFactor = conversionFactor #1.


def expand_image(image, shape):
    """ Expand image from original shape to requested shape. Output shape
    must be an integer multiple of input image shape for each axis. """
    if (shape[0] % image.shape[0]) or (shape[1] % image.shape[1]):
        raise ValueError("Output shape must be an integer multiple of input "
                         "image shape.")
    sx = shape[1] // image.shape[1]
    sy = shape[0] // image.shape[0]

    ox = (sx - 1.0) / (2.0 * sx)
    oy = (sy - 1.0) / (2.0 * sy)

    # generate output coordinates:
    y, x = np.indices(shape, dtype=np.float)
    x = x / sx - ox
    y = y / sy - oy

    # interpolate:
    return bilinear_interp(image, x, y)


def bilinear_interp(data, x, y):
    """ Interpolate input ``data`` at "pixel" coordinates ``x`` and ``y``. """
    x = np.asarray(x)
    y = np.asarray(y)
    if x.shape != y.shape:
        raise ValueError("X- and Y-coordinates must have identical shapes.")

    out_shape = x.shape
    out_size = x.size
    x = x.ravel()
    y = y.ravel()

    x0 = np.empty(out_size, dtype=np.int)
    y0 = np.empty(out_size, dtype=np.int)
    np.clip(x, 0, data.shape[1] - 2, out=x0)
    np.clip(y, 0, data.shape[0] - 2, out=y0)
    x1 = x0 + 1
    y1 = y0 + 1

    f00 = data[(y0, x0)]
    f10 = data[(y1, x0)]
    f01 = data[(y0, x1)]
    f11 = data[(y1, x1)]

    w00 = (x1 - x) * (y1 - y)
    w10 = (x1 - x) * (y - y0)
    w01 = (x - x0) * (y1 - y)
    w11 = (x - x0) * (y - y0)

    interp = w00 * f00 + w10 * f10 + w01 * f01 + w11 * f11

    return interp.reshape(out_shape).astype(data.dtype.type)
