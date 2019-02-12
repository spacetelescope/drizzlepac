"""
`wfpc2Data` module provides classes used to import WFPC2 specific instrument data.

:Authors: Warren Hack, Ivo Busko, Christopher Hanley

:License: :doc:`LICENSE`

"""
import os

from astropy.io import fits
import numpy as np

from stsci.tools import fileutil, readgeis

from .imageObject import imageObject
from . import buildmask

# Translation table for any image that does not use the DQ extension of the MEF
# for the DQ array.
DQ_EXTNS = {'c0h': 'sdq', 'c0f': 'sci', 'c0m': 'sci'}

#### Calibrated gain and readnoise values for each chip
WFPC2_GAINS = {
    1: {7:[7.12,5.24],15:[13.99,7.02]},
    2: {7:[7.12,5.51],15:[14.50,7.84]},
    3: {7:[6.90,5.22],15:[13.95,6.99]},
    4: {7:[7.10,5.19],15:[13.95,8.32]}
}
WFPC2_DETECTOR_NAMES = {1: "PC", 2: "WF2", 3: "WF3", 4: "WF4"}

class WFPC2InputImage (imageObject):

    SEPARATOR = '_'

    def __init__(self, filename, group=None):
        super().__init__(filename, group=group)
        # define the cosmic ray bits value to use in the dq array
        self.cr_bits_value = 4096
        self._instrument=self._image["PRIMARY"].header["INSTRUME"]
        self._effGain = 1
        self.errExt = None

        # Attribute defining the pixel dimensions of WFPC2 chips.
        self.full_shape = (800, 800)
        self.native_units = "COUNTS"
        self.flatkey = 'FLATFILE'

        # Reference Plate Scale used for updates to MDRIZSKY, we should get this from the wcs class
        #self.refplatescale = 0.0996 # arcsec / pixel
        for chip in range(1,self._numchips+1,1):
            self._assignSignature(chip) #this is used in the static mask
            self._image[self.scienceExt,chip].cte_dir = -1 # independent of amp, chip
            det=int(self._image[self.scienceExt,chip].header["DETECTOR"])
            self._image[self.scienceExt,chip]._detector=WFPC2_DETECTOR_NAMES[det]
            self._image[self.scienceExt,chip].darkcurrent = self.getdarkcurrent(chip)

    def find_DQ_extension(self):
        """ Return the suffix for the data quality extension and the name of
        the file which that DQ extension should be read from.

        """
        dqfile = None
        # Look for additional file with DQ array, primarily for WFPC2 data

        indx = self._filename.find('.fits')

        if indx > 3:
            suffix = self._filename[indx-4:indx]
            dqfile = self._filename.replace(suffix[:3],'_c1')

        elif indx < 0 and len(self._filename) > 3 and \
             self._filename[-4] == os.extsep and \
             self._filename[-1].lower() == 'h':

            # assume we've got a GEIS file
            dqfile  = self._filename[:-2]+'1'+self._filename[-1]
            hdulist = readgeis.readgeis(dqfile)
            prih    = hdulist[0].header
            if 'FILETYPE' in prih:
                dq_suffix = prih['FILETYPE'].strip().upper()
            else:
                # assume extension name is 'SDQ' for WFPC2 GEIS files
                dq_suffix = 'SDQ'
            hdulist.close()
            return dqfile,dq_suffix

        else:
            raise ValueError("Input file {} does not appear to be neither " \
                        "a FITS file nor a GEIS file.".format(self._filename))

        if os.path.exists(dqfile):
            dq_suffix = fits.getval(dqfile, "EXTNAME", ext=1, memmap=False)
        else:
            dq_suffix = "SCI"

        return dqfile, dq_suffix

    def getEffGain(self):
        """
        Method used to return the effective gain of a instrument's
        detector.

        Returns
        -------
        gain : float
            The effective gain.
        """
        return self._effGain

    def setInstrumentParameters(self, instrpars):
        """ This method overrides the superclass to set default values into
            the parameter dictionary, in case empty entries are provided.
        """
        pri_header = self._image[0].header
        self.proc_unit = instrpars['proc_unit']
        instrpars['gnkeyword'] = 'ATODGAIN'  # hard-code for WFPC2 data
        instrpars['rnkeyword'] = None

        if self._isNotValid (instrpars['exptime'], instrpars['expkeyword']):
            instrpars['expkeyword'] = 'EXPTIME'

        for chip in self.returnAllChips(extname=self.scienceExt):
            chip._headergain = self.getInstrParameter(
                instrpars['gain'], pri_header, instrpars['gnkeyword']
            )
            chip._exptime = self.getInstrParameter(
                instrpars['exptime'], pri_header, instrpars['expkeyword']
            )
            # We need to treat Read Noise as a special case since it is
            # not populated in the WFPC2 primary header
            if instrpars['rnkeyword'] is None:
                chip._rdnoise = None
            else:
                chip._rdnoise = self.getInstrParameter(
                    instrpars['rdnoise'], pri_header, instrpars['rnkeyword']
                )

            if chip._headergain is None or chip._exptime is None:
                print('ERROR: invalid instrument task parameter')
                raise ValueError

        # We need to determine if the user has used the default readnoise/gain value
        # since if not, they will need to supply a gain/readnoise value as well
        usingDefaultGain = instrpars['gnkeyword'] == 'ATODGAIN'
        usingDefaultReadnoise = instrpars['rnkeyword'] in [None, 'None']

        # If the user has specified either the readnoise or the gain, we need to make sure
        # that they have actually specified both values.  In the default case, the readnoise
        # of the system depends on what the gain

        if usingDefaultReadnoise and usingDefaultGain:
            self._setchippars()
        elif usingDefaultReadnoise and not usingDefaultGain:
            raise ValueError("ERROR: You need to supply readnoise information\n when not using the default gain for WFPC2.")
        elif not usingDefaultReadnoise and usingDefaultGain:
            raise ValueError("ERROR: You need to supply gain information when\n not using the default readnoise for WFPC2.")
        else:
            # In this case, the user has specified both a gain and readnoise values.  Just use them as is.
            for chip in self.returnAllChips(extname=self.scienceExt):
                chip._gain = chip._headergain
            print("Using user defined values for gain and readnoise")

        # Convert the science data to electrons
        self.doUnitConversions()

    def getflat(self,chip):
        """
        Method for retrieving a detector's flat field.

        Returns
        -------
        flat : array
            The flat-field array in the same shape as the input image.

        """
        # For the WFPC2 flat we need to invert
        # for use in Multidrizzle
        flat = 1.0 / super().getflat(chip)
        return flat

    def doUnitConversions(self):
        """ Apply unit conversions to all the chips, ignoring the group parameter.
            This insures that all the chips get the same conversions when this
            gets done, even if only 1 chip was specified to be processed.
        """
         # Image information
        _handle = fileutil.openImage(self._filename, mode='readonly', memmap=False)

        # Now convert the SCI array(s) units
        for det in range(1,self._numchips+1):
            chip=self._image[self.scienceExt,det]
            conversionFactor = 1.0
            # add D2IMFILE to outputNames for removal by 'clean()' method later
            if 'D2IMFILE' in _handle[0].header and _handle[0].header['D2IMFILE'] not in ["","N/A"]:
                chip.outputNames['d2imfile'] = _handle[0].header['D2IMFILE']

            if chip._gain is not None:
                """
                # Multiply the values of the sci extension pixels by the gain.
                print "Converting %s[%d] from COUNTS to ELECTRONS"%(self._filename,det)

                # If the exptime is 0 the science image will be zeroed out.
                np.multiply(_handle[self.scienceExt,det].data,chip._gain,_handle[self.scienceExt,det].data)
                chip.data=_handle[self.scienceExt,det].data

                # Set the BUNIT keyword to 'electrons'
                chip._bunit = 'ELECTRONS'
                chip.header.update('BUNIT','ELECTRONS')
                _handle[self.scienceExt,det].header.update('BUNIT','ELECTRONS')

                # Update the PHOTFLAM value
                photflam = _handle[self.scienceExt,det].header['PHOTFLAM']
                _handle[self.scienceExt,det].header.update('PHOTFLAM',(photflam/chip._gain))
                """
                conversionFactor = chip._gain
                chip._effGain = chip._gain #1.
                chip._conversionFactor = conversionFactor #1.

            else:
                msg = "Invalid gain value for data, no conversion done"
                print(msg)
                raise ValueError(msg)

        # Close the files and clean-up
        _handle.close()
        self._effGain = conversionFactor # 1.

    def getdarkcurrent(self,exten):
        """
        Return the dark current for the WFPC2 detector.  This value
        will be contained within an instrument specific keyword.
        The value in the image header will be converted to units
        of electrons.

        Returns
        -------
        darkcurrent : float
            Dark current for the WFPC3 detector in **units of counts/electrons**.

        """
        darkrate = 0.005 # electrons / s
        if self.proc_unit == 'native':
            darkrate = darkrate / self.getGain(exten) #count/s

        try:
            chip = self._image[0]
            darkcurrent = chip.header['DARKTIME'] * darkrate

        except:
            msg =  "#############################################\n"
            msg += "#                                           #\n"
            msg += "# Error:                                    #\n"
            msg += "#   Cannot find the value for 'DARKTIME'    #\n"
            msg += "#   in the image header.  WFPC2 input       #\n"
            msg += "#   images are expected to have this header #\n"
            msg += "#   keyword.                                #\n"
            msg += "#                                           #\n"
            msg += "# Error occured in the WFPC2InputImage class#\n"
            msg += "#                                           #\n"
            msg += "#############################################\n"
            raise ValueError(msg)

        return darkcurrent

    def getReadNoise(self, exten):
        """
        Method for returning the readnoise of a detector (in counts).

        Returns
        -------
        readnoise : float
            The readnoise of the detector in **units of counts/electrons**.

        """
        rn = self._image[exten]._rdnoise
        if self.proc_unit == 'native':
            rn = self._rdnoise / self.getGain(exten)
        return rn

    def buildMask(self, chip, bits=0, write=False):
        """ Build masks as specified in the user parameters found in the
            configObj object.
        """
        sci_chip = self._image[self.scienceExt,chip]
        ### For WFPC2 Data, build mask files using:
        maskname = sci_chip.dqrootname+'_dqmask.fits'
        dqmask_name = buildmask.buildShadowMaskImage(sci_chip.dqfile,sci_chip.detnum,sci_chip.extnum,maskname,bitvalue=bits,binned=sci_chip.binned)
        sci_chip.dqmaskname = dqmask_name
        sci_chip.outputNames['dqmask'] = dqmask_name
        sci_chip.outputNames['tmpmask'] = 'wfpc2_inmask%d.fits'%(sci_chip.detnum)
        dqmask = fits.getdata(dqmask_name, ext=0, memmap=False)
        return dqmask

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
        sig = (instr + WFPC2_DETECTOR_NAMES[detnum], (nx, ny), chip) #signature is a tuple
        sci_chip.signature=sig #signature is a tuple

    def _setchippars(self):
        for chip in self.returnAllChips(extname=self.scienceExt):
            try:
                chip._gain,chip._rdnoise = WFPC2_GAINS[chip.detnum][chip._headergain]
            except KeyError:
                raise ValueError("! Header gain value is not valid for WFPC2")
