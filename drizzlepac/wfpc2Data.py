"""
``wfpc2Data`` module provides classes used to import WFPC2 specific instrument data.

:Authors: Warren Hack, Ivo Busko, Christopher Hanley

:License: :doc:`/LICENSE`

"""
import copy
import os
import shutil
import glob

from astropy.io import fits
import numpy as np

from stsci.tools import fileutil, readgeis

from .imageObject import imageObject
from . import buildmask

# Define default public CRDS server URL to use in case user does not specify one in ``os.environ``
PUBLIC_CRDS_SERVER_URL = "https://hst-crds.stsci.edu"

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

    flat_file_map = {}

    def __init__(self, filename, output=None, group=None):
        super().__init__(filename, output=output, group=group)
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
        _flt_file = False

        if self._filename.endswith('_flt.fits'):
            dqfile = self._filename
            _flt_file = True

        elif indx > 3:
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
            if _flt_file:
                dq_suffix = 'DQ'
            else:
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

    def getflat(self, chip, flat_file=None, flat_ext=None):
        """
        Method for retrieving a detector's flat field.

        Parameters
        ----------
        chip : int
            Chip number. Same as FITS ``EXTVER``.

        flat_file : str, None
            Flat field file name. If not specified, it will be determined
            automatically from image header.

        flat_ext : str, None
            Flat field extension name (same as FITS ``EXTNAME``). Specifies
            extension name containing flat field data.

        Returns
        -------
        flat : numpy.ndarray
            The flat-field array in the same shape as the input image.

        """
        # For the WFPC2 flat we need to invert
        # for use in Multidrizzle
        if flat_file is None:
            filename = fileutil.osfn(self._image["PRIMARY"].header[self.flatkey])
            if filename in WFPC2InputImage.flat_file_map:
                flat_file, mef_flat_ext = WFPC2InputImage.flat_file_map[filename]
            else:
                h = fileutil.openImage(filename, mode='readonly', memmap=False)
                flat_file = h.filename()
                mef_flat_ext = h[0].header.get('FILETYPE', '')
                mef_flat_ext = h[1].header.get('EXTNAME', mef_flat_ext)
                h.close()
                WFPC2InputImage.flat_file_map[filename] = (flat_file, mef_flat_ext)
            if flat_ext is None:
                flat_ext = mef_flat_ext

        elif flat_ext is None:
            h = fileutil.openImage(flat_file, mode='readonly', memmap=False,
                                   writefits=False)
            flat_ext = h[0].header.get('FILETYPE', '')
            flat_ext = h[1].header.get('EXTNAME', flat_ext)
            h.close()

        flat = 1.0 / super().getflat(chip, flat_file, flat_ext)
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
        sig = (self.outroot, (nx, ny), chip) #signature is a tuple
        sci_chip.signature=sig #signature is a tuple

    def _setchippars(self):
        for chip in self.returnAllChips(extname=self.scienceExt):
            try:
                chip._gain,chip._rdnoise = WFPC2_GAINS[chip.detnum][chip._headergain]
            except KeyError:
                raise ValueError("! Header gain value is not valid for WFPC2")


# ----------------------------------------------------------------------------
# Functions to convert WFPC2 C0M/C1M files into a single MEF FLT file
#
# ----------------------------------------------------------------------------
def wfpc2_to_flt(imgname):
    """Convert separate GEIS-based FITS files into single FLT file

    Parameters
    ----------
    imgname : str
        Filename of calibrated WFPC2 SCI (*_c0m.fits) image

    Returns
    -------
    flt_filename : str
        Filename of WFPC2 MEF _*flt.fits file that was written out

    """
    is_mef = 'c0m' in imgname
    if not is_mef:
        raise TypeError("MEF C0M file needed as input.")

    dq_file = imgname.replace('c0m', 'c1m')
    is_dq = os.path.exists(dq_file)
    flt_filename = imgname.replace('c0m', 'flt')

    # Read in input SCI file
    in_sci = fits.open(imgname)

    # Add keywords to be more compatible with ACS and WFC3 data
    num_sci = fileutil.countExtn(imgname)
    det_name = 'PC'
    in_sci[0].header['DETECTOR'] = det_name
    in_sci[0].header['PRIMESI'] = det_name

    if is_dq:
        # Read in existing input DQ file
        in_dq = fits.open(dq_file)
        dq_extns = [extn for extn in in_dq[1:]]
    else:
        # Could not find a DQ file, so create empty DQ arrays
        # based on SCI arrays
        dq_extns = [extn for extn in copy.deepcopy(in_sci[1:])]
        for extn in dq_extns:
            extn.data = np.zeros(extn.data.shape, dtype=np.int32)

    # Update EXTNAME to be consistent with ACS and WFC3 DQ extname
    for i,extn in enumerate(dq_extns):
        extn.header['extname'] = 'DQ'
        extn.header['extver'] = i+1

    # Now create ERR arrays as well...
    err_extns =[extn for extn in copy.deepcopy(in_sci[1:])]
    for i, extn in enumerate(err_extns):
        # Initialize using Poisson error estimate
        extn.data = np.sqrt(extn.data)
        extn.header['extname'] = 'ERR'
        extn.header['extver'] = i+1

    # Create output FLT file now to avoid having astropy
    # create a tmp* file that doesn't always get cleaned up...
    out_hdu = copy.deepcopy(in_sci)
    fname_kw = out_hdu[0].header['filename']
    out_hdu[0].header['filename'] = f"{fname_kw[:-8]}flt.fits"
    for dq_extn, err_extn in zip(dq_extns, err_extns):
        out_hdu.append(dq_extn)
        out_hdu.append(err_extn)

    print(f"Writing out {flt_filename}")
    out_file = open(flt_filename, 'wb')
    out_hdu.writeto(out_file, overwrite=True)
    in_sci.close()
    del in_sci
    if is_dq:
        in_dq.close()

    return flt_filename


# ------------------------------------------------------
# Function for updating headers to latest
# reference files from CRDS
# ------------------------------------------------------
def apply_bestrefs(filename=None, dirname=None,
                   uref_path=None, crds_path=None,
                   reftypes=['idctab', 'dgeofile', 'offtab']):
    """Update WFPC2 data to use the latest reference files from CRDS

    .. note::
        See `https://hst-crds.stsci.edu/docs/cmdline_bestrefs/
        <https://hst-crds.stsci.edu/docs/cmdline_bestrefs/>`_
        for details on how to configure CRDS for your local system and for definitions
        of all the environment variables used by CRDS.

    Parameters
    ----------
    filename : str, optional
        Filename of input file to be updated.
        If not specified, **all RAW and calibrated files** from the
        current directory, or ``dirname`` directory if given, will
        be updated.

    dirname : str, optional
        Name of directory containing WFPC2 data to be updated.
        If not specified, current directory will be checked.

    uref_path : str, optional
        Path for ``uref`` directory on local system.
        If not provided, the one defined in `os.environ` will be used.

    crds_path : str
        Path for the ``CRDS_PATH`` directory on local system.
        .. warning:: If no value is specified and ``CRDS_PATH`` is not
        already defined locally, this function will create a temporary
        ``CRDS_PATH`` directory tree under the directory with the files
        to be processed, populate it with the latest reference files
        and mappings needed by CRDS, then delete it when done.

    reftypes : list, optional
        List of reference files to be updated.  If None or an empty list,
        all reference files will be updated.

    """
    if reftypes is None:
        reftypes = []

    starting_dir = os.getcwd()

    wfpc2_dir = dirname if dirname else starting_dir
    os.chdir(wfpc2_dir)

    c0m = []
    if filename:
        if os.path.exists(filename):
            # User specified a single filename to process (default case)
            d0m = [filename]
        elif '*' in filename:
            # User specified a wild-carded filename to use as input
            d0m = sorted(glob.glob(filename))
        else:
            # Single input filename provided that could not be found
            os.chdir(starting_dir)
            raise ValueError(f"WFPC2 image {filename} not found in {wfpc2_dir}")
    else:
        # Get the list of ALL WFPC2 images from the specified directory
        c0m = sorted(glob.glob("*c0m.fits"))
        d0m = sorted(glob.glob("*d0m.fits"))

    # Make sure we are in a directory with WFPC2 images
    if len(d0m) == 0:
        # Return to original starting directory
        os.chdir(starting_dir)
        raise ValueError(f"ERROR: No WFPC2 data in {wfpc2_dir}")

    orig_crds = {'CRDS_PATH': os.environ.get('CRDS_PATH'),
                 'uref': os.environ.get('uref'),
                 'CRDS_SERVER_URL': os.environ.get('CRDS_SERVER_URL'),
                 'CRDS_OBSERVATORY': os.environ.get('CRDS_OBSERVATORY')}

    # Now, define what CRDS directories will be used for this update...
    remove_local_cache = False
    sync_refs = False if os.environ.get('CRDS_READONLY_CACHE') == '1' else True
    if not orig_crds['CRDS_PATH']:
        # User has not set up any local CRDS cache, so
        # we need to define one under the current working directory
        crds_cache = os.path.join(wfpc2_dir, 'crds_cache')
        crds_hst_cache = os.path.join(crds_cache, 'references', 'hst')
        crds_map_cache = os.path.join(crds_cache, 'mappings', 'hst')
        crds_uref_path = uref_path if uref_path else os.path.join(crds_hst_cache, 'wfpc2', os.sep)
        if not os.path.exists(crds_hst_cache):
            os.makedirs(crds_hst_cache)
            print(f"Creating temporary CRDS cache directory: {crds_hst_cache}")
        if not os.path.exists(crds_map_cache):
            os.makedirs(crds_map_cache)
            print(f"Creating temporary CRDS MAP cache directory: {crds_map_cache}")
        remove_local_cache = True
    else:
        # User has already defined a local CRDS cache, *and* it already exists, so use it.
        if crds_path:
            crds_cache = crds_path
        else:
            crds_cache = orig_crds['CRDS_PATH'] if os.path.exists(orig_crds['CRDS_PATH']) else None

        crds_hst_cache = os.path.join(crds_cache, 'references', 'hst')
        if os.path.exists(orig_crds['uref']):
            # Use the `uref` defined by the user in their `os.environ`
            crds_uref_path = orig_crds['uref']
        else:
            # User wants to explicitly use a CRDS cache they specify on input
            # and default to a path based on CRDS_PATH if `uref_path` is not
            # actually passed in through this function
            crds_uref_path = uref_path if uref_path else os.path.join(crds_hst_cache, 'wfpc2', os.sep)

        # `orig_crds["CRDS_PATH"]` has been defined, now verify it is valid...
        if not os.path.exists(crds_cache):
            # If CRDS_PATH was defined, but does not exists, raise an Exception
            # so that the user can finish setting up their CRDS cache properly.
            os.chdir(starting_dir)
            raise EnvironmentError(f"CRDS_PATH was defined as {crds_cache}, but does not exist!")
        if not os.path.exists(crds_uref_path):
            # If CRDS_PATH was defined, but does not exists, raise an Exception
            # so that the user can finish setting up their CRDS cache properly.
            os.chdir(starting_dir)
            raise EnvironmentError(f"CRDS path to `uref` was defined as {crds_uref_path}, but does not exist!")

    # Now that we have confirmed we have images to update...
    # configure CRDS for use in updating the WFPC2 data
    os.environ['CRDS_SERVER_URL'] = orig_crds['CRDS_SERVER_URL'] if orig_crds['CRDS_SERVER_URL'] else PUBLIC_CRDS_SERVER_URL
    os.environ['CRDS_OBSERVATORY'] = "hst"
    os.environ['CRDS_PATH'] = crds_cache
    os.environ['uref'] = crds_uref_path

    #    os.environ['CRDS_PATH'] = "D:\data\crds_cache"
    #    os.environ['uref'] = "D:\\data\\crds_cache\\references\\hst\\wfpc2\\"

    # Only import this package if there is data to be updated
    import crds

    print(f"Running CRDS.assign_bestrefs on: {d0m} for reftypes={reftypes}")
    # Apply bestrefs to images after downloading references to local CRDS cache
    crds.assign_bestrefs(d0m, reftypes=reftypes, sync_references=sync_refs)
    if len(c0m) > 0:
        crds.assign_bestrefs(c0m, reftypes=reftypes, sync_references=sync_refs)

    # clean up temp crds_cache dir, if created
    if remove_local_cache:
        shutil.rmtree(crds_cache)

    # Now, revert os.environ to original state prior to running this function
    for crds_key, crds_val in orig_crds.items():
        if crds_val is None:
            # If it was originally None, then
            # there was no key defined originally, so delete it.
            del os.environ[crds_key]
        else:
            os.environ[crds_key] = crds_val

    # Return to original starting directory
    os.chdir(starting_dir)
