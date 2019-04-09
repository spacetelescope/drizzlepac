"""
A class which makes image objects for each input filename.

:Authors: Warren Hack

:License: :doc:`LICENSE`

"""
import copy, os, re, sys

import numpy as np
from stwcs import distortion

from stsci.tools import fileutil, logutil, textutil
from astropy.io import fits
from . import util
from . import wcs_functions
from . import buildmask
from .version import *

__all__ = ['baseImageObject', 'imageObject', 'WCSObject']


log = logutil.create_logger(__name__, level=logutil.logging.NOTSET)


_NUMPY_TO_IRAF_DTYPES = {'float64': -64, 'float32': -32, 'uint8': 8,
                         'int16': 16, 'int32': 32, 'int64': 64}

_IRAF_DTYPES_TO_NUMPY = {-64: 'float64', -32: 'float32', 8: 'uint8',
                         16: 'int16', 32: 'int32', 64: 'int64'}


class baseImageObject:
    """ Base ImageObject which defines the primary set of methods. """
    def __init__(self,filename):
        """
        """
        self.scienceExt= "SCI" # the extension the science image is stored in
        self.maskExt="DQ" #the extension with the mask image in it
        self.errExt = "ERR"  # the extension the ERR array can be found in
        self._filename = filename
        self._original_file_name = filename
        self.native_units='ELECTRONS'

        self.flatkey = None  # keyword which points to flat-field reference file

        self._image = None
        self._instrument=None
        self._rootname=None
        self.outputNames={}
        self.outputValues = {}
        self.createContext = True

        self.inmemory = False # flag for all in-memory operations
        #this is the number of science chips to be processed in the file
        self._numchips=1
        self._nextend=0
        # this is the number of chip which will be combined based on 'group' parameter
        self._nmembers = 0

    def __getitem__(self,exten):
        """ Overload  getitem to return the data and header
            these only work on the HDU list already in memory
            once the data has been zero's in self._image you should
            use getData or getHeader to re-read the file.
        """
        return fileutil.getExtn(self._image,extn=exten)

    def __cmp__(self, other):
        """ Overload the comparison operator
            just to check the filename of the object?
        """
        return (isinstance(other, imageObject) and self._filename == other._filename)

    def _isNotValid(self, par1, par2):
        """ Method used to determine if a value or keyword is
            supplied as input for instrument specific parameters.
        """
        invalidValues = [None,'None','INDEF','']
        return (par1 in invalidValues and par2 in invalidValues)

    def info(self):
        """ Return fits information on the _image.
        """
        #if the file hasn't been closed yet then we can
        #use the fits info which looks at the extensions
        #if(self._isSimpleFits):
        #    print self._filename," is a simple fits image"
        #else:
        self._image.info()

    def close(self):
        """ Close the object nicely and release all the data
            arrays from memory YOU CANT GET IT BACK, the pointers
            and data are gone so use the getData method to get
            the data array returned for future use. You can use
            putData to reattach a new data array to the imageObject.
        """
        if self._image is None:
            return

        # mcara: I think the code below is not necessary but in order to
        #        preserve the same functionality as the code removed below,
        #        I make an empty copy of the image object:
        empty_image = fits.HDUList()
        for u in self._image:
            empty_image.append(u.__class__(data=None, header=None))
        # mcara: END unnecessary code

        self._image.close()  #calls fits.close()

        self._image = empty_image

        #we actuallly want to make sure that all the
        #data extensions have been closed and deleted
        #since we could have the DQ,ERR and others read in
        #at this point, but I'd like there to be something
        #valid there afterwards that I can play with

        # mcara: REMOVED unnecessary code:
        #
        # if not self._isSimpleFits:
        #     for ext,hdu in enumerate(self._image):
        #         #use the datatype for the extension
        #         #dtype=self.getNumpyType(hdu.header["BITPIX"])
        #         hdu.data = None #np.array(0,dtype=dtype)  #so we dont get io errors on stuff that wasn't read in yet
        # else:
        #     self._image.data= None # np.array(0,dtype=self.getNumpyType(self._image.header["BITPIX"]))

    def clean(self):
        """ Deletes intermediate products generated for this imageObject.
        """
        clean_files = ['blotImage','crmaskImage','finalMask',
                        'staticMask','singleDrizMask','outSky',
                        'outSContext','outSWeight','outSingle',
                        'outMedian','dqmask','tmpmask',
                        'skyMatchMask']

        log.info('Removing intermediate files for %s' % self._filename)
        # We need to remove the combined products first; namely, median image
        util.removeFileSafely(self.outputNames['outMedian'])
        # Now remove chip-specific intermediate files, if any were created.
        for chip in self.returnAllChips(extname='SCI'):
            for fname in clean_files:
                if fname in chip.outputNames:
                    util.removeFileSafely(chip.outputNames[fname])

    def getData(self,exten=None):
        """ Return just the data array from the specified extension
            fileutil is used instead of fits to account for non-
            FITS input images. openImage returns a fits object.
        """
        if exten.lower().find('sci') > -1:
            # For SCI extensions, the current file will have the data
            fname = self._filename
        else:
            # otherwise, the data being requested may need to come from a
            # separate file, as is the case with WFPC2 DQ data.
            #
            # convert exten to 'sci',extver to get the DQ info for that chip
            extn = exten.split(',')
            sci_chip = self._image[self.scienceExt,int(extn[1])]
            fname = sci_chip.dqfile

        extnum = self._interpretExten(exten)
        if self._image[extnum].data is None:
            if os.path.exists(fname):
                _image=fileutil.openImage(fname, clobber=False, memmap=False)
                _data=fileutil.getExtn(_image, extn=exten).data
                _image.close()
                del _image
                self._image[extnum].data = _data
            else:
                _data = None
        else:
            _data = self._image[extnum].data

        return _data

    def getHeader(self,exten=None):
        """ Return just the specified header extension fileutil
            is used instead of fits to account for non-FITS
            input images. openImage returns a fits object.
        """
        _image=fileutil.openImage(self._filename, clobber=False, memmap=False)
        _header=fileutil.getExtn(_image,extn=exten).header
        _image.close()
        del _image
        return _header

    def _interpretExten(self,exten):
        #check if the exten is a string or number and translate to the correct chip
        _extnum=0

        if ',' in str(exten): #assume a string like "sci,1" has been given
            _extensplit=exten.split(',')
            _extname=_extensplit[0]
            _extver=int(_extensplit[1])
            _extnum=self.findExtNum(_extname,_extver)
        else:
            #assume that a direct extnum has been given
            _extnum=int(exten)

        if _extnum is None:
            msg = "no extension number found"
            log.error(msg)
            raise ValueError(msg)

        return _extnum

    def updateData(self,exten,data):
        """ Write out updated data and header to
            the original input file for this object.
        """
        _extnum=self._interpretExten(exten)
        fimg = fileutil.openImage(self._filename, mode='update', memmap=False)
        fimg[_extnum].data = data
        fimg[_extnum].header = self._image[_extnum].header
        fimg.close()

    def putData(self,data=None,exten=None):
        """ Now that we are removing the data from the object to save memory,
            we need something that cleanly puts the data array back into
            the object so that we can write out everything together  using
            something like fits.writeto....this method is an attempt to
            make sure that when you add an array back to the .data section
            of the hdu it still matches the header information for that
            section ( ie. update the bitpix to reflect the datatype of the
            array you are adding). The other header stuff is up to you to verify.

            Data should be the data array exten is where you want to stick it,
            either extension number or a string like 'sci,1'
        """

        if data is None:
            log.warning("No data supplied")
        else:
            extnum = _interpretExten(exten)
            ext = self._image[extnum]
            # update the bitpix to the current datatype, this aint fancy and
            # ignores bscale
            ext.header['BITPIX'] = _NUMPY_TO_IRAF_DTYPES[data.dtype.name]
            ext.data = data

    def getAllData(self,extname=None,exclude=None):
        """ This function is meant to make it easier to attach ALL the data
            extensions of the image object so that we can write out copies of
            the original image nicer.

            If no extname is given, the it retrieves all data from the original
            file and attaches it. Otherwise, give the name of the extensions
            you want and all of those will be restored.

            Ok, I added another option. If you want to get all the data
            extensions EXCEPT a particular one, leave extname=NONE and
            set exclude=EXTNAME. This is helpfull cause you might not know
            all the extnames the image has, this will find out and exclude
            the one you do not want overwritten.
        """

        extensions = self._findExtnames(extname=extname,exclude=exclude)

        for i in range(1,self._nextend+1,1):
            if hasattr(self._image[i],'_extension') and \
                "IMAGE" in self._image[i]._extension:
                extver = self._image[i].header['extver']
                if (self._image[i].extname in extensions) and self._image[self.scienceExt,extver].group_member:
                    self._image[i].data=self.getData(self._image[i].extname + ','+str(self._image[i].extver))

    def returnAllChips(self,extname=None,exclude=None):
        """ Returns a list containing all the chips which match the
            extname given minus those specified for exclusion (if any).
        """
        extensions = self._findExtnames(extname=extname,exclude=exclude)
        chiplist = []
        for i in range(1,self._nextend+1,1):
            if 'extver' in self._image[i].header:
                extver = self._image[i].header['extver']
            else:
                extver = 1
            if hasattr(self._image[i],'_extension') and \
                "IMAGE" in self._image[i]._extension:
                if (self._image[i].extname in extensions) and self._image[self.scienceExt,extver].group_member:
                    chiplist.append(self._image[i])
        return chiplist

    def _findExtnames(self, extname=None, exclude=None):
        """ This method builds a list of all extensions which have 'EXTNAME'==extname
            and do not include any extensions with 'EXTNAME'==exclude, if any are
            specified for exclusion at all.
        """
        #make a list of the available extension names for the object
        extensions=[]
        if extname is not None:
            if not isinstance(extname,list): extname=[extname]
            for extn in extname:
                extensions.append(extn.upper())
        else:
        #restore all the extensions data from the original file, be careful here
        #if you've altered data in memory you want to keep!
            for i in range(1,self._nextend+1,1):
                if hasattr(self._image[i],'_extension') and \
                    "IMAGE" in self._image[i]._extension:
                    if self._image[i].extname.upper() not in extensions:
                        extensions.append(self._image[i].extname)
        #remove this extension from the list
        if exclude is not None:
            exclude.upper()
            if exclude in extensions:
                newExt=[]
                for item in extensions:
                    if item != exclude:
                        newExt.append(item)
            extensions=newExt
            del newExt
        return extensions

    def findExtNum(self, extname=None, extver=1):
        """Find the extension number of the give extname and extver."""
        extnum = None
        extname = extname.upper()

        if not self._isSimpleFits:
            for ext in self._image:
                if (hasattr(ext,'_extension') and 'IMAGE' in ext._extension and
                    (ext.extname == extname) and (ext.extver == extver)):
                    extnum = ext.extnum
        else:
            log.info("Image is simple fits")

        return extnum

    def _assignRootname(self, chip):
        """ Assign a unique rootname for the image based in the expname. """
        extname=self._image[self.scienceExt,chip].header["EXTNAME"].lower()
        extver=self._image[self.scienceExt,chip].header["EXTVER"]
        expname = self._rootname

        # record extension-based name to reflect what extension a mask file corresponds to
        self._image[self.scienceExt,chip].rootname=expname + "_" + extname + str(extver)
        self._image[self.scienceExt,chip].sciname=self._filename + "[" + extname +","+str(extver)+"]"
        self._image[self.scienceExt,chip].dqrootname=self._rootname + "_" + extname + str(extver)
        # Needed to keep EXPNAMEs associated properly (1 EXPNAME for all chips)
        self._image[self.scienceExt,chip]._expname=expname
        self._image[self.scienceExt,chip]._chip =chip

    def _setOutputNames(self,rootname,suffix='_drz'):
        """ Define the default output filenames for drizzle products,
            these are based on the original rootname of the image
            filename should be just 1 filename, so call this in a loop
            for chip names contained inside a file.
        """
        # Define FITS output filenames for intermediate products

        # Build names based on final DRIZZLE output name
        # where 'output' normally would have been created
        #   by 'process_input()'
        #
        outFinal = rootname+suffix+'.fits'
        outSci = rootname+suffix+'_sci.fits'
        outWeight = rootname+suffix+'_wht.fits'
        outContext = rootname+suffix+'_ctx.fits'
        outMedian = rootname+'_med.fits'

        # Build names based on input name
        origFilename = self._filename.replace('.fits','_OrIg.fits')
        outSky = rootname + '_sky.fits'
        outSingle = rootname+'_single_sci.fits'
        outSWeight = rootname+'_single_wht.fits'
        crCorImage = rootname+'_crclean.fits'

        # Build outputNames dictionary
        fnames={
            'origFilename': origFilename,
            'outFinal': outFinal,
            'outMedian': outMedian,
            'outSci': outSci,
            'outWeight': outWeight,
            'outContext': outContext,
            'outSingle': outSingle,
            'outSWeight': outSWeight,
            'outSContext': None,
            'outSky': outSky,
            'crcorImage': crCorImage,
            'ivmFile': None
        }

        return fnames

    def _setChipOutputNames(self,rootname,chip):
        blotImage = rootname + '_blt.fits'
        crmaskImage = rootname + '_crmask.fits'

        # Start with global names
        fnames = self.outputNames

        # Now add chip-specific entries
        fnames['blotImage'] = blotImage
        fnames['crmaskImage'] = crmaskImage
        sci_chip = self._image[self.scienceExt,chip]
        # Define mask names as additional entries into outputNames dictionary
        fnames['finalMask']=sci_chip.dqrootname+'_final_mask.fits' # used by final_drizzle
        fnames['singleDrizMask']=fnames['finalMask'].replace('final','single')
        fnames['staticMask']=None

        # Add the following entries for use in creating outputImage object
        fnames['data'] = sci_chip.sciname
        return fnames

    ##############################################################
    #
    # Methods related to managing virtual intermediate output products
    # as opposed to writing them out as files on disk
    #
    ##############################################################
    def _initVirtualOutputs(self):
        """ Sets up the structure to hold all the output data arrays for
        this image in memory.
        """
        self.virtualOutputs = {}
        for product in self.outputNames:
            self.virtualOutputs[product] = None

    def saveVirtualOutputs(self,outdict):
        """ Assign in-memory versions of generated products for this
        ``imageObject`` based on dictionary 'outdict'.
        """
        if not self.inmemory:
            return
        for outname in outdict:
            self.virtualOutputs[outname] = outdict[outname]

    def getOutputName(self,name):
        """ Return the name of the file or PyFITS object associated with that
        name, depending on the setting of self.inmemory.
        """
        val = self.outputNames[name]
        if self.inmemory: # if inmemory was turned on...
            # return virtualOutput object saved with that name
            val = self.virtualOutputs[val]
        return val

    ##############################################################
    # Methods for managing output values associated with this input
    ##############################################################
    def updateOutputValues(self,output_wcs):
        """ Copy info from output WCSObject into outputnames for each chip
        for use in creating outputimage object.
        """
        outputvals = self.outputValues

        outputvals['output'] = output_wcs.outputNames['outFinal']
        outputvals['outnx'], outputvals['outny'] = output_wcs.wcs.pixel_shape
        outputvals['texptime'] = output_wcs._exptime
        outputvals['texpstart'] = output_wcs._expstart
        outputvals['texpend'] = output_wcs._expend
        outputvals['nimages'] = output_wcs.nimages

        outputvals['scale'] = output_wcs.wcs.pscale #/ self._image[self.scienceExt,1].wcs.pscale
        outputvals['exptime'] = self._exptime

        outnames = self.outputNames
        outnames['outMedian'] = output_wcs.outputNames['outMedian']
        outnames['outFinal'] = output_wcs.outputNames['outFinal']
        outnames['outSci'] = output_wcs.outputNames['outSci']
        outnames['outWeight'] = output_wcs.outputNames['outWeight']
        outnames['outContext'] = output_wcs.outputNames['outContext']

    def updateContextImage(self, contextpar):
        """ Reset the name of the context image to `None` if parameter
        ``context`` is `False`.

        """
        self.createContext = contextpar
        if not contextpar:
            log.info('No context image will be created for %s' %
                     self._filename)
            self.outputNames['outContext'] = None

    def find_DQ_extension(self):
        """ Return the suffix for the data quality extension and the name of
        the file which that DQ extension should be read from.
        """
        dqfile = None
        dq_suffix=None
        if(self.maskExt is not None):
            for hdu in self._image:
                # Look for DQ extension in input file
                if 'extname' in hdu.header and hdu.header['extname'].lower() == self.maskExt.lower():
                    dqfile = self._filename
                    dq_suffix=self.maskExt
                    break

        return dqfile,dq_suffix

    def getKeywordList(self, kw):
        """
        Return lists of all attribute values for all active chips in the
        ``imageObject``.

        """
        kwlist = []
        for chip in range(1,self._numchips+1,1):
            sci_chip = self._image[self.scienceExt,chip]
            if sci_chip.group_member:
                kwlist.append(sci_chip.__dict__[kw])

        return kwlist

    def getGain(self, exten):
        return self._image[exten]._gain

    def getflat(self, chip):
        """
        Method for retrieving a detector's flat field.

        Returns
        -------
        flat: array
            This method will return an array the same shape as the image in
            **units of electrons**.

        """
        sci_chip = self._image[self.scienceExt, chip]
        # The keyword for ACS flat fields in the primary header of the flt
        # file is pfltfile.  This flat file is already in the required
        # units of electrons.

        # The use of fileutil.osfn interprets any environment variable, such as
        # jref$, used in the specification of the reference filename
        filename = fileutil.osfn(self._image["PRIMARY"].header[self.flatkey])
        hdulist = None
        try:
            hdulist = fileutil.openImage(filename, mode='readonly',
                                         memmap=False)
            data = hdulist[(self.scienceExt, chip)].data

            if data.shape[0] != sci_chip.image_shape[0]:
                ltv2 = int(np.round(sci_chip.ltv2))
            else:
                ltv2 = 0
            size2 = sci_chip.image_shape[0] + ltv2

            if data.shape[1] != sci_chip.image_shape[1]:
                ltv1 = int(np.round(sci_chip.ltv1))
            else:
                ltv1 = 0
            size1 = sci_chip.image_shape[1] + ltv1

            flat = data[ltv2:size2, ltv1:size1]

        except FileNotFoundError:
            flat = np.ones(sci_chip.image_shape, dtype=sci_chip.image_dtype)
            log.warning("Cannot find flat field file '{}'".format(filename))
            log.warning("Treating flatfield as a constant value of '1'.")

        finally:
            if hdulist is not None:
                hdulist.close()

        return flat

    def getReadNoiseImage(self, chip):
        """
        Notes
        =====
        Method for returning the readnoise image of a detector
        (in electrons).

        The method will return an array of the same shape as the image.

        :units: electrons

        """
        sci_chip = self._image[self.scienceExt,chip]
        return np.ones(sci_chip.image_shape,dtype=sci_chip.image_dtype) * sci_chip._rdnoise

    def getexptimeimg(self,chip):
        """
        Notes
        =====
        Return an array representing the exposure time per pixel for the detector.
        This method will be overloaded for IR detectors which have their own
        EXP arrays, namely, WFC3/IR and NICMOS images.

        :units:
          None

        Returns
        =======
        exptimeimg : numpy array
            The method will return an array of the same shape as the image.

        """
        sci_chip = self._image[self.scienceExt,chip]
        if sci_chip._wtscl_par == 'expsq':
            wtscl = sci_chip._exptime*sci_chip._exptime
        else:
            wtscl = sci_chip._exptime

        return np.ones(sci_chip.image_shape,dtype=sci_chip.image_dtype)*wtscl

    def getdarkimg(self,chip):
        """
        Notes
        =====
        Return an array representing the dark image for the detector.

        The method will return an array of the same shape as the image.

        :units: electrons
        """
        sci_chip = self._image[self.scienceExt,chip]
        return np.ones(sci_chip.image_shape,dtype=sci_chip.image_dtype)*sci_chip.darkcurrent

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
        return np.ones(sci_chip.image_shape,dtype=sci_chip.image_dtype)*sci_chip.subtractedSky

    def getdarkcurrent(self):
        """
        Notes
        =====
        Return the dark current for the detector.  This value
        will be contained within an instrument specific keyword.
        The value in the image header will be converted to units
        of electrons.

        :units: electrons

        """
        pass

    # the following two functions are basically doing the same thing,
    # how are they used differently in the code?
    def getExtensions(self, extname='SCI', section=None):
        """ Return the list of EXTVER values for extensions with name specified
        in extname.

        """
        if section is None:
            numext = 0
            section = []
            for hdu in self._image:
                if 'extname' in hdu.header and hdu.header['extname'] == extname:
                    section.append(hdu.header['extver'])
        else:
            if not isinstance(section,list):
                section = [section]

        return section

    def _countEXT(self,extname="SCI"):
        """ Count the number of extensions in the file with the given name
        (``EXTNAME``).

        """
        count=0 #simple fits image

        if (self._image['PRIMARY'].header["EXTEND"]):
            for i,hdu in enumerate(self._image):
                if i > 0:
                    hduExtname = False
                    if 'EXTNAME' in hdu.header:
                        self._image[i].extnum=i
                        self._image[i].extname=hdu.header["EXTNAME"]
                        hduExtname = True
                    if 'EXTVER' in hdu.header:
                        self._image[i].extver=hdu.header["EXTVER"]
                    else:
                        self._image[i].extver = 1

                    if ((extname is not None) and \
                            (hduExtname and (hdu.header["EXTNAME"] == extname))) \
                            or extname is None:
                        count=count+1
        return count

    def getNumpyType(self, irafType):
        """ Return the corresponding numpy data type. """
        return _IRAF_DTYPES_TO_NUMPY[irafType]

    def buildMask(self,chip,bits=0,write=False):
        """
        Build masks as specified in the user parameters found in the
        configObj object.

        We should overload this function in the instrument specific
        implementations so that we can add other stuff to the badpixel
        mask? Like vignetting areas and chip boundries in nicmos which
        are camera dependent? these are not defined in the DQ masks, but
        should be masked out to get the best results in multidrizzle.
        """
        dqarr = self.getData(exten=self.maskExt+','+str(chip))
        dqmask = buildmask.buildMask(dqarr,bits)

        if write:
            phdu = fits.PrimaryHDU(data=dqmask,header=self._image[self.maskExt,chip].header)
            dqmask_name = self._image[self.scienceExt,chip].dqrootname+'_dqmask.fits'
            log.info('Writing out DQ/weight mask: %s' % dqmask_name)
            if os.path.exists(dqmask_name): os.remove(dqmask_name)
            phdu.writeto(dqmask_name)
            del phdu
            self._image[self.scienceExt,chip].dqmaskname = dqmask_name
            # record the name of this mask file that was created for later
            # removal by the 'clean()' method
            self._image[self.scienceExt,chip].outputNames['dqmask'] = dqmask_name
        del dqarr
        return dqmask

    def buildEXPmask(self, chip, dqarr):
        """ Builds a weight mask from an input DQ array and the exposure time
        per pixel for this chip.
        """
        log.info("Applying EXPTIME weighting to DQ mask for chip %s" %
                 chip)
        #exparr = self.getexptimeimg(chip)
        exparr = self._image[self.scienceExt,chip]._exptime
        expmask = exparr*dqarr

        return expmask.astype(np.float32)

    def buildIVMmask(self ,chip, dqarr, scale):
        """ Builds a weight mask from an input DQ array and either an IVM array
        provided by the user or a self-generated IVM array derived from the
        flat-field reference file associated with the input image.
        """
        sci_chip = self._image[self.scienceExt,chip]
        ivmname = self.outputNames['ivmFile']

        if ivmname is not None:
            log.info("Applying user supplied IVM files for chip %s" % chip)
            #Parse the input file name to get the extension we are working on
            extn = "IVM,{}".format(chip)

            #Open the mask image for updating and the IVM image
            ivm =  fileutil.openImage(ivmname, mode='readonly', memmap=False)
            ivmfile = fileutil.getExtn(ivm, extn)

            # Multiply the IVM file by the input mask in place.
            ivmarr = ivmfile.data * dqarr

            ivm.close()

        else:
            log.info("Automatically creating IVM files for chip %s" % chip)
            # If no IVM files were provided by the user we will
            # need to automatically generate them based upon
            # instrument specific information.

            flat = self.getflat(chip)
            RN = self.getReadNoiseImage(chip)
            darkimg = self.getdarkimg(chip)
            skyimg = self.getskyimg(chip)

            #exptime = self.getexptimeimg(chip)
            #exptime = sci_chip._exptime
            #ivm = (flat*exptime)**2/(darkimg+(skyimg*flat)+RN**2)
            ivm = (flat)**2/(darkimg+(skyimg*flat)+RN**2)

           # Multiply the IVM file by the input mask in place.
            ivmarr = ivm * dqarr

        # Update 'wt_scl' parameter to match use of IVM file
        sci_chip._wtscl = pow(sci_chip._exptime,2)/pow(scale,4)
        #sci_chip._wtscl = 1.0/pow(scale,4)

        return ivmarr.astype(np.float32)

    def buildERRmask(self,chip,dqarr,scale):
        """
        Builds a weight mask from an input DQ array and an ERR array
        associated with the input image.
        """
        sci_chip = self._image[self.scienceExt,chip]

        # Set default value in case of error, or lack of ERR array
        errmask = dqarr

        if self.errExt is not None:
            try:
                # Attempt to open the ERR image.
                err = self.getData(exten=self.errExt+','+str(chip))

                log.info("Applying ERR weighting to DQ mask for chip %s" %
                         chip)

                # Multiply the scaled ERR file by the input mask in place.
                #exptime = self.getexptimeimg(chip)
                exptime = sci_chip._exptime
                errmask = (exptime/err)**2 * dqarr

                # Update 'wt_scl' parameter to match use of IVM file
                #sci_chip._wtscl = pow(sci_chip._exptime,2)/pow(scale,4)
                sci_chip._wtscl = 1.0/pow(scale,4)

                del err

            except:
                # We cannot find an 'ERR' extension and the data isn't WFPC2.
                # Print a generic warning message and continue on with the
                # final drizzle step.

                print(textutil.textbox(
                    'WARNING: No ERR weighting will be applied to the mask '
                    'used in the final drizzle step!  Weighting will be only '
                    'by exposure time.\n\nThe data provided as input does not '
                    'contain an ERR extension'), file=sys.stderr)
                print('\n Continue with final drizzle step...', sys.stderr)
        else:
            # If we were unable to find an 'ERR' extension to apply, one
            # possible reason was that the input was a 'standard' WFPC2 data
            # file that does not actually contain an error array.  Test for
            # this condition and issue a Warning to the user and continue on to
            # the final drizzle.

            print(textutil.textbox(
                "WARNING: No ERR weighting will be applied to the mask used "
                "in the final drizzle step!  Weighting will be only by "
                "exposure time.\n\nThe WFPC2 data provided as input does not "
                "contain ERR arrays.  WFPC2 data is not supported by this "
                "weighting type.\n\nA workaround would be to create inverse "
                "variance maps and use 'IVM' as the final_wht_type.  See the "
                "HELP file for more details on using inverse variance maps."),
                file=sys.stderr)
            print("\n Continue with final drizzle step...", file=sys.stderr)

        return errmask.astype(np.float32)

    def updateIVMName(self,ivmname):
        """ Update outputNames for image with user-supplied IVM filename.
        """
        self.outputNames['ivmFile'] = ivmname

    def set_mt_wcs(self, image):
        """ Reset the WCS for this image based on the WCS information from
            another imageObject.
        """
        for chip in range(1,self._numchips+1,1):
            sci_chip = self._image[self.scienceExt,chip]
            ref_chip = image._image[image.scienceExt,chip]
            # Do we want to keep track of original WCS or not? No reason now...
            sci_chip.wcs = ref_chip.wcs.copy()

    def set_wtscl(self, chip, wtscl_par):
        """ Sets the value of the wt_scl parameter as needed for drizzling.
        """
        sci_chip = self._image[self.scienceExt,chip]

        exptime = 1 #sci_chip._exptime
        _parval = 'unity'
        if wtscl_par is not None:
            if type(wtscl_par) == type(''):
                if not wtscl_par.isdigit():
                    # String passed in as value, check for 'exptime' or 'expsq'
                    _wtscl_float = None
                    try:
                        _wtscl_float = float(wtscl_par)
                    except ValueError:
                        _wtscl_float = None
                    if _wtscl_float is not None:
                        _wtscl = _wtscl_float
                    elif wtscl_par == 'expsq':
                        _wtscl = exptime*exptime
                        _parval = 'expsq'
                    else:
                        # Default to the case of 'exptime', if
                        #   not explicitly specified as 'expsq'
                        _wtscl = exptime
                else:
                    # int value passed in as a string, convert to float
                    _wtscl = float(wtscl_par)
            else:
                # We have a non-string value passed in...
                _wtscl = float(wtscl_par)
        else:
            # Default case: wt_scl = exptime
            _wtscl = exptime

        sci_chip._wtscl_par = _parval
        sci_chip._wtscl = _wtscl

    def set_units(self):
        """ Record the units for this image, both BUNITS from header and
            in_units as needed internally. This method will be defined
            specifically for each instrument.
        """
        pass

    def getInstrParameter(self, value, header, keyword):
        """ This method gets a instrument parameter from a
            pair of task parameters: a value, and a header keyword.

            The default behavior is:
              - if the value and header keyword are given, raise an exception.
              - if the value is given, use it.
              - if the value is blank and the header keyword is given, use
                the header keyword.
              - if both are blank, or if the header keyword is not
                found, return None.
        """
        if isinstance(value, str) and value in ['None', '', ' ', 'INDEF']:
            value = None

        if value and (keyword is not None and keyword.strip() != ''):
            exceptionMessage = "ERROR: Your input is ambiguous!  Please specify either a value or a keyword.\n  You specifed both " + str(value) + " and " + str(keyword)
            raise ValueError(exceptionMessage)

        elif value is not None and value != '':
            return self._averageFromList(value)

        elif keyword is not None and keyword.strip() != '':
            return self._averageFromHeader(header, keyword)

        else:
            return None

    def _averageFromHeader(self, header, keyword):
        """ Averages out values taken from header. The keywords where
            to read values from are passed as a comma-separated list.
        """
        _list = ''
        for _kw in keyword.split(','):
            if _kw in header:
                _list = _list + ',' + str(header[_kw])
            else:
                return None
        return self._averageFromList(_list)

    def _averageFromList(self, param):
        """ Averages out values passed as a comma-separated
            list, disregarding the zero-valued entries.
        """
        _result = 0.0
        _count = 0

        for _param in param.split(','):
            if _param != '' and float(_param) != 0.0:
                _result = _result + float(_param)
                _count  += 1

        if _count >= 1:
            _result = _result / _count
        return _result

class imageObject(baseImageObject):
    """
        This returns an imageObject that contains all the
        necessary information to run the image file through
        any multidrizzle function. It is essentially a
        PyFits object with extra attributes.

        There will be generic keywords which are good for
        the entire image file, and some that might pertain
        only to the specific chip.
    """

    def __init__(self,filename,group=None,inmemory=False):
        super().__init__(filename)

        #filutil open returns a fits object
        try:
            self._image=fileutil.openImage(filename, clobber=False, memmap=False)

        except IOError:
            raise IOError("Unable to open file: %s" % filename)

        #populate the global attributes which are good for all the chips in the file
        #self._rootname=self._image['PRIMARY'].header["ROOTNAME"]
        self._rootname=fileutil.buildNewRootname(filename)
        self.outputNames=self._setOutputNames(self._rootname)

        # flag to indicate whether or not to write out intermediate products
        # to disk (default) or keep everything in memory
        self.inmemory = inmemory
        self._initVirtualOutputs()

        #self._exptime=self._image["PRIMARY"].header["EXPTIME"]
        #exptime should be set in the image subclass code since it's kept in different places
#        if(self._exptime == 0):
        self._exptime =1. #to avoid divide by zero
 #           print "Setting exposure time to 1. to avoid div/0!"

        #this is the number of science chips to be processed in the file
        self._numchips=self._countEXT(extname=self.scienceExt)

        self.proc_unit = None

        #self._nextend=self._image["PRIMARY"].header["NEXTEND"]
        self._nextend = self._countEXT(extname=None)

        if (self._numchips == 0):
            #the simple fits image contains the data in the primary extension,
            #this will help us deal with the rest of the code that looks
            #and acts on chips :)
            #self._nextend=1
            self._numchips=1
            self.scienceExt="PRIMARY"
            self.maskExt=None
            self._image["PRIMARY"].header["EXTNAME"] = "PRIMARY"
            self._image["PRIMARY"].header["EXTVER"] = 1
            self._image["PRIMARY"].extnum = 0

        self._isSimpleFits = False

        # Clean out any stray MDRIZSKY keywords from PRIMARY headers
        fimg = fileutil.openImage(filename, mode='update', memmap=False)
        if 'MDRIZSKY' in fimg['PRIMARY'].header:
            del fimg['PRIMARY'].header['MDRIZSKY']
        fimg.close()
        del fimg

        if group not in [None,'']:
            # Only use selected chip
            if ',' in group:
                group_id = group.split(',')
                if group_id[0].isalpha(): # user specified a specific extname,extver
                    self.group = [int(group_id[1])]
                else: # user specified a list of extension numbers to process
                    self.group = []
                    for grp in group_id:
                        # find extname/extver which corresponds to this extension number
                        group_extname = self._image[int(grp)].header['EXTNAME']
                        group_extver = self._image[int(grp)].header['EXTVER']
                        self.group.append(group_extver)
            else:
                # find extname/extver which corresponds to this extension number
                group_extver = self._image[int(group)].header['EXTVER']
                self.group = [int(group_extver)]
        else:
            # Use all chips
            self.group = None

        if not self._isSimpleFits:

            #assign chip specific information
            for chip in range(1,self._numchips+1,1):

                self._assignRootname(chip)
                sci_chip = self._image[self.scienceExt,chip]

                # Set a flag to indicate whether this chip should be included
                # or not, based on user input from the 'group' parameter.
                if self.group is None or (self.group is not None and chip in self.group):
                    sci_chip.group_member = True
                    self._nmembers += 1
                else:
                    sci_chip.group_member = False

                sci_chip.signature = None

                sci_chip.dqname = None
                sci_chip.dqmaskname = None

                sci_chip.dqfile,sci_chip.dq_extn = self.find_DQ_extension()
                #self.maskExt = sci_chip.dq_extn
                if(sci_chip.dqfile is not None):
                    sci_chip.dqname = sci_chip.dqfile +'['+sci_chip.dq_extn+','+str(chip)+']'

                # build up HSTWCS object for each chip, which will be necessary for drizzling operations
                sci_chip.wcs=wcs_functions.get_hstwcs(self._filename,self._image,sci_chip.extnum)
                sci_chip.detnum,sci_chip.binned = util.get_detnum(sci_chip.wcs,self._filename,chip)
                sci_chip.wcslin_pscale = 1.0

                #assuming all the chips don't have the same dimensions in the file
                sci_chip._naxis1=sci_chip.header["NAXIS1"]
                sci_chip._naxis2=sci_chip.header["NAXIS2"]

                # record the exptime values for this chip so that it can be
                # easily used to generate the composite value for the final output image
                sci_chip._expstart,sci_chip._expend = util.get_expstart(sci_chip.header,self._image['PRIMARY'].header)

                sci_chip.outputNames=self._setChipOutputNames(sci_chip.rootname,chip).copy() #this is a dictionary
                # Set the units: both bunit and in_units
                self.set_units(chip)

                #initialize gain, readnoise, and exptime attributes
                # the actual values will be set by each instrument based on
                # keyword names specific to that instrument by 'setInstrumentParamters()'
                sci_chip._headergain = 1 # gain value read from header
                sci_chip._gain = 1.0     # calibrated gain value
                sci_chip._rdnoise = 1.0  # calibrated readnoise
                sci_chip._exptime = 1.0
                sci_chip._effGain = 1.0
                sci_chip._conversionFactor = 1.0
                sci_chip._wtscl = 1.0

                # Keep track of the sky value that should be subtracted from this chip
                # Read in value from image header, in case user has already
                # determined the sky level
                #
                # .computedSky:   value to be applied by the
                #                 adrizzle/ablot steps.
                # .subtractedSky: value already (or will be by adrizzle/ablot)
                #                 subtracted from the image
                #
                if "MDRIZSKY" in sci_chip.header:
                    subsky = sci_chip.header['MDRIZSKY']
                    log.info('Reading in MDRIZSKY of %s' % subsky)
                    sci_chip.subtractedSky = subsky
                    sci_chip.computedSky = subsky
                else:
                    sci_chip.subtractedSky = 0.0
                    sci_chip.computedSky = None

                sci_chip.darkcurrent = 0.0

                # The following attributes are used when working with sub-arrays
                # and get reference file arrays for auto-generation of IVM masks
                try:
                    sci_chip.ltv1 = sci_chip.header['LTV1'] * -1
                    sci_chip.ltv2 = sci_chip.header['LTV2'] * -1
                except KeyError:
                    sci_chip.ltv1 = 0
                    sci_chip.ltv2 = 0
                if sci_chip.ltv1 < 0:
                    sci_chip.ltv1 = 0
                if sci_chip.ltv2 < 0:
                    sci_chip.ltv2 = 0
                sci_chip.size1 = sci_chip.header['NAXIS1'] + np.round(sci_chip.ltv1)
                sci_chip.size2 = sci_chip.header['NAXIS2'] + np.round(sci_chip.ltv2)
                #sci_chip.image_shape = (sci_chip.size2,sci_chip.size1)
                sci_chip.image_shape = (sci_chip.header['NAXIS2'],sci_chip.header['NAXIS1'])

                # Interpret the array dtype by translating the IRAF BITPIX value
                if sci_chip.header['BITPIX'] in _IRAF_DTYPES_TO_NUMPY:
                    sci_chip.image_dtype = _IRAF_DTYPES_TO_NUMPY[sci_chip.header['BITPIX']]

                if self.inmemory:
                    # read image data array into memory
                    shape = sci_chip.data.shape

    def setInstrumentParameters(self,instrpars):
        """ Define instrument-specific parameters for use in the code.
            By definition, this definition will need to be overridden by
            methods defined in each instrument's sub-class.
        """
        pass

    def compute_wcslin(self,undistort=True):
        """ Compute the undistorted WCS based solely on the known distortion
            model information associated with the WCS.
        """
        for chip in range(1,self._numchips+1,1):
            sci_chip = self._image[self.scienceExt,chip]
            chip_wcs = sci_chip.wcs.copy()

            if chip_wcs.sip is None or not undistort or chip_wcs.instrument=='DEFAULT':
                chip_wcs.sip = None
                chip_wcs.cpdis1 = None
                chip_wcs.cpdis2 = None
                chip_wcs.det2im = None
                undistort=False

            # compute the undistorted 'natural' plate scale for this chip
            wcslin = distortion.utils.output_wcs([chip_wcs],undistort=undistort)
            sci_chip.wcslin_pscale = wcslin.pscale

    def set_units(self,chip):
        """ Define units for this image.
        """
        # Determine output value of BUNITS
        # and make sure it is not specified as 'ergs/cm...'
        sci_chip = self._image[self.scienceExt,chip]

        _bunit = None
        if 'BUNIT' in sci_chip.header and sci_chip.header['BUNIT'].find('ergs') < 0:
            _bunit = sci_chip.header['BUNIT']
        else:
            _bunit = 'ELECTRONS/S'
        sci_chip._bunit = _bunit
        #
        if '/s' in _bunit.lower():
            _in_units = 'cps'
        else:
            _in_units = 'counts'
        sci_chip.in_units = _in_units


class WCSObject(baseImageObject):
    def __init__(self,filename,suffix='_drz'):
        super().__init__(filename)

        self._image = fits.HDUList()
        self._image.append(fits.PrimaryHDU())

        # Build rootname, but guard against the rootname being given without
        # the '_drz.fits' suffix
        #patt = re.compile(r"_dr[zc]\w*.fits$")
        drz_extn = suffix
        patt = re.compile(r"_dr[zc]")
        m = patt.search(filename)
        if m:
            self._rootname = filename[:m.start()]
            drz_extn = m.group()
        else:
            # Guard against having .fits in the rootname
            if '.fits' in filename:

                self._rootname = filename[:filename.find('.fits')]
                drz_extn = ''
            else:
                self._rootname = filename

        self.outputNames = self._setOutputNames(self._rootname,suffix=drz_extn)
        self.nimages = 1

        self._bunit = 'ELECTRONS/S'
        self.default_wcs = None
        self.final_wcs = None
        self.single_wcs = None

    def restore_wcs(self):
        self.wcs = copy.copy(self.default_wcs)
