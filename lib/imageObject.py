#!/usr/bin/env python
"""
A class which makes image objects for 
each input filename

"""

import sys,copy
from pytools import fileutil
import pyfits
import util,wcs_functions
import buildmask
import numpy as np

# Translation table for any image that does not use the DQ extension of the MEF
# for the DQ array.
DQ_EXTNS = {'WFPC2':{'c0h':'sdq','c0f':'sci'}}

__version__ = '0.1dev1'

class baseImageObject:
    def __init__(self,filename):

        self.scienceExt= "SCI" # the extension the science image is stored in
        self.maskExt="DQ" #the extension with the mask image in it
        self._filename = filename

        self._image = None
        self._instrument=None
        self._rootname=None
        self.outputNames={}
        self.outputValues = {}
         
        #this is the number of science chips to be processed in the file
        self._numchips=1
        self._nextend=0
        

    def __getitem__(self,exten):
        """overload  getitem to return the data and header
            these only work on the HDU list already in memory
            once the data has been zero's in self._image you should
            use getData or getHeader to re-read the file
        """
        return fileutil.getExtn(self._image,extn=exten)
    
    
    def __cmp__(self, other):
        """overload the comparison operator
            just to check the filename of the object?
         """
        if isinstance(other,imageObject):
            if (self._filename == other._filename):
                return True            
        return False

    def _isNotValid(self, par1, par2):
        """ Method used to determine if a value or keyword is supplied as 
            input for instrument specific parameters.
        """
        if (par1 == None or par1 == '') and (par2 == None or par2 == ''):
            return True
        else:
            return False
    
    def info(self):
        """return fits information on the _image"""
        #if the file hasn't been closed yet then we can
        #use the pyfits info which looks at the extensions
        if(self._isSimpleFits):
            print self._filename," is a simple fits image"
        else:
            self._image.info()    
 
    def close(self):
        """close the object nicely
           and release all the data arrays from memory
           YOU CANT GET IT BACK, the pointers and data are gone
           so use the getData method to get the data array
           returned for future use. You can use putData to 
           reattach a new data array to the imageObject
        """
        self._image.close()  #calls pyfits.close()
        
        #we actuallly want to make sure that all the
        #data extensions have been closed and deleted
        #since we could have the DQ,ERR and others read in
        #at this point, but I'd like there to be something
        #valid there afterwards that I can play with
        
        if not self._isSimpleFits: 
            for ext in range(1,self._nextend+1,1):
                #use the datatype for the extension
                dtype=self.getNumpyType(self._image[ext].header["BITPIX"])
                self._image[ext].data = np.array(0,dtype=dtype)  #so we dont get io errors on stuff that wasn't read in yet     
        else:            
            self._image.data=np.array(0,dtype=self.getNumpyType(self._image.header["BITPIX"]))
            
            
    def getData(self,exten=None):
        """return just the data array from the specified extension 
        
            fileutil is used instead of pyfits to account for
            non FITS input images. openImage returns a pyfits object
        
        """
        _image=fileutil.openImage(self._filename,clobber=False,memmap=0)
        _data=fileutil.getExtn(_image,extn=exten).data
        _image.close()
        del _image
        return _data
                
    def getHeader(self,exten=None):
        """return just the specified header extension
           
        fileutil is used instead of pyfits to account for
        non FITS input images. openImage returns a pyfits object        
        """
        _image=fileutil.openImage(self._filename,clobber=False,memmap=0)
        _header=fileutil.getExtn(_image,extn=exten).header
        _image.close()
        del _image
        return _header

    def putData(self,data=None,exten=None):
        """Now that we are removing the data from the object to save memory,
            we need something that cleanly puts the data array back into
            the object so that we can write out everything together  using
            something like pyfits.writeto....this method is an attempt to
            make sure that when you add an array back to the .data section
            of the hdu it still matches the header information for that
            section ( ie. update the bitpix to reflect the datatype of the
            array you are adding). The other header stuff is  up to you to verify...
            
            data should be the data array
            exten is where you want to stick it, either extension number or
                a string like 'sci,1'
        """
        if (data == None):
            print "No data supplied"
        else:   
        
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
                
            if(_extnum == None):
                print "no extension number found"
                return ValueError
                
            iraf={'float64':-64,'float32':-32,'uint8':8,'int16':16,'int32':32}
                    
            #update the bitpix to the current datatype, this aint fancy and ignores bscale
            self._image[_extnum].header["BITPIX"]=iraf[data.dtype.name]
            self._image[_extnum].data=data

    def getAllData(self,extname=None,exclude=None):
        """ this function is meant to make it easier to attach ALL the data
        extensions of the image object so that we can write out copies of the
        original image nicer.
        
        if no extname is given, the it retrieves all data from the original
        file and attaches it. Otherwise, give the name of the extensions
        you want and all of those will be restored.
        
        ok, I added another option. If you want to get all the data
        extensions EXCEPT a particular one, leave extname=NONE and
        set exclude=EXTNAME. This is helpfull cause you might not know
        all the extnames the image has, this will find out and exclude
        the one you do not want overwritten.
        """
        
        extensions = self._findExtnames(extname=extname,exclude=exclude)
                   
        for i in range(1,self._nextend+1,1):
            if (self._image[i].extname in extensions) and self._image[i].group_member:
                self._image[i].data=self.getData(self._image[i].extname + ','+str(self._image[i].extver))

    def returnAllChips(self,extname=None,exclude=None):
        """ Returns a list containing all the chips which match the extname given
            minus those specified for exclusion (if any). 
        """
        extensions = self._findExtnames(extname=extname,exclude=exclude)
        chiplist = []
        for i in range(1,self._nextend+1,1):
            if (self._image[i].extname in extensions) and self._image[i].group_member:
                chiplist.append(self._image[i])
        return chiplist
        
    def _findExtnames(self,extname=None,exclude=None):
        """ This method builds a list of all extensions which have 'EXTNAME'==extname
            and do not include any extensions with 'EXTNAME'==exclude, if any are 
            specified for exclusion at all.
        """
        #make a list of the available extension names for the object
        extensions=[]
        if extname != None:
            if not isinstance(extname,list): extname=[extname]
            for extn in extname:
                extensions.append(extn.upper())
        else:
        #restore all the extensions data from the original file, be careful here
        #if you've altered data in memory you want to keep!
            for i in range(1,self._nextend+1,1):
                if self._image[i].extname.upper() not in extensions:
                    extensions.append(self._image[i].extname)
        #remove this extension from the list
        if exclude != None:
            exclude.upper()
            if exclude in extensions:
                newExt=[]
                for item in extensions:
                    if item != exclude:
                        newExt.append(item)
            extensions=newExt
            del newExt
        return extensions
    
    def findExtNum(self,extname=None,extver=1):
        """find the extension number of the give extname and extver"""      
        extnum=None
        _extname=extname.upper()
         
        if not self._isSimpleFits:
            for ext in range(1,self._nextend+1,1):
                if (self._image[ext].extname == _extname):
                    if (self._image[ext].extver == extver):
                        extnum=self._image[ext].extnum
        else:
            print "Image is simple fits"
            
        return extnum        
        
    def _assignRootname(self, chip):
        """assign a unique rootname for the image based in the expname"""
        
        extname=self._image[self.scienceExt,chip].header["EXTNAME"].lower()
        extver=self._image[self.scienceExt,chip].header["EXTVER"]
        expname=self._image[self.scienceExt,chip].header["EXPNAME"].lower()

        # record extension-based name to reflect what extension a mask file corresponds to
        self._image[self.scienceExt,chip].rootname=expname + "_" + extname + str(extver)
        self._image[self.scienceExt,chip].sciname=self._filename + "[" + extname +","+str(extver)+"]"
        self._image[self.scienceExt,chip].dqrootname=self._rootname + "_" + extname + str(extver)
        # Needed to keep EXPTIMEs associated properly (1 EXPTIME for all chips)
        self._image[self.scienceExt,chip]._expname=expname
        self._image[self.scienceExt,chip]._chip =chip
        

    def _setOutputNames(self,rootname):
        """
        Define the default output filenames for drizzle products,
        these are based on the original rootname of the image 

        filename should be just 1 filename, so call this in a loop
        for chip names contained inside a file

        """
        # Define FITS output filenames for intermediate products
        
        # Build names based on final DRIZZLE output name
        # where 'output' normally would have been created 
        #   by 'process_input()'
        #

        outFinal = rootname+'_drz.fits'
        outSci = rootname+'_drz_sci.fits'
        outWeight = rootname+'_drz_weight.fits'
        outContext = rootname+'_drz_context.fits'
        outMedian = rootname+'_med.fits'
        
        # Build names based on input name
        indx = self._filename.find('.fits')
        origFilename = self._filename[:indx]+'_OrIg.fits'
        outSky = rootname + '_sky.fits'
        outSingle = rootname+'_single_sci.fits'
        outSWeight = rootname+'_single_wht.fits'
        
        # Build outputNames dictionary
        fnames={
            'origFilename':origFilename,
            'outFinal':outFinal,
            'outMedian':outMedian,
            'outSci':outSci,
            'outWeight':outWeight,
            'outContext':outContext,
            'outSingle':outSingle,
            'outSWeight':outSWeight,
            'outSContext':None,
            'outSky':outSky,
            'ivmFile':None}
        

        return fnames

    def _setChipOutputNames(self,rootname,chip):
        blotImage = rootname + '_blt.fits'
        crmaskImage = rootname + '_crmask.fits'
        crcorImage = rootname + '_cor.fits'


        # Start with global names
        fnames = self.outputNames

        # Now add chip-specific entries
        fnames['blotImage'] = blotImage
        fnames['crcorImage'] = crcorImage
        fnames['crmaskImage'] = crmaskImage
        sci_chip = self._image[self.scienceExt,chip]
        # Define mask names as additional entries into outputNames dictionary
        fnames['finalMask']=sci_chip.dqrootname+'_final_mask.fits' # used by final_drizzle
        fnames['singleDrizMask']=fnames['finalMask'].replace('final','single')
        fnames['staticMask']=None
        
        # Add the following entries for use in creating outputImage object
        fnames['data'] = sci_chip.sciname
        return fnames

    def updateOutputValues(self,output_wcs):
        """Copy info from output WCSObject into outputnames for each chip
           for use in creating outputimage object. 
        """
        
        outputvals = self.outputValues
        
        outputvals['output'] = output_wcs.outputNames['outFinal']
        outputvals['outnx'] = output_wcs.wcs.naxis1
        outputvals['outny'] = output_wcs.wcs.naxis2
        outputvals['texptime'] = output_wcs._exptime
        outputvals['texpstart'] = output_wcs._expstart
        outputvals['texpend'] = output_wcs._expend
        outputvals['nimages'] = output_wcs.nimages
        # Required for blot?
        outputvals['scale'] = output_wcs.wcs.pscale / self._image[self.scienceExt,1].wcs.pscale
        
        outnames = self.outputNames
        outnames['outMedian'] = output_wcs.outputNames['outMedian']
        outnames['outFinal'] = output_wcs.outputNames['outFinal']
        outnames['outSci'] = output_wcs.outputNames['outSci']
        outnames['outWeight'] = output_wcs.outputNames['outWeight']
        outnames['outContext'] = output_wcs.outputNames['outContext']
        
        
    def _find_DQ_extension(self):
        ''' Return the suffix for the data quality extension and the name of the file
            which that DQ extension should be read from.
        '''
        dqfile = None
        for hdu in self._image:
            # Look for DQ extension in input file
            if hdu.header.has_key('extname') and hdu.header['extname'].lower() == self.maskExt.lower():
                dqfile = self._filename
                dq_suffix=self.maskExt
                break
        # This should be moved to a WFPC2-specific version of the imageObject class
        if dqfile == None:
            # Look for additional file with DQ array, primarily for WFPC2 data
            indx = self._filename.find('.fits')
            suffix = self._filename[indx-4:indx]
            dqfile = self._filename.replace(suffix[:3],'_c1')
            dq_suffix = DQ_EXTNS[self._instrument][suffix[1:]]

        return dqfile,dq_suffix
            
    
    def getKeywordList(self,kw):
        """return lists of all attribute values 
           for all active chips in the imageObject
        """
        kwlist = []
        for chip in range(1,self._numchips+1,1):
            sci_chip = self._image[self.scienceExt,chip]
            if sci_chip.group_member:
                kwlist.append(sci_chip.__dict__[kw])
            
        return kwlist

#the following two functions are basically doing the same thing,
#how are they used differently in the code?                
    def getExtensions(self,extname='SCI',section=None):
        ''' Return the list of EXTVER values for extensions with name specified in extname.
        '''
        if section == None:
            numext = 0
            section = []
            for hdu in self._image:
               if hdu.header.has_key('extname') and hdu.header['extname'] == extname:
                    section.append(hdu.header['extver'])
        else:
            if not isinstance(section,list):
                section = [section]

        return section
        
        
         
    def _countEXT(self,extname="SCI"):

        """
            count the number of extensions in the file
            with the given name (EXTNAME)
        """

        count=0 #simple fits image
        
        if (self._image['PRIMARY'].header["EXTEND"]):
            nextend=int(self._image['PRIMARY'].header["NEXTEND"])
            for i in range (1,nextend+1,1):
                self._image[i].extnum=i
                self._image[i].extname=self._image[i].header["EXTNAME"]
                self._image[i].extver=self._image[i].header["EXTVER"]
                
                if (self._image[i].header["EXTNAME"] == extname):
                    count=count+1    
                          
        return count
    
    def getNumpyType(self,irafType):
        """return the corresponding numpy data type"""
        
        iraf={-64:'float64',-32:'float32',8:'uint8',16:'int16',32:'int32'}
        
        return iraf[irafType]
        
    def buildMask(self,chip,bits=0,write=False):
        """ Build masks as specified in the user parameters found in the 
            configObj object.
            
            we should overload this function in the instrument specific
            implementations so that we can add other stuff to the badpixel
            mask? Like vignetting areas and chip boundries in nicmos which
            are camera dependent? these are not defined in the DQ masks, but
            should be masked out to get the best results in multidrizzle
        """
        dqarr = self.getData(exten=self.maskExt+','+str(chip))
        dqmask = self._buildMask(dqarr,bits)
        if write:
            phdu = pyfits.PrimaryHDU(data=dqmask,header=self._image[self.maskExt,chip].header)
            dqmask_name = self._image[self.scienceExt,chip].dqrootname+'_dqmask.fits'
            print 'Writing out DQ mask: ',dqmask_name
            phdu.writeto(dqmask_name)
            del phdu
        del dqarr            
        return dqmask
        """
        ### For WFPC2 Data, build mask files using:
        buildShadowMaskImage(sci_chip.dqfile,sci_chip.detnum,sci_chip.extnum,maskname,bitvalue=bits,binned=sci_chip.binned)
        """

    def _buildMask(self,dqarr,bitvalue):
        """ Builds a bit-mask from an input DQ array and a bitvalue flag"""
        if bitvalue == None:
            return (dqarr * 0.0) + 1.0
        _maskarr = np.bitwise_or(dqarr,np.array([bitvalue]))
        return np.choose(np.greater(_maskarr,bitvalue),(1,0)).astype(np.uint8)

    def updateIVMName(self,ivmname):
        """ Update outputNames for image with user-supplied IVM filename."""
        self.outputNames['ivmFile'] = ivmname

    def set_units(self):
        """ Record the units for this image, both BUNITS from header and 
            in_units as needed internally.
            This method will be defined specifically for each instrument.
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
        if (value != None and value != '')  and (keyword != None and keyword.strip() != ''):
            exceptionMessage = "ERROR: Your input is ambiguous!  Please specify either a value or a keyword.\n  You specifed both " + str(value) + " and " + str(keyword) 
            raise ValueError, exceptionMessage
        elif value != None and value != '':
            return self._averageFromList(value)
        elif keyword != None and keyword.strip() != '':
            return self._averageFromHeader(header, keyword)
        else:
            return None

    def _averageFromHeader(self, header, keyword):
        """ Averages out values taken from header. The keywords where
            to read values from are passed as a comma-separated list.
        """
        _list = ''
        for _kw in keyword.split(','):
            if header.has_key(_kw):
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
    PyFits object with extra attributes
    
    There will be generic keywords which are good for
    the entire image file, and some that might pertain
    only to the specific chip. 
    
    """
    
    def __init__(self,filename,group=None):
        baseImageObject.__init__(self,filename)
        
        #filutil open returns a pyfits object
        try:
            self._image=fileutil.openImage(filename,clobber=False,memmap=0)
            
        except IOError:
            print "\nUnable to open file:",filename
            raise IOError

        #populate the global attributes which are good for all the chips in the file
        self._rootname=self._image['PRIMARY'].header["ROOTNAME"]
        self.outputNames=self._setOutputNames(self._rootname)
        
        self._exptime=self._image["PRIMARY"].header["EXPTIME"]
        if(self._exptime == 0): self._exptime =1. #to avoid divide by zero
         
        #this is the number of science chips to be processed in the file
        self._numchips=self._countEXT(extname=self.scienceExt)

        self.proc_unit = None

        if (self._numchips == 0):
            self._isSimpleFits = True
            self._nextend=0
        else:
            self._isSimpleFits = False
        
        if group not in [None,'']:
            # Only use selected chip(s?)
            group_id = fileutil.parseExtn(str(group))
            if group_id[0] == '':
                # find extname/extver which corresponds to this extension number
                group_extname = self._image[group_id[1]].header['EXTNAME']
                group_extver = self._image[group_id[1]].header['EXTVER']
                self.group = [group_extname,group_extver]
            else:
                self.group = group_id
        else:
            # Use all chips
            self.group = None
        
        if not self._isSimpleFits:
            self._nextend=self._image["PRIMARY"].header["NEXTEND"]

            #assign chip specific information
            for chip in range(1,self._numchips+1,1):
                
                self._assignRootname(chip)
                sci_chip = self._image[self.scienceExt,chip]

                # Set a flag to indicate whether this chip should be included
                # or not, based on user input from the 'group' parameter.
                if self.group is None or (self.group is not None and self.group[1] == chip):
                    sci_chip.group_member = True
                else:
                    sci_chip.group_member = False

                sci_chip.signature = None
                
                sci_chip.dqfile,sci_chip.dq_extn = self._find_DQ_extension()               
                sci_chip.dqname = sci_chip.dqfile+'['+sci_chip.dq_extn+','+str(chip)+']'

                # build up HSTWCS object for each chip, which will be necessary for drizzling operations
                sci_chip.wcs=wcs_functions.get_hstwcs(self._filename,self._image,sci_chip.extnum)
                sci_chip.detnum,sci_chip.binned = util.get_detnum(sci_chip.wcs,self._filename,chip)

                #assuming all the chips don't have the same dimensions in the file
                sci_chip._naxis1=sci_chip.header["NAXIS1"]
                sci_chip._naxis2=sci_chip.header["NAXIS2"]            

                # record the exptime values for this chip so that it can be
                # easily used to generate the composite value for the final output image
                sci_chip._exptime,sci_chip._expstart,sci_chip._expend = util.get_exptime(sci_chip.header,self._image['PRIMARY'].header)
                            
                sci_chip.outputNames=self._setChipOutputNames(sci_chip.rootname,chip).copy() #this is a dictionary
                # Set the units: both bunit and in_units
                self.set_units(chip)

    def setInstrumentParameters(self,instrpars):
        """ Define instrument-specific parameters for use in the code. 
            By definition, this definition will need to be overridden by 
            methods defined in each instrument's sub-class.
        """
        pass
                                    
    def set_units(self,chip):
        """ Define units for this image."""
        # Determine output value of BUNITS
        # and make sure it is not specified as 'ergs/cm...'
        sci_chip = self._image[self.scienceExt,chip]

        _bunit = None
        if sci_chip.header.has_key('BUNIT') and sci_chip.header['BUNIT'].find('ergs') < 0:
            _bunit = sci_chip.header['BUNIT']
        else:
            _bunit = 'ELECTRONS/S'
        sci_chip._bunit = _bunit
        #
        sci_chip.in_units = 'counts'
                            

class WCSObject(baseImageObject):
    def __init__(self,filename,suffix='_drz.fits'):
        baseImageObject.__init__(self,filename)
                
        self._image = pyfits.HDUList()
        self._image.append(pyfits.PrimaryHDU())
        self._rootname = filename[:filename.find(suffix)]
        self.outputNames = self._setOutputNames(self._rootname)
        self.nimages = 1
    
        self._bunit = 'ELECTRONS/S'
        self.default_wcs = None
        self.final_wcs = None
        self.single_wcs = None

    def restore_wcs(self):
        self.wcs = copy.copy(self.default_wcs)

