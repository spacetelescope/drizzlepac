#!/usr/bin/env python
"""
A class which makes image objects for 
each input filename

"""

import sys
import util,wcs_functions
from pytools import fileutil
import buildmask
import numpy as np

# Translation table for any image that does not use the DQ extension of the MEF
# for the DQ array.
DQ_EXTNS = {'WFPC2':{'c0h':'sdq','c0f':'sci'}}

__version__ = '0.1dev1'

class imageObject():
    """
    This returns an imageObject that contains all the
    necessary information to run the image file through
    any multidrizzle function. It is essentially a 
    PyFits object with extra attributes
    
    There will be generic keywords which are good for
    the entire image file, and some that might pertain
    only to the specific chip. 
    
    """
    
    def __init__(self,filename):
        
        #filutil open returns a pyfits object
        try:
            self._image=fileutil.openImage(filename,clobber=False,memmap=0)
            
        except IOError:
            print "\nUnable to open file:",filename
            raise IOError
            

        #populate the global attributes which are good for all the chips in the file
        self._instrument=self._image['PRIMARY'].header["INSTRUME"]
        self.scienceExt= "SCI" # the extension the science image is stored in
        self.maskExt="DQ" #the extension with the mask image in it
        self._filename = filename
        self._rootname=self._image['PRIMARY'].header["ROOTNAME"]
        self.outputNames=self._setOutputNames(self._rootname)
         
        #this is the number of science chips to be processed in the file
        self._numchips=self._countEXT(extname=self.scienceExt)
        
        if (self._numchips == 0):
            self._isSimpleFits = True
        else:
            self._isSimpleFits = False
            
        self._nextend=0
        
        if not self._isSimpleFits:
            self._nextend=self._image["PRIMARY"].header["NEXTEND"]

            #assign chip specific information
            for chip in range(1,self._numchips+1,1):
                self._assignRootname(chip)
                sci_chip = self._image[self.scienceExt,chip]
                sci_chip._staticmask=None #this will be replaced with a  pointer to a StaticMask object

                sci_chip.dqfile,sci_chip.dq_extn = self._find_DQ_extension()               
                sci_chip.dqname = sci_chip.dqfile+'['+sci_chip.dq_extn+','+str(chip)+']'

                # build up HSTWCS object for each chip, which will be necessary for drizzling operations
                sci_chip.wcs=wcs_functions.get_hstwcs(self._filename,self._image,sci_chip.header,self._image['PRIMARY'].header)
                sci_chip.detnum,sci_chip.binned = util.get_detnum(sci_chip.wcs,self._filename,chip)

                #assuming all the chips don't have the same dimensions in the file
                sci_chip._naxis1=sci_chip.header["NAXIS1"]
                sci_chip._naxis2=sci_chip.header["NAXIS2"]            
                self._assignSignature(chip) #this is used in the static mask, static mask name also defined here, must be done after outputNames

                # record the exptime values for this chip so that it can be
                # easily used to generate the composite value for the final output image
                sci_chip._exptime,sci_chip._expstart,sci_chip._expend = util.get_exptime(sci_chip.header,self._image['PRIMARY'].header)

                sci_chip.outputNames=self._setChipOutputNames(sci_chip.rootname,chip).copy() #this is a dictionary

                # Determine output value of BUNITS
                # and make sure it is not specified as 'ergs/cm...'
                _bunit = None
                if sci_chip.header.has_key('BUNIT') and sci_chip.header['BUNIT'].find('ergs') < 0:
                    _bunit = sci_chip.header['BUNIT']
                else:
                    _bunit = 'ELECTRONS/S'
                sci_chip._bunit = _bunit


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
    
    def info(self):
        """return fits information on the _image"""
        #if the file hasn't been closed yet then we can
        #use the pyfits info which looks at the extensions
        if(self._isSimpleFits):
            print self._filename," is a simple fits image"
        else:
            if(len(self._image[1].data.shape) !=0 ): 
                self._image.info()    
            
            #otherwise, we need to do something else
            else:  
                print "Data sections have already been closed\n"  
    
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
            for ext in range(1,self._nextend,1):
                self._image[1].data = np.array(0)  #so we dont get io errors on stuff that wasn't read in yet        

        else:
            self._image.data=np.array(0)
            
    def getData(self,exten=None):
        """return just the data array from the specified extension 
            this method should be 
        
            fileutil is used instead of pyfits to account for
            non FITS input images. openImage returns a pyfits object
        
        """
        _image=fileutil.openImage(self._filename,clobber=False,memmap=0)
        _data=fileutil.getExtn(_image,extn=exten).data
        
        del _image
        return _data
                
    def getHeader(self,exten=None):
        """return just the specified header extension
           
        fileutil is used instead of pyfits to account for
        non FITS input images. openImage returns a pyfits object        
        """
        image=fileutil.openImage(self._filename,clobber=False,memmap=0)
        header=fileutil.getExtn(self._image,extn=exten).header
        image.close()
        del image
        return header

    def putData(self,exten=None):
        """Now that we are removing the data from the object to save memory,
            we need something that cleanly puts the data array back into
            the object so that we can write out everything together  using
            something like pyfits.writeto....this method is an attempt to
            make sure that when you add an array back to the .data section
            of the hdu it still matches the header information for that
            section ( ie. update the bitpix to reflect the datatype of the
            array you are adding). The other header stuff is  up to you to verify...
            
        """
        

    def _assignRootname(self, chip):
        """assign a unique rootname for the image based in the expname"""
        extname=self._image[self.scienceExt,chip].header["EXTNAME"].lower()
        extver=self._image[self.scienceExt,chip].header["EXTVER"]
        expname=self._image[self.scienceExt,chip].header["EXPNAME"].lower()

        # record extension-based name to reflect what extension a mask file corresponds to
        self._image[self.scienceExt,chip].rootname=expname + "_" + extname + str(extver)
        self._image[self.scienceExt,chip].sciname=self._filename + "[" + extname +","+str(extver)+"]"
        self._image[self.scienceExt,chip].dqrootname=self._rootname + "_" + extname + str(extver)
        self._image[self.scienceExt,chip]._expname=expname
        self._image[self.scienceExt,chip]._chip = chip

        
    def _assignSignature(self, chip):
        """assign a unique signature for the image based 
           on the  instrument, detector, chip, and size
           this will be used to uniquely identify the appropriate
           static mask for the image
           
           this also records the filename for the static mask to the outputNames dictionary
           
        """
        instr=self._instrument
        detector=self._image['PRIMARY'].header["DETECTOR"]
        ny=self._image[self.scienceExt,chip]._naxis1
        nx=self._image[self.scienceExt,chip]._naxis2
        detnum = self._image[self.scienceExt,chip].detnum
        
        self._image[self.scienceExt,chip].signature=(instr+detector,(nx,ny),detnum) #signature is a tuple

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
            'outMedian':outMedian,
            'outFinal':outFinal,
            'outSci':outSci,
            'outWeight':outWeight,
            'outContext':outContext,
            'outSingle':outSingle,
            'outSWeight':outSWeight,
            'outSContext':None,
            'outSky':outSky}
        

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

        # Define mask names as additional entries into outputNames dictionary
        fnames['drizMask']=self._image[self.scienceExt,chip].dqrootname+'_final_mask.fits'
        fnames['singleDrizMask']=fnames['drizMask'].replace('final','single')

        # Add entries to outputNames for use by 'drizzle'
        fnames['inData']=self._image[self.scienceExt,chip].sciname
        fnames['extver']=chip

        return fnames
        
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
           for all chips in the imageObject
        """
        kwlist = []
        for chip in range(1,self._numchips+1,1):
            sci_chip = self._image[self.scienceExt,chip]
            kwlist.append(sci_chip.__dict__[kw])
            
        return kwlist
                
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
            for i in range (1,nextend,1):
                if (self._image[i].header["EXTNAME"] == extname):
                    count=count+1    
            
        return count
    
    def _averageFromHeader(self, header, keyword):
        """ Averages out values taken from header. The keywords from which
            to read values are passed as a comma-separated list.
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

