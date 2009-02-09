#!/usr/bin/env python
"""
A class which makes image objects for 
each input filename

"""

import sys
import util,wcs_functions
from pytools import fileutil
import buildmask

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
    
    def __init__(self,filename=None):
        
        #filutil open returns a pyfits object
        try:
            self._image=fileutil.openImage(filename,clobber=False,memmap=0)
            
        except IOError:
            print "\nUnable to open file:",filename
            raise IOError
            

        #populate the global attributes which are good for all the chips in the file
        self._instrument=self._image['PRIMARY'].header["INSTRUME"]
        self.scienceExt= 'SCI' # the extension the science image is stored in
        self.maskExt='DQ' #the extension with the mask image in it
        self._filename=self._image['PRIMARY'].header["FILENAME"] 
        self._rootname=self._image['PRIMARY'].header["ROOTNAME"]
        self.outputNames=self._setOutputNames(self._rootname)
         
        #this is the number of science chips to be processed in the file
        self._numchips=self._countEXT(extname=self.scienceExt)

        #assign chip specific information
        for chip in range(1,self._numchips+1,1):
            self._assignRootname(chip)
            self._staticmask=None #this will be replaced with a  pointer to a StaticMask object
            sci_chip = self._image[self.scienceExt,chip]
            
            sci_chip.dqfile,sci_chip.dq_extn = self._find_DQ_extension()               
            sci_chip.dqname = sci_chip.dqfile+'['+sci_chip.dq_extn+','+str(chip)+']'

            # build up HSTWCS object for each chip, which will be necessary for drizzling operations
            sci_chip.wcs=wcs_functions.get_hstwcs(self._filename,self._image,sci_chip.header,self._image['PRIMARY'].header)
            sci_chip.detnum,sci_chip.binned = util.get_detnum(sci_chip.wcs,self._filename,chip)

            #assuming all the chips don't have the same dimensions in the file
            sci_chip._naxis1=sci_chip.header["NAXIS1"]
            sci_chip._naxis2=sci_chip.header["NAXIS2"]            
            sci_chip.outputNames=self._setOutputNames(sci_chip.rootname) #this is a dictionary
            self._assignSignature(chip) #this is used in the static mask, static mask name also defined here, must be done after outputNames

            # record the exptime values for this chip so that it can be
            # easily used to generate the composite value for the final output image
            sci_chip._exptime,sci_chip._expstart,sci_chip._expend = util.get_exptime(sci_chip.header,self._image['PRIMARY'].header)
                        
            # Define mask names as additional entries into outputNames dictionary
            sci_chip.outputNames['driz_mask']=sci_chip.dqrootname+'_final_mask.fits'
            sci_chip.outputNames['single_driz_mask']=sci_chip.outputNames['driz_mask'].replace('final','single')

            #
            # SYNTAX for creating mask files:
            #for final mask
            #util.build_mask(sci_chip.dqfile,sci_chip.outputNames['driz_mask'],sci_chip.detnum,sci_chip.dq_extn,bits_final,chip,self._instrument)
            #for single mask
            #util.build_mask(sci_chip.dqfile,sci_chip.outputNames['single_driz_mask'],sci_chip.detnum,sci_chip.dq_extn,bits_single,chip,self._instrument)
            # static_mask needs to be applied to single_driz_mask upon creation


    def __getitem__(self,exten):
        """overload  getitem to return the data and header"""
        return fileutil.getExtn(self._image,extn=exten)
    
    def __setitem__(self,kw,value):
        """overload setitem to update information, not sure this is right yet"""
        # This operation only updates keyword values in the primary header
        self._image.header.update(kw,value)
    
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
        self._image.info()    
        
    
    def close(self):
        """close the object nicely"""
        self._image.close()  
        #do we want to  del self._image here?     

    def getData(self,exten=None):
        """return just the specified data extension """
        return fileutil.getExtn(self._image,extn=exten).data
        
    def getHeader(self,exten=None):
        """return just the specified header extension"""
        return fileutil.getExtn(self._image,extn=exten).header


    def _assignRootname(self, chip):
        """assign a unique rootname for the image based in the expname"""
        extname=self._image[self.scienceExt,chip].header["EXTNAME"].lower()
        extver=self._image[self.scienceExt,chip].header["EXTVER"]
        expname=self._image[self.scienceExt,chip].header["EXPNAME"].lower()

        # record extension-based name to reflect what extension a mask file corresponds to
        self._image[self.scienceExt,chip].rootname=expname + "_" + extname + str(extver)
        self._image[self.scienceExt,chip].dqrootname=self._rootname + "_" + extname + str(extver)
        self._image[self.scienceExt,chip]._expname=expname


    def _assignDQName(self,chip):
        """assign a unique rootname for the DQ image associated with this chip """
        self._image[self.scienceExt,chip].dqname = self.dqfile+'['+self.dq_extn+','+str(chip)+']'
        
    def _assignSignature(self, chip):
        """assign a unique signature for the image based 
           on the  instrument, detector, chip, and size
           this will be used to uniquely identify the appropriate
           static mask for the image
           
           this also records the filename for the static mask to the outputNames dictionary
           
        """
        instr=self._instrument
        detector=self._image['PRIMARY'].header["DETECTOR"]
        nx=self._image[self.scienceExt,chip]._naxis1
        ny=self._image[self.scienceExt,chip]._naxis2
        detnum = self._image[self.scienceExt,chip].detnum
        
        self._image[self.scienceExt,chip].signature=[instr+detector,(nx,ny),detnum]
        self._image[self.scienceExt,chip].outputNames["staticMask"]=instr+"_"+detector+"_"+str(nx)+"x"+str(ny)+"_"+str(detnum)+"_staticMask.fits"
            

    def _setOutputNames(self,rootname):
        """
        Define the default output filenames for drizzle products,
        these are based on the original rootname of the image 

        filename should be just 1 filename, so call this in a loop
        for chip names contained inside a file

        """

        # Define FITS output filenames for intermediate products
        outFinal = rootname+'_drz.fits'
        outSci = rootname+'_sci.fits'
        outWeight = rootname+'_weight.fits'
        outContext = rootname+'_context.fits'
        outSky = rootname + '_sky.fits'
        blotImage = rootname + '_blot.fits'
        crImage = rootname + '_cr.fits'
        outSingle = rootname+'_single.fits'
        outSWeight = rootname+'_wht.fits'
        outSContext = None

        fnames={'outFinal':outFinal,
                'outSci':outSci,
                'outWeight':outWeight,
                'outContext':outContext,
                'outSingle':outSingle,
                'outSWeight':outSWeight,
                'outSContext':outSContext,
                'blotImage':blotImage,
                'crImage':crImage,
                'outSingle':outSingle,
                'outSWeight':outSWeight,
                'outSContext':outSContext,
                'outSky':outSky}

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
        """return lists of all exptime,expstart,expend values 
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
           
        
         
    def _countEXT(self,extname='SCI'):

        """
            count the number of extensions in the file
            with the given name (EXTNAME)
        """

        _sciext="SCI"
        count=0
        nextend=self._image['PRIMARY'].header["NEXTEND"]

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

