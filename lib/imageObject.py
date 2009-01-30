#!/usr/bin/env python
"""
A class which makes image objects for 
each input filename

"""

import sys
import util
from pytools import fileutil

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
            print "\Unable to open file:",filename
            raise IOError
            

        #populate the global attributes which are good for all the chips in the file
        self._instrument=self._image[0].header["INSTRUME"]
        self.scienceExt= 'SCI' # the extension the science image is stored in
        self.maskExt='DQ' #the extension with the mask image in it
        self._filename=self._image[0].header["filename"] #can we make this unchangeable?
        
        #assuming all the chips have the same dimensions in the file
        self._naxis1=self._image[self.scienceExt,1].header["NAXIS1"]
        self._naxis2=self._image[self.scienceExt,1].header["NAXIS2"]
        
        
        #this is the number of science chips to be processed in the file
        self._numchips=self._countEXT(extname=self.scienceExt)
        
        #get the rootnames for the chip and add output filename information
        for chip in range(1,self._numchips+1,1):
            self._assignRootname(chip)
            self._image[self.scienceExt,chip].outputNames=util.setOutputNames(self._image[self.scienceExt,chip].rootname) #this is a dictionary
           
    def _assignRootname(self, chip):
        """assign a unique rootname for the image based in the expname"""
        extname=self._image[self.scienceExt,chip].header["EXTNAME"].lower()
        extver=self._image[self.scienceExt,chip].header["EXTVER"]
        expname=self._image[self.scienceExt,chip].header["EXPNAME"]
        self._image[self.scienceExt,chip].rootname=expname + "_" + extname + str(extver)
        
        
    def getData(self,exten=None):
        """return just the specified data extension """
        return fileutil.getExtn(self._image,extn=exten).data
        
    def getHeader(self,exten=None):
        """return just the specified header extension"""
        return fileutil.getExtn(self._image,extn=exten).header
        
        
    def __getitem__(self,exten):
        """overload  getitem to return the data and header"""
        return fileutil.getExtn(self._image,extn=exten)
    
    def __setitem__(self,kw,value):
        """overload setitem to update information, not sure this is right yet"""
        self._image.header.update[kw] = value
    
    def __cmp__(self, other):
        """overload the comparison operator???
            just to check the filename of the object
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
        
         
    def _countEXT(self,extname='SCI'):

        """
            count the number of extensions in the file
            with the given name (EXTNAME)
        """

        _sciext="SCI"
        count=0
        nextend=self._image[0].header["NEXTEND"]
        
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

