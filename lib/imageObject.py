#!/usr/bin/env python
"""
A class which makes image objects for 
each input filename

"""

import sys
import pytools
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
    only to the specific chip. Does it make sense to
    have an object with attributes which are general 
    and then dictionaries for each chip? Like a list
    of dictionaries maybe?
    
 
    
    """
    
    def __init__(self,filename=None):
        
        #filutil open returns a pyfits object
        #am I using the self syntax correctly here? Or should it just be _image
        try:
            self.image=fileutil.openImage(filename,clobber=False,memmap=0)
            
        except IOError:
            print "\Unable to open file:",filename
            raise IOError
            

        #populate the global attributes which are good for all the chips in the file
        self.instrument=self.image[0].header["INSTRUME"]
        self.scienceExt= 'SCI' # the extension the science image is stored in
        self.filename=self.image[0].header["filename"]
        
        #assuming all the chips have the same dimensions in the file
        self.naxis1=self.image[self.scienceExt,1].header["NAXIS1"]
        self.naxis2=self.image[self.scienceExt,1].header["NAXIS2"]
        
        
        #this is the number of science chips to be processed in the file
        self.numchips=self._countEXT(extname=self.scienceExt)
        
        #get the rootnames for the chip
        for chip in range(1,self.numchips,1):
            self.image[chip].rootname=self.image[self.scienceExt,chip].header["EXPNAME"]
               
            
    def _getHeader(self,extname,extver=1):
        """return the header for the specified extension """
        if(extname==''): #ask for specific information
            print "Please specify a header extension to return"
            raise ValueError
            
        return self.image[extname,extver].header
        
    def _getData(self,extname,extver=1):
        """return the data for the specified extension"""
        if(extname==''): #ask for specific information
            print "Please specify a header extension to return"        
        return self.image[extname,extver].data
    
    
    def __cmp__(self, other):
        """overload the comparison operator??? """
        if isinstance(other,imageObject):
            if (self.filename == other.filename):
                return True            
        return False
        
    
    def close(self):
        """close the object nicely"""
        self.image.close()       
        
         
    def _countEXT(self,extname='SCI'):

        """
            count the number of extensions in the file
            with the given name (EXTNAME)
        """

        _sciext="SCI"
        count=0
        nextend=self.image[0].header["NEXTEND"]
        
        for i in range (1,nextend,1):
            if (self.image[i].header["EXTNAME"] == extname):
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

