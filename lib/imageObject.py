#!/usr/bin/env python
"""
functions associated with HST instruments


"""

import sys
import pytools
import util
from pytools import fileutil

class imageObject:
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
            self._image=fileutil.openImage(filename,clobber=False,memmap=0)
        except:
            print "\Unable to open file:",filename
            raise ValueError
            

        #populate the global attributes
        self.instrument=self._image[0].header["INSTRUME"]
        self.scienceExt= 'SCI' # the extension the science image is stored in
        print self.scienceExt
        print self._image[0].header["NEXTEND"]
        
        #this is the number of science chips to be processed in the file
        self.numchips=self._countEXT(extname=self.scienceExt)
        
        #get the rootnames for the chip
        for chip in numchips:
            thisImage[chip].rootname=self._image[self.scienceExt,chip].header["EXPNAME"]
                
        return self._image
        
    
    #close the object nicely, this should be calling pyfits.close() I think
    def close(self):
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
            if (image[i].header["EXTNAME"] == extname):
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

