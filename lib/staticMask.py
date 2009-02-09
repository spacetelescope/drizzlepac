"""
   Authors:    Ivo Busko, Christopher Hanley, Warren Hack, Megan Sosey
   Program:    staticMask.py
   Purpose:    Class that manages the creation of a global static
               mask which is used to mask pixels that are some
               sigma BELOW the mode computed for the image.

    This class manages the creation of the global static mask which
    masks pixels which are negative in the SCI array.
    A static mask numarray object gets created for each global
    mask needed, one for each chip from each instrument/detector.
    Each static mask array has type Int16, and resides in memory.

"""

import numpy as np
from pytools import fileutil
import pyfits

class staticMask:
    """
    This class manages the creation of the global static mask which
    masks pixels which are negative in the SCI array.
    A static mask numarray object gets created for each global
    mask needed, one for each chip from each instrument/detector.
    Each static mask array has type Int16, and resides in memory.

    The parameter 'badval' defaults to 64 and represents the
    DQ value used to mark pixels flagged as bad by this mask.

    The parameter 'goodval' defaults to 1.0 and represents the
    pixel value of the good pixels.


    """
    
    def __init__ (self, chipImage=None, configObj={}, saveFiles=True): 

        # For now, we don't use badval. It is supposed to
        # be used to flag back the DQ array of the input
        # images. This may lead to confusion when running
        # the task repeatedly over a set of images, whenever
        # additional images are included in the set each
        # time.
        #
        # chipImage is a pointer back to the header+plus data of the chip in imageObject
        # signature is created in the imageObject class
        #
        
        if (chipImage == None):
            print "No image data supplied"
        else:
            self.signature=chipImage.signature    

        self.parameters=self._setDefaults(configObj)
        self.maskPtr=None #points back to the imageObject
        self.saveFiles=saveFiles

    def addMember(self):
        """
        Combines the input image with the static mask that
        has the same signature.  The signature parameter
        consists of the tuple:
        (instrument/detector, (nx,ny), chip_id)
        
        signature is defined in the image object now
        """
        # If this is a new signature, create a new Static Mask file
        if not self.masklist.has_key(self.signature):
            self.masklist[self.signature] = self._buildMaskArray(self.signature)

        # Operate on input image DQ array to flag 'bad' pixels in the
        # global static mask
        stats = ImageStats(chipImage.data,nclip=3,fields='mode')
        mode = stats.mode
        rms  = stats.stddev
        del stats

        print('  mode = %9f;   rms = %7f')  %  (mode,rms)
        #
        # The scale value (3.0) could potentially become a useful 
        # user settable parameter.  Consider for future revisions.
        # 29-April-2004 WJH/CJH/AMK
        #
        sky_rms_diff = mode - (self.static_sig*rms)

        np.bitwise_and(self.masklist[self.signature],np.logical_not(np.less( sci_arr, sky_rms_diff)),self.masklist[self.signature])

    def _buildMaskArray(self):
        """ Creates numarray array for static mask array signature. """
        return np.ones(self.signature[1],dtype=np.int16)

    def getMask(self):
        """ Returns the appropriate StaticMask array for the image. """
        if self.masklist.has_key(self.signature):
            mask =  self.masklist[self.signature]
        else:
            mask = None
        return mask

    def delete(self):
        """ Deletes all static mask objects. """

        for key in self.masklist.keys():
            self.masklist[key] = None
        self.masklist = {}
        
    def saveToFile(self, filename):
        """ saves the static mask to a file
        
        
        """
        
        #check to see if the file already exists on disk
        if not (fileutil.checkFileExists(filename)):
            #create a new fits image with the mask array and a standard header
            for mask in self.masklist:
                #open a new header and data unit
                newHDU = pyfits.PrimaryHDU()
                newHDU.data = masklist[mask]                
            try:
                newHDU.writeto(filename)
                print "Saving static mask to disk:",filename
                
            except IOError:
                print "Problem saving static mask file: ",filename," to disk!\n"
                raise IOError
                
            

    def _setDefaults(configObj={}):
        """set the default parameters for the class"""
 
        static_sig = 4.0
        static_badval = 64
        static_goodval = 1.0
        masklist = {}
        
        paramDict={"static_sig":static_sig,
                    "static_badval":static_badval,
                    "static_goodval":static_goodval,
                    "masklist":masklist}
                    
        if(len(configObj) != 0):
            for key in configObj:
                paramDict[key]=configObj[key] 
                
        self.masklist=paramDict["masklist"]
                  
        return paramDict
