"""
   Authors:    Ivo Busko, Christopher Hanley, Warren Hack, Megan Sosey
   Program:    staticMask.py
   Purpose:    Class that manages the creation of a global static
               mask which is used to mask pixels that are some
               sigma BELOW the mode computed for the image.

    This class manages the creation of the global static mask which
    masks pixels which are negative in the SCI array.
    A static mask numpy object gets created for each global
    mask needed, one for each chip from each instrument/detector.
    Each static mask array has type Int16, and resides in memory.

"""

import numpy as np
from pytools import fileutil
import pyfits
from imagestats import ImageStats

class staticMask:
    """
    This class manages the creation of the global static mask which
    masks pixels that are negative in the SCI array.
    A static mask numarray object gets created for each global
    mask needed, one for each chip from each instrument/detector.
    Each static mask array has type Int16, and resides in memory.

    The parameter 'badval' defaults to 64 and represents the
    DQ value used to mark pixels flagged as bad by this mask.

    The parameter 'goodval' defaults to 1.0 and represents the
    pixel value of the good pixels.


    """
    
    def __init__ (self, configObj={}): 

        # For now, we don't use badval. It is supposed to
        # be used to flag back the DQ array of the input
        # images. This may lead to confusion when running
        # the task repeatedly over a set of images, whenever
        # additional images are included in the set each
        # time.
        #
        # the signature is created in the imageObject class
        #
        
        self._setDefaults(configObj)
            
            
    def addMember(self, imagePtr=None):
        """
        Combines the input image with the static mask that
        has the same signature.  The signature parameter
        consists of the tuple:
        (instrument/detector, (nx,ny), chip_id)
       
        signature is defined in the image object for each chip
        
        imagePtr is an imageObject reference
        """
        
        numchips=imagePtr._numchips
        
        print "Computing static masks:\n"
        for chip in range(1,numchips+1,1):
            chipid=imagePtr.scienceExt + ','+ str(chip)
            chipimage=imagePtr.getData(chipid)
            signature=imagePtr[chipid].signature

            # If this is a new signature, create a new Static Mask file which is empty
            # only create a new mask if one doesn't already exist
            if ((not self.masklist.has_key(signature)) or (len(self.masklist) == 0)):
                self.masklist[signature] = self._buildMaskArray(signature)

            # Operate on input image DQ array to flag 'bad' pixels in the
            # global static mask
            stats = ImageStats(chipimage,nclip=3,fields='mode')
            mode = stats.mode
            rms  = stats.stddev
            del stats
            
            print('  mode = %9f;   rms = %7f')  %  (mode,rms)

            sky_rms_diff = mode - (self.static_sig*rms)
            np.bitwise_and(self.masklist[signature],np.logical_not(np.less( chipimage, sky_rms_diff)),self.masklist[signature])
            del chipimage

                
    def _buildMaskArray(self,signature):
        """ Creates empty  numpy array for static mask array signature. """
        return np.ones(signature[1],dtype=np.int16)

    def getMaskArray(self, signature):
        """ Returns the appropriate StaticMask array for the image. """
        if self.masklist.has_key(signature):
            mask =  self.masklist[signature]
        else:
            mask = None
        return mask

    def getMaskName(self,signature):
        """returns the name of the output mask file that
            should reside on disk for the given signature """
             
        filename=self.constructFilename(signature)

        if(fileutil.checkFileExists(filename)):
            return filename
        else:
            print "\nmMask file for ",str(signature)," does not exist on disk"
            return None
            
    def constructFilename(self,signature):
        """construct an output filename for the given signature
             signature=[instr+detector,(nx,ny),detnum]
        """
        filename=signature[0]+"_"+str(signature[1][0])+"x"+str(signature[1][1])+"_"+str(signature[2])+"_staticMask.fits"
        return filename

        
    def delete(self):
        """ Deletes all static mask objects. """

        for key in self.masklist.keys():
            self.masklist[key] = None
        self.masklist = {}
        
    def deleteMask(self,signature):
        """delete just the mask that matches the signature given"""
        if self.masklist.has_key(signature):
            self.masklist[signature] = None
        else:
            print "No matching mask"
        
    def saveToFile(self):
        """ saves the static mask to a file
            it uses the signatures associated with each
            mask to contruct the filename for the output mask image
        """
        
        for key in self.masklist.keys():
            #check to see if the file already exists on disk
            filename=self.constructFilename(key)
            
            if not(fileutil.checkFileExists(filename)):
                #create a new fits image with the mask array and a standard header
                #open a new header and data unit
                newHDU = pyfits.PrimaryHDU()
                newHDU.data = self.masklist[key]     
                           
                try:
                    newHDU.writeto(filename)
                    print "Saving static mask to disk:",filename

                except IOError:
                    print "Problem saving static mask file: ",filename," to disk!\n"
                    raise IOError

    def close(self):
        """close out the static mask cleanly"""
        for mask in self.masklist.keys():
            self.masklist[mask]=0.
        self.masklist={}
                  

    def _setDefaults(self,configObj={}):
        """set the default parameters for the class"""
 
        static_sig = 4.0
        static_badval = 64
        static_goodval = 1.0
        masklist = {}
        
        paramDict={"static_sig":static_sig,
                    "static_badval":static_badval,
                    "static_goodval":static_goodval,
                    "masklist":masklist}
                    
        #if a masklist is passed in, it should follow the convention
        #of using the signatures to define the key, otherwise new ones
        #will be created and used automatically
                    
        if(len(configObj) != 0):
            for key in configObj:
                paramDict[key]=configObj[key] 
                
        self.masklist=paramDict["masklist"]
        self.static_sig=paramDict["static_sig"]
        self.static_badval=paramDict["static_badval"]
        self.static_goodval=paramDict["static_goodval"]
        
