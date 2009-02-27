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
import util


__taskname__ = "BigBlackBox.staticMask"
_step_num_ = 1


#help information that TEAL will look for
def getHelpAsString():
    return "Static Mask Help will eventually be here"


#this is called by the user
def staticMask(imageList=None,static_sig=None,editpars=False,**inputDict):
    """the user can input a list of images if they like to create static masks
       as well as optional values for static_sig and inputDict
       
       the configObj.cfg file will set the defaults and then override them
       with the user options
    """
    if not isinstance(imageList,list):
        imageList=[imageList]
        
    if(static_sig != None):
        inputDict["static_sig"]=static_sig
        inputDict["input"]=imageList

    #this accounts for a user called init where config is not defined yet

    configObj = util.getDefaultConfigObj(__taskname__,configObj,inputDict,loadOnly=loadOnly(not editpars))
        
    run(configObj,inputDict)
    
#this is called by the TEAL interface
def run(configObj=None):
    imgObjList,outwcs = processInput.setCommonInput(configObj,createOutwcs=False) #outwcs is not neaded here

    createStaticMask(imgObjList,configObj)


#this is the workhorse function
def createStaticMask(imageObjectList=[],configObj=None):

    if (not isinstance(imageObjectList,list) or (len(imageObjectList) ==0)):
        print "Invalid image object list given to static mask"
        return ValueError
    
    #create a static mask object
    staticMask=staticMask.staticMask(configObj)
    
    for image in imageObjectList:
        staticMask.addMember(image)
        
    #save the masks to disk for later access  
    staticMask.saveToFile()
    staticMask.close()

class staticMask:
    """
    This class manages the creation of the global static mask which
    masks pixels that are unwanted in the SCI array.
    A static mask  object gets created for each global
    mask needed, one for each chip from each instrument/detector.
    Each static mask array has type Int16, and resides in memory

    """
    
    def __init__ (self, configObj=None): 

        # the signature is created in the imageObject class
        
        self.static_sig=4. #just a reasonable number
        self.masklist={}        
        self.step_name=util.getSectionName(configObj,_step_num_)    
                               
            
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
        
        print "Computing static mask:\n"
        for chip in range(1,numchips+1,1):
            chipid=imagePtr.scienceExt + ','+ str(chip)
            chipimage=imagePtr.getData(chipid)
            signature=imagePtr[chipid].signature

            # If this is a new signature, create a new Static Mask file which is empty
            # only create a new mask if one doesn't already exist
            if ((not self.masklist.has_key(signature)) or (len(self.masklist) == 0)):
                self.masklist[signature] = self._buildMaskArray(signature)

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

    def getFilename(self,signature):
        """returns the name of the output mask file that
            should reside on disk for the given signature """
             
        filename=self.constructFilename(signature)

        if(fileutil.checkFileExists(filename)):
            return filename
        else:
            print "\nmMask file for ",str(signature)," does not exist on disk"
            return None
            
    def getMaskname(self,chipid):
        """construct an output filename for the given signature
             signature=[instr+detector,(nx,ny),detnum]
             
             the signature is in the image object and the
             name of the static mask file is saved as sci_chip.outputNames["staticMask"]
        """
        
        return self._image[chipid].outputNames["staticMask"]

    def constructFilename(signature):
        """construct an output filename for the given signature
             signature=[instr+detector,(nx,ny),detnum]
             
             the signature is in the image object 
        """
        filename=signature[0]+"_"+str(signature[1][0])+"x"+str(signature[1][1])+"_"+str(signature[2])+"_staticMask.fits"
        return filename        
    
    
       
    def close(self):
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
               
                   
           
