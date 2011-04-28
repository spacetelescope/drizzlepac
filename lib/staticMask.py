"""
This class manages the creation of the global static mask which
masks pixels which are negative in the SCI array.
A static mask numpy object gets created for each global
mask needed, one for each chip from each instrument/detector.
Each static mask array has type Int16, and resides in memory.

:Authors:    
    Ivo Busko, Christopher Hanley, Warren Hack, Megan Sosey
:Program:
    staticMask.py
:Purpose:    
    Class that manages the creation of a global static
    mask which is used to mask pixels that are some
    sigma BELOW the mode computed for the image.
"""
from __future__ import division # confidence high

import numpy as np
from pytools import fileutil, teal
import pyfits
from imagestats import ImageStats
import util
import os
import processInput



__taskname__ = "betadrizzle.staticMask"
_step_num_ = 1


def help():
    print getHelpAsString()
    
#help information that TEAL will look for
def getHelpAsString():
    """ 
    return useful help from a file in the script directory called module.help
    """
    helpString = teal.getHelpFileAsString(__taskname__,__file__)

    return helpString


#this is called by the user
def createMask(input=None, static_sig=4.0, group=None, editpars=False, configObj=None, **inputDict):
    """the user can input a list of images if they like to create static masks
       as well as optional values for static_sig and inputDict
       
       the configObj.cfg file will set the defaults and then override them
       with the user options
    """
        
    if input is not None:
        inputDict["static_sig"]=static_sig
        inputDict["group"]=group
        inputDict["updatewcs"]=False
        inputDict["input"]=input
    else:
        print "Please supply an input image\n"
        raise ValueError
      
    #this accounts for a user-called init where config is not defined yet
    configObj = util.getDefaultConfigObj(__taskname__,configObj,inputDict,loadOnly=(not editpars))
    if configObj is None:
        return

    if editpars == False:
        run(configObj)
    
#this is called by the TEAL interface
def run(configObj):

    #now we really just need the imageObject list created for the dataset
    filelist,output,ivmlist,oldasndict=processInput.processFilenames(configObj['input'],None)

    imageObjList=processInput.createImageObjectList(filelist,instrpars={},group=configObj['group'])  
    createStaticMask(imageObjList,configObj)


#this is the workhorse function called by MultiDrizzle
def createStaticMask(imageObjectList=[],configObj=None,procSteps=None):
    if procSteps is not None:
        procSteps.addStep('Static Mask')
        
    step_name = util.getSectionName(configObj,_step_num_)
    
    if not configObj[step_name]['static']:
        print 'Static Mask step not performed.'
        procSteps.endStep('Static Mask')
        return
    
    if (not isinstance(imageObjectList,list) or (len(imageObjectList) ==0)):
        print "Invalid image object list given to static mask"
        return ValueError
    
    print "\nUSER INPUT PARAMETERS for Static Mask Step:"
    util.printParams(configObj[step_name])        

    #create a static mask object
    myMask=staticMask(configObj)
    
    for image in imageObjectList:
        myMask.addMember(image)
        
    #save the masks to disk for later access  
    myMask.saveToFile()
    myMask.close()

    if procSteps is not None:
        procSteps.endStep('Static Mask')

def constructFilename(signature):
    """construct an output filename for the given signature::
    
         signature=[instr+detector,(nx,ny),detnum]
         
    The signature is in the image object. 
    """
    filename=signature[0]+"_"+str(signature[1][0])+"x"+str(signature[1][1])+"_"+str(signature[2])+"_staticMask.fits"
    return filename        

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
        
        self.masklist={}   
        self.step_name=util.getSectionName(configObj,_step_num_)    
        if configObj is not None:
            self.static_sig = configObj[self.step_name]['static_sig']
        else:
            self.static_sig = 4. # define a reasonable number
            print 'WARNING:  Using default of 4. for static mask sigma.'
                               
            
    def addMember(self, imagePtr=None):
        """
        Combines the input image with the static mask that
        has the same signature.  

        Parameters
        ----------
        imagePtr: object
            An imageObject reference

        Notes
        -----
        The signature parameter consists of the tuple::
        
            (instrument/detector, (nx,ny), chip_id)
       
        The signature is defined in the image object for each chip
        
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
            
            print('  mode = %9f;   rms = %7f;   static_sig = %0.2f')  %  (mode,rms,self.static_sig)

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
             
        filename=constructFilename(signature)

        if(fileutil.checkFileExists(filename)):
            return filename
        else:
            print "\nmMask file for ",str(signature)," does not exist on disk"
            return None
            
    def getMaskname(self,chipid):
        """construct an output filename for the given signature::
        
             signature=[instr+detector,(nx,ny),detnum]
             
        The signature is in the image object and the
        name of the static mask file is saved as sci_chip.outputNames["staticMask"]
        """
        
        return self._image[chipid].outputNames["staticMask"]    
    
       
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
            filename=constructFilename(key)
            
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
               
                   
           
