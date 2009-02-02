#!/usr/bin/env python

"""
A class to read in the users list of image files and control
the running of multidrizzle using those object

Files can be in GEIS or MEF format (but not waiver fits).

"""
from pytools import parseinput, fileutil, irafglob
import os
import util
from sky import subtractSky
import staticMask
from optparse import OptionParser

class MultiDrizzle:
"""

Create an object which contains the basic
information for multidrizzle including the
steps which have been performed and
a run function to actually "multidrizzle"
all the images together.

This includes the list of image files that
the user supplied, and the steps that
have been run on the images.

"""	

	def __init__(self,inputImageList=[],configObj={},saveFiles=True):
    
        """ inputImageList is a list of filenames supplied by the user
            configObj are the optional user overrides for the parameters
            savefFiles will write output files for every step
        """
    
    	"""Check to see what kind of file input was given"""
        self.ivmList=[] #just to open the list, hmm, should each image have an ivm attribute?
   		self.objectList=[] #the list of imageObject object references       
        self.fileList=[] #the list of filenames given by the user, and parsed by the code
                         #there should be one imageObject for each file in the list
		
        #Keep track of steps to perform on the object
        self.staticMaskDone=False
        self.skySubtractionDone=False
        self.drizzleSeperateDone=False
        self.medianImageDone=False
        self.blotDone=False
        self.derivCRDone=False
        self.drizFinalDone=False
        
        #setup default parameters including user overrides       
        self.parameters=_setDefaults(configObj)
                
        #create the list of inputImage object pointers
        if(len(inputImageList) == 0):
            print "No input images were specified!"
            return ValueError
        else:
       		for filename in self.fileList:
            	self.objectList.append(imageObject(filename))
        
    
	def run(self):
    	"""step through all the functions to perform full drizzling """
        
        if (self.doStaticMask && !(self.staticMaskDone)):
            for imageSet in self.objectList:
                numchips=imageSet._numchips
                for chip in range(1,numchips+1,1):
                    image=imageSet._image[imageSet.scienceExt,chip]
                	image.staticMask=StaticMask(image)
                    if(saveFiles):
                        image.staticMask.saveToFile(image.outputNames["staticMask"])                    
            
       	if (self.doSkySubtraction && (!(self.skySubtractionDone)):
        	subtractSky(self.objectList)
            
        if (self.doDrizzleSeparate && !(self.drizzleSeperateDone)):
        	drizzleSeperate(self.objectList)
            
        if (self.doMakeMedian && !(self.medianImagedone)):
        	medianImage(self.objectList)
            
        if (self.doBlot && !self.blotDone && self.medianImageDone):
        	blot(self.objectList)
            
        if (self.doDerivCr && !(self.derivCRDrone)):
        	derivCR(self.objectList)
            
        if (self.doFinalDrizzle && !(self.drizFinalDone)):
        	drizFinal(self.objectList)    
         

    def _setDefaults(configObj={}):
        """ set the defaults for the user input section"""

        #what steps shall we perform? I'm turning them on as they are completed
        self.doStaticMask=True
        self.doSkySubtraction=True
        self.doDrizzleSeparate=False
        self.doMakeMedian=False
        self.doBlot=False
        self.doDerivCr=False
        self.doFinalDrizzle=False
        
        params={'output':'',
                'mdriztab':'',
                'refimage':'',
                'runfile':'',
                'workinplace':False,
                'updatewcs':True,
                'proc_unit':'native',
                'coeffs'='header',
                'output=''
                'mdriztab=False
                'refimage=''
                'runfile=''
                'workinplace'=False
                'updatewcs'=True
                'proc_unit'="native" 
                'coeffs'='header'
                'context'=False
                'clean'=False
                'group'=''
                'ra'='' 
                'dec'='' 
                'build'=True
                'shiftfile'='' 
                'staticfile'='' }
        
        #override defaults        
        if(len(configObj) !=0 ):
            for key in configObj:
                params[key]=configObj[key]
                
       return params

if __name__ == __main__:

    #parse the command line options
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-c", "--comments", dest="comments",default="none", type="string",
                      help="File that contains comments to add", metavar="COMMENTS")

    parser.add_option("-i", "--image", dest="image",default="none", type="string",
                      help="FITS image to update", metavar="IMAGE")

    parser.add_option("-o", "--overwrite", dest="overimage",action="store_true",
                      default=0,help="Overwrite original FITS image", metavar="OVERIMAGE")

    parser.add_option("-k", "--keeplog", dest="keeplog", action="store_true",default=0,
                      help="Keep logfile of results", metavar="KEEPLOG")

    parser.add_option("-l", "--logfile", dest="logfile",default="fixOpticalHdr.log",type="string",
                      help="file to store log information to",metavar="LOGFILE")

    parser.add_option("-r", "--report", dest="report", default=0, action="store_true",
                      help="report keywords to log that have N/A values",metavar="REPORT")

    (options, args) = parser.parse_args()
