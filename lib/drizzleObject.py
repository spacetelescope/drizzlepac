#!/usr/bin/env python

"""
A class to read in the users list of image files and control
the running of multidrizzle using those object

Files can be in GEIS or MEF format (but not waiver fits).

"""
from pytools import parseinput, fileutil, irafglob
import os
import util


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

    def __init__(self,inputImageList=[],configObj={}):

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

        if (self.doStaticMask and (not self.staticMaskDone)):
            staticMask(self.objectList)
      
        if (self.doSkySubtraction and (not self.skySubtractionDone)):
            subtractSky(self.objectList)
            
        if (self.doDrizzleSeparate and (not self.drizSeperateDone)):
        	drizSeparate(self.objectList)

        if (self.doMakeMedian and (not self.medianImagedone)):
        	medianImage(self.objectList)

        if (self.doBlot and not self.blotDone and self.medianImageDone):
        	blot(self.objectList)

        if (self.doDerivCr and (not self.derivCRDrone)):
        	derivCR(self.objectList)
       
        if (self.doFinalDrizzle and (not self.drizFinalDone)):
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
                'coeffs':'header',
                'output':'',
                'mdriztab':False,
                'refimage':'',
                'runfile':'',
                'workinplace':False,
                'updatewcs':True,
                'proc_unit':"native", 
                'coeffs':'header',
                'context':False,
                'clean':False,
                'group':'',
                'ra':'', 
                'dec':'', 
                'build':True,
                'shiftfile':'', 
                'staticfile':'' }
        
        #override defaults        
        if(len(configObj) !=0 ):
            for key in configObj:
                params[key]=configObj[key]
                
        return params
