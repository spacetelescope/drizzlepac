#!/usr/bin/env python

"""
A class to read in the users list of image files and control
the running of multidrizzle using those object

Files can be in GEIS or MEF format (but not waiver fits).

I'm writing this to help me work through running the
whole package as I'm coding it
"""

from pytools import fileutil
import os
import sky
from staticMask import staticMask
from imageObject import imageObject

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
        self._ivmList=[] #just to open the list 
        self._objectList=[] #pointers to the imageObjects                                              
        self._fileList=[] #the list of filenames given by the user, and parsed by the code    
                         #there should be one imageObject for each file in the list           
        self._drizSepList=[] #this is a list of the singly drizzled images                     
        self._medianImage='' #the name of the output median images                             
        self._blotImlist=[] #the list of blotted image names                                   

        #Keep track of steps to perform on the object                                         
        self.staticMaskDone=False                                                             
        self.skySubtractionDone=False                                                         
        self.drizzleSeparateDone=False                                                        
        self.medianImageDone=False                                                            
        self.blotDone=False                                                                   
        self.derivCRDone=False                                                                
        self.drizFinalDone=False                                                              

        #setup default parameters including user overrides                                    
        self.parameters=self._setDefaults(configObj)                                               
        self.saveFiles=saveFiles
        
        #create the list of inputImage object pointers                                        
        if(len(inputImageList) == 0):                                                         
            print "No input images were specified!"                                           
            return ValueError                                                                 

        else:                                                                                 
            self._fileList=inputImageList                                                     

        #create a static mask object for the image set
        self.staticMask=staticMask(configObj)  #this creates a pointer back to the static mask object
                                               #but no masks have been added to it yet
                                               
        self._populateObjects() #read the images into memory (data is lazy instantiation)            


    def run(self):
        """step through all the functions to perform full drizzling """

        #These can be run on individual images, 
        #they dont have to be in memory together or submitted as a list               
        for imageSet in self._objectList:    
            self.computeStaticMask(imageSet)
            self.computeSky(imageSet)
            self.createDrizSep(imageSet)

            imageSet.close()   #the output images have been saved to larger separate files

        self.computeMedian() #this step needs a list of all the separately drizled images   
        self.createBlotImages()
        self.calcDerivCr()

        self.runFinalDrizzle() #give it the list of images

        print "MultiDrizzle finished!"
        

    def _populateObjects(self):
        """read the images into memory and populate the objectlist """
        if(len(self._objectList) == 0):
            for image in self._fileList:
                imageSet=(imageObject(image))
                self._objectList.append(imageSet)  
        else:
            print "Object list already poplulated\nNot adding new ones"
             

    def computeStaticMask(self,imageSet):
        """run static mask step, where imageSet is a single imageObject"""   
                                                                   
        if (self.doStaticMask and not(self.staticMaskDone)):                                       
            try:                                                                                
                self.staticMask.addMember(imageSet)     
                self.staticMaskDone=True  
            except:                                                                             
                print "Problem occured during static mask step"                                 
                return ValueError     
        
        else:
            print "Step is not properly setup"
            return ValueError
                                          
                                                                              
        if(self.saveFiles):                                                                     
            self.staticMask.saveToFile()                        

    def computeSky(self,imageSet):
        """ run sky subtraction """

        if (self.doSkySubtraction and (not(self.skySubtractionDone))):
            try:
                sky.subtractSky(imageSet,self.parameters, self.saveFiles)
                self.skySubtractionDone=True
            except:
                print "Problem occured during sky subtraction step"
                return ValueError

    def createDrizSep(self,imageSet):
        """ drizzle separate images """
        if (self.doDrizzleSeparate and (not(self.drizzleSeparateDone))):
            try:
                self._drizSepList.append(drizzleSeparate(imageSet, self.parameters, self.saveFiles))
                self.drizzleSeparateDone=True
            except:
                print "Problem running driz separate step"
                return ValueError

    def computeMedian(self):
        """ create a median image from the separately drizzled images """
        
        if(self.drizzleSeparateDone):
            try:
                self.medianImage=mkMedian(self._drizSepList, self.parameters,self.saveFiles)
                self.medianImageDone=True
            except:
                print "Problem running median combinations step"
                return ValueError

    def createBlotImages(self):
        """ create blotted images from the median image """
        
        if (self.doBlot and (not(self.blotDone) and self.medianImageDone)):
            try:
                blot(self.medianImage, self._fileList, self.parameters,self.saveFiles)
                self.blotDone=True
            except:
                print "problem running blot image step"
                return ValueError

    def calcDerivCr(self):
        """ run deriv_cr to look for cosmic rays """

        if (self.doDerivCr and (not(self.derivCRDrone)) ):
            try:
                derivCR(self._objectList,self.parameters,self.saveFiles)
                self.dericCRDone=True
            except:
                print "Problem running deriv cr step"
                return ValueError

    def runFinalDrizzle(self):
        """ run through the final drizzle process """
        if (self.doFinalDrizzle and not(self.drizFinalDone)):
            try:
                drizFinal(self._objectList,self.parameters,self.saveFiles)    
            except:
                print "Problem running final drizzle"
                return ValueError


    def _setDefaults(self,configObj={}):
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
                'dec':'' ,
                'build':True,
                'shiftfile':'' ,
                'staticfile':'' }

        #override defaults        
        if(len(configObj) !=0 ):
            for key in configObj:
                params[key]=configObj[key]

        return params

