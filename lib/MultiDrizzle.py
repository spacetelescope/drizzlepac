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
import staticMask
import imageObject

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
        self._objectList=[] #the list of imageObject object references                         
        self._fileList=[] #the list of filenames given by the user, and parsed by the code    
                         #there should be one imageObject for each file in the list           
        self._drizSepList=[] #this is a list of the singly drizzled images                     
        self._medianImage='' #the name of the output median images                             
        self._blotImlist=[] #the list of blotted image names                                   

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
            self._fileList=inputImageList                                                     


    def run(self):
        """step through all the functions to perform full drizzling """

        #These can be run on individual images, 
        #they dont have to  be in memory together        
        for filename in self.fileList:
            imageSet=(imageObject(filename))

            _computeStaticMask(imageSet)
            _computeSky(imageSet)
            _createDrizSep(imageSet)

            imageSet.close()   #the output images have been saved to larger separate files

        _computeMedian() #this step needs a list of all the separately drizled images   
        _createBlotImages()
        _calcDerivCr()

        _runFinalDrizzle() #give it the list of images

        print "MultiDrizzle finished!"
        
        

    def _computeStaticMask(self,imageSet):
        """run static mask step"""   
                                                                   
        if (self.doStaticMask and not(self.staticMaskDone)):                                       
            try:                                                                                
                staticMask(imageSet, self.parameters,self.saveFiles)       
            except:                                                                             
                print "Problem occured during static mask step"                                 
                return ValueError     
                                                                          
        if(self.saveFiles):                                                                     
            imageSet.staticMask.saveToFile(image.outputNames["staticMask"])                        

    def _computeSky(self,imageSet):
        """ run sky subtraction """

        if (self.doSkySubtraction and (not(self.skySubtractionDone))):
            try:
                sky.subtractSky(imageSet,self.parameters, self.saveFiles)
            except:
                print "Problem occured during sky subtraction step"
                return ValueError

    def _createDrizSep(self,imageSet):
        """ drizzle seperate images """
        if (self.doDrizzleSeparate and (not(self.drizzleSeperateDone))):
            try:
                self._drizSepList.append(drizzleSeperate(imageSet, self.parameters, self.saveFiles))
            except:
                print "Problem running driz seperate step"
                return ValueError

    def _computeMedian(self):
        """ create a median image from the separately drizzled images """
        try:
            self.medianImage=mkMedian(self._drizSepList, self.parameters,self.saveFiles)
        except:
            print "Problem running median combinations step"
            return ValueError

    def _createBlotImages(self):
        """ create blotted images from the median image """
        
        if (self.doBlot and (not(self.blotDone) and self.medianImageDone)):
            try:
                blot(self.medianImage, self._filelist, self.parameters,self.saveFiles)
            except:
                print "problem running blot image step"
                return ValueError

    def _calcDerivCr():
        """ run deriv_cr to look for cosmic rays """

        if (self.doDerivCr and (not(self.derivCRDrone)) ):
            try:
                derivCR(self._objectList,self.parameters,self.saveFiles)
            except:
                print "Problem running deriv cr step"
                return ValueError

    def runFinalDrizzle():
        """ run through the final drizzle process """
        if (self.doFinalDrizzle and not(self.drizFinalDone)):
            try:
                drizFinal(self._objectList,self.parameters,self.saveFiles)    
            except:
                print "Problem running final drizzle"
                return ValueError


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
                'dec':'' ,
                'build':True,
                'shiftfile':'' ,
                'staticfile':'' }

        #override defaults        
        if(len(configObj) !=0 ):
            for key in configObj:
                params[key]=configObj[key]

        return params

