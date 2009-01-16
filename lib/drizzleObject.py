#!/usr/bin/env python

"""
A class to read in the users list of image files and create objects
which can be passsed on to other classes

Files can be in GEIS or MEF format (but not waiver fits).

Runs some sanity checks on the input files.
If necessary converts files to MEF format (this should not be left to makewcs 
because 'updatewcs' may be False).
Runs makewcs.
Returns the list of created objects which have been
initialized



"""
from pytools import parseinput, fileutil, readgeis, asnutil, irafglob, check_files
import pyfits
import os
import util
from hstwcs import updatewcs
import instrumentData


class drizzleObject
Object:
"""

Create an object which contains the basic
information for multidrizzle including the
steps which have been performed and
a run function to actually "multidrizzle"
all the images together.

This includes the list of image files that
the user supplied, and the steps that
have been run on the images.

a python list of pointers to instrumentData
objects is part of its attribute set.


"""	

	def __init__(self,*args,**keywords):
    
    	"""Check to see what kind of file input was given"""
        self.ivmList=[] #just to open the list, hmm, should each image have an ivm attribute?
   		self.objectList=[] #the list of InputImage object references       
        self.fileList=[] #the list of filenames given by the user, and parsed by the code
		
        #what if they ask for the object themselves and want to supply a dictionary of the options?
        #we still need to get the defaults from the config if no all the params are specified
        self.paramDict=utils.populateParameters() #find and populate the user set params for this step
        
        
        keys=keywords.keys()
        keys.sort()
        
        #Keep track of steps to perform on the object
        self.staticMaskDone=False
        self.skySubtractionDone=False
        self.drizzleSeperateDone=False
        self.medianImageDone=False
        self.blotDone=False
        self.derivCRDone=False
        self.drizFinalDone=False
        
       #set up the list of science image files
        self.fileList, self.outputFilename=getInputAsList(inputFilelist, output=None, ivmlist=None, prodonly=False))
        
        #create the list of inputImage objects
   		for filename in self.fileList:
        	self.objectList.append(instrumentData(filename))
        
    
	def run(self):
    	"""step through all the functions to perform full drizzling """
        
        if !(self.staticMaskDone):
        	staticMask(self.objectList)
       	if !(self.skySubtractionDone):
        	subtractSky(self.objectList)
        if !(self.drizzleSeperateDone):
        	drizzleSeperate(self.objectList)
        if !(self.medianImagedone):
        	medianImage(self.objectList)
        if !self.blotDone && self.medianImage:
        	blot(self.objectList)
        if !(self.derivCRDrone):
        	derivCR(self.objectList)
        if !(self.drizFinalDone):
        	drizFinal(self.objectList)    
         
