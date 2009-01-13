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


class imageObject:
"""
create a list of image objects from input
given by the user which can be:

- a python list of files
- a comma separated string of filenames (including wild card characters)
- an association table
- an @file (can have a second column with names of ivm files)

a python list of pointers to objects is returned (may contain only 1 object)
each of the objects is an InputImage

This image object class contains all the parameters that pertain
to all the input images to make the final output image
"""	

	def __init__(self, inputFilelist,outputFilename):
    
    	"""Check to see what kind of file input was given"""
		self.ivmList=[] #just to open the list, hmm, should each image have an ivm attribute?
        self.objectList=[] #the list of InputImage object references       
		self.fileList=[] #the list of filenames given by the user
        self.outputFilename=outputFilename
		utils.populateParameters(self.name) #find and populate the user set params for this step
        
        #Keep tracks of steps to perform on the object
        self.staticMaskDone=False
        self.skySubtractionDone=False
        self.drizzleSeperateDone=False
        self.medianImageDone=False
        self.blotDone=False
        self.derivCRDone=False
        self.drizFinalDone=False
        
       #set up the list of science image files
        if(isSingleFile(inputFilelist)):  #a single filename was given
        	fileList=.append(inputFilelist)
        else:
        	self.fileList, self.outputFilename=getInputAsList(inputFilelist, output=None, ivmlist=None, prodonly=False))
        
        #create the list of inputImage objects
   		for filename in self.fileList:
        	self.objectList.append(InputImage(filename))
        
    
	def run(self):
    	"""step through all the functions to perform full drizzling """
        
        if !(staticMaskDone):
        	staticMask(self.objectList)
       	if !(skySubtractionDone):
        	skySubtract(self.objectList)
        if !(drizzleSeperateDone):
        	drizzleSeperate(self.objectList)
        if !(medianImagedone):
        	medianImage(self.objectList)
        if !(blotDone):
        	blot(self.objectList)
        if !(derivCRDrone):
        	derivCR(self.objectList)
        if !(drizFinalDone):
        	drizFinal(self.objectList)    
         
