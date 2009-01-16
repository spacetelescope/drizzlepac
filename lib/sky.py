#!/usr/bin/env python
"""

Function for computing and subtracting the backgroud of an image.  The
algorithm employed here uses a sigma clipped median of  each *sci* image in a data file.  
Then the sky value for each detector is compared and the lowest value is 
subtracted from all chips in the detector.  Finally, the MDRIZSKY keyword 
is updated in the header of the input files.

:author: Christopher Hanley
:author: Megan Sosey

inputImageList is a python list of image filename
configObj is there as a placeholder to pass in the full parameter set
	for now, pass in configObj as a dictionary of all the necessary params

These are parameters that the configObj should contain:


params that should be in configobj
---------------------------------------
skyuser		'KEYWORD indicating a sky subtraction value if done by user.
skysub		'Perform sky subtraction?'
skywidth	'Bin width for sampling sky statistics (in sigma)'
skystat	 	'Sky correction statistics parameter'
skylower	'Lower limit of usable data for sky (always in electrons)'
skyupper	'Upper limit of usable data for sky (always in electrons)'
skyclip		'Number of clipping iterations'
skylsigma	'Lower side clipping factor (in sigma)'
skyusigma	'Upper side clipping factor (in sigma)'


stuff I got rid of
-------------------
#can't we do something to get rid of this parameter? I've always hated it
#or do we need it for some reason I haven't figured out yet?
group		'the group number for the fits image' 


"""

import util
import instrumentData
from pytools import cfgpars

def subtractSky(inputImageList,configObj=configObj):

	"""inputImageList is a python list of image filenames
       configObj is a dictionary of additional parameters needed to run the sky subtraction
    """
	            
    #inputImageList here is assumed to be a python list of filenames
    if len(inputImageList) == 0:
    	print "Empty image list given to Sky routine"
        return ValueError
             
   # Print out the parameters provided by the user
    print "USER INPUT PARAMETERS for SKY SUBTRACTION:"
	util.printParams(configObj)        
        
    _computedSky = 0.0
    _subtractedSky = 0.0
    _skyValue=0.0  
        
    """ Processes sky in input images."""
    
    # User Subtraction Case, User has done own sky subtraction, we use 
	# image header value for subtractedsky value
    if configObj["skyuser"] != '':
        print "User has done their own sky subtraction, updating MDRIZSKY with supplied value..."
        for image in inputImageList:
            try:
                imageHandle = fileutil.openImage(image,mode='update',memmap=0)
                _skyValue = imageHandle[0].header[configObj["skyuser"]]
            except:
                print "**************************************************************"
                print "*"
                print "*  Cannot find keyword ",configObj["skyuser"]," in ",image," to update"
                print "*"
                print "**************************************************************\n\n\n"
				imageHandle.close()
                raise KeyError
            
            imageHandle[0].header[configObj["skyuser"]] = _skyValue
            imageHandle.close()
            
                           
    elif (configObj["skysub"]):
        # Compute our own sky values and subtract them from the image copies.
        print "Subtracting sky..."
        imageMinDict = {}
        currentImageName = "no match"
        
        for image in inputImageList:
        	imageRootname=util.findrootname(image)
            _computedSky = computeSky(image, configObj, extn="[sci,1]",memmap=0)
            computedImageSky = _computedSky 
            if (imageRootname != currentImageName):
                currentMinSky = computedImageSky
                currentImageName = imageRootname #cause we are just passing in file names right now
                imageMinDict[currentImageName]=currentMinSky
            else:
                if (computedImageSky < imageMinDict[imageRootname]):
                    imageMinDict[imageRootname]=computedImageSky

        for image in inputImageList:
            # We need to account for the fact that STIS associations contain
            # separate exposures for the same chip within the same file.
            # In those cases, we do not want to use the minimum sky value
            # for the set of chips, but rather use each chip's own value.
            # NOTE: This can be generalized later with changes to PyDrizzle
            #       to provide an attribute that specifies whether each member
            #       associated with file is a separate exposure or not.
            #   WJH/CJH 
            if (configObj['exposure'].header['INSTRUME'] != 'STIS'):
                configObj['image'].setSubtractedSky(imageMinDict[configObj['rootname']])
            else:
                _computedSky=(image,getComputedSky(image))
                subtractSky(image,_computedSky)

    else:
        # Default Case, sky subtraction is turned off.  No sky subtraction done to image.
        print "No sky subtraction requested, MDRIZSKY set to a value of 0..."
        for image in imageFileList:
            setSubtractedSky(image,sky=0.0)


    #Update the MDRIZSKY KEYWORD
    for image in imageFileList:
    	updateMDRIZSKY(image, skyValue=_skyValue,memmap=0)
            
           
            
    def computeSky(imageName,  configSky, extn, memmap=0):

        """ Compute the sky value based upon the sci array of the chip
        	 
        
        """

        # Open input image and get pointer to SCI data
        #
 		#
        try:
            _handle = fileutil.openImage(imageName,mode='update',memmap=0)
            _sciext = fileutil.getExtn(_handle) #will grab the first extension with data if none specified
        except:
            raise IOError, "Unable to open %s for sky level computation"%imageName

		#this object contains the returned values from the sky stats routine
        _tmp = ImageStats(_sciext.data,
                fields      = skypars['skystat'],
                lower       = skypars['skylower'],
                upper       = skypars['skyupper'],
                nclip       = skypars['skyclip'],
                lsig        = skypars['skylsigma'],
                usig        = skypars['skyusigma'],
                binwidth    = skypars['skywidth']
                )

        _computedsky = _extractSkyValue(_tmp,skypars['skystat'].lower())
        print "Computed sky value for ",imageName," : ",self._computedsky

        # Close input image filehandle
        _handle.close()
        del _sciext,_handle


    def _extractSkyValue(self,ImageStatsObject,skystat):
        if (skystat =="mode"):
            return ImageStatsObject.mode
        elif (skystat == "mean"):
            return ImageStatsObject.mean
        else:
            return ImageStatsObject.median


    def _subtractSky(self):
        try:
            try:
                _handle = fileutil.openImage(self.name,mode='update',memmap=self.memmap)
                _sciext = fileutil.getExtn(_handle,extn=self.extn)
                print "%s (computed sky,subtracted sky) : (%f,%f)"%(self.name,self.getComputedSky(),self.getSubtractedSky())
                np.subtract(_sciext.data,self.getSubtractedSky(),_sciext.data)
            except:
                raise IOError, "Unable to open %s for sky subtraction"%self.name
        finally:
            _handle.close()
            del _sciext,_handle


    def updateMDRIZSKY(image, skyValue=0.0,memmap=0):
    """ update the header value of MDRIZSKY with the 
    	subtracted skyValue
    """ 
                
        try:
            _handle = fileutil.openImage(image,mode='update',memmap=0)
        except:
            raise IOError, "Unable to open %s for sky level computation"%image
        try:
            try:
                # Assume MDRIZSKY lives in primary header
                print "Updating MDRIZSKY in %s with %f"%(image,skyValue)
                _handle[0].header['MDRIZSKY'] = skyValue
            except:
                print "Cannot find keyword MDRIZSKY in %s to update"%image
                print "Adding MDRIZSKY keyword to primary header with value %f"%skyValue
                _handle[0].header.update('MDRIZSKY',self.getSubtractedSky(), 
                    comment="Sky value subtracted by Multidrizzle")
        finally:
            _handle.close()
    
