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

old params
-----------
skysub		'Perform sky subtraction?'
skywidth	'Bin width for sampling sky statistics (in sigma)'
skystat	 	'Sky correction statistics parameter'
skylower	'Lower limit of usable data for sky (always in electrons)'
skyupper	'Upper limit of usable data for sky (always in electrons)'
skyclip		'Number of clipping iterations'
skylsigma	'Lower side clipping factor (in sigma)'
skyusigma	'Upper side clipping factor (in sigma)'

new params
----------
skyuser		'KEYWORD indicating a sky subtraction value if done by user.
group		'the group number for the fits image'



"""

import util
from pytools import cfgpars

def subtractSky(inputImageList,configObj):
	            
    #inputImageList here is assumed to be a python list of filenames
    if len(inputImageList) == 0:
    	print "silly user, I need images to work on"
        return ValueError
             
   # Print out the parameters provided by the user
    print "USER INPUT PARAMETERS for SKY SUBTRACTION:"
	util.printParams(configObj)        
        
        
    """ Processes sky in input images."""
    
    # User Subtraction Case, User has done own sky subtraction, we use 
	# image header value for subtractedsky value
    if configObj["skyuser"] != '':
        print "User has done their own sky subtraction, updating MDRIZSKY with supplied value..."
        for image in inputImageList:
            if int(configObj['group']) == 1:
                try:
                    imageHandle = fileutil.openImage(image,mode='update',memmap=0)
                    userSkyValue = imageHandle[0].header[configObj['skyuser']]
                    imageHandle.close()
                except:
                    print "**************************************************************"
                    print "*"
                    print "*  Cannot find keyword ",skypars['skyuser']," in ",configObj['image'].datafile," to update"
                    print "*"
                    print "**************************************************************\n\n\n"
                    raise KeyError
            image.setSubtractedSky(userSkyValue)
            else:
            	print "Sky updates only implemented for group=1 right now"
                raise KeyError
                
    elif (skysub):
        # Compute our own sky values and subtract them from the image copies.
        print "Subtracting sky..."
        imageMinDict = {}
        currentImageName = "no match"
        for p in self.assoc.parlist:
            configObj['image'].computeSky(skypars)
            computedImageSky = configObj['image'].getComputedSky()
            if (configObj['rootname'] != currentImageName):
                currentMinSky = computedImageSky
                currentImageName = configObj['rootname']
                imageMinDict[currentImageName]=currentMinSky
            else:
                if (computedImageSky < imageMinDict[configObj['rootname']]):
                    imageMinDict[configObj['rootname']]=computedImageSky

        for p in self.assoc.parlist:
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
                configObj['image'].setSubtractedSky(configObj['image'].getComputedSky())
            configObj['image'].subtractSky()

    else:
        # Default Case, sky subtraction is turned off.  No sky subtraction done to image.
        print "No sky subtraction requested, MDRIZSKY set to a value of 0..."
        for p in self.assoc.parlist:
            configObj['image'].setSubtractedSky(0)

    #Update the MDRIZSKY KEYWORD
    for p in self.assoc.parlist:
        if int(configObj['group']) == 1:
            configObj['image'].updateMDRIZSKY(configObj['orig_filename'])
            #configObj['image'].updateMDRIZSKY(configObj['image'].datafile)
