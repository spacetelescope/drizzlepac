#!/usr/bin/env python
"""

    Function for computing and subtracting the backgroud of
    an image.  The algorithm employed here uses a sigma 
    clipped median of  each *sci* image in a data file.   
    Then the sky value for each detector is compared 
    and the lowest value is  subtracted from all chips
    in the detector.  Finally, the MDRIZSKY keyword 
    is updated in the header of the input files.

    :author: Christopher Hanley
    :author: Megan Sosey


    stuff I got rid of
    -------------------
    #can't we do something to get rid of this parameter? I've always hated it
    #or do we need it for some reason I haven't figured out yet?
    group		'the group number for the fits image' 


"""
import util
import instrumentData
from pytools import cfgpars

def subtractSky(configObj={},inputImageList=[], skyuser="", skysub=True, skywidth=0.1,
		skystat="median", skylower="INDEF", skyupper="INDEF", skyclip=5, skysligma=4.,skyusigma=4.,):

    """
    inputImageList is a python list of image filename
    configObj is there as a placeholder to pass in the 
        full parameter set for now, pass in configObj 
        as a dictionary right now

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

    """

    #inputImageList here is assumed to be a python list of filenames
    #though I think this might be part of configObj too
    if len(inputImageList) == 0:
        print "Empty image list given to Sky routine"
        return ValueError

    #make up a dictionary of the task parameter values
    paramDict={"skyuser":skyuser,"skysub":skysub,"skywidth":skywidth,
    			"skystat":skystat,"skylower":skylower,"skyupper":skyupper,
                "skyclip":skyclip,"skylsigma":skylsigma,"skyusigma":skyusigma}


    #if configobj has been provided, then use those parameters instead
    if (len(configObj) != 0):
    # Print out the parameters provided by the user
        print "USER INPUT PARAMETERS for SKY SUBTRACTION:"
        util.printParams(configObj)        

        #now replace the defaults with the user specs
        #this assumes only matching keys, though with this
        #syntax configObj could contain more than the standard keys
        #and they will get appended onto the paramDict
        for key in configObj:
            paramDict[key]=configObj[key]


    """ Processes sky in input images."""
    #this will be the final sky value computed for the exposure                  
    _skyValue=0.0                                                                  

    # User Subtraction Case, User has done own sky subtraction,  
	# so use the image header value for subtractedsky value
    if paramDict["skyuser"] != '':
        print "User has done their own sky subtraction, updating MDRIZSKY with supplied value..."
        for image in inputImageList:


            try:
                imageHandle = fileutil.openImage(image,mode='update',memmap=0)
                _skyValue = imageHandle[0].header[paramDict["skyuser"]]

            except:
                print "**************************************************************"
                print "*"
                print "*  Cannot find keyword ",paramDict["skyuser"]," in ",image," to update"
                print "*"
                print "**************************************************************\n\n\n"
                imageHandle.close()
                raise KeyError

            imageHandle[0].header[paramDict["MDRIZSKY"]] = _skyValue
            imageHandle.close()


    elif (paramDict["skysub"]):
        # Compute our own sky values and subtract them from the image copies.
        # for all instruments other than ACS the minimum sky value from all the
        # science chips in the exposure is used as the reference sky for each chip

        print "Subtracting sky..."

        for image in inputImageList:
            instDat=instrumentData.getBasicInstrData(image)
            numsci=instDat["NUMCHIPS"]
            minSky=[]

            for chip in numsci:
            	myext="[SCI,"+str(chip)+"]"
                _computedSky= computeSky(image, paramDict, extn=myext, memmap=0)
                minSky.append(_computedSky)
           
            _skyValue=min(minSky) #what if the sky is negative?

            #now subtract that value from all the chips in the exposure
            subtractSky(image,_skyValue,numsci)

        for image in inputImageList:
            # We need to account for the fact that STIS associations contain
            # separate exposures for the same chip within the same file.
            # In those cases, we do not want to use the minimum sky value
            # for the set of chips, but rather use each chip's own value.
            # NOTE: This can be generalized later with changes to PyDrizzle
            #       to provide an attribute that specifies whether each member
            #       associated with file is a separate exposure or not.
            #   WJH/CJH 
            if (paramDict['exposure'].header['INSTRUME'] != 'STIS'):
                paramDict['image'].setSubtractedSky(imageMinDict[paramDict['rootname']])
            else:
                _computedSky=(image,getComputedSky(image))


    else:
        # Default Case, sky subtraction is turned off.  No sky subtraction done to image.
        print "No sky subtraction requested, MDRIZSKY set to a value of 0..."
        for image in imageFileList:
            setSubtractedSky(image,sky=_skyValue)


    #Update the MDRIZSKY KEYWORD
    for image in imageFileList:
    	updateMDRIZSKY(image, skyValue=_skyValue,memmap=0)


####Helper functions follow ####### 


    def computeSky(imageName,  configSky, extn, memmap=0):

        """ 
        Compute the sky value based upon the sci array of the chip	         
        """

        # Open input image and get pointer to SCI data
        #
 		#
        try:
            _handle = fileutil.openImage(imageName,mode='update',memmap=0)
            _sciext = fileutil.getExtn(_handle, extn=extn) #will grab the first extension with data if none specified
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
        print "Computed sky value for ",imageName," : ",_computedsky

        # Close input image filehandle
        _handle.close()
        del _sciext,_handle

        return _computedsky



    def _extractSkyValue(ImageStatsObject,skystat):
        if (skystat =="mode"):
            return ImageStatsObject.mode
        elif (skystat == "mean"):
            return ImageStatsObject.mean
        else:
            return ImageStatsObject.median



    def _subtractSky(imageName,skyValue,numsci,memmap):
        """
        subtract the given sky value from each chip in the exposure
        """

        try:
            _handle = fileutil.openImage(imageName,mode='update',memmap=memmap)

            for chip in numsci: 
                myext="[SCI,"+str(chip)+"]"   
                _sciext = fileutil.getExtn(_handle,extn=myext)
                np.subtract(_sciext.data,skyValue,_sciext.data)

        except:
            raise IOError, "Unable to open %s for sky subtraction"%imageName

        finally:
            _handle.close()
            del _sciext,_handle



    def updateMDRIZSKY(image, skyValue=0.0,memmap=0):
        """ update the global header value of MDRIZSKY with the 
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

