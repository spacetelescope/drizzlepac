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
import assert

#this is the main function that takes an imageObject and the config obj runs with it
#this is eventually what the user function calls as well

def subtractSky(imageObject,paramDict={}):
    """
    subtract the sky from all the chips in the imagefile that imageObject represents
    
    imageObject contains all the information about the chips in the image file
    configObj is represented as a dict for now, but will prolly be an actual config object
    
    """
   """ Processes sky in input images."""
   
                   
    _skyValue=0.0    #this will be the final sky value computed for the exposure                                                                  

    #just making sure, tricky users and all, these are things that will be used
    #by the sky function so we want them defined at least
    try:
        assert (imageObject._numchips > 0), "invalid value for number of chips"
        assert (imageObject._filename != ''), "image object filename is empty!, doh!"
        assert (imageObject.scienceExt !=''), "image object science extension is empty!"
        assert (imageObject._instrument !=''), "image object instrument name is empty!"
    except:
        raise ValueError
        
    numchips=imageObject._numchips
    sciExt=imageObject.scienceExt
    
    # User Subtraction Case, User has done own sky subtraction,  
	# so use the image header value for subtractedsky value
    
    if paramDict["skyuser"] != '':
        print "User has done their own sky subtraction, updating MDRIZSKY with supplied value..."
       
        for chip in numchips:
            try:
                image=imageObject[sciExt+','+str(chip)]                
                _skyValue = image.header[paramDict["skyuser"]]

            except:
                print "**************************************************************"
                print "*"
                print "*  Cannot find keyword ",paramDict["skyuser"]," in ",imageObject._filename
                print "*"
                print "**************************************************************\n\n\n"
                imageHandle.close()
                raise KeyError

            image.header[paramDict["MDRIZSKY"]] = _skyValue

    elif (paramDict["skysub"]):
        # Compute our own sky values and subtract them from the image copies.
        # for all instruments other than ACS the minimum sky value from all the
        # science chips in the exposure is used as the reference sky for each chip

        print "Computing minimum sky ..."

        image=imageObject[sciExt+','+str(chip)] 
        minSky=[] #store the sky for each chip
         
           
        if ("STIS" in imageObject._instrument):
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
            for chip in numchips:
            	    myext="[SCI,"+str(chip)+"]"
                    image=imageObject[myext]
                    _computedSky= computeSky(image.data, paramDict, memmap=0)
                    minSky.append(_computedSky)

                _skyValue=min(minSky) #what if the sky is negative?

            #now subtract that value from all the chips in the exposure
            subtractSky(image,_skyValue,numsci)



    else:
        # Default Case, sky subtraction is turned off.  No sky subtraction done to image.
        print "No sky subtraction requested, MDRIZSKY set to a value of 0..."
        for chip in numchips:
            imageObject[sciExt+','+str(chip)].header["MDRIZSKY"] = _skyValue
    
       
   
#this function can be called by users and will create an imageObject to send to
#the official function. I dunno, what's really a good name for this that the users
#can easily differentiate from the call we want? I chose "my" as the prefix cause
#it would be easy to add that to all the user independent calls and make it
#somewhat uniform

def mySubtractSky(configObj={},inputImageList=[], skyuser="", skysub=True, skywidth=0.1,
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
        
    #create image object    
    #create a configObject here with the paramDict settings?
    #call the real sky subtraction routine
    subtractSky(imageObject,paramDict)
    


###############################
##  Helper functions follow  ## 
###############################

def computeSky(dataArray, skypars, memmap=0):

    """ 
    Compute the sky value for the data array passed to the function        
    """

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
    print "Computed sky value for data array" : ",_computedsky

    return _computedsky



def _extractSkyValue(ImageStatsObject,skystat):
    if (skystat =="mode"):
        return ImageStatsObject.mode
    elif (skystat == "mean"):
        return ImageStatsObject.mean
    else:
        return ImageStatsObject.median



def _subtractSky(dataArray,skyValue,memmap=0):
    """
    subtract the given sky value from each the data array
    that has been passed
    """
    try:
        np.subtract(dataArray,skyValue,dataArray)

    except:
        raise IOError, "Unable to perform sky subtraction on data array"


#this is really related to each individual chip
def getreferencesky(imageObject,extension=1):

    _subtractedSky=imageObject[extension].header["MDRIZSKY"]
    _refplatescale=imageObject[extension].header["REFPLTSCL"]
    _platescale=imageObject[extension].header["PLATESCL"]
    
    return (_subtractedsky * (_refplatescale / _platescale)**2 )                















#this is the original function
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


