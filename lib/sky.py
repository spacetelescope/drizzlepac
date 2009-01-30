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
import imageObject
from pytools import cfgpars
import imagestats
import numpy as N

#this is the main function that takes an imageObject and a parameter dictions
#made from the config obj. This is what the user function calls as well
def subtractSky(imageObject,paramDict={},saveFile=True):
    """
    subtract the sky from all the chips in the imagefile that imageObject represents
    
    imageObject contains all the information about the chips in the image file
    configObj is represented as a dict for now, but will prolly be an actual config object
    if saveFile=True, then images that have been sky subtracted are saved to a predetermined output name
    """
   
                   
    _skyValue=0.0    #this will be the sky value computed for the exposure                                                                  
    skyKW="MDRIZSKY" #header keyword that contains the sky that's been subtracted
    
    #just making sure, tricky users and all, these are things that will be used
    #by the sky function so we want them defined at least
    try:
        assert imageObject._numchips > 0, "invalid value for number of chips"
        assert imageObject._filename != '', "image object filename is empty!, doh!"
        assert imageObject.scienceExt !='', "image object science extension is empty!"
        assert imageObject._instrument !='', "image object instrument name is empty!"

    except AssertionError:
        raise AssertionError
        
    numchips=imageObject._numchips
    sciExt=imageObject.scienceExt
    
    # User Subtraction Case, User has done own sky subtraction,  
	# so use the image header value for subtractedsky value
    
    if paramDict["skyuser"] != '':
        print "User has done their own sky subtraction, updating MDRIZSKY with supplied value..."
       
        for chip in range(1,numchips+1,1):
            try:
                image=imageObject[sciExt+','+str(chip)]                
                _skyValue = image.header[paramDict["skyuser"]]

            except:
                print "**************************************************************"
                print "*"
                print "*  Cannot find keyword ",paramDict["skyuser"]," in ",imageObject._filename
                print "*"
                print "**************************************************************\n\n\n"
                raise KeyError

            image.header[paramDict[skyKW]] = _skyValue

    else:
        # Compute our own sky values and subtract them from the image copies.
        # for all instruments other than ACS the minimum sky value from all the
        # science chips in the exposure is used as the reference sky for each chip

        print "Computing minimum sky ..."
        minSky=[] #store the sky for each chip
                    
        #####FIX THIS#######            
        if ("STIS" in imageObject._instrument):
            for chip in range(1,numchips+1,1):
                # We need to account for the fact that STIS associations contain
                # separate exposures for the same chip within the same file.
                # In those cases, we do not want to use the minimum sky value
                # for the set of chips, but rather use each chip's own value.
                # NOTE: This can be generalized later with changes to PyDrizzle
                #       to provide an attribute that specifies whether each member
                #       associated with file is a separate exposure or not.
                #   WJH/CJH 
                image=imageObject["SCI,"+str(chip)]
                _skyValue=_computeSky(image,paramDict)
                _subtractSky(image,_skyValue)
                _updateKW(image,skyKW,_skyValue)
        
        else:
            for chip in range(1,numchips+1,1):
            	    myext=imageObject.scienceExt+","+str(chip)
                    image=imageObject[myext]
                    _skyValue= computeSky(image, paramDict, memmap=0)
                    minSky.append(_skyValue)
                    image.computedSky=_skyValue #update the keyword in the actual header here as well?

            _skyValue=min(minSky) #what if the sky is negative?

            #now subtract that value from all the chips in the exposure
            #and update the chips header keyword with the sub
            for chip in range(1,numchips+1,1):
                _subtractSky(image,_skyValue)
                _updateKW(image,skyKW,_skyValue)
            
        #update the value of MDRIZSKY in the global header
        _updateKW(imageObject[0],skyKW,_skyValue)
   
        if(saveFile):
            print "Saving output sky subtracted images....\n"
            for chip in range(1,numchips+1,1):
                myext="SCI,"+str(chip)
                image=imageObject[myext]
                print image.outputNames['outSky']
                image.writeto(image.outputNames['outSky'])
            
   
#this function can be called by users and will create an imageObject to send to
#the official function. I dunno, what's really a good name for this that the users
#can easily differentiate from the call we want? I chose "my" as the prefix cause
#it would be easy to add that to all the user independent calls and make it
#somewhat uniform

def mySubtractSky(configObj={},inputImageList=[], skyuser="", skysub=True, skywidth=0.1,skystat="median", 
    skylower="INDEF", skyupper="INDEF",skyclip=5, skysligma=4.,skyusigma=4.,saveFile=True):

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
        print "Empty input image list given to Sky routine, checking configObj"
        if(len(configObj["imageList"]) == 0):
            print "no image list in configObject either"
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
    subtractSky(imageObject,paramDict,saveFile)
    


###############################
##  Helper functions follow  ## 
###############################

def _computeSky(image, skypars, memmap=0):

    """ 
    Compute the sky value for the data array passed to the function 
    image is a pyfits object which contains the data and the header
    for one image extension
    
    skypars is passed in as paramDict
    """

	#this object contains the returned values from the image stats routine
    _tmp = imagestats.ImageStats(image.data,
            fields      = skypars['skystat'],
            lower       = skypars['skylower'],
            upper       = skypars['skyupper'],
            nclip       = skypars['skyclip'],
            lsig        = skypars['skylsigma'],
            usig        = skypars['skyusigma'],
            binwidth    = skypars['skywidth']
            )

    _skyValue = _extractSkyValue(_tmp,skypars['skystat'].lower())
    print "Computed sky value for %s: "%image.rootname, _skyValue

    return _skyValue



def _extractSkyValue(imstatObject,skystat):
    if (skystat =="mode"):
        return imstatObject.mode 
    elif (skystat == "mean"):
        return imstatObject.mean
    else:
        return imstatObject.median 



def _subtractSky(image,skyValue,memmap=0):
    """
    subtract the given sky value from each the data array
    that has been passed. image is a pyfits object that 
    contains the data and header for one image extension
    """
    try:
        N.subtract(image.data,skyValue,image.data)

    except IOError:
        print "Unable to perform sky subtraction on data array"
        raise IOError 


def _updateKW(image, skyKW, _skyValue):
    """update the header with the kw,value"""

    image.header.update(skyKW,_skyValue)
    
    
#this is really related to each individual chip
#so pass in the image for that chip, image contains header and data
def getreferencesky(image,keyval):

    _subtractedSky=image.header[keyval]
    _refplatescale=image.header["REFPLTSCL"]
    _platescale=image.header["PLATESCL"]
    
    return (_subtractedsky * (_refplatescale / _platescale)**2 )                
















