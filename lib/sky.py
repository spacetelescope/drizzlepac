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
"""
import util
from imageObject import imageObject
from pytools import cfgpars
import imagestats
import numpy as np

__taskname__= "BigBlackBox.sky" #looks in BigBlackBox for sky.cfg
_step_num_ = 2  #this relates directly to the syntax in the cfg file

def getHelpAsString():
    """ I'm thinking we could just make a file called sky.help
    then use this function to read it into an array or list and return that?
    """

    helpString="  Since the minimum sky is calculated between all chips,"
    helpString+= "it is possible that the chips have a different platescale"
    helpString+= "so the minimum value needs to be compared on the sky." 
    helpString+= "Each sky minimum is ratioed with the platescale and THAT"
    helpString+= "value is stored in the MDRIZSKY keyword in the header. It"
    helpString+= "is also assumed that when the user has subtracted the sky themselves,"
    helpString+= "that has already taken this into account, so no extra scaling is"
    helpString+= "done later on in the code to account for it. "

    return helpString

#this is the user access function
def sky(imageList=None,configObj=None, editpars=False, **inputDict):
    """
    imageList is a python list of image filenames, or just a single filename
    configObj is an instance of configObject
    inputDict is an optional list of parameters specified by the user
    saveFile decided whether to save intermediate files, set to False to delete them
    
    These are parameters that the configObj should contain by default,
    they can be altered on the fly using the inputDict


    params that should be in configobj
    ---------------------------------------
    skyuser		'KEYWORD in header which indicates a sky subtraction value to use'.
    skysub		'Perform sky subtraction?'
    skywidth	'Bin width for sampling sky statistics (in sigma)'
    skystat	 	'Sky correction statistics parameter'
    skylower	'Lower limit of usable data for sky (always in electrons)'
    skyupper	'Upper limit of usable data for sky (always in electrons)'
    skyclip		'Number of clipping iterations'
    skylsigma	'Lower side clipping factor (in sigma)'
    skyusigma	'Upper side clipping factor (in sigma)'


    the output from sky subtraction is a copy of the original input file
    where all the science data extensions have been sky subtracted
    
    """
    inputDict["input"]=imageList        
    configObj = util.getDefaultConfigObj(__taskname__,configObj,inputDict,loadOnly=(not editpars))
    if configObj is None:
        return

    run(configObj)
     

#this is the function that will be called from TEAL
def run(configObj):
 
    imgObjList,outwcs = processInput.setCommonInput(configObj,createOutwcs=False) #outwcs is not neaded here

    subtractSky(imgObjList,configObj)


#this is the workhorse function
def subtractSky(imageObjList,configObj):
    print "image list:",imageObjList
    print
    print "configobj:",configObj
    print
    for image in imageObjList:
        print "Working on sky for: ",image
        _subtractSky(configObj,imageSet=image,saveFile=configObj["clean"])
    

#this is the main function that does all the real work
def _subtractSky(configObj=None,imageSet=None,saveFile=True):
    """
    subtract the sky from all the chips in the imagefile that imageSet represents
    imageSet is a single imageObject reference
    configObj should be an actual config object by now
    if saveFile=True, then images that have been sky subtracted are saved to a predetermined output name

    the output from sky subtraction is a copy of the original input file
    where all the science data extensions have been sky subtracted

    """
   
    #General values to use    
    step_name=util.getSectionName(configObj,_step_num_)  
    
    #get the sub-dictionary of values for this step alone
    paramDict=configObj[step_name]         

    print "\nUSER INPUT PARAMETERS for SKY SUBTRACTION:"
    util.printParams(paramDict)        
             
    _skyValue=0.0    #this will be the sky value computed for the exposure                                                                  
    skyKW="MDRIZSKY" #header keyword that contains the sky that's been subtracted
    
    #just making sure, tricky users and all, these are things that will be used
    #by the sky function so we want them defined at least
    try:
        assert imageSet._numchips > 0, "invalid value for number of chips"
        assert imageSet._filename != '', "image object filename is empty!, doh!"
        assert imageSet._rootname != '', "image rootname is empty!, doh!"
        assert imageSet.scienceExt !='', "image object science extension is empty!"
        
    except AssertionError:
        raise AssertionError
        
    numchips=imageSet._numchips
    sciExt=imageSet.scienceExt
     
    # User Subtraction Case, User has done own sky subtraction,  
    # so use the image header value for subtractedsky value    
    skyuser=paramDict["skyuser"]
    
    if skyuser != '':
        print "User has done their own sky subtraction, updating MDRIZSKY with supplied value..."
       
        for chip in range(1,numchips+1,1):
            try:
                _skyValue = imageSet._image["PRIMARY"].header[skyuser]

            except:
                print "**************************************************************"
                print "*"
                print "*  Cannot find keyword ",skyuser," in ",imageSet._filename
                print "*"
                print "**************************************************************\n\n\n"
                raise KeyError
                
            _updateKW(imageSet[sciExt+','+str(chip)],skyKW,_skyValue)
                        
        #update the value of MDRIZSKY in the global header
        _updateKW(imageSet["PRIMARY"],skyKW,_skyValue)
        print skyKW,"=",_skyValue

    else:
        # Compute our own sky values and subtract them from the image copies.
        # The minimum sky value from all the  science chips in the exposure
        # is used as the reference sky for each chip

        print "Computing minimum sky ..."
        minSky=[] #store the sky for each chip
        
        for chip in range(1,numchips+1,1):
            myext=sciExt+","+str(chip)
            
            #add the data back into the chip, leave it there til the end of this function
            imageSet[myext].data=imageSet.getData(myext)
            
            image=imageSet[myext]
            _skyValue= _computeSky(image, paramDict, memmap=0)
            #scale the sky value by the area on sky
            pscale=imageSet[myext].wcs.pscale
            scaledSky=_skyValue / (pscale**2)
            _skyValue=scaledSky
            minSky.append(_skyValue)
            
            #update the keyword in the actual header here as well
            image.computedSky=_skyValue #this is the scaled sky value

        _skyValue = min(minSky)
        print "Minimum sky value for all chips ",_skyValue

        #now subtract that value from all the chips in the exposure
        #and update the chips header keyword with the sub
        for chip in range(1,numchips+1,1):
            image=imageSet._image[sciExt,chip]
            _subtractSky(image,(_skyValue * image.wcs.pscale))
            _updateKW(image,skyKW,_skyValue)
            
        #update the value of MDRIZSKY in the global header
        # This does not make sense for STIS ASN files that
        #haven't been chunked up into separate fits files already
        _updateKW(imageSet[0],skyKW,_skyValue)
   
    if(saveFile):
        print "Saving output sky subtracted image: ",imageSet.outputNames["outSky"]
        #get the rest of the data extensions
        imageSet.getAllData(exclude="SCI")
        try:
            imageSet._image.writeto(imageSet.outputNames['outSky'])
        except IOError:
            print "Image already exists on disk!"
            return IOError
                    
    imageSet.close() #remove the data from memory



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
    
    del _tmp
    
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
        np.subtract(image.data,skyValue,image.data)

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













