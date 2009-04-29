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
from pytools import fileutil
import processInput
import imagestats
import os
import numpy as np

__taskname__= "betadrizzle.sky" #looks in betadrizzle for sky.cfg
_step_num_ = 2  #this relates directly to the syntax in the cfg file


def help():
    print getHelpAsString()
        
def getHelpAsString():
    """ 
    return useful help from a file in the script directory called module.help
    """
    #get the local library directory where the code is stored
    localDir=os.path.split(__file__)
    helpfile=__taskname__.split(".")
    helpfile=localDir[0]+"/"+helpfile[1]+".help"
    
    if os.access(helpfile,os.R_OK):
        fh=open(helpfile,'r')
        ss=fh.readlines()
        fh.close()
        helpString=""
        for line in ss:
            helpString+=line
    else:    
        helpString=__doc__

    return helpString

#this is the user access function
def sky(input=None,outExt=None,configObj=None, group=None, editpars=False, **inputDict):
    """
    input is a python list of image filenames, or just a single filename
    configObj is an instance of configObject
    inputDict is an optional list of parameters specified by the user
    outExt is the extension of the output image. If the output already exists
      then the input image is overwritten
    
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
    
    
    if input is not None:
        inputDict['input']=input  
        inputDict['output']=None
        inputDict['updatewcs']=False
        inputDict['workinplace']=True   
        inputDict['group']=group   
    else:
        print "Please supply an input image"
        raise ValueError

    configObj = util.getDefaultConfigObj(__taskname__,configObj,inputDict,loadOnly=(not editpars))
    if configObj is None:
        return
    
    if not editpars:
        run(configObj,outExt=outExt)
         

#this is the function that will be called from TEAL
def run(configObj,outExt=None):
 
    #now we really just need the imageObject list created for the dataset
    filelist,output,ivmlist,oldasndict=processInput.processFilenames(configObj['input'],None)

    imageObjList=processInput.createImageObjectList(filelist,instrpars={},group=configObj['group'])
    
    #set up the output names, if no extension given the default will be used
    #otherwise, the user extension is used and if the file already exists it's overwritten
    saveFile = False
    if(outExt not in [None,'','None']):    
        saveFile = True
        for image in imageObjList:
            outsky=image.outputNames['outSky']
            if outExt not in outsky:
                outsky=outsky.replace("sky",outExt)
                image.outputNames['outSky']=outsky
                print outsky

    subtractSky(imageObjList,configObj,saveFile=saveFile)


#this is the workhorse looping function
def subtractSky(imageObjList,configObj,saveFile=False):
    if not util.getConfigObjPar(configObj,'skysub'):
        print 'Sky Subtraction step not performed.'
        return

    for image in imageObjList:
        print "Working on sky for: ",image._filename
        _skySub(configObj,image,saveFile=saveFile)
    

#this is the main function that does all the real work
def _skySub(configObj=None,imageSet=None,saveFile=False):
    """
    subtract the sky from all the chips in the imagefile that imageSet represents
    imageSet is a single imageObject reference
    configObj should be an actual config object by now
    if saveFile=True, then images that have been sky subtracted are saved to a predetermined output name
    else, overwrite the input images with the sky-subtracted results

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
                
            _updateKW(imageSet[sciExt+','+str(chip)],imageSet._filename,(sciExt,chip),skyKW,_skyValue)
                        
        #update the value of MDRIZSKY in the global header
        _updateKW(imageSet["PRIMARY"],imageSet._filename,"PRIMARY",skyKW,_skyValue)
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
            pscale=imageSet[myext].wcs.idcscale
            print 'pixel area on sky=',pscale
            _scaledSky=_skyValue / (pscale**2)
            print 'scaledSky=',_scaledSky
            #_skyValue=_scaledSky
            minSky.append(_scaledSky)
            
            #update the keyword in the actual header here as well
            image.computedSky=_scaledSky #this is the scaled sky value

        _skyValue = min(minSky)
        print "Minimum sky value for all chips ",_skyValue

        #now subtract that value from all the chips in the exposure
        #and update the chips header keyword with the sub
        for chip in range(1,numchips+1,1):
            image=imageSet._image[sciExt,chip]
            myext = sciExt+","+str(chip)
            _scaledSky=_skyValue * (image.wcs.idcscale**2)
            image.subtractedSky = _scaledSky
            print "\nsubtracting scaled sky from chip %d: %f\n"%(chip,_scaledSky)
            _subtractSky(image,(_scaledSky))
            _updateKW(image,imageSet._filename,(sciExt,chip),skyKW,_scaledSky) #I updated this so that the keyword in the image is 
                                            #the sky value actually subtracted from the image

            if not saveFile:
                print "**Updating input image with sky subtracted image"
                # Write out the sky-subtracted array back to the input image
                imageSet.updateData(sciExt+","+str(chip),image.data)
            
        #update the value of MDRIZSKY in the global header
        # This does not make sense for STIS ASN files that
        #haven't been chunked up into separate fits files already
        _updateKW(imageSet["PRIMARY"],imageSet._filename,'PRIMARY',skyKW,_skyValue)
   
    if(saveFile):
        print "Saving output sky subtracted image: ",imageSet.outputNames["outSky"]
        #get the rest of the data extensions
        imageSet.getAllData(exclude="SCI")
        if os.access(imageSet.outputNames['outSky'],os.F_OK):
            os.remove(imageSet.outputNames['outSky'])
            
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
    print "Computed unscaled sky value for %s: "%image.rootname, _skyValue
    
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


def _updateKW(image, filename, exten, skyKW, Value):
    """update the header with the kw,value"""
    # Update the value in memory
    image.header.update(skyKW,Value)

    # Now update the value on disk
    if isinstance(exten,tuple):
        strexten = '[%s,%s]'%(exten[0],str(exten[1]))
    else:
        strexten = '[%s]'%(exten)
    print 'Updating keyword ',skyKW,' in ',filename+strexten
    fobj = fileutil.openImage(filename,mode='update')
    fobj[exten].header.update(skyKW,Value)
    fobj.close()
    
#this is really related to each individual chip
#so pass in the image for that chip, image contains header and data
def getreferencesky(image,keyval):

    _subtractedSky=image.header[keyval]
    _refplatescale=image.header["REFPLTSCL"]
    _platescale=image.header["PLATESCL"]
    
    return (_subtractedsky * (_refplatescale / _platescale)**2 )                













