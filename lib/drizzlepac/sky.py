#!/usr/bin/env python
"""

Function for computing and subtracting the backgroud of
an image.  The algorithm employed here uses a sigma
clipped median of  each *sci* image in a data file.
Then the sky value for each detector is compared
and the lowest value is  subtracted from all chips
in the detector.  Finally, the MDRIZSKY keyword
is updated in the header of the input files.

:Authors:
    Christopher Hanley, Megan Sosey
"""


from __future__ import division  # confidence medium

import os
import sys

import util
from imageObject import imageObject
from stsci.tools import fileutil, teal, logutil
import processInput
import stsci.imagestats as imagestats
import numpy as np


__taskname__= "drizzlepac.sky" #looks in drizzlepac for sky.cfg
_step_num_ = 2  #this relates directly to the syntax in the cfg file


log = logutil.create_logger(__name__)


def help():
    print getHelpAsString()


def getHelpAsString():
    """
    return useful help from a file in the script directory called module.help
    """

    helpString = teal.getHelpFileAsString(__taskname__, __file__)

    return helpString


#this is the user access function
def sky(input=None,outExt=None,configObj=None, group=None, editpars=False, **inputDict):
    """
    Perform sky subtraction on input list of images

    Parameters
    ----------
    input : str or list of str
        a python list of image filenames, or just a single filename
    configObj : configObject
        an instance of configObject
    inputDict : dict, optional
        an optional list of parameters specified by the user
    outExt : str
        The extension of the output image. If the output already exists
        then the input image is overwritten

    Notes
    -----
    These are parameters that the configObj should contain by default,
    they can be altered on the fly using the inputDict

    Parameters that should be in configobj:

    ==========  ===================================================================
    Name        Definition
    ==========  ===================================================================
    skyuser		'KEYWORD in header which indicates a sky subtraction value to use'.
    skysub		'Perform sky subtraction?'
    skywidth	'Bin width for sampling sky statistics (in sigma)'
    skystat	 	'Sky correction statistics parameter'
    skylower	'Lower limit of usable data for sky (always in electrons)'
    skyupper	'Upper limit of usable data for sky (always in electrons)'
    skyclip		'Number of clipping iterations'
    skylsigma	'Lower side clipping factor (in sigma)'
    skyusigma	'Upper side clipping factor (in sigma)'
    ==========  ===================================================================

    The output from sky subtraction is a copy of the original input file
    where all the science data extensions have been sky subtracted.

    """


    if input is not None:
        inputDict['input']=input
        inputDict['output']=None
        inputDict['updatewcs']=False
        inputDict['group']=group
    else:
        print >> sys.stderr, "Please supply an input image"
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
            outsky = image.outputNames['outSky']
            if outExt not in outsky:
                outsky = outsky.replace("sky",outExt)
                image.outputNames['outSky']=outsky
                log.info(outsky)

    subtractSky(imageObjList,configObj,saveFile=saveFile)


#this is the workhorse looping function
def subtractSky(imageObjList,configObj,saveFile=False,procSteps=None):
    if procSteps is not None:
        procSteps.addStep('Subtract Sky')

    if not util.getConfigObjPar(configObj, 'skysub'):
        log.info('Sky Subtraction step not performed.')
        _addDefaultSkyKW(imageObjList)
        procSteps.endStep('Subtract Sky')
        return

    #General values to use
    step_name=util.getSectionName(configObj,_step_num_)
    paramDict = configObj[step_name]
    #get the sub-dictionary of values for this step alone and print them out
    log.info('USER INPUT PARAMETERS for Sky Subtraction Step:')
    util.printParams(paramDict, log=log)
    if 'skyfile' in paramDict and not util.is_blank(paramDict['skyfile']):
        _skyUserFromFile(imageObjList,paramDict['skyfile'])
    else:
        for image in imageObjList:
            log.info('Working on sky for: %s' % image._filename)
            _skySub(image,paramDict,saveFile=saveFile)

    if procSteps is not None:
        procSteps.endStep('Subtract Sky')

# this function applies user supplied sky values from an input file
def _skyUserFromFile(imageObjList,skyFile):
    """
    Apply sky value as read in from a user-supplied input file
    """
    skyKW="MDRIZSKY" #header keyword that contains the sky that's been subtracted

    # create dict of fname=sky pairs
    skyvals = {}
    skyapplied = False # flag whether sky has already been applied to images
    for line in open(skyFile):
        if line[0] == '#' and 'applied' in line:
            if '=' in line: linesep = '='
            if ':' in line: linesep = ':'
            appliedstr = line.split(linesep)[1].strip()
            if appliedstr.lower() in ['yes','true','y','t']:
                skyapplied = True
                print '...Sky values already applied by user...'

        if not util.is_blank(line) and line[0] != '#':
            lspl = line.split()
            svals = []
            for lvals in lspl[1:]:
                svals.append(float(lvals))
            skyvals[lspl[0]] = svals

    # Apply user values to appropriate input images
    for imageSet in imageObjList:
        fname = imageSet._filename
        numchips=imageSet._numchips
        sciExt=imageSet.scienceExt
        if fname in skyvals:
            print "    ...updating MDRIZSKY with user-supplied value."
            for chip in range(1,numchips+1,1):
                if len(skyvals[fname]) == 1:
                    _skyValue = skyvals[fname][0]
                else:
                    _skyValue = skyvals[fname][chip-1]

                chipext = '%s,%d'%(sciExt,chip)
                _updateKW(imageSet[chipext],fname,(sciExt,chip),skyKW,_skyValue)

                # Update internal record with subtracted sky value
                #
                # .computedSky:   value to be applied by the
                #                 adrizzle/ablot steps.
                # .subtractedSky: value already (or will be by adrizzle/ablot)
                #                 subtracted from the image
                if skyapplied:
                    imageSet[chipext].computedSky = 0.0 # used by adrizzle/ablot
                else:
                    imageSet[chipext].computedSky = _skyValue
                imageSet[chipext].subtractedSky = _skyValue
                print "Setting ",skyKW,"=",_skyValue
        else:
            print "*"*40
            print "*"
            print "WARNING:"
            print "    .... NO user-supplied sky value found for ",fname
            print "    .... Setting sky to a value of 0.0! "
            print "*"
            print "*"*40

#this is the main function that does all the real work in computing the
# statistical sky value for each image (set of chips)
def _skySub(imageSet,paramDict,saveFile=False):
    """
    subtract the sky from all the chips in the imagefile that imageSet represents

    imageSet is a single imageObject reference
    paramDict should be the subset from an actual config object
    if saveFile=True, then images that have been sky subtracted are saved to a predetermined output name
    else, overwrite the input images with the sky-subtracted results

    the output from sky subtraction is a copy of the original input file
    where all the science data extensions have been sky subtracted

    """


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
        print "User has computed their own sky values..."

        if skyuser != skyKW:
            print "    ...updating MDRIZSKY with supplied value."
            for chip in range(1,numchips+1,1):
                try:
                    chipext = '%s,%d'%(sciExt,chip)
                    _skyValue = imageSet[chipext].header[skyuser]

                except:
                    print "**************************************************************"
                    print "*"
                    print "*  Cannot find keyword ",skyuser," in ",imageSet._filename
                    print "*"
                    print "**************************************************************\n\n\n"
                    raise KeyError

                _updateKW(imageSet[sciExt+','+str(chip)],imageSet._filename,(sciExt,chip),skyKW,_skyValue)

                # Update internal record with subtracted sky value
                imageSet[chipext].subtractedSky = _skyValue
                imageSet[chipext].computedSky = 0.0
                print "Setting ",skyKW,"=",_skyValue

    else:
        # Compute our own sky values and record the values for use later.
        # The minimum sky value from all the  science chips in the exposure
        # is used as the reference sky for each chip

        log.info("Computing minimum sky ...")
        minSky=[] #store the sky for each chip
        minpscale = []

        for chip in range(1,numchips+1,1):
            myext=sciExt+","+str(chip)

            #add the data back into the chip, leave it there til the end of this function
            imageSet[myext].data=imageSet.getData(myext)

            image=imageSet[myext]
            _skyValue= _computeSky(image, paramDict, memmap=0)
            #scale the sky value by the area on sky
            # account for the case where no IDCSCALE has been set, due to a
            # lack of IDCTAB or to 'coeffs=False'.
            pscale=imageSet[myext].wcs.idcscale
            if pscale is None:
                log.warning("No Distortion coefficients available...using "
                            "default plate scale.")
                pscale = imageSet[myext].wcs.pscale
            _scaledSky=_skyValue / (pscale**2)
            #_skyValue=_scaledSky
            minSky.append(_scaledSky)
            minpscale.append(pscale)

        _skyValue = min(minSky)

        _reportedSky = _skyValue*(minpscale[minSky.index(_skyValue)]**2)
        log.info("Minimum sky value for all chips %s" % _reportedSky)

        #now subtract that value from all the chips in the exposure
        #and update the chips header keyword with the sub
        for chip in range(1,numchips+1,1):
            image=imageSet[sciExt,chip]
            myext = sciExt+","+str(chip)
            # account for the case where no IDCSCALE has been set, due to a
            # lack of IDCTAB or to 'coeffs=False'.
            idcscale = image.wcs.idcscale
            if idcscale is None: idcscale = image.wcs.pscale
            _scaledSky=_skyValue * (idcscale**2)
            image.subtractedSky = _scaledSky
            image.computedSky = _scaledSky
            log.info("Using sky from chip %d: %f\n" % (chip,_scaledSky))
            ###_subtractSky(image,(_scaledSky))
            # Update the header so that the keyword in the image is
            #the sky value which should be subtracted from the image
            _updateKW(image,imageSet._filename,(sciExt,chip),skyKW,_scaledSky)


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
    log.info("    Computed sky value/pixel for %s: %s "%
             (image.rootname, _skyValue))

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
    log.info('Updating keyword %s in %s' % (skyKW, filename + strexten))
    fobj = fileutil.openImage(filename, mode='update')
    fobj[exten].header.update(skyKW, Value,
                    comment='Sky value computed by AstroDrizzle')
    fobj.close()

def _addDefaultSkyKW(imageObjList):
    """Add MDRIZSKY keyword to SCI headers of all input images,
        if that keyword does not already exist.
    """
    skyKW = "MDRIZSKY"
    Value = 0.0
    for imageSet in imageObjList:
        fname = imageSet._filename
        numchips=imageSet._numchips
        sciExt=imageSet.scienceExt
        fobj = fileutil.openImage(fname, mode='update')
        for chip in range(1,numchips+1,1):
            exten = (sciExt,chip)
            if skyKW not in fobj[exten].header:
                fobj[exten].header.update(skyKW, Value,
                                comment='Sky value computed by AstroDrizzle')
                log.info("MDRIZSKY keyword not found in the %s[%s,%d] header."%(
                            fname,sciExt,chip))
                log.info("    Adding MDRIZSKY to header with default value of 0.")
        fobj.close()

#this is really related to each individual chip
#so pass in the image for that chip, image contains header and data
def getreferencesky(image,keyval):

    _subtractedSky=image.header[keyval]
    _refplatescale=image.header["REFPLTSCL"]
    _platescale=image.header["PLATESCL"]

    return (_subtractedsky * (_refplatescale / _platescale)**2 )
