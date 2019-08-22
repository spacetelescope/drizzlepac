#!/usr/bin/env python
"""
This step measures, subtracts and/or equalizes the sky from each
input image while recording the subtracted value in the image header.

:Authors:
    Christopher Hanley, Megan Sosey, Mihai Cara

:License: :doc:`LICENSE`

"""
import os, sys

from .imageObject import imageObject
from stsci.tools import fileutil, teal, logutil
try:
    from stsci.tools.bitmask import interpret_bit_flags
except ImportError:
    from stsci.tools.bitmask import interpret_bits_value as interpret_bit_flags


from stsci.skypac.skymatch import skymatch
from stsci.skypac.utils import MultiFileLog, ResourceRefCount, ext2str, \
     file_name_components, in_memory_mask, temp_mask_file, openImageEx
from stsci.skypac.parseat import FileExtMaskInfo, parse_at_file

from . import processInput
import stsci.imagestats as imagestats
import numpy as np

from . import util
from .version import *


__taskname__= "drizzlepac.sky" #looks in drizzlepac for sky.cfg
_step_num_ = 2  #this relates directly to the syntax in the cfg file


log = logutil.create_logger(__name__, level=logutil.logging.NOTSET)


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
    skymethod   'Sky computation method'
    skysub		'Perform sky subtraction?'
    skywidth	'Bin width of histogram for sampling sky statistics (in sigma)'
    skystat	 	'Sky correction statistics parameter'
    skylower	'Lower limit of usable data for sky (always in electrons)'
    skyupper	'Upper limit of usable data for sky (always in electrons)'
    skyclip		'Number of clipping iterations'
    skylsigma	'Lower side clipping factor (in sigma)'
    skyusigma	'Upper side clipping factor (in sigma)'
    skymask_cat 'Catalog file listing image masks'
    use_static  'Use static mask for skymatch computations?'
    sky_bits    'Integer mask bit values considered good pixels in DQ array'
    skyfile     'Name of file with user-computed sky values'
    skyuser     'KEYWORD indicating a sky subtraction value if done by user'
    in_memory   'Optimize for speed or for memory use'

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
        print("Please supply an input image", file=sys.stderr)
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
    # if neither 'skyfile' nor 'skyuser' are specified, subtractSky will
    # call _skymatch to perform "sky background matching". When 'skyuser'
    # is specified, subtractSky will call the old _skysub.

    if procSteps is not None:
        procSteps.addStep('Subtract Sky')

    #General values to use
    step_name=util.getSectionName(configObj,_step_num_)
    paramDict = configObj[step_name]

    if not util.getConfigObjPar(configObj, 'skysub'):
        log.info('Sky Subtraction step not performed.')
        _addDefaultSkyKW(imageObjList)
        if 'skyuser' in paramDict and not util.is_blank(paramDict['skyuser']):
            kwd = paramDict['skyuser'].lstrip()
            if kwd[0] == '@':
                # user's sky values are in a file:
                log.info("Retrieving user computed sky values from file '{}'"
                         .format(kwd[1:]))
                _skyUserFromFile(imageObjList, kwd[1:],apply_sky=False)
            else:
                # user's sky values are stored in a header keyword:
                log.info("Retrieving user computed sky values from image "
                         "headers ")
                log.info("recorded in the '{:s}' header keywords."
                         .format(paramDict['skyuser']))
                for image in imageObjList:
                    log.info('Working on sky for: %s' % image._filename)
                    _skyUserFromHeaderKwd(image, paramDict)
        else:
            # reset "computedSky" chip's attribute:
            for image in imageObjList:
                numchips    = image._numchips
                extname     = image.scienceExt
                for extver in range(1, numchips + 1, 1):
                    chip = image[extname, extver]
                    if not chip.group_member:
                        continue
                    chip.computedSky = None

        if procSteps is not None:
            procSteps.endStep('Subtract Sky')
        return

    #get the sub-dictionary of values for this step alone and print them out
    log.info('USER INPUT PARAMETERS for Sky Subtraction Step:')
    util.printParams(paramDict, log=log)
    if 'skyfile' in paramDict and not util.is_blank(paramDict['skyfile']):
        _skyUserFromFile(imageObjList,paramDict['skyfile'])
    else:
        # in_memory:
        if 'in_memory' in configObj:
            inmemory = configObj['in_memory']
        elif len(imageObjList) > 0 and imageObjList[0].inmemory is not None:
            inmemory = imageObjList[0].inmemory
        else:
            inmemory = False
        # clean:
        if 'STATE OF INPUT FILES' in configObj and \
           'clean' in configObj['STATE OF INPUT FILES']:
            clean = configObj['STATE OF INPUT FILES']['clean']
        else:
            clean = True

        _skymatch(imageObjList, paramDict, inmemory, clean, log)

    if procSteps is not None:
        procSteps.endStep('Subtract Sky')


def _skymatch(imageList, paramDict, in_memory, clean, logfile):
    # '_skymatch' converts input imageList and other parameters to
    # data structures accepted by the "skymatch" package.
    # It also creates a temporary mask by combining 'static' mask,
    # DQ image, and user-supplied mask. The combined mask is then
    # passed to 'skymatch' to be used for excluding "bad" pixels.

    #header keyword that contains the sky that's been subtracted
    skyKW = "MDRIZSKY"

    nimg = len(imageList)
    if nimg == 0:
        ml.logentry("Skymatch needs at least one images to perform{0}" \
                    "sky matching. Nothing to be done.",os.linesep)
        return

    # create a list of input file names as provided by the user:
    user_fnames   = []
    loaded_fnames = []
    filemaskinfos = nimg * [ None ]
    for img in imageList:
        user_fnames.append(img._original_file_name)
        loaded_fnames.append(img._filename)

    # parse sky mask catalog file (if any):
    catfile = paramDict['skymask_cat']
    if catfile:
        #extname = imageList[0].scienceExt
        #assert(extname is not None and extname != '')
        catfile = catfile.strip()
        mfindx = parse_at_file(fname = catfile,
                               default_ext = ('SCI','*'),
                               default_mask_ext = 0,
                               clobber      = False,
                               fnamesOnly   = True,
                               doNotOpenDQ  = True,
                               match2Images = user_fnames,
                               im_fmode     = 'update',
                               dq_fmode     = 'readonly',
                               msk_fmode    = 'readonly',
                               logfile      = MultiFileLog(console=True),
                               verbose      = True)
        for p in mfindx:
            filemaskinfos[p[1]] = p[0]

    # step through the list of input images and create
    # combined (static + DQ + user supplied, if any) mask, and
    # create a list of FileExtMaskInfo objects to be passed
    # to 'skymatch' function.
    #
    # This needs to be done in several steps, mostly due to the fact that
    # the mask catalogs use "original" (e.g., GEIS, WAIVER FITS) file names
    # while ultimately we want to open the version converted to MEF. Second
    # reason is that we want to combine user supplied masks with DQ+static
    # masks provided by astrodrizzle.
    new_fi = []
    sky_bits = interpret_bit_flags(paramDict['sky_bits'])
    for i in range(nimg):
        # extract extension information:
        extname = imageList[i].scienceExt
        extver  = imageList[i].group
        if extver is None:
            extver = imageList[i].getExtensions()
        assert(extname is not None and extname != '')
        assert(extver)

        # create a new FileExtMaskInfo object
        fi = FileExtMaskInfo(default_ext=(extname,'*'),
                             default_mask_ext=0,
                             clobber=False,
                             doNotOpenDQ=True,
                             fnamesOnly=False,
                             im_fmode='update',
                             dq_fmode='readonly',
                             msk_fmode='readonly')

        # set image file and extensions:
        fi.image = loaded_fnames[i]
        extlist  = [ (extname,ev) for ev in extver ]
        fi.append_ext(extlist)

        # set user masks if any (this will open the files for a later use):
        fi0 = filemaskinfos[i]
        if fi0 is not None:
            nmask = len(fi0.mask_images)
            for m in range(nmask):
                mask = fi0.mask_images[m]
                ext  = fi0.maskext[m]
                fi.append_mask(mask, ext)
        fi.finalize()

        # combine user masks with static masks:
        assert(len(extlist) == fi.count) #TODO: <-- remove after thorough testing

        masklist = []
        mextlist = []

        for k in range(fi.count):
            if fi.mask_images[k].closed:
                umask = None
            else:
                umask = fi.mask_images[k].hdu[fi.maskext[k]].data
            (mask, mext) = _buildStaticDQUserMask(imageList[i], extlist[k],
                               sky_bits, paramDict['use_static'],
                               fi.mask_images[k], fi.maskext[k], in_memory)

            masklist.append(mask)
            mextlist.append(mext)

        # replace the original user-supplied masks with the
        # newly computed combined static+DQ+user masks:
        fi.clear_masks()
        for k in range(fi.count):
            if in_memory and mask is not None:
                # os.stat() on the "original_fname" of the mask will fail
                # since this is a "virtual" mask. Therefore we need to compute
                # mask_stat ourselves. We will simply use id(data) for this:
                mstat = os.stat_result((0,id(mask.hdu)) + 8*(0,))
                fi.append_mask(masklist[k], mextlist[k], mask_stat=mstat)
            else:
                fi.append_mask(masklist[k], mextlist[k])
            if masklist[k]:
                masklist[k].release()
        fi.finalize()

        new_fi.append(fi)

    # Run skymatch algorithm:
    skymatch(new_fi,
             skymethod   = paramDict['skymethod'],
             skystat     = paramDict['skystat'],
             lower       = paramDict['skylower'],
             upper       = paramDict['skyupper'],
             nclip       = paramDict['skyclip'],
             lsigma      = paramDict['skylsigma'],
             usigma      = paramDict['skyusigma'],
             binwidth    = paramDict['skywidth'],
             skyuser_kwd = skyKW,
             units_kwd   = 'BUNIT',
             readonly    = not paramDict['skysub'],
             dq_bits     = None,
             optimize    = 'inmemory' if in_memory else 'balanced',
             clobber     = True,
             clean       = clean,
             verbose     = True,
             flog        = MultiFileLog(console = True, enableBold = False),
             _taskname4history = 'AstroDrizzle')

    # Populate 'subtractedSky' and 'computedSky' of input image objects:
    for i in range(nimg):
        assert(not new_fi[i].fnamesOnly and not new_fi[i].image.closed)
        image = imageList[i]
        skysubimage = new_fi[i].image.hdu
        numchips = image._numchips
        extname = image.scienceExt
        assert(os.path.samefile(image._filename, skysubimage.filename()))

        for extver in range(1, numchips + 1, 1):
            chip = image[extname, extver]
            if not chip.group_member:
                continue
            subtracted_sky = skysubimage[extname, extver].header.get(skyKW, 0.)
            chip.subtractedSky = subtracted_sky
            chip.computedSky = subtracted_sky

    # clean-up:
    for fi in new_fi:
        fi.release_all_images()

def _buildStaticDQUserMask(img, ext, sky_bits, use_static, umask,
                           umaskext, in_memory):
    # creates a temporary mask by combining 'static' mask,
    # DQ image, and user-supplied mask.

    def merge_masks(m1, m2):
        if m1 is None: return m2
        if m2 is None: return m1
        return np.logical_and(m1, m2).astype(np.uint8)

    mask = None

    # build DQ mask
    if sky_bits is not None:
        mask = img.buildMask(img[ext]._chip,bits=sky_bits)

    # get correct static mask mask filenames/objects
    staticMaskName = img[ext].outputNames['staticMask']
    smask = None
    if use_static:
        if img.inmemory:
            if staticMaskName in img.virtualOutputs:
                smask = img.virtualOutputs[staticMaskName].data
        else:
            if staticMaskName is not None and os.path.isfile(staticMaskName):
                sm, dq = openImageEx(
                    staticMaskName,
                    mode='readonly',
                    memmap=False,
                    saveAsMEF=False,
                    clobber=False,
                    imageOnly=True,
                    openImageHDU=True,
                    openDQHDU=False,
                    preferMEF=False,
                    verbose=False
                )
                if sm.hdu is not None:
                    smask = sm.hdu[0].data
                    sm.release()
            else:
                log.warning("Static mask for file \'{}\', ext={} NOT FOUND." \
                            .format(img._filename, ext))
        # combine DQ and static masks:
        mask = merge_masks(mask, smask)

    # combine user mask with the previously computed mask:
    if umask is not None and not umask.closed:
        if mask is None:
            # return user-supplied mask:
            umask.hold()
            return (umask, umaskext)
        else:
            # combine user mask with the previously computed mask:
            dtm  = umask.hdu[umaskext].data
            mask = merge_masks(mask, dtm)

    if mask is None:
        return (None, None)
    elif mask.sum() == 0:
        log.warning("All pixels masked out when applying DQ, " \
                    "static, and user masks!")

    # save mask to a temporary file:
    (root,suffix,fext) = file_name_components(img._filename)
    if in_memory:
        tmpmask = in_memory_mask(mask)
        strext = ext2str(ext, compact=True, default_extver=None)
        tmpmask.original_fname = "{1:s}{0:s}{2:s}{0:s}{3:s}" \
            .format('_', root, suffix, 'in-memory_skymatch_mask')
    else:
        (tmpfname, tmpmask) = temp_mask_file(mask, root,
            prefix='', suffix='skymatch_mask', ext=ext,
            randomize_prefix=False)
        img[ext].outputNames['skyMatchMask'] = tmpfname

    return (tmpmask, 0)

# this function applies user supplied sky values from an input file
def _skyUserFromFile(imageObjList, skyFile, apply_sky=None):
    """
    Apply sky value as read in from a user-supplied input file.

    """
    skyKW="MDRIZSKY" #header keyword that contains the sky that's been subtracted

    # create dict of fname=sky pairs
    skyvals = {}
    if apply_sky is None:
        skyapplied = False # flag whether sky has already been applied to images
    else:
        skyapplied = apply_sky

    for line in open(skyFile):
        if apply_sky is None and line[0] == '#' and 'applied' in line:
            if '=' in line: linesep = '='
            if ':' in line: linesep = ':'
            appliedstr = line.split(linesep)[1].strip()
            if appliedstr.lower() in ['yes','true','y','t']:
                skyapplied = True
                print('...Sky values already applied by user...')

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
            print("    ...updating MDRIZSKY with user-supplied value.")
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
                    imageSet[chipext].computedSky = None # used by adrizzle/ablot
                else:
                    imageSet[chipext].computedSky = _skyValue
                imageSet[chipext].subtractedSky = _skyValue
                print("Setting ",skyKW,"=",_skyValue)
        else:
            print("*"*40)
            print("*")
            print("WARNING:")
            print("    .... NO user-supplied sky value found for ",fname)
            print("    .... Setting sky to a value of 0.0! ")
            print("*")
            print("*"*40)

def _skyUserFromHeaderKwd(imageSet,paramDict):
    """
    subtract the sky from all the chips in the imagefile that imageSet represents

    imageSet is a single imageObject reference
    paramDict should be the subset from an actual config object

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
        print("User has computed their own sky values...")

        if skyuser != skyKW:
            print("    ...updating MDRIZSKY with supplied value.")
            for chip in range(1,numchips+1,1):
                chipext = '%s,%d'%(sciExt,chip)
                if not imageSet[chipext].group_member:
                    # skip extensions/chips that will not be processed
                    continue
                try:
                    _skyValue = imageSet[chipext].header[skyuser]
                except:
                    print("**************************************************************")
                    print("*")
                    print("*  Cannot find keyword ",skyuser," in ",imageSet._filename)
                    print("*")
                    print("**************************************************************\n\n\n")
                    raise KeyError

                _updateKW(imageSet[sciExt+','+str(chip)],
                          imageSet._filename,(sciExt,chip),skyKW,_skyValue)

                # Update internal record with subtracted sky value
                imageSet[chipext].subtractedSky = _skyValue
                imageSet[chipext].computedSky = None
                print("Setting ",skyKW,"=",_skyValue)

#this is the main function that does all the real work in computing the
# statistical sky value for each image (set of chips)
# mcara: '_skySub' is obsolete now:
#        was replaced with '_skyUserFromHeaderKwd' and '_skymatch'
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
        print("User has computed their own sky values...")

        if skyuser != skyKW:
            print("    ...updating MDRIZSKY with supplied value.")
            for chip in range(1,numchips+1,1):
                try:
                    chipext = '%s,%d'%(sciExt,chip)
                    _skyValue = imageSet[chipext].header[skyuser]

                except:
                    print("**************************************************************")
                    print("*")
                    print("*  Cannot find keyword ",skyuser," in ",imageSet._filename)
                    print("*")
                    print("**************************************************************\n\n\n")
                    raise KeyError

                _updateKW(imageSet[sciExt+','+str(chip)],imageSet._filename,(sciExt,chip),skyKW,_skyValue)

                # Update internal record with subtracted sky value
                imageSet[chipext].subtractedSky = _skyValue
                imageSet[chipext].computedSky = None
                print("Setting ",skyKW,"=",_skyValue)

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
            _skyValue= _computeSky(image, paramDict, memmap=False)
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

def _computeSky(image, skypars, memmap=False):

    """
    Compute the sky value for the data array passed to the function
    image is a fits object which contains the data and the header
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



def _subtractSky(image,skyValue,memmap=False):
    """
    subtract the given sky value from each the data array
    that has been passed. image is a fits object that
    contains the data and header for one image extension
    """
    try:
        np.subtract(image.data,skyValue,image.data)

    except IOError:
        print("Unable to perform sky subtraction on data array")
        raise IOError


def _updateKW(image, filename, exten, skyKW, Value):
    """update the header with the kw,value"""
    # Update the value in memory
    image.header[skyKW] = Value

    # Now update the value on disk
    if isinstance(exten,tuple):
        strexten = '[%s,%s]'%(exten[0],str(exten[1]))
    else:
        strexten = '[%s]'%(exten)
    log.info('Updating keyword %s in %s' % (skyKW, filename + strexten))
    fobj = fileutil.openImage(filename, mode='update', memmap=False)
    fobj[exten].header[skyKW] = (Value, 'Sky value computed by AstroDrizzle')
    fobj.close()

def _addDefaultSkyKW(imageObjList):
    """Add MDRIZSKY keyword to "commanded" SCI headers of all input images,
        if that keyword does not already exist.
    """
    skyKW = "MDRIZSKY"
    Value = 0.0
    for imageSet in imageObjList:
        fname = imageSet._filename
        numchips=imageSet._numchips
        sciExt=imageSet.scienceExt
        fobj = fileutil.openImage(fname, mode='update', memmap=False)
        for chip in range(1,numchips+1,1):
            ext = (sciExt,chip)
            if not imageSet[ext].group_member:
                # skip over extensions not used in processing
                continue
            if skyKW not in fobj[ext].header:
                fobj[ext].header[skyKW] = (Value, 'Sky value computed by AstroDrizzle')
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


def help(file=None):
    """
    Print out syntax help for running astrodrizzle

    Parameters
    ----------
    file : str (Default = None)
        If given, write out help to the filename specified by this parameter
        Any previously existing file with this name will be deleted before
        writing out the help.

    """
    helpstr = getHelpAsString(docstring=True, show_ver = True)
    if file is None:
        print(helpstr)
    else:
        if os.path.exists(file): os.remove(file)
        f = open(file, mode = 'w')
        f.write(helpstr)
        f.close()


def getHelpAsString(docstring = False, show_ver = True):
    """
    return useful help from a file in the script directory called
    __taskname__.help

    """
    install_dir = os.path.dirname(__file__)
    taskname = util.base_taskname(__taskname__, __package__)
    htmlfile = os.path.join(install_dir, 'htmlhelp', taskname + '.html')
    helpfile = os.path.join(install_dir, taskname + '.help')

    if docstring or (not docstring and not os.path.exists(htmlfile)):
        if show_ver:
            helpString = os.linesep + \
                ' '.join([__taskname__, 'Version', __version__,
                ' updated on ', __version_date__]) + 2*os.linesep
        else:
            helpString = ''
        if os.path.exists(helpfile):
            helpString += teal.getHelpFileAsString(taskname, __file__)
        else:
            if __doc__ is not None:
                helpString += __doc__ + os.linesep
    else:
        helpString = 'file://' + htmlfile

    return helpString


sky.__doc__ = getHelpAsString(docstring = True, show_ver = False)
