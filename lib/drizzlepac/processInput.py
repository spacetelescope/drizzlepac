"""
Process input to MultiDrizzle/PyDrizzle.

The input can be one of:

    * a python list of files
    * a comma separated string of filenames (including wild card characters)
    * an association table
    * an @file (can have a second column with names of ivm files)

No mixture of instruments is allowed.
No mixture of association tables, @files and regular fits files is allowed.
Files can be in GEIS or MEF format (but not waiver fits).

Runs some sanity checks on the input files.
If necessary converts files to MEF format (this should not be left to makewcs
because 'updatewcs' may be False).
Runs makewcs.
The function 'process_input' returns an association table, ivmlist, output name

The common interface interpreter for MultiDrizzle tasks, 'processCommonInput()',
not only runs 'process_input()' but 'createImageObject()' and 'defineOutput()'
as well to fully setup all inputs for use with the rest of the MultiDrizzle
steps either as stand-alone tasks or internally to MultiDrizzle itself.

"""

from __future__ import division  # confidence high
import datetime
import os
import shutil
import string
import sys

import numpy as np
import pyfits

from stsci.tools import (cfgpars, parseinput, fileutil, asnutil, irafglob,
                         check_files, logutil, mputil, textutil)
from stwcs import updatewcs as uw
from stwcs.wcsutil import altwcs, wcscorr

import wcs_functions
import util
import resetbits
import mdzhandler
import imageObject

log = logutil.create_logger(__name__)

# list parameters which correspond to steps where multiprocessing can be used
parallel_steps = [(3,'driz_separate'),(6,'driz_cr')]

if util.can_parallel:
    import multiprocessing


def setCommonInput(configObj, createOutwcs=True):
    """
    The common interface interpreter for MultiDrizzle tasks which not only runs
    'process_input()' but 'createImageObject()' and 'defineOutput()' as well to
    fully setup all inputs for use with the rest of the MultiDrizzle steps either
    as stand-alone tasks or internally to MultiDrizzle itself.

    Parameters
    ----------
    configObj : object
        configObj instance or simple dictionary of input parameters
    imageObjectList : list of imageObject objects
        list of imageObject instances, 1 for each input exposure
    outwcs : object
        imageObject instance defining the final output frame

    Notes
    -----
    At a minimum, the configObj instance (dictionary) should contain:
        configObj = {'input':None,'output':None,
                    'updatewcs':None}

    If provided, the configObj should contain the values of all the multidrizzle parameters
    as set by the user with TEAL. If no configObj is given, it will retrieve
    the default values automatically.  In either case, the values from the input_dict
    will be merged in with the configObj before being used by the rest of the
    code.

    Examples
    --------
    You can set *createOutwcs=False* for the cases where you only want the
    images processed and no output wcs information in necessary; as in:

    >>>imageObjectList,outwcs = processInput.processCommonInput(configObj)


    """

    if not createOutwcs or not configObj['coeffs']:
        # we're probably just working on single images here
        configObj['updatewcs']=False

    # maybe we can chunk this part up some more so that we can call just the
    # parts we want

    # Interpret input, read and convert and update input files, then return
    # list of input filenames and derived output filename
    asndict, ivmlist, output = process_input(
            configObj['input'], configObj['output'],
            updatewcs=configObj['updatewcs'], wcskey=configObj['wcskey'],
            **configObj['STATE OF INPUT FILES'])

    if not asndict:
        return None, None
    # convert the filenames from asndict into a list of full filenames
    files = [fileutil.buildRootname(f) for f in asndict['order']]

    # interpret MDRIZTAB, if specified, and update configObj accordingly
    # This can be done here because MDRIZTAB does not include values for
    # input, output, or updatewcs.
    if 'mdriztab' in configObj and configObj['mdriztab']:
        mdriztab_dict = mdzhandler.getMdriztabParameters(files)

        # Update configObj with values from mpars
        cfgpars.mergeConfigObj(configObj, mdriztab_dict)

    # Convert interpreted list of input files from process_input into a list
    # of imageObject instances for use by the MultiDrizzle tasks.
    instrpars = configObj['INSTRUMENT PARAMETERS']
    # pass in 'proc_unit' to initialize unit conversions as necessary
    instrpars['proc_unit'] = configObj['proc_unit']

    undistort = True
    if not configObj['coeffs']:
        undistort = False

    # determine whether parallel processing will be performed
    use_parallel = False
    if util.can_parallel:
        # look to see whether steps which can be run using multiprocessing
        # have been turned on
        for stepnum in parallel_steps:
            sname = util.getSectionName(configObj,stepnum[0])
            if configObj[sname][stepnum[1]]:
                use_parallel = True
                break

    # interpret all 'bits' related parameters and convert them to integers
    configObj['resetbits'] = util.interpret_bits_value(configObj['resetbits'])
    step3name = util.getSectionName(configObj,3)
    configObj[step3name]['driz_sep_bits'] = util.interpret_bits_value(
                                        configObj[step3name]['driz_sep_bits'])
    step7name = util.getSectionName(configObj,7)
    configObj[step7name]['final_bits'] = util.interpret_bits_value(
                                        configObj[step7name]['final_bits'])
    # Verify any refimage parameters to be used
    step3aname = util.getSectionName(configObj,'3a')
    if not util.verifyRefimage(configObj[step3aname]['driz_sep_refimage']):
        msg = 'No refimage with WCS found!\n '+\
        ' This could be caused by one of 2 problems:\n'+\
        '   * filename does not specify an extension with a valid WCS.\n'+\
        '   * can not find the file.\n'+\
        'Please check the filename specified in the "refimage" parameter.'
        print textutil.textbox(msg)
        return None,None
    step7aname = util.getSectionName(configObj,'7a')
    if not util.verifyRefimage(configObj[step7aname]['final_refimage']):
        msg = 'No refimage with WCS found!\n '+\
        ' This could be caused by one of 2 problems:\n'+\
        '   * filename does not specify an extension with a valid WCS.\n'+\
        '   * can not find the file.\n'+\
        'Please check the filename specified in the "refimage" parameter.'
        print textutil.textbox(msg)
        return None,None


    # Build imageObject list for all the valid, shift-updated input files
    log.info('-Creating imageObject List as input for processing steps.')
    if 'in_memory' in configObj:
        virtual = configObj['in_memory']
    else:
        virtual = False
    imageObjectList = createImageObjectList(files, instrpars,
                                            group=configObj['group'],
                                            undistort=undistort,
                                            inmemory=virtual)

    # apply context parameter
    applyContextPar(imageObjectList, configObj['context'])

    # reset DQ bits if requested by user
    resetDQBits(imageObjectList, cr_bits_value=configObj['resetbits'])

    # Add info about input IVM files at this point to the imageObjectList
    addIVMInputs(imageObjectList, ivmlist)

    if createOutwcs:
        log.info('-Creating output WCS.')

        # Build output WCS and update imageObjectList with output WCS info
        outwcs = wcs_functions.make_outputwcs(imageObjectList, output,
                                              configObj=configObj)
        outwcs.final_wcs.printwcs()
    else:
        outwcs = None

    try:
        # Provide user with some information on resource usage for this run
        # raises ValueError Exception in interactive mode and user quits
        num_cores = configObj.get('num_cores') if use_parallel else 1

        reportResourceUsage(imageObjectList, outwcs, num_cores)
    except ValueError:
        imageObjectList = None

    return imageObjectList, outwcs


def reportResourceUsage(imageObjectList, outwcs, num_cores,
                        interactive=False):
    """ Provide some information to the user on the estimated resource
    usage (primarily memory) for this run.
    """

    if outwcs is None:
        output_mem = 0
    else:
        if isinstance(outwcs,imageObject.WCSObject):
            owcs = outwcs.final_wcs
        else:
            owcs = outwcs
        output_mem = owcs.naxis1*owcs.naxis2*4*3 # bytes used for output arrays
    img1 = imageObjectList[0]
    numchips = 0
    input_mem = 0
    for img in imageObjectList:
        numchips += img._nmembers # account for group parameter set by user

    # if we have the cpus and s/w, ok, but still allow user to set pool size
    pool_size = util.get_pool_size(num_cores)
    pool_size = pool_size if (numchips >= pool_size) else numchips

    inimg = 0
    chip_mem = 0
    for img in imageObjectList:
        for chip in range(img._numchips):
            cmem = img[chip].image_shape[0]*img[chip].image_shape[1]*4
            inimg += 1
            if inimg < pool_size:
                input_mem += cmem*2
            if chip_mem == 0:
                chip_mem = cmem
    max_mem = (input_mem + output_mem*pool_size + chip_mem*2)//(1024*1024)

    print '*'*80
    print '*'
    print '*  Estimated memory usage:  up to %d Mb.'%(max_mem)
    print '*  Output image size:       %d X %d pixels. '%(owcs.naxis1,owcs.naxis2)
    print '*  Output image file:       ~ %d Mb. '%(output_mem//(1024*1024))
    print '*  Cores used by task:      %d'%(pool_size)
    print '*'
    print '*'*80

    if interactive:
        print 'Continue with processing?'
        k = raw_input("(y)es or (n)o")
        if 'n' in k.lower():
            raise ValueError


def getMdriztabPars(input):
    """ High-level function for getting the parameters from MDRIZTAB

    Used primarily for TEAL interface.
    """
    filelist,output,ivmlist,oldasndict=processFilenames(input,None)

    try:
        mdrizdict = mdzhandler.getMdriztabParameters(filelist)
    except KeyError:
        print 'No MDRIZTAB found for "%s". Parameters remain unchanged.'%(filelist[0])
        mdrizdict = {}

    return mdrizdict

def addIVMInputs(imageObjectList,ivmlist):
    """ Add IVM filenames provided by user to outputNames dictionary for each input imageObject.
    """
    if ivmlist is None:
        return

    for img,ivmname in zip(imageObjectList,ivmlist):
        img.updateIVMName(ivmname)

def checkMultipleFiles(input):
    """ Evaluates the input to determine whether there is 1 or more than 1 valid input file.
    """
    f,i,o,a=buildFileList(input)
    return len(f) > 1

def createImageObjectList(files,instrpars,group=None,
                            undistort=True, inmemory=False):
    """ Returns a list of imageObject instances, 1 for each input image in the list of input filenames.
    """
    imageObjList = []
    mtflag = False
    mt_refimg = None
    for img in files:
        image = _getInputImage(img,group=group)
        image.setInstrumentParameters(instrpars)
        image.compute_wcslin(undistort=undistort)
        if 'MTFLAG' in image._image['PRIMARY'].header:
            # check to see whether we are dealing with moving target observations...
            _keyval = image._image['PRIMARY'].header['MTFLAG']
            if not util.is_blank(_keyval):
                if isinstance(_keyval,bool):
                    mtflag = _keyval
                else:
                    if 'T' in _keyval:
                        mtflag = True
                    else:
                        mtflag = False
        else:
            mtflag = False
        if mtflag is True:
            print "#####\nProcessing Moving Target Observations using reference image as WCS for all inputs!\n#####\n"
            if mt_refimg is None:
                mt_refimg = image
            else:
                image.set_mt_wcs(mt_refimg)
        image.inmemory = inmemory # set flag for inmemory processing
        # Now add (possibly updated) image object to list
        imageObjList.append(image)
    return imageObjList

def applyContextPar(imageObjectList,contextpar):
    """ Apply the value of the parameter `context` to the input list, setting
        the name of the output context image to None if `context` is False
    """
    for img in imageObjectList:
        img.updateContextImage(contextpar)

def _getInputImage (input,group=None):
    """ Factory function to return appropriate imageObject class instance"""
    # extract primary header and SCI,1 header from input image
    sci_ext = 'SCI'
    if group in [None,'']:
        exten = '[sci,1]'
        phdu = pyfits.getheader(input)
    else:
        # change to use pyfits more directly here?
        if group.find(',') > 0:
            grp = group.split(',')
            if grp[0].isalpha():
                grp = (grp[0],int(grp[1]))
            else:
                grp = int(grp[0])
        else:
            grp = int(group)
        phdu = pyfits.getheader(input)
        phdu.extend(pyfits.getheader(input,ext=grp))

    # Extract the instrument name for the data that is being processed by Multidrizzle
    _instrument = phdu['INSTRUME']

    # Determine the instrument detector in use.  NICMOS is a special case because it does
    # not use the 'DETECTOR' keyword.  It instead used 'CAMERA' to identify which of it's
    # 3 camera's is in use.  All other instruments support the 'DETECTOR' keyword.
    if _instrument == 'NICMOS':
        _detector = phdu['CAMERA']
    else:
        try:
            _detector = phdu['DETECTOR']
        except KeyError:
            # using the phdu as set above (pyfits.getheader) is MUCH faster and
            # works for the majority of data; but fileutil handles waivered fits
            phdu = fileutil.getHeader(input+exten)
            _detector = phdu['DETECTOR'] # if this fails, let it throw

    del phdu # just to keep clean
    # Match up the instrument and detector with the right class
    # only importing the instrument modules as needed.
    try:
        if _instrument == 'ACS':
            import acsData
            if _detector == 'HRC': return acsData.HRCInputImage(input,group=group)
            if _detector == 'WFC': return acsData.WFCInputImage(input,group=group)
            if _detector == 'SBC': return acsData.SBCInputImage(input,group=group)
        if _instrument == 'NICMOS':
            import nicmosData
            if _detector == 1: return nicmosData.NIC1InputImage(input)
            if _detector == 2: return nicmosData.NIC2InputImage(input)
            if _detector == 3: return nicmosData.NIC3InputImage(input)


        if _instrument == 'WFPC2':
            import wfpc2Data
            return wfpc2Data.WFPC2InputImage(input,group=group)
        """
            if _detector == 1: return wfpc2Data.PCInputImage(input)
            if _detector == 2: return wfpc2Data.WF2InputImage(input)
            if _detector == 3: return wfpc2Data.WF3InputImage(input)
            if _detector == 4: return wfpc2Data.WF4InputImage(input)
        """
        if _instrument == 'STIS':
            import stisData
            if _detector == 'CCD': return stisData.CCDInputImage(input)
            if _detector == 'FUV-MAMA': return stisData.FUVInputImage(input)
            if _detector == 'NUV-MAMA': return stisData.NUVInputImage(input)
        if _instrument == 'WFC3':
            import wfc3Data
            if _detector == 'UVIS': return wfc3Data.WFC3UVISInputImage(input)
            if _detector == 'IR': return wfc3Data.WFC3IRInputImage(input)

    except ImportError:
        msg = 'No module implemented for '+str(_instrument)+'!'
        raise ValueError,msg
    # If a supported instrument is not detected, print the following error message
    # and raise an exception.
    msg = 'Instrument: ' + str(_instrument) + '/' + str(_detector) + ' not yet supported!'
    raise ValueError, msg

#### Remaining functions support process_input()
def processFilenames(input=None,output=None,infilesOnly=False):
    """Process the input string which contains the input file information and
       return a filelist,output
    """
    ivmlist = None
    oldasndict = None

    if input is None:
        print "No input files provided to processInput"
        raise ValueError

    if (isinstance(input, list) == False) and \
       ('_asn' in input or '_asc' in input) :
        # Input is an association table
        # Get the input files, and run makewcs on them
        oldasndict = asnutil.readASNTable(input, prodonly=infilesOnly)

        if not infilesOnly:
            if output in ["",None,"None"]:
                output = oldasndict['output']

        #filelist = [fileutil.buildRootname(fname) for fname in oldasndict['order']]
        filelist = buildASNList(oldasndict['order'],input)

    elif (not isinstance(input, list)) and \
       (input[0] == '@') :
        # input is an @ file
        f = open(input[1:])
        # Read the first line in order to determine whether
        # IVM files have been specified in a second column...
        line = f.readline()
        f.close()
        # Parse the @-file with irafglob to extract the input filename
        filelist = irafglob.irafglob(input, atfile=util.atfile_sci)
        # If there is a second column...
        if len(line.split()) == 2:
            # ...parse out the names of the IVM files as well
            ivmlist = irafglob.irafglob(input, atfile=util.atfile_ivm)
        if output in ['',None,"None"]:
            if len(filelist) == 1:
                output = fileutil.buildNewRootname(filelist[0])
            else:
                output = 'final'
    else:
        #input is a string or a python list
        try:
            filelist, output = parseinput.parseinput(input, outputname=output)
            if output in ['',None,"None"]:
                if len(filelist) == 1:
                    output = fileutil.buildNewRootname(filelist[0])
                else:
                    output = 'final'
            if not isinstance(input, list):
                filelist.sort()
        except IOError: raise

    # sort the list of input files
    # this ensures the list of input files has the same order on all platforms
    # it can have ifferent order because listdir() uses inode order, not unix type order
    #filelist.sort()



    return filelist, output, ivmlist, oldasndict


def process_input(input, output=None, ivmlist=None, updatewcs=True,
                  prodonly=False,  wcskey=None, **workinplace):
    """
    Create the full input list of filenames after verifying and converting
    files as needed.
    """

    newfilelist, ivmlist, output, oldasndict = buildFileList(
            input, output=output, ivmlist=ivmlist, wcskey=wcskey,
            updatewcs=updatewcs, **workinplace)

    if not newfilelist:
        buildEmptyDRZ(input, output)
        return None, None, output

    # run all WCS updating -- Now done in buildFileList
    #pydr_input = _process_input_wcs(newfilelist, wcskey, updatewcs)
    pydr_input = newfilelist

    # AsnTable will handle the case when output==None
    if not oldasndict:# and output is not None:
        oldasndict = asnutil.ASNTable(pydr_input, output=output)
        oldasndict.create()

    asndict = update_member_names(oldasndict, pydr_input)

    # Build output filename
    drz_extn = '_drz.fits'
    for img in newfilelist:
        # special case logic to automatically recognize when _flc.fits files
        # are provided as input and produce a _drc.fits file instead
        if '_flc.fits' in img:
            drz_extn = '_drc.fits'
            break

    if output in [None,'']:
        output = fileutil.buildNewRootname(asndict['output'],
                                           extn=drz_extn)
    else:
        if '.fits' in output.lower():
            pass
        elif drz_extn[:4] not in output.lower():
            output = fileutil.buildNewRootname(output, extn=drz_extn)

    log.info('Setting up output name: %s' % output)

    return asndict, ivmlist, output


def _process_input_wcs(infiles, wcskey, updatewcs):
    """
    This is a subset of process_input(), for internal use only.  This is the
    portion of input handling which sets/updates WCS data, and is a performance
    hit - a target for parallelization. Returns the expanded list of filenames.
    """

    # Run parseinput though it's likely already been done in processFilenames
    outfiles = parseinput.parseinput(infiles)[0]

    # Disable parallel processing here for now until hardware I/O gets "wider".
    # Since this part is IO bound, parallelizing doesn't help more than a little
    # in most cases, and may actually slow this down on some desktop nodes.
#   cfgval_num_cores = None # get this from paramDict
#   pool_size = util.get_pool_size(cfgval_num_cores, len(outfiles))
    pool_size = 1

    # do the WCS updating
    if wcskey in ['', ' ', 'INDEF', None]:
        if updatewcs:
            log.info('Updating input WCS using "updatewcs"')
    else:
        log.info('Resetting input WCS to be based on WCS key = %s' % wcskey)

    if pool_size > 1:
        log.info('Executing %d parallel workers' % pool_size)
        subprocs = []
        for fname in outfiles:
            p = multiprocessing.Process(target=_process_input_wcs_single,
                name='processInput._process_input_wcs()', # for err msgs
                args=(fname, wcskey, updatewcs) )
            subprocs.append(p)
        mputil.launch_and_wait(subprocs, pool_size) # blocks till all done
    else:
        log.info('Executing serially')
        for fname in outfiles:
            _process_input_wcs_single(fname, wcskey, updatewcs)

    return outfiles


def _process_input_wcs_single(fname, wcskey, updatewcs):
    """
    See docs for _process_input_wcs.
    This is separated to be spawned in parallel.
    """
    if wcskey in ['', ' ', 'INDEF', None]:
        if updatewcs:
            uw.updatewcs(fname, checkfiles=False)
    else:
        numext = fileutil.countExtn(fname)
        extlist = []
        for extn in xrange(1, numext + 1):
            extlist.append(('SCI', extn))
        if wcskey in string.uppercase:
            wkey = wcskey
            wname = ' '
        else:
            wname = wcskey
            wkey = ' '
        altwcs.restoreWCS(fname, extlist, wcskey=wkey, wcsname=wname)
    # make an asn table at the end
    # Make sure there is a WCSCORR table for each input image
    if wcskey not in ['', ' ', 'INDEF', None] or updatewcs:
        wcscorr.init_wcscorr(fname)


def buildFileList(input, output=None, ivmlist=None,
                wcskey=None, updatewcs=True, **workinplace):
    """
    Builds a file list which has undergone various instrument-specific
    checks for input to MultiDrizzle, including splitting STIS associations.
    """
    filelist,output,ivmlist,oldasndict=processFilenames(input,output)

    # verify that all input images specified can be updated as needed
    filelist = util.verifyFilePermissions(filelist)
    if filelist is None or len(filelist) == 0:
        return None, None, None, None

    manageInputCopies(filelist,**workinplace)

    # Check format of FITS files - convert Waiver/GEIS to MEF if necessary
    filelist, ivmlist = check_files.checkFITSFormat(filelist, ivmlist)

    # check for non-polynomial distortion correction
    filelist = checkDGEOFile(filelist)

    # run all WCS updating
    updated_input = _process_input_wcs(filelist, wcskey, updatewcs)

    newfilelist, ivmlist = check_files.checkFiles(updated_input, ivmlist)

    return newfilelist,ivmlist,output,oldasndict


def buildASNList(rootnames, asnname):
    """
    Return the list of filenames for a given set of rootnames
    """

    # Recognize when multiple valid inputs with the same rootname are present
    # this would happen when both CTE-corrected (_flc) and non-CTE-corrected (_flt)
    # products are in the same directory as an ASN table
    filelist, duplicates = checkForDuplicateInputs(rootnames)

    if duplicates:
        # Build new ASN tables for each set of input files
        origasn = changeSuffixinASN(asnname, 'flt')
        dupasn = changeSuffixinASN(asnname, 'flc')

        errstr = 'ERROR:\nMultiple valid input files found:\n'
        for fname, dname in zip(filelist, duplicates):
            errstr += '    %s    %s\n' % (fname, dname)
        errstr += ('\nNew association files have been generated for each '
                   'version of these files.\n    %s\n    %s\n\nPlease '
                   're-start astrodrizzle using of these new ASN files or '
                   'use widlcards for the input to only select one type of '
                   'input file.' % (dupasn, origasn))

        print >> sys.stderr, textutil.textbox(errstr)

        # generate new ASN files for each case,
        # report this case of duplicate inputs to the user then quit
        raise ValueError

    return filelist


def changeSuffixinASN(asnfile, suffix):
    """
    Create a copy of the original asn file and change the name of all members
    to include the suffix.
    """
    # Start by creating a new name for the ASN table
    _new_asn = asnfile.replace('_asn.fits','_'+suffix+'_asn.fits')
    if os.path.exists(_new_asn) == True:
        os.remove(_new_asn)
    # copy original ASN table to new table
    shutil.copy(asnfile,_new_asn)

    # Open up the new copy and convert all MEMNAME's to lower-case
    fasn = pyfits.open(_new_asn,'update')
    for i in xrange(len(fasn[1].data)):
        if "prod" not in fasn[1].data[i].field('MEMTYPE').lower():
            fasn[1].data[i].setfield('MEMNAME',fasn[1].data[i].field('MEMNAME')+'_'+suffix)
    fasn.close()

    return _new_asn


def checkForDuplicateInputs(rootnames):
    """
    Check input files specified in ASN table for duplicate versions with
    multiple valid suffixes (_flt and _flc, for example).
    """

    flist = []
    duplist = []

    for fname in rootnames:
        # Look for any recognized CTE-corrected products
        f1 = fileutil.buildRootname(fname,ext=['_flc.fits'])
        f2 = fileutil.buildRootname(fname)
        flist.append(f2)
        if os.path.exists(f1) and f1 != f2:
            # More than 1 valid input found for this rootname
            duplist.append(f1)

    return flist,duplist


def runmakewcs(input):
    """
    Runs make wcs and recomputes the WCS keywords

    Parameters
    ----------
    input : str or list of str
        a list of files

    Returns
    -------
    output : list of str
        returns a list of names of the modified files
        (For GEIS files returns the translated names.)
    """
    newNames = uw.updatewcs(input, checkfiles=False)
    #newNames = makewcs.run(input)
    return newNames


def resetDQBits(imageObjectList, cr_bits_value=4096):
    """Reset the CR bit in each input image's DQ array"""

    if cr_bits_value > 0:
        for img in imageObjectList:
            for chip in range(1,img._numchips+1,1):
                sci_chip = img._image[img.scienceExt,chip]
                resetbits.reset_dq_bits(sci_chip.dqfile, cr_bits_value,
                                        extver=chip, extname=sci_chip.dq_extn)


def update_member_names(oldasndict, pydr_input):
    """
    Update names in a member dictionary.

    Given an association dictionary with rootnames and a list of full
    file names, it will update the names in the member dictionary to
    contain '_*' extension. For example a rootname of 'u9600201m' will
    be replaced by 'u9600201m_c0h' making sure that a MEf file is passed
    as an input and not the corresponding GEIS file.
    """

    omembers = oldasndict['members'].copy()
    nmembers = {}
    translated_names = [f.split('.fits')[0] for f in pydr_input]

    newkeys = [fileutil.buildNewRootname(file) for file in pydr_input]
    keys_map = zip(newkeys, pydr_input)

    iter = omembers.iteritems()
    while True:
        try:
            okey, oval = iter.next()
            if okey in newkeys:
                nkey = pydr_input[newkeys.index(okey)]
                nmembers[nkey.split('.fits')[0]] = oval
        except StopIteration:
            break
    oldasndict.pop('members')
    # replace should be always True to cover the case when flt files were removed
    # and the case when names were translated

    oldasndict.update(members=nmembers, replace=True)
    oldasndict['order'] = translated_names
    return oldasndict


def manageInputCopies(filelist, **workinplace):
    """
    Creates copies of all input images in a sub-directory.

    The copies are made prior to any processing being done to the images at all,
    including updating the WCS keywords. If there are already copies present,
    they will NOT be overwritten, but instead will be used to over-write the
    current working copies.
    """

    # Find out what directory is being used for processing
    workingdir = os.getcwd()
    # Only create sub-directory for copies of inputs, if copies are requested
    # Create name of sub-directory for copies
    origdir = os.path.join(workingdir,'OrIg_files')
    if workinplace['overwrite'] or workinplace['preserve']:
        # if sub-directory does not exist yet, create it
        if not os.path.exists(origdir):
            os.mkdir(origdir)

    printMsg = True
    # check to see if copies already exist for each file
    for fname in filelist:
        copymade = False # If a copy is made, no need to restore
        copyname = os.path.join(origdir,fname)
        short_copyname = os.path.join('OrIg_files',fname)
        if workinplace['overwrite']:
            print 'Forcibly archiving original of: ',fname, 'as ',short_copyname
            # make a copy of the file in the sub-directory
            if os.path.exists(copyname): os.chmod(copyname, 0666)
            shutil.copy(fname,copyname)
            os.chmod(copyname,0444) # make files read-only
            if printMsg:
                print '\nTurning OFF "preserve" and "restore" actions...\n'
                printMsg = False # We only need to print this one time...
            copymade = True

        if (workinplace['preserve'] and not os.path.exists(copyname)) \
                and not workinplace['overwrite']:
            # Preserving a copy of the input, but only if not already archived
            print 'Preserving original of: ',fname, 'as ',short_copyname
            # make a copy of the file in the sub-directory
            shutil.copy(fname,copyname)
            os.chmod(copyname,0444) # make files read-only
            copymade = True

        if 'restore' in workinplace and not copymade:
            if (os.path.exists(copyname) and workinplace['restore']) and not workinplace['overwrite']:
                print 'Restoring original input for ',fname,' from ',short_copyname
                # replace current files with original version
                os.chmod(fname, 0666)
                shutil.copy(copyname, fname)
                os.chmod(fname, 0666)


def buildEmptyDRZ(input, output):
    """
    Create an empty DRZ file.

    This module creates an empty DRZ file in a valid FITS format so that the HST
    pipeline can handle the Multidrizzle zero expossure time exception
    where all data has been excluded from processing.

    Parameters
    ----------
    input : str
        filename of the initial input to process_input
    output : str
        filename of the default empty _drz.fits file to be generated

    """

    # Identify the first input image
    inputfile = parseinput.parseinput(input)[0]
    if not inputfile:
        print >> sys.stderr, '\n******* ERROR *******'
        print >> sys.stderr, (
               'No input file found!  Check specification of parameter '
               '"input". ')
        print >> sys.stderr, 'Quitting...'
        print >> sys.stderr, '******* ***** *******\n'
        return # raise IOError, "No input file found!"

    # Set up output file here...
    if output is None:
        if len(input) == 1:
            oname = fileutil.buildNewRootname(input[0])
        else:
            oname = 'final'
        output = fileutil.buildNewRootname(oname, extn='_drz.fits')
    else:
        if 'drz' not in output:
            output = fileutil.buildNewRootname(output, extn='_drz.fits')

    log.info('Setting up output name: %s' % output)

    # Open the first image (of the excludedFileList?) to use as a template to build
    # the DRZ file.
    try :
        log.info('Building empty DRZ file from %s' % inputfile[0])
        img = pyfits.open(inputfile[0])
    except:
        raise IOError('Unable to open file %s \n' % inputfile)

    # Create the fitsobject
    fitsobj = pyfits.HDUList()
    # Copy the primary header
    hdu = img[0].copy()
    fitsobj.append(hdu)

    # Modify the 'NEXTEND' keyword of the primary header to 3 for the
    #'sci, wht, and ctx' extensions of the newly created file.
    fitsobj[0].header['NEXTEND'] = 3

    # Create the 'SCI' extension
    hdu = pyfits.ImageHDU(header=img['sci', 1].header.copy())
    hdu.header['EXTNAME'] = 'SCI'
    fitsobj.append(hdu)

    # Create the 'WHT' extension
    hdu = pyfits.ImageHDU(header=img['sci', 1].header.copy())
    hdu.header['EXTNAME'] = 'WHT'
    fitsobj.append(hdu)

    # Create the 'CTX' extension
    hdu = pyfits.ImageHDU(header=img['sci', 1].header.copy())
    hdu.header['EXTNAME'] = 'CTX'
    fitsobj.append(hdu)

    # Add HISTORY comments explaining the creation of this file.
    fitsobj[0].header.add_history("** AstroDrizzle has created this empty "
                                  "DRZ product because**")
    fitsobj[0].header.add_history("** all input images were excluded from "
                                  "processing.**")


    # Change the filename in the primary header to reflect the name of the output
    # filename.
    fitsobj[0].header['FILENAME'] = str(output)  # +"_drz.fits"

    # Change the ROOTNAME keyword to the ROOTNAME of the output PRODUCT
    fitsobj[0].header['ROOTNAME'] = str(output.split('_drz.fits')[0])
    # Modify the ASN_MTYP keyword to contain "PROD-DTH" so it can be properly
    # ingested into the archive catalog.

    # stis has this keyword in the [1] header, so I am directing the code
    #t o first look in the primary, then the 1
    try:
        fitsobj[0].header['ASN_MTYP'] = 'PROD-DTH'
    except:
        fitsobj[1].header['ASN_MTYP'] = 'PROD-DTH'

    # If the file is already on disk delete it and replace it with the
    # new file
    dirfiles = os.listdir(os.curdir)
    if dirfiles.count(output) > 0:
        os.remove(output)
        log.info("       Replacing %s..." % output)

    # Write out the empty DRZ file
    fitsobj.writeto(output)

    print >> sys.stderr, textutil.textbox(
        'ERROR:\nAstroDrizzle has created an empty DRZ product because all '
        'input images were excluded from processing or a user requested the '
        'program to stop.') + '\n'

    return


def checkDGEOFile(filenames):
    """
    Verify that input file has been updated with NPOLFILE

    This function checks for the presence of 'NPOLFILE' kw in the primary header
    when 'DGEOFILE' kw is present and valid (i.e. 'DGEOFILE' is not blank or 'N/A').
    It handles the case of science files downloaded from the archive before the new
    software was installed there.
    If 'DGEOFILE' is present and 'NPOLFILE' is missing, print a message and let the user
    choose whether to (q)uit and update the headers or (c)ontinue and run astrodrizzle
    without the non-polynomial correction.
    'NPOLFILE' will be populated in the pipeline before astrodrizzle is run.

    In the case of WFPC2 the old style dgeo files are used to create detector to image
    correction at runtime.

    Parameters
    ----------
    filenames : list of str
        file names of all images to be checked

    """

    msg = """
            A 'DGEOFILE' keyword is present in the primary header but 'NPOLFILE' keyword was not found.
            This version of the software uses a new format for the residual distortion DGEO files.
            Please consult the instrument web pages for which reference files to download.
            A small (new style) dgeofile is needed ('_npl.fits' extension) and possibly a
            detector to image correction file ('_d2i.fits' extension).
            The names of these files must be added to the primary header either using the task XXXX
            or manually, for example:

            hedit %s[0] npolfile fname_npl.fits add+
            hedit %s[0] d2imfile fname_d2i.fits add+

            where fname_npl.fits is the name of the new style dgeo file and fname_d2i.fits is
            the name of the detector to image correction. After adding these keywords to the
            primary header, updatewcs must be run to update the science files:

            from stwcs import updatewcs
            updatewcs.updatewcs(%s)

            Alternatively you may choose to run astrodrizzle without DGEO and detector to image correction.

            To stop astrodrizzle and update the dgeo files, type 'q'.
            To continue running astrodrizzle without the non-polynomial distortion correction, type 'c':
            """

    for inputfile in filenames:
        if pyfits.getval(inputfile, 'INSTRUME') == 'WFPC2':
            update_wfpc2_d2geofile(inputfile)
        else:
            try:
                dgeofile = pyfits.getval(inputfile, 'DGEOFILE')
            except KeyError:
                continue
            if dgeofile not in ["N/A", "n/a", ""]:
                message = msg % (inputfile, inputfile, inputfile)
                try:
                    npolfile = pyfits.getval(inputfile, 'NPOLFILE')
                except KeyError:
                    ustop = userStop(message)
                    while ustop == None:
                        ustop = userStop(message)
                    if ustop == True:
                        return None
                    elif ustop == False:
                        pass
    return filenames

def userStop(message):
    user_input = raw_input(message)
    if user_input == 'q':
        return True
    elif user_input == 'c':
        return False
    else:
        return None


def update_wfpc2_d2geofile(filename, fhdu=None):
    """
    Creates a D2IMFILE from the DGEOFILE for a WFPC2 image (input), and
    modifies the header to reflect the new usage.
    """

    close_fhdu = False
    if fhdu is None:
        fhdu = fileutil.openImage(filename, mode='update')
        close_fhdu = True

    dgeofile = fhdu['PRIMARY'].header.get('DGEOFILE', None)
    if dgeofile not in [None, "N/A", "", " "]:
        log.info('Converting DGEOFILE %s into D2IMFILE...' % dgeofile)
        rootname = filename[:filename.find('.fits')]
        d2imfile = convert_dgeo_to_d2im(dgeofile,rootname)
        fhdu['PRIMARY'].header.update('ODGEOFIL', dgeofile)
        fhdu['PRIMARY'].header.update('DGEOFILE', 'N/A')
        fhdu['PRIMARY'].header.update('D2IMFILE', d2imfile)
    else:
        d2imfile = None
        fhdu['PRIMARY'].header.update('DGEOFILE', 'N/A')
        fhdu['PRIMARY'].header.update('D2IMFILE', 'N/A')

    # Only close the file handle if opened in this function
    if close_fhdu:
        fhdu.close()

    # return the d2imfile name so that calling routine can keep
    # track of the new file created and delete it later if necessary
    # (multidrizzle clean=True mode of operation)
    return d2imfile


def convert_dgeo_to_d2im_OLD(dgeofile,output,clobber=True):
    """ Routine that converts the WFPC2 DGEOFILE into a D2IMFILE.
    """
    dgeo = fileutil.openImage(dgeofile)
    outname = output+'_d2im.fits'

    util.removeFileSafely(outname)

    scihdu = pyfits.ImageHDU(data=dgeo['dy',1].data[:,0])
    dgeo.close()
    # add required keywords for D2IM header
    scihdu.header.update('EXTNAME','DY',comment='Extension name')
    scihdu.header.update('EXTVER',1,comment='Extension version')
    pyfits_str = 'PYFITS Version '+str(pyfits.__version__)
    scihdu.header.update('ORIGIN',pyfits_str,comment='FITS file originator')
    scihdu.header.update('INHERIT',False,comment='Inherits global header')

    dnow = datetime.datetime.now()
    scihdu.header.update('DATE',str(dnow).replace(' ','T'),comment='Date FITS file was generated')

    scihdu.header.update('CRPIX1',0,comment='Distortion array reference pixel')
    scihdu.header.update('CDELT1',0,comment='Grid step size in first coordinate')
    scihdu.header.update('CRVAL1',0,comment='Image array pixel coordinate')
    scihdu.header.update('CRPIX2',0,comment='Distortion array reference pixel')
    scihdu.header.update('CDELT2',0,comment='Grid step size in second coordinate')
    scihdu.header.update('CRVAL2',0,comment='Image array pixel coordinate')
    
    d2imhdu = pyfits.HDUList()
    d2imhdu.append(pyfits.PrimaryHDU())
    d2imhdu.append(scihdu)
    d2imhdu.writeto(outname)
    d2imhdu.close()

    return outname

def convert_dgeo_to_d2im(dgeofile,output,clobber=True):
    """ Routine that converts the WFPC2 DGEOFILE into a D2IMFILE.
    """
    dgeo = fileutil.openImage(dgeofile)
    outname = output+'_d2im.fits'

    util.removeFileSafely(outname)
    data = np.array([dgeo['dy',1].data[:,0]])
    scihdu = pyfits.ImageHDU(data=data)
    dgeo.close()
    # add required keywords for D2IM header
    scihdu.header.update('EXTNAME','DY',comment='Extension name')
    scihdu.header.update('EXTVER',1,comment='Extension version')
    pyfits_str = 'PYFITS Version '+str(pyfits.__version__)
    scihdu.header.update('ORIGIN',pyfits_str,comment='FITS file originator')
    scihdu.header.update('INHERIT',False,comment='Inherits global header')

    dnow = datetime.datetime.now()
    scihdu.header.update('DATE',str(dnow).replace(' ','T'),comment='Date FITS file was generated')

    scihdu.header.update('CRPIX1',0,comment='Distortion array reference pixel')
    scihdu.header.update('CDELT1',1,comment='Grid step size in first coordinate')
    scihdu.header.update('CRVAL1',0,comment='Image array pixel coordinate')
    scihdu.header.update('CRPIX2',0,comment='Distortion array reference pixel')
    scihdu.header.update('CDELT2',1,comment='Grid step size in second coordinate')
    scihdu.header.update('CRVAL2',0,comment='Image array pixel coordinate')

    phdu = pyfits.PrimaryHDU()
    phdu.header.update('INSTRUME', 'WFPC2')
    d2imhdu = pyfits.HDUList()
    d2imhdu.append(phdu)
    scihdu.header.update('DETECTOR', 1, comment='CCD number of the detector: PC 1, WFC 2-4 ')
    d2imhdu.append(scihdu.copy())
    scihdu.header.update('EXTVER', 2, comment='Extension version')
    scihdu.header.update('DETECTOR', 2, comment='CCD number of the detector: PC 1, WFC 2-4 ')
    d2imhdu.append(scihdu.copy())
    scihdu.header.update('EXTVER', 3, comment='Extension version')
    scihdu.header.update('DETECTOR', 3, comment='CCD number of the detector: PC 1, WFC 2-4 ')
    d2imhdu.append(scihdu.copy())
    scihdu.header.update('EXTVER', 4, comment='Extension version')
    scihdu.header.update('DETECTOR', 4, comment='CCD number of the detector: PC 1, WFC 2-4 ')
    d2imhdu.append(scihdu.copy())
    d2imhdu.writeto(outname)
    d2imhdu.close()

    return outname

def _setDefaults(input_dict={}):
    """ Define full set of default values for unit-testing this module.[OBSOLETE]"""
    paramDict = {
        'input':'*flt.fits',
        'output':None,
        'mdriztab':None,
        'refimage':None,
        'runfile':None,
        'workinplace':False,
        'updatewcs':True,
        'proc_unit':'native',
        'coeffs':True,
        'context':False,
        'clean':True,
        'group':None,
        'ra':None,
        'dec':None,
        'build':True,
        'gain':None,
        'gnkeyword':None,
        'readnoise':None,
        'rnkeyword':None,
        'exptime':None,
        'expkeyword':None,
        'crbitval':4096,
        'static':True,
        'static_sig':4.0,
        'skysub':True,
        'skywidth':0.1,
        'skystat':"median",
        'skylower':None,
        'skyupper':None,
        'skyclip':5,
        'skylsigma':4.0,
        'skyusigma':4.0,
        'driz_separate':True,
        'driz_sep_outnx':None,
        'driz_sep_outny':None,
        'driz_sep_kernel':'turbo',
        'driz_sep_scale':None,
        'driz_sep_pixfrac':1.0,
        'driz_sep_rot':None,
        'driz_sep_fillval':None,
        'driz_sep_bits':0,
        'median':True,
        'median_newmasks':True,
        'combine_type':"minmed",
        'combine_nsigma':"4 3",
        'combine_nlow':0,
        'combine_nhigh':1,
        'combine_lthresh':None,
        'combine_hthresh':None,
        'combine_grow':1,
        'blot':True,
        'blot_interp':'poly5',
        'blot_sinscl':1.0,
        'driz_cr':True,
        'driz_cr_corr':False,
        'driz_cr_snr':"3.5 3.0",
        'driz_cr_scale':"1.2 0.7",
        'driz_cr_cteg':0,
        'driz_cr_grow':1,
        'driz_combine':True,
        'final_wht_type':"EXP",
        'final_outnx':None,
        'final_outny':None,
        'final_kernel':'square',
        'final_scale':None,
        'final_pixfrac':1.0,
        'final_rot':None,
        'final_fillval':None,
        'final_bits':0}

    paramDict.update(input_dict)

    print '\nUser Input Parameters for Init Step:'
    util.printParams(paramDict)

    return paramDict
