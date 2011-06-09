from __future__ import division # confidence high
import datetime
import os
import shutil

import pyfits

from pytools import parseinput, fileutil, asnutil, irafglob, check_files
from stwcs import updatewcs
from stwcs.wcsutil import altwcs, wcscorr

import wcs_functions
import util
import resetbits
import mdzhandler

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

def setCommonInput(configObj,createOutwcs=True):
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
    if not createOutwcs:
        configObj['updatewcs']=False #we're probably just working on single images here

    #maybe we can chunk this part up some more so that we can call just the parts we want

    # Interpret input, read and convert and update input files, then return
    # list of input filenames and derived output filename
    asndict,ivmlist,output = process_input(configObj['input'], configObj['output'],
            updatewcs=configObj['updatewcs'], wcskey=configObj['wcskey'],
            **configObj['STATE OF INPUT FILES'])

    if not asndict:
        return None, None
    # convert the filenames from asndict into a list of full filenames
    files = [fileutil.buildRootname(f) for f in asndict['order']]

    # interpret MDRIZTAB, if specified, and update configObj accordingly
    # This can be done here because MDRIZTAB does not include values for
    # input, output, or updatewcs.
    if configObj['mdriztab']:
        mdriztab_dict = mdzhandler.getMdriztabParameters(files)
        # Update configObj with values from mpars
        util.mergeConfigObj(configObj,mdriztab_dict)

    # Convert interpreted list of input files from process_input into a list
    # of imageObject instances for use by the MultiDrizzle tasks.
    instrpars = configObj['INSTRUMENT PARAMETERS']
    # pass in 'proc_unit' to initialize unit conversions as necessary
    instrpars['proc_unit'] = configObj['proc_unit']

    # Build imageObject list for all the valid, shift-updated input files
    print '\n-Creating imageObject List as input for processing steps.\n'
    imageObjectList = createImageObjectList(files,instrpars,group=configObj['group'])

    # apply context parameter
    applyContextPar(imageObjectList,configObj['context'])

    # reset DQ bits if requested by user
    resetDQBits(imageObjectList,cr_bits_value=configObj['resetbits'])

    # Add info about input IVM files at this point to the imageObjectList
    addIVMInputs(imageObjectList,ivmlist)

    if(createOutwcs):
        print '\n-Creating output WCS.\n'
        # Build output WCS and update imageObjectList with output WCS info
        outwcs = wcs_functions.make_outputwcs(imageObjectList,output,configObj=configObj)
        return imageObjectList,outwcs
    else:
        return imageObjectList,None

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

def createImageObjectList(files,instrpars,group=None):
    """ Returns a list of imageObject instances, 1 for each input image in the list of input filenames.
    """
    imageObjList = []
    mtflag = False
    mt_refimg = None
    for img in files:
        image = _getInputImage(img,group=group)
        image.setInstrumentParameters(instrpars)
        if image._image['PRIMARY'].header.has_key('MTFLAG'):
            # check to see whether we are dealing with moving target observations...
            _keyval = image._image['PRIMARY'].header['MTFLAG']
            if _keyval not in [""," ", None, "INDEF"] and "T" in _keyval: mtflag = True
        else:
            mtflag = False
        if mtflag is True:
            print "#####\nProcessing Moving Target Observations using reference image as WCS for all inputs!\n#####\n"
            if mt_refimg is None:
                mt_refimg = image
            else:
                image.set_mt_wcs(mt_refimg)

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
    if group in [None,'']:
        exten = '[sci,1]'
    else:
        if group.find(',') > 0:
            grp = group.split(',')
        else:
            grp = ['SCI',int(group)]
        fimg = fileutil.openImage(input)
        exten = '['+str(fileutil.findExtname(fimg,extname=grp[0],extver=grp[1]))+']'
        fimg.close()

    phdu = fileutil.getHeader(input+exten)

    # Extract the instrument name for the data that is being processed by Multidrizzle
    _instrument = phdu['INSTRUME']

    # Determine the instrument detector in use.  NICMOS is a special case because it does
    # not use the 'DETECTOR' keyword.  It instead used 'CAMERA' to identify which of it's
    # 3 camera's is in use.  All other instruments support the 'DETECTOR' keyword.
    if (_instrument == 'NICMOS'):
        _detector = phdu['CAMERA']
    else:
        _detector = phdu['DETECTOR']

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

        filelist = [fileutil.buildRootname(fname) for fname in oldasndict['order']]

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


    return filelist,output,ivmlist,oldasndict

def process_input(input, output=None, ivmlist=None, updatewcs=True, prodonly=False,  wcskey=None, **workinplace):
    """ Create the full input list of filenames after verifying and converting
        files as needed.
    """
    newfilelist,ivmlist,output,oldasndict = buildFileList(input,output=output,ivmlist=ivmlist,**workinplace)

    if not newfilelist or len(newfilelist) == 0:
        buildEmptyDRZ(input,output)
        return None, None, output

    # check for non-polynomial distortion correction
    newfilelist = checkDGEOFile(newfilelist)
    if newfilelist == None:

        return None, None, None

    #make an asn table at the end
    if wcskey not in ['',' ','INDEF',None] or updatewcs:
        # Make sure there is a WCSCORR table for each input image
        for img in newfilelist:
            wcscorr.init_wcscorr(img)

    if wcskey in ['',' ','INDEF',None]:
        if updatewcs:
            print 'Updating input WCS using "updatewcs"'
            pydr_input = runmakewcs(newfilelist)
            # update WCSCORR table with newly updated WCS
            for imgname in newfilelist:
                idcname = os.path.split(fileutil.osfn(pyfits.getval(imgname,'idctab')))[1]
                idcname = idcname[:idcname.find('_idc.fits')]
                wcscorr.archive_wcs_file(imgname,wcs_id='IDC_'+idcname)
        else:
            pydr_input = newfilelist
    else:
        print 'Resetting input WCS to be based on WCS key = ',wcskey
        for fname in newfilelist:
            numext = fileutil.countExtn(fname)
            extlist = []
            for extn in xrange(1,numext+1):
                extlist.append(('SCI',extn))
            altwcs.restoreWCS(fname,extlist,wcskey=wcskey,clobber=True)
        pydr_input = newfilelist

    # AsnTable will handle the case when output==None
    if not oldasndict:# and output is not None:
        oldasndict = asnutil.ASNTable(pydr_input, output=output)
        oldasndict.create()

    asndict = update_member_names(oldasndict, pydr_input)

    # Build output filename
    if output in [None,'']:
        output = fileutil.buildNewRootname(asndict['output'],extn='_drz.fits')
    else:
        if 'drz' not in output:
            output = fileutil.buildNewRootname(output,extn='_drz.fits')

    print 'Setting up output name: ',output

    return asndict, ivmlist, output

def buildFileList(input, output=None, ivmlist=None,**workinplace):
    """
    Builds a file list which has undergone various instrument-specific
    checks for input to MultiDrizzle, including splitting STIS associations.
    """
    filelist,output,ivmlist,oldasndict=processFilenames(input,output)

    manageInputCopies(filelist,**workinplace)

    newfilelist, ivmlist = check_files.checkFiles(filelist, ivmlist)

    return newfilelist,ivmlist,output,oldasndict

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
    newNames = updatewcs.updatewcs(input,checkfiles=False)
    #newNames = makewcs.run(input)
    return newNames


def resetDQBits(imageObjectList,cr_bits_value=4096):
    """ Reset the CR bit in each input image's DQ array
    """
    if cr_bits_value > 0:
        for img in imageObjectList:
            for chip in range(1,img._numchips+1,1):
                sci_chip = img._image[img.scienceExt,chip]
                resetbits.reset_dq_bits(sci_chip.dqfile,cr_bits_value,extver=chip,extname=sci_chip.dq_extn)

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
            okey,oval = iter.next()
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
    # Create name of sub-directory for copies
    origdir = os.path.join(workingdir,'OrIg_files')
    # if sub-directory does not exist yet, create it
    if not os.path.exists(origdir):
        os.mkdir(origdir)

    printMsg = True
    # check to see if copies already exist for each file
    for fname in filelist:
        copyname = os.path.join(origdir,fname)
        if workinplace['overwrite']:
            print 'Forcibly archiving original of: ',fname, 'as ',copyname
            # make a copy of the file in the sub-directory
            if os.path.exists(copyname): os.chmod(copyname, 0666)
            shutil.copy(fname,copyname)
            os.chmod(copyname,0444) # make files read-only
            if printMsg:
                print '\nTurning OFF "preserve" and "restore" actions...\n'
                printMsg = False # We only need to print this one time...
                
        if (workinplace['preserve'] and not os.path.exists(copyname)) \
                and not workinplace['overwrite']:
            # Preserving a copy of the input, but only if not already archived            
            print 'Preserving original of: ',fname, 'as ',copyname
            # make a copy of the file in the sub-directory
            shutil.copy(fname,copyname)
            os.chmod(copyname,0444) # make files read-only
            
        if 'restore' in workinplace:
            if (os.path.exists(copyname) and workinplace['restore']) and not workinplace['overwrite']:
                print 'Restoring original input for ',fname,' from ',copyname
                # replace current files with original version
                os.chmod(fname,0666)
                shutil.copy(copyname,fname)
                os.chmod(fname,0666)

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
    if len(inputfile) == 0:
        print '\n******* ERROR *******'
        print 'No input file found!  Check specification of parameter "input". '
        print 'Quitting...'
        print '******* ***** *******\n'
        return #raise IOError, "No input file found!"

    # Set up output file here...
    if output == None:
        if len(input) == 1:
            oname = fu.buildNewRootname(input[0])
        else:
            oname = 'final'
        output = fileutil.buildNewRootname(oname,extn='_drz.fits')
    else:
        if 'drz' not in output:
            output = fileutil.buildNewRootname(output,extn='_drz.fits')

    print 'Setting up output name: ',output

    # Open the first image (of the excludedFileList?) to use as a template to build
    # the DRZ file.
    try :
        img = pyfits.open(inputfile[0])
    except:
        raise IOError, 'Unable to open file %s \n' %inputfile


    # Create the fitsobject
    fitsobj = pyfits.HDUList()
    # Copy the primary header
    hdu = img[0].copy()
    fitsobj.append(hdu)

    # Modify the 'NEXTEND' keyword of the primary header to 3 for the
    #'sci, wht, and ctx' extensions of the newly created file.
    fitsobj[0].header['NEXTEND'] = 3

    # Create the 'SCI' extension
    hdu = pyfits.ImageHDU(header=img['sci',1].header.copy(),data=None)
    hdu.header['EXTNAME'] = 'SCI'
    fitsobj.append(hdu)

    # Create the 'WHT' extension
    hdu = pyfits.ImageHDU(header=img['sci',1].header.copy(),data=None)
    hdu.header['EXTNAME'] = 'WHT'
    fitsobj.append(hdu)

    # Create the 'CTX' extension
    hdu = pyfits.ImageHDU(header=img['sci',1].header.copy(),data=None)
    hdu.header['EXTNAME'] = 'CTX'
    fitsobj.append(hdu)

    # Add HISTORY comments explaining the creation of this file.
    fitsobj[0].header.add_history("** Multidrizzle has created this empty DRZ product because**")
    fitsobj[0].header.add_history("** all input images were excluded from processing.**")


    # Change the filename in the primary header to reflect the name of the output
    # filename.
    fitsobj[0].header['FILENAME'] = str(output) #+"_drz.fits"

    # Change the ROOTNAME keyword to the ROOTNAME of the output PRODUCT
    fitsobj[0].header['ROOTNAME'] = str(output.split('_drz.fits')[0])
    # Modify the ASN_MTYP keyword to contain "PROD-DTH" so it can be properly
    # ingested into the archive catalog.

    #stis has this keyword in the [1] header, so I am directing the code
    #to first look in the primary, then the 1
    try:
        fitsobj[0].header['ASN_MTYP'] = 'PROD-DTH'
    except:
        fitsobj[1].header['ASN_MTYP'] = 'PROD-DTH'


    errstr =  "###############################################################################\n"
    errstr += "#                                                                             #\n"
    errstr += "# ERROR:                                                                      #\n"
    errstr += "# Multidrizzle has created an empty DRZ product because all input images were #\n"
    errstr += "# excluded from processing or a user requested the program to stop.           #\n"
    errstr += "#                                                                             #\n"
    errstr += "###############################################################################\n\n"
    print errstr

    # If the file is already on disk delete it and replace it with the
    # new file
    dirfiles = os.listdir(os.curdir)
    if (dirfiles.count(output) > 0):
        os.remove(output)
        print "       Replacing "+output+"..."

    # Write out the empty DRZ file
    fitsobj.writeto(output)
    return

def checkDGEOFile(filenames):
    """
    Verify that input file has been updated with NPOLFILE

    This function checks for the presence of 'NPOLFILE' kw in the primary header
    when 'DGEOFILE' kw is present and valid (i.e. 'DGEOFILE' is not blank or 'N/A').
    It handles the case of science files downloaded from the archive before the new
    software was installed there.
    If 'DGEOFILE' is present and 'NPOLFILE' is missing, print a message and let the user
    choose whether to (q)uit and update the headers or (c)ontinue and run betadrizzle
    without the non-polynomial correction.
    'NPOLFILE' will be populated in the pipeline before betadrizzle is run.

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

            Alternatively you may choose to run betadrizzle without DGEO and detector to image correction.

            To stop betadrizzle and update the dgeo files, type 'q'.
            To continue running betadrizzle without the non-polynomial distortion correction, type 'c':
            """
    for inputfile in filenames:
        if (pyfits.getval(inputfile,'INSTRUME') == 'WFPC2'):
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

def update_wfpc2_d2geofile(filename,fhdu=None):
    """ Creates a D2IMFILE from the DGEOFILE for a WFPC2 image (input),
        and modifies the header to reflect the new usage.
    """
    close_fhdu = False
    if fhdu is None:
        fhdu = fileutil.openImage(filename,mode='update')
        close_fhdu = True

    dgeofile = fhdu['PRIMARY'].header.get('DGEOFILE',None)
    if dgeofile not in [None, "N/A", "", " "]:
        print 'Converting DGEOFILE ',dgeofile,' into D2IMFILE...'
        rootname = filename[:filename.find('.fits')]
        d2imfile = convert_dgeo_to_d2im(dgeofile,rootname)
        fhdu['PRIMARY'].header.update('ODGEOFIL',dgeofile)
        fhdu['PRIMARY'].header.update('DGEOFILE','N/A')
        fhdu['PRIMARY'].header.update('D2IMFILE',d2imfile)
    else:
        d2imfile = None
        fhdu['PRIMARY'].header.update('DGEOFILE','N/A')
        fhdu['PRIMARY'].header.update('D2IMFILE','N/A')

    # Only close the file handle if opened in this function
    if close_fhdu:
        fhdu.close()

    # return the d2imfile name so that calling routine can keep
    # track of the new file created and delete it later if necessary
    # (multidrizzle clean=True mode of operation)
    return d2imfile

def convert_dgeo_to_d2im(dgeofile,output,clobber=True):
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
        'staticfile':None,
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
