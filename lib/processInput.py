from pytools import parseinput, fileutil, readgeis, asnutil,irafglob,check_files
import pyfits
import os,shutil

import wcs_functions,util
import mdzhandler

from updatewcs import hstwcs

"""
Process input to MultiDrizzle/PyDrizzle.
Input can be one of 

- a python list of files
- a comma separated string of filenames (including wild card characters)
- an association table
- an @file (can have a second column with names of ivm files)

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

    Syntax:
        imageObjectList,outwcs = processInput.processCommonInput(configObj)

        where,
        configObj: configObj instance or simple dictionary of input parameters        
        imageObjectList: list of imageObject instances, 1 for each input exposure
        outwcs: imageObject instance defining the final output frame

        you can set createOutwcs=False for the cases where you only want the
        images processed and no output wcs information in necessary
        
    At a minimum, the configObj instance (dictionary) should contain:
        configObj = {'input':None,'output':None,
                    'updatewcs':None,'shiftfile':None}

    If provided, the configObj should contain the values of all the multidrizzle parameters 
    as set by the user with TEAL. If no configObj is given, it will retrieve
    the default values automatically.  In either case, the values from the input_dict
    will be merged in with the configObj before being used by the rest of the 
    code. 

    """
    if not createOutwcs:
        configObj['updatewcs']=False #we're probably just working on single images here
        
        
    #maybe we can chunk this part up some more so that we can call just the parts we want
    
        
    # Interpret input, read and convert and update input files, then return
    # list of input filenames and derived output filename
    asndict,ivmlist,output = process_input(configObj['input'], configObj['output'], 
            updatewcs=configObj['updatewcs'], workinplace=configObj['workinplace'])

    # convert the filenames from asndict into a list of full filenames
    files = [fileutil.buildRootname(f) for f in asndict['order']]

    # interpret MDRIZTAB, if specified, and update configObj accordingly
    # This can be done here because MDRIZTAB does not include values for 
    # input, output, updatewcs, or shiftfile.
    if configObj['mdriztab']:
        mdriztab_dict = mdzhandler.getMdriztabParameters(files)
        # Update configObj with values from mpars
        configObj.update(mdriztab_dict)
        
    # Convert interpreted list of input files from process_input into a list
    # of imageObject instances for use by the MultiDrizzle tasks.
    instrpars = configObj['INSTRUMENT PARAMETERS']
    # pass in 'proc_unit' to initialize unit conversions as necessary
    instrpars['proc_unit'] = configObj['proc_unit']

    if configObj['shiftfile'] not in [None,""]:
        print '\n-Applying shiftfile ',configObj['shiftfile'],' to input images...\n'
        # Update all input images with shifts from shiftfile
        wcs_functions.createHeaderlets(configObj['shiftfile'])

    # Build imageObject list for all the valid, shift-updated input files
    print '\n-Creating imageObject List as input for processing steps.\n'
    imageObjectList = createImageObjectList(files,instrpars,group=configObj['group'])

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
    """ Add IVM filenames provided by user to outputNames dictionary for each
        input imageObject.
    """
    if ivmlist is None:
        return

    for img,ivmname in zip(imageObjectList,ivmlist):
        img.updateIVMName(ivmname)

def checkMultipleFiles(input):
    a,i,o = process_input(input,updatewcs=False)
    return len(a['members']) > 1

def createImageObjectList(files,instrpars,group=None):
    """ Returns a list of imageObject instances, 1 for each input image in the
        list of input filenames.
    """
    imageObjList = []
    for img in files:
        image = _getInputImage(img,group=group)
        image.setInstrumentParameters(instrpars)
        imageObjList.append(image)

    return imageObjList

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
def atfile_sci(f):
    return f.split()[0]
    
def atfile_ivm(f):
    return f.split()[1]    

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
        oldasndict = asnutil.readASNTable(input, prodonly=prodonly)
        
        if not infilesOnly:
            if not output:
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
        filelist = irafglob.irafglob(input, atfile=atfile_sci)
        # If there is a second column...
        if len(line.split()) == 2:
            # ...parse out the names of the IVM files as well 
            ivmlist = irafglob.irafglob(input, atfile=atfile_ivm)        
    else:
        #input is a string or a python list
        try:
            filelist, output = parseinput.parseinput(input, outputname=output)
            if output in ['',None]: 
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
    
def process_input(input, output=None, ivmlist=None, updatewcs=True, prodonly=False, workinplace=True):
    
    filelist,output,ivmlist,oldasndict=processFilenames(input,output)
    if not workinplace:
        createInputCopies(filelist)
    newfilelist, ivmlist = check_files.checkFiles(filelist, ivmlist)
    
    if not newfilelist:
        buildEmptyDRZ(input,output)
        return None, None, output 
    
    #make an asn table at the end
    if updatewcs:
        pydr_input = runmakewcs(newfilelist)  
    else:
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

def runmakewcs(input):
    """
    Runs make wcs and recomputes the WCS keywords
    input: a list of files
    output: returns a list of names of the modified files
            (For GEIS files returns the translated names.)
    """
    newNames = hstwcs.updatewcs(input,checkfiles=False)
    #newNames = makewcs.run(input)
    return newNames

def update_member_names(oldasndict, pydr_input):
    """
    Purpose
    =======
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


def createInputCopies(filelist):
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
    
    # check to see if copies already exist for each file
    for fname in filelist:
        copyname = os.path.join(origdir,fname)
        if not os.path.exists(copyname):
            print 'Preserving original of: ',fname, 'as ',copyname
            # make a copy of the file in the sub-directory
            shutil.copy(fname,copyname)
            os.chmod(copyname,0444) # make files read-only
        else:
            print 'Restoring original input for ',fname,' from ',copyname
            # replace current files with original version
            os.chmod(fname,0666)
            shutil.copy(copyname,fname)
            os.chmod(fname,0666)

def buildEmptyDRZ(input, output):
    """
    
    FUNCTION: buildEmptyDRZ
    PURPOSE : Create an empty DRZ file in a valid FITS format so that the HST
              pipeline can handle the Multidrizzle zero expossure time exception
              where all data has been excluded from processing.
    INPUT   : None
    OUTPUT  : DRZ file on disk
     
    """
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
            
    # Open the first image of the excludedFileList to use as a template to build
    # the DRZ file.
    inputfile = parseinput.parseinput(input)[0]
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
    fitsobj[0].header.add_history("** Multidrizzle has created this empty DRZ **")
    fitsobj[0].header.add_history("** product because all input images were   **")
    fitsobj[0].header.add_history("** excluded from processing because their  **")
    fitsobj[0].header.add_history("** header EXPTIME values were 0.0.  If you **")
    fitsobj[0].header.add_history("** still wish to use this data make the    **")
    fitsobj[0].header.add_history("** EXPTIME values in the header non-zero.  **")
    
    # Change the filename in the primary header to reflect the name of the output
    # filename.
    fitsobj[0].header['FILENAME'] = str(output) #+"_drz.fits"
    
    # Change the ROOTNAME keyword to the ROOTNAME of the output PRODUCT
    fitsobj[0].header['ROOTNAME'] = str(output.split('_drz.fits')[0])
    print 'self.output', output
    # Modify the ASN_MTYP keyword to contain "PROD-DTH" so it can be properly
    # ingested into the archive catalog.
    
    #stis has this keyword in the [1] header, so I am directing the code
    #to first look in the primary, then the 1
    try:
        fitsobj[0].header['ASN_MTYP'] = 'PROD-DTH'
    except:
        fitsobj[1].header['ASN_MTYP'] = 'PROD-DTH'
 
        
    errstr =  "#############################################\n"
    errstr += "#                                           #\n"
    errstr += "# ERROR:                                    #\n"
    errstr += "#  Multidrizzle has created this empty DRZ  #\n"
    errstr += "#  product because all input images were    #\n"
    errstr += "#  excluded from processing because their   #\n"
    errstr += "#  header EXPTIME values were 0.0.  If you  #\n"
    errstr += "#  still wish to use this data make the     #\n"
    errstr += "#  EXPTIME values in the header non-zero.   #\n"
    errstr += "#                                           #\n"
    errstr += "#############################################\n\n"
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

def _setDefaults(input_dict={}):
    """ Define full set of default values for unit-testing this module."""
    paramDict = {
        'input':'*flt.fits',
        'output':None,
        'mdriztab':None,
        'refimage':None,
        'runfile':None,
        'workinplace':False,
        'updatewcs':True,
        'proc_unit':'native',
        'coeffs':'header',
        'context':False,
        'clean':True,
        'group':None,
        'ra':None,
        'dec':None,
        'build':True,
        'shiftfile':None,
        'gain':None,
        'gnkeyword':None,
        'readnoise':None,
        'rnkeyword':None,
        'exptime':None,
        'expkeyword':None,
        'crbitval':4096,
        'shiftfile':None,
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
