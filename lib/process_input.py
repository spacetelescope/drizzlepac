from pytools import parseinput, fileutil, readgeis, makewcs, asnutil,irafglob
import pyfits
import os 

"""
Process input to pydrizzle.
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
Returns an association table, ivmlist, output name

"""

def atfile_sci(f):
    return f.split()[0]
    
def atfile_ivm(f):
    return f.split()[1]    


def process_input(input, output=None, ivmlist=None, updatewcs=True, prodonly=False, shiftfile=None):
    
    ivmlist = None
    oldasndict = None
    
    if (isinstance(input, list) == False) and \
       ('_asn' in input or '_asc' in input) :
        # Input is an association table
        # Get the input files, and run makewcs on them
        oldasndict = asnutil.readASNTable(input, prodonly=prodonly)
        if not output:
            output = oldasndict['output']

        filelist = [fileutil.buildRootname(fname) for fname in oldasndict['order']]
        
    elif (isinstance(input, list) == False) and \
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
            #filelist.sort()
        except IOError: raise
    
    # sort the list of input files
    # this ensures the list of input files has the same order on all platforms
    # it can have ifferent order because listdir() uses inode order, not unix type order 
    filelist.sort()
    newfilelist, ivmlist = checkFiles(filelist, ivmlist)
    
   
    if not newfilelist:
        buildEmptyDRZ(input,output)
        return None, None, output 
    
    #make an asn table at the end
    if updatewcs:
        pydr_input = runmakewcs(newfilelist)  
    else:
        pydr_input = newfilelist

    # AsnTable will handle the case when output==None
    if not oldasndict:        
        oldasndict = asnutil.ASNTable(pydr_input, output=output)
        oldasndict.create()
                
    if shiftfile:
        oldasndict.update(shiftfile=shiftfile)

    asndict = update_member_names(oldasndict, pydr_input)
    
    # Build output filename
    if output == None:
        output = fileutil.buildNewRootname(asndict['output'],extn='_drz.fits')
    else:
        if 'drz' not in output:
            output = fileutil.buildNewRootname(output,extn='_drz.fits')
        
    print 'Setting up output name: ',output

    return asndict, ivmlist, output


def checkFiles(filelist, ivmlist = None):
    """
    1. Converts waiver fits sciece and data quality files to MEF format
    2. Converts all GEIS science and data quality files to MEF format
    3. Checks for stis association tables 
    4. Checks if kw idctab exists, if not tries to populate it 
        based on the spt file
    5. Removes files with EXPTIME=0 and the corresponding ivm files
    6. Removes files with NGOODPIX == 0 (to exclude saturated images)
    """
    #newfilelist = []
    removed_files = []
    translated_names = []
    newivmlist = []
    
    if ivmlist == None:
        ivmlist = [None for l in filelist]

    sci_ivm = zip(filelist, ivmlist)
    
    for file in sci_ivm:
        #find out what the input is
        # if science file is not found on disk, add it to removed_files for removal
        try:
            imgfits,imgtype = fileutil.isFits(file[0])
        except IOError:
            print "Warning:  File %s could not be found\n" %file[0]
            print "Removing file %s from input list" %file[0]
            removed_files.append(file)
            continue
        if file[1] != None:
            #If an ivm file is not found on disk
            # Remove the corresponding science file
            try:
                ivmfits,ivmtype = fileutil.isFits(file[1])
            except IOError:
                print "Warning:  File %s could not be found\n" %file[1]
                print "Removing file %s from input list" %file[0]
                removed_files.append(file)
        # Check for existence of waiver FITS input, and quit if found.
        # Or should we print a warning and continue but not use that file
        if imgfits and imgtype == 'waiver':
            newfilename = waiver2mef(file[0], convert_dq=True)
            if newfilename == None:
                print "Removing file %s from input list - could not convert waiver to mef" %file[0]
                removed_files.append(file[0])
            else:
                translated_names.append(newfilename)

        # If a GEIS image is provided as input, create a new MEF file with 
        # a name generated using 'buildFITSName()'
        # Convert the corresponding data quality file if present    
        if not imgfits:
            newfilename = geis2mef(file[0], convert_dq=True)
            if newfilename == None:
                print "Removing file %s from input list - could not convert geis to mef" %file[0]
                removed_files.append(file[0])
            else:
                translated_names.append(newfilename)
        if file[1] != None:
            if ivmfits and ivmtype == 'waiver':
                print "Warning: PyDrizzle does not support waiver fits format.\n"
                print "Convert the input files to GEIS or multiextension FITS.\n"
                print "File %s appears to be in waiver fits format \n" %file[1]
                print "Removing file %s from input list" %file[0] 
                removed_files.append(file[0])
  
            if not ivmfits:
                newfilename = geis2mef(file[1], convert_dq=False)
                if newfilename == None:
                    print "Removing file %s from input list" %file[0]
                    removed_files.append(file[0])
                else:
                    newivmlist.append(newfilename)

    newfilelist, ivmlist = update_input(filelist, ivmlist, removed_files)

    if newfilelist == []:
        return [], []
    """
    errormsg = "\n No valid input was found. Quitting ...\n"
    raise IOError, errormsg
    """
    if translated_names != []:
        # Since we don't allow input from different instruments
        # we can abandon the original input list and provide as 
        # input only the translated names
        removed_expt_files = check_exptime(translated_names)
        newfilelist, ivmlist = update_input(translated_names, newivmlist, removed_expt_files)
    else:
        # check for STIS association files. This must be done before 
        # the check for EXPTIME in order to handle correctly stis 
        # assoc files
        if pyfits.getval(newfilelist[0], 'INSTRUME') == 'STIS':
            newfilelist, ivmlist = checkStisFiles(newfilelist, ivmlist)
            #removed_files = check_exptime(newflist)
        
        removed_expt_files = check_exptime(newfilelist)
        newfilelist, ivmlist = update_input(newfilelist, ivmlist, removed_expt_files)
    if removed_expt_files:
        errorstr =  "#############################################\n"
        errorstr += "#                                           #\n"
        errorstr += "# ERROR:                                    #\n"
        errorstr += "#                                           #\n"
        errorstr += "#  The following files were excluded from   #\n"
        errorstr += "#  Multidrizzle processing because their    #\n"
        errorstr += "#  header keyword EXPTIME values were 0.0:  #\n"
        for name in removed_expt_files:
            errorstr += "         "+ str(name) + "\n" 
        errorstr += "#                                           #\n"
        errorstr += "#############################################\n\n"
        print errorstr
        
    removed_ngood_files = checkNGOODPIX(newfilelist)
    newfilelist, ivmlist = update_input(newfilelist, ivmlist, removed_ngood_files)
    if removed_ngood_files:
        msgstr =  "####################################\n"
        msgstr += "#                                  #\n"
        msgstr += "# WARNING:                         #\n"
        msgstr += "#  NGOODPIX keyword value of 0 in  #\n"
        for name in removed_ngood_files:
            msgstr += "         "+ str(name) + "\n" 
        msgstr += "#  has been detected.  Images with #\n"
        msgstr += "#  no valid pixels will not be     #\n"
        msgstr += "#  used during processing.  If you #\n"
        msgstr += "#  wish this file to be used in    #\n"
        msgstr += "#  processing, please check its DQ #\n"
        msgstr += "#  array and reset driz_sep_bits   #\n"
        msgstr += "#  and final_bits parameters       #\n"
        msgstr += "#  to accept flagged pixels.       #\n"
        msgstr += "#                                  #\n"
        msgstr += "####################################\n"
        print msgstr   
                
    return newfilelist, ivmlist
    
def waiver2mef(sciname, newname=None, convert_dq=True):
    """
    Converts a GEIS science file and its corresponding 
    data quality file (if present) to MEF format
    Writes out both files to disk.
    Returns the new name of the science image.
    """
    
    def convert(file):
        newfilename = fileutil.buildNewRootname(file, extn='_c0h.fits')
        try:
            newimage = fileutil.openImage(file,writefits=True,
                                          fitsname=newfilename,clobber=True)
            del newimage
            return newfilename
        except IOError:
            print 'Warning: File %s could not be found' % file     
            return None
        
    newsciname = convert(sciname)
    if convert_dq:
        dq_name = convert(fileutil.buildNewRootname(sciname, extn='_c1h.fits'))
        
    return newsciname   



def geis2mef(sciname, convert_dq=True):
    """
    Converts a GEIS science file and its corresponding 
    data quality file (if present) to MEF format
    Writes out both files to disk.
    Returns the new name of the science image.
    """
        
    def convert(file):
        newfilename = fileutil.buildFITSName(file)
        try:
            newimage = fileutil.openImage(file,writefits=True,
                fitsname=newfilename, clobber=True)            
            del newimage
            return newfilename
        except IOError:
            print 'Warning: File %s could not be found' % file     
            return None

    newsciname = convert(sciname)
    if convert_dq:
        dq_name = convert(sciname.split('.')[0] + '.c1h')
        
    return newsciname


def checkStisFiles(filelist, ivmlist=None):
    newflist = []
    newilist = []
    
    if len(filelist) != len(ivmlist):
        errormsg = "Input file list and ivm list have different lenghts\n"
        errormsg += "Quitting ...\n"
        raise ValueError, errormsg
        
    for t in zip(filelist, ivmlist):
        sci_count = stisObsCount(t[0])
        if sci_count >1:
            newfilenames = splitStis(t[0], sci_count)
            newflist.extend(newfilenames)
            if t[1] != None:
                newivmnames = splitStis(t[1], sci_count)
                newilist.extend(newivmnames)
            else:
                newilist.append(None)
        elif sci_count == 1:
            newflist.append(t[0])
            newilist.append(t[1])
        else:
            errormesg = "No valid 'SCI extension in STIS file\n"
            raise ValueError, errormsg

    return newflist, newilist


def runmakewcs(input):
    """
    Runs make wcs and recomputes the WCS keywords
    input: a list of files
    output: returns a list of names of the modified files
            (For GEIS files returns the translated names.)
    """
    newNames = makewcs.run(input)
    return newNames

def check_exptime(filelist):
    """
    Removes files with EXPTIME==0 from filelist.
    """
    removed_files = []
    
    for f in filelist:
        if fileutil.getKeyword(f, 'EXPTIME') <= 0: 
            removed_files.append(f)
            
    return removed_files

def checkNGOODPIX(filelist):
    """
    Only for ACS, and STIS, check NGOODPIX
    If all pixels are 'bad' on all chips, exclude this image
    from further processing. 
    Similar checks requiring comparing 'driz_sep_bits' against
    WFPC2 c1f.fits arrays and NICMOS DQ arrays will need to be
    done separately (and later).
    """
    removed_files = []
    for inputfile in filelist:
        if (fileutil.getKeyword(inputfile,'instrume') == 'ACS') \
           or fileutil.getKeyword(inputfile,'instrume') == 'STIS': 
            _file = fileutil.openImage(inputfile)
            _ngood = 0
            for extn in _file:
                if extn.header.has_key('EXTNAME') and extn.header['EXTNAME'] == 'SCI':
                    _ngood += extn.header['NGOODPIX']
            _file.close()
            
            if (_ngood == 0):
                removed_files.append(inputfile)
    return removed_files

def update_input(filelist, ivmlist=None, removed_files=None):
    """
    Removes files flagged to be removed from the input filelist.
    Removes the corresponding ivm files if present.
    """
    newfilelist = []

    if removed_files == []:
        return filelist, ivmlist
    else:
        sci_ivm = zip(filelist, ivmlist)
        for f in removed_files:
            result=[sci_ivm.remove(t) for t in sci_ivm if t[0] == f ]
        ivmlist = [el[1] for el in sci_ivm] 
        newfilelist = [el[0] for el in sci_ivm] 
        return newfilelist, ivmlist 
  

def stisObsCount(input):
    """
    Input: A stis multiextension file
    Output: Number of stis science extensions in input
    """
    count = 0
    f = pyfits.open(input)
    for ext in f:
        if ext.header.has_key('extname'):
            if (ext.header['extname'].upper() == 'SCI'):
                count += 1
    f.close()
    return count

def splitStis(stisfile, sci_count):
    """
    Purpose
    =======
    
    Split a STIS association file into multiple imset MEF files.
    Split the corresponding spt file if present into single spt files.
    If an spt file can't be split or is missing a Warning is printed.
    
    Output: a list with the names of the new flt files.
    """
    newfiles = []
    
    f = pyfits.open(stisfile)
    hdu0 = f[0].copy()


    for count in range(1,sci_count+1):
        #newfilename = rootname+str(count)+'.fits'
        fitsobj = pyfits.HDUList()            
        fitsobj.append(hdu0)
        hdu = f['sci',count].copy()
        fitsobj.append(hdu)
        rootname = hdu.header['EXPNAME']
        newfilename = fileutil.buildNewRootname(rootname, extn='_flt.fits')
        try:
            # Verify error array exists
            if f['err',count].data == None:
                raise ValueError
            # Verify dq array exists
            if f['dq',count].data == None:
                raise ValueError
            # Copy the err extension
            hdu = f['err',count].copy()
            fitsobj.append(hdu)
            # Copy the dq extension
            hdu = f['dq',count].copy()
            fitsobj.append(hdu)
        except:
            errorstr =  "\n###############################\n"
            errorstr += "#                             #\n"
            errorstr += "# ERROR:                      #\n"
            errorstr += "#  The input image:           #\n"
            errorstr += "      " + str(stisfile) +"\n"
            errorstr += "#  does not contain required  #\n"
            errorstr += "#  image extensions.  Each    #\n"
            errorstr += "#  must contain populated sci,#\n"
            errorstr += "#  dq, and err arrays.        #\n"
            errorstr += "#                             #\n"
            errorstr += "###############################\n"
            raise ValueError, errorstr
        
        
        # Update the 'EXTNER' keyword to indicate the new extnesion number
        # for the single exposure files.
        fitsobj[1].header['EXTVER'] = 1
        fitsobj[2].header['EXTVER'] = 1
        fitsobj[3].header['EXTVER'] = 1
        
        # Determine if the file you wish to create already exists on the disk.
        # If the file does exist, replace it.
        if (os.path.exists(newfilename)):
            os.remove(newfilename)
            print "       Replacing "+newfilename+"..."
            
            # Write out the new file
        fitsobj.writeto(newfilename)
        newfiles.append(newfilename)
    f.close()

    sptfilename = fileutil.buildNewRootname(stisfile, extn='_spt.fits')
    try:
        sptfile = pyfits.open(sptfilename)
    except IOError:
        print 'SPT file not found %s \n' % sptfilename

    if sptfile:
        hdu0 = sptfile[0].copy()
        try:
            for count in range(1,sci_count+1):
                fitsobj = pyfits.HDUList()            
                fitsobj.append(hdu0)
                hdu = sptfile[count].copy()
                fitsobj.append(hdu)
                rootname = hdu.header['EXPNAME']
                newfilename = fileutil.buildNewRootname(rootname, extn='_spt.fits')
                fitsobj[1].header['EXTVER'] = 1
                if (os.path.exists(newfilename)):
                    os.remove(newfilename)
                    print "       Replacing "+newfilename+"..."
            
                # Write out the new file
                fitsobj.writeto(newfilename)
        except:
            print "Warning: Unable to split spt file %s " % sptfilename
        sptfile.close()
    
    return newfiles 


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



def buildEmptyDRZ(input, output):
    """
    
    METHOD  : _buildEmptyDRZ
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
    fitsobj[0].header['ASN_MTYP'] = 'PROD-DTH'
    
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

