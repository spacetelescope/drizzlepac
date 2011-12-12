#!/usr/bin/env python
"""
A library of utility functions

"""
from __future__ import division # confidence medium
import os
import logging
import traceback
import sys
import string

import numpy as np
import pyfits
from stsci.tools import asnutil, fileutil, teal, cfgpars

__version__ = "0.1.0tng1"
__pyfits_version__ = pyfits.__version__
__numpy_version__ = np.__version__


multiprocessing = None
can_parallel = False
_cpu_count = -1
if 'ASTRODRIZ_NO_PARALLEL' not in os.environ:
    try:
        import multiprocessing
        try:
            # sanity check - do we even have the hardware?
            _cpu_count = multiprocessing.cpu_count()
            can_parallel = _cpu_count > 1
        except:
            can_parallel = False
    except ImportError:
        print '\nCould not import multiprocessing, will only take advantage of a single CPU core'

def get_pool_size(usr_config_value = None):
    """ Use the suggested pool size (from cpu_count) and the config
    object value, to get the right pool size to use. Consolidate all
    such logic here, not in the caller. """
    if not can_parallel:
        return 0
    if usr_config_value != None:
        return usr_config_value
    return _cpu_count


DEFAULT_LOGNAME = 'astrodrizzle.log'
blank_list = [None, '', ' ',"None","INDEF"]

def is_blank(val):
    """ Determines whether or not a value is considered 'blank'.
    """
    blank = False
    if val in blank_list: blank = True
    return blank

def check_blank(cvar):
    """ Converts blank value (from configObj?) into a value of None.
    """
    if cvar in blank_list: val = None
    else: val = cvar
    return val

"""
Logging routines
"""
class StreamLogger(object):
    """ Class to manage trapping of STDOUT and STDERR messages to a trailer file
    """

    def __init__(self, stream, logfile, mode='w', prefix=''):
        self.stream = stream
        if is_blank(prefix):
            self.prefix = ''
        else:
            self.prefix = '['+prefix+'] '
        self.data = ''

        # set up logfile
        if not is_blank(logfile):
            self.log = open(logfile,mode)
            self.filename = logfile
            # clear out any previous exceptions, so that only those generated
            # by this code will be picked up in the trailer file
            sys.exc_clear()
            print '[astrodrizzle] Trailer file will be written out to: ',self.filename
        else:
            self.log = None
            self.filename = None
            print '[astrodrizzle] No trailer file will be created...'

    def write(self, data):
        self.stream.write(data)
        self.stream.flush()

        if self.log is not None:
            self.data += data
            tmp = str(self.data)
            if '\x0a' in tmp or '\x0d' in tmp:
                tmp = tmp.rstrip('\x0a\x0d')
                self.log.write('%s%s\n' % (self.prefix,tmp))
                self.data = ''
    def flush(self):
        self.stream.flush()

def init_logging(logfile=DEFAULT_LOGNAME,default=None):
    """ Set up logfile for capturing stdout/stderr messages.
        Must be called prior to writing any messages that you want to log.
    """
    if logfile == "INDEF":
        if not is_blank(default):
            logname = fileutil.buildNewRootname(default) +'.log'
        else:
            logname = DEFAULT_LOGNAME
    elif logfile not in [None,""," "]:
        if '.log' in logfile:
            logname = logfile
        else:
            logname = logfile+'.log'
    else:
        logname = None
    if logname is not None:
        # redirect logging of stdout to logfile
        sys.stdout = StreamLogger(sys.stdout, logname)
    else:
        print '[astrodrizzle] No trailer file created...'

def end_logging():
    """ Close log file and restore stdout/stderr to system defaults.
    """
    if hasattr(sys.stdout,'log'): # only need to close the log if one was created
        if sys.stdout.log is not None:
            print '[astrodrizzle] Trailer file written to: ',sys.stdout.filename
            sys.stdout.log.flush()
            sys.stdout.log.close()

            # Add any traceback information to the trailer file to document
            # the error that caused the code to stop processing
            if sys.exc_info()[0] is not None:
                errfile = open(sys.stdout.filename,mode="a")
                traceback.print_exc(None,errfile)
                errfile.close()
        else:
            print '[astrodrizzle] No trailer file saved...'

        sys.stdout = sys.__stdout__

def print_pkg_versions(packages=None,svn=False):
    pkgs = ['numpy','pyfits','stwcs','pywcs']
    if packages is not None:
        if not isinstance(packages,list):
            packages = [packages]
        pkgs.extend(packages)
        
    print 'Version Information'
    print '-'*20
    print 'Python Version '+sys.version
    for software in pkgs:
        try:
            package = __import__(software)
            vstr = "%s "%(software)
            try:
                vstr+= "Version -> "+package.__version__+" "
            except:
                vstr += "No version defined.  "
            if svn:
                try:
                    vstr += "\n    SVN version -> "+package.__svn_version__.rstrip()
                except:
                    vstr += " "
        except:
            vstr = software+" not found in path..."
        print vstr

class ProcSteps:
    """ This class allows MultiDrizzle to keep track of the
        start and end times of each processing step that gets run
        as well as computing/reporting the elapsed time for each step.

        The code for each processing step must call the 'addStep()'
        method to initialize the information for that step, then
        the 'endStep()' method to record the end and elapsed times.

        The 'reportTimes()' method can then be used to provide a summary
        of all the elapsed times and total run time.
    """
    __report_header = '\n   %20s          %s\n'%('-'*20,'-'*20)
    __report_header += '   %20s          %s\n'%('Step','Elapsed time')
    __report_header += '   %20s          %s\n'%('-'*20,'-'*20)

    def __init__(self):
        self.steps = {}
        self.order = []
        self.start = _ptime()
        self.end = None

    def addStep(self,key):
        """
        Add information about a new step to the dict of steps
        The value 'ptime' is the output from '_ptime()' containing
        both the formatted and unformatted time for the start of the
        step.
        """
        ptime = _ptime()
        print '==== Processing Step ',key,' started at ',ptime[0]
        self.steps[key] = {'start':ptime}
        self.order.append(key)

    def endStep(self,key):
        """
        Record the end time for the step.

        If key==None, simply record ptime as end time for class to represent
        the overall runtime since the initialization of the class.
        """
        ptime = _ptime()
        if key is not None:
            self.steps[key]['end'] = ptime
            self.steps[key]['elapsed'] = ptime[1] - self.steps[key]['start'][1]
        self.end = ptime

        print'==== Processing Step ',key,' finished at ',ptime[0]

    def reportTimes(self):
        """ Print out a formatted summary of the elapsed times for all
            the performed steps.
        """
        self.end = _ptime()
        total_time = 0
        print ProcSteps.__report_header
        for step in self.order:
            if self.steps[step].has_key('elapsed'):
                _time = self.steps[step]['elapsed']
            else:
                _time = 0.0
            total_time += _time
            print '   %20s          %0.4f sec.'%(step,_time)

        print '   %20s          %s'%('='*20,'='*20)
        print '   %20s          %0.4f sec.'%('Total',total_time)

        # Compute overall runtime of entire program, including overhead
        #total = self.end[1] - self.start[1]
        #print '   %20s          %0.4f sec.'%('Total Runtime',total)

def _ptime():
    import time
    try:
        import datetime as dtime
    except ImportError:
        dtime = None
    ftime = time.time()
    if dtime:
        # This time stamp includes sub-second timing...
        _ltime = dtime.datetime.fromtimestamp(ftime)
        tlm_str = _ltime.strftime("%H:%M:%S")+str(_ltime.microsecond/1e+6)[1:-3]+_ltime.strftime(" (%d/%m/%Y)")
    else:
        # Basic time stamp which only includes integer seconds
        # Format time values for keywords IRAF-TLM, and DATE
        _ltime = time.localtime(ftime)
        tlm_str = time.strftime('%H:%M:%S (%d/%m/%Y)',_ltime)
        #date_str = time.strftime('%Y-%m-%dT%H:%M:%S',_ltime)
    return tlm_str,ftime

def findrootname(filename):
    """
    Return the rootname of the given file.
    """

    puncloc = [filename.find(char) for char in string.punctuation]
    val = sys.maxint
    for num in puncloc:
        if num !=-1 and num < val:
            val = num
    return filename[0:val]

def removeFileSafely(filename,clobber=True):
    """ Delete the file specified, but only if it exists and clobber is True.
    """
    if filename is not None and filename.strip() != '':
        if os.path.exists(filename) and clobber: os.remove(filename)

def getDefaultConfigObj(taskname,configObj,input_dict={},loadOnly=True):
    """ Return default configObj instance for task updated
        with user-specified values from input_dict.

        Parameters
        ----------
        taskname : string
            Name of task to load into TEAL

        configObj : string
            The valid values for 'configObj' would be::

                None                      - loads last saved user .cfg file
                'defaults'                - loads task default .cfg file
                name of .cfg file (string)- loads user-specified .cfg file

        input_dict : dict
            Set of parameters and values specified by user to be different from
            what gets loaded in from the .cfg file for the task

        loadOnly : bool
            Setting 'loadOnly' to False causes the TEAL GUI to start allowing the
            user to edit the values further and then run the task if desired.

    """
    if configObj is None:
        # Start by grabbing the default values without using the GUI
        # This insures that all subsequent use of the configObj includes
        # all parameters and their last saved values
        configObj = teal.load(taskname)
    elif isinstance(configObj,str):
        if configObj.lower().strip() == 'defaults':
            # Load task default .cfg file with all default values
            configObj = teal.load(taskname,defaults=True)
            # define default filename for configObj
            configObj.filename = taskname.lower()+'.cfg'
        else:
            # Load user-specified .cfg file with its special default values
            # we need to call 'fileutil.osfn()' to insure all environment
            # variables specified by the user in the configObj filename are
            # expanded to the full path
            configObj = teal.load(fileutil.osfn(configObj))

    # merge in the user values for this run
    # this, though, does not save the results for use later
    if input_dict not in [None,{}]:# and configObj not in [None, {}]:
        cfgpars.mergeConfigObj(configObj, input_dict)
        # Update the input .cfg file with the updated parameter values
        #configObj.filename = os.path.join(cfgpars.getAppDir(),os.path.basename(configObj.filename))
        #configObj.write()

    if not loadOnly:
    # We want to run the GUI AFTER merging in any parameters
    # specified by the user on the command-line and provided in
    # input_dict
        configObj = teal.teal(configObj,loadOnly=False)

    return configObj

def getSectionName(configObj,stepnum):
    """ Return section label based on step number.
    """
    for key in configObj.keys():
        if key.find('STEP '+str(stepnum)) >= 0:
            return key

def getConfigObjPar(configObj, parname):
    """ Return parameter value without having to specify which section
        holds the parameter.
    """
    return cfgpars.findFirstPar(configObj, parname)[1]

def displayMakewcsWarningBox(display=True, parent=None):
    """ Displays a warning box for the 'makewcs' parameter.
    """
    import tkMessageBox

    ans = {'yes':True,'no':False}
    if ans[display]:
        msg = 'Setting "updatewcs=yes" will result '+ \
              'in all input WCS values to be recomputed '+ \
              'using the original distortion model and alignment.'
        tkMessageBox.showwarning(parent=parent,message=msg, title="WCS will be overwritten!")
    return True

"""
These two functions are for reading in an 'at file' which contains
two columns of filenames, the first column is assumed to
be the science image and the second column is assumed
to be the IVM file that is associated with it.
"""

def atfile_sci(filename):
    """
    Return the filename of the science image
    which is assumed to be the first word
    in the atfile the user gave.
    """
    return filename.split()[0]


def atfile_ivm(filename):
    """
    Return the filename of the IVM file
    which is assumed to be the second word
    in the atfile the user gave.
    """
    return filename.split()[1]


def printParams(paramDictionary,all=False):
    """ Print nicely the parameters from the dictionary.
    """

    if (len(paramDictionary) == 0):
        print "\nNo parameters were supplied\n"
    else:
        keys=paramDictionary.keys()
        keys.sort()
        for key in keys:
            if all or (not isinstance(paramDictionary[key],dict)):
                print "\t",key,":\t",paramDictionary[key]
        print '\n'
    sys.stdout.flush()

def isASNTable(inputFilelist):
    """Return TRUE if inputFilelist is a fits ASN file."""
    if ("_asn"  or "_asc") in inputFilelist:
        return True
    return False

def isCommaList(inputFilelist):
    """Return True if the input is a comma separated list of names."""
    if "," in inputFilelist:
        return True
    return False

def loadFileList(inputFilelist):
    """Open up the '@ file' and read in the science and possible
       ivm filenames from the first two columns.
    """
    f = open(inputFilelist[1:])
    # check the first line in order to determine whether
    # IVM files have been specified in a second column...
    lines = f.readline()
    f.close()

    # If there is a second column...
    if len(line.split()) == 2:
        # ...parse out the names of the IVM files as well
        ivmlist = irafglob.irafglob(input, atfile=atfile_ivm)

    # Parse the @-file with irafglob to extract the input filename
    filelist = irafglob.irafglob(input, atfile=atfile_sci)
    return filelist


def readCommaList(fileList):
    """ Return a list of the files with the commas removed. """
    names=fileList.split(',')
    fileList=[]
    for item in names:
        fileList.append(item)
    return fileList

def runmakewcs(input):
    """
    Runs 'updatewcs' to recompute the WCS keywords for the input image.

    Parameters
    ----------
    input : list of str
        A list of file names.

    Returns
    -------
    output : list of str
        Returns a list of names of the modified files
        (For GEIS files returns the translated names).

    """
    newNames = updatewcs.updatewcs(input, checkfiles=False)

    return newNames


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


####
#
# The following functions were required for use with the drizzling code
# and were copied in from pydrizzle_tng.py.
#
####

def countImages(imageObjectList):
    expnames = []
    for img in imageObjectList:
        expnames += img.getKeywordList('_expname')
    imgnames = []

    nimages = 0
    for e in expnames:
        if e not in imgnames:
            imgnames.append(e)
            nimages += 1
    return nimages

def get_detnum(hstwcs,filename,extnum):
    detnum = None
    binned = None
    if hstwcs.filename == filename and hstwcs.extver == extnum:
        detnum = hstwcs.chip
        binned = hstwcs.binned

    return detnum,binned

def get_expstart(header,primary_hdr):
    """shouldn't this just be defined in the instrument subclass of imageobject?"""

    if primary_hdr.has_key('expstart'):
        exphdr = primary_hdr
    else:
        exphdr = header

    if exphdr.has_key('EXPSTART'):
        expstart = float(exphdr['EXPSTART'])
        expend = float(exphdr['EXPEND'])
    else:
        expstart = 0.
        expend = 0.0

    return (expstart,expend)

def compute_texptime(imageObjectList):
    """
    Add up the exposure time for all the members in
    the pattern, since 'drizzle' doesn't have the necessary
    information to correctly set this itself.
    """
    expnames = []
    exptimes = []
    start = []
    end = []
    for img in imageObjectList:
        expnames += img.getKeywordList('_expname')
        exptimes += img.getKeywordList('_exptime')
        start += img.getKeywordList('_expstart')
        end += img.getKeywordList('_expend')

    exptime = 0.
    expstart = min(start)
    expend = max(end)
    exposure = None
    for n in range(len(expnames)):
        if expnames[n] != exposure:
            exposure = expnames[n]
            exptime += exptimes[n]

    return (exptime,expstart,expend)

def computeRange(corners):
    """ Determine the range spanned by an array of pixel positions. """
    _xrange = (np.minimum.reduce(corners[:,0]),np.maximum.reduce(corners[:,0]))
    _yrange = (np.minimum.reduce(corners[:,1]),np.maximum.reduce(corners[:,1]))
    return _xrange,_yrange

def getRotatedSize(corners,angle):
    """ Determine the size of a rotated (meta)image."""
    # If there is no rotation, simply return original values
    if angle == 0.:
        _corners = corners
    else:
        # Find center
        #_xr,_yr = computeRange(corners)
        #_cen = ( ((_xr[1] - _xr[0])/2.)+_xr[0],((_yr[1]-_yr[0])/2.)+_yr[0])
        _rotm = fileutil.buildRotMatrix(angle)
        # Rotate about the center
        #_corners = N.dot(corners - _cen,_rotm)
        _corners = np.dot(corners,_rotm)

    return computeRange(_corners)

def readcols(infile,cols=[0,1,2,3],hms=False):
    """
    Read the columns from an ASCII file as numpy arrays.

    Parameters
    ----------
    infile : str
        Filename of ASCII file with array data as columns.

    cols : list of int
        List of 0-indexed column numbers for columns to be turned into numpy arrays
        (DEFAULT- [0,1,2,3]).

    Returns
    -------
    outarr : list of numpy arrays
        Simple list of numpy arrays in the order as specifed in the 'cols' parameter.

    """

    fin = open(infile,'r')
    outarr = []
    for l in fin.readlines():
        l = l.strip()
        if len(l) == 0 or len(l.split()) < len(cols) or (len(l) > 0 and l[0] == '#' or (l.find("INDEF") > -1)): continue
        for i in range(10):
            lnew = l.replace("  "," ")
            if lnew == l: break
            else: l = lnew
        lspl = lnew.split(" ")

        if len(outarr) == 0:
            for c in range(len(cols)): outarr.append([])

        for c,n in zip(cols,range(len(cols))):
            if not hms:
                val = float(lspl[c])
            else:
                val = lspl[c]
            outarr[n].append(val)
    fin.close()
    for n in range(len(cols)):
        outarr[n] = np.array(outarr[n])
    return outarr

def parse_colnames(colnames,coords=None):
    """ Convert colnames input into list of column numbers.
    """
    cols = []
    if not isinstance(colnames,list):
        colnames = colnames.split(',')
    # parse column names from coords file and match to input values
    if coords is not None and fileutil.isFits(coords)[0]:
        # Open FITS file with table
        ftab = pyfits.open(coords)
        # determine which extension has the table
        for extn in ftab:
            if isinstance(extn,pyfits.BinTableHDU):
                # parse column names from table and match to inputs
                cnames = extn.columns.names
                if colnames is not None:
                    for c in colnames:
                        for name,i in zip(cnames,xrange(len(cnames))):
                            if c == name.lower(): cols.append(i)
                    if len(cols) < len(colnames):
                        errmsg = "Not all input columns found in table..."
                        ftab.close()
                        raise ValueError, errmsg
                else:
                    cols = cnames[:2]
                break
        ftab.close()
    else:
        for c in colnames:
            if isinstance(c, str):
                if c[0].lower() == 'c': cols.append(int(c[1:])-1)
                else:
                    cols.append(int(c))
            else:
                if isinstance(c, int):
                    cols.append(c)
                else:
                    errmsg = "Unsupported column names..."
                    raise ValueError, errmsg
    return cols


def createFile(dataArray=None, outfile=None, header=None):
    """Create a simple fits file for the given data array and header."""

    try:
        assert(outfile != None), "Please supply an output filename for createFile"
        assert(dataArray != None), "Please supply a data array for createFiles"
    except AssertionError:
        raise AssertionError

    print 'Creating output : ',outfile

    try:
        # Create the output file
        fitsobj = pyfits.HDUList()
        if (header != None):
            del(header['NAXIS1'])
            del(header['NAXIS2'])
            if header.has_key('XTENSION'):
                del(header['XTENSION'])
            if header.has_key('EXTNAME'):
                del(header['EXTNAME'])
            if header.has_key('EXTVER'):
                del(header['EXTVER'])

            if header.has_key('NEXTEND'):
                header['NEXTEND'] = 0

            hdu = pyfits.PrimaryHDU(data=dataArray,header=header)
            del hdu.header['PCOUNT']
            del hdu.header['GCOUNT']

        else:
            hdu = pyfits.PrimaryHDU(data=dataArray)

        fitsobj.append(hdu)
        fitsobj.writeto(outfile)

    finally:
        # CLOSE THE IMAGE FILES
        fitsobj.close()
        del fitsobj
