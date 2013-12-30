#!/usr/bin/env python

"""
A library of utility functions

"""

from __future__ import division  # confidence medium
import logging
import functools
import os
import sys
import string
import errno

import numpy as np
import pyfits
from stsci.tools import asnutil, fileutil, teal, cfgpars, logutil
from stsci.tools import check_files
from stsci.tools import configobj

from stwcs import wcsutil
from stwcs.wcsutil import altwcs

from .version import *

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


def get_pool_size(usr_config_value=None, num_tasks=None):
    """ Determine size of thread/process-pool for parallel processing.
    This examines the cpu_count to decide and return the right pool
    size to use.  Also take into account the user's wishes via the config
    object value, if specified.  On top of that, don't allow the pool size
    returned to be any higher than the number of parallel tasks, if specified.
    Only use what we need (mp.Pool starts pool_size processes, needed or not).
    If number of tasks is unknown, call this with "num_tasks" set to None.
    Returns 1 when indicating that parallel processing should not be used.
    Consolidate all such logic here, not in the caller. """

    if not can_parallel:
        return 1
    # Give priority to their specified cfg value, over the actual cpu count
    if usr_config_value != None:
        if num_tasks == None:
            return usr_config_value
        else:
            # usr_config_value may be needlessly high
            return min(usr_config_value, num_tasks)
    # they haven't specified a cfg value, so go with the cpu_count
    if num_tasks == None:
        return _cpu_count
    else:
        # run no more workers than tasks
        return min(_cpu_count, num_tasks)


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

#
# Logging routines
#

_log_file_handler = None

def init_logging(logfile=DEFAULT_LOGNAME, default=None, level=logging.INFO):
    """
    Set up logger for capturing stdout/stderr messages.

    Must be called prior to writing any messages that you want to log.
    """

    if logfile == "INDEF":
        if not is_blank(default):
            logname = fileutil.buildNewRootname(default, '.log')
        else:
            logname = DEFAULT_LOGNAME
    elif logfile not in [None, "" , " "]:
        if logfile.endswith('.log'):
            logname = logfile
        else:
            logname = logfile + '.log'
    else:
        logname = None

    if logname is not None:
        logutil.setup_global_logging()
        # Don't use logging.basicConfig since it can only be called once in a
        # session
        # TODO: Would be fine to use logging.config.dictConfig, but it's not
        # available in Python 2.5
        global _log_file_handler
        root_logger = logging.getLogger()
        if _log_file_handler:
            root_logger.removeHandler(_log_file_handler)
        # Default mode is 'a' which is fine
        _log_file_handler = logging.FileHandler(logname)
        # TODO: Make the default level configurable in the task parameters
        _log_file_handler.setLevel(level)
        _log_file_handler.setFormatter(
            logging.Formatter('[%(levelname)-8s] %(message)s'))
        root_logger.addHandler(_log_file_handler)

        print 'Setting up logfile : ', logname

        stdout_logger = logging.getLogger('stsci.tools.logutil.stdout')
        # Disable display of prints to stdout from all packages except
        # drizzlepac
        stdout_logger.addFilter(logutil.EchoFilter(include=['drizzlepac']))
    else:
        print 'No trailer file created...'


def end_logging(filename=None):
    """
    Close log file and restore system defaults.
    """

    if logutil.global_logging_started:
        if filename:
            print 'Trailer file written to: ', filename
        else:
            # This generally shouldn't happen if logging was started with
            # init_logging and a filename was given...
            print 'No trailer file saved...'

        logutil.teardown_global_logging()
    else:
        print 'No trailer file saved...'


class WithLogging(object):
    def __init__(self):
        self.depth = 0

    def __call__(self, func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            from .processInput import processFilenames

            errorobj = None
            filename = None
            if self.depth == 0:
                # Setup logging for this task; otherwise we're being called as
                # a subtask and will use the existing log
                # The first arg must be the configobj
                try:
                    configobj = args[0]
                    input = processFilenames(configobj['input'])
                    inputs, output, _, _ = input
                    if output is not None:
                        default = output
                    elif inputs:
                        default = inputs[0]
                    else:
                        default = None

                    if default is not None:
                        # astrodrizzle and tweakreg have this parameter.
                        # Could we do away with it altogether?
                        if 'runfile' in configobj:
                            filename = configobj['runfile']
                        else:
                            filename = default
                    verbose_level=logging.INFO
                    if 'verbose' in configobj and configobj['verbose']:
                        verbose_level=logging.DEBUG
                    init_logging(filename,level=verbose_level)
                except (KeyError, IndexError, TypeError):
                    pass

            self.depth += 1

            # This looks utterly bizarre, but it seems to be the only way I can
            # ensure that any exceptions that occur in the wrapped function are
            # logged before teardown_global_logging() is called.  Unless the
            # except clause is explicitly included here, even with just 'pass'
            # under it, Python discards the sys.exc_info() data by the time the
            # finally clause is reached.

            try:
                func(*args, **kwargs)
            except Exception as errorobj:
                pass
            finally:
                self.depth -= 1
                if self.depth == 0:
                    end_logging(filename)
                    # Insure that any exception raised by the code gets passed on
                    # (hope that end_logging didn't change the last exception raised)
                    if errorobj:
                        raise

        return wrapper


with_logging = WithLogging()


def print_pkg_versions(packages=None, svn=False, log=None):
    if log is not None:
        def output(msg):
            log.info(msg)
    else:
        def output(msg):
            print msg

    pkgs = ['numpy', 'pyfits', 'stwcs', 'pywcs']
    if packages is not None:
        if not isinstance(packages, list):
            packages = [packages]
        pkgs.extend(packages)

    output('Version Information')
    output('-' * 20)
    output('Python Version %s' % sys.version)
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
        output(vstr)


class ProcSteps(object):
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
        """
        Print out a formatted summary of the elapsed times for all the
        performed steps.
        """

        self.end = _ptime()
        total_time = 0
        print ProcSteps.__report_header
        for step in self.order:
            if 'elapsed' in self.steps[step]:
                _time = self.steps[step]['elapsed']
            else:
                _time = 0.0
            total_time += _time
            print '   %20s          %0.4f sec.' % (step, _time)

        print '   %20s          %s' % ('=' * 20, '=' * 20)
        print '   %20s          %0.4f sec.' % ('Total', total_time)

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

def displayEmptyInputWarningBox(display=True, parent=None):
    """ Displays a warning box for the 'input' parameter.
    """
    import tkMessageBox

    if display:
        msg = 'No valid input files found! '+\
        'Please check the value for the "input" parameter.'
        tkMessageBox.showwarning(parent=parent,message=msg, title="No valid inputs!")
    return "yes"

def displayBadRefimageWarningBox(display=True, parent=None):
    """ Displays a warning box for the 'input' parameter.
    """
    import tkMessageBox

    if display:
        msg = 'No refimage with WCS found!\n '+\
        ' This could be caused by one of 2 problems:\n'+\
        '   * filename does not specify an extension with a valid WCS.\n'+\
        '   * can not find the file.\n'+\
        'Please check the filename specified in the "refimage" parameter.'

        tkMessageBox.showwarning(parent=parent,message=msg, title="No valid inputs!")
    return "yes"

def updateNEXTENDKw(fobj):
    """ Update NEXTEND keyword in PRIMARY header (if present) to accurately
        reflect the number of extensions in the MEF file.
    """
    if 'nextend' in fobj[0].header:
        fobj[0].header['nextend'] = len(fobj)-1

def count_sci_extensions(filename):
    """ Return the number of SCI extensions and the EXTNAME from a input MEF file.
    """
    num_sci = 0
    extname = 'SCI'
    num_ext = 0
    for extn in fileutil.openImage(filename):
        num_ext += 1
        if 'extname' in extn.header and extn.header['extname'] == extname:
            num_sci += 1
    if num_sci == 0:
        extname = 'PRIMARY'
        num_sci = 1

    return num_sci,extname

def verifyUniqueWcsname(fname,wcsname):
    """
    Report whether or not the specified WCSNAME already exists in the file
    """
    uniq = True
    numsci,extname = count_sci_extensions(fname)
    wnames = altwcs.wcsnames(fname,ext=(extname,1))

    if wcsname in wnames.values():
        uniq = False

    return uniq

def verifyUpdatewcs(fname):
    """
    Verify the existence of WCSNAME in the file.  If it is not present,
    report this to the user and raise an exception.  Returns True if WCSNAME
    was found in all SCI extensions.
    """
    updated = True
    numsci,extname = count_sci_extensions(fname)
    for n in range(1,numsci+1):
        hdr = pyfits.getheader(fname,extname=extname,extver=n)
        if 'wcsname' not in hdr:
            updated = False
            break
    return updated

def verifyRefimage(refimage):
    """
    Verify that the value of refimage specified by the user points to an
    extension with a proper WCS defined. It starts by making sure an extension gets
    specified by the user when using a MEF file. The final check comes by looking
    for a CD matrix in the WCS object itself. If either test fails, it returns
    a value of False.
    """
    valid = True

    # start by trying to see whether the code can even find the file
    if is_blank(refimage):
        valid=True
        return valid

    refroot = fileutil.parseFilename(refimage)[0]
    if not os.path.exists(refroot):
        valid = False
        return valid

    # start by checking to make sure user specified an extension specified
    # when using an MEF as refimage
    ftype = fileutil.isFits(refimage)
    if ftype[1] == 'mef' and '[' not in refimage:
        valid = False
    # if a MEF has been specified, make sure extension contains a valid WCS
    if valid:
        # check for CD matrix in WCS object
        refwcs = wcsutil.HSTWCS(refimage)
        if not refwcs.wcs.has_cd():
            valid = False
        else:
            valid = True
        del refwcs

    return valid

def verifyFilePermissions(filelist, chmod=True):
    """ Verify that images specified in 'filelist' can be updated.

    A message will be printed reporting the names of any images which
    do not have write-permission, then quit.
    """
    badfiles = []
    archive_dir = False
    for img in filelist:
        fname = fileutil.osfn(img)
        if 'OrIg_files' in os.path.split(fname)[0]:
            archive_dir = True
        try:
            fp = open(fname,mode='a')
            fp.close()
        except IOError as e:
            if e.errno == errno.EACCES:
                badfiles.append(img)
            # Not a permission error.
            pass

    num_bad = len(badfiles)
    if num_bad > 0:
        if archive_dir:
            print '\n'
            print '#'*40
            print '    Working in "OrIg_files" (archive) directory. '
            print '    This directory has been created to serve as an archive'
            print '    for the original input images. '
            print '\n    These files should be copied into another directory'
            print '     for processing. '
            print '#'*40

        print '\n'
        print '#'*40
        print 'Found %d files which can not be updated!'%(num_bad)
        for img in badfiles:
            print '    %s'%(img)
        print '\nPlease reset permissions for these files and restart...'
        print '#'*40
        print '\n'
        filelist = None

    return filelist

def getFullParList(configObj):
    """
    Return a single list of all parameter names included in the configObj
    regardless of which section the parameter was stored
    """
    plist = []
    for par in configObj.iterkeys():
        if isinstance(configObj[par],configobj.Section):
            plist.extend(getFullParList(configObj[par]))
        else:
            plist.append(par)
    return plist

def validateUserPars(configObj,input_dict):
    """ Compares input parameter names specified by user with those already
        recognized by the task.

        Any parameters provided by the user that does not match a known
        task parameter will be reported and a ValueError exception will be
        raised.
    """
    # check to see whether any input parameters are unexpected.
    # Any unexpected parameters provided on input should be reported and
    # the code should stop
    plist = getFullParList(configObj)
    extra_pars = []
    for kw in input_dict:
        if kw not in plist:
            extra_pars.append(kw)
    if len(extra_pars) > 0:
        print '='*40
        print 'The following input parameters were not recognized as valid inputs:'
        for p in extra_pars:
            print "    %s"%(p)
        print '\nPlease check the spelling of the parameter(s) and try again...'
        print '='*40
        raise ValueError

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
        # check to see whether any input parameters are unexpected.
        # Any unexpected parameters provided on input should be reported and
        # the code should stop
        validateUserPars(configObj,input_dict)

        # If everything looks good, merge user inputs with configObj and continue
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
        if key.find('STEP '+str(stepnum)+':') >= 0:
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


def printParams(paramDictionary, all=False, log=None):
    """
    Print nicely the parameters from the dictionary.
    """

    if log is not None:
        def output(msg):
            log.info(msg)
    else:
        def output(msg):
            print msg

    if not paramDictionary:
        output('No parameters were supplied')
    else:
        for key in sorted(paramDictionary):
            if all or (not isinstance(paramDictionary[key], dict)) \
            and key[0] != '_':
                output('\t' + '\t'.join([str(key) + ' :',
                                         str(paramDictionary[key])]))
        if log is None:
            output('\n')


def isASNTable(inputFilelist):
    """Return TRUE if inputFilelist is a fits ASN file."""
    if ("_asn"  or "_asc") in inputFilelist:
        return True
    return False

def isCommaList(inputFilelist):
    """Return True if the input is a comma separated list of names."""
    if isinstance(inputFilelist, int) or isinstance(inputFilelist, np.int32):
        ilist = str(inputFilelist)
    else:
        ilist = inputFilelist
    if "," in ilist:
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


def interpret_bits_value(val):
    """ Converts input bits value from string to a single integer value or None.
    If a comma-separated set of values are provided, they are summed.
    """
    if isinstance(val,int):
        intval = val
    else:
        val = str(val)
        intval = 0
        if isCommaList(val):
            valspl = val.split(',')
            for v in valspl:
                intval += int(v)
        elif '+' in val:
            valspl = val.split('+')
            for v in valspl:
                intval += int(v)
        elif is_blank(val):
            intval = None
        else:
            intval = int(val)

    return intval

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

    if 'expstart' in primary_hdr:
        exphdr = primary_hdr
    else:
        exphdr = header

    if 'EXPSTART' in exphdr:
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
    """Create a simple fits file for the given data array and header.
        Returns either the FITS object in-membory when outfile==None or
            None when the FITS file was written out to a file.
    """

    # Insure that at least a data-array has been provided to create the file
    assert(dataArray != None), "Please supply a data array for createFiles"

    try:
        # Create the output file
        fitsobj = pyfits.HDUList()
        if (header != None):
            del(header['NAXIS1'])
            del(header['NAXIS2'])
            if 'XTENSION' in header:
                del(header['XTENSION'])
            if 'EXTNAME' in header:
                del(header['EXTNAME'])
            if 'EXTVER' in header:
                del(header['EXTVER'])

            if 'NEXTEND' in header:
                header['NEXTEND'] = 0

            hdu = pyfits.PrimaryHDU(data=dataArray,header=header)
            del hdu.header['PCOUNT']
            del hdu.header['GCOUNT']

        else:
            hdu = pyfits.PrimaryHDU(data=dataArray)

        fitsobj.append(hdu)
        if outfile:
            fitsobj.writeto(outfile)
    finally:
        # CLOSE THE IMAGE FILES
        fitsobj.close()

        if outfile:
            del fitsobj
            fitsobj = None
    return fitsobj
