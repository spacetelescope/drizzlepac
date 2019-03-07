#!/usr/bin/env python

# $Id: updatenpol.py 8609 2010-01-19 16:22:48Z hack $

"""
`updatenpol`: Update the header of ACS file(s) with the names of new
``NPOLFILE`` and ``D2IMFILE`` reference files for use with the
C version of MultiDrizzle (astrodrizzle).

:Authors: Warren Hack

:License: :doc:`LICENSE`


:Usage:
    This task can be run from the operating system command line with::

        updatenpol [options] input [refdir]

:Command-line Options:
    `input`
        The specification of the files to be updated, either as a single filename,
        an ASN table name, or wild-card specification
        of a list of files.
    `refdir`
        The name of the directory containing all the new reference files
        (``*_npl.fits`` and ``*_d2i.fits`` files).
        If no directory is given, it will look in `jref$` by default.

    ``-h``
        Print the help (this text).

    ``-l``
        If specified, copy NPOLFILEs and D2IMFILEs to local directory
        for use with the input files.

    ``-i``
        If specified, the program will interactively request the exact
        names of the ``NPOLFILE`` and ``D2IMFILE`` reference files to be used
        for updating the header of each file. The value of 'refdir'
        will be ignored in interactive mode.


.. warning:: It will ask for the names of the ``NPOLFILE`` and ``D2IMFILE`` for
             EACH separate INPUT file when the option `-i` has been specified.

:Example:
    1. This command will update all the FLT files in the current directory
    with the new ``NPOLFILE`` and ``D2IMFILE`` reference files found in the 'myjref'
    directory as defined in the environment::

        updatenpol *flt.fits myjref$


:Compatability with MultiDrizzle:
    The new version of ``MultiDrizzle`` (``AstroDrizzle``) and `updatewcs`
    only work with the new ``NPOLFILE`` reference file for the ``DGEO`` correction
    (to replace the use of DGEOFILE).
    In fact, ``AstroDrizzle`` has been extensively modified to
    prompt the user with a very lengthy explanation on whether it should
    stop and allow the user to update the header or continue without
    applying the ``DGEO`` correction under circumstances when the ``NPOLFILE``
    keyword can not be found for ACS.

"""
__docformat__ = 'restructuredtext'

__taskname__ = "updatenpol"

# This is specifically NOT intended to match the package-wide version information.
__version__ = '1.1.0'
__version_date__ = '16-Aug-2011'

import os,sys,shutil

from astropy.io import fits
from stsci.tools import fileutil as fu
from stsci.tools import parseinput
from stsci.tools import teal

from stwcs import updatewcs
from . import util


def update(input,refdir="jref$",local=None,interactive=False,wcsupdate=True):
    """
    Updates headers of files given as input to point to the new reference files
    NPOLFILE and D2IMFILE required with the new C version of MultiDrizzle.

    Parameters
    -----------
    input : string or list
                Name of input file or files acceptable forms:
                  - single filename with or without directory
                  - @-file
                  - association table
                  - python list of filenames
                  - wildcard specification of filenames

    refdir : string
                Path to directory containing new reference files, either
                environment variable or full path.

    local : boolean
                Specifies whether or not to copy new reference files to local
                directory for use with the input files.

    interactive : boolean
                Specifies whether or not to interactively ask the user for the
                exact names of the new reference files instead of automatically
                searching a directory for them.

    updatewcs : boolean
                Specifies whether or not to update the WCS information in this
                file to use the new reference files.

    Examples
    --------
    1. A set of associated images specified by an ASN file can be updated to use
       the NPOLFILEs and D2IMFILE found in the local directory defined using
       the `myjref$` environment variable under PyRAF using::

            >>> import updatenpol
            >>> updatenpol.update('j8bt06010_asn.fits', 'myref$')

    2. Another use under Python would be to feed it a specific list of files
       to be updated using::

          >>> updatenpol.update(['file1_flt.fits','file2_flt.fits'],'myjref$')

    3. Files in another directory can also be processed using::

          >>> updatenpol.update('data$*flt.fits','../new/ref/')

    Notes
    -----
    .. warning::
        This program requires access to the `jref$` directory in order
        to evaluate the DGEOFILE specified in the input image header.
        This evaluation allows the program to get the information it
        needs to identify the correct NPOLFILE.

    The use of this program now requires that a directory be set up with
    all the new NPOLFILE and D2IMFILE reference files for ACS (a single
    directory for all files for all ACS detectors will be fine, much like
    jref).  Currently, all the files generated by the ACS team has initially
    been made available at::

        /grp/hst/acs/lucas/new-npl/

    The one known limitation to how this program works comes from
    confusion if more than 1 file could possibly be used as the new
    reference file. This would only happen when NPOLFILE reference files
    have been checked into CDBS multiple times, and there are several
    versions that apply to the same detector/filter combination.  However,
    that can be sorted out later if we get into that situation at all.

    """
    print('UPDATENPOL Version',__version__+'('+__version_date__+')')
    # expand (as needed) the list of input files
    files,fcol = parseinput.parseinput(input)

    if not interactive:
        # expand reference directory name (if necessary) to
        # interpret IRAF or environment variable names
        rdir = fu.osfn(refdir)
        ngeofiles,ngcol = parseinput.parseinput(os.path.join(rdir,'*npl.fits'))
        # Find D2IMFILE in refdir for updating input file header as well
        d2ifiles,d2col = parseinput.parseinput(os.path.join(rdir,"*d2i.fits"))

    # Now, build a matched list of input files and DGEOFILE reference files
    # to use for selecting the appropriate new reference file from the
    # refdir directory.
    for f in files:
        print('Updating: ',f)
        fdir = os.path.split(f)[0]
        # Open each file...
        fimg = fits.open(f, mode='update', memmap=False)
        phdr = fimg['PRIMARY'].header
        fdet = phdr['detector']
        # get header of DGEOFILE
        dfile = phdr.get('DGEOFILE','')
        if dfile in ['N/A','',' ',None]:
            npolname = ''
        else:
            dhdr = fits.getheader(fu.osfn(dfile), memmap=False)
            if not interactive:
                # search all new NPOLFILEs for one that matches current DGEOFILE config
                npol = find_npolfile(ngeofiles,fdet,[phdr['filter1'],phdr['filter2']])
            else:
                if sys.version_info[0] >= 3:
                    npol = input("Enter name of NPOLFILE for %s:"%f)
                else:
                    npol = raw_input("Enter name of NPOLFILE for %s:"%f)
                if npol == "": npol = None

            if npol is None:
                errstr =  "No valid NPOLFILE found in "+rdir+" for detector="+fdet+"\n"
                errstr += " filters = "+phdr['filter1']+","+phdr['filter2']
                raise ValueError(errstr)

            npolname = os.path.split(npol)[1]
            if local:
                npolname = os.path.join(fdir,npolname)
                # clobber any previous copies of this reference file
                if os.path.exists(npolname): os.remove(npolname)
                shutil.copy(npol,npolname)
            else:
                if '$' in refdir:
                    npolname = refdir+npolname
                else:
                    npolname = os.path.join(refdir,npolname)
        phdr.set('NPOLFILE', value=npolname,
                 comment="Non-polynomial corrections in Paper IV LUT",
                 after='DGEOFILE')

        # Now find correct D2IFILE
        if not interactive:
            d2i = find_d2ifile(d2ifiles,fdet)
        else:
            if sys.version_info[0] >= 3:
                d2i = input("Enter name of D2IMFILE for %s:"%f)
            else:
                d2i = raw_input("Enter name of D2IMFILE for %s:"%f)
            if d2i == "": d2i = None

        if d2i is None:
            print('=============\nWARNING:')
            print("    No valid D2IMFILE found in "+rdir+" for detector ="+fdet)
            print("    D2IMFILE correction will not be applied.")
            print('=============\n')
            d2iname = ""
        else:
            d2iname = os.path.split(d2i)[1]
            if local:
                # Copy D2IMFILE to local data directory alongside input file as well
                d2iname = os.path.join(fdir,d2iname)
                # clobber any previous copies of this reference file
                if os.path.exists(d2iname): os.remove(d2iname)
                shutil.copy(d2i,d2iname)
            else:
                if '$' in refdir:
                    d2iname = refdir+d2iname
                else:
                    d2iname = os.path.join(refdir,d2iname)

        phdr.set('D2IMFILE', value=d2iname,
                 comment="Column correction table",
                 after='DGEOFILE')

        # Close this input file header and go on to the next
        fimg.close()

        if wcsupdate:
            updatewcs.updatewcs(f)

def find_d2ifile(flist,detector):
    """ Search a list of files for one that matches the detector specified.
    """
    d2ifile = None
    for f in flist:
        fdet = fits.getval(f, 'detector', memmap=False)
        if fdet == detector:
            d2ifile = f
    return d2ifile

def find_npolfile(flist,detector,filters):
    """ Search a list of files for one that matches the configuration
        of detector and filters used.
    """
    npolfile = None
    for f in flist:
        fdet = fits.getval(f, 'detector', memmap=False)
        if fdet == detector:
            filt1 = fits.getval(f, 'filter1', memmap=False)
            filt2 = fits.getval(f, 'filter2', memmap=False)
            fdate = fits.getval(f, 'date', memmap=False)
            if filt1 == 'ANY' or \
             (filt1 == filters[0] and filt2 == filters[1]):
                npolfile = f
    return npolfile

#
#### Interfaces used by TEAL
#
def run(configobj=None,editpars=False):
    """ Teal interface for running this code.
    """

    if configobj is None:
        configobj =teal.teal(__taskname__,loadOnly=(not editpars))

    update(configobj['input'],configobj['refdir'],
        local=configobj['local'],interactive=configobj['interactive'],
        wcsupdate=configobj['wcsupdate'])


def main():

    import getopt

    try:
        optlist, args = getopt.getopt(sys.argv[1:], 'hli')
    except getopt.error as e:
        print(str(e))
        print(__doc__)
        print("\t", __version__)

    # initialize default values
    help = 0
    local = False
    interactive = False

    # read options
    for opt, value in optlist:
        if opt == "-h":
            help = 1
        if opt == "-l":
            local = True
        if opt == "-i":
            interactive = True
    if len(args) < 2:
        args.append('jref$')
    if (help):
        print(__doc__)
        print("\t", __version__+'('+__version_date__+')')
    else:
        update(args[:-1],args[-1],local=local,interactive=interactive)

if __name__ == "__main__":
    main()
"""
    Copyright (C) 2003 Association of Universities for Research in Astronomy (AURA)

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

        1. Redistributions of source code must retain the above copyright
          notice, this list of conditions and the following disclaimer.

        2. Redistributions in binary form must reproduce the above
          copyright notice, this list of conditions and the following
          disclaimer in the documentation and/or other materials provided
          with the distribution.

        3. The name of AURA and its representatives may not be used to
          endorse or promote products derived from this software without
          specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
    WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
    INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
    BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
    OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
    TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
    USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
    DAMAGE.
"""

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
    taskname = util.base_taskname(__taskname__, '')
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
