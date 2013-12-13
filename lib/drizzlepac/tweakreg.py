""" ImageReg - A replacement for IRAF-based 'tweakshifts'
"""
import os

from stsci.tools import parseinput, teal
from stsci.tools import logutil, textutil
from stwcs import updatewcs

import util

# __version__ and __vdate__ are defined here, prior to the importing
# of the modules below, so that those modules can use the values
# from these variable definitions, allowing the values to be designated
# in one location only.
#
# This is specifically NOT intended to match the package-wide version information.
__version__ = '1.2.3'
__vdate__ = '11-Dec-2013'

import tweakutils
import imgclasses
import catalogs
import imagefindpars

__taskname__ = 'tweakreg' # unless someone comes up with anything better

PSET_SECTION = '_SOURCE FINDING PARS_'
# !! Use pre and post underscores to hide this added section in TEAL, so that
# TEAL allows its data to stay alongside the expected data during a call to
# TweakReg().  All of this needs to be revisited.

log = logutil.create_logger(__name__)

#
# Interfaces used by TEAL
#
def getHelpAsString(docstring=False):
    """
    return useful help from a file in the script directory called __taskname__.help
    """
    install_dir = os.path.dirname(__file__)
    htmlfile = os.path.join(install_dir,'htmlhelp',__taskname__+'.html')
    helpfile = os.path.join(install_dir,__taskname__+'.help')
    if docstring or (not docstring and not os.path.exists(htmlfile)):
        helpString = __taskname__+' Version '+__version__+' updated on '+__vdate__+'\n\n'
        if os.path.exists(helpfile):
            helpString += teal.getHelpFileAsString(__taskname__,__file__)
    else:
        helpString = 'file://'+htmlfile

    return helpString


def _managePsets(configobj,iparsobj=None):
    """ Read in parameter values from PSET-like configobj tasks defined for
        source-finding algorithms, and any other PSET-like tasks under this task,
        and merge those values into the input configobj dictionary.
    """
    # Merge all configobj instances into a single object
    configobj[PSET_SECTION] = {}

    if iparsobj is None:
        iparsobj = teal.load(imagefindpars.__taskname__)
        del iparsobj['_task_name_']

    # merge these parameters into full set
    configobj[PSET_SECTION].merge(iparsobj)

    # clean up configobj a little to make it easier for later...
#   if '_RULES_' in configobj:
#       del configobj['_RULES_']

def edit_imagefindpars():
    """ Allows the user to edit the imagefindpars configObj in a TEAL GUI
    """
    teal.teal(imagefindpars.__taskname__, returnAs=None,
              autoClose=True, loadOnly=False, canExecute=False)


@util.with_logging
def run(configobj):
    """ Primary Python interface for image registration code
        This task replaces 'tweakshifts'
    """
    print 'TweakReg Version %s(%s) started at: %s \n'%(
                    __version__,__vdate__,util._ptime()[0])
    util.print_pkg_versions()

    # Check to see whether or not the imagefindpars parameters have
    # already been loaded, as done through the python interface
    if PSET_SECTION not in configobj:
        # Manage PSETs for source finding algorithms
        _managePsets(configobj)

    # print out user set input parameter values for running this task
    log.info("USER INPUT PARAMETERS common to all Processing Steps:")
    util.printParams(configobj, log=log)

    # start interpretation of input parameters
    input = configobj['input']
    # Start by interpreting the inputs
    use_catfile = True
    filenames,catnames = tweakutils.parse_input(input)

    if not filenames:
        print 'No filenames matching input %r were found.' % input
        raise IOError

    # Verify that files are writable (based on file permissions) so that
    #    they can be updated if either 'updatewcs' or 'updatehdr' have
    #    been turned on (2 cases which require updating the input files)
    if configobj['updatewcs'] or configobj['UPDATE HEADER']['updatehdr']:
        filenames = util.verifyFilePermissions(filenames)
        if filenames is None or len(filenames) == 0:
            raise IOError

    if configobj['UPDATE HEADER']['updatehdr']:
        wname = configobj['UPDATE HEADER']['wcsname']
        # verify that a unique WCSNAME has been specified by the user
        for fname in filenames:
            uniq = util.verifyUniqueWcsname(fname,wname)
            if not uniq:
                errstr = 'WCSNAME "%s" already present in "%s".  '%(wname,fname)+\
                'A unique value for the "wcsname" parameter needs to be ' + \
                'specified. \n\nQuitting!'
                print textutil.textbox(errstr,width=60)
                raise IOError

    if configobj['updatewcs']:
        print '\nRestoring WCS solutions to original state using updatewcs...\n'
        updatewcs.updatewcs(filenames)

    if catnames in [None,'',' ','INDEF'] or len(catnames) == 0:
        catfile_par = configobj['COORDINATE FILE DESCRIPTION']['catfile']
        # check to see whether the user specified input catalogs through other parameters
        if catfile_par not in [None,'',' ','INDEF']:
            # read in catalog file list provided by user
            catnames = tweakutils.parse_atfile_cat('@'+catfile_par)
        else:
            use_catfile = False

    if 'exclusions' in configobj and \
        configobj['exclusions'] not in [None,'',' ','INDEF']:
        if os.path.exists(configobj['exclusions']):
            exclusion_files = tweakutils.parse_atfile_cat(
                '@'+configobj['exclusions'])
        else:
            print 'Could not find specified exclusions file "%s"'%(configobj['exclusions'])
            raise IOError
    else:
        exclusion_files = [None]*len(filenames)

    # Verify that we have the same number of catalog files as input images
    if catnames is not None and (len(catnames) > 0):
        rcat = configobj['REFERENCE CATALOG DESCRIPTION']['refcat']
        if (len(catnames) != len(filenames)):
            print 'The number of input catalogs does not match the number of input images'
            print 'Catalog files specified were:'
            print catnames
            print 'Input images specified were:'
            print filenames
            raise IOError
    else:
        # setup array of None values as input to catalog parameter for Image class
        catnames = [None]*len(filenames)
        use_catfile = False

    # convert input images and any provided catalog file names into Image objects
    input_images = []
    # copy out only those parameters needed for Image class
    catfile_kwargs = tweakutils.get_configobj_root(configobj)
    # define default value for 'xyunits' assuming sources to be derived from image directly
    catfile_kwargs['xyunits'] = 'pixels' # initialized here, required by Image class
    del catfile_kwargs['exclusions']

    if use_catfile:
        # reset parameters based on parameter settings in this section
        catfile_kwargs.update(configobj['COORDINATE FILE DESCRIPTION'])
        for sort_par in imgclasses.sortKeys:
            catfile_kwargs['sort_'+sort_par] = catfile_kwargs[sort_par]
    # Update parameter set with 'SOURCE FINDING PARS' now
    catfile_kwargs.update(configobj[PSET_SECTION])
    uphdr_par = configobj['UPDATE HEADER']
    hdrlet_par = configobj['HEADERLET CREATION']
    objmatch_par = configobj['OBJECT MATCHING PARAMETERS']
    catfit_pars = configobj['CATALOG FITTING PARAMETERS']

    catfit_pars['minobj'] = objmatch_par['minobj']
    objmatch_par['residplot'] = catfit_pars['residplot']

    hdrlet_par.update(uphdr_par) # default hdrlet name
    catfile_kwargs['updatehdr'] = uphdr_par['updatehdr']

    shiftpars = configobj['OPTIONAL SHIFTFILE OUTPUT']

    # verify a valid hdrname was provided, if headerlet was set to True
    imgclasses.verify_hdrname(**hdrlet_par)

    print '\nFinding shifts for: '
    for f in filenames:
        print '    ',f
    print '\n'

    log.info("USER INPUT PARAMETERS for finding sources for each input:")
    util.printParams(catfile_kwargs, log=log)

    try:
        for imgnum in xrange(len(filenames)):
            # Create Image instances for all input images
            input_images.append(imgclasses.Image(filenames[imgnum],
                                input_catalogs=catnames[imgnum],
                                exclusions=exclusion_files[imgnum],
                                **catfile_kwargs))
    except KeyboardInterrupt:
        for img in input_images:
            img.close()
        print 'Quitting as a result of user request (Ctrl-C)...'
        return
    # create set of parameters to pass to RefImage class
    kwargs = tweakutils.get_configobj_root(configobj)
    # Determine a reference image or catalog and
    #    return the full list of RA/Dec positions
    # Determine what WCS needs to be used for reference tangent plane
    refcat_par = configobj['REFERENCE CATALOG DESCRIPTION']
    if refcat_par['refcat'] not in [None,'',' ','INDEF']: # User specified a catalog to use
        # Update kwargs with reference catalog parameters
        kwargs.update(refcat_par)

    # otherwise, extract the catalog from the first input image source list
    if configobj['refimage'] not in [None, '',' ','INDEF']: # User specified an image to use
        #refimg = imgclasses.Image(configobj['refimage'],**catfile_kwargs)
        # Check to see whether the user specified a separate catalog
        #    of reference source positions and replace default source list with it
        if refcat_par['refcat'] not in [None,'',' ','INDEF']: # User specified a catalog to use
            ref_source = refcat_par['refcat']
        else:
            refimg = imgclasses.Image(configobj['refimage'],**catfile_kwargs)
            ref_source = refimg.all_radec

        try:
            refimage = imgclasses.RefImage(configobj['refimage'],ref_source,**kwargs)
            refwcs = refimage.wcs
            ref_source = refimage.all_radec
            refwcs_fname = refwcs.filename
        except KeyboardInterrupt:
            refimage.close()
            for img in input_images:
                img.close()
            print 'Quitting as a result of user request (Ctrl-C)...'
            return
    else:
        refwcs = []
        for i in input_images:
            refwcs.extend(i.get_wcs())
        try:
            refimg = input_images[0]
            ref_source = refimg.all_radec
            refimage = imgclasses.RefImage(refwcs,ref_source,
                                        xycatalog=refimg.xy_catalog,**kwargs)
            refwcs_fname = refimg.name
        except KeyboardInterrupt:
            refimage.close()
            for img in input_images:
                img.close()
            print 'Quitting as a result of user request (Ctrl-C)...'
            return

    print '\n'+'='*20+'\n'
    print 'Aligning all input images to WCS defined by ',refwcs_fname
    print '\n'+'='*20+'\n'

    if refimage.outxy is not None:
        try:
            log.info("USER INPUT PARAMETERS for matching sources:")
            util.printParams(objmatch_par, log=log)

            log.info("USER INPUT PARAMETERS for fitting source lists:")
            util.printParams(configobj['CATALOG FITTING PARAMETERS'], log=log)

            if hdrlet_par['headerlet']:
                log.info("USER INPUT PARAMETERS for creating headerlets:")
                util.printParams(hdrlet_par, log=log)

            if shiftpars['shiftfile']:
                log.info("USER INPUT PARAMETERS for creating a shiftfile:")
                util.printParams(shiftpars, log=log)

            # Now, apply reference WCS to each image's sky positions as well as the
            # reference catalog sky positions,
            # then perform the fit between the reference catalog positions and
            #    each image's positions
            quit_immediately = False
            xycat_lines = ''
            xycat_filename = None
            for img in input_images:
                if xycat_filename is None:
                    xycat_filename = img.rootname+'_xy_catfile.list'
                # Keep a record of all the generated input_xy catalogs
                xycat_lines += img.get_xy_catnames()

                print '\n'+'='*20
                print 'Performing fit for: ',img.name,'\n'
                img.match(refimage, **objmatch_par)

                img.performFit(**catfit_pars)
                if img.quit_immediately:
                    quit_immediately = True
                    img.close()
                    break
                img.updateHeader(wcsname=uphdr_par['wcsname'])
                if hdrlet_par['headerlet']:
                    img.writeHeaderlet(**hdrlet_par)
                if configobj['clean']:
                    img.clean()
                img.close()

            if not quit_immediately:
                if configobj['writecat'] and not configobj['clean']:
                    # Write out catalog file recording input XY catalogs used
                    # This file will be suitable for use as input to 'tweakreg'
                    # as the 'catfile' parameter
                    if os.path.exists(xycat_filename): os.remove(xycat_filename)
                    f=open(xycat_filename,mode='w')
                    f.writelines(xycat_lines)
                    f.close()

                # write out shiftfile (if specified)
                if shiftpars['shiftfile']:
                    tweakutils.write_shiftfile(input_images,shiftpars['outshifts'],outwcs=shiftpars['outwcs'])
        except KeyboardInterrupt:
            refimage.close()
            for img in input_images:
                img.close()
                del img
            print 'Quitting as a result of user request (Ctrl-C)...'
            return
    else:
        print 'No valid sources in reference frame. Quitting...'
        return

#
# Primary interface for running this task from Python
#
def TweakReg(files=None, editpars=False, configobj=None, imagefindcfg=None,
                **input_dict):
    """
    """
    # support input of filenames from command-line without a parameter name
    # then copy this into input_dict for merging with TEAL ConfigObj parameters
    if files is None and configobj is None:
        raise TypeError('TweakReg() needs either "files" or "configobj" arg')

    if files and not util.is_blank(files):
        if input_dict is None:
            input_dict = {}
        input_dict['input'] = files

    # Get default or user-specified configobj for primary task
    if isinstance(configobj, str):
        if not os.path.exists(configobj):
            print 'Cannot find .cfg file: '+configobj
            return
        configobj = teal.load(configobj, strict=False)

    if configobj is None:
        configobj = teal.load(__taskname__)

    # Merge PSET configobj with full task configobj
    _managePsets(configobj,iparsobj=imagefindcfg)
    # !! NOTE - the above line needs to be done so that getDefaultConfigObj()
    # can merge in input_dict, however the TEAL GUI is not going to understand
    # the extra section, (or use it), so work needs to be done here - some
    # more thinking about how to accomplish what we want to with PSETS.
    # !! For now, warn them that imagefindpars args will be ignored in GUI.
    if editpars and input_dict:
        idkeys = input_dict.keys()
        for i in idkeys:
            if i in configobj[PSET_SECTION]:
                print 'WARNING: ignoring imagefindpars setting "'+i+ \
                     '='+str(input_dict[i])+'", for now please enter directly into TEAL.'
                input_dict.pop(i)
        del configobj[PSET_SECTION] # force run() to pull it again after GUI use

    # If called from interactive user-interface, configObj will not be
    # defined yet, so get defaults using EPAR/TEAL.
    #
    # Also insure that the input_dict (user-specified values) are folded in
    # with a fully populated configObj instance.
    try:
        configObj = util.getDefaultConfigObj(__taskname__, configobj,
                                            input_dict,
                                            loadOnly=(not editpars))
    except ValueError:
        print "Problem with input parameters. Quitting..."
        return

    if configObj is None:
        return
    # If 'editpars' was set to True, util.getDefaultConfigObj() will have already
    # called 'run()'.
    if editpars == False:
        # Pass full set of parameters on to the task
        run(configObj)

def help(file=None):
    """
    Print out syntax help for running tweakreg

    Parameters
    ----------
    file : str (Default = None)
        If given, write out help to the filename specified by this parameter
        Any previously existing file with this name will be deleted before
        writing out the help.

    """
    helpstr = getHelpAsString(docstring=True)
    if file is None:
        print helpstr
    else:
        if os.path.exists(file): os.remove(file)
        f = open(file,mode='w')
        f.write(helpstr)
        f.close()

# Append help file as docstring for use in Sphinx-generated documentation/web pages
TweakReg.__doc__ = getHelpAsString(docstring=True)
