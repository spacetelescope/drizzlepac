""" ImageReg - A replacement for IRAF-based 'tweakshifts'
"""
import os

from stsci.tools import parseinput, teal
from stwcs import updatewcs

import util

# __version__ and __vdate__ are defined here, prior to the importing
# of the modules below, so that those modules can use the values 
# from these variable definitions, allowing the values to be designated 
# in one location only.
__version__ = '0.6.10'
__vdate__ = '10-Jan-2012'

import tweakutils
import imgclasses
import catalogs
import imagefindpars
    
__taskname__ = 'tweakreg' # unless someone comes up with anything better

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


def _managePsets(configobj):
    """ Read in parameter values from PSET-like configobj tasks defined for
        source-finding algorithms, and any other PSET-like tasks under this task,
        and merge those values into the input configobj dictionary.
    """
    # Merge all configobj instances into a single object
    configobj['SOURCE FINDING PARS'] = {}
        
    iparsobj = teal.teal(imagefindpars.__taskname__,loadOnly=True)

    # merge these parameters into full set
    configobj['SOURCE FINDING PARS'].merge(iparsobj)
        
    # clean up configobj a little to make it easier for later...
    if '_RULES_' in configobj:
        del configobj['_RULES_']

    
def edit_imagefindpars():
    """ Allows the user to edit the imagefindpars configObj in a TEAL GUI
    """
    iparsobj = teal.teal(imagefindpars.__taskname__, returnDict=False, loadOnly=False, canExecute=False)
        
        
@util.with_logging
def run(configobj):
    """ Primary Python interface for image registration code
        This task replaces 'tweakshifts'
    """
    print 'TweakReg Version %s(%s) started at: %s \n'%(
                    __version__,__vdate__,util._ptime()[0])
    util.print_pkg_versions()
    
    # Manage PSETs for source finding algorithms
    _managePsets(configobj)

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

    if configobj['updatewcs']:
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
    # Update parameter set with 'SOURCE FINDING PARS' now
    catfile_kwargs.update(configobj['SOURCE FINDING PARS'])
    uphdr_par = configobj['UPDATE HEADER']
    hdrlet_par = configobj['HEADERLET CREATION']
    hdrlet_par.update(uphdr_par) # default hdrlet name
    catfile_kwargs['updatehdr'] = uphdr_par['updatehdr']
    
    print '\n'+__taskname__+': Finding shifts for ',filenames,'\n'

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

    # otherwise, extract the catalog from the first input image source list
    if configobj['refimage'] not in [None, '',' ','INDEF']: # User specified an image to use
        refimg = imgclasses.Image(configobj['refimage'],**catfile_kwargs)
        refwcs = refimg.get_wcs()
        ref_source = refimg.all_radec
        refwcs_fname = refwcs[0].filename
    else:
        refwcs = []
        for i in input_images:
            refwcs.extend(i.get_wcs())
        ref_source = input_images[0].all_radec
        refwcs_fname = input_images[0].name

    print '\n'+'='*20+'\n'
    print 'Aligning all input images to WCS defined by ',refwcs_fname
    print '\n'+'='*20+'\n'
    
    # Check to see whether the user specified a separate catalog 
    #    of reference source positions and replace default source list with it
    refcat_par = configobj['REFERENCE CATALOG DESCRIPTION']
    if refcat_par['refcat'] not in [None,'',' ','INDEF']: # User specified a catalog to use
        ref_source = refcat_par['refcat']
        # Update kwargs with reference catalog parameters
        kwargs.update(refcat_par)

    try:
        # Create Reference Catalog object
        refimage = imgclasses.RefImage(refwcs,ref_source,**kwargs)
    except KeyboardInterrupt:
        refimage.close()
        for img in input_images:
            img.close()
        print 'Quitting as a result of user request (Ctrl-C)...'
        return

    if refimage.outxy is not None:    
        try:
            # Now, apply reference WCS to each image's sky positions as well as the
            # reference catalog sky positions,
            # then perform the fit between the reference catalog positions and 
            #    each image's positions    
            quit_immediately = False
            for img in input_images:
                print '\n'+'='*20
                print 'Performing fit for: ',img.name,'\n'
                img.match(refimage.outxy, refimage.wcs,refimage.name,
                            **configobj['OBJECT MATCHING PARAMETERS'])
                configobj['CATALOG FITTING PARAMETERS']['minobj'] = \
                                configobj['OBJECT MATCHING PARAMETERS']['minobj']
                img.performFit(**configobj['CATALOG FITTING PARAMETERS'])
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
                # write out shiftfile (if specified)
                shiftpars = configobj['OPTIONAL SHIFTFILE OUTPUT']
                if shiftpars['shiftfile']:
                    tweakutils.write_shiftfile(input_images,shiftpars['outshifts'],outwcs=shiftpars['outwcs'])
        except KeyboardInterrupt:
            refimage.close()
            for img in input_images:
                img.close()
            print 'Quitting as a result of user request (Ctrl-C)...'
            return
    else:
        print 'No valid sources in reference frame. Quitting...'
        return
# 
# Primary interface for running this task from Python
#
def TweakReg(files, editpars=False, configObj=None, **input_dict):
    """
    """
    # support input of filenames from command-line without a parameter name
    # then copy this into input_dict for merging with TEAL ConfigObj parameters
    if not util.is_blank(files):
        if input_dict is None:
            input_dict = {}
        input_dict['input'] = files
    # If called from interactive user-interface, configObj will not be 
    # defined yet, so get defaults using EPAR/TEAL.
    #
    # Also insure that the input_dict (user-specified values) are folded in
    # with a fully populated configObj instance.
    try:
        configObj = util.getDefaultConfigObj(__taskname__, configObj,
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

def help():
    print getHelpAsString(docstring=True)
    
# Append help file as docstring for use in Sphinx-generated documentation/web pages
TweakReg.__doc__ = getHelpAsString(docstring=True)
