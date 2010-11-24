""" ImageReg - A replacement for IRAF-based 'tweakshifts'
"""
import os

from pytools import parseinput, teal

import betadrizzle # for util.getDefaultConfigObj()

import tweakutils
import imgclasses
import catalogs
import sextractorpars
import imagefindpars

try:
    import sextractor
except:
    sextractor = None
    
__version__ = '0.1'
__vdate__ = '20-Sep-2010'

__taskname__ = 'tweakreg' # unless someone comes up with anything better


# 
# Interfaces used by TEAL
#
def getHelpAsString():
    helpString = 'TweakReg Version %s(%s)\n'%(__version__,__vdate__)
    helpString += teal.getHelpFileAsString(__taskname__,__file__)
    
    return helpString

def _managePsets(configobj):
    """ Read in parameter values from PSET-like configobj tasks defined for
        source-finding algorithms, and any other PSET-like tasks under this task,
        and merge those values into the input configobj dictionary.
    """
    # Merge all configobj instances into a single object
    configobj['SOURCE FINDING PARS'] = {}
    if configobj['findmode'] == 'imagefind' or (configobj['findmode'] == 'sextractor' and sextractor is None):
        if configobj['findmode'] == 'sextractor' and sextractor is None:
            print '**********\nWARNING\n**********'
            print 'No version of SExtractor found!'
            print ' Defaulting to basic source finding algorithm...'
            print '**********'
            # set the edit parameter to be the same as used for SExtractor 
            configobj['editipars'] = configobj['editspars']
            
        # Load parameters for this source finding mode
        iloadonly = not(configobj['editipars'])
        
        iparsobj = teal.teal(imagefindpars.__taskname__,loadOnly=iloadonly,canExecute=False)
        # Catch the case when the user hit 'CANCEL' for the PSET, we still the parameter values
        if iparsobj is None:
            iparsobj = teal.teal(imagefindpars.__taskname__,loadOnly=True)

        # merge these parameters into full set
        configobj['SOURCE FINDING PARS'].merge(iparsobj)

    elif configobj['findmode'] == 'sextractor':
        # Load parameters for this source finding mode
        sloadonly = not(configobj['editspars'])

        sparsobj = teal.teal(sextractorpars.__taskname__,loadOnly=sloadonly,canExecute=False)
        # Catch the case when the user hit 'CANCEL' for the PSET, we still the parameter values
        if sparsobj is None:
            sparsobj = teal.teal(sextractorpars.__taskname__,loadOnly=True)

        # merge these parameters into full set
        configobj['SOURCE FINDING PARS'].merge(sparsobj['USER SUPPLIED PARAMETERS'])
        
    # clean up configobj a little to make it easier for later...
    if '_RULES_' in configobj:
        del configobj['_RULES_']
    
def run(configobj):
    """ Primary Python interface for image registration code
        This task replaces 'tweakshifts'
    """
    # Manage PSETs for source finding algorithms
    _managePsets(configobj)

    # start interpretation of input parameters
    input = configobj['input']
    # Start by interpreting the inputs
    use_catfile = True
    filenames,catnames = tweakutils.parse_input(input)
    if catnames in [None,'',' ','INDEF'] or len(catnames) == 0:
        catfile_par = configobj['COORDINATE FILE DESCRIPTION']['catfile']
        # check to see whether the user specified input catalogs through other parameters
        if catfile_par not in [None,'',' ','INDEF']:
            # read in catalog file list provided by user
            catnames = tweakutils.parse_atfile_cat('@'+catfile_par)
        else:
            use_catfile = False
        
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
    kwargs = tweakutils.get_configobj_root(configobj)
    # define default value for 'xyunits' assuming sources to be derived from image directly
    kwargs['xyunits'] = 'pixels' # initialized here, required by Image class
    if use_catfile:
        # reset parameters based on parameter settings in this section
        kwargs.update(configobj['COORDINATE FILE DESCRIPTION'])
    # Update parameter set with 'SOURCE FINDING PARS' now
    kwargs.update(configobj['SOURCE FINDING PARS'])

    print '\n'+__taskname__+': Finding shifts for ',filenames,'\n'
    
    for imgnum in xrange(len(filenames)):
        # Create Image instances for all input images
        input_images.append(imgclasses.Image(filenames[imgnum],input_catalogs=catnames[imgnum],**kwargs))
    
    # create set of parameters to pass to RefImage class
    kwargs = tweakutils.get_configobj_root(configobj)
    # Determine a reference image or catalog and 
    #    return the full list of RA/Dec positions
    refcat_par = configobj['REFERENCE CATALOG DESCRIPTION']
    if refcat_par['refcat'] not in [None,'',' ','INDEF']: # User specified a catalog to use
        ref_source = refcat_par['refcat']
        # Update kwargs with reference catalog parameters
        kwargs.update(refcat_par)
    else: # otherwise, extract the catalog from the first input image source list
        ref_source = input_images[0].all_radec
        
    # Determine what WCS needs to be used for reference tangent plane
    if configobj['refimage'] not in [None, '',' ','INDEF']: # User specified an image to use
        refwcs = configobj['refimage']
    else:
        refwcs = []
        for i in input_images:
            refwcs.extend(i.get_wcs())

    # Create Reference Catalog object
    refimage = imgclasses.RefImage(refwcs,ref_source,**kwargs)

    # Now, apply reference WCS to each image's sky positions as well as the
    # reference catalog sky positions,
    # then perform the fit between the reference catalog positions and 
    #    each image's positions    
    for img in input_images:
        img.match(refimage.outxy, refimage.wcs,**configobj['OBJECT MATCHING PARAMETERS'])
        img.performFit(**configobj['CATALOG FITTING PARAMETERS'])
        if configobj['updatehdr']:
            img.updateHeader()
        if configobj['clean']:
            img.clean()

    # write out shiftfile (if specified)
    if configobj['output'] not in [None,'',' ','INDEF']:
        tweakutils.write_shiftfile(input_images,configobj['output'],outwcs=configobj['outwcs'])

# 
# Primary interface for running this task from Python
#
def TweakReg(files, editpars=False, configObj=None, **input_dict):
    """
    """
    # support input of filenames from command-line without a parameter name
    # then copy this into input_dict for merging with TEAL ConfigObj parameters
    if files not in ['',' ','INDEF', None]:
        if input_dict is None:
            input_dict = {}
        input_dict['input'] = files
    # If called from interactive user-interface, configObj will not be 
    # defined yet, so get defaults using EPAR/TEAL.
    #
    # Also insure that the input_dict (user-specified values) are folded in
    # with a fully populated configObj instance.
    configObj = betadrizzle.util.getDefaultConfigObj(__taskname__,configObj,input_dict,loadOnly=(not editpars))
    if configObj is None:
        return
    # If 'editpars' was set to True, util.getDefaultConfigObj() will have already
    # called 'run()'.
    if editpars == False:
        # Pass full set of parameters on to the task
        run(configObj)

# Append help file as docstring for use in Sphinx-generated documentation/web pages
TweakReg.__doc__ = getHelpAsString()
