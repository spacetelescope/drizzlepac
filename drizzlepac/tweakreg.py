""" ``TweakReg`` - A replacement for IRAF-based ``tweakshifts``

:Authors: Warren Hack, Mihai Cara

:License: :doc:`LICENSE`

"""
import os
import sys
import numpy as np
from copy import copy

from stsci.tools import parseinput, teal
from stsci.tools import logutil, textutil
from stsci.tools.cfgpars import DuplicateKeyError
from stwcs import updatewcs

from . import util

# __version__ and __version_date__ are defined here, prior to the importing
# of the modules below, so that those modules can use the values
# from these variable definitions, allowing the values to be designated
# in one location only.
#
# This is specifically NOT intended to match the package-wide version information.
__version__ = '1.4.7'
__version_date__ = '18-April-2018'

from . import tweakutils
from . import imgclasses
from . import catalogs
from . import imagefindpars
from . import refimagefindpars

__taskname__ = 'tweakreg' # unless someone comes up with anything better

PSET_SECTION        = '_SOURCE FINDING PARS_'
PSET_SECTION_REFIMG = '_REF IMAGE SOURCE FINDING PARS_'
# !! Use pre and post underscores to hide this added section in TEAL, so that
# TEAL allows its data to stay alongside the expected data during a call to
# TweakReg().  All of this needs to be revisited.

log = logutil.create_logger(__name__, level=logutil.logging.NOTSET)


def _managePsets(configobj, section_name, task_name, iparsobj=None, input_dict=None):
    """ Read in parameter values from PSET-like configobj tasks defined for
        source-finding algorithms, and any other PSET-like tasks under this task,
        and merge those values into the input configobj dictionary.
    """
    # Merge all configobj instances into a single object
    configobj[section_name] = {}

    # Load the default full set of configuration parameters for the PSET:
    iparsobj_cfg = teal.load(task_name)

    # Identify optional parameters in input_dicts that are from this
    # PSET and add it to iparsobj:
    if input_dict is not None:
        for key in list(input_dict.keys()):
            if key in iparsobj_cfg:
                if iparsobj is not None and key in iparsobj:
                    raise DuplicateKeyError("Duplicate parameter '{:s}' "
                        "provided for task {:s}".format(key, task_name))
                iparsobj_cfg[key] = input_dict[key]
                del input_dict[key]

    if iparsobj is not None:
        iparsobj_cfg.update(iparsobj)

    del iparsobj_cfg['_task_name_']

    # merge these parameters into full set
    configobj[section_name].merge(iparsobj_cfg)

    # clean up configobj a little to make it easier for later...
#   if '_RULES_' in configobj:
#       del configobj['_RULES_']


def edit_imagefindpars():
    """ Allows the user to edit the imagefindpars configObj in a TEAL GUI
        """
    teal.teal(imagefindpars.__taskname__, returnAs=None,
              autoClose=True, loadOnly=False, canExecute=False)

def edit_refimagefindpars():
    """ Allows the user to edit the refimagefindpars configObj in a TEAL GUI
        """
    teal.teal(refimagefindpars.__taskname__, returnAs=None,
              autoClose=True, loadOnly=False, canExecute=False)


@util.with_logging
def run(configobj):
    """ Primary Python interface for image registration code
        This task replaces 'tweakshifts'
    """
    print('TweakReg Version %s(%s) started at: %s \n'%(
                    __version__,__version_date__,util._ptime()[0]))
    util.print_pkg_versions()

    # make sure 'updatewcs' is set to False when running from GUI or if missing
    # from configObj:
    if 'updatewcs' not in configobj:
        configobj['updatewcs'] = False

    # Check to see whether or not the imagefindpars parameters have
    # already been loaded, as done through the python interface.
    # Repeat for refimagefindpars
    if PSET_SECTION not in configobj:
        # Manage PSETs for source finding algorithms
        _managePsets(configobj, PSET_SECTION, imagefindpars.__taskname__)
    #print configobj[PSET_SECTION]
    if PSET_SECTION_REFIMG not in configobj:
        # Manage PSETs for source finding algorithms in reference image
        _managePsets(configobj, PSET_SECTION_REFIMG,
                     refimagefindpars.__taskname__)

    log.debug('')
    log.debug("==== TweakReg was invoked with the following parameters: ====")
    log.debug('')
    util.print_cfg(configobj, log.debug)

    # print out user set input parameter values for running this task
    log.info('')
    log.info("USER INPUT PARAMETERS common to all Processing Steps:")
    util.printParams(configobj, log=log)

    # start interpretation of input parameters
    input_files = configobj['input']
    # Start by interpreting the inputs
    use_catfile = True
    expand_refcat = configobj['expand_refcat']
    enforce_user_order = configobj['enforce_user_order']

    filenames, catnames = tweakutils.parse_input(
        input_files, sort_wildcards=not enforce_user_order
    )

    catdict = {}
    for indx,f in enumerate(filenames):
        if catnames is not None and len(catnames) > 0:
            catdict[f] = catnames[indx]
        else:
            catdict[f] = None

    if not filenames:
        print('No filenames matching input %r were found.' % input_files)
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
        if not configobj['UPDATE HEADER']['reusename']:
            for fname in filenames:
                uniq = util.verifyUniqueWcsname(fname,wname)
                if not uniq:
                    errstr = 'WCSNAME "%s" already present in "%s".  '%(wname,fname)+\
                    'A unique value for the "wcsname" parameter needs to be ' + \
                    'specified. \n\nQuitting!'
                    print(textutil.textbox(errstr,width=60))
                    raise IOError

    if configobj['updatewcs']:
        print('\nRestoring WCS solutions to original state using updatewcs...\n')
        updatewcs.updatewcs(filenames)

    if catnames in [None,'',' ','INDEF'] or len(catnames) == 0:
        catfile_par = configobj['COORDINATE FILE DESCRIPTION']['catfile']
        # check to see whether the user specified input catalogs through other parameters
        if catfile_par not in [None,'',' ','INDEF']:
            # read in catalog file list provided by user
            catnames,catdict = tweakutils.parse_atfile_cat('@'+catfile_par)
        else:
            use_catfile = False

        if 'exclusions' in configobj and \
            configobj['exclusions'] not in [None,'',' ','INDEF']:
            if os.path.exists(configobj['exclusions']):
                excl_files, excl_dict = tweakutils.parse_atfile_cat(
                    '@'+configobj['exclusions'])

                # make sure the dictionary is well formed and that keys are base
                # file names and that exclusion files have been expanded:
                exclusion_files = []
                exclusion_dict = {}
                rootpath = os.path.abspath(
                    os.path.split(configobj['exclusions'])[0]
                )

                for f in excl_dict.keys():
                    print(f)
                    bf = os.path.basename(f)
                    exclusion_files.append(bf)
                    reglist = excl_dict[f]

                    if reglist is None:
                        exclusion_dict[bf] = None
                        continue

                    new_reglist = []
                    for regfile in reglist:
                        if regfile in [ None, 'None', '', ' ', 'INDEF' ]:
                            new_reglist.append(None)
                        else:
                            abs_regfile = os.path.normpath(
                                os.path.join(rootpath, regfile)
                            )
                            new_reglist.append(abs_regfile)

                    exclusion_dict[bf] = new_reglist

            else:
                raise IOError('Could not find specified exclusions file "{:s}"'
                      .format(configobj['exclusions']))

        else:
            exclusion_files = [None]*len(filenames)
            exclusion_dict = {}

            for f in filenames:
                exclusion_dict[os.path.basename(f)] = None

    # Verify that we have the same number of catalog files as input images
    if catnames is not None and (len(catnames) > 0):
        missed_files = []

        for f in filenames:
            if f not in catdict:
                missed_files.append(f)
            if len(missed_files) > 0:
                print('The input catalogs does not contain entries for the following images:')
                print(missed_files)
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

    print('')
    print('Finding shifts for: ')
    for f in filenames:
        print('    {}'.format(f))
    print('')

    log.info("USER INPUT PARAMETERS for finding sources for each input image:")
    util.printParams(catfile_kwargs, log=log)
    log.info('')

    try:
        minsources = max(1, catfit_pars['minobj'])
        omitted_images = []
        all_input_images = []
        for imgnum in range(len(filenames)):
            # Create Image instances for all input images
            try:
                regexcl = exclusion_dict[os.path.basename(filenames[imgnum])]
            except KeyError:
                regexcl = None
                pass

            img = imgclasses.Image(filenames[imgnum],
                                   input_catalogs=catdict[filenames[imgnum]],
                                   exclusions=regexcl,
                                   **catfile_kwargs)

            all_input_images.append(img)
            if img.num_sources < minsources:
                warn_str = "Image '{}' will not be aligned " \
                           "since it contains fewer than {} sources." \
                           .format(img.name, minsources)
                print('\nWARNING: {}\n'.format(warn_str))
                log.warning(warn_str)
                omitted_images.append(img)
                continue
            input_images.append(img)

    except KeyboardInterrupt:
        for img in input_images:
            img.close()
        print('Quitting as a result of user request (Ctrl-C)...')
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

    # input_images list can be modified below.
    # So, make a copy of the original:
    input_images_orig_copy = copy(input_images)
    do_match_refimg = False

    # otherwise, extract the catalog from the first input image source list
    if configobj['refimage'] not in [None, '',' ','INDEF']: # User specified an image to use
        # A hack to allow different source finding parameters for
        # the reference image:
        ref_sourcefind_pars = \
            tweakutils.get_configobj_root(configobj[PSET_SECTION_REFIMG])
        ref_catfile_kwargs = catfile_kwargs.copy()
        ref_catfile_kwargs.update(ref_sourcefind_pars)
        ref_catfile_kwargs['updatehdr'] = False

        log.info('')
        log.info("USER INPUT PARAMETERS for finding sources for "
                 "the reference image:")
        util.printParams(ref_catfile_kwargs, log=log)

        #refimg = imgclasses.Image(configobj['refimage'],**catfile_kwargs)
        # Check to see whether the user specified a separate catalog
        #    of reference source positions and replace default source list with it
        if refcat_par['refcat'] not in [None,'',' ','INDEF']: # User specified a catalog to use
            ref_source = refcat_par['refcat']
            cat_src = ref_source
            xycat = None
            cat_src_type = 'catalog'
        else:
            try:
                regexcl = exclusion_dict[configobj['refimage']]
            except KeyError:
                regexcl = None
                pass

            refimg = imgclasses.Image(configobj['refimage'],
                                      exclusions=regexcl,
                                      **ref_catfile_kwargs)
            ref_source = refimg.all_radec
            cat_src = None
            xycat = refimg.xy_catalog
            cat_src_type = 'image'

        try:
            if 'use_sharp_round' in ref_catfile_kwargs:
                kwargs['use_sharp_round'] = ref_catfile_kwargs['use_sharp_round']
            refimage = imgclasses.RefImage(configobj['refimage'],
                            ref_source, xycatalog=xycat,
                            cat_origin=cat_src, **kwargs)
            refwcs_fname = refimage.wcs.filename
            if cat_src is not None:
                refimage.name = cat_src

        except KeyboardInterrupt:
            refimage.close()
            for img in input_images:
                img.close()
            print('Quitting as a result of user request (Ctrl-C)...')
            return

        if len(input_images) < 1:
            warn_str = "Fewer than two images are available for alignment. " \
                       "Quitting..."
            print('\nWARNING: {}\n'.format(warn_str))
            log.warning(warn_str)
            for img in input_images:
                img.close()
            return

        image = _max_overlap_image(refimage, input_images, expand_refcat,
                                   enforce_user_order)

    elif refcat_par['refcat'] not in [None,'',' ','INDEF']:
        # a reference catalog is provided but not the reference image/wcs
        if len(input_images) < 1:
            warn_str = "No images available for alignment. Quitting..."
            print('\nWARNING: {}\n'.format(warn_str))
            log.warning(warn_str)
            for img in input_images:
                img.close()
            return

        if len(input_images) == 1:
            image = input_images.pop(0)
        else:
            image, image2 = _max_overlap_pair(input_images, expand_refcat,
                                              enforce_user_order)
            input_images.insert(0, image2)

        # Workaround the defect described in ticket:
        # http://redink.stsci.edu/trac/ssb/stsci_python/ticket/1151
        refwcs = []
        for i in all_input_images:
            refwcs.extend(i.get_wcs())
        kwargs['ref_wcs_name'] = image.get_wcs()[0].filename

        # A hack to allow different source finding parameters for
        # the reference image:
        ref_sourcefind_pars = \
            tweakutils.get_configobj_root(configobj[PSET_SECTION_REFIMG])
        ref_catfile_kwargs = catfile_kwargs.copy()
        ref_catfile_kwargs.update(ref_sourcefind_pars)
        ref_catfile_kwargs['updatehdr'] = False

        log.info('')
        log.info("USER INPUT PARAMETERS for finding sources for "
                 "the reference image (not used):")
        util.printParams(ref_catfile_kwargs, log=log)

        ref_source = refcat_par['refcat']
        cat_src = ref_source
        xycat = None

        try:
            if 'use_sharp_round' in ref_catfile_kwargs:
                kwargs['use_sharp_round'] = ref_catfile_kwargs['use_sharp_round']
            kwargs['find_bounding_polygon'] = True
            refimage = imgclasses.RefImage(refwcs,
                                           ref_source, xycatalog=xycat,
                                           cat_origin=cat_src, **kwargs)
            refwcs_fname = refimage.wcs.filename
            refimage.name = cat_src
            cat_src_type = 'catalog'
        except KeyboardInterrupt:
            refimage.close()
            for img in input_images:
                img.close()
            print('Quitting as a result of user request (Ctrl-C)...')
            return

    else:
        if len(input_images) < 2:
            warn_str = "Fewer than two images available for alignment. " \
                "Quitting..."
            print('\nWARNING: {}\n'.format(warn_str))
            log.warning(warn_str)
            for img in input_images:
                img.close()
            return

        kwargs['use_sharp_round'] = catfile_kwargs['use_sharp_round']

        cat_src = None

        refimg, image = _max_overlap_pair(input_images, expand_refcat,
                                          enforce_user_order)

        refwcs = []
        #refwcs.extend(refimg.get_wcs())
        #refwcs.extend(image.get_wcs())
        #for i in input_images:
            #refwcs.extend(i.get_wcs())
        # Workaround the defect described in ticket:
        # http://redink.stsci.edu/trac/ssb/stsci_python/ticket/1151
        for i in all_input_images:
            refwcs.extend(i.get_wcs())
        kwargs['ref_wcs_name'] = refimg.get_wcs()[0].filename

        try:
            ref_source = refimg.all_radec
            refimage = imgclasses.RefImage(refwcs, ref_source,
                                        xycatalog=refimg.xy_catalog, **kwargs)
            refwcs_fname = refimg.name
            cat_src_type = 'image'
        except KeyboardInterrupt:
            refimage.close()
            for img in input_images:
                img.close()
            print('Quitting as a result of user request (Ctrl-C)...')
            return

        omitted_images.insert(0, refimg) # refimage *must* be first
        do_match_refimg = True

    print("\n{0}\nPerforming alignment in the projection plane defined by the "
          "WCS\nderived from '{1}'\n{0}\n".format('='*63, refwcs_fname))

    if refimage.outxy is not None:
        if cat_src is None:
            cat_src = refimage.name

        try:
            log.info("USER INPUT PARAMETERS for matching sources:")
            util.printParams(objmatch_par, log=log)

            log.info('')
            log.info("USER INPUT PARAMETERS for fitting source lists:")
            util.printParams(configobj['CATALOG FITTING PARAMETERS'], log=log)

            if hdrlet_par['headerlet']:
                log.info('')
                log.info("USER INPUT PARAMETERS for creating headerlets:")
                util.printParams(hdrlet_par, log=log)

            if shiftpars['shiftfile']:
                log.info('')
                log.info("USER INPUT PARAMETERS for creating a shiftfile:")
                util.printParams(shiftpars, log=log)

            # Now, apply reference WCS to each image's sky positions as well as the
            # reference catalog sky positions,
            # then perform the fit between the reference catalog positions and
            #    each image's positions
            quit_immediately = False
            xycat_lines = ''
            xycat_filename = None

            for img in input_images_orig_copy:
                if xycat_filename is None:
                    xycat_filename = img.rootname+'_xy_catfile.list'
                # Keep a record of all the generated input_xy catalogs
                xycat_lines += img.get_xy_catnames()

            retry_flags = len(input_images)*[0]
            objmatch_par['cat_src_type'] = cat_src_type
            while image is not None:
                print ('\n'+'='*20)
                print ('Performing fit for: {}\n'.format(image.name))
                image.match(refimage, quiet_identity=False, **objmatch_par)
                assert(len(retry_flags) == len(input_images))

                if not image.goodmatch:
                    # we will try to match it again once reference catalog
                    # has expanded with new sources:
                    #if expand_refcat:
                    input_images.append(image)
                    retry_flags.append(1)
                    if len(retry_flags) > 0 and retry_flags[0] == 0:
                        retry_flags.pop(0)
                        image = input_images.pop(0)
                        # try to match next image in the list
                        continue
                    else:
                        # no new sources have been added to the reference
                        # catalog and we have already tried to match
                        # images to the existing reference catalog

                        #input_images.append(image) # <- add it back for later reporting
                        #retry_flags.append(1)
                        break

                image.performFit(**catfit_pars)
                if image.quit_immediately:
                    quit_immediately = True
                    image.close()
                    break

                # add unmatched sources to the reference catalog
                # (to expand it):
                if expand_refcat:
                    refimage.append_not_matched_sources(image)

                image.updateHeader(wcsname=uphdr_par['wcsname'],
                                    reusename=uphdr_par['reusename'])
                if hdrlet_par['headerlet']:
                    image.writeHeaderlet(**hdrlet_par)
                if configobj['clean']:
                    image.clean()
                image.close()

                if refimage.dirty and len(input_images) > 0:
                    # The reference catalog has been updated with new sources.
                    # Clear retry flags and get next image:
                    image = _max_overlap_image(
                        refimage, input_images, expand_refcat,
                        enforce_user_order
                    )
                    retry_flags = len(input_images)*[0]
                    refimage.clear_dirty_flag()
                elif len(input_images) > 0 and retry_flags[0] == 0:
                    retry_flags.pop(0)
                    image = input_images.pop(0)
                else:
                    break

            assert(len(retry_flags) == len(input_images))

            if not quit_immediately:
                # process images that have not been matched in order to
                # update their headers:
                si = 0
                if do_match_refimg:
                    image = omitted_images[0]
                    image.match(refimage, quiet_identity=True, **objmatch_par)
                    si = 1

                # process omitted (from start) images separately:
                for image in omitted_images[si:]:
                    image.match(refimage, quiet_identity=False, **objmatch_par)

                # add to the list of omitted images, images that could not
                # be matched:
                omitted_images.extend(input_images)

                if len(input_images) > 0:
                    print("\nUnable to match the following images:")
                    print("-------------------------------------")
                    for image in input_images:
                        print(image.name)
                    print("")

                # update headers:
                for image in omitted_images:
                    image.performFit(**catfit_pars)
                    if image.quit_immediately:
                        quit_immediately = True
                        image.close()
                        break
                    image.updateHeader(wcsname=uphdr_par['wcsname'],
                                    reusename=uphdr_par['reusename'])
                    if hdrlet_par['headerlet']:
                        image.writeHeaderlet(**hdrlet_par)
                    if configobj['clean']:
                        image.clean()
                    image.close()

                if configobj['writecat'] and not configobj['clean']:
                    # Write out catalog file recording input XY catalogs used
                    # This file will be suitable for use as input to 'tweakreg'
                    # as the 'catfile' parameter
                    if os.path.exists(xycat_filename):
                        os.remove(xycat_filename)
                    f = open(xycat_filename, mode='w')
                    f.writelines(xycat_lines)
                    f.close()

                    if expand_refcat:
                        base_reg_name = os.path.splitext(
                            os.path.basename(cat_src))[0]
                        refimage.write_skycatalog(
                            'cumulative_sky_refcat_{:s}.coo' \
                            .format(base_reg_name),
                            show_flux=True, show_id=True
                        )

                # write out shiftfile (if specified)
                if shiftpars['shiftfile']:
                    tweakutils.write_shiftfile(input_images_orig_copy,
                                               shiftpars['outshifts'],
                                               outwcs=shiftpars['outwcs'])

        except KeyboardInterrupt:
            refimage.close()
            for img in input_images_orig_copy:
                img.close()
                del img
            print ('Quitting as a result of user request (Ctrl-C)...')
            return
    else:
        print ('No valid sources in reference frame. Quitting...')
        return


def _overlap_matrix(images):
    nimg = len(images)
    m = np.zeros((nimg,nimg), dtype=np.float)
    for i in range(nimg):
        for j in range(i+1,nimg):
            p = images[i].skyline.intersection(images[j].skyline)
            area = np.fabs(p.area())
            m[j,i] = area
            m[i,j] = area
    return m


def _max_overlap_pair(images, expand_refcat, enforce_user_order):
    assert(len(images) > 1)
    if len(images) == 2 or not expand_refcat or enforce_user_order:
        # for the special case when only two images are provided
        # return (refimage, image) in the same order as provided in 'images'.
        # Also, when ref. catalog is static - revert to old tweakreg behavior
        im1 = images.pop(0) # reference image
        im2 = images.pop(0)
        return (im1, im2)

    m = _overlap_matrix(images)
    imgs = [f.name for f in images]
    n = m.shape[0]
    index = m.argmax()
    i = index // n
    j = index % n
    si = np.sum(m[i])
    sj = np.sum(m[:,j])

    if si < sj:
        c = j
        j = i
        i = c

    if i < j:
        j -= 1

    im1 = images.pop(i) # reference image
    im2 = images.pop(j)

    # Sort the remaining of the input list of images by overlap area
    # with the reference image (in decreasing order):
    row = m[i]
    row = np.delete(row, i)
    row = np.delete(row, j)
    sorting_indices = np.argsort(row)[::-1]
    images_arr = np.asarray(images)[sorting_indices]
    while len(images) > 0:
        del images[0]
    for k in range(images_arr.shape[0]):
        images.append(images_arr[k])

    return (im1, im2)


def _max_overlap_image(refimage, images, expand_refcat, enforce_user_order):
    nimg = len(images)
    assert(nimg > 0)
    if not expand_refcat or enforce_user_order:
        # revert to old tweakreg behavior
        return images.pop(0)

    area = np.zeros(nimg, dtype=np.float)
    for i in range(nimg):
        area[i] = np.fabs(
            refimage.skyline.intersection(images[i].skyline).area()
        )

    # Sort the remaining of the input list of images by overlap area
    # with the reference image (in decreasing order):
    sorting_indices = np.argsort(area)[::-1]
    images_arr = np.asarray(images)[sorting_indices]
    while len(images) > 0:
        del images[0]
    for k in range(images_arr.shape[0]):
        images.append(images_arr[k])

    return images.pop(0)


#
# Primary interface for running this task from Python
#
def TweakReg(files=None, editpars=False, configobj=None, imagefindcfg=None,
             refimagefindcfg=None, **input_dict):
    """
    """
    # support input of filenames from command-line without a parameter name
    # then copy this into input_dict for merging with TEAL ConfigObj parameters

    # Get default or user-specified configobj for primary task
    if isinstance(configobj, (str, bytes)):
        if configobj == 'defaults':
            # load "TEAL"-defaults (from ~/.teal/):
            configobj = teal.load(__taskname__)
        else:
            if not os.path.exists(configobj):
                raise RuntimeError('Cannot find .cfg file: '+configobj)
            configobj = teal.load(configobj, strict=False)
    elif configobj is None:
        # load 'astrodrizzle' parameter defaults as described in the docs:
        configobj = teal.load(__taskname__, defaults=True)

    if files and not util.is_blank(files):
        input_dict['input'] = files
    elif configobj is None:
        raise TypeError("TweakReg() needs either 'files' or "
                        "'configobj' arguments")

    if 'updatewcs' in input_dict: # user trying to explicitly turn on updatewcs
        configobj['updatewcs'] = input_dict['updatewcs']
        del input_dict['updatewcs']

    # Merge PSET configobj with full task configobj
    _managePsets(configobj, PSET_SECTION, imagefindpars.__taskname__,
                 iparsobj=imagefindcfg, input_dict=input_dict)
    _managePsets(configobj, PSET_SECTION_REFIMG,
                 refimagefindpars.__taskname__, iparsobj=refimagefindcfg)
    # !! NOTE - the above line needs to be done so that getDefaultConfigObj()
    # can merge in input_dict, however the TEAL GUI is not going to understand
    # the extra section, (or use it), so work needs to be done here - some
    # more thinking about how to accomplish what we want to with PSETS.
    # !! For now, warn them that imagefindpars args will be ignored in GUI.
    if editpars and input_dict:
        idkeys = input_dict.keys()
        for i in idkeys:
            if i in configobj[PSET_SECTION]:
                print('WARNING: ignoring imagefindpars setting "'+i+
                     '='+str(input_dict[i])+'", for now please enter directly into TEAL.')
                input_dict.pop(i)
        del configobj[PSET_SECTION] # force run() to pull it again after GUI use
        del configobj[PSET_SECTION_REFIMG] # force run() to pull it again after GUI use

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
        print("Problem with input parameters. Quitting...")
        return

    if configObj is None:
        return
    # If 'editpars' was set to True, util.getDefaultConfigObj() will have
    # already called 'run()'.
    if not editpars:
        # Pass full set of parameters on to the task
        run(configObj)


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


TweakReg.__doc__ = getHelpAsString(docstring = True, show_ver = False)
