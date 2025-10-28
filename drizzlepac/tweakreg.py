""" ``TweakReg`` - A replacement for IRAF-based ``tweakshifts``

:License: :doc:`/LICENSE`

"""
import os
import numpy as np
from copy import copy

from stsci.tools import teal
from stsci.tools import logutil, textutil
from stsci.tools.cfgpars import DuplicateKeyError
from stwcs import updatewcs

from . import util

# __version__ is defined here, prior to the importing
# of the modules below, so that those modules can use the value
# from this variable definition, allowing the value to be designated
# in one location only.
#

from . import tweakutils
from . import imgclasses
from . import imagefindpars
from . import refimagefindpars
from . import __version__

__taskname__ = 'tweakreg'
__all__ = ['TweakReg']

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

@util.with_logging
def run(configobj):
    """ Primary Python interface for image registration code
        This task replaces 'tweakshifts'
    """
    print('TweakReg Version %s started at: %s \n'%(
                    __version__, util._ptime()[0]))
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
        for fname in filenames:
            unique_wcsname = util.verifyUniqueWcsname(
                fname,
                wname,
                include_primary = not configobj['UPDATE HEADER']['reusename']
            )
            if not unique_wcsname:
                errstr = (f"WCSNAME '{wname}' already present in '{fname}'. "
                          "A unique value for the 'wcsname' parameter needs "
                          "to be specified.")
                print(textutil.textbox(errstr + "\n\nQuitting!", width=80))
                raise ValueError(errstr)

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
    m = np.zeros((nimg,nimg), dtype=float)
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

    area = np.zeros(nimg, dtype=float)
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
    Tweakreg provides an automated interface for computing residual shifts
    between input exposures being combined using ``AstroDrizzle``. The offsets
    computed by Tweakreg correspond to pointing differences after applying the WCS
    information from the input image's headers.  Such errors would, for example,
    be due to errors in guide-star positions when combining observations from
    different observing visits or from slight offsets introduced upon re-acquiring
    the guide stars in a slightly different position.

    Parameters
    ----------
    file : str or list of str  (Default = ``'*flt.fits'``)
        Input files (passed in from *files* parameter)
        This paramater can be provided in any of several forms:

        - filename of a single image
        
        - filename of an association (ASN)table
        
        - wild-card specification for files in directory (using ``\*``, ``?`` etc.)
        
        - comma-separated list of filenames
        
        - filelist containing list of desired input filenames with one filename 
          on each line of the file.

    editpars : bool (Default = False)
        A parameter that allows user to edit input parameters by hand in the GUI.
        ``True`` to use the GUI to edit parameters.

    configobj : ConfigObjPars, ConfigObj, dict (Default = None)
        An instance of ``stsci.tools.cfgpars.ConfigObjPars`` or
        ``stsci.tools.configobj.ConfigObj`` which overrides default parameter
        settings. When ``configobj`` is ``defaults``, default parameter values are
        loaded from the user local configuration file usually located in
        ``~/.teal/tweakreg.cfg`` or a matching configuration file in the
        current directory. This configuration file stores most recent
        settings that an user used when running ``TweakReg`` through the
        `TEAL <https://stscitools.readthedocs.io/en/latest/teal_guide.html>`_
        interface. When ``configobj`` is ``None``, ``TweakReg``
        parameters not provided explicitly will be initialized with their
        default values as described in the "Other Parameters" section.

    imagefindcfg : dict, configObject (Default = None)
        An instance of ``dict`` or ``configObject`` which overrides default source
        finding (for input images) parameter settings. See help for
        ``imagefindpars`` PSET for a list of available parameters. **Only** the
        parameters that are different from default values  **need** to be specified
        here.

    refimagefindcfg : dict, configObject (Default = None)
        An instance of ``dict`` or ``configObject`` which overrides default source
        finding (for input reference image) parameter settings. See help for
        ``refimagefindpars`` PSET for a list of available parameters. **Only** the
        parameters that are different from default values **need** to be specified
        here.

    input_dict : dict, optional
        An optional list of parameters specified by the user, which can also
        be used to override the defaults. This list of parameters **can** include
        the ``updatewcs`` parameter, even though this parameter no longer can be 
        set through the TEAL GUI. This list of parameters **can** also contain parameters
        specific to the ``TweakReg`` task itself described here in the "Other Parameters"
        section and **may not** contain parameters from the ``refimagefindpars``
        PSET. For compatibility purpose with previous ``TweakReg`` versions,
        ``input_dict`` may contain parameters from the the ``imagefindpars``
        PSET. However, if ``imagefindcfg`` is not ``None``, then ``imagefindpars``
        parameters specified through ``input_dict`` may not duplicate
        parameters specified through ``imagefindcfg``.


    Notes
    -----

    Additional parameters that can be set using the configobj, or input_dict
    arguments.

    refimage : str (Default = '')
        Filename of reference image. Sources derived from this image will be
        used as the reference for matching with sources from all input images
        unless a separate catalog is provided through the ``refcat`` parameter.
        In addition, this image file must contain a valid not distorted WCS that
        will define the projection plane in which image alignment is performed
        ("reference WCS"). When ``refimage`` is not provided, a reference WCS
        will be derived from input images.

    expand_refcat : bool (Default = False)
        Specifies whether to add new sources from just matched images to
        the reference catalog to allow next image to be matched against an
        expanded reference catalog.

    enforce_user_order : bool (Default = True)
        Specifies whether images should be aligned in the order specified in
        the ``file`` input parameter or ``TweakReg`` should optimize the order
        of alignment by intersection area of the images. Default value (`True`)
        will align images in the user specified order, except when some images
        cannot be aligned in which case ``TweakReg`` will optimize the image
        alignment order. Alignment order optimization is available *only*
        when ``expand_refcat`` = ``True``.

    exclusions: string (Default = '')
        This parameter allows the user to specify a set of files which contain
        regions in the image to ignore (or to include) when finding sources.
        This file MUST have 1 line for each input image with the name of the input
        file in the first column.  Subsequent columns would be used to specify an
        inclusion and/or exclusion file for each chip in ``'SCI,<n>'`` index order.
        If a chip does not require an exclusion file, the string ``None`` or ``INDEF``
        can be used as a placeholder for that chip. Each exclusion file can be
        either a mask provided as a simple FITS file or a region file in DS9-format.

        When a mask file is provided, ``TweakReg`` will look for the first
        image-like extension with image data of the same dimensions as the input
        image. Zeros in the mask will be interpreted as "bad"
        (excluded from search) pixels while non-zero pixels will be interpreted
        as "good" pixels. It is recommended that mask files be FITS files without
        extensions and mask data (preferably of integer type) reside in the
        primary HDU.

        If a region file is provided then it should conform to the 'region' file
        format generated by DS9. The region files can contain both regular
        ("include") regions as well as "exclude" regions. Regular ("include")
        regions indicate the regions of the image that should be searched for
        sources while "exclude" regions indicate parts of the image that should
        not be used for source detection. The "ruler", "compass", and "projection"
        regions are not supported (ignored). When **all regions**
        in a region file are "exclude" regions, then it will be assumed that the
        entire image is "good" before the exclude regions are processed. In other
        words, an "include" region corresponding to the entire image will be
        *prepended* to the list of exclude regions.

        .. note::

            Regions in a region file are processed in the order they appear in the
            region file. Thus, when region files contain *both* "include" and
            "exclude" regions, the order in which these regions appear may affect
            the results.

        .. warning::

            ``TweakReg`` relies on ``pyregion`` package for work with region files.
            At the time of writing, ``pyregion`` uses a different algorithm from DS9
            for converting regions from sky coordinates to image coordinate
            (this conversion is performed before regions are converted to masks).
            For these reasons, regions provided in sky coordinates may not produce
            the expected (from DS9) results. While in most instances these
            discrepancies should be tolerable, it is important to keep this in mind.

        During testing it was observed that conversion to image coordinates is
        most accurate for polygonal regions and less accurate for other regions
        Therefore, if one must provide regions in sky coordinates, it is
        types. recommended to use polygonal and circular regions and to avoid
        elliptical and rectangular regions as their conversion to image
        coordinates is less accurate. One may use ``mapreg`` task in the
        ``drizzlepac`` package to convert region files from sky coordinates to image
        coordinates. This will allow one to see the actual regions that will be
        used by source finding routine in ``TweakReg``.

    updatewcs : bool (Default = No)
        **NOT available through TEAL GUI interface.**
        This parameter can only be set through the Python interface to Tweakreg by
        passing it in as part of the input_dict in order to insure that running
        ``updatewcs`` **does not overwrite** a previously determined solution written out
        to the input file headers.

    writecat : bool (Default = Yes)
        Specify whether or not to write out the source catalogs generated for
        each input image by the built-in source extraction algorithm.

    clean : bool (Default = No)
        Specify whether or not to remove the temporary files created by
        ``TweakReg``, including any catalog files generated for the shift
        determination.

    interactive : bool (Default = Yes)
        This switch controls whether the program stops and waits for the user to
        examine any generated plots before continuing on to the next image.  If
        turned off, plots will still be displayed, but they will also be saved to
        disk automatically as a PNG image with an autogenerated name without
        requiring any user input.

    verbose : bool (Default = No)
        Specify whether or not to print extra messages during
        processing.

    runfile : string (Default = 'tweakreg.log')
        Specify the filename of the processing log.


    *UPDATE HEADER*

    updatehdr : bool (Default = No)
        Specify whether or not to update the headers of each input image
        directly with the shifts that were determined. This will allow the
        input images to be combined by ``AstroDrizzle`` without having to provide
        the shiftfile as well.

    wcsname : str (Default = 'TWEAK')
        Name of updated primary WCS.

    reusename : bool (Default = False)
        Allows overwriting of an existing primary WCS with the same name as
        specified by ``wcsname`` parameter.


    *HEADERLET CREATION*

    headerlet: bool (Default = No)
        Specify whether or not to generate a headerlet from the images at
        the end of the task? If turned on, this will create a headerlet
        from the images regardless of the value of the ``updatehdr``
        parameter.

    attach: bool (Default = Yes)
        If creating a headerlet, choose whether or not to attach the new
        headerlet to the input image as a new extension.


    hdrfile: string (Default = '')
        Filename to use for writing out headerlet to a separate file. If the name
        does not contain ``.fits``, it will create a filename from the
        rootname of the input image, the value of this string, and it will end
        in ``'_hlet.fits'``. For example, if only ``'hdrlet1'`` is given, the
        full filename created will be ``'j99da1f2q_hdrlet1_hlet.fits'`` when
        creating a headerlet for image ``'j99da1f2q_flt.fits'``.

    clobber: bool (Default = No)
        If a headerlet with 'hdrfile' already exists on disk, specify
        whether or not to overwrite that previous file.

    hdrname: string (Default = '')
        Unique name to give to headerlet solution. This name will be used to
        identify this specific WCS alignment solution contained in the headerlet.

    author: string, optional (Default = '')
        Name of the creator of the headerlet.

    descrip: string, optional (Default = '')
        Short (1-line) description to be included in headerlet as ``DESCRIP`` keyword.
        This can be used to provide a quick look description of the WCS alignment
        contained in the headerlet.

    catalog: string, optional (Default = '')
        Name of reference catalog used as the basis for the image alignment.

    history: string, optional (Default = '')
        Filename of a file containing detailed information regarding the history
        of the WCS solution contained in the headerlet. This can include information
        on the catalog used for the alignment, or notes on processing that went
        into finalizing the WCS alignment stored in this headerlet. This information
        will be reformatted as 70-character wide FITS HISTORY keyword section.


    *OPTIONAL SHIFTFILE OUTPUT*

    shiftfile : bool (Default = No)
        Create output shiftfile?

    outshifts : str (Default = 'shifts.txt')
        The name for the output shift file created by ``TweakReg``.  This
        shiftfile will be formatted for use as direct input to ``AstroDrizzle``.


    outwcs : str  (Default = 'shifts_wcs.fits')
        Filename to be given to the OUTPUT reference WCS file created
        by ``TweakReg``. This reference WCS defines the WCS from which the
        shifts get measured, and will be used by ``AstroDrizzle`` to interpret
        those shifts. This reference WCS file will be a FITS file that
        only contains the WCS keywords in a Primary header with no image
        data itself. The values will be derived from the FIRST input image
        specified.


    *COORDINATE FILE DESCRIPTION*

    catfile : str (Default = '')
        Name of file that contains a list of input images and associated
        catalog files generated by the user. Each line of this file will
        contain the name of an input image in the first column. The remaining
        columns will provide the names of the source catalogs for each chip
        in order of the science extension numbers ((SCI,1), (SCI,2), ...).

        A sample catfile, with one line per image would look like::

            image1_flt.fts  cat1_sci1.coo  cat1_sci2.coo
            image2_flt.fts  cat2_sci1.coo  cat2_sci2.coo

        .. note::

            Catalog files themselves must be text files containing
            "white space"-separated list of values (``xcol``, ``ycol``, etc.)

    xcol : int (Default = 1)
        Column number of X position from the user-generated catalog files
        specified in the catfile.

    ycol : int (Default = 2)
        Column number of Y position from the user-generated catalog files
        specified in the catfile.

    fluxcol : int (Default = None)
        Column number for the flux values from the user-generated catalog
        files specified in the catfile. These values will only be used if
        a flux limit has been specified by the user using the ``maxflux`` or
        ``minflux`` parameters.

    maxflux : float (Default = None)
        Limiting flux value for selecting valid objects in the input image's
        catalog.  If specified, this flux will serve as the upper limit of a
        range for selecting objects to be used in matching with objects
        identified in the reference image. If the value is set to ``None``, all
        objects with fluxes brighter than the minimum specified in ``minflux``
        will be used. If both values are set to ``None``, all objects will be used.

    minflux : float (Default = None)
        Limiting flux value for selecting valid objects in the input image's
        catalog. If specified, this flux value will serve as the lower limit
        of a range for selecting objects to be used in matching with objects
        identified in the reference image. If the value is set to ``None``, all
        objects fainter than the limit specified by ``maxflux`` will be used.
        If both values are set to ``None``, all objects will be used.

    fluxunits : str {'counts', 'cps', 'mag'} (Default = 'counts')
        This allows the task to correctly interpret the flux limits specified
        by ``maxflux`` and ``minflux`` when sorting the object list for trimming
        of fainter objects.

    xyunits : str {'pixels', 'degrees'} (Default = 'pixels')
        Specifies whether the positions in this catalog are already sky pixel
        positions, or whether they need to be transformed to the sky.

    nbright : int (Default = None)
        The number of brightest objects to keep after sorting the full object
        list. If nbright is set equal to ``None``, all objects will be used.


    *REFERENCE CATALOG DESCRIPTION*

    refcat : str (Default = '')
        Name of the external reference catalog file to be used in place of the
        catalog extracted from one of the input images. When ``refimage`` is not
        specified, reference WCS to be used with reference catalog will be
        derived from input images.

        .. note::

            Reference catalog must be text file containing
            "white space"-separated list of values (``xcol``, ``ycol``, etc.)

    refxcol : int (Default = 1)
        Column number of RA in the external catalog file specified by the
        refcat.

    refycol : int (Default = 2)
        Column number of Dec in the external catalog file specified by the
        refcat.

    refxyunits : str {'pixels','degrees'} (Default = 'degrees')
        Units of sky positions.

    rfluxcol : int (Default = None)
        Column number of flux/magnitude values in the external catalog file
        specified by the refcat.

    rmaxflux : float (Default = None)
        Limiting flux value used to select valid objects in the external
        catalog. If specified, the flux value will serve as the upper limit
        of a range for selecting objects to be used in matching with objects
        identified in the reference image. If the value is set to ``None``,
        all objects with fluxes brighter than the minimum specified in
        ``rminflux`` will be used. If both values are set to ``None``, all
        objects will be used.

    rminflux : float (Default = None)
        Limiting flux value used to select valid objects in the external
        catalog. If specified, the flux will serve as the lower limit
        of a range for selecting objects to be used in matching with objects
        identified in the reference image. If the value is set to ``None``,
        all objects fainter than the limit specified by ``rmaxflux`` will be
        used. If both values are set to ``None``, all objects will be used.

    rfluxunits : {'counts', 'cps', 'mag'} (Default = 'mag')
        This allows the task to correctly interpret the flux limits specified
        by ``rmaxflux`` and ``rminflux`` when sorting the object list for trimming
        of fainter objects.

    refnbright : int (Default = None)
        Number of brightest objects to keep after sorting the full object
        list. If refnbright is set to ``None``, all objects will be used. Used in
        conjunction with refcat.


    *OBJECT MATCHING PARAMETERS*
    

    minobj : int (Default = 15)
        Minimum number of identified objects from each input image to use
        in matching objects from other images.

    searchrad : float (Default = 1.0)
        The search radius for a match.

    searchunits : str (Default = 'arcseconds')
        Units for search radius.

    use2dhist : bool (Default = Yes)
        Use 2d histogram to find initial offset?

    see2dplot : bool (Default = Yes)
        See 2d histogram for initial offset?

    tolerance : float (Default = 1.0)
        The matching tolerance in pixels after applying an initial solution
        derived from the 'triangles' algorithm.  This parameter gets passed
        directly to ``xyxymatch`` for use in matching the object lists from each
        image with the reference image's object list.

    separation : float (Default = 0.0)
        The  minimum  separation for objects in the input and reference
        coordinate lists. Objects closer together than 'separation' pixels
        are removed from the input and reference coordinate lists prior
        to matching. This parameter gets passed directly to ``xyxymatch`` for
        use in matching the object lists from each image with the reference
        image's object list.

    xoffset : float (Default = 0.0)
        Initial estimate for the offset in X between the images and the
        reference frame. This offset will be used for all input images
        provided. If the parameter value is set to ``None``, no offset will
        be assumed in matching sources in ``xyxymatch``.

    yoffset : float (Default = 0.0)
        Initial estimate for the offset in Y between the images and the
        reference frame. This offset will be used for all input images
        provided.If the parameter value is set to None, no offset will
        be assumed in matching sources in ``xyxymatch``.


    *CATALOG FITTING PARAMETERS*

    fitgeometry : str {'shift', 'rscale', 'general'} (Default = 'rscale')
        The fitting geometry to be used in fitting the matched object lists.
        This parameter is used in fitting the offsets, rotations and/or scale
        changes from the matched object lists. The 'general' fit geometry
        allows for independent scale and rotation for each axis.

    residplot : str {'No plot', 'vector', 'residuals', 'both'}  (Default = 'both')
        Plot residuals from fit? If 'both' is selected, the 'vector'
        and 'residuals' plots will be displayed in separate plotting windows at
        the same time.

    nclip : int (Default = 3)
        Number of clipping iterations in fit.

    sigma : float (Default = 3.0)
        Clipping limit in sigma units.


    *ADVANCED PARAMETERS AVAILABLE FROM COMMAND LINE*

    updatewcs : bool  (Default = No)
        This parameter specifies whether the WCS keywords are to be updated by
        running updatewcs on the input data, or left alone. The update performed
        by updatewcs not only recomputes the WCS based on the currently
        used ``IDCTAB``, but also populates the header with
        the ``SIP`` coefficients. For ``ACS/WFC`` images, the time-dependence
        correction will also be applied to the ``WCS`` and ``SIP`` keywords.
        This parameter should be set to 'No' (`False`) when the WCS keywords
        have been carefully set by some other method, and need to be passed
        through to drizzle 'as is', otherwise those updates will be over-written
        by this update.

        .. note::

            This parameter was preserved in the API for compatibility purposes with
            existing user processing pipe-lines. However, it has been removed from
            the ``TEAL`` interface because it is easy to have it set to 'yes'
            (especially between consecutive runs of ``AstroDrizzle``) with
            potentially disastrous effects on input image WCS (for example it
            could wipe-out previously aligned WCS).

    .. note::

        Tweakreg supports the use of calibrated, distorted images (such as FLT
        images for ACS and WFC3, or ``_c0m.fits`` images for WFPC2) as input images.
        All coordinates for sources derived from these images (either by this task
        or as provided by the user directly) will be corrected for distortion using
        the distortion model information specified in each image's header. This
        eliminates the need to run ``AstroDrizzle`` on the input images prior to
        running ``TweakReg``.

    .. note::

        All calibrated input images must have been updated using
        ``updatewcs`` from the ``STWCS`` package, to include the full
        distortion model in the header. Alternatively, one can set
        ``updatewcs`` parameter to ``True`` when running either ``TweakReg``
        or ``AstroDrizzle`` from command line (Python interpreter)
        **the first time** on such images.

    This task will use catalogs, and catalog-matching, based on the ``xyxymatch``
    algorithm to determine the offset between the input images. The primary
    mode of operation will be to extract a catalog of source positions from
    each input image using either a 'DAOFIND-like' algorithm or ``SExtractor`` (if
    the user has ``SExtractor`` installed). Alternatively, the user can provide
    their catalogs of source positions derived from **each input chip**.

    .. note::

        Catalog files must be text files containing
        "white space"-separated list of values (``xcol``, ``ycol``, etc.)

    The reference frame will be defined either by:

        * the image with the largest overlap with another input image AND with 
          the largest total overlap with the rest of the input images,
          
        * a catalog derived from a reference image specified by the user, or
        
        * a catalog of undistorted sky positions (RA/Dec) and fluxes provided by 
          the user.

    For a given observation, the distortion model is applied to all distorted
    input positions, and the sources from each chip are then combined into a
    single catalog of undistorted positions.

    The undistorted positions for each observation then get passed to
    ``xyxymatch`` for matching to objects from the reference catalog.

    The source lists from each image will generally include cosmic-rays as
    detected sources, which can at times significantly confuse object
    identification between images. Observations that include long exposures
    often have more cosmic-ray events than source objects. As such, isolating
    the cosmic-ray events in those cases would significantly improve the
    efficiency of common source identification between images. One such method
    for trimming potential false detections from each source list would be
    to set a flux limit to exclude detections below that limit. As the fluxes
    reported in the default source object lists are provided as magnitude
    values, setting the ``maxflux`` or ``minflux`` parameter value to a magnitude-
    based limit, and then setting the ``ascend`` parameter to ``True``, will allow
    for the creations of catalogs trimmed of all sources fainter than the
    provided limit. The trimmed source list can then be used in matching
    sources between images and in establishing the final fitting for the shifts.

    A fit can then be performed on the matched set of positions between the
    input and the reference to produce the 'shiftfile'. If the user is
    confident that the solution will be correct, the header of each input image
    can be updated directly with the fit derived for that image. Otherwise,
    the 'shiftfile' can be passed to AstroDrizzle for aligning the images.

    .. note::

        Because of the nature of the used algorithm it may be necessary
        to run this task multiple time until new shifts, rotations,
        and/or scales are small enough for the required precision.

    New sources (that are not in the reference catalog) from the matched images
    are added to the reference catalog in order to allow next image to be
    matched to a larger reference catalog. This allows alignment of images
    that do not overlap directly with the reference image and/or catalog and
    it is particularly useful in image registration of large mosaics. Addition
    of new sources to the reference catalog can be turned off by
    setting ``expand_refcat`` to ``False`` when using an external reference catalog.
    When an external catalog is not provided (``refcat``='') or when
    using an external reference catalog with ``expand_refcat`` set to ``True``
    (assuming ``writecat`` = ``True`` and ``clean`` = ``False``),
    the list of all sources in the expanded reference catalog is saved
    in a catalog file named ``cumulative_sky_refcat_###.coo`` where ### is the
    base file name derived from either the external catalog (if provided) or
    the name of the image used as the reference image.

    When ``enforce_user_order`` is ``False``, image catalogs are matched to the
    reference catalog in order of decreasing overlap area with the reference
    catalog, otherwise user order of files specified in the ``file`` parameter
    is used.


    **Format of Exclusion Catalog**

    The format for the exclusions catalog requires 1 line in the file for
    every input image, regardless of whether or not that image has
    any defined exclusion regions.  A sample file would look like::

        j99da1emq_flt.fits
        j99da1f2q_flt.fits test_exclusion.reg

    This file specifies no exclusion files for the first image, and only
    an regions file for SCI,1 of the second image.  NOTE: The first file can be
    dropped completely from the exclusion catalog file.

    In the above example, should an exclusion regions file only be needed for the
    second chip in the second image, the file would need to look like::

        j99da1emq_flt.fits
        j99da1f2q_flt.fits None test_sci2_exclusion.reg

    The value ``None`` could also be replaced by ``INDEF`` if desired, but either
    string needs to be present to signify no regions file for that chip while
    the code continues parsing the line to find a file for the second chip.

    
    **Format of Region Files**

    The format of the exclusions catalogs referenced in the 'exclusions'
    file defaults to the format written out by DS9 using the 'DS9/Funtools'
    region file format.  A sample file with circle() regions will look like::

        # Region file format: DS9 version 4.1
        # Filename: j99da1f2q_flt.fits[SCI]
        global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman"
        select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
        image
        circle(3170,198,20)
        ellipse(3269,428,30,10,45) # a rotated ellipse
        box(3241.1146,219.78132,20,20,15) # a rotated box
        circle(200,200,50)  # outer circle
        -circle(200,200,30) # inner circle

    This region file will be interpreted as "find all sources in the image that
    **are inside** the four regions above but **not inside** the
    region -circle(200,200,30)". Effectively we will instruct ``TweakReg``
    to find all the sources *inside* the following regions::

        circle(3170,198,20)
        ellipse(3269,428,30,10,45) # a rotated ellipse
        box(3241.1146,219.78132,20,20,15) # a rotated box
        annulus(200,200,30,50)  # outer circle(r=50) - inner circle(r=30)

    Examples
    --------
    The tweakreg task can be run from either the TEAL GUI or from the command-line
    using Python.
    These examples illustrate the various syntax options available.

    **Example 1:**  Align a set of calibrated (``_flt.fits``) images
    using ``IMAGEFIND``, a built-in source
    finding algorithm based on ``DAOPHOT``. Auto-detect the sky sigma value and select sources > 200
    sigma.   (Auto-sigma is computed from the first input exposure
    as:  ``1.5*imstat(image,nclip=3,fields='stddev')``. )
    Set the convolution kernel width to ``~2x`` the value of the PSF FWHM.
    Save the residual offsets (``dx``, ``dy``, ``rot``, ``scale``, ``xfit_rms``, ``yfit_rms``) to a text file.

     1. Run the task from Python using the command line while individually
         specifying source finding parameters for the reference image and
         input images:

        >>> import drizzlepac
        >>> from drizzlepac import tweakreg
        >>> tweakreg.TweakReg('*flt.fits',
        ...       imagefindcfg={'threshold' : 200, 'conv_width' : 3.5},
        ...       refimagefindcfg={'threshold' : 400, 'conv_width' : 2.5},
        ...       updatehdr=False, shiftfile=True, outshifts='shift.txt')

         or, using ``dict`` constructor,

        >>> import drizzlepac
        >>> from drizzlepac import tweakreg
        >>> tweakreg.TweakReg('*flt.fits',
        ...       imagefindcfg=dict(threshold=200, conv_width=3.5),
        ...       refimagefindcfg=dict(threshold=400, conv_width=2.5),
        ...       updatehdr=False, shiftfile=True, outshifts='shift.txt')

         Or, run the same task from the Python command line, but specify all parameters in
         a config file named "myparam.cfg":

        >>> tweakreg.TweakReg('*flt.fits', configobj='myparam.cfg')

         Alternately, edit the imagefind parameters in a TEAL GUI window
         prior to running the task:

        >>> tweakreg.edit_imagefindpars()

     2. Help can be accessed via the "Help" pulldown menu in the TEAL GUI.  It can also
         be accessed from the Python command-line and saved to a text file:

        >>> from drizzlepac import tweakreg
        >>> tweakreg.help()

         or

        >>> tweakreg.help(file='help.txt')

    See Also
    --------
    drizzlepac.astrodrizzle

    """
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
