"""
``tweakback`` - propagate the "tweaked" solutions back to the original
input files.

Version 0.4.0 - replaced previous algorithm that used fitting of WCS
footprints to reconstruct the transformation that was applied to the
old drizzled image (to align it with another image) to obtain the new
drizzled image WCS with an algorithm that is based on linearization of
the exact compound operator that transforms current image coordinates to
the "aligned" (to the new drizzled WCS) image coordinates.

:Authors: Warren Hack, Mihai Cara

:License: :doc:`/LICENSE`

"""
import os
import string
from datetime import date

import numpy as np
from astropy.io import fits
from astropy.utils.decorators import deprecated

from stwcs import wcsutil
from stwcs.wcsutil import altwcs
from stsci.tools import parseinput, logutil
from stsci.skypac.utils import get_ext_list, ext2str
from stsci.skypac.parseat import parse_cs_line


from . import updatehdr
from . import linearfit
from . import util
from . import __version__


__taskname__ = 'tweakback'


log = logutil.create_logger(__name__, level=logutil.logging.NOTSET)

if hasattr(np, 'float128'):
    ndfloat128 = np.float128
elif hasattr(np, 'float96'):
    ndfloat128 = np.float96
else:
    ndfloat128 = np.float64


def _wcs_key_name(wcs_key, wcs_name, fi, ext, default_key, param_name):
    wnames = altwcs.wcsnames(fi.image.hdu, ext=ext)
    if wcs_key is None:
        if wcs_name is None:
            wcs_key = sorted(wnames)[default_key]
            wcs_name = wnames[wcs_key]
            # report our guess:
            print(f"Using WCS with WCSNAME '{wcs_name}' (WCS key '{wcs_key}') "
                  f"for '{param_name}'")
        else:
            # find associated wcs_key
            wcs_name_u = wcs_name.upper()
            for k, name in wnames.items():
                if name.upper() == wcs_name_u:
                    wcs_key = k
                    break
            else:
                fi.release_all_images()
                raise KeyError(f"WCS with wcsname '{wcs_name}' not found.")
    else:
        wcs_name = wnames[wcs_key]

    return wcs_key, wcs_name


def _process_wcs_key_par(par_name, kwargs):
    if par_name in kwargs:
        wcs_key = kwargs[par_name]
        if len(wcs_key) != 1 or wcs_key.upper().strip() not in string.ascii_uppercase:
            raise ValueError(
                f"Parameter '{par_name}' must be a character - "
                "one of 'A'-'Z' or ' '."
            )
    else:
        wcs_key = None
    return wcs_key


def apply_tweak(drz_file, orig_wcs_name, output_wcs_name=None, input_files=None,
                default_extname='SCI', **kwargs):
    """
    Apply WCS solution recorded in drizzled file to distorted input images
    (``_flt.fits`` files) used to create the drizzled file.

    It is assumed that if input images given by ``input_files`` are drizzled
    together, they would produce the drizzled image given by ``drz_file`` image
    and with the "original" primary WCS. It is also assumed that drizzled image
    was aligned using ``tweakreg`` either to another image or to an external
    reference catalog. We will refer to the primary WCS in the drizzled image
    **before** ``tweakreg`` was run as the "original" WCS and the WCS **after**
    ``tweakreg`` was run as "tweaked" WCS.

    By comparing both "original" and "tweaked" WCS, ``apply_wcs`` computes
    the correction that was applied by ``tweakreg`` to the "original" WCS
    and converts this correction in the drizzled image frame into a correction
    in the input image's (``input_files``) frame that will be applied to the
    primary WCS of input images. If updated input images are now resampled
    again, they would produce an image very close to ``drz_file`` but with
    a primary WCS very similar to the "tweaked" WCS instead of the "original"
    WCS.

    Parameters
    ----------
    drz_file : str
        File name of the drizzled image that contains both the "original" and
        "tweaked" WCS. Even though wildcards are allowed in the file name,
        their expansion must resolve to a single image file. By default,
        ``apply_tweak`` looks for the first image-like HDU in the drizzled
        image. To specify a particular extension from which to load WCS,
        append extension specification after the file name, for example:

            - ``'image_drz.fits[sci,1]'`` for first "sci" extension
            - ``'image_drz.fits[1]'`` for the first extension
            - ``'image_drz.fits[0]'`` for the primary HDU

    orig_wcs_name : str
        Name (provided by the ``WCSNAME?`` header keyword where ``?``
        respesents a letter A-Z) of the "original" WCS. This is the WCS of
        the resampled image (obtained by drizzling all input images)  **before**
        this resampled image was aligned ("tweaked") to another image/catalog.

        If ``orig_wcs_name`` is ``None``, the the original WCS **must be
        specified** using ``orig_wcs_key``. When ``orig_wcs_key`` is provided,
        ``orig_wcs_name`` is ignored altogether.

    output_wcs_name : str, None
        Value of ``WCSNAME`` to be used to label the updated solution in the
        input (e.g., ``_flt.fits``) files.  If left blank or ``None``, it will
        default to using either the current (primary) ``WCSNAME`` value from
        the ``drz_file`` or from the alternate WCS given by the
        ``tweaked_wcs_name`` or ``tweaked_wcs_key`` parameters.

    input_files : str, None
        Filenames of distorted images whose primary WCS is to be updated
        with the same transformation as used in the "tweaked" drizzled image.
        Default value of ``None`` indicates that input image filenames will be
        derived from the ``D*DATA`` keywords written out by the ``AstroDrizzle``.
        If they can not be found, the task will quit.

        ``input_files`` string can contain one of the following:

            * a comma-separated list of valid science image file names
              (see note below) and (optionally) extension specifications,
              e.g.: ``'j1234567q_flt.fits[1], j1234568q_flt.fits[sci,2]'``;

            * an @-file name, e.g., ``'@files_to_match.txt'``.

        .. note::
            **Valid** **science** **image** **file** **names** are:

            * file names of existing FITS, GEIS, or WAIVER FITS files;

            * partial file names containing wildcard characters, e.g.,
              ``'*_flt.fits'``;

            * Association (ASN) tables (must have ``_asn``, or ``_asc``
              suffix), e.g., ``'j12345670_asn.fits'``.

        .. warning::
            @-file names **MAY** **NOT** be followed by an extension
            specification.

        .. warning::
            If an association table or a partial file name with wildcard
            characters is followed by an extension specification, it will be
            considered that this extension specification applies to **each**
            file name in the association table or **each** file name
            obtained after wildcard expansion of the partial file name.

    default_extname : str
        Extension name of extensions in input images whose primary WCS
        should be updated. This value is used only when file names provided in
        ``input_files`` do not contain extension specifications.

    Other Parameters
    ----------------
    tweaked_wcs_name : str
        Name of the "tweaked" WCS. This is the WCS of
        the resampled image (obtained by drizzling all input images)  _after_
        this resampled image was aligned ("tweaked") to another image/catalog.

        When neither ``tweaked_wcs_name`` nor ````tweaked_wcs_key`` are not
        provided, ``apply_tweak`` will take the current primary WCS in the
        drizzled image as a "tweaked" WCS. ``tweaked_wcs_name`` is ignored
        when ``tweaked_wcs_key`` is provided.

    tweaked_wcs_key : {' ', 'A'-'Z'}
        Same as ``tweaked_wcs_name`` except it specifies a WCS by key instead
        of name. When provided, ``tweaked_wcs_name`` is ignored.

    orig_wcs_key : {' ', 'A'-'Z'}
        Same as ``orig_wcs_name`` except it specifies a WCS by key instead
        of name. When provided, ``orig_wcs_name`` is ignored.

    Notes
    -----
    The algorithm used by this function is based on linearization of
    the exact compound operator that converts input image coordinates
    to the coordinates (in the input image) that would result in
    alignment with the new drizzled image WCS.

    .. warning::
        Parameters ``orig_wcs_name`` and ``tweaked_wcs_name`` (or their "key"
        equivalents) allow computation of transformation between *any two
        WCS* in the drizzled image and application of this transformation to the
        primary WCS of the input images. This will produce an
        expected result **only if** the WCS pointed to by ``orig_wcs_name`` was
        obtained by drizzling input images with their current primary WCS.


    EXAMPLES
    --------
    A drizzled image named ``acswfc_mos2_drz.fits`` was created from 4 images
    using ``AstroDrizzle``. The primary WCS of this drizzled image was named
    ``'INITIAL_GUESS'``. This drizzled image was then aligned to some other
    image using ``TweakReg`` and the updated ("tweaked") primary WCS was named
    ``'BEST_WCS'`` while the previous primary WCS - the WCS named
    ``'INITIAL_GUESS'`` - was archived by ``TweakReg`` under WCS key ``'C'``.
    We will refer to this archived WCS as the "original" WCS.
    ``apply_tweak`` can now be used to compute the
    transformation between the original and the tweaked WCS and apply this
    transformation to the WCS of each of the input images that were
    drizzle-combined to produce the resampled image ``acswfc_mos2_drz.fits``.

    The simplest way to accomplish this would be to run ``apply_tweak()`` using
    default parameters:

    >>> from drizzlepac import tweakback
    >>> tweakback.apply_tweak('acswfc_mos2_drz.fits', orig_wcs_name='INITIAL_GUESS')

    or

    >>> tweakback.apply_tweak('acswfc_mos2_drz.fits', orig_wcs_key='C')

    If the same WCS should be applied to a specific set of images or
    extensions in those images, then we can explicitly specify input files:

    >>> tweakback.apply_tweak(
    ...     'acswfc_mos2_drz.fits',
    ...     input='img_mos2a_flt.fits,img_mos2c_flt.fits[1],img_mos2d_flt.fits[sci,1]'
    ... )

    In the examples above, current primary WCS of the input
    ``'img_mos2?_flt.fits'`` files will be archived and the primary WCS will
    be replaced with a "tweaked" WCS obtained by applying relevant
    transformations to the current primary WCS. Because we did not specify
    ``output_wcs_name``, the name of this tweaked primary WCS in the
    input images will be set to ``'BEST_WCS'``.

    See Also
    --------
    stwcs.wcsutil.altwcs

    """
    print(f"\n*** 'apply_tweak' version {__version__:s} started "
          f"at {util._ptime()[0]:s}: ***\n")

    tweaked_wcs_name = kwargs.get('tweaked_wcs_name', None)
    tweaked_wcs_key = kwargs.get('tweaked_wcs_key', None)
    orig_wcs_key = kwargs.get('orig_wcs_key', None)

    tweaked_wcs_key = _process_wcs_key_par('tweaked_wcs_key', kwargs)
    orig_wcs_key = _process_wcs_key_par('orig_wcs_key', kwargs)

    # load drizzled image and extract input file names (if needed) and
    # load specified WCS:

    fis = parse_cs_line(
        drz_file, default_ext='*', fnamesOnly=False, doNotOpenDQ=True,
        im_fmode="readonly"
    )
    if len(fis) == 0:
        raise FileNotFoundError(f"Drizzled file '{drz_file}' not found.")
    elif len(fis) > 1:
        for f in fis:
            f.release_all_images()
        raise ValueError("When expanded, 'drz_file' should correspond to a "
                         "single file.")

    fi = fis[0]
    hdul = fi.image.hdu
    if len(fi.fext) == 1:
        drz_sciext = fi.fext[0]

    elif fi.fext:
        fi.release_all_images()
        raise ValueError(
            "Input drizzled image contains multiple image-like extensions. "
            "Please explicitly specify a single extension containing desired "
            "WCS."
        )

    else:
        fi.release_all_images()
        raise ValueError(
            "Specified extension was not found in the input drizzled image."
        )

    # check that there are at least two WCS in the drizzled image header:
    wkeys = altwcs.wcskeys(hdul, ext=drz_sciext)
    if len(wkeys) < 2:
        fi.release_all_images()
        raise ValueError(f"'{fi.image}[{ext2str(drz_sciext)}]' must "
                         "contain at least two valid WCS: original and "
                         "updated by tweakreg.")

    # load "tweaked" WCS
    tweaked_wcs_key, tweaked_wcs_name = _wcs_key_name(
        tweaked_wcs_key,
        tweaked_wcs_name,
        fi=fi,
        ext=drz_sciext,
        default_key=0,
        param_name='tweaked_wcs_name'
    )
    tweaked_wcs = wcsutil.HSTWCS(hdul, ext=drz_sciext, wcskey=tweaked_wcs_key)

    # load "original" WCS
    if orig_wcs_key is None and orig_wcs_name is None:
        fi.release_all_images()
        raise ValueError("Either 'orig_wcs_name' or 'orig_wcs_key' must be specified.")

    # default_key=-1 below is useless since we require that either
    # orig_wcs_key or orig_wcs_name be specified. However, in the future,
    # if we allow both to be None, we can use the last WCS key as the default
    # for the "original" WCS (last WCS key in the list).
    orig_wcs_key, orig_wcs_name = _wcs_key_name(
        orig_wcs_key,
        orig_wcs_name,
        fi=fi,
        ext=drz_sciext,
        default_key=-1,
        param_name='orig_wcs_name'
    )
    orig_wcs = wcsutil.HSTWCS(hdul, ext=drz_sciext, wcskey=orig_wcs_key)

    # get RMS values reported for new solution
    crderr1 = fi.image.hdu[drz_sciext].header.get('CRDER1' + orig_wcs_key, 0.0)
    crderr2 = fi.image.hdu[drz_sciext].header.get('CRDER2' + orig_wcs_key, 0.0)

    fi.release_all_images()  # done with the resampled image

    # Process the list of input files:
    if not isinstance(default_extname, str):
        raise TypeError("Argument 'default_extname' must be a string")
    default_extname = default_extname.strip()
    if default_extname.upper() == 'PRIMARY':
        ext2get = ('PRIMARY', 1)
    else:
        ext2get = (default_extname, '*')

    if input_files is None:
        # get input (FLT) file names from the drizzled image. This information
        # is recorded in the primary header of the drizzled image.
        input_files = ",".join(hdul[0].header["D???DATA"].values())

    # Build a list of input files and extensions
    fnames_ext = {}
    fis = parse_cs_line(
        input_files, default_ext=ext2get, fnamesOnly=False,
        doNotOpenDQ=True, im_fmode="readonly"
    )
    for f in fis:
        f.release_all_images()
        if f.image in fnames_ext:
            fnames_ext[f.image] |= set(f.fext)
        else:
            fnames_ext[f.image] = set(f.fext)

    if output_wcs_name is None:
        output_wcs_name = tweaked_wcs_name
        print(f"\n* Setting 'output_wcs_name' to '{output_wcs_name}'")
        auto_output_name = True
    else:
        auto_output_name = False

    output_wcs_name_u = output_wcs_name.strip().upper()

    # Compute tweakback transformation to each extension of each input file.
    # This is the main part of this function.
    # Before applying new WCS solution, make sure we can use the same
    # output WCS name for the updated WCS in all input images.
    # Also, this gives us opportunity to remove duplicate extensions, if any.

    final_wcs_info = []

    for fname, extlist in fnames_ext.items():
        print(f"\n* Working on input image {fname:s} ...")

        fis = parse_cs_line(
            f"{fname}",
            default_ext=ext2get, fnamesOnly=False,
            doNotOpenDQ=True, im_fmode="readonly"
        )
        if len(fis) != 1:
            fi.release_all_images()
            raise AssertionError("The algorithm should not open more than one file.")
        fi = fis[0]

        if not fi.fext:
            fi.release_all_images()
            print(f"  No valid input image extension found. Skipping image {fname}.\n")
            continue

        current_wcs_info = {
            'fname': fname,
            'extlist': [],
            'archived_wcs_name': [],
            'updated_primary_wcs': []
        }
        final_wcs_info.append(current_wcs_info)

        # Process extensions
        hdu_list = []  # to avoid processing duplicate hdus
        try:
            for ext in extlist:
                imhdulist = fi.image.hdu
                hdu = imhdulist[ext]
                if hdu in hdu_list:
                    continue
                hdu_list.append(hdu)

                current_wcs_info['extlist'].append(ext)

                # Find the name under which to archive current WCS:
                all_wcs_names = [
                    v.upper() for v in
                    altwcs.wcsnames(imhdulist, ext, include_primary=False).values()
                ]

                if output_wcs_name_u in all_wcs_names:
                    if auto_output_name:
                        raise ValueError(
                            "Current value of 'output_wcs_name' was set to "
                            f"'{tweaked_wcs_name}' by default. However, this "
                            f"WCS name value was already used in {fname:s}[{ext2str(ext)}]. "
                            "Please re-run 'apply_tweak' again and explicitly "
                            "provide a unique value for the output WCS name."
                        )
                    else:
                        raise ValueError(
                            "Provided value of 'output_wcs_name' - '{tweaked_wcs_name}' - "
                            f"was already used in {fname:s}[{ext2str(ext)}]. "
                            "Please re-run 'apply_tweak' again and explicitly "
                            "provide a unique value for the output WCS name."
                        )

                if 'WCSNAME' in imhdulist[ext].header:
                    pri_wcs_name = imhdulist[ext].header['WCSNAME'].strip()
                else:
                    pri_wcs_name = 'NONAME'

                # add current output WCS name to the list so that archived
                # primary WCS will be archived under a different name:
                all_wcs_names.append(output_wcs_name)

                archived_name = altwcs._auto_increment_wcsname(pri_wcs_name, all_wcs_names)
                current_wcs_info['archived_wcs_name'].append(archived_name)

                # compute updated WCS:
                new_wcs = wcsutil.HSTWCS(imhdulist, ext=ext)

                update_chip_wcs(new_wcs, orig_wcs, tweaked_wcs,
                                xrms=crderr1, yrms=crderr2)
                new_wcs.setOrient()
                current_wcs_info['updated_primary_wcs'].append(new_wcs)

                print(f"  - Computed new WCS solution for {fname:s}[{ext2str(ext)}]:")
                repr_wcs = repr(new_wcs)
                print('\n'.join(['      ' + l.strip() for l in repr_wcs.split('\n')]))

        finally:
            fi.release_all_images()

    print("\n* Saving updated WCS to image headers:")

    for fwi in final_wcs_info:
        if not fwi['extlist']:
            continue

        fname = fwi['fname']

        fis = parse_cs_line(
            f"{fname}",
            default_ext=ext2get, fnamesOnly=False,
            doNotOpenDQ=True, im_fmode="update"
        )
        fi = fis[0]

        # Process extensions
        try:
            for ext, archived_name, new_wcs in zip(fwi['extlist'],
                                                   fwi['archived_wcs_name'],
                                                   fwi['updated_primary_wcs']):
                imhdulist = fi.image.hdu
                hdu = imhdulist[ext]

                hdu.header['HISTORY'] = (
                    f"apply_tweak version: {__version__} ({date.today().isoformat():s})"
                )

                # Archive current primary WCS:
                awcs_key, awcs_name = altwcs.archive_wcs(
                    imhdulist, ext, wcsname=archived_name,
                    mode=altwcs.ArchiveMode.NO_CONFLICT
                )
                hdu.header['HISTORY'] = (
                    "apply_tweak: Archived Primary WCS under key "
                    f"'{awcs_key}' on {date.today().isoformat():s}"
                )
                hdu.header['HISTORY'] = (
                    f"apply_tweak: WCSNAME{awcs_key}='{awcs_name}'"
                )

                # Update primary WCS of this extension:
                wcs_hdr = new_wcs.wcs2header(idc2hdr=new_wcs.idcscale is not None, relax=True)
                wcs_hdr.set('WCSNAME', output_wcs_name, before=0)
                wcs_hdr.set(
                    'WCSTYPE',
                    updatehdr.interpret_wcsname_type(output_wcs_name),
                    after=0
                )
                wcs_hdr.set('ORIENTAT', new_wcs.orientat, after=len(wcs_hdr))
                hdu.header.update(wcs_hdr)
                hdu.header['HISTORY'] = (
                    f"apply_tweak: Applied Primary WCS correction on {date.today().isoformat():s}"
                )

            str_extlist = '; '.join(map(ext2str, fwi['extlist']))
            print(f"  - Updated '{fname:s}', extensions: {str_extlist}")

        finally:
            util.updateNEXTENDKw(imhdulist)
            fi.release_all_images()


@deprecated(since='3.4.2', name='tweakback', alternative='apply_tweak')
def tweakback(drzfile, input=None,  origwcs = None,
              newname = None, wcsname = None,
              extname='SCI', force=False, verbose=False):
    """
    Apply WCS solution recorded in drizzled file to distorted input images
    (``_flt.fits`` files) used to create the drizzled file.  This task relies on
    the original WCS and updated WCS to be recorded in the drizzled image's
    header as the last 2 alternate WCSs.

    Parameters
    ----------
    drzfile : str (Default = '')
        filename of undistorted image which contains the new WCS
        and WCS prior to being updated

    newname : str (Default = None)
        Value of ``WCSNAME`` to be used to label the updated solution in the
        output (e.g., ``_flt.fits``) files.  If left blank or ``None``, it will
        default to using the current ``WCSNAME`` value from the input drzfile.

    input : str (Default = '')
        filenames of distorted images to be updated using new WCS
        from 'drzfile'.  These can be provided either as an ``@-file``,
        a comma-separated list of filenames or using wildcards.

        .. note:: A blank value will indicate that the task should derive the
           filenames from the 'drzfile' itself, if possible. The filenames will be
           derived from the ``D*DATA`` keywords written out by
           ``AstroDrizzle``. If they can not be found, the task will quit.

    origwcs : str (Default = None)
        Value of ``WCSNAME`` keyword prior to the drzfile image being updated
        by ``TweakReg``.  If left blank or None, it will default to using the
        second to last ``WCSNAME*`` keyword value found in the header.

    wcsname : str (Default = None)
        Value of WCSNAME for updated solution written out by ``TweakReg`` as
        specified by the ``wcsname`` parameter from ``TweakReg``.  If this is
        left blank or ``None``, it will default to the current ``WCSNAME``
        value from the input drzfile.

    extname : str (Default = 'SCI')
        Name of extension in ``input`` files to be updated with new WCS

    force : bool  (Default = False)
        This parameters specified whether or not to force an update of the WCS
        even though WCS already exists with this solution or ``wcsname``?

    verbose : bool (Default = False)
        This parameter specifies whether or not to print out additional
        messages during processing.


    Notes
    -----
    The algorithm used by this function is based on linearization of
    the exact compound operator that converts input image coordinates
    to the coordinates (in the input image) that would result in
    alignment with the new drizzled image WCS.

    If no input distorted files are specified as input, this task will attempt
    to generate the list of filenames from the drizzled input file's own
    header.


    EXAMPLES
    --------
    An image named ``acswfc_mos2_drz.fits`` was created from 4 images using
    astrodrizzle. This drizzled image was then aligned to another image using
    tweakreg and the header was updated using the ``WCSNAME`` = ``TWEAK_DRZ``.
    The new WCS can then be used to update each of the 4 images that were
    combined to make up this drizzled image using:

    >>> from drizzlepac import tweakback
    >>> tweakback.tweakback('acswfc_mos2_drz.fits')

    If the same WCS should be applied to a specific set of images, those images
    can be updated using:

    >>> tweakback.tweakback('acswfc_mos2_drz.fits',
    ...                     input='img_mos2a_flt.fits,img_mos2e_flt.fits')

    See Also
    --------
    stwcs.wcsutil.altwcs

    """
    print("TweakBack Version {:s} started at: {:s}\n"
          .format(__version__, util._ptime()[0]))

    # Interpret input list/string into list of filename(s)
    fltfiles = parseinput.parseinput(input)[0]

    if fltfiles is None or len(fltfiles) == 0:
        # try to extract the filenames from the drizzled file's header
        fltfiles = extract_input_filenames(drzfile)
        if fltfiles is None:
            print('*'*60)
            print('*')
            print('* ERROR:')
            print('*    No input filenames found! ')
            print('*    Please specify "fltfiles" or insure that input drizzled')
            print('*    image contains D*DATA keywords. ')
            print('*')
            print('*'*60)
            raise ValueError

    if not isinstance(fltfiles,list):
        fltfiles = [fltfiles]

    sciext = determine_extnum(drzfile, extname='SCI')
    scihdr = fits.getheader(drzfile, ext=sciext, memmap=False)

    ### Step 1: Read in updated and original WCS solutions
    # determine keys for all alternate WCS solutions in drizzled image header
    wkeys = wcsutil.altwcs.wcskeys(drzfile, ext=sciext)
    if len(wkeys) < 2:
        raise ValueError(f"'{drzfile}' must contain at least two valid WCS: original and updated.")
    wnames = wcsutil.altwcs.wcsnames(drzfile, ext=sciext)
    if not util.is_blank(newname):
        final_name = newname
    else:
        final_name = wnames[wkeys[-1]]

    # Read in HSTWCS objects for final,updated WCS and previous WCS from
    # from drizzled image header
    # The final solution also serves as reference WCS when using updatehdr

    if not util.is_blank(wcsname):
        for wkey, wname in wnames.items():
            if wname == wcsname:
                wcskey = wkey
                break
        else:
            raise ValueError(f"WCS with name '{wcsname}' not found in '{drzfile}'")
    else:
        wcskey = wkeys[-1]

    final_wcs = wcsutil.HSTWCS(drzfile, ext=sciext, wcskey=wcskey)

    if not util.is_blank(origwcs):
        for wkey, wname in wnames.items():
            if wname == origwcs:
                orig_wcskey = wkey
                break
        else:
            raise ValueError(f"WCS with name '{origwcs}' not found in '{drzfile}'")
    else:
        _, orig_wcskey = determine_orig_wcsname(scihdr, wnames, wkeys)

    orig_wcs = wcsutil.HSTWCS(drzfile, ext=sciext, wcskey=orig_wcskey)

    # read in RMS values reported for new solution
    crderr1kw = 'CRDER1'+wkeys[-1]
    crderr2kw = 'CRDER2'+wkeys[-1]

    if crderr1kw in scihdr:
        crderr1 = fits.getval(drzfile, crderr1kw, ext=sciext, memmap=False)
    else:
        crderr1 = 0.0

    if crderr2kw in scihdr:
        crderr2 = fits.getval(drzfile, crderr2kw, ext=sciext, memmap=False)
    else:
        crderr2 = 0.0
    del scihdr

    ### Step 2: Apply solution to input file headers
    for fname in fltfiles:
        logstr = "....Updating header for {:s}...".format(fname)
        if verbose:
            print("\n{:s}\n".format(logstr))
        else:
            log.info(logstr)

        # reset header WCS keywords to original (OPUS generated) values
        imhdulist = fits.open(fname, mode='update', memmap=False)
        extlist = get_ext_list(imhdulist, extname='SCI')
        if not extlist:
            extlist = [0]

        # Process MEF images...
        for ext in extlist:
            logstr = "Processing {:s}[{:s}]".format(imhdulist.filename(),
                                                    ext2str(ext))
            if verbose:
                print("\n{:s}\n".format(logstr))
            else:
                log.info(logstr)
            chip_wcs = wcsutil.HSTWCS(imhdulist, ext=ext)

            update_chip_wcs(chip_wcs, orig_wcs, final_wcs,
                            xrms=crderr1, yrms = crderr2)

            # Update FITS file with newly updated WCS for this chip
            extnum = imhdulist.index(imhdulist[ext])
            updatehdr.update_wcs(imhdulist, extnum, chip_wcs,
                                 wcsname=final_name, reusename=False,
                                 verbose=verbose)

        imhdulist.close()


def linearize(wcsim, wcsima, wcs_olddrz, wcs_newdrz, imcrpix, hx=1.0, hy=1.0):
    # linearization using 5-point formula for first order derivative
    x0 = imcrpix[0]
    y0 = imcrpix[1]
    p = np.asarray([[x0, y0],
                    [x0 - hx, y0],
                    [x0 - hx * 0.5, y0],
                    [x0 + hx * 0.5, y0],
                    [x0 + hx, y0],
                    [x0, y0 - hy],
                    [x0, y0 - hy * 0.5],
                    [x0, y0 + hy * 0.5],
                    [x0, y0 + hy]],
                   dtype=np.float64)
    # convert image coordinates to old drizzled image coordinates:
    p = wcs_olddrz.wcs_world2pix(wcsim.wcs_pix2world(p, 1), 1)
    # convert to sky coordinates using the new drizzled image's WCS:
    p = wcs_newdrz.wcs_pix2world(p, 1)
    # convert back to image coordinate system using partially (CRVAL only)
    # aligned image's WCS:
    p = wcsima.wcs_world2pix(p, 1).astype(ndfloat128)

    # derivative with regard to x:
    u1 = ((p[1] - p[4]) + 8 * (p[3] - p[2])) / (6*hx)
    # derivative with regard to y:
    u2 = ((p[5] - p[8]) + 8 * (p[7] - p[6])) / (6*hy)

    return (np.asarray([u1, u2]).T, p[0])


def update_chip_wcs(chip_wcs, drz_old_wcs, drz_new_wcs,
                    xrms=None, yrms=None):
    cd_eye = np.eye(chip_wcs.wcs.cd.shape[0])

    # estimate precision necessary for iterative processes:
    naxis1, naxis2 = chip_wcs.pixel_shape
    maxiter = 100
    crpix2corners = np.dstack([i.flatten() for i in np.meshgrid(
        [1, naxis1], [1, naxis2])])[0] - chip_wcs.wcs.crpix
    maxUerr = 1.0e-5 / np.amax(np.linalg.norm(crpix2corners, axis=1))

    # estimate step for numerical differentiation. We need a step
    # large enough to avoid rounding errors and small enough to get a
    # better precision for numerical differentiation.
    # TODO: The logic below should be revised at a later time so that it
    # better takes into account the two competing requirements.
    hx = max(1.0, min(20.0, (chip_wcs.wcs.crpix[0] - 1.0)/100.0,
                      (naxis1 - chip_wcs.wcs.crpix[0])/100.0))
    hy = max(1.0, min(20.0, (chip_wcs.wcs.crpix[1] - 1.0)/100.0,
                      (naxis2 - chip_wcs.wcs.crpix[1])/100.0))

    # compute new CRVAL for the image WCS:
    chip_wcs_orig = chip_wcs.deepcopy()
    crpix_in_old_drz = drz_old_wcs.wcs_world2pix([chip_wcs.wcs.crval], 1)
    chip_wcs.wcs.crval = drz_new_wcs.wcs_pix2world(crpix_in_old_drz, 1)[0]
    chip_wcs.wcs.set()

    # initial approximation for CD matrix of the image WCS:
    (U, u) = linearize(chip_wcs_orig, chip_wcs, drz_old_wcs, drz_new_wcs,
                       chip_wcs_orig.wcs.crpix, hx=hx, hy=hy)
    err0 = np.amax(np.abs(U-cd_eye)).astype(np.float64)
    chip_wcs.wcs.cd = np.dot(chip_wcs.wcs.cd.astype(ndfloat128), U).astype(np.float64)
    chip_wcs.wcs.set()

    # NOTE: initial solution is the exact mathematical solution (modulo numeric
    # differentiation). However, e.g., due to rounding errors, approximate
    # numerical differentiation, the solution may be improved by performing
    # several iterations. The next step will try to perform
    # fixed-point iterations to "improve" the solution
    # but this is not really required.

    # Perform fixed-point iterations to improve the approximation
    # for CD matrix of the image WCS (actually for the U matrix).
    for i in range(maxiter):
        (U, u) = linearize(chip_wcs_orig, chip_wcs, drz_old_wcs, drz_new_wcs,
                           chip_wcs_orig.wcs.crpix, hx=hx, hy=hy)
        err = np.amax(np.abs(U-cd_eye)).astype(np.float64)
        if err > err0:
            break
        chip_wcs.wcs.cd = np.dot(chip_wcs.wcs.cd, U).astype(np.float64)
        chip_wcs.wcs.set()
        if err < maxUerr:
            break
        err0 = err

    if xrms is not None:
        chip_wcs.wcs.crder = np.array([xrms,yrms])


#--------------------------------
# TEAL Interface functions
# (these functions are deprecated)
#---------------------------------
def run(configobj):
    # Interpret user-input from TEAL GUI and call function
    tweakback(configobj['drzfile'], newname = configobj['newname'],
            input=configobj['input'], origwcs = configobj['origwcs'],
            wcsname = configobj['wcsname'],
            extname=configobj['extname'],verbose=configobj['verbose'],
            force=configobj['force'])


#### Utility functions
#
def extract_input_filenames(drzfile):
    """
    Generate a list of filenames from a drizzled image's header
    """
    data_kws = fits.getval(drzfile, 'd*data', ext=0, memmap=False)
    if len(data_kws) == 0:
        return None
    fnames = []
    for kw in data_kws.cards:
        f = kw.value.split('[')[0]
        if f not in fnames:
            fnames.append(f)

    return fnames

def determine_extnum(drzfile, extname='SCI'):
    # Determine what kind of drizzled file input has been provided: MEF or single
    hdulist = fits.open(drzfile, memmap=False)
    numext = len(hdulist)
    sciext = 0
    for e,i in zip(hdulist,list(range(numext))):
        if 'extname' in e.header and e.header['extname'] == extname:
            sciext = i
            break
    hdulist.close()

    return sciext

def determine_orig_wcsname(header, wnames, wkeys):
    """
    Determine the name of the original, unmodified WCS solution
    """
    orig_wcsname = None
    orig_key = " "
    if orig_wcsname is None:
        for k,w in wnames.items():
            if w[:4] == 'IDC_':
                orig_wcsname = w
                orig_key = k
                break
    if orig_wcsname is None:
        # No IDC_ wcsname found... revert to second to last if available
        if len(wnames) > 1:
            orig_key = wkeys[-2]
            orig_wcsname = wnames[orig_key]
    return orig_wcsname,orig_key


util._def_help_functions(
    locals(), module_file=__file__, task_name=__taskname__, module_doc=__doc__
)
