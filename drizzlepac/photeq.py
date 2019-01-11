"""
A tool to adjust data values of images by equalizing each chip's PHOTFLAM value
to a single common value so that all chips can be treated equally
by ``AstroDrizzle``.

:Authors: Mihai Cara

:License: :doc:`LICENSE`

"""

__all__ = ['photeq']
__taskname__ = 'drizzlepac.photeq'
__version__ = '0.2'
__version_date__ = '06-Nov-2015'
__author__ = 'Mihai Cara'

# HISTORY:
# v0.1 - 14-Aug-2015 - Initial release
# v0.2 - 06-Nov-2015 - Added doctstrings, logging, TEAL interface


# STDLIB
import os
from datetime import datetime
import logging

# THIRD PARTY
import numpy as np
from astropy.io import fits

try:
    from stsci.tools import teal
except ImportError:
    teal = None

# LOCAL
from stsci.skypac import parseat, utils
from . import util

# create logger
logging.captureWarnings(capture=True)
_log = logging.getLogger(__name__)
_log.setLevel(logging.DEBUG)
_log_formatter = logging.Formatter('%(message)s')

# create console handler and set level to debug
_sh_log = logging.StreamHandler()
_sh_log.setLevel(logging.DEBUG)
_sh_log.setFormatter(_log_formatter)

# add _sh_log to logger
_log.addHandler(_sh_log)


def _mlinfo(msg, *args, **kwargs):
    lines = msg.splitlines()
    for line in lines:
        _log.info(line, *args, **kwargs)


def _mlwarn(msg, *args, **kwargs):
    lines = msg.splitlines()
    for line in lines:
        _log.warning(line, *args, **kwargs)


def photeq(files='*_flt.fits', sciext='SCI', errext='ERR',
           ref_phot=None, ref_phot_ext=None,
           phot_kwd='PHOTFLAM', aux_phot_kwd='PHOTFNU',
           search_primary=True,
           readonly=True, clobber=False, logfile='photeq.log'):
    """
    Adjust data values of images by equalizing each chip's PHOTFLAM value
    to a single common value so that all chips can be treated equally
    by ``AstroDrizzle``.

    Parameters
    ----------

    files : str (Default = ``'*_flt.fits'``)

        A string containing one of the following:

            * a comma-separated list of valid science image file names,
              e.g.: ``'j1234567q_flt.fits, j1234568q_flt.fits'``;

            * an @-file name, e.g., ``'@files_to_match.txt'``. See notes
              section for details on the format of the @-files.

        .. note::

            **Valid science image file names** are:

            * file names of existing FITS, GEIS, or WAIVER FITS files;

            * partial file names containing wildcard characters, e.g.,
              ``'*_flt.fits'``;

            * Association (ASN) tables (must have ``_asn``, or ``_asc``
              suffix), e.g., ``'j12345670_asn.fits'``.

    sciext : str (Default = 'SCI')
        Extension *name* of extensions whose data and/or headers should
        be corrected.

    errext : str (Default = 'ERR')
        Extension *name* of the extensions containing corresponding error
        arrays. Error arrays are corrected in the same way as science data.

    ref_phot : float, None (Default = None)
        A number indicating the new value of PHOTFLAM or PHOTFNU
        (set by 'phot_kwd') to which the data should be adjusted.

    ref_phot_ext : int, str, tuple, None (Default = None)
        Extension from which the `photeq` should get the reference photometric
        value specified by the `phot_kwd` parameter. This parameter is ignored
        if `ref_phot` **is not** `None`. When `ref_phot_ext` is `None`, then
        the reference inverse sensitivity value will be picked from the
        first `sciext` of the first input image containing `phot_kwd`.

    phot_kwd : str (Default = 'PHOTFLAM')
        Specifies the primary keyword which contains inverse sensitivity
        (e.g., PHOTFLAM). It is used to compute conversion factors by
        which data should be rescaled.

    aux_phot_kwd : str, None, list of str (Default = 'PHOTFNU')
        Same as `phot_kwd` but describes *other* photometric keyword(s)
        that should be corrected by inverse of the scale factor used to correct
        data. These keywords are *not* used to compute conversion factors.
        Multiple keywords can be specified as a Python list of strings:
        ``['PHOTFNU', 'PHOTOHMY']``.

        .. note::

            If specifying multiple secondary photometric keywords in the TEAL
            interface, use a comma-separated list of keywords.

    search_primary : bool (Default = True)
        Specifies whether to first search the primary header for the
        presence of `phot_kwd` keyword and compute conversion factor based on
        that value. This is (partially) ignored when `ref_phot` is not `None` in
        the sense that the value specified by `ref_phot` will be used as the
        reference *but* in all images primary will be searched for `phot_kwd`
        and `aux_phot_kwd` and those values will be corrected
        (if ``search_primary=True``).

    readonly : bool (Default = True)
        If `True`, `photeq` will not modify input files (nevertheless, it will
        convert input GEIS or WAVERED FITS files to MEF and could overwrite
        existing MEF files if `clobber` is set to `True`).
        The (console or log file) output however will be identical to the case
        when ``readonly=False`` and it can be examined before applying these
        changes to input files.

    clobber : bool (Default = False)
        Overwrite existing MEF files when converting input WAVERED FITS or GEIS
        to MEF.

    logfile : str, None (Default = 'photeq.log')
        File name of the log file.

    Notes
    -----

    By default, `photeq` will search for the first inverse sensitivity
    value (given by the header keyword specified by the `phot_kwd` parameter,
    e.g., PHOTFLAM or PHOTFNU) found in the input images and it will equalize
    all other images to this reference value.

    It is possible to tell `photeq` to look for the reference inverse
    sensitivity value only in a specific extension of input images, e.g.: 3,
    ('sci',3), etc. This can be done by setting `ref_phot_ext` to a specific
    extension. This may be useful, for example, for WFPC2 images: WF3 chip was
    one of the better calibrated chips, and so, if one prefers to have
    inverse sensitivities equalized to the inverse sensitivity of the WF3 chip,
    one can set ``ref_phot_ext=3``.

    Alternatively, one can provide their own reference inverse sensitivity
    value to which all other images should be "equalized" through the
    parameter `ref_phot`.

    .. note::

       Default parameter values (except for `files`, `readonly`, and `clobber`)
       should be acceptable for most HST images.

    .. warning::

       If images are intended to be used with ``AstroDrizzle``, it is
       recommended that sky background measurement be performed on "equalized"
       images as the `photeq` is not aware of sky user keyword in the image
       headers and thus it cannot correct sky values already recorded in the
       headers.

    Examples
    --------

    #. In most cases the default parameters should suffice:

           >>> from drizzlepac import photeq
           >>> photeq.photeq(files='*_flt.fits', readonly=False)

    #. If the re-calibration needs to be done on PHOTFNU rather than
       PHOTFLAM, then:

           >>> photeq.photeq(files='*_flt.fits', ref_phot='PHOTFNU',
           ... aux_phot_kwd='PHOTFLAM')

    #. If for WFPC2 data one desires that PHOTFLAM from WF3 be used as the
       reference in WFPC2 images, then:

           >>> photeq.photeq(files='*_flt.fits', ref_phot_ext=3) # or ('sci',3)

    """

    # Time it
    runtime_begin = datetime.now()

    # check that input file name is a string:
    if not isinstance(files, str):
        raise TypeError("Argument 'files' must be a comma-separated list of "
                        " file names")

    # Set-up log files:
    if isinstance(logfile, str):
        # first, in case there are any "leftover" file handlers,
        # close and remove them:
        for h in _log.handlers:
            if h is not _sh_log and isinstance(h, logging.FileHandler):
                h.close()
                _log.removeHandler(h)
        # create file handler:
        log_formatter = logging.Formatter('[%(levelname)s:] %(message)s')
        log_file_handler = logging.FileHandler(logfile)
        log_file_handler.setFormatter(log_formatter)
        # add log_file_handler to logger
        _log.addHandler(log_file_handler)

    elif logfile is not None:
        raise TypeError("Unsupported 'logfile' type")

    #  BEGIN:
    _mlinfo("***** {0} started on {1}".format(__taskname__, runtime_begin))
    _mlinfo("      Version {0} ({1})".format(__version__, __version_date__))

    # check that extension names are strings (or None for error ext):
    if sciext is None:
        sci_ext4parse = '*'
        ext2get = None
    else:
        if not isinstance(sciext, str):
            raise TypeError("Argument 'sciext' must be a string or None")
        sciext = sciext.strip()
        if sciext.upper() == 'PRIMARY':
            sciext = sciext.upper()
            ext2get = (sciext, 1)
        else:
            ext2get = (sciext, '*')

        sci_ext4parse = ext2get

    if errext is not None and not isinstance(errext, str):
        raise TypeError("Argument 'errext' must be a string or None")

    # check that phot_kwd is supported:
    if not isinstance(phot_kwd, str):
        raise TypeError("Argument 'phot_kwd' must be a string")
    phot_kwd = phot_kwd.strip().upper()

    # check that ref_phot_ext has correct type:
    if ref_phot_ext is not None and not \
       (isinstance(ref_phot_ext, int) or isinstance(ref_phot_ext, str) \
        or (isinstance(ref_phot_ext, tuple) and len(ref_phot_ext) == 2 \
            and isinstance(ref_phot_ext[0], str) and \
            isinstance(ref_phot_ext[1], int))):
        raise TypeError("Unsupported 'ref_phot_ext' type")
    if isinstance(ref_phot_ext, str):
        ref_phot_ext = (ref_phot_ext, 1)

    if aux_phot_kwd is None:
        aux_phot_kwd = []

    elif isinstance(aux_phot_kwd, str):
        aux_phot_kwd = [aux_phot_kwd.strip().upper()]
        if phot_kwd == aux_phot_kwd:
            raise ValueError("Auxiliary photometric keyword must be different "
                             "from the main photometric keyword 'phot_kwd'.")

    elif hasattr(aux_phot_kwd, '__iter__'):
        if not all([isinstance(phot, str) for phot in aux_phot_kwd]):
            raise TypeError("Argument 'aux_phot_kwd' must be a string, list of "
                        "strings, or None")
        aux_phot_kwd = [phot.strip().upper() for phot in aux_phot_kwd]
        if ref_phot in aux_phot_kwd:
            raise ValueError("Auxiliary photometric keyword(s) must be "
                             "different from the main photometric keyword "
                             "'phot_kwd'.")

    else:
        raise TypeError("Argument 'aux_phot_kwd' must be a string, list of "
                        "strings, or None")

    # read input file list:
    fl = parseat.parse_cs_line(csline=files, default_ext=sci_ext4parse,
                               im_fmode='readonly' if readonly else 'update',
                               clobber=clobber, fnamesOnly=True,
                               doNotOpenDQ=True)

    # check if user supplied file extensions, set them to the sciext,
    # and warn that they will be ignored:
    for f in fl:
        if f.count > 1 or f.fext[0] != sci_ext4parse:
            _mlwarn("WARNING: Extension specifications for file {:s} "
                    "will be ignored. Using all {:s} extensions instead."
                    .format(f.image,  'image-like' if sciext is None else \
                            "{:s}".format(utils.ext2str(sciext,
                                                        default_extver=None))))

    # find the reference PHOTFLAM/PHOTNU:
    flc = fl[:]
    ref_hdu = None
    ref_ext = None
    ref_user = True

    if ref_phot is None:
        ref_user = False
        for f in flc:
            f.convert2ImageRef()

            # get primary hdu:
            pri_hdu = f.image.hdu[0]

            # find all valid extensions:
            if ref_phot_ext is None:
                if sciext == 'PRIMARY':
                    extnum = [0]
                else:
                    extnum = utils.get_ext_list(f.image, sciext)

                is_pri_hdu = [f.image.hdu[ext] is pri_hdu for ext in extnum]

                # if necessary, add primary header to the hdu list:
                if search_primary:
                    try:
                        pri_index = is_pri_hdu.index(True)
                        extnum.insert(0, extnum.pop(pri_index))
                    except ValueError:
                        extnum.insert(0, 0)

            else:
                extnum = [ref_phot_ext]

            for ext in extnum:
                hdu = f.image.hdu[ext]
                if phot_kwd in hdu.header:
                    ref_phot = hdu.header[phot_kwd]
                    ref_ext = ext
                    ref_hdu = hdu
                    break

            if ref_phot is None:
                _mlwarn("WARNING: Could not find specified inverse "
                        "         sensitivity keyword '{:s}'\n"
                        "         in any of the {} extensions of file '{}'.\n"
                        "         This input file will be ignored."
                        .format(phot_kwd, 'image-like' if sciext is None else \
                                "{:s}".format(utils.ext2str(sciext,
                                                            default_extver=None)),
                                os.path.basename(f.image.original_fname)))
                f.release_all_images()
                fl.remove(f)

            else:
                break

    if ref_phot is None:
        raise RuntimeError("Could not find the inverse sensitivity keyword "
                           "'{:s}' in the specified headers of "
                           "the input image(s).\nCannot continue."
                           .format(phot_kwd))

    aux_phot_kwd_list = ','.join(aux_phot_kwd)

    _mlinfo("\nPRIMARY PHOTOMETRIC KEYWORD: {:s}".format(phot_kwd))
    _mlinfo("SECONDARY PHOTOMETRIC KEYWORD(S): {:s}"
              .format(aux_phot_kwd_list if aux_phot_kwd_list else 'None'))
    if ref_user:
        _mlinfo("REFERENCE VALUE PROVIDED BY USER: '{:s}'={}\n"
                .format(phot_kwd, ref_phot))
    else:
        _mlinfo("REFERENCE VALUE FROM FILE: '{:s}[{:s}]'\n"
                .format(os.path.basename(f.image.original_fname),
                          utils.ext2str(ref_ext)))
        _mlinfo("REFERENCE '{:s}' VALUE IS: {}".format(phot_kwd, ref_phot))

    # equalize PHOTFLAM/PHOTNU
    for f in fl:
        # open the file if necessary:
        if f.fnamesOnly:
            _mlinfo("\nProcessing file '{:s}'".format(f.image))
            f.convert2ImageRef()
        else:
            _mlinfo("\nProcessing file '{:s}'".format(f.image.original_fname))

        # first, see if photflam is in the primary header and save this value:
        pri_conv = None
        if search_primary:
            whdu = f.image.hdu[0]
            if phot_kwd in whdu.header:
                _mlinfo("   * Primary header:")
                if whdu is ref_hdu:
                    pri_conv = 1.0
                    _mlinfo("     - '{}' = {} found in the primary header."
                            .format(phot_kwd, whdu.header[phot_kwd]))
                    _mlinfo("     - Data conversion factor based on primary "
                              "header: {}".format(pri_conv))
                else:
                    _mlinfo("     - '{}' found in the primary header."
                            .format(phot_kwd))
                    pri_conv = whdu.header[phot_kwd] / ref_phot
                    _mlinfo("     - Setting {:s} in the primary header to {} "
                              "(old value was {})"
                            .format(phot_kwd, ref_phot, whdu.header[phot_kwd]))
                    _mlinfo("     - Data conversion factor based on primary "
                            "header: {}".format(pri_conv))
                    whdu.header[phot_kwd] = ref_phot

            # correct the "other" photometric keyword, if present:
            if pri_conv is not None and whdu is not ref_hdu:
                for aux_kwd in aux_phot_kwd:
                    if aux_kwd in whdu.header:
                        old_aux_phot = whdu.header[aux_kwd]
                        new_aux_phot = old_aux_phot / pri_conv
                        whdu.header[aux_kwd] = new_aux_phot
                        _mlinfo("     - Setting {:s} in the primary header "
                                "to {} (old value was {})"
                                .format(aux_kwd, new_aux_phot, old_aux_phot))

            # process data and error arrays when 'sciext' was specifically set to
            # 'PRIMARY':
            if sciext == 'PRIMARY' and pri_conv is not None:
                has_data = (hasattr(whdu, 'data') and
                            whdu.data is not None)

                # correct data:
                if has_data:
                    if np.issubdtype(whdu.data.dtype, np.floating):
                        whdu.data *= pri_conv
                        _mlinfo("     - Data have been multiplied by {}"
                                .format(pri_conv))
                    else:
                        _mlwarn("WARNING: Data not converted because it is of "
                                "non-floating point type.")

                # correct error array:
                if errext is not None:
                    eext = (errext, 1)
                    try:
                        whdu = f.image.hdu[eext]
                    except KeyError:
                        _mlwarn("     - WARNING: Error extension {:s} not found."
                                .format(utils.ext2str(eext)))

                        f.release_all_images()
                        continue

                    if hasattr(whdu, 'data') and whdu.data is not None:
                        if np.issubdtype(whdu.data.dtype, np.floating):
                            whdu.data *= pri_conv
                            _mlinfo("     - Error array (ext={}) has been "
                                    "multiplied by {}".format(eext, pri_conv))
                        else:
                            _mlinfo("     - Error array in extension {:s} "
                                    "contains non-floating point data.\n"
                                    "       Skipping this extension"
                                    .format(utils.ext2str(ext)))

                f.release_all_images()
                continue

        # find all valid extensions:
        extnum = utils.get_ext_list(f.image, sciext)

        for ext in extnum:
            whdu = f.image.hdu[ext]
            conv = None

            if whdu is ref_hdu:
                _mlinfo("   * EXT: {} - This is the \"reference\" extension.\n"
                        "          Nothing to do. Skipping this extension..."
                        .format(ext))
                continue

            has_data = (hasattr(whdu, 'data') and
                        whdu.data is not None)

            if has_data and not np.issubdtype(whdu.data.dtype, np.floating):
                _mlinfo("   * EXT: {} contains non-floating point data. "
                        "Skipping this extension".format(ext))

            # find all auxiliary photometric keywords present in the header:
            paux = [aux_kwd for aux_kwd in aux_phot_kwd if aux_kwd \
                    in whdu.header]

            if phot_kwd in whdu.header:
                _mlinfo("   * EXT: {}".format(ext))
                old_phot = whdu.header[phot_kwd]
                conv = old_phot / ref_phot
                _mlinfo("     - Setting {:s} to {} (old value was {})"
                        .format(phot_kwd, ref_phot, old_phot))
                whdu.header[phot_kwd] = ref_phot
                _mlinfo("     - Computed conversion factor for data: {}"
                        .format(conv))

            elif pri_conv is None:
                _mlinfo("   * EXT: {}".format(ext))
                _mlinfo("     - '{:s} not found. Skipping this extension..."
                        .format(phot_kwd))
                continue

            else:
                _mlinfo("   * EXT: {}".format(ext))

                # if paux:
                    # print("ERROR: Primary photometric keyword ('{:s}') is "
                          # "missing but\n       the secondary keywords ('{:s}') "
                          # "are present. This extension cannot be processed."
                          # .format(phot_kwd, ','.join(paux)))
                    # continue

                _mlinfo("     - '{:s} not found. Using conversion factor "
                        "based\n       on the primary header: {}"
                        .format(phot_kwd, pri_conv))
                conv = pri_conv

            # correct the "other" photometric keyword, if present:
            if conv is not None:
                for aux_kwd in paux:
                    old_aux_phot = whdu.header[aux_kwd]
                    new_aux_phot = old_aux_phot / conv
                    whdu.header[aux_kwd] = new_aux_phot
                    _mlinfo("     - Setting {:s} to {} (old value was {})"
                            .format(aux_kwd, new_aux_phot, old_aux_phot))

            # correct data:
            if has_data:
                if conv is None:
                    _mlinfo("   * EXT: {}".format(ext))

                if np.issubdtype(whdu.data.dtype, np.floating):
                    whdu.data *= conv
                    _mlinfo("     - Data have been multiplied by {}"
                            .format(conv))
                else:
                    _mlinfo("WARNING: Non-floating point data. Data cannot "
                            "be re-scaled.")

            # correct error array:
            if errext is not None and isinstance(ext, tuple) and len(ext) == 2:
                eext = (errext, ext[1])
                try:
                    whdu = f.image.hdu[eext]
                except KeyError:
                    continue

                if hasattr(whdu, 'data') and whdu.data is not None:
                    if np.issubdtype(whdu.data.dtype, np.floating):
                        whdu.data *= conv
                        _mlinfo("     - Error array (ext={}) has been "
                                "multiplied by {}".format(eext, conv))
                    else:
                        _mlinfo("     - Error array in extension {:s} "
                                "contains non-floating point data.\n"
                                "       Skipping this extension"
                                .format(utils.ext2str(ext)))

        f.release_all_images()

    _mlinfo("\nDone.")

    if readonly:
        _mlinfo("\nNOTE: '{:s}' was run in READONLY mode\n"
                "       and input image(s)' content WAS NOT MODIFIED."
                .format(__taskname__))

    # close all log file handlers:
    for h in _log.handlers:
        if h is not _sh_log and isinstance(h, logging.FileHandler):
            h.close()
            _log.removeHandler(h)


#--------------------------
# TEAL Interface functions
#--------------------------
def run(configObj):
    logfile = configObj['logfile'] if len(configObj['logfile'].strip()) > 0 \
        else None
    photeq(files=configObj['files'],
           sciext=configObj['sciext'],
           errext=configObj['errext'],
           ref_phot=configObj['ref_phot'],
           ref_phot_ext=util.check_blank(configObj['ref_phot_ext']),
           phot_kwd=configObj['phot_kwd'],
           aux_phot_kwd=_split_kwd_list(configObj['aux_phot_kwd']),
           search_primary=configObj['search_primary'],
           readonly=configObj['readonly'],
           clobber=configObj['clobber'],
           logfile=logfile)


def _split_kwd_list(kwd_list):
    kwdl = kwd_list.split(',')
    newl = []
    for kwd in kwdl:
        k = kwd.strip()
        if k:
            newl.append(k)
    nkwd = len(newl)
    if nkwd == 0:
        return None
    elif nkwd == 1:
        return newl[0]
    else:
        return newl


def getHelpAsString(docstring = False, show_ver = True):
    """
    return useful help from a file in the script directory called
    __taskname__.help

    """
    install_dir = os.path.dirname(__file__)
    taskname = util.base_taskname(__taskname__, __package__)
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
                helpString += __doc__ + os.linesep + photeq.__doc__
    else:
        helpString = 'file://' + htmlfile

    return helpString


def help(file=None):
    """
    Print out syntax help for running skymatch

    Parameters
    ----------
    file : str (Default = None)
        If given, write out help to the filename specified by this parameter
        Any previously existing file with this name will be deleted before
        writing out the help.
    """
    helpstr = getHelpAsString(docstring=True)
    if file is None:
        print(helpstr)
    else:
        if os.path.exists(file): os.remove(file)
        f = open(file,mode='w')
        f.write(helpstr)
        f.close()
