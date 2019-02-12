"""
This module provides functions for mapping DS9 region files given in sky
coordinates to DS9 region files specified in image coordinates
of multiple images using the WCS information from the images.

:Authors: Mihai Cara

:License: :doc:`LICENSE`

"""
from astropy.io import fits
import stregion as pyregion
import stwcs
import os
from stsci.tools.fileutil import findExtname
from stsci.tools import teal
from . import util
from .regfilter import fast_filter_outer_regions


# This is specifically NOT intended to match the package-wide version information.
__version__ = '0.1'
__version_date__ = '11-Nov-2013'
__taskname__ = 'mapreg'
__author__ = 'Mihai Cara'


class _AuxSTWCS:
    def __init__(self, fobj=None, ext=None, minerr=0.0, wcskey=" "):
        if isinstance(fobj, stwcs.wcsutil.HSTWCS):
            self._stwcs = fobj
        else:
            self._stwcs = stwcs.wcsutil.HSTWCS(fobj=fobj, ext=ext,
                                               minerr=minerr, wcskey=wcskey)

    @property
    def wcs(self):
        return self._stwcs.wcs

    def sub(self, axes):
        return self._stwcs.sub(axes=axes)

    def wcs_sky2pix(self, *args, **kwargs):
        ar = list(args[:])
        if 'origin' in kwargs:
            ar.append(kwargs['origin'])
        return self._stwcs.all_world2pix( *tuple(ar) )

    def wcs_pix2sky(self, *args, **kwargs):
        ar = list(args[:])
        if 'origin' in kwargs:
            ar.append(kwargs['origin'])
        return self._stwcs.all_pix2world( *tuple(ar) )


def MapReg(input_reg, images, img_wcs_ext='sci', refimg='', ref_wcs_ext='sci',
           chip_reg='', outpath='./regions', filter='', catfname='', iteractive=False,
           append=False, verbose=True):
    from .util import check_blank
    from .tweakutils import parse_input

    input_reg       = _simple_parse_teal_fname( check_blank(input_reg),
                                                parse_at = True )
    images_par, cat = parse_input( check_blank(images) )
    img_wcs_ext_par = _simple_parse_teal_extn( check_blank(img_wcs_ext) )
    refimg_par      = check_blank(refimg) #TODO: may need something similar to
                                          # _simple_parse_teal_fname when we will
                                          # support it
    ref_wcs_ext_par = _simple_parse_teal_extn( check_blank(ref_wcs_ext) )
    chip_reg_par    = _simple_parse_teal_fname( check_blank(chip_reg),
                                                parse_at = False )
    outpath_par     = check_blank(outpath)
    filter_par      = check_blank(filter)

    catfname_par    =  check_blank(catfname)

    map_region_files( input_reg, images=images_par, img_wcs_ext=img_wcs_ext_par,
                      refimg=refimg_par, ref_wcs_ext=ref_wcs_ext_par,
                      chip_reg=chip_reg_par,
                      outpath=outpath_par, filter=filter_par,
                      catfname = catfname_par, iteractive=iteractive,
                      append=append, verbose=verbose )


def _simple_parse_teal_fname(fnamestr, parse_at=False):
    if not fnamestr: return None
    fnames = [fname if fname else None \
              for fname in map(str.strip,fnamestr.split(','))]
    if not parse_at:
        if ',' in fnamestr:
            return fnames
        else:
            return fnames[0]
    fnamelst = []
    for fname in fnames:
        if not fname:
            continue
        if fname[0] == '@':
            fh = open(fname[1:],'r')
            lines = fh.readlines()
            fh.close()
            for line in lines:
                f = _simple_parse_teal_fname(
                            map(str.strip,line.replace(',',' ').split())[0],
                            parse_at = parse_at )
                if isinstance(f, str):
                    fnamelst.append(f)
                else:
                    fnamelst += f
        else:
            fnamelst.append(fname)
    if ',' in fnamestr or len(fnamelst) > 1:
        return fnamelst
    else:
        return fnamelst[0]


def _simple_parse_teal_extn(extnstr):
    if extnstr and isinstance(extnstr, str) and extnstr.strip():
        extlst = list(map(str.strip,extnstr.split(';')))
        if len(extlst) == 1:
            return extlst[0]
        else:
            return extlst
    else:
        return None


def map_region_files(input_reg, images, img_wcs_ext='sci',
                     refimg=None, ref_wcs_ext='sci', chip_reg=None,
                     outpath='./regions', filter=None, catfname=None,
                     iteractive=False, append=False, verbose=True):
    # Check that output directory exists:
    if outpath in [None, ""]:
        outpath = os.path.curdir + os.path.sep
    elif not os.path.isdir(outpath):
        raise IOError("The output directory \'%s\' does not exist." % outpath)

    if filter is not None and filter.lower() not in ["fast","precise"]:
        raise TypeError("The 'filter' argument can be None, 'fast', " \
                        "or 'precise'.")

    #TODO: Implement the "precise" checking of region intersection
    # with the image's bounding box.
    if filter is not None and filter.lower() == "precise":
        _print_warning("\"precise\" filter option is not yet supported. "\
                       "Reverting to filter = \"fast\" instead.")
        filter = "fast"

    # create a 2D list of region - header (with WCS) pairs:
    # [[reg1,ref_wcs_header1],[reg2,ref_wcs_header2],...]
    regheaders = build_reg_refwcs_header_list(input_reg, refimg, ref_wcs_ext,
                                              verbose)

    # create a list of input files *without* extensions (if present) and
    # a 2D list of "chip" region - image extension pairs:
    # [[reg1, ext1], [reg2,ext2],...]
    imgfnames, cregext = build_img_ext_reg_list(images, chip_reg, img_wcs_ext,
                                                verbose)

    # build a "master" list of regions in sky coordinates:
    all_sky_regions = []
    for p in regheaders:
        if p[1] is None: # Already in sky WCS
            all_sky_regions += p[0]
        else:
            all_sky_regions += shapes_img2sky(p[0], p[1]) #TODO: implement shapes_img2sky

    all_sky_regions = pyregion.ShapeList(list(all_sky_regions))

    cattb = []
    # create a region file (with regions in image coordinates) for each
    # image extension from the input_reg and chip_reg regions
    for fname in imgfnames:
        imghdu = None
        imghdu = fits.open(fname, memmap=False)
        catreg = []
        try:
            for extp in cregext:
                ext = extp[1]

                # check that the extension is of supported type:
                if not _is_supported_hdu(imghdu[ext].header):
                    raise RuntimeError("Extension {} is of unsupported " \
                                       "type.".format(ext))

                # Remove references to non-SIP distortions from the header
                #TODO: remove this line of code as well as the remove_header_tdd
                # function once publicly available release of pyregion uses
                # all_world2pix and all_pix2world functions and does this header
                # cleanup itself.
#                remove_header_tdd(imghdu[ext].header)

                ####### added to pass hstwcs instead of header to pyregion
#                if isinstance(ext, str):
#                    ext = findExtname(imghdu, extname = ext, extver = 1)
#                elif isinstance(ext, tuple):
#                    ext = findExtname(imghdu, extname = ext[0], extver = ext[1])

#                wcs = stwcs.wcsutil.HSTWCS(imghdu, ext=ext)
                wcs    = _AuxSTWCS(imghdu, ext=ext)
                extreg = all_sky_regions.as_imagecoord(wcs, rot_wrt_axis=2)
                #######

                #extreg = all_sky_regions.as_imagecoord(imghdu[ext].header, rot_wrt_axis=2)

                if extp[0]: # add "fixed" regions if any
                    extreg = pyregion.ShapeList(extreg+extp[0])

                # filter out regions outside the image:
                if filter and filter.lower() == 'fast':
                    fast_filter_outer_regions( extreg,
                                               imghdu[ext].header['NAXIS1'],
                                               imghdu[ext].header['NAXIS2'],
                                               origin = 1 )

                if len(extreg) < 1: # do not create empty region files
                    catreg.append('None')
                    continue

                # generate output region file name:
                extsuffix = _ext2str_suffix(ext)
                basefname, fext = os.path.splitext(os.path.basename(fname))
                regfname = basefname + extsuffix + os.extsep + "reg"
                fullregfname = os.path.join(outpath, regfname)

                catreg.append(regfname)

                # save regions to a file:
                if append and os.path.isfile(fullregfname):
                    old_extreg = pyregion.open(fullregfname)
                    extreg = pyregion.ShapeList(old_extreg + extreg)

                #TODO: pyregion.write does not work well. For now use _regwite
                # (until the bug fixes get to the pyregion project).
                #TODO: we replaced pyregion.write with our implementation
                # of the write (code from _regwrite) in the locally maintained
                # pyregion until these changes get to be implemented in the
                # publicly available release of pyregion.
                #
                _regwrite(extreg, fullregfname)
                #extreg.write(fullregfname) # <- use this instead of _regwrite
                                            # once the pyregion bugs are fixed.
            cattb.append([fname, catreg])
        except:
            if imghdu:
                imghdu.close()
            raise

    # create exclusions catalog file:
    if catfname:
        catfh = open(os.path.join(outpath, catfname), 'w')
        for catentry in cattb:
            catfh.write(catentry[0]) # image file name
            for reg in catentry[1]:
                catfh.write(' ' + reg) # region file name
            catfh.write('\n')
        catfh.close()


def remove_header_tdd(hdr):
    # a workaround to pyregion using FITS 'header' (instead of 'WCS')...
    # For some images header alone is not enough...
    # Removes offending corrections from the header...
    #
    # Code below is taken on 'remove_distortion_keywords' from fitsblender/blendheaders.py
    #
    distortion_kws = ['TDDALPHA','TDDBETA','D2IMEXT','D2IMERR',
                      'DGEOEXT','NPOLEXT']

    # Remove any reference to TDD correction from
    #    distortion-corrected products
    # We also need to remove the D2IM* keywords so that HSTWCS/PyWCS
    # does not try to look for non-existent extensions
    for kw in distortion_kws:
        if kw in hdr:
            try:
                del hdr[kw]
            except KeyError:
                pass

    # Remove paper IV related keywords related to the
    #   DGEO correction here
    for k in list(hdr.items()):
        try:
            if (k[0][:2] == 'DP'):
                del hdr[k[0]+'*']
                del hdr[k[0]+'.*']
                del hdr[k[0]+'.*.*']
            if (k[0][:2] == 'CP'):
                del hdr[k[0]]
        except KeyError:
            pass


def _regwrite(shapelist,outfile):
    """ Writes the current shape list out as a region file """
    # This function corrects bugs and provides improvements over the pyregion's
    # ShapeList.write method in the following:
    #
    # 1. ShapeList.write crashes if regions have no comments;
    # 2. ShapeList.write converts 'exclude' ("-") regions to normal regions ("+");
    # 3. ShapeList.write does not support mixed coordinate systems in a
    #    region list.
    #
    # NOTE: This function is provided as a temoprary workaround for the above
    #    listed problems of the ShapeList.write. We hope that a future version
    #    of pyregion will address all these issues.
    #
    #TODO: Push these changes to pyregion.

    if len(shapelist) < 1:
        _print_warning("The region list is empty. The region file \"%s\" "\
                       "will be empty." % outfile)
        try:
            outf = open(outfile,'w')
            outf.close()
            return
        except IOError as e:
            cmsg = "Unable to create region file \'%s\'." % outfile
            if e.args:
                e.args = (e.args[0] + "\n" + cmsg,) + e.args[1:]
            else:
                e.args=(cmsg,)
            raise e
        except:
            raise

    prev_cs = shapelist[0].coord_format

    outf = None
    try:
        outf = open(outfile,'w')

        attr0 = shapelist[0].attr[1]
        defaultline = " ".join(["%s=%s" % (a,attr0[a]) for a in attr0 \
                                if a!='text'])

        # first line is globals
        print("global", defaultline, file=outf)
        # second line must be a coordinate format
        print(prev_cs, file=outf)

        for shape in shapelist:
            shape_attr = '' if prev_cs == shape.coord_format \
                            else shape.coord_format+"; "
            shape_excl = '-' if shape.exclude else ''
            text_coordlist = ["%f" % f for f in shape.coord_list]
            shape_coords = "(" + ",".join(text_coordlist) + ")"
            shape_comment = " # " + shape.comment if shape.comment else ''

            shape_str = shape_attr + shape_excl + shape.name + shape_coords + \
                        shape_comment

            print(shape_str, file=outf)

    except IOError as e:
        cmsg = "Unable to create region file \'%s\'." % outfile
        if e.args:
            e.args = (e.args[0] + "\n" + cmsg,) + e.args[1:]
        else:
            e.args=(cmsg,)
        if outf: outf.close()
        raise e
    except:
        if outf: outf.close()
        raise

    outf.close()


def _is_supported_hdu(header):
    # check that NAXIS == 2:
    if not 'NAXIS' in header or header['NAXIS'] != 2:
        return False
    if header['NAXIS1'] < 1 or header['NAXIS2'] < 1:
        return False
    # check that this is either a primary HDU or an 'IMAGE' extension
    return ('SIMPLE' in header) or \
           ('XTENSION' in header and header['XTENSION'] == 'IMAGE')


def _ext2str_suffix(ext):
    if isinstance(ext, tuple):
        return "_{}{}_twreg".format(ext[0],ext[1])
    elif isinstance(ext, int):
        return "_extn{}_twreg".format(ext)
    else:
        return "_{}_twreg".format(ext) # <- we should not get here...


def shapes_img2sky(reglist, wcs):
    #TODO: implement shapes_img2sky to do something useful
    # (like real CS transformation)
    from pyregion.wcs_helper import image_like_coordformats
    return [r for r in reglist
            if r.coord_format not in image_like_coordformats]


def build_reg_refwcs_header_list(input_reg, refimg, ref_wcs_ext, verbose):
    # check input region file names and initialize region lists:
    region_lists = []
    single_inp_reg = False

    if isinstance(input_reg, str) and len(input_reg) > 0: # a single region file
        single_inp_reg  = True
        input_regfnames = [input_reg]
        ref_wcs_exts    = [0] # for a single input region we assume that the
                              # default location of the WCS is in the primary
                              # header (if nothing else specified)

    elif isinstance(input_reg, list) and input_reg: # a non-empty list of files
        # select non-zero length file names;replace all other elements with None
        input_regfnames = [fname if isinstance(fname,str) and \
                           len(fname)>0 else None for fname in input_reg]

        # Populate ref_wcs_exts with "default" FITS extension numbers according
        # to the input region position in the input_regfnames list
        # (starting with 1). That is, if multiple regions are given, then we
        # assume that (if no extension name is provided) the WCS is located in
        # "extensions" of the input reference FITS file (as opposite to the
        # primary header).
        ref_wcs_exts = list(range(1, len(input_regfnames)+1))

        # Filter out elements of ref_wcs_exts and input_regfnames that
        # correspond to None elements in input_regfnames:
        ref_wcs_exts = [ext for reg,ext in zip(input_regfnames,ref_wcs_exts) \
                        if reg is not None]
        input_regfnames = [reg for reg in input_regfnames if reg is not None]

    else:
        input_regfnames = []
        ref_wcs_exts    = [] # unnecessary, really, because of the next "if"

    # Check that we have at least one "valid" input region file name
    if len(input_regfnames) < 1:
        raise RuntimeError("Parameter 'input_reg' must be either a non-zero " \
            "length string or a non-empty list with at least "               \
            "one non-zero length string file name.")

    # Check that the region files exist and try to open/read them:
    ireg = 0
    for fname in input_regfnames:
        # check region file existence:
        if not os.path.isfile(fname):
            raise IOError("The input region file \'%s\' does not exist." % \
                fname)
        # try to read regions:
        try:
            reglist = pyregion.open(fname)
            region_lists.append(reglist)

        except IOError:
            raise IOError("Unable to open the input region file \'%s\'." % \
                fname)

        except:
            raise RuntimeError("Unexpected error parsing input region " \
                "file \'%s\'. Suspecting either a corrupt file or "     \
                "an unrecognizable region file format." % fname)

        # check if WCS is needed to convert from image to sky CS
        if not _needs_ref_WCS(reglist):
            ref_wcs_exts[ireg] = None

        ireg += 1

    # If WCS is needed for conversion to sky coordinates, find the correct
    # FITS extension based either on extension specified in the input reference
    # image 'refimg', optional input argument 'ref_wcs_ext', and/or position of
    # the region file in the input region list 'input_reg'
    if [True for ireg in ref_wcs_exts if ireg is not None]:
        # A reference WCS is required. Check that 'refimg' is an existing file:
        if refimg is None:
            raise ValueError("Argument 'refimg' cannot be None when some " \
                             "input regions are given in logical coordinates.")

        if not isinstance(refimg, str):
            raise TypeError("Argument 'refimg' must be a string with the " \
                            "name of an existing FITS file and, optionally, " \
                            "followed by an extension specification.")

        (refimg_fname, frefext) = extension_from_filename(refimg)

        if not os.path.isfile(refimg_fname):
            raise IOError("The reference FITS file \'%s\' does not exist." % \
                          refimg_fname)

        if frefext:
            # Try to determine the extension name from the reference file name
            refext = parse_ext(frefext)

            if single_inp_reg:

                if isinstance(refext, tuple) or isinstance(refext, int):
                    ref_wcs_exts = [refext]
                elif isinstance(refext, str): # it is a string:
                    ref_wcs_exts = [(refext, 1)]
                else:
                    raise RuntimeError("Logical error in the code.")

            else:
                # check that FITS extension specified in the reference image file
                # name is not too restrictive (for multiple regions it cannot
                # contain a specific extver):
                if isinstance(refext, tuple):
                    raise RuntimeError("Extension version ('EXTVER') in the "  \
                    "reference file name should not be present when multiple " \
                    "region files (in logical CS) are provided as input. "     \
                    "Only extension name ('EXTNAME') is allowed (e.g., "       \
                    "[SCI], [DQ], etc.)")

                if isinstance(refext, int):
                    raise RuntimeError("Extension number (e.g., [0], [2], "    \
                    "etc.) in the reference file name should not be present "  \
                    "when multiple region files (in logical CS) are "          \
                    "provided as input. Only extension name ('EXTNAME') is "   \
                    "allowed (e.g., [SCI], [DQ], etc.)")

                ref_wcs_exts = [None if extn is None else (refext, extn) \
                                for extn in ref_wcs_exts]

        elif ref_wcs_ext:
            # Try to determine the extension name from the 'ref_wcs_ext'
            # argument:
            refext = parse_ext(ref_wcs_ext)

            if single_inp_reg:
                if isinstance(refext, tuple) or isinstance(refext, int):
                    ref_wcs_exts = [refext]
                elif isinstance(refext, str): # it is a string:
                    ref_wcs_exts = [(refext, 1)]
                else:
                    raise RuntimeError("Logical error in the code.")

            else:
                # check that FITS extension specified in the reference image
                # file name is not too restrictive (for multiple regions it
                # cannot contain a specific extver):
                if isinstance(refext, tuple):
                    raise RuntimeError("Extension version ('EXTVER') in the "  \
                    "'ref_wcs_ext' argument should not be present when "       \
                    "multiple region files (in logical CS) are provided "      \
                    "as input. Only extension name ('EXTNAME') is allowed "    \
                    "(e.g., [SCI], [DQ], etc.)")

                if isinstance(refext, int):
                    raise RuntimeError("Extension number (e.g., [0], [2], "    \
                    "etc.) in the 'ref_wcs_ext' argument should not be "       \
                    "present when multiple region files (in logical CS) are "  \
                    "provided as input. Only extension name ('EXTNAME') is "   \
                    "allowed (e.g., [SCI], [DQ], etc.)")

                ref_wcs_exts = [None if extn is None else (refext, extn) \
                                for extn in ref_wcs_exts]
        # check if the extensions found above are present
        # in the reference WCS FITS file:
        if (refimg is not None and
            not _check_FITS_extensions(refimg, ref_wcs_exts)):
            raise RuntimeError("Not all FITS extensions derived based on the " \
                "input region(s) and extension names have been found in the "  \
                "input reference WCS file. Unable to proceed.")

        # load headers containing WCS from the reference input FITS file:
        ref_wcs_headers = [None if extn is None \
                           else fits.getheader(refimg_fname, ext=extn, memmap=False) \
                           for extn in ref_wcs_exts ] #TODO: return WCS instead of header

    else:
        # no WCS is needed (regions are already in sky coordinates):
        ref_wcs_headers = [None for reg in region_lists]

    # Return a list of pairs of the form [region, header] with missing regions
    # (i.e., region = None) filtered out. (We do not need to keep order anymore)
    return [p for p in map(list, list(zip(region_lists, ref_wcs_headers))) \
            if p[0] is not None]


def build_img_ext_reg_list(images, chip_reg=None, img_wcs_ext='sci',
                           verbose=True):
    # Check that the 'chip_reg' argument is of correct type/form:
    if chip_reg is None:
        nreg     = 0
        multireg = False
        chip_reg = []
    elif isinstance(chip_reg, str):
        nreg = 1
        chip_reg = [chip_reg]
        multireg = False
    elif isinstance(chip_reg, list):
        chip_reg = chip_reg[:]
        nreg     = len(chip_reg)
        multireg = True
        # check that all elements of the list are either strings or None:
        if [True for reg in chip_reg
            if reg is not None and not isinstance(reg, str)]:
            raise TypeError("Argument 'chip_reg' can be either None, " \
                            "a string, or a list of stings and/or None.")
    else:
        raise TypeError("Argument 'chip_reg' can be either None, " \
                        "a string, or a list of stings and/or None.")

    # from the 'images' argument built a list of file *names*
    # ignoring extensions(!) and check that the files exist:
    imgfnames = []
    if not images:
        raise TypeError("Argument 'images' cannot be an empty list or None.")

    if not isinstance(images, list):
        images = [images]

    for fname in images:
        if not isinstance(fname, str):
            raise RuntimeError("Argument 'images' must be either a string " \
                    "or a non-empty list of strings of valid file names.")
        (fn, ext) = extension_from_filename(fname)
        if not os.path.isfile(fn):
            raise IOError("The image file \'%s\' does not exist." % fn)
        imgfnames.append(fn)

    # Get the HDU list of the first file in the list of images. This will be
    # re-used to check available extensions.
    try:
        hdulist = fits.open(imgfnames[0], memmap=False)
        hdulist.close()
    except IOError as e:
        cmsg = "Unable to open the image file \'%s\'." % imgfnames[0]
        if e.args:
            e.args = (e.args[0] + "\n" + cmsg,) + e.args[1:]
        else:
            e.args=(cmsg,)
        raise e

    except:
        raise

    # initial processing of extensions (get a list of integers, str, or tuples):
    exts = parse_ext(img_wcs_ext, default_extver=None)
    if not isinstance(exts, list):
        exts = [exts]

    # for each string extension in the extension list 'exts', find out how many
    # extension versions exist and create a new list of extensions in which the
    # string (e.g., 'sci') gets replaced with tuples for each available
    # extension version: [('sci',1),('sci',2), etc.]
    extensions = []
    for ext in exts:
        if isinstance(ext,int) or isinstance(ext,tuple):
            extensions.append(ext)
            continue
        extv = get_extver_list(hdulist, extname=ext)
        for v in extv:
            extensions.append((ext,v))

    if not _check_FITS_extensions(hdulist, extensions):
        raise RuntimeError("Not all extensions computed based on the "         \
        "provided 'img_wcs_ext' could be found in the HDU list of the "        \
        "image(s) specified by the argument 'images'. Make sure that all "     \
        "extensions in the 'img_wcs_ext' list are present in input 'images'.")

    nexts = len(extensions) # <- number of "valid" extensions

    # Warn users if some regions will be dropped out:
    if nreg > nexts:
        print("")
        _print_warning("The number of region files provided through " \
                       "'chip_reg' is larger than the number of extensions " \
                       "derived based on 'img_wcs_ext' that were found " \
                       "in 'images'.")
        _print_warning("The following \"fixed\" region files will " \
                       "be dropped out:")
        for i in range(nexts,nreg):
            _print_warning("chip_reg[%d]: %s" % (i,str(chip_reg[i])))
        chip_reg = chip_reg[:nexts]

    # Warn users if there are more extensions than "chip" regions:
    if nreg == 0:
        chip_reg = nexts * [None]

    elif nreg < nexts:
        print("")
        _print_warning("The number of region files provided through " \
                       "'chip_reg' is smaller than the number of extensions " \
                       "derived based on 'img_wcs_ext' that were found " \
                       "in 'images'.")
        _print_warning("The following extensions have not been assigned any "  \
                       "\"fixed\" (chip related) regions:")
        extlist = "   "
        for ext in extensions[nreg:-1]:
            chip_reg.append(None)
            extlist += str(ext)+", "
        extlist += str(extensions[-1])
        chip_reg.append(None)
        _print_warning(extlist)

    # Check that the region files exist and try to open/read them:
    region_lists = []

    for fname in chip_reg:
        if fname is None:
            region_lists.append(None)
            continue
        # check region file existence:
        if not os.path.isfile(fname):
            raise IOError("The input \"chip\" region file \'%s\' does not "    \
                          "exist." % fname)
        # try to read regions:
        try:
            reglist = pyregion.open(fname)
            region_lists.append(reglist)
            #TODO: Check that regions are all in *image* coordinates

        except IOError:
            raise IOError("Unable to open the input region file \'%s\'." %     \
                          fname)

        except:
            raise RuntimeError("Unexpected error parsing input region "        \
                               "file \'%s\'. Suspecting either a corrupt "     \
                               "file or an unrecognizable region file format." \
                               % fname)

    # Print the list of pairs extension->chip region:
    if verbose and nreg > 0:
        print("")
        _print_important("Please verify that the following association "  \
              "between the commanded FITS extensions and provided \"fixed\" " \
              "region files is correct:")
        print("-----------------------------")
        print("EXTENSION:   -->   REGION FILE:")
        for i in range(nexts):
            print("{!s:<18.18s} {}".format(extensions[i], chip_reg[i]))

    # Return a list of image file names (with any extensions removed!) to which
    # regions will be mapped AND a list of pairs of the form
    # [region, extension].
    return imgfnames, list(map(list, list(zip(region_lists, extensions))))


def _print_warning(msg):
    print("\033[1m"+"WARNING: "+"\33[0m"+msg)

def _print_important(msg):
    print("\033[1mIMPORTANT: \33[0m{}".format(msg))

def _needs_ref_WCS(reglist):
    """ Check if the region list contains shapes in image-like coordinates
    """
    from pyregion.wcs_helper import image_like_coordformats

    for r in reglist:
        if r.coord_format in image_like_coordformats:
            return True
    return False


def _split_sky_img_regions(reglist):
    # Given an input ShapeList object, _split_sky_img_regions separates the
    # regions in the ShapeList object into the regions defined in sky coordinates
    # and regions defined in image coordinates.
    #
    # Return Value: a tuple of lists of regions (image, sky)
    from pyregion.wcs_helper import image_like_coordformats

    img = []
    sky = []

    for r in reglist:
        if r.coord_format in image_like_coordformats:
            img.append(r)
        else:
            sky.append(r)

    return (img, sky)


def extension_from_filename(filename):
    """
    Parse out filename from any specified extensions.
    Returns rootname and string version of extension name.
    """
    # Parse out any extension specified in filename
    _indx1 = filename.find('[')
    _indx2 = filename.find(']')
    if _indx1 > 0:
        # check for closing square bracket:
        if _indx2 < _indx1:
            raise RuntimeError("Incorrect extension specification in file " \
                                "name \'%s\'." % filename)
        # Read extension name provided
        _fname = filename[:_indx1]
        _extn = filename[_indx1+1:_indx2].strip()
    else:
        _fname = filename
        _extn = None

    return _fname, _extn


def parse_ext(extn, default_extver=None):

    if isinstance(extn, list):
        return [parse_ext(ext, default_extver) for ext in extn]

    if default_extver is not None and not isinstance(default_extver, int):
        raise TypeError("Argument 'default_extver' must be an integer or None.")

    if extn is None:
        return 0

    if isinstance(extn, int):
        return extn

    if isinstance(extn, tuple):
        if len(extn) < 1:
            return 0

        extname = extn[0]
        extver  = default_extver

        if not isinstance(extname, str):
            raise TypeError("The first element of a tuple must be a string " \
                "with the extension name ('EXTNAME').")

        if len(extn) == 1:
            return (extname, default_extver) if default_extver is not None \
                                             else extname

        if len(extn) > 1:
            if not (isinstance(extn[1], int) or extn[1] is None):
                raise TypeError("The second element of a tuple must be an " \
                    "integer number corresponding to the extension " \
                    "version ('EXTVER').")
            extver = default_extver if extn[1] is None else extn[1]

        return (extname if extver is None else (extname, extver))

    if not isinstance(extn, str):
        raise TypeError("Input extension name/number is not of any of the "\
                        "supported types: integer number, string, or tuple.")

    # check that extn is not a "forced" 'extname':
    forced_string_delim = ['\'','\"']
    extn = extn.strip()

    if extn and (extn[0]  in forced_string_delim) and \
                (extn[-1] in forced_string_delim):
        if default_extver is None:
            return extn[1:-1]
        else:
            return (extn[1:-1], default_extver)

    #TODO: Add support for strings like '(SCI,1)' that *include* parentheses?

    # check if extension version is present:
    commapos = extn.rfind(',')
    if commapos >= 0:
        # see if extension version is an integer:
        try:
            extver_str = extn[commapos+1:].strip()
            extver = int(extver_str) if extver_str else default_extver
        except ValueError as e:
            cmsg = "FITS extension version must be a valid integer."
            if e.args:
                e.args = (e.args[0] + "\n" + cmsg,) + e.args[1:]
            else:
                e.args=(cmsg,)
            raise e

        except:
            raise

        extname = extn[:commapos]

        return (extname if extver is None else (extname, extver))

    else:
        # no extver present. However, check if the 'extn' string can be
        # interpreted as an integer. If so, interpret it as extension version.
        try:
            extver = int(extn)
            return extver
        except:
            pass

        # assume default extver if present:
        if default_extver is None:
            return extn[:]
        else:
            return (extn[:], default_extver)
#
# Tests of 'parse_ext':
#
#--> parse_ext(2)
#2
#--> parse_ext(None)
#0
#--> parse_ext(('sci',2))
#('sci', 2)
#--> parse_ext(('sci'))
#'sci'
#--> parse_ext(['file1','file2'])
#['file1', 'file2']
#--> parse_ext('sci')
#'sci'
#--> parse_ext('sci,2')
#('sci', 2)
#--> parse_ext('sci,2.3')
#Traceback (innermost last):
#    File "<console>", line 1, in <module>
#    File "./tweakregtools.py", line 164, in parse_ext
#ValueError: invalid literal for int() with base 10: '2.3'
#FITS extension version must be a valid integer.


def count_extensions(img, extname='SCI'):
    """ Return the number of 'extname' extensions. 'img' can be either a file
    name, an HDU List object (from fits), or None (to get the number of all
    HDU headers.
    """
    if isinstance(img, str):
        img = fits.open(img, memmap=False)
        img.close()
    elif not isinstance(img, fits.HDUList):
        raise TypeError("Argument 'img' must be either a file name (string) " \
                        "or a `astropy.io.fits.HDUList` object.")

    if extname is None:
        return len(img)

    if not isinstance(extname, str):
        raise TypeError("Argument 'extname' must be either a string " \
        "indicating the value of the 'EXTNAME' keyword of the extensions " \
        "to be counted or None to return the count of all HDUs in the " \
        "'img' FITS file.")

    extname = extname.upper()

    n = 0
    for e in img:
        #if isinstance(e, fits.ImageHDU): continue
        if 'EXTNAME' in list(map(str.upper, list(e.header.keys()))) \
            and e.header['extname'].upper() == extname:
            n += 1

    return n


def get_extver_list(img, extname='SCI'):
    """ Return a list of all extension versions of 'extname' extensions.
    'img' can be either a file name or a HDU List object (from fits).
    """
    if isinstance(img, str):
        img = fits.open(img, memmap=False)
        img.close()
    elif not isinstance(img, fits.HDUList):
        raise TypeError("Argument 'img' must be either a file name (string) "  \
                        "or a fits.HDUList object.")

    # when extver is None - return the range of all FITS extensions
    if extname is None:
        extver = list(range(len(img)))
        return extver

    if not isinstance(extname, str):
        raise TypeError("Argument 'extname' must be either a string "          \
            "indicating the value of the 'EXTNAME' keyword of the extensions " \
            "whose versions are to be returned or None to return "             \
            "extension numbers of all HDUs in the 'img' FITS file.")

    extname = extname.upper()

    extver = []
    for e in img:
        #if not isinstance(e, fits.ImageHDU): continue
        hkeys = list(map(str.upper, list(e.header.keys())))
        if 'EXTNAME' in hkeys and e.header['EXTNAME'].upper() == extname:
            extver.append(e.header['EXTVER'] if 'EXTVER' in hkeys else 1)

    return extver


def _check_FITS_extvers(img, extname, extvers):
    """Returns True if all (except None) extension versions specified by the
    argument 'extvers' and that are of the type specified by the argument
    'extname' are present in the 'img' FITS file. Returns False if some of the
    extension versions for a given EXTNAME cannot be found in the FITS image.
    """
    default_extn = 1 if isinstance(extname, str) else 0

    if isinstance(extvers, list):
        extv = [default_extn if ext is None else ext for ext in extvers]
    else:
        extv = [default_extn if extvers is None else extvers]

    extv_in_fits = get_extver_list(img, extname)

    return set(extv).issubset(set(extv_in_fits))


def _check_FITS_extensions(img, extensions):
    """
    """
    if not isinstance(extensions, list):
        extensions = [extensions]

    int_ext  = []
    ext_dict = { }

    # sort 'extensions' elements into a list of integers or a dictionary having
    # as keys the first elements of the tuples and the values being lists of
    # the second elements of the tuples corresponding to a given key:
    for extn in extensions:

        if isinstance(extn, int):
            int_ext.append(extn)

        elif isinstance(extn, tuple):
            if len(extn) < 2 or not isinstance(extn[0], str) or \
                                not isinstance(extn[1], int):
                raise TypeError("Argument 'extensions' must be a list "     \
                    "containing a mixture of None, integers (FITS "         \
                    "extension numbers), strings (FITS EXTNAME), or tuple " \
                    "types (str, int) with the first element being the "    \
                    "EXTNAME and the second element being EXTVER.")
            if extn[0] in ext_dict:
                ext_dict[extn[0]].append(extn[1])
            else:
                ext_dict.update({extn[0]: [extn[1]]})

        elif isinstance(extn, str):
            if extn in ext_dict:
                ext_dict[extn].append(1)
            else:
                ext_dict.update({extn: [1]})

        elif extn is None:
            int_ext.append(0)

        else:
            raise TypeError("Argument 'extensions' must be a list "     \
                "containing a mixture of None, integers (FITS "         \
                "extension numbers), strings (FITS EXTNAME), or tuple " \
                "types (str, int) with the first element being the "    \
                "EXTNAME and the second element being EXTVER.")

    # if 'img' is a file name - open the FITS file:
    if isinstance(img, str):
        img = fits.open(img, memmap=False)
        img.close()
    elif not isinstance(img, fits.HDUList):
        raise TypeError("Argument 'img' must be either a file name (string) " \
                        "or a fits.HDUList object.")

    all_present = True

    if len(int_ext) > 0:
        if max(int_ext) >= len(img):
            all_present = False

    for extname in ext_dict.keys():
        if not _check_FITS_extvers(img, extname, ext_dict[extname]):
            all_present = False
            break

    return all_present


#--------------------------
# TEAL Interface functions
#--------------------------
def run(configObj):
    MapReg(input_reg   = configObj['input_reg'],
           images      = configObj['images'],
           img_wcs_ext = configObj['img_wcs_ext'],
           refimg      = '', #configObj['refimg'],
           ref_wcs_ext = 'sci', #configObj['ref_wcs_ext'],
           chip_reg    = configObj['chip_reg'],
           outpath     = configObj['outpath'],
           filter      = configObj['filter'],
           catfname    = configObj['catfname'],
           iteractive  = False, #configObj['iteractive'],
           append      = configObj['append'],
           verbose     = configObj['verbose'])


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


MapReg.__doc__ = getHelpAsString(docstring = True, show_ver = False)
