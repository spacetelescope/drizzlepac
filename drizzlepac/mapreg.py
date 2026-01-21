"""
This module provides functions for mapping DS9 region files given in sky
coordinates to DS9 region files specified in image coordinates
of multiple images using the WCS information from the images.

:py:func:`~drizzlepac.mapreg.MapReg` provides an automated interface for converting
a region file to the image coordinate system (CS) of multiple images (and
their extensions) using WCS information from the image(s) header(s).
This conversion does not take into account pointing errors and, therefore,
an examination and adjustment (if required) of output region files is highly
recommended. This task is designed to simplify the creation of the exclusions
and/or inclusions region files used with :py:func:`~drizzlepac.tweakreg.TweakReg`
task for sources finding.

"""

import glob
import os
from astropy.io import fits
from astropy.wcs import WCS
from regions import Regions
from astropy.utils import deprecated

__all__ = ["MapReg", "map_region_files"]


def _warn_unused_parameter(name, value, default, since, message=None):
    if value == default:
        return

    note = message or "The '%s' parameter is deprecated and ignored." % name

    @deprecated(since=since, name=name, message=note, obj_type="parameter")
    def _deprecated_param():
        return None

    _deprecated_param()


def MapReg(input_reg, images, img_wcs_ext="sci", refimg="", ref_wcs_ext="sci", 
           chip_reg="", outpath="./regions", filter="", catfname="", 
           iteractive=False, append=False, verbose=True):
    """Primary interface to map DS9 region files given in sky coordinates.

    Parameters
    ----------
    input_reg : string or list of strings (Default = '')
        Input region files that need to be mapped to image CS using WCS information
        from ``images`` (see below). Only region files saved in sky CS are allowed
        in this release. Regions specified in image-like coordinates (e.g., image,
        physical) will be ignored.

        This paramater can be provided in any of several forms:

        * filename of a single image
        * comma-separated list of filenames
        * ``@-file`` filelist containing list of desired input region filenames

        The ``@-file`` filelist needs to be provided as an ASCII text file
        containing a list of filenames for all input region files with one
        filename on each line of the file.

    images : string or list of strings (Default = ``*.fits``)
        FITS images onto which the region files ``input_reg`` will be mapped. These
        image files must contain WCS information in their headers in order to
        convert ``input_reg`` from sky coordinates to correct image coordinates.

        This paramater can be provided in any of several forms:

        * filename of a single image
        * filename of an association (ASN)table
        * wild-card specification for files in a directory (using ``*``, ``?`` etc.)
        * comma-separated list of filenames
        * ``@-file`` filelist containing list of desired input filenames
          (and optional inverse variance map filenames)

        The ``@-file`` filelist needs to be provided as an ASCII text file
        containing a list of filenames for all input images (to which ``input_reg``
        regions should be mapped) with one filename on each line of the file.

    img_wcs_ext : string or list of strings (Default = ``SCI``)
        Extension name, extension name and version, or extension number of FITS
        extensions in the ``images`` to which the input regions ``input_reg`` should be
        mapped. The header of each extension must contain WCS information that will
        be used to convert ``input_reg`` from sky CS to image-like CS. Multiple
        extensions must be separated by semicolon while extension name and version
        (if present) must be separated by comma, e.g., ``'SCI;DQ,1;0'``.
        When specifying the extension name only, internally it will be expanded
        into a list of extension names and versions for each version of that
        extension name present in the input ``images``. For example, if a FITS file
        has four SCI and four DQ extensions, then ``'SCI;DQ,1;0'`` will be expanded into
        ``'SCI,1;SCI,2;SCI,3;SCI,4;DQ,1;0'``.
    
    outpath : string (Default = ``./regions``)
        The directory to which the transformed regions should be saved.

    filter : string {'None', 'fast', or 'precise' } (Default = 'None')
        Specify whether or not to remove the regions in the transformed region files
        that are outside the image array. With the 'fast' mehod only intersection
        of the bounding boxes is being checked.

    refimg : string (Default = '', deprecated)
        **Reserved for future use.**
        Filename of the reference image. May contain extension specifier:
        [extname,extver], [extname], or [extnumber].

    chip_reg : string or list of strings (Default = '', deprecated)
        Input region files in image CS associated with each extension specified by
        the ``img_wcs_ext`` parameter above. These regions will be added directly
        (without any transformation) to the ``input_reg`` regions mapped to each
        extension of the input ``images``. These regions must be specified in
        image-like coordinates. Typically, these regions should contain "exclude"
        regions to exclude parts of the image specific to the detector **chip**
        (e.g., vignetted regions due to used filters, or occulting finger in ACS/HRC
        images) from being used for source finding.

        This paramater can be provided in one of the following forms:

        * filename of a single image (if ``img_wcs_ext`` specifies a single FITS
          extension);
        * comma-separated list of filenames (if ``img_wcs_ext`` specifies more than
          one extension) or ``None`` for extensions that do not need any
          chip-specific regions to be excluded/included;
        * '' (empty string) or None if no chip-specific region files are provided.

        The number of regions ideally must be equal to the number of extensions
        specified by the ``img_wcs_ext`` parameter. If the number of chip-specific
        regions is less than the number of ``img_wcs_ext`` extensions then 'chip_reg'
        regions will be assigned to the first extensions from ``img_wcs_ext``
        (after internal expansion described in help for the ``img_wcs_ext``
        parameter above). If the number of 'chip_reg' is larger than the number of
        ``img_wcs_ext`` extensions then extra regions will be ignored.

    catfname : string (Default = ``exclusions_cat.txt``, deprecated)
        The file name of the output exclusions catalog file to be created from the
        supplied image and region file names. This file can be passed as an input
        to TweakReg task. Verify that the created file is correct!

    append : bool (Default = False, deprecated)
        Specify whether or not to append the transformed regions to the existing
        region files with the same name.

    iteractive : bool (Default = False, deprecated)
        **Reserved for future use.** (This switch controls whether the program stops
        and waits for the user to examine any generated region files before
        continuing on to the next image.)

    verbose : bool (Default = False, deprecated)
        Specify whether or not to print extra messages during processing.


    Notes
    -----
    **NOTE 1:** This task takes a region file (or multiple files) that describe(s)
    what regions of sky should be used for source finding (*include* regions) and
    what regions should be avoided (*exclude* regions) and transforms/maps
    this region file onto a number of image files that need to be aligned.

    The idea behind this task is automate the creation of region files that then can
    be passed to *exclusions* parameter of the ``TweakReg`` task.

    The same thing can be achieved manually using, for example, external FITS
    viewers, e.g., SAO DS9. For example, based on some image ``refimg.fits`` we can
    select a few small regions of sky that contain several good (bright, not
    saturated) point-like sources that could be used for image alignment of other
    images (say ``img1.fits``, ``img2.fits``, etc.). We can save this region file in
    sky coordinates (e.g., ``fk5``), e.g., under the name ``input_reg.reg``. We can then
    load a specific extension of each of the images ``img1.fits``, ``img2.fits``, etc.
    one by one into DS9 and then load onto those images the previously created
    include/exclude region file ``input_reg.reg``. Now we can save the regions using
    *image* coordinates. To do conversion from the sky coordinates to image
    coordinates, DS9 will use the WCS info from the image onto which the region file
    was loaded. The :py:func:`~drizzlepac.mapreg.MapReg` task tries to automate this process.

    **NOTE 2:** :py:func:`~drizzlepac.mapreg.MapReg` does not take into account pointing errors and thus the
    produced region files can be somewhat misaligned compared to their intended
    position around the sources identified in the "reference" image. Threfore, it is
    highly recommended that the produced region files be loaded into DS9 and their
    position be adjusted manually to include the sources of interest (or to avoid
    the regions that need to be avoided).
    If possible, the *include* or *exclude* regions should be large enough as to
    allow for most pointing errors.

    Examples
    --------
    Let's say that based on some image ``refimg.fits`` we have produced a "master"
    reference image (``master.reg``) that includes regions around sources that we want
    to use for image alignment in task :py:func:`~drizzlepac.tweakreg.TweakReg` and excludes regions that we
    want to avoid being used for image alignment (e.g, diffraction spikes, saturated
    quasars, stars, etc.). We save the file ``master.reg`` in sky CS (e.g., ``fk5``).

    Also, let's assume that we have a set of images ``img1.fits``, ``img2.fits``, etc.
    with four FITS extensions named 'SCI' and 'DQ'. For some of the extensions,
    after analyzing the ``img*.fits`` images we have identified parts of the chips that
    cannot be used for image alignment. We create region files for those extensions
    and save the files in image CS as, e.g., ``img1_chip_sci2.reg`` (in our example
    this will be the only chip that needs "special" treatment).

    Finally, let's say we want to "replicate" the "master" region file to all SCI
    extensions of the ``img*.fits`` images as well as to the 2nd DQ extension and to
    the 8th extension of the ``img*.fits`` images.

    To do this we run:

        >>> mapreg(input_reg = 'master.reg', images='img*.fits',
        ...        img_wcs_ext='sci;dq,2;8')

    This will produce six region files in the ./regions subdirectory for *each*
    input image::

        img1_sci1_twreg.reg,    img1_sci2_twreg.reg,    img1_sci3_twreg.reg,
        img1_sci4_twreg.reg,    img1_dq2_twreg.reg,     img1_extn8_twreg.reg
        ...

    ::

        img2_sci1_twreg.reg,    img2_sci2_twreg.reg,    img2_sci3_twreg.reg,
        img2_sci4_twreg.reg,    img2_dq2_twreg.reg,     img2_extn8_twreg.reg
        ...
    """

    _warn_unused_parameter("refimg", refimg, "", "3.11.0")
    _warn_unused_parameter("ref_wcs_ext", ref_wcs_ext, "sci", "3.11.0")
    _warn_unused_parameter("chip_reg", chip_reg, "", "3.11.0")
    _warn_unused_parameter("catfname", catfname, "", "3.11.0")
    _warn_unused_parameter("iteractive", iteractive, False, "3.11.0")
    _warn_unused_parameter("append", append, False, "3.11.0")
    _warn_unused_parameter("verbose", verbose, True, "3.11.0")

    map_region_files(
        input_reg,
        images=images,
        img_wcs_ext=img_wcs_ext,
        outpath=outpath,
        filter=filter,
    )


def _validate_input_reg(input_reg):
    """Normalize input region specifications into concrete file paths."""

    seen_lists = set()

    def _expand_string_spec(spec, context_dir=None):
        value = spec.strip()
        if not value:
            raise ValueError("input_reg entries must not be empty")

        if value.startswith("@"):
            list_path = value[1:].strip()
            if not list_path:
                raise ValueError("input_reg list specifications must name a file")
            if context_dir and not os.path.isabs(list_path):
                list_path = os.path.join(context_dir, list_path)
            list_path = os.path.normpath(list_path)
            if not os.path.isfile(list_path):
                raise IOError("The input region list '%s' does not exist." % list_path)

            key = os.path.abspath(list_path)
            if key in seen_lists:
                raise ValueError("input_reg list '%s' references itself" % list_path)

            seen_lists.add(key)
            entries = []
            with open(list_path, "r", encoding="utf-8") as handle:
                for raw_line in handle:
                    line = raw_line.strip()
                    if not line or line.startswith("#"):
                        continue
                    entries.extend(
                        _expand_string_spec(line, os.path.dirname(list_path))
                    )
            seen_lists.remove(key)

            if not entries:
                raise ValueError(
                    "input_reg list '%s' did not yield any filenames" % list_path
                )
            return entries

        if "," in value:
            parts = [part.strip() for part in value.split(",") if part.strip()]
            if not parts:
                raise ValueError(
                    "input_reg comma-separated specifications must list filenames"
                )
            expanded = []
            for part in parts:
                expanded.extend(_expand_string_spec(part, context_dir))
            return expanded

        if context_dir and not os.path.isabs(value):
            value = os.path.normpath(os.path.join(context_dir, value))

        return [value]

    def _expand_entry(entry):
        if isinstance(entry, str):
            return _expand_string_spec(entry)
        raise TypeError("input_reg entries must be strings")

    if isinstance(input_reg, str):
        candidates = _expand_string_spec(input_reg)
    elif isinstance(input_reg, (list, tuple)):
        candidates = []
        for item in input_reg:
            candidates.extend(_expand_entry(item))
    else:
        raise TypeError("input_reg must be a string or a sequence of strings")

    if not candidates:
        raise ValueError("input_reg must specify at least one region file")

    normalized = []
    for candidate in candidates:
        path = candidate.strip()
        if not path:
            raise ValueError("input_reg entries must not be empty")
        if not os.path.isfile(path):
            raise IOError("The input region file '%s' does not exist." % path)
        normalized.append(path)

    return normalized


def _validate_images(images):
    """Expand image specifications into a list of existing files."""

    seen_lists = set()

    def _normalize_path(path, base_dir):
        expanded = os.path.expanduser(path)
        if base_dir and not os.path.isabs(expanded):
            expanded = os.path.join(base_dir, expanded)
        return os.path.abspath(expanded)

    def _expand_string_spec(spec, base_dir=None):
        value = spec.strip()
        if not value:
            raise ValueError("images entries must not be empty")

        if value.startswith("@"):
            list_name = value[1:].strip()
            if not list_name:
                raise ValueError("images list specifications must name a file")
            list_path = _normalize_path(list_name, base_dir)
            if not os.path.isfile(list_path):
                raise IOError("The image list '%s' does not exist." % list_path)

            key = os.path.abspath(list_path)
            if key in seen_lists:
                raise ValueError("images list '%s' references itself" % list_path)

            seen_lists.add(key)
            collected = []
            with open(list_path, "r", encoding="utf-8") as handle:
                for raw_line in handle:
                    line = raw_line.strip()
                    if not line or line.startswith("#"):
                        continue
                    collected.extend(
                        _expand_string_spec(line, os.path.dirname(list_path))
                    )
            seen_lists.remove(key)

            if not collected:
                raise ValueError(
                    "images list '%s' did not yield any filenames" % list_path
                )
            return collected

        if "," in value:
            parts = [part.strip() for part in value.split(",") if part.strip()]
            if not parts:
                raise ValueError(
                    "images comma-separated specifications must list filenames"
                )
            expanded = []
            for part in parts:
                expanded.extend(_expand_string_spec(part, base_dir))
            return expanded

        normalized = _normalize_path(value, base_dir)

        if any(char in value for char in "*?[]"):
            matches = sorted(glob.glob(normalized))
            if not matches:
                raise IOError(
                    "The image wildcard '%s' did not match any files." % value
                )
            return matches

        if not os.path.isfile(normalized):
            raise IOError("The image file '%s' does not exist." % normalized)

        return [normalized]

    def _expand_entry(entry):
        if isinstance(entry, str):
            return _expand_string_spec(entry)
        raise TypeError("images entries must be strings")

    if isinstance(images, str):
        normalized = _expand_string_spec(images)
    elif isinstance(images, (list, tuple)):
        normalized = []
        for item in images:
            normalized.extend(_expand_entry(item))
    else:
        raise TypeError("images must be a string or a sequence of strings")

    if not normalized:
        raise ValueError("images must contain at least one file")

    return normalized


def _validate_img_wcs_ext(img_wcs_ext):
    if img_wcs_ext in (None, ""):
        return None

    def _normalize_single(ext):
        if isinstance(ext, str):
            value = ext.strip()
            if not value:
                raise ValueError("img_wcs_ext entries must not be empty strings")
            return value

        if isinstance(ext, int):
            if ext < 0:
                raise ValueError("img_wcs_ext integer entries must be non-negative")
            return ext

        if isinstance(ext, tuple):
            if len(ext) != 2:
                raise ValueError("img_wcs_ext tuple entries must be (extname, extver)")
            name, ver = ext
            if not isinstance(name, str) or not name.strip():
                raise TypeError(
                    "img_wcs_ext tuple entries must start with a non-empty string"
                )
            if ver is None or not isinstance(ver, int):
                raise TypeError("img_wcs_ext tuple extver must be an integer")
            if ver < 0:
                raise ValueError("img_wcs_ext tuple extver must be non-negative")
            return (name.strip(), ver)

        raise TypeError(
            "img_wcs_ext must contain strings, integers, or (str, int) tuples"
        )

    if isinstance(img_wcs_ext, (list, tuple)):
        normalized = [_normalize_single(ext) for ext in img_wcs_ext]
        if not normalized:
            raise ValueError("img_wcs_ext sequences must not be empty")
        return normalized

    return _normalize_single(img_wcs_ext)


def _validate_outpath(outpath):
    if outpath in (None, ""):
        return os.path.curdir + os.path.sep

    if not isinstance(outpath, str):
        raise TypeError("outpath must be a string")

    path = outpath.strip()
    if not path:
        path = os.path.curdir + os.path.sep

    if not os.path.isdir(path):
        raise IOError("The output directory '%s' does not exist." % path)

    return path


def _validate_filter(filter_option):
    if filter_option is None:
        return None

    if isinstance(filter_option, str):
        value = filter_option.strip().lower()
        if not value or value == "none":
            return None
        if value == "precise":
            print("\"precise\" filter option is not supported. Using 'fast' instead.")
            return "fast"
        if value == "fast":
            return "fast"
        raise ValueError("filter must be None, 'fast', 'precise', or 'none'")

    raise TypeError("filter must be None or a string value")


def _region_in_image(region, wcs, shape, mode="fast"):
    """
    Check whether a region is within an image.

    Parameters
    ----------
    region : regions.Region
        SkyRegion or PixelRegion.
    wcs : astropy.wcs.WCS
        WCS for the image extension.
    shape : tuple
        Image shape as (ny, nx).
    mode : {"fast", "precise"}
        - "fast": bounding-box containment (approximate, very fast)
        - "precise": pixel-mask containment (exact, slower)

    Returns
    -------
    inside : bool
        True if region is fully within the image.
    """

    # Convert to pixel coordinates if needed
    pix = region.to_pixel(wcs) if hasattr(region, "to_pixel") else region

    ny, nx = shape

    if mode == "fast":
        bb = pix.bounding_box
        return bb.ixmin >= 0 and bb.iymin >= 0 and bb.ixmax <= nx and bb.iymax <= ny

    elif mode == "precise":
        mask = pix.to_mask()
        img = mask.to_image((ny, nx))
        return img is not None and img.all()

    else:
        raise ValueError("mode must be 'fast' or 'precise'")


def _ext2str_suffix(ext):
    if isinstance(ext, tuple):
        return "_{}{}_twreg".format(ext[0], ext[1])
    elif isinstance(ext, int):
        return "_extn{}_twreg".format(ext)
    else:
        return "_{}_twreg".format(ext)  # <- we should not get here...


def map_region_files(input_reg, images, img_wcs_ext=None, outpath="./regions", 
                     filter="fast"):
    """Map DS9 sky-coordinate regions onto image coordinates.

    Parameters
    ----------
    input_reg : str or sequence of str
        Region file specification(s) written in a sky-based coordinate system.
        Supports single filenames, comma-separated filenames, wildcard
        patterns, and iterable collections of filenames.
    images : str or sequence of str
        FITS image specification(s) providing the WCS transformations. Strings
        may contain single filenames, comma-separated filenames, wildcard
        patterns, or ``@`` list files. Iterable inputs can mix these forms.
    img_wcs_ext : str, int, tuple, or sequence, optional
        Extension selector(s) identifying which HDUs receive mapped regions.
        Strings reference EXTNAME values, integers reference extension
        numbers, and ``(extname, extver)`` tuples point to specific versions.
        The default ``'sci'`` expands to every science extension found in each
        image.
    outpath : str, optional
        Destination directory for the generated region files. Defaults to
        ``./regions`` and must exist ahead of time.
    filter : {None, 'fast', 'precise', 'none'}, optional
        Region filtering strategy. ``'fast'`` performs bounding-box checks,
        ``'precise'`` is accepted but coerced to ``'fast'``, and ``None``/``'none'``
        skips filtering entirely.

    Raises
    ------
    IOError
        If the input region file, any image file, or the output directory is
        missing.
    TypeError
        If an argument is not the expected type.
    ValueError
        If an argument is syntactically valid but empty or otherwise invalid.
    """

    region_paths = _validate_input_reg(input_reg)
    image_list = _validate_images(images)
    ext_spec = _validate_img_wcs_ext(img_wcs_ext)
    output_dir = _validate_outpath(outpath)
    filter_mode = _validate_filter(filter)

    sky_regions = []
    for region_path in region_paths:
        sky_regions.extend(Regions.read(region_path))

    all_sky_regions = Regions(sky_regions)
    for image_fname in image_list:
        current_exts = ext_spec

        if current_exts is None or current_exts == "sci":
            from drizzlepac.util import count_sci_extensions

            current_exts = count_sci_extensions(image_fname, return_ind=True)

        if not isinstance(current_exts, list):
            current_exts = [current_exts]

        with fits.open(image_fname, memmap=False) as hdu:
            for ext in current_exts:
                all_pix_region = []
                data_hdu = hdu[ext]
                wcs = WCS(data_hdu.header, hdu)
                shape = data_hdu.data.shape if data_hdu.data is not None else None

                for region in all_sky_regions:
                    if filter_mode and shape is not None:
                        if not _region_in_image(
                            region, wcs, shape, mode=filter_mode
                        ):
                            continue
                    pix_reg = region.to_pixel(wcs)
                    all_pix_region.append(pix_reg)

                full_pix_region = Regions(all_pix_region)

                # get output region file name
                extsuffix = _ext2str_suffix(ext)
                basefname, _ = os.path.splitext(os.path.basename(image_fname))
                regfname = basefname + extsuffix + os.extsep + "reg"
                fullregfname = os.path.join(output_dir, regfname)

                full_pix_region.write(fullregfname, overwrite=True)
