"""
This module provides functions for mapping DS9 region files given in sky
coordinates to DS9 region files specified in image coordinates of multiple
images using the WCS information from the images.

:py:func:`~drizzlepac.mapreg.MapReg` provides an automated interface for
converting a region file to the image coordinate system (CS) of multiple
images (and their extensions) using WCS information from the image(s)
header(s). This conversion does not take into account pointing errors and,
therefore, an examination and adjustment (if required) of output region
files is highly recommended. This task is designed to simplify the creation
of the exclusions and/or inclusions region files used with
:py:func:`~drizzlepac.tweakreg.TweakReg` task for sources finding.

"""

import os
import sys
import glob
import logging
import warnings
from astropy.io import fits
from astropy.wcs import WCS
from astropy.utils.exceptions import AstropyDeprecationWarning
from regions import Regions
from drizzlepac import __version__
from drizzlepac.util import count_sci_extensions

__all__ = ["MapReg"]
package_level_logger = logging.getLogger('drizzlepac')
log = logging.getLogger(f'drizzlepac.mapreg')


def MapReg(input_reg, images, img_wcs_ext="sci", outpath="./regions", filter=None,
           refimg="", chip_reg="", catfname="", append=False, ref_wcs_ext=None,
           iteractive=None, verbose=None):
    """Primary interface to map DS9 region files given in sky coordinates.

    Parameters
    ----------
    input_reg : string or list of strings (Default = '')
        Input region files that need to be mapped to image CS using WCS
        information from ``images`` (see below). Only region files saved in
        sky CS are allowed in this release. Regions specified in image-like
        coordinates (e.g., image, physical) will be ignored.

        This parameter can be provided in any of several forms:

        * filename of a single image
        * comma-separated list of filenames
        * ``@-file`` filelist containing list of desired input region filenames

        The ``@-file`` filelist needs to be provided as an ASCII text file
        containing a list of filenames for all input region files with one
        filename on each line of the file.

    images : string or list of strings (Default = ``*.fits``)
        FITS images onto which the region files ``input_reg`` will be mapped.
        These image files must contain WCS information in their headers to
        convert ``input_reg`` from sky coordinates to correct image coordinates.

        This parameter can be provided in any of several forms:

        * filename of a single image
        * filename of an association (ASN) table
        * wild-card specification for files in a directory (using ``*``, ``?``
            etc.)
        * comma-separated list of filenames
        * ``@-file`` filelist containing list of desired input filenames (and
          optional inverse variance map filenames)

        The ``@-file`` filelist needs to be provided as an ASCII text file
        containing a list of filenames for all input images (to which
        ``input_reg`` regions should be mapped) with one filename on each
        line of the file.

    img_wcs_ext : string or list of integers (Default = ``sci``)
        Extension selection for mapping. A string specifies an EXTNAME to use
        for every matching extension in each image (for example ``"SCI"``). A
        list of integers specifies explicit extension numbers to process. The
        associated headers must contain valid WCS information for the
        transformation.

        *Note:* Previously a tuple of strings was accepted to specify multiple
        EXTNAMEs; this usage is no longer supported.

    outpath : string (Default = ``./regions``)
        The directory to which the transformed regions should be saved. If the
        directory does not exist it will be created.

    filter : string {None, 'fast', 'precise', 'none'} (Default = None)
        Specify whether or not to remove the regions in the transformed region
        files that are outside the image array. With the 'fast' method only
        intersection of the bounding boxes is being checked.

    chip_reg : string or list of strings (Default = '', deprecated)
        Input region files in image CS associated with each extension specified
        by the ``img_wcs_ext`` parameter above. These regions will be added
        directly (without any transformation) to the ``input_reg`` regions
        mapped to each extension of the input ``images``. These regions must be
        specified in image-like coordinates. Typically, these regions should
        contain "exclude" regions to exclude parts of the image specific to the
        detector **chip** (e.g., vignetted regions due to used filters, or
        occulting finger in ACS/HRC images) from being used for source finding.

        This parameter can be provided in one of the following forms:

        * filename of a single image (if ``img_wcs_ext`` specifies a single
            FITS extension)
        * comma-separated list of filenames (if ``img_wcs_ext`` specifies more
            than one extension) or ``None`` for extensions that do not need any
            chip-specific regions to be excluded/included
        * '' (empty string) or None if no chip-specific region files are
            provided

        The number of regions ideally must be equal to the number of
        extensions specified by the ``img_wcs_ext`` parameter. If the number
        of chip-specific regions is less than the number of ``img_wcs_ext``
        extensions then ``chip_reg`` regions will be assigned to the first
        extensions from ``img_wcs_ext`` (after internal expansion described in
        help for the ``img_wcs_ext`` parameter above). If the number of
        ``chip_reg`` entries is larger than the number of ``img_wcs_ext``
        extensions then extra regions will be ignored.

    catfname : string (Default = ``exclusions_cat.txt``)
        The file name of the output exclusions catalog file to be created
        from the supplied image and region file names. This file can be
        passed as an input to TweakReg task. Verify that the created file is
        correct!

    append : bool (Default = False)
        Specify whether or not to append the transformed regions to the
        existing region files with the same name.

    verbose : bool (Default = False)
        Specify whether or not to print extra messages during processing.

    refimg : string (Default = '', deprecated)
        **Reserved for future use.** Filename of the reference image. May
        contain extension specifier: [extname,extver], [extname], or
        [extnumber].

    ref_wcs_ext : string (Default = 'sci', deprecated)
        **Reserved for future use.** Extension name and/or version of FITS
        extensions in the ``refimg`` that contain WCS information that will
        be used to convert ``input_reg`` from image-like CS to sky CS. NOTE:
        Only extension name is allowed when ``input_reg`` is a list of region
        files that contain regions in image-like CS. In this case, the number
        of regions in ``input_reg`` must agree with the number of extensions
        with name specified by ``ref_wcs_ext`` present in the ``refimg`` FITS
        image.

    iteractive : bool (Default = False, deprecated)
        **Reserved for future use.** This switch controls whether the
        program stops and waits for the user to examine any generated region
        files before continuing on to the next image.

    Notes
    -----
    **NOTE 1:** This task takes a region file (or multiple files) that
    describe(s) what regions of sky should be used for source finding
    (*include* regions) and what regions should be avoided (*exclude*
    regions) and transforms/maps this region file onto a number of image
    files that need to be aligned.

    The idea behind this task is automate the creation of region files that
    then can be passed to *exclusions* parameter of the ``TweakReg`` task.

    The same thing can be achieved manually using, for example, external
    FITS viewers, e.g., SAO DS9. For example, based on some image
    ``refimg.fits`` we can select a few small regions of sky that contain
    several good (bright, not saturated) point-like sources that could be
    used for image alignment of other images (say ``img1.fits``,
    ``img2.fits``, etc.). We can save this region file in sky coordinates
    (e.g., ``fk5``) under the name ``input_reg.reg``. We can then load a
    specific extension of each of the images ``img1.fits``, ``img2.fits``,
    etc. one by one into DS9 and then load onto those images the previously
    created include/exclude region file ``input_reg.reg``. Now we can save the
    regions using *image* coordinates. To do conversion from the sky
    coordinates to image coordinates, DS9 will use the WCS info from the
    image onto which the region file was loaded. The
    :py:func:`~drizzlepac.mapreg.MapReg` task tries to automate this process.

    **NOTE 2:** :py:func:`~drizzlepac.mapreg.MapReg` does not take into account
    pointing errors and thus the produced region files can be somewhat
    misaligned compared to their intended position around the sources
    identified in the "reference" image. Therefore, it is highly
    recommended that the produced region files be loaded into DS9 and their
    position be adjusted manually to include the sources of interest (or to
    avoid the regions that need to be avoided). If possible, the *include*
    or *exclude* regions should be large enough as to allow for most pointing
    errors.

    Examples
    --------
    Let's say that based on some image ``refimg.fits`` we have produced a
    "master" reference image (``master.reg``) that includes regions around
    sources that we want to use for image alignment in task
    :py:func:`~drizzlepac.tweakreg.TweakReg` and excludes regions that we
    want to avoid being used for image alignment (e.g, diffraction spikes,
    saturated quasars, stars, etc.). We save the file ``master.reg`` in sky
    CS (e.g., ``fk5``).

    Also, let's assume that we have a set of images ``img1.fits``,
    ``img2.fits``, etc. with four FITS extensions named 'SCI' and 'DQ'. For
    some of the extensions, after analyzing the ``img*.fits`` images we have
    identified parts of the chips that cannot be used for image alignment.
    We create region files for those extensions and save the files in image
    CS as, e.g., ``img1_chip_sci2.reg`` (in our example this will be the only
    chip that needs "special" treatment).

    Finally, let's say we want to "replicate" the "master" region file to all
    SCI extensions of the ``img*.fits`` images.

    To do this we run:

        >>> mapreg(input_reg='master.reg', images='img*.fits', img_wcs_ext=[2,8])

    This will produce region files in the ``./regions`` subdirectory for *each*
    input image::

        img1_sci2_twreg.reg,    img1_sci8_twreg.reg,

    ::

        img2_sci2_twreg.reg,    img2_sci8_twreg.reg,

    """
    
    # add logging
    if verbose:
        default_log_level = logging.DEBUG
        formatter = logging.Formatter('[%(levelname)s:%(name)s] %(message)s')
    else:
        default_log_level = logging.INFO
        formatter = logging.Formatter('[%(levelname)-8s] %(message)s')

    file_handler = logging.FileHandler('mapreg.log', mode='w', encoding='utf-8')
    stream_handler = logging.StreamHandler(sys.stdout)

    file_handler.setLevel(default_log_level)
    stream_handler.setLevel(default_log_level)
    file_handler.setFormatter(formatter)
    stream_handler.setFormatter(formatter)
    package_level_logger.addHandler(file_handler)
    package_level_logger.addHandler(stream_handler)
    package_level_logger.setLevel(default_log_level)

    # deprecation warnings for removed parameters
    if refimg not in ("", None):
        warnings.warn("The 'refimg' parameter is deprecated and ignored.",
            AstropyDeprecationWarning)

    if ref_wcs_ext not in ("sci", "SCI", None, ""):
        warnings.warn("The 'ref_wcs_ext' parameter is deprecated and ignored.",
            AstropyDeprecationWarning)

    if chip_reg not in ("", None):
        warnings.warn("The 'chip_reg' parameter is no longer supported; please"
            "make individual MapReg calls for each image/extension pair.",
            AstropyDeprecationWarning)

    if iteractive:
        warnings.warn("The 'iteractive' parameter is deprecated and ignored.",
            AstropyDeprecationWarning)

    # log input parameters
    log.info("Starting MapReg...")
    log.info("Mapping region files to image coordinate systems...")
    log.info(f"Using drizzlepac version: {__version__}")
    log.debug(f"Input region file(s): {input_reg}")
    log.debug(f"Input image file(s): {images}")
    log.debug(f"Image WCS extension(s): {img_wcs_ext}")
    log.debug(f"Output directory: {outpath}")
    log.debug(f"Output exclusions catalog file: {catfname}")
    log.debug(f"Region filtering method: {filter if filter else 'None'}")
    log.debug(f"Append mode: {append}")
    log.debug(f"Verbose mode: {verbose}")
    
    # validate and normalize inputs
    region_paths = tuple(_validate_input_reg(input_reg))
    image_list = tuple(_validate_images(images))
    ext_spec = _validate_img_wcs_ext(img_wcs_ext)
    normalized_catfname = _validate_catfname(catfname)
    output_dir = _validate_outpath(outpath)
    filter_mode = _validate_filter(filter)
    append_flag = _validate_append(append)

    # log updated parameters
    log.debug(f"Expanded region paths: {region_paths}")
    log.debug(f"Expanded image list: {image_list}")
    log.debug(f"Normalized img_wcs_ext: {ext_spec}")
    log.debug(f"Normalized outpath: {output_dir}")
    log.debug(f"Normalized catfname: {normalized_catfname}")
    log.debug(f"Normalized filter: {filter_mode}")
    log.debug(f"Normalized append: {append_flag}")

    # normalize extension specification
    if isinstance(ext_spec, list):
        normalized_ext_spec = ext_spec
    elif isinstance(ext_spec, tuple):
        normalized_ext_spec = tuple(ext_spec)
    else:
        normalized_ext_spec = ext_spec

    # read all input region files and combine into a single Regions object
    sky_regions = []
    for region_path in region_paths:
        sky_regions.extend(Regions.read(region_path))
    all_sky_regions = Regions(sky_regions)
    log.debug("Total number of sky regions: %d", len(all_sky_regions))

    # perform the mapping
    map_region_files(
        sky_regions=all_sky_regions,
        images=image_list,
        img_wcs_ext=normalized_ext_spec,
        catfname=normalized_catfname,
        outpath=output_dir,
        filter=filter_mode,
        append=append_flag,
    )

    # remove logging handlers
    package_level_logger.removeHandler(file_handler)
    package_level_logger.removeHandler(stream_handler)
    file_handler.close()
    stream_handler.close()


def _validate_input_reg(input_reg):
    """Normalize input region specifications into concrete file paths."""
    log.debug(f"Validating input_reg: {input_reg}")
    seen_lists = set()
    
    def _expand_string_spec(spec, context_dir=None):
        value = spec.strip()
        if not value:
            raise ValueError("input_reg entries must not be empty")

        if value.startswith("@-"):
            list_name = value[2:].strip()
            if not list_name:
                raise ValueError("input_reg list specifications must name a file")
            if context_dir and not os.path.isabs(list_name):
                list_name = os.path.join(context_dir, list_name)
            list_path = os.path.normpath(os.path.expanduser(list_name))
            if not os.path.isfile(list_path):
                raise IOError(f"The input region list '{list_path}' does not exist.")

            key = os.path.abspath(list_path)
            if key in seen_lists:
                raise ValueError(f"input_reg list '{list_path}' references itself")

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
                    f"input_reg list '{list_path}' did not yield any filenames"
                )
            return entries

        if value.startswith("@"):
            raise ValueError(
                "input_reg list specifications must use the '@-' prefix"
            )

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

        normalized = os.path.expanduser(value)
        if context_dir and not os.path.isabs(normalized):
            normalized = os.path.join(context_dir, normalized)
        normalized = os.path.normpath(normalized)

        if any(char in normalized for char in "*?[]"):
            matches = sorted(glob.glob(normalized))
            if not matches:
                raise IOError(
                    f"The input region wildcard '{value}' did not match any files."
                )
            return matches

        return [normalized]

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
            raise IOError(f"The input region file '{path}' does not exist.")
        normalized.append(path)

    return normalized


def _validate_images(images):
    """Expand image specifications into a list of existing files."""
    log.debug(f"Validating images: {images}")
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

        if value.startswith("@-"):
            list_name = value[2:].strip()
            if not list_name:
                raise ValueError("images list specifications must name a file")
            list_path = _normalize_path(list_name, base_dir)
            if not os.path.isfile(list_path):
                raise IOError(f"The image list '{list_path}' does not exist.")

            key = os.path.abspath(list_path)
            if key in seen_lists:
                raise ValueError(f"images list '{list_path}' references itself")

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
                    f"images list '{list_path}' did not yield any filenames"
                )
            return collected

        if value.startswith("@"):
            raise ValueError("images list specifications must use the '@-' prefix")

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
                    f"The image wildcard '{value}' did not match any files."
                )
            return matches

        if not os.path.isfile(normalized):
            raise IOError(f"The image file '{normalized}' does not exist.")

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
    log.debug(f"Validating img_wcs_ext: {img_wcs_ext}")
    
    if img_wcs_ext in (None, ""):
        return None

    def _normalize_single(ext):
        if not isinstance(ext, int):
            raise TypeError("img_wcs_ext sequences must contain integers")
        if ext < 0:
            raise ValueError("img_wcs_ext integers must be non-negative")
        return ext

    if isinstance(img_wcs_ext, str):
        value = img_wcs_ext.strip()
        if not value:
            raise ValueError("img_wcs_ext string must not be empty")
        if any(sep in value for sep in ";,"):
            raise ValueError("img_wcs_ext string must name a single"+
                             " EXTNAME without separators")
        return value

    if isinstance(img_wcs_ext, (list, tuple)):
        if not img_wcs_ext:
            raise ValueError("img_wcs_ext sequences must not be empty")
        return [_normalize_single(ext) for ext in img_wcs_ext]

    raise TypeError("img_wcs_ext must be a string or a sequence of integers")


def _validate_outpath(outpath):
    log.debug(f"Validating outpath: {outpath}")
    
    if outpath in (None, ""):
        return os.path.curdir + os.path.sep

    if not isinstance(outpath, str):
        raise TypeError("outpath must be a string")

    path = outpath.strip()
    if not path:
        path = os.path.curdir + os.path.sep

    if not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)
        log.debug(f"Created output directory: {path}")

    return path

def _validate_catfname(catfname):
    log.debug(f"Validating catfname: {catfname}")
    
    if catfname in (None, ""):
        return "exclusions_cat.txt"

    if not isinstance(catfname, str):
        raise TypeError("catfname must be a string")

    path = catfname.strip()
    if not path:
        path = "exclusions_cat.txt"

    return path


def _validate_filter(filter_option):
    log.debug(f"Validating filter_option: {filter_option}")
    
    if filter_option is None:
        return None

    if isinstance(filter_option, str):
        value = filter_option.strip().lower()
        if not value or value == "none" or value == "":
            return None
        if value == "fast":
            return "fast"
        if value == "precise":
            return "precise"
        raise ValueError("filter must be None, 'fast', 'precise', or 'none'")

    raise TypeError("filter must be None or a string value")


def _validate_append(append_option):
    log.debug(f"Validating append_option: {append_option}")
    if isinstance(append_option, bool):
        return append_option
    raise TypeError("append must be a boolean value")


def map_region_files(sky_regions, images, img_wcs_ext=None, outpath="./regions",
                     catfname="exclusions_cat.txt", filter=None, append=False):
    """Map DS9 sky-coordinate regions onto image coordinates.

    Parameters
    ----------
    sky_regions : regions.Regions
        Regions defined in sky coordinates to project into each output image.
        Callers must supply a populated :class:`regions.Regions` object.
    images : tuple[str, ...]
        Tuple of FITS image paths providing WCS transformations. Each file
        must exist on disk prior to invocation.
    img_wcs_ext : str | tuple[int, ...] | None, optional
        Normalized extension selector from :func:`_validate_img_wcs_ext`.
        ``"sci"`` expands to all science extensions, a tuple selects explicit
        HDU indices, and ``None`` implies automatic SCI discovery.
    outpath : str, optional
        Destination directory for the generated region files. The directory
        must already exist and should be provided by :func:`_validate_outpath`.
    catfname : str, optional
        Destination path for the exclusion catalog output. Supply the value
        returned by :func:`_validate_catfname`.
    filter : None | {'fast', 'precise'}, optional
        Filtering mode produced by :func:`_validate_filter`. ``None`` disables
        filtering.
    append : bool, optional
        Whether to append to existing outputs. Use the result from
        :func:`_validate_append`.
    """

    created_region_files = []

    # iterate of images and extensions
    for image_fname in images:
        image_label = os.path.basename(image_fname)
        current_exts = img_wcs_ext

        # expand "sci" to all science extensions
        if current_exts is None or (isinstance(current_exts, str) and 
                                    current_exts.lower() == "sci"):
            current_exts = count_sci_extensions(image_fname, return_ind=True)
            log.debug(f"Expanded 'sci' to extensions: {current_exts}")

        # make sure that we are working with a list of extensions and not a single value
        if not isinstance(current_exts, list):
            current_exts = [current_exts]

        with fits.open(image_fname, memmap=False) as hdu:
            for ext in current_exts:
                all_pix_region = []
                data_hdu = hdu[ext]
                wcs = WCS(data_hdu.header, hdu)
                shape = data_hdu.data.shape if data_hdu.data is not None else None

                # Project each sky region into pixel coordinates and filter if needed
                for region in sky_regions:
                    if filter and shape is not None:
                        # skip regions outside the image
                        if not _check_if_region_in_image(region, wcs, shape, mode=filter):
                            log.debug(f"Skipping region outside image: {region}")
                            continue
                    pix_reg = region.to_pixel(wcs)
                    all_pix_region.append(pix_reg)
                full_pix_region = Regions(all_pix_region)

                # do not create empty region files
                if not full_pix_region:
                    log.warning(f"No overlap between image '{image_label}' and "
                        f"extension {ext}; skipping file creation.")
                    continue
                
                # get output region file name
                extsuffix = _ext2str_suffix(ext)
                basefname, _ = os.path.splitext(os.path.basename(image_fname))
                regfname = basefname + extsuffix + os.extsep + "reg"
                fullregfname = os.path.join(outpath, regfname)
                if append:
                    _save_region_file(full_pix_region, fullregfname, 
                                      exclusion_file=False, append=append)
                else:
                    full_pix_region.write(fullregfname, format="ds9", overwrite=True)
                log.debug(f"Wrote region file: {fullregfname}")
                created_region_files.append(os.path.basename(fullregfname))
                
                # write exclusion catalog entry
                _save_region_file(full_pix_region, catfname, exclusion_file=True, 
                                  append=append)
                log.debug(f"Wrote exclusion catalog entry: {catfname}")
    if created_region_files:
        log.info(f"Output region files: {', '.join(created_region_files)}")
    else:
        log.info("Output region files: none")
    log.info("MapReg processing complete.")


def _check_if_region_in_image(region, wcs, shape, mode="fast"):
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


def _save_region_file(region, filename, exclusion_file, append):
    """Export a ds9 region file with optional appending and exclusion markers.

    When ``append`` is True, new shapes append to an existing file; otherwise
    the file is overwritten. When ``exclusion_file`` is True, each shape
    gains the '-' prefix required for tweakreg exclusion regions. For simple
    saves prefer the ``Regions.write`` method.

    Parameters
    ----------
    region : regions.Regions
        The regions to be converted.
    filename : str
        The output exclusion catalog filename.
    exclusion_file : bool
        Essentially, adds a '-' before each region shape type if True.
    append : bool
        Whether to append to the output file if it exists.
    """
    
    allowable_ds9_region_types = ["polygon", "circle", "ellipse", "box"]
    if len(region) == 0:
        raise ValueError("No regions were provided for exclusion conversion")

    header_lines = []
    header_captured = False
    converted_lines = []
    for reg in region:
        reg_serial = reg.serialize(format="ds9")
        shape_lines = []
        candidate_header = []
        saw_shape = False

        for raw_line in reg_serial.splitlines():
            stripped = raw_line.strip()

            if not stripped:
                continue

            if stripped.startswith("#"):
                if not header_captured and not saw_shape:
                    candidate_header.append(stripped)
                continue

            token = stripped.lstrip("-")
            lowered = token.lower()

            is_shape = False
            for region_type in allowable_ds9_region_types:
                if lowered.startswith(region_type):
                    is_shape = True
                    break

            if is_shape:
                if exclusion_file and not stripped.startswith("-"):
                    stripped = "-" + stripped
                shape_lines.append(stripped)
                saw_shape = True
                continue

            if not header_captured and not saw_shape:
                candidate_header.append(stripped)

        if not shape_lines:
            raise ValueError(f"No DS9 shapes found in region: {reg}")

        converted_lines.extend(shape_lines)
        if not header_captured and candidate_header:
            header_lines = candidate_header
            header_captured = True

    if filename:
        writing_mode = "a" if append else "w"
        write_header = True

        if append and os.path.exists(filename):
            write_header = os.path.getsize(filename) == 0

        with open(filename, writing_mode, encoding="utf-8") as handle:
            if write_header and header_lines:
                for header_line in header_lines:
                    handle.write(header_line + "\n")
            for entry in converted_lines:
                handle.write(entry + "\n")
        log.debug(f"Wrote exclusion catalog file: {filename}")

