""" DrizzlePac - Python tasks to 'dither' combine images

This package provides the tools to align, manage and combine images using
algorithms originally developed as part of IRAF's 'dither' package. Tasks
supported by this package include:
  * astrodrizzle
  * tweakreg
  * resetbits
  * updatenpol
  * pixtosky, skytopix and pixtopix
  * tweakback

**Output**: The primary output from this task is the distortion-corrected,
cosmic-ray cleaned, and combined image as a FITS file.

"""
import contextlib
import importlib
from typing import TYPE_CHECKING

from .version import __version__

_OPTIONAL_IMPORTS = ("stsci.skypac",)

for _name in _OPTIONAL_IMPORTS:
    with contextlib.suppress(Exception):
        with contextlib.redirect_stdout(None):
            importlib.import_module(_name)

# Lazily import task modules so package import does not require compiled extensions.
_LAZY_SUBMODULES = {
        "ablot": "drizzlepac.ablot",
        "adrizzle": "drizzlepac.adrizzle",
        "align": "drizzlepac.align",
        "astrodrizzle": "drizzlepac.astrodrizzle",
        "buildmask": "drizzlepac.buildmask",
        "buildwcs": "drizzlepac.buildwcs",
        "catalogs": "drizzlepac.catalogs",
        "createMedian": "drizzlepac.createMedian",
        "drizCR": "drizzlepac.drizCR",
        "haputils": "drizzlepac.haputils",
        "imageObject": "drizzlepac.imageObject",
        "imagefindpars": "drizzlepac.imagefindpars",
        "imgclasses": "drizzlepac.imgclasses",
        "make_custom_mosaic": "drizzlepac.make_custom_mosaic",
        "mapreg": "drizzlepac.mapreg",
        "mdzhandler": "drizzlepac.mdzhandler",
        "outputimage": "drizzlepac.outputimage",
        "photeq": "drizzlepac.photeq",
        "pixreplace": "drizzlepac.pixreplace",
        "pixtopix": "drizzlepac.pixtopix",
        "pixtosky": "drizzlepac.pixtosky",
        "processInput": "drizzlepac.processInput",
        "refimagefindpars": "drizzlepac.refimagefindpars",
        "resetbits": "drizzlepac.resetbits",
        "runastrodriz": "drizzlepac.runastrodriz",
        "sky": "drizzlepac.sky",
        "skytopix": "drizzlepac.skytopix",
        "staticMask": "drizzlepac.staticMask",
        "tweakback": "drizzlepac.tweakback",
        "tweakreg": "drizzlepac.tweakreg",
        "tweakutils": "drizzlepac.tweakutils",
        "updatenpol": "drizzlepac.updatenpol",
        "util": "drizzlepac.util",
        "wcs_functions": "drizzlepac.wcs_functions",
}

if TYPE_CHECKING:  # pragma: no cover - aid static analyzers without runtime import
    from . import ablot as ablot  # noqa: F401
    from . import adrizzle as adrizzle  # noqa: F401
    from . import align as align  # noqa: F401
    from . import astrodrizzle as astrodrizzle  # noqa: F401
    from . import buildmask as buildmask  # noqa: F401
    from . import buildwcs as buildwcs  # noqa: F401
    from . import catalogs as catalogs  # noqa: F401
    from . import createMedian as createMedian  # noqa: F401
    from . import drizCR as drizCR  # noqa: F401
    from . import haputils as haputils  # noqa: F401
    from . import imageObject as imageObject  # noqa: F401
    from . import imagefindpars as imagefindpars  # noqa: F401
    from . import imgclasses as imgclasses  # noqa: F401
    from . import make_custom_mosaic as make_custom_mosaic  # noqa: F401
    from . import mapreg as mapreg  # noqa: F401
    from . import mdzhandler as mdzhandler  # noqa: F401
    from . import outputimage as outputimage  # noqa: F401
    from . import photeq as photeq  # noqa: F401
    from . import pixreplace as pixreplace  # noqa: F401
    from . import pixtopix as pixtopix  # noqa: F401
    from . import pixtosky as pixtosky  # noqa: F401
    from . import processInput as processInput  # noqa: F401
    from . import refimagefindpars as refimagefindpars  # noqa: F401
    from . import resetbits as resetbits  # noqa: F401
    from . import runastrodriz as runastrodriz  # noqa: F401
    from . import sky as sky  # noqa: F401
    from . import skytopix as skytopix  # noqa: F401
    from . import staticMask as staticMask  # noqa: F401
    from . import tweakback as tweakback  # noqa: F401
    from . import tweakreg as tweakreg  # noqa: F401
    from . import tweakutils as tweakutils  # noqa: F401
    from . import updatenpol as updatenpol  # noqa: F401
    from . import util as util  # noqa: F401
    from . import wcs_functions as wcs_functions  # noqa: F401


def __getattr__(name):
    if name in _LAZY_SUBMODULES:
        module = importlib.import_module(_LAZY_SUBMODULES[name])
        globals()[name] = module
        return module
    raise AttributeError(f"module 'drizzlepac' has no attribute '{name}'")


def __dir__():
    return sorted(set(globals()) | set(_LAZY_SUBMODULES))


def help():
    msg = \
""" The DrizzlePac package contains a suite of tasks that allow users to align HST images, combine them, and perform coordinate transformations on source positions.

drizzlepac:
       astrodrizzle - primary task for combining images, removing cosmic rays, and removing distortion
           tweakreg - task to compute offsets in WCS between images and a reference image or reference frame
      imagefindpars - sub-task containing parameters to find point sources used by tweakreg to build source catalogs for each tweakreg input image
          tweakback - apply an updated WCS solution created by tweakreg for a drizzled image to the constituent distorted (flt.fits) images
             mapreg - task to map a DS9 region file to multiple images based on the WCS information of each image.
           pixreplace - task to replace pixel values such as NaNs in images with another value
           pixtopix - task to convert pixel positions from an input image to pixel positions in an output WCS or image
           pixtosky - task to convert pixel positions from an input image to sky coordinates with full distortion correction as appropriate
           photeq   - task to equalize sensitivities of images across chips and epochs
           skytopix - task to convert sky positions to pixel positions in an image
          resetbits - sub-task to reset specified flt.fits data quality (DQ) values to 0
         updatenpol - task to add the names of the new ACS distortion reference files NPOLFILE and D2IMFILE then update headers to include residual distortion corrections as image extensions

fitsblender:
       blendheaders - task to merge the keywords from all images used to create a drizzle product into a single header with a table extension using rules defined for each instrument

stwcs:
    apply_headerlet - apply a headerlet to a file
  archive_headerlet - save a WCS solution as a headerlet extension and write it out as a headerlet FITS file
   attach_headerlet - attach a headerlet as an extension to a file
   delete_headerlet - delete a headerlet extension from a file
  extract_headerlet - write out a headerlet extension as a separate FITS file
  headerlet_summary - print a summary of all headerlet extensions in a file
  restore_headerlet - replace current WCS solution with the WCS solution from a headerlet extension
    write_headerlet - save a WCS solution as a separate headerlet FITS file
          updatewcs - recompute the WCS keywords and import the distortion model from the reference files
"""
    print(msg)
