.. _skytopix:

********************************************************
skytopix: Coordinate transformation from sky coordinates
********************************************************
This task allows a user to perform coordinate transformations
with the full WCS and distortion model on source positions
from sky coordinates to the WCS defined by an image.
This task serves as a replacement for
the ``IRAF.STSDAS`` task `rd2xy`, albeit with the
added capability of understanding the full distortion model
provided in the headers of images updated for use with astrodrizzle and
tweakreg.

.. moduleauthor:: Warren Hack <help@stsci.edu>

.. automodule:: drizzlepac.skytopix

.. autofunction:: rd2xy
