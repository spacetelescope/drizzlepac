.. _pixtosky:

******************************************************
pixtosky: Coordinate transformation to sky coordinates
******************************************************
This task allows a user to perform coordinate transformations
with the full WCS and distortion model on source positions
from an input image to sky coordinates.  This task serves as
a replacement for the IRAF.STSDAS task `xy2rd`, albeit with the
added capability of understanding the full distortion model
provided in the headers of images updated for use with astrodrizzle and
tweakreg.

.. moduleauthor:: Warren Hack <help@stsci.edu>

.. automodule:: drizzlepac.pixtosky

.. autofunction:: xy2rd
