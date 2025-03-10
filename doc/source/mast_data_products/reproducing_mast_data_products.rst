.. _running-astrodrizzle:

******************************
Reproducing MAST DATA Products
******************************

``runastrodriz`` is a module to control operation of astrodrizzle which removes distortion and combines HST images in the pipeline.


Typical Usage
=============

    >>> runastrodriz.py [-fhibn] inputFilename [newpath]


Alternative Usage
=================

    >>> python
    >>> from wfc3tools import runastrodriz
    >>> runastrodriz.process(inputFilename,force=False,newpath=None,inmemory=False)


GUI Usage under Python (no longer supported)
============================================

    >>> python
    >>> from stsci.tools import teal
    >>> import wfc3tools
    >>> cfg = teal.teal('runastrodriz')


Options
=======

If the '-i' option gets specified, no intermediate products will be written out
to disk. These products, instead, will be kept in memory. This includes all
single drizzle products (*single_sci* and *single_wht*), median image,
blot images, and crmask images.  The use of this option will therefore require
significantly more memory than usual to process the data.

If a value has been provided for the newpath parameter, all processing will be
performed in that directory/ramdisk.  The steps involved are:

    * create a temporary directory under that directory named after the input file
    * copy all files related to the input to that new directory
    * change to that new directory and run astrodrizzle
    * change back to original directory
    * move (not copy) ALL files from temp directory to original directory
    * delete temp sub-directory

The '-b' option will run this task in BASIC mode without creating headerlets
for each input image.

The '-n' option allows the user to specify the number of cores to be used in
running AstroDrizzle.

.. note:: This value will be forced to a value of '1' (one) on Windows systems due to
  exceptions caused by threaded logging under Windows.  Future versions will
  lift this enforced restriction on Windows systems once issues with logging are
  resolved.