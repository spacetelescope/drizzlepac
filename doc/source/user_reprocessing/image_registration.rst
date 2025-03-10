.. _imagereg:

========================
Image Registration Tasks
========================

A number of tasks have been developed to support the registration
of images.  This includes documentation for the replacement task
for IRAF's ``tweakshifts``, currently named ``TweakReg``, along
with tasks for updating the WCS in HST images and performing
photometry equalization for WFPC2 data.

The current implementation of this
code relies on a very basic source finding algorithm loosely patterned
after the DAOFIND algorithm and does not provide all the same features
or outputs found in DAOFIND. The fitting algorithm also reproduces the
fitting performed by IRAF's ``geomap`` in a limited fashion; primarily,
it only performs fits equivalent to ``geomap``'s 'shift' and 'rscale'
solutions. These algorithms will be upgraded as soon as replacements
are available.

.. toctree::
   :maxdepth: 2

   tweakreg
   tweakutils
   tweakback
   updatehdr
   UPDATEWCS (STWCS) <https://stwcs.readthedocs.io/en/latest/updatewcs.html>
   wcscorr
