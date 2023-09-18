.. _imagereg:

==========================
Image Registration Tasks
==========================

A number of tasks have been developed to support the registration
of images.  This includes documentation for the replacement task
for IRAF's ``tweakshifts``, currently named ``TweakReg``, along
with tasks for updating the WCS in HST images and performing
photometry equalization for WFPC2 data.

These pages describe how to run the new ``TEAL``-enabled task,
as well as use the classes in the task to generate catalogs interactively
for any chip and work with that catalog. The current implementation of this
code relies on a very basic source finding algorithm loosely patterned
after the DAOFIND algorithm and does not provide all the same features
or outputs found in DAOFIND. The fitting algorithm also reproduces the
fitting performed by IRAF's ``geomap`` in a limited fashion; primarily,
it only performs fits equivalent to ``geomap``'s 'shift' and 'rscale'
solutions. These algorithms will be upgraded as soon as replacements
are available.

.. toctree::
   :maxdepth: 2

   ../reference/tweakreg
   ../reference/refimagefindpars
   ../reference/imagefindpars
   ../reference/image
   ../reference/wcscorr
   ../reference/tweakutils
   ../reference/updatehdr
   ../reference/mapreg
   ../reference/photeq
   ../reference/pixreplace


ACS Header Update Task
----------------------
A task, 'updatenpol', has been written to automate the updating of ACS image headers
with the filename of the appropriate NPOLFILE based on the DGEOFILE specified in the image header.
This task should be used to update all ACS images prior to processing them with 'astrodrizzle'.

.. toctree::
   :maxdepth: 2

   ../reference/updatenpol


================================
Coordinate Transformation Tasks
================================

Several tasks have also been developed to perform
coordinate transformations that take into account the full
distortion model specified in HST image headers.

.. toctree::
   :maxdepth: 1

   ../reference/pixtopix
   ../reference/pixtosky
   ../reference/skytopix
