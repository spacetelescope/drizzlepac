.. _multivisit:

=============================
Multi-visit Mosaic Processing
=============================

Multi-Visit Mosaic (MVM)
    A Multi-Visit Mosaic (MVM) is a single image (product) made by combining all observations taken of the same part of the sky.

Observations taken of the same part of the sky over all the years that HST has been operational can enable unique science
due to the high-resolution of the HST cameras.  Generating useful mosaics from these observations, though, requires
solving a number of key problems; namely,

  * aligning all the images to the same coordinate system
  * define the size on the sky of each mosaic
  * defining what exposures should go into each mosaic

MVM processing implemented as part of the Hubble Advanced Products (HAP) pipeline generates new products based on
solutions implemented for these critical issues.


Alignment
==========
Generating MVM products relies entirely on the definition of the WCS of each input image in order to define where each
exposure will land in the output mosaic.  Errors in an input exposure's WCS will lead to misalignment of overlapping
images resulting in smeared sources or even multiple images of the overlapping sources in the output mosaic.  Errors
in the WCS primarily stem from images being aligned to different astrometric catalogs.  For example,
one exposure may be aligned successfully to the GAIADR2 catalog, yet another overlapping exposure can only be aligned
to the 2MASS catalog due to the lack of GAIA sources in the offset exposure.

MVM processing relies entirely on the alignment that can be performed on the exposures during standard pipeline processing
and subsequently as part of the Single-Visit Mosaic (SVM) processing.  This allows each visit to be aligned to the most
accurate catalog available while taking into account the proper motions of the astrometric sources in the field based on
the date the exposures were taken in the visit.

By default, the most current GAIA catalog will
be used as the default astrometric catalog for aligning all exposures which typically results in images that are aligned
well enough to avoid serious misalignment artifacts in the regions of overlap. Unfortunately, due to the nature of some
of the fields observed by HST, not all exposures can be aligned to a GAIA-based
catalog using GAIA astrometric sources.  As a result, the sources used for alignment will have much larger errors on the
sky compared to the errors of GAIA astrometric sources which can lead to offsets from the GAIA coordinate system of a
significant fraction of an HST pixel (or even multiple HST pixels).  This will lead to artifacts in MVM products where
such exposures overlap other exposures aligned to other catalogs, especially GAIA-based catalogs.


Defining the Mosaic on the Sky
==============================
MVM products need to be defined in order to understand what input exposures will contribute to each mosaic.  The solution
implemented for HAP MVM products relies on tesselation of the entire sky using the same basic tiles defined by the
PanSTARRS project as described at `PanSTARRS Sky tessellation patterns
<https://outerspace.stsci.edu/display/PANSTARRS/PS1+Sky+tessellation+patterns>`_.

.. figure::
  :width: 537 pixels
  :target: images/figure_aitoff.png
  :alt: Aitoff plot of all 2,009 PS1 projection cells for the 3PI survey.

  Aitoff plot of all 2,009 PS1 projection cells for the 3PI survey.  The coverage extends from declination −30° to the
  north celestial pole.


.. image::
  :width: 433 pixels
  :target: images/figure_pole.png
  :alt: PS1 projection cells near the north celestial pole.

  PS1 projection cells near the north celestial pole, where the image overlap is greatest due to convergence of the RA grid.
  The projection cells are 4°x4° in size and are on rings spaced by 4° in declination.

MVM processing uses these pre-defined projection cells after extending them to cover the entire sky
with cells that are actually 4.2°x4.2° in size and are on rings spaced by 4° apart in declination.  This provides
sufficient overlap to allow data from one projection cell to be merged with data from a neighboring cell if necessary.

All these definitions, including finely each ProjectionCell should be divided to define the SkyCells, are encoded in a
FITS table installed as part of the `drizzlepac` package:

    drizzlepac/pars/allsky_cells.fits


Defining each SkyCell
----------------------
Each projection cell is split into a grid of 21 x 21 'sky cells' which serves as the most basic MVM product generated
during MVM processing.  Sky cells are approximately 0.2degrees x 0.2degrees in size (~21500 x ~21500 pixels) and
they have the same WCS as the ‘projection cell’.  Each skycell gets identified by its position within the projection cell
as shown in this figure:

.. image::
  :width: 620 pixels
  :target: images/SkyCell_numbering.png
  :alt: Numbering convention for SkyCells within a Projection Cell.

  Numbering convention for SkyCells within a Projection Cell used for naming the SkyCell.

This provides a way to uniquely identify any position on the sky that can be used as the basis for a unique filename for
all products generated from all the expsosures that overlap each SkyCell.


The WCS for each SkyCell gets defined as a subarray of the Projection cell's WCS.  This allows data across SkyCells in
the same projection cell to be combined into larger mosaics as part of the same tangent plane without performing any
additional resampling.


Defining SkyCell Image Exposures
---------------------------------
Defining the SkyCell for a region on the sky allows for the identification of all exposures that overlap that WCS.
However, creating a single mosaic from data taken with different detectors and filters would not result in a
meaningful result.  Therefore, the exposures that overlap each SkyCell get grouped based on the detector and filter used
to take the exposure to define a 'layer' of the SkyCell.  Each layer can then be generated as the primary basic image
product for each SkyCell.  Exposures taken with spectroscopic elements, like grisms and prisms, and exposures taken of
moving targets can not be used to create layers due to the inability to align them with the rest of the observations.
Therefore, only images taken with standard filters (like the WFC3/UVIS F275W filter) will be used to define SkyCell
mosaics (layers).

The default plate scale for all MVM image products for each SkyCell has been defined as 0.04"/pixel to match the higher
resolution imaging performed by the WFC3/UVIS detector.  However, WFC3/IR data suffers from serious resampling artifacts
when drizzling IR data to that plate scale. So in addition to creating IR mosaics at the 0.04"/pixel 'fine' plate scale,
IR mosaics are also generated at a 'coarse' plate scale of 0.12"/pixel to minimize the resampling artifacts while also being easily
scaled to the 'fine' plate scale mosaics.

SkyCell Example
'''''''''''''''
For example, observations have been taken of NGC5474 with both the ACS and WFC3 cameras.  The ACS observations were taken
with the ACS/WFC detector using the F814W and F606W filters, while the WFC3 observations were taken using the IR detector
using the F110W and F160W filters.  All these observations fall within the SkyCell at position X06y11 in the
projection cell p2381, but given the dramatic plate scale differences, these observations can not be used to create a
single mosaic.  Instead, 6 separate layers get defined for this SkyCell; namely,

  * acs_wfc_f814w  (0.04"/pixel)
  * acs_wfc_f606w  (0.04"/pixel)
  * wfc3_ir_f110w  (0.04"/pixel)
  * wfc3_ir_f160w  (0.04"/pixel)
  * wfc3_ir_f110w_coarse  (0.12"/pixel)
  * wfc3_ir_f160w_coarse  (0.12"/pixel)

Since they all have the same WCS, modulo the plate scale differences, they can be overlaid directly with each other for
analysis.

