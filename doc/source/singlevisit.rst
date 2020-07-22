.. _singlevisit:

==============================
Single-visit Mosaic Processing
==============================
Standard calibration pipeline processing focuses on aligning data taken as 
associations or as single exposures as closely to a standard astrometric coordinate
frame as possible (see :ref:`runastrodriz-description` for full details).  
This involves alignment of very few images all taken under
as identical conditions as possible; namely, detector, filter, guide stars, 
guiding mode, and so on.  These observations were intended to by the user to 
represent a single view of the object of interest in a way that allows for 
removal of as many calibration effects as possible. 

Observations taken within a single visit, on the other hand, represent data 
intended to produce a multi-wavelength view of the objects in the desired
field-of-view. Most of the time, those observations are taken using pairs of guide 
stars to provide very stable field-of-view throughout the entire visit resulting 
in images which overlap almost perfectly.  Unfortunately, this is not
always possible due to increasing limitations of the aging telescope systems.
The result is that an increasing number of visits are taken where observations 
drift and/or roll during the course of the visit or re-acquire at slightly 
different pointings from one orbit to the next.  Whatever the reasons, data across
each visit can not automatically be assumed to align.  Single-visit mosaic (SVM) 
processing attempts to correct these relative
alignment errors and to align the data to an absolute astrometric frame so that
all the data across all the filters used can be drizzled onto the same pixel grid.

Understanding the quality of the SVM products requires knowing what processing
took place to generate that product, and more importantly, what limitations that
processing may have had.  The following sections provide the description of the
single-visit processing.  Definitions of a number of processing-specific terms used 
in this description can be found in the :ref:`hap-glossary`.

Processing Steps
================
Single-visit processing relies on the results of the standard astrometric 
processing of the individual exposures and associations as the starting point
for alignment. This processing then follows these steps to create the final products::

  * interpret the list of filenames for all exposures taken as part of a single visit
  * copy the pipeline-calibrated (FLT/FLC) files to the current directory for processing
  * rename the input files to conform to the single-visit naming conventions
    
    * This step insures that the original pipeline results remain available in 
      the archive unchanged
  
  * Define what output products can be generated 
  * Align all exposures in a relative sense (all to each other)
  * Create a composite source catalog from all aligned input exposures
  * Cross-match and fit this composite catalog to GAIA to determine new WCS solution
  * Update renamed input exposures with results of alignment to GAIA
  * Create each of the output products using the updated WCS solutions
 

Single Visit Naming Convention
==============================



Input Data
===========



Output Products
===============

Unique SVM Keywords
--------------------
