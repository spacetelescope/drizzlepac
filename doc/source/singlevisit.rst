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
All files processed as part of a single visit get renamed from the standard
pipeline filenames into something which describes the data more clearly.  The 
convention used for the names of these input and output files uses these 
components:

  * ``<propid>`` : the proposal ID for this visit
  * ``<obsetid>`` : the 2-digit visit ID from this proposal 
  * ``<instr>`` : 3 letter designation of the instrument used for the observations
  * ``<detector>`` : name of the detector used for the observations
  * ``<filter>`` : hyphen-separated list of filter names used for the observations
  * ``<ipppssoo>`` : standard 8 character rootname for a single exposure defined by the pipeline
  * ``<ipppss>`` : standard 6 character designation of the ``<instr>``/``<propid>``/``<obsetid>`` for this visit
  * ``dr[cz].fits`` : suffix for drizzled products
  * ``fl[ct].fits`` : suffix for pipeline-calibrated files
  * ``asn.fits`` : suffix for the association table
  * ``hlet.fits`` : suffix for headerlet files containing the WCS solution used to create the final drizzle products/mosaics
  
These components get combined to create filenames specific to each type of file being
processed.  The following table provides a complete list of all the products 
created as a result of single-visit processing.

.. list-table:: Single-visit product filenames
  :widths: 8 25 83
  :header-rows: 1
  
  * - Product
    - File Type
    - Filename for Files Produced
  * - Exposure
    - drizzle product
    - hst_<propid>_<obsetid>_<instr>_<detector>_<filter>_<ipppssoo>_dr[cz].fits
  * -
    - flat-field product
    - hst_<propid>_<obsetid>_<instr>_<detector>_<filter>_<ipppssoo>_fl[ct].fits
  * - 
    - headerlet file
    - hst_<propid>_<obsetid>_<instr>_<detector>_<filter>_<ipppssoo>_hlet.fits
  * -
    - trailer file
    - hst_<propid>_<obsetid>_<instr>_<detector>_<filter>_<ipppssoo>_trl.txt
  * -
    - preview (full size)
    - hst_<propid>_<obsetid>_<instr>_<detector>_<filter>_<ipppssoo>_dr[cz].jpg
  * - 
    - preview (thumbnail)
    - hst_<propid>_<obsetid>_<instr>_<detector>_<filter>_<ipppssoo>_dr[cz]_thumb.jpg
  * - Filter
    - drizzle product
    - hst_<propid>_<obsetid>_<instr>_<detector>_<filter>_<ipppss>_dr[cz].fits
  * -
    - point-source catalog
    - hst_<propid>_<obsetid>_<instr>_<detector>_<filter>_<ipppss>_point-cat.ecsv
  * -
    - segment-source catalog
    - hst_<propid>_<obsetid>_<instr>_<detector>_<filter>_<ipppss>_segment-cat.ecsv
  * - 
    - trailer file
    - hst_<propid>_<obsetid>_<instr>_<detector>_<filter>_<ipppss>_trl.txt
  * -
    - preview (full size)
    - hst_<propid>_<obsetid>_<instr>_<detector>_<filter>_<ipppss>_dr[cz].jpg
  * -
    - preview (thumbnail)
    - hst_<propid>_<obsetid>_<instr>_<detector>_<filter>_<ipppss>_dr[cz]_thumb.jpg
  * - Total
    - drizzle product
    - hst_<propid>_<obsetid>_<instr>_<detector>_total_<ipppss>_dr[cz].fits
  * -
    - point-source catalog
    - hst_<propid>_<obsetid>_<instr>_<detector>_total_<ipppss>_point-cat.ecsv
  * -
    - segment-source catalog
    - hst_<propid>_<obsetid>_<instr>_<detector>_total_<ipppss>_segment-cat.ecsv
  * - 
    - trailer file
    - hst_<propid>_<obsetid>_<instr>_<detector>_total_<ipppss>_trl.txt
  * -
    - preview (full size)
    - hst_<propid>_<obsetid>_<instr>_<detector>_total_<ipppss>_dr[cz].jpg
  * - 
    - preview (thumbnail)
    - hst_<propid>_<obsetid>_<instr>_<detector>_total_<ipppss>_dr[cz]_thumb.jpg
  * - 
    - color preview (full size)
    - hst_<propid>_<obsetid>_<instr>_<detector>_total_<ipppss>_<filters>_dr[cz].jpg
  * - 
    - color preview (thumbnail)
    - hst_<propid>_<obsetid>_<instr>_<detector>_total_<ipppss>_<filters>_dr[cz]_thumb.jpg
    

Input Data
===========



Output Products
===============

Unique SVM Keywords
--------------------
