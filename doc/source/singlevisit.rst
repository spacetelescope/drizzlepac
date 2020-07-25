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


Primary User-Interface
=======================
One task has been written to perform the single-visit processing: ``runsinglehap``. 
It gets used by STScI to generate the single-visit products which
can be found in the Mikulsik Archive for Space Telescopes (MAST) archive. This task 
can also be run from the operating system command-line or from within a
Python session to reproduce those results, or with modification of the input 
parameters, perhaps improve on the standard archived results.  Full details on 
how to run this task can be found in the description of the task at :ref:`runsinglehap_api`.


Processing Steps
================
Single-visit processing performed by ``runsinglehap`` 
relies on the results of the standard astrometric 
processing of the individual exposures and associations as the starting point
for alignment. This processing then follows these steps to create the final products:

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
 

.. _svm_naming_convention:

Single Visit Naming Convention
==============================
All files processed as part of a single visit get renamed from the standard
pipeline filenames into something which describes the data more clearly.  The 
convention used for the names of these input and output files uses these 
components:

  * **<propid>** : the proposal ID for this visit
  * **<obsetid>** : the 2-digit visit ID from this proposal 
  * **<instr>** : 3 letter designation of the instrument used for the observations
  * **<detector>** : name of the detector used for the observations
  * **<filter>** : hyphen-separated list of filter names used for the observations
  * **<ipppssoo>** : standard 8 character rootname for a single exposure defined by the pipeline
  * **<ipppss>** : standard 6 character designation of the **<instr>**/**<propid>**/**<obsetid>** for this visit
  * **dr[cz].fits** : suffix for drizzled products
  * **fl[ct].fits** : suffix for pipeline-calibrated files
  * **asn.fits** : suffix for the association table
  * **hlet.fits** : suffix for headerlet files containing the WCS solution used to create the final drizzle products/mosaics
  
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
    

Processing the Input Data
=========================
SVM processing starts with a list of all the single exposures 
which were taken as part of a visit.  Any associations which were defined by the
proposal are ignored, since the visit itself gets treated, in essence, as a new 
association.  The input files can be specified either using the **poller** file format
used by the STScI automated processing or a file with a simple list of filenames.

Automated poller input file format
----------------------------------
The automated processing performed to populate the MAST archive at 
STScI provides a file with the following format::

    ic0s17h4q_flt.fits,12861,C0S,17,602.937317,F160W,IR,ic0s/ic0s17h4q/ic0s17h4q_flt.fits
    ic0s17h5q_flt.fits,12861,C0S,17,602.937317,F160W,IR,ic0s/ic0s17h5q/ic0s17h5q_flt.fits
    ic0s17h7q_flt.fits,12861,C0S,17,602.937317,F160W,IR,ic0s/ic0s17h7q/ic0s17h7q_flt.fits
    ic0s17hhq_flt.fits,12861,C0S,17,602.937317,F160W,IR,ic0s/ic0s17hhq/ic0s17hhq_flt.fits

This example comes from the 'ic0s1' visit where the columns are:

  #. exposure filename
  #. proposal ID (numeric value)
  #. program ID - ppp value from exposure filename
  #. obset_id - visit number from proposal 
  #. exposure time of the exposure
  #. filters used for the exposure
  #. detector used to take the exposure
  #. location of the exposure in a local cache


Filtering the input data
--------------------------
Not all HST imaging observations can be aligned using SVM processing.  Observations
taken with the GRISM or in SPATIAL SCAN mode result in sources which can not be 
aligned, for example.  The :ref:`analyze_api` module evaluates all
input exposures using these header keywords for the stated rejection criteria.

.. list-table:: Single-visit product filenames
  :widths: 26 27 60
  :header-rows: 1
  
  * - Header Keyword
    - Values Which Trigger Rejection
    - Explanation
  * - OBSTYPE
    - (not IMAGING)
    - Only Imaging mode data processed
  * - MTFLAG
    - T 
    - No moving targets, WCS and background sources vary
  * - SCAN_TYP
    - C or D (or not N)
    - Can not align streaked sources
  * - FILTER or FILTER1, FILTER2
    - C*, PR*
    - G=Grism and PR=Prism, Can not align streaked sources
  * - EXPTIME
    - 0 
    - no exposure time, no data to align
  * - TARGNAME
    - DARK, TUNGSTEN, BIAS, FLAT, 
    - No alignable external sources in these calibration modes 
  * - 
    - EARTH-CALIB, DEUTERIUM
    - No alignable external sources in these calibration modes 
  * - CHINJECT
    - not NONE
    - No alignable external sources in these calibration modes 


Any observation which meets any of these criteria are flagged to be ignored (not
processed).  All observations which are alignable based on these criteria are then
passed along as a table to create the SVM products.  Those inputs which can be
processed are then copied and renamed using the :ref:`svm_naming_convention`.  This 
insures that no SVM processing will affect or otherwise modify the original 
pipeline-processed input files.  Only the SVM named input files will be updated
with new SVM-aligned WCS solutions and then used to produce the drizzle products.  


Defining the Output Products
=============================
The table with the set of observations which can be processed now gets interpreted 
in order to identify what exposures can be combined to create unique products.  
This interpretation gets performed using the code in :ref:`poller_utils_api` by
grouping similar observations.    The rules used for grouping the inputs into output
products result in outputs which have the same detector and filter.  These products
are referred to as **filter** products.  

All exposures for a single detector are also identified and grouped to 
define a **total** product.  This **total** product provides the deepest available 
view of the field-of-view from this visit which will be used to produce the master
catalog of sources for this visit.  The master catalog of source positions will
be used to perform photometry on each exposure, whether the source can be identified
in the exposure at that position or not.  This **forced photometry** results in
limits for the photometry in cases where the sources are not bright enough to be
identified in a given filter.

Two separate source catalogs for each filter are also pre-defined; namely, 

  * a point-source catalog derived using ``photutils`` ``DAOStarFinder`` 
  * a segmentation-based catalog derived using ``photutils`` segmentation code

These two catalogs provide complimentary views of each field-of-view to try to
highlight all types of compact sources found in the exposures. 


Example Visit
--------------
For example, a relatively simple visit of a fairly bright and crowded field with 
6 F555W exposures (two 15-second and four 30-second exposures) and 
6 F814W exposures (two 5-second and four 15-second exposures)
would result in the definition of these output products: 

  * a drizzled image for each separate exposure
  * a single F555W product
  * a single F814W product, and
  * a single **total** product
  * a point-source catalog for the F555W product
  * a segmentation-based source catalog for the F555W product
  * a point-source catalog for the F814W product
  * a segmentation-based source catalog for the F814W product
  
The function ``haputils.poller_utils.interpret_obset_input`` serves as the sole interface 
for this interpretation. A basic tree gets defined (as a dictionary of dictionaries) 
by this function where the
output exposures are identified along with all the names of the input exposures.
This tree then serves as the basis for organizing the rest of the SVM processing.

In addition to defining what output products need to be generated, all the SVM
products names are defined using the :ref:`svm_naming_convention`.  This insures
that all the output products have filenames which are not only unique but also 
understandable (if a bit long) that are easily grouped on disk.  


Aligning the Input Data
=======================
All input exposures should have already been aligned either individually or by 
association table as close to GAIA as possible during standard pipeline calibration
processing.  However, each exposure or association (of exposures) can be aligned
to slightly different fits or catalogs due to differences in the source objects 
which can be identified in each separate exposure.  The primary goal of SVM 
processing is to refine this alignment so that all exposures in the visit for 
the same detector (those exposures which contribute to each **total** product)
share the same WCS (pixels on the sky).  

Alignment of all the exposures for a **total** product uses the same alignment
code as the standard calibration pipeline.  The basic steps it follows is:

  * generate a source catalog for each exposure (using :ref:`amutils_api`)
  * obtain the WCS from each exposure
  * perform a relative fit between the exposures using ``tweakwcs``
  * obtain an astrometric catalog for the field-of-view 
  * perform a final fit of all the exposures at once to the astrometric catalog
  * update each WCS with the final corrected WCS generated by ``tweakwcs``
  
The limits for performing the relative alignment and absolute fit to the astrometric
catalog (defaults to **GAIADR2**) are lower under the expectation that large 
offsets (> 0.5 arcseconds) have already been removed in the pipeline.  This makes
the SVM alignment more robust across a wider range of types of fields-of-view.
The final updated WCS will be provided with a name that reflects this cross-filter
alignment using **-FIT-SVM-<catalog name>** as the final half of the **WCSNAME** 
keyword.  More details on the WCS naming conventions can be found in the
:ref:`wcsname-conventions` section.



Creating the Output Products
============================



Unique SVM Keywords
-------------------
A small set of keywords have been added to the standard drizzle headers to reflect
the unique characteristics of the SVM products.  These keywords are:

.. glossary::

  NPIXFRAC
    Fraction of pixels with data
    
  MEANEXPT
    Mean exposure time per pixel with data
  
  MEDEXPT
    Median exposure time per pixel with data
  
  MEANNEXP
    Mean number of exposures per pixel with data
  
  MEDNEXP
    Median number of exposures per pixel with data
    



