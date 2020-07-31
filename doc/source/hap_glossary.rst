.. _hap-glossary:

==================================
Hubble Advanced Products Glossary
==================================

.. glossary::

    *a priori* WCS
      A WCS derived using pre-existing information about how the observation
      was taken.  The primary example of an *a priori* WCS would be the **GSC240**
      WCS solutions from the astrometry database which were derived by correcting
      the coordinates that were used for guiding the observation with updated
      coordinates for those guide stars derived from the GAIA catalog.

    *a posteriori* WCS
      A WCS derived using information measured in the observations themselves and 
      fit to some external reference. A primary example would be the **HSC30**
      WCS solutions from the astrometry database.  These were computed using the
      `Hubble Source Catalog (HSC) <https://archive.stsci.edu/hst/hsc/>`_ 
      positions of sources from the observation, then
      cross-matching and fitting them to sources from the GAIA catalog.  Another 
      example would be the **-FIT** WCS solutions computed during standard pipeline
      calibration processing.  These solutions get derived by identifying sources
      from each exposure, then cross-matching those sources with GAIA catalog 
      sources and applying the fit to the WCS. 
       
    single visit
      This term refers to the set of observations taken by a single instrument
      using the same guide stars potentially spanning multiple orbits.  All 
      observations taken by all detectors of that instrument which were requested
      by the PI in the original proposal are included as this 'single visit' set.    
    
    single-visit mosaic (SVM)
       A set of drizzled products for all the observations taken in a single visit.
       All observations in this set of products should be aligned to each other 
       with the exposures for each detector/filter combination being combined into
       it's own product.  All products have the same WCS definition and therefore
       the same pixel grid to enable direct pixel-by-pixel comparison of the 
       data across all the detectors/filters. 
       
    singleton
       A single exposure taken as part of a visit that is not part of any 
       pre-defined association based on the proposal.
       
    association
       A group of exposures taken together in the same visit, typically defined
       as a single defined observation in the proposal.  For example, a proposal
       may specify use of a DITHER-LINE pattern and request that it be used to
       take 15second observations with the F555W filter.  All those exposures would be 
       'associated' by the proposal and result in defining an association table 
       that specifies the names of each of the 15second exposures and the name 
       of the image created by combining all those exposures.  
       
    association table
       A FITS table listing the filenames of all the input exposures that should
       be combined together, along with the name or names of the products to be
       created by combining the input images. 
       
    FLT/FLC image
       This term refers to the pipeline calibrated version of the input exposures.
       These files have either **_flt.fits(FLT)** or **_flc.fits(FLC)** in their filename, thus
       the term FLT/FLC.  The FLC files are the CTE-corrected versions of the FLT
       files, while both copies have identical WCS solutions.
       
    Active WCS
    Primary WCS
       The WCS solution in the same image header as the science array which 
       defines the transformation from pixel coordinates to world coordinates. 
       The set of keywords that make up this WCS are::
        
         * CRVAL1, CRVAL2
         * CRPIX1, CRPIX2
         * CD1_1, CD1_2, CD2_1, CD2_2
         * CTYPE1, CTYPE2
         * A_\*_\*, B_\*_\*
         
    Alternate WCS
       This term refers to any set of WCS keywords which have a single alphabetic
       character (A-Z) appended to the end, such as CRVAL1A, as defined by
       FITS WCS Paper I.

