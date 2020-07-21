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
       
     single-visit mosaic (SVM)
       A set of drizzled products for all the observations taken in a single visit.
       All observations in this set of products should be aligned to each other 
       with the exposures for each detector/filter combination being combined into
       it's own product.  All products have the same WCS definition and therefore
       the same pixel grid to enable direct pixel-by-pixel comparison of the 
       data across all the detectors/filters. 
       
     
