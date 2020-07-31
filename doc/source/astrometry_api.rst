.. _advanced_products_api:

*******************************
Hubble Advanced Products API
*******************************
These modules provide the basic functionality used to process automatically
data using this package to apply the distortion models to the WCS of HST 
observations and to verify the alignment of the observations. 

.. _amutils_api:

haputils.astrometric_utils
--------------------------

.. automodule:: drizzlepac.haputils.astrometric_utils
.. autofunction:: drizzlepac.haputils.astrometric_utils.create_astrometric_catalog
.. autofunction:: drizzlepac.haputils.astrometric_utils.get_catalog
.. autofunction:: drizzlepac.haputils.astrometric_utils.find_gsc_offset
.. autofunction:: drizzlepac.haputils.astrometric_utils.extract_sources
.. autofunction:: drizzlepac.haputils.astrometric_utils.classify_sources
.. autofunction:: drizzlepac.haputils.astrometric_utils.generate_source_catalog
.. autofunction:: drizzlepac.haputils.astrometric_utils.generate_sky_catalog
.. autofunction:: drizzlepac.haputils.astrometric_utils.compute_photometry
.. autofunction:: drizzlepac.haputils.astrometric_utils.filter_catalog
.. autofunction:: drizzlepac.haputils.astrometric_utils.build_self_reference
.. autofunction:: drizzlepac.haputils.astrometric_utils.within_footprint
.. autofunction:: drizzlepac.haputils.astrometric_utils.find_hist2d_offset
.. autofunction:: drizzlepac.haputils.astrometric_utils.build_wcscat
.. autofunction:: drizzlepac.haputils.astrometric_utils.compute_similarity
.. autofunction:: drizzlepac.haputils.astrometric_utils.determine_focus_index
.. autofunction:: drizzlepac.haputils.astrometric_utils.max_overlap_diff
.. autofunction:: drizzlepac.haputils.astrometric_utils.detect_point_sources


.. _analyze_api:

haputils.analyze
-----------------
.. automodule:: drizzlepac.haputils.analyze
.. autofunction:: drizzlepac.haputils.analyze.analyze_data


.. _align_utils_api:

haputils.align_utils
---------------------
.. automodule:: drizzlepac.haputils.align_utils
.. autoclass:: drizzlepac.haputils.align_utils.AlignmentTable
.. autofunction:: drizzlepac.haputils.align_utils.match_relative_fit
.. autofunction:: drizzlepac.haputils.align_utils.match_default_fit
.. autofunction:: drizzlepac.haputils.align_utils.match_2dhist_fit


