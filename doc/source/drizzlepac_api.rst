DrizzlePac Interface
====================
.. _drizzlepac_api:

DrizzlePac is written in a hierarchical manner with top-level code in python, and some lower-level code in C. The primary user interface is the astrodrizzle.py ``Astrodrizzle()`` function. Single- and Multi-Visit Mosaic creation is handled by hapsequencer.py ``run_hap_processing()`` and hapmultisequencer.py ``run_mvm_processing()``, respectively, which perform additional setup before making calls to ``Astrodrizzle()``. Examples of how to use these functions can be found in the DrizzlePac `Notebook Tutorials <https://github.com/spacetelescope/notebooks/tree/master/notebooks/DrizzlePac>`_. 

The main steps of the drizzling process are listed below. For more information on these steps please refer to the DrizzlePac `Handbook <https://www.stsci.edu/scientific-community/software/drizzlepac.html>`_. 

.. toctree:: 
    :maxdepth: 1

    drizzlepac_api/astrodrizzle
    drizzlepac_api/process
    drizzlepac_api/static
    drizzlepac_api/sky
    drizzlepac_api/adrizzle
    drizzlepac_api/median
    drizzlepac_api/ablot
    drizzlepac_api/drizcr
    drizzlepac_api/outimage


Index
-----

* :ref:`genindex`
* :ref:`modindex`