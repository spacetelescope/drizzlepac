DrizzlePac Interface
====================
.. _drizzlepac_api:

DrizzlePac is written in a hierarchical manner with top-level code in python, and some lower-level code in C. The primary user interface for drizzling individual frames is the astrodrizzle.py ``Astrodrizzle()`` function, which accepts input parameters as variables or as a user-specified configuration (.cfg) file. Examples of how to use these functions can be found in the DrizzlePac `Notebook Tutorials <https://github.com/spacetelescope/notebooks/tree/master/notebooks/DrizzlePac>`_. 

The main steps of the drizzling process are listed below. For more information on these steps please refer to the DrizzlePac `Handbook <https://www.stsci.edu/scientific-community/software/drizzlepac.html>`_. Additional information on the Hubble Advanced Product (HAP) Single- and Multi-visit mosaics can be found in the MAST Data Products section.

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