.. _multivisit_api:

******************************************
Python Interface for MVM Processing Code
******************************************
The ``hapmultisequencer.py`` module is part of the drizzlepac package and calls
utilities which are intended to be used with
Multi-visit Mosaic (MVM) data and to be imported into other Python programs or
a Python session for interactive data analysis.

.. toctree::
   :maxdepth: 1

.. _hapmulti :

hapmultisequencer
---------------------
This module provides the primary user interface for performing MVM processing.

.. automodule:: drizzlepac.hapmultisequencer
.. autofunction:: drizzlepac.hapmultisequencer.run_mvm_processing

haputils.cell_utils
-------------------
Module for defining ProjectionCells and SkyCells for a set of exposures.

.. automodule:: drizzlepac.haputils.cell_utils
.. autofunction:: drizzlepac.haputils.cell_utils.get_sky_cells
.. autofunction:: drizzlepac.haputils.cell_utils.interpret_scells
.. autoclass:: drizzlepac.haputils.cell_utils.ProjectionCell
.. autoclass:: drizzlepac.haputils.cell_utils.SkyCell


haputils.make_poller_files
---------------------------
Module for creating CSV-formatted input 'poller' files for MVM processing.

.. automodule:: drizzlepac.haputils.make_poller_files
.. autofunction:: drizzlepac.haputils.make_poller_files.generate_poller_file

haputils.poller_utils
---------------------
Module for interpreting input 'poller' files for MVM processing.

.. automodule:: drizzlepac.haputils.poller_utils
.. autofunction:: drizzlepac.haputils.poller_utils.interpret_mvm_input
.. autofunction:: drizzlepac.haputils.poller_utils.build_poller_table


haputils.product
-----------------
Module for defining classes for MVM processing.

.. automodule:: drizzlepac.haputils.product
.. autofunction:: drizzlepac.haputils.product.HAPProduct
.. autofunction:: drizzlepac.haputils.product.FilterProduct
.. autofunction:: drizzlepac.haputils.product.SkyCellExposure
.. autofunction:: drizzlepac.haputils.product.SkyCellProduct



