.. _multivisit_api:

******************************************
Python Interface for MVM Processing Code
******************************************
The `~drizzlepac.hapmultisequencer` module is part of the drizzlepac package and calls
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


.. _make_poller_files:

haputils.make_poller_files
---------------------------
Module for creating CSV-formatted input 'poller' files for MVM processing.

.. automodule:: drizzlepac.haputils.make_poller_files
.. autofunction:: drizzlepac.haputils.make_poller_files.generate_poller_file


haputils.poller_utils
---------------------
Module for interpreting input 'poller' files for MVM processing.

.. autofunction:: drizzlepac.haputils.poller_utils.interpret_mvm_input


haputils.product
-----------------
Module for defining classes for MVM processing.

.. autofunction:: drizzlepac.haputils.product.SkyCellExposure
.. autofunction:: drizzlepac.haputils.product.SkyCellProduct


.. _generate_custom_svm_mvm_param_file:

haputils.generate_custom_svm_mvm_param_file
---------------------------------------------
Module for defining custom sets of parameters for SVM and MVM processing.

.. automodule:: drizzlepac.haputils.generate_custom_svm_mvm_param_file
.. autofunction:: drizzlepac.haputils.generate_custom_svm_mvm_param_file.make_svm_input_file


