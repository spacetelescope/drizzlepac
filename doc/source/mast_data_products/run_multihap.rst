.. _runmultihap_api:

The MVM processing can be performed using one of two interfaces:
  - command-line interface: :ref:`runmultihap<runmultihap_api>`
  - Python interface: :doc:`hapmultisequencer<mvm_api>`

==========================================
Command-line Interface for MVM Processing
==========================================
.. toctree::
   :maxdepth: 2

The task ``runmultihap`` serves as the primary interface for processing data
from multiple visits that share a single skycell into a uniform set of images.

.. automodule:: drizzlepac.runmultihap
.. autofunction:: drizzlepac.runmultihap.perform

Associated Helper Code
======================
These modules and functions assist with some of the logistics associated with multi-visit processing.

.. _which_skycell:

drizzlepac.haputils.which_skycell
---------------------------------
.. automodule:: drizzlepac.haputils.which_skycell
.. autofunction:: drizzlepac.haputils.which_skycell.report_skycells