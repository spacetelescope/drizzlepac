.. _runmultihap_api:

===================
API for runmultihap
===================
.. toctree::
   :maxdepth: 2

The task ``runmultihap`` serves as the primary interface for processing data
from multiple visits that share a single skycell into a uniform set of images.

.. automodule:: drizzlepac.runmultihap
.. autofunction:: drizzlepac.runmultihap.perform

Associated Helper Code
======================
These modules and functions assist with some of the logistics associated with multi-visit processing.

.. _make_poller_files:

drizzlepac.haputils.make_poller_files
-------------------------------------
.. automodule:: drizzlepac.haputils.make_poller_files
.. autofunction:: drizzlepac.haputils.make_poller_files.generate_poller_file

.. _which_skycell:

drizzlepac.haputils.which_skycell
---------------------------------
.. automodule:: drizzlepac.haputils.which_skycell
.. autofunction:: drizzlepac.haputils.which_skycell.report_skycells