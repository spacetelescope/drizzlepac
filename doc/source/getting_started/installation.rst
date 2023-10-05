Installation
------------

Conda (Recommended)
===================

``Drizzlepac`` is installed when you install the ``stenv`` conda environment, formerly (astropy). Select your desired release through the `stenv readthedocs page <https://stenv.readthedocs.io/en/latest/getting_started.html>`_ and follow the instructions for installation. 


From Source
===========

Clone this repository
*********************
.. code-block:: shell

    git clone https://github.com/spacetelescope/drizzlepac
    cd drizzlepac


Build the documentation
=======================

*Note:* If you intend to use ``drizzlepac``'s embedded help feature from within
an interactive ``python`` or ``ipython`` session, we recommend you do not skip
this step.


.. code-block:: shell

    cd doc/
    make html


Install drizzlepac
==================

.. code-block:: shell

    pip install .

The option ``--no-use-pep517`` MAY be required in order to correctly build 
the C extensions with ``pip`` versions up to 22.2, after commenting out 
the ``build-backend`` from the ``pyproject.toml`` config file.

**Support for installing using `pip` is still evolving, so use of this 
command is provided on an experimental basis for now.**
