.. _processinput:

*************
Process Input
*************

This module supports the interpretation and initial verification of all the
input files specified by the user. These functions:

    * reads in parameter values from MDRIZTAB reference file and merges those merges those values in with the rest of the parameters from the GUI/configObj, if use of MDRIZTAB was specified
    * insure that all input files are multi-extension FITS files and converts them if they are not
    * updates all input WCS's to be consistent with IDCTAB, if specified
    * generates the ``ImageObject`` instances for each input file
    * resets the DQ bits if specified by the user
    * adds info about any user-provided IVM files to the ImageObjects
    * generates the output WCS based on user inputs

.. moduleauthor:: Warren Hack <help@stsci.edu>

.. automodule:: drizzlepac.processInput
   :members:
   :undoc-members:


ResetBits Update of Input
=========================

This module provides the capability to set a specific set of bit values in the input DQ arrays to zero. This allows a user to reset pixels previously erroneously flagged as cosmic-rays to good for reprocessing with improved alignment or detection parameters.

.. automodule:: drizzlepac.resetbits
   :members:
