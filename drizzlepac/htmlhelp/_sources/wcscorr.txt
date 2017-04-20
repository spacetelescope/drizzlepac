***************************************
Functions to Manage WCS Table Extension
***************************************
These functions provide the basic support for initializing,
creating and updating the WCS table extension which serves as
the archive of updates made to the WCS information in the image headers.

.. moduleauthor:: Warren Hack <help@stsci.edu>

.. automodule:: stwcs.wcsutil.wcscorr
   :members:
   :undoc-members:


*************************************************************
Functions to Manage Legacy OPUS WCS Keywords in the WCS Table
*************************************************************
The previously released versions of ``makewcs`` provided
with `MultiDrizzle` *archives* the original OPUS generated WCS
keywords using header keywords which have a prefix of "O",
such as "OCRPIX1".  In order to avoid overwriting or ignoring
these original values, these functions can be used to convert
the prefixed OPUS WCS keywords into WCS table entries compatible
with the new code.

Strictly to provide complete support for these OPUS keywords,
the code will also create, if the user desires, prefix "O" WCS
keywords from the alternate WCS FITS conventions OPUS keywords.
This would allow images processed using the new code only can
then be used with older versions of `MultiDrizzle`, if the user
needs such compatibility.

.. moduleauthor:: Warren Hack <help@stsci.edu>

.. automodule:: stwcs.wcsutil.convertwcs
   :members:
   :undoc-members:
