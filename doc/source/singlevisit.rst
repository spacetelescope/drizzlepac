.. _singlevisit:

==============================
Single-visit Mosaic Processing
==============================

Observations taken within a single visit represent data intended to produce a
multi-wavelength view of the objects in the desired field-of-view. Most of the
time, those observations are taken using pairs of guide stars to provide very
stable field-of-view throughout the entire visit resulting in images which overlap
almost perfectly.  Unfortunately, this is not
always possible due to increasing limitations of the aging telescope systems.
The result is that an increasing number of visits are taken where observations 
drift and/or roll during the course of the visit or re-acquire at slightly 
different pointings from one orbit to the next.  Whatever the reasons, data across
each visit can not automatically be assumed to align.  Single-visit mosaic (SVM) 
processing attempts to correct these relative
alignment errors and to align the data to an absolute astrometric frame so that
all the data across all the filters used can be drizzled onto the same pixel grid.

 
