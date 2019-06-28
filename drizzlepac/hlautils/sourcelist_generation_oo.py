#!/usr/bin/env python

"""This script contains code to support creation of photometric sourcelists using two techniques: aperture photometry
segmentation-map based photometry.
"""
import pdb
__taskname__ = 'sourcelist_generation_oo'


# ----------------------------------------------------------------------------------------------------------------------

class Point_source_photometry(object):
    """Using aperture photometry, generate photometric sourcelist for specified image(s).
    """
    def __init__(self):
        self.label="Point_source_photometry"
        self.description="A set of routines to generate photometric sourcelists using aperture photometry"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def identify_sources(self):
        """Create a master coordinate list of sources identified in the specified total detection product image"""
        print("IDENTIFY_SOURCES")
        pass
        return()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def perform_photometry(self,foo):
        """Perform aperture photometry on identified sources"""
        print("PERFORM_PHOTOMETRY"+foo)
        pass
        return()


# ----------------------------------------------------------------------------------------------------------------------


class Segment_photometry(object):
    """Create the Sextractor-like source catalog using PhotUtils"""

    def __init__(self):
        self.label = "Segment_photometry"
        self.description ="A set of routines to generate photometric sourcelists using segment maps"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def placeholder(self):
        pass

# ----------------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    blarg = Point_source_photometry()
    pdb.set_trace()