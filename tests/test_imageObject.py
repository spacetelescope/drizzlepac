#!/usr/bin/env python

import pytest
from drizzlepac import imageObject
from drizzlepac.haputils import astroquery_utils as aqutils

#from http://blog.moertel.com/articles/2008/03/19/property-checking-with-pythons-nose-testing-framework
def forall_cases(cases):
    def decorate(testfn):
        def gen():
            for case in cases:
                yield testfn, case
        gen.__name__ = "test_%s_for_a_case" % testfn.__name__
        return gen
    return decorate


class TestimageObject():
    """test the implementation of imageObject by creating and examining an object
       from an image file thats on disk
    """

    def noFilename(self, filename):
        image=imageObject.imageObject(filename)
    def test_NoFilename(self,filename=''):
        with pytest.raises(IOError):
            self.noFilename(filename)

    def test_Attributes(self, filename="./j6ll01yiq_flt.fits"):
        input_files = []
        input_files = aqutils.retrieve_observation('j6ll01yiq*', suffix=["FLT"], product_type="pipeline")
        image=imageObject.imageObject(input_files[0])

        #just check to make sure the global attributes are not empty
        assert(image._filename != '')
        assert(image._rootname != '')
        assert(image._instrument != '')
