#!/usr/bin/env python

import nose
import imageObject

#from http://blog.moertel.com/articles/2008/03/19/property-checking-with-pythons-nose-testing-framework
def forall_cases(cases):
    def decorate(testfn):
        def gen():
            for case in cases:
                yield testfn, case
        gen.__name__ = "test_%s_for_a_case" % testfn.__name__
        return gen
    return decorate



class test_imageObject():
    """test the implementation of imageObject by creating and examining an object
       from an image file thats on disk
    """

    #not sure if this is the correct way to do it
    def testNoFilename(self,filename=''):
        self.assertRaises(IOError,imageObject.imageObect,filename)

    def testAttributes(self, filename="./j8uq10lbq_flt.fits"):
        image=imageObject.imageObject(filename)

        #just check to make sure the global attributes are not empty
        assert(image._filename != '')
        assert(image._naxis1 > 0)
        assert(image._naxis2 > 0)
        assert(image._instrument != '')

if __name__ == "__main__":
    nose.run()
