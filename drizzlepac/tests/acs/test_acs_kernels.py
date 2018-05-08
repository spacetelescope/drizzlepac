import os

from stsci.tools import teal
import drizzlepac
from drizzlepac import astrodrizzle
from drizzlepac.helpers.mark import require_bigdata, remote_data
from drizzlepac.helpers.io import get_bigdata

from ..helpers import BaseACS, download_file_cgi, raw_from_asn

def pytest_generate_tests(metafunc):
    # called once per each test function
    funcarglist = metafunc.cls.params[metafunc.function.__name__]
    argnames = sorted(funcarglist[0])
    metafunc.parametrize(argnames, [[funcargs[name] for name in argnames]
            for funcargs in funcarglist])


class TestAcsKernels(BaseACS):
    
    params = {
        'test_kernels' : [dict(output='acs_square', kernel='square'),
                          dict(output='acs_turbo', kernel='turbo'),
                          dict(output='acs_tophat', kernel='tophat'),
                          dict(output='acs_point', kernel='point'),
                          dict(output='acs_lanczos3', kernel='lanczos3'),
                          dict(output='acs_gaussian', kernel='gaussian'),],
        }
        
    def test_kernels(self, output, kernel):
        
        input_file = 'j8bt06nyq_flt.fits'
        print("Running the test for kernel={}".format(kernel))

        # Prepare input files.
        input_file = get_bigdata(self.input_loc, input_file)
        #download_file_cgi(self.tree, self.input_loc, input_file,
        #                  timeout=self.timeout)
        print("Testing input: {} \nkernel: {}".format(input_file, kernel))
        # run astrodrizzle now...
        parObj = teal.load('astrodrizzle', defaults=True)  # get all default values
        parObj['output'] = output
        parObj['build'] = True
        parObj['in_memory'] = True
        parObj['runfile'] = 'drizzlepac.run'
        parObj['STATE OF INPUT FILES']['preserve'] = False
        parObj['STATE OF INPUT FILES']['clean'] = True
        parObj['STEP 1: STATIC MASK']['static'] = False
        parObj['STEP 2: SKY SUBTRACTION']['skysub'] = False
        parObj['STEP 3: DRIZZLE SEPARATE IMAGES']['driz_separate'] = False
        parObj['STEP 4: CREATE MEDIAN IMAGE']['median'] = False
        parObj['STEP 5: BLOT BACK THE MEDIAN IMAGE']['blot'] = False
        parObj['STEP 6: REMOVE COSMIC RAYS WITH DERIV, DRIZ_CR']['driz_cr'] = False
        parObj['STEP 7: DRIZZLE FINAL COMBINED IMAGE']['final_kernel'] = kernel
        parObj['STEP 7: DRIZZLE FINAL COMBINED IMAGE']['final_units'] = 'counts'
        
        astrodrizzle.AstroDrizzle(input_file, configobj=parObj)

        # Compare results
        outfile = '{}_drz.fits'.format(output)
        reffile = 'reference_{}.fits'.format(kernel)
        outputs = [(outfile, reffile)]
        self.compare_outputs(outputs)
