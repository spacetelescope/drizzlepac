import os

from stsci.tools import teal
from drizzlepac import astrodrizzle
from ci_watson.artifactory_helpers import get_bigdata
from ..resources import BaseACS


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
        input_file = self.get_input_file('input', input_file)

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
        self.ignore_keywords += ['D001DATA', 'D001MASK']
        self.compare_outputs(outputs)
