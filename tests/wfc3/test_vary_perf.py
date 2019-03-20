import os

from stsci.tools import teal
from drizzlepac import astrodrizzle
from stwcs import updatewcs
from ..resources import BaseWFC3

def pytest_generate_tests(metafunc):
    # called once per each test function
    funcarglist = metafunc.cls.params[metafunc.function.__name__]
    argnames = sorted(funcarglist[0])
    metafunc.parametrize(argnames, [[funcargs[name] for name in argnames]
            for funcargs in funcarglist])


class TestVaryPerf(BaseWFC3):

    params = {
        'test_perf' : [dict(output='wfc3_default_perf', in_memory=False,
                            num_cores=None, use_static=True, skywidth=0.1),
                       dict(output='wfc3_parallel_hi', in_memory=False,
                            num_cores=12, use_static=True, skywidth=0.1),
                       dict(output='wfc3_parallel_off', in_memory=False,
                            num_cores=1, use_static=True, skywidth=0.1),
                       dict(output='wfc3_parallel_on', in_memory=False,
                            num_cores=3, use_static=True, skywidth=0.1),
                       dict(output='wfc3_virtual', in_memory=True,
                            num_cores=1, use_static=True, skywidth=0.1),
                       dict(output='wfc3_virtual_parallel', in_memory=True,
                            num_cores=3, use_static=True, skywidth=0.1),
                       dict(output='wfc3_binned_asn', in_memory=True,
                            num_cores=None, use_static=False, skywidth=0.3),
                      ],
             }

    ignore_keywords = ['origin', 'filename', 'date', 'iraftlm', 'fitsdate', 'upwtim',
           'wcscdate', 'upwcsver', 'pywcsver', 'prod_ver', 'rulefile',
           'history']

    def test_perf(self, output, in_memory, num_cores, use_static, skywidth):

        # Prepare input files.
        raw_inputs = ["ib6m02d9q_flt.fits", "ib6m02daq_flt.fits"]
        inputs = [os.path.basename(self.get_input_file('input', i))
                      for i in raw_inputs]

        # Merge common parameter settings with test-specific settings
        input_pars = {'build':True, 'preserve':False,
                       'clean':True, 'sky_bits':None}
        input_pars['output'] = output
        input_pars['in_memory'] = in_memory
        input_pars['num_cores'] = num_cores
        input_pars['use_static'] = use_static
        input_pars['skywidth'] = skywidth
        run_file = '{}.log'.format(output)
        input_pars['runfile'] = run_file

        # Update WCS for all inputs
        driz_inputs = updatewcs.updatewcs(inputs, use_db=False)

        # run astrodrizzle now...
        parObj = teal.load('astrodrizzle', defaults=True)  # get all default values
        astrodrizzle.AstroDrizzle(driz_inputs, configobj=parObj, **input_pars)

        # Compare results
        outfile = '{}_drz.fits'.format(output)
        reffile = 'reference_{}.fits'.format(output)
        outputs = [(outfile, reffile)]
        self.compare_outputs(outputs)
