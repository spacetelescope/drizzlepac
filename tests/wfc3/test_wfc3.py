import os

from stsci.tools import teal
import drizzlepac
from drizzlepac import astrodrizzle

from stwcs import updatewcs

from ..resources import BaseWFC3
from ci_watson.hst_helpers import raw_from_asn


class TestWFC3(BaseWFC3):

    def test_binned_single(self):
        rootname = 'iacs01t9q'
        input_name = '{}_flt.fits'.format(rootname)
        output = '{}_drz.fits'.format(rootname)
        ref_file = 'reference_wfc3_binned_single.fits'
        runfile = ref_file.replace('.fits','.log')

        # Prepare input files.
        input_file = os.path.basename(self.get_input_file('input', input_name))

        input_pars = {'output': output, 'runfile':runfile,
                      'build':True, 'in_memory':True,
                      'preserve': False, 'clean':True, 'static':False,
                      'skysub':False, 'driz_separate':False, 'median':False,
                      'blot':False, 'driz_cr':False}

        # Update wcs in input file
        driz_input = updatewcs.updatewcs(input_file)

        # run astrodrizzle now...
        parObj = teal.load('astrodrizzle', defaults=True)  # get all default values
        astrodrizzle.AstroDrizzle(driz_input, configobj=parObj, **input_pars)

        # Compare results
        outputs = [(output, ref_file)]
        self.compare_outputs(outputs)


    def test_uvis_single(self):
        rootname = 'iacr51ohq'
        input_name = '{}_flt.fits'.format(rootname)
        output = '{}_drz.fits'.format(rootname)
        ref_file = 'reference_wfc3_uvis_single.fits'
        runfile = ref_file.replace('.fits','.log')

        # Prepare input files.
        input_file = os.path.basename(self.get_input_file('input', input_name))

        input_pars = {'output': output, 'runfile':runfile,
                      'build':True, 'in_memory':True,
                      'preserve': False, 'clean':True, 'static':True,
                      'skysub':True, 'skywidth':0.3, 'use_static':False,
                      'sky_bits':None,
                      'driz_separate':False, 'median':False,
                      'blot':False, 'driz_cr':False}

        # Update wcs in input file
        driz_input = updatewcs.updatewcs(input_file)

        # run astrodrizzle now...
        parObj = teal.load('astrodrizzle', defaults=True)  # get all default values
        astrodrizzle.AstroDrizzle(driz_input, mdriztab=False,
                                  configobj=parObj, **input_pars)

        # Compare results
        outputs = [(output, ref_file)]
        self.compare_outputs(outputs)


    def test_uvis_asn(self):
        rootname = 'iacr51010'
        asn_file = '{}_asn.fits'.format(rootname)
        output = '{}_drz.fits'.format(rootname)

        ref_file = 'reference_wfc3_uvis_asn.fits'
        runfile = ref_file.replace('.fits','.log')

        # Prepare input files.
        input_files = []
        input_asn = self.get_data('input', asn_file)
        for raw_file in raw_from_asn(asn_file, suffix='_flt.fits'):
            input_files.append(os.path.basename(self.get_input_file('input', raw_file)))

        # define input parameters to be used
        input_pars = {'output': output, 'runfile':runfile,
                      'build':True, 'in_memory':False,
                      'preserve': False, 'clean':True, 'static':True,
                      'skysub':True, 'skywidth':0.3, 'use_static':False,
                      'sky_bits':None}

        # run astrodrizzle now...
        parObj = teal.load('astrodrizzle', defaults=True)  # get all default values
        astrodrizzle.AstroDrizzle(asn_file, mdriztab=False,
                                  configobj=parObj, **input_pars)

        # Compare results
        outputs = [(output, ref_file)]
        self.compare_outputs(outputs)


    def test_wfc3_ir_saturated(self):
        # Prepare input files.
        raw_inputs = ['ib5w02n0q_flt.fits', 'ib5w02naq_flt.fits',
                      'ib5w02njq_flt.fits']
        inputs = [os.path.basename(self.get_input_file('input', i))
                  for i in raw_inputs]

        outname = 'wfc3_ir_sat'
        output = '{}_drz.fits'.format(outname)

        ref_file = 'reference_{}.fits'.format(outname)
        runfile = ref_file.replace('.fits','.log')

        # define input parameters to be used
        input_pars = {'output': outname, 'runfile':runfile,
                      'build':True, 'in_memory':False,
                      'preserve': False, 'clean':True, 'static':True,
                      'skysub':True, 'skystat':'mode', 'skylower':-100.0,
                      'use_static':False, 'sky_bits':None, 'driz_sep_bits':832,
                      'driz_cr_snr':'5.0 4.0', 'final_bits': 832}

        # run astrodrizzle now...
        parObj = teal.load('astrodrizzle', defaults=True)  # get all default values
        astrodrizzle.AstroDrizzle(inputs, mdriztab=False,
                                  configobj=parObj, **input_pars)

        # Compare results
        outputs = [(output, ref_file)]
        self.compare_outputs(outputs)
