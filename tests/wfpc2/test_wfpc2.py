import os
import pytest

from stsci.tools import teal
from drizzlepac import astrodrizzle
from stwcs import updatewcs
from ..resources import BaseWFPC2


class TestWFPC2(BaseWFPC2):

    @pytest.mark.skip(reason="disable until truth files can be updated")
    def test_waiver_single(self):
        """ This test confirms that drizzlepac can correcly process .
        """
        # Prepare input files.
        raw_inputs = ["u40x010hm_c0f.fits", "u40x010hm_c1f.fits"]
        inputs = [os.path.basename(self.get_input_file('input', i))
                  for i in raw_inputs]

        output = 'wfpc2_single_waiver'
        outfile = '{}_drz.fits'.format(output)
        reffile = 'reference_single_waiver.fits'

        # Update WCS for all inputs
        driz_inputs = updatewcs.updatewcs(inputs[0], use_db=False)

        # run astrodrizzle now...
        adriz_parobj = teal.load('astrodrizzle', defaults=True)
        adriz_parobj['output'] = output
        adriz_parobj['build'] = True
        adriz_parobj['in_memory'] = False
        adriz_parobj['runfile'] = 'wfpc2_single_waiver.log'
        adriz_parobj['STATE OF INPUT FILES']['preserve'] = False
        adriz_parobj['STATE OF INPUT FILES']['clean'] = True
        adriz_parobj['STEP 1: STATIC MASK']['static'] = False
        adriz_parobj['STEP 2: SKY SUBTRACTION']['skysub'] = True
        adriz_parobj['STEP 2: SKY SUBTRACTION']['skystat'] = 'mode'
        adriz_parobj['STEP 2: SKY SUBTRACTION']['skywidth'] = 0.3
        adriz_parobj['STEP 2: SKY SUBTRACTION']['skylower'] = -100.0
        adriz_parobj['STEP 2: SKY SUBTRACTION']['use_static'] = False
        adriz_parobj['STEP 2: SKY SUBTRACTION']['sky_bits'] = None
        adriz_parobj['STEP 3: DRIZZLE SEPARATE IMAGES']['driz_separate'] = False
        adriz_parobj['STEP 4: CREATE MEDIAN IMAGE']['median'] = False
        adriz_parobj['STEP 5: BLOT BACK THE MEDIAN IMAGE']['blot'] = False
        adriz_parobj['STEP 6: REMOVE COSMIC RAYS WITH DERIV, DRIZ_CR']['driz_cr'] = False

        astrodrizzle.AstroDrizzle(driz_inputs, configobj=adriz_parobj)

        # Compare results
        outputs = [(outfile, reffile)]
        self.compare_outputs(outputs)

    @pytest.mark.skip(reason="disable until truth files can be updated")
    def test_waiver_asn(self):
        """ This test confirms that drizzlepac can correcly process input
        WFPC2 data stored in WAIVER fits format.
        """

        # Prepare input files.
        raw_inputs = ['u40x010hm_c0f.fits', 'u40x010im_c0f.fits',
                      'u40x010jm_c0f.fits', 'u40x010km_c0f.fits',
                      'u40x010hm_c1f.fits', 'u40x010im_c1f.fits',
                      'u40x010jm_c1f.fits', 'u40x010km_c1f.fits']
        inputs = [os.path.basename(self.get_input_file('input', i))
                  for i in raw_inputs]

        output = 'wfpc2_waiver'
        outfile = '{}_drz.fits'.format(output)
        reffile = 'reference_wfpc2_asn_waiver.fits'

        # Update WCS for all inputs
        driz_inputs = updatewcs.updatewcs(inputs[:4], use_db=False)

        # run astrodrizzle now...
        adriz_parobj = teal.load('astrodrizzle', defaults=True)
        adriz_parobj['output'] = output
        adriz_parobj['build'] = True
        adriz_parobj['in_memory'] = True
        adriz_parobj['runfile'] = 'wfpc2_asn_waiver.log'
        adriz_parobj['STATE OF INPUT FILES']['preserve'] = False
        adriz_parobj['STATE OF INPUT FILES']['clean'] = True
        adriz_parobj['STEP 1: STATIC MASK']['static_sig'] = 3.0
        adriz_parobj['STEP 2: SKY SUBTRACTION']['skysub'] = False

        astrodrizzle.AstroDrizzle(driz_inputs, configobj=adriz_parobj)

        # Compare results
        outputs = [(outfile, reffile)]
        self.compare_outputs(outputs)

    def test_wfpc2_single(self):
        """ This test confirms that drizzlepac can correcly process single
        WFPC2 exposures.
        """

        # Prepare input files.
        raw_inputs = ["u9yq0703m_c0m.fits", "u9yq0703m_c1m.fits"]
        inputs = [os.path.basename(self.get_input_file('input', i))
                  for i in raw_inputs]

        output = 'wfpc2_single_mef'
        outfile = '{}_drz.fits'.format(output)
        reffile = 'reference_single_mef.fits'

        # Update WCS for all inputs
        driz_inputs = updatewcs.updatewcs(inputs[0], use_db=False)

        # run astrodrizzle now...
        adriz_parobj = teal.load('astrodrizzle', defaults=True)
        adriz_parobj['output'] = output
        adriz_parobj['build'] = True
        adriz_parobj['in_memory'] = False
        adriz_parobj['runfile'] = 'wfpc2_single_mef.log'
        adriz_parobj['STATE OF INPUT FILES']['preserve'] = False
        adriz_parobj['STATE OF INPUT FILES']['clean'] = True
        adriz_parobj['STEP 1: STATIC MASK']['static'] = False
        adriz_parobj['STEP 2: SKY SUBTRACTION']['skysub'] = False
        adriz_parobj['STEP 3: DRIZZLE SEPARATE IMAGES']['driz_separate'] = False
        adriz_parobj['STEP 4: CREATE MEDIAN IMAGE']['median'] = False
        adriz_parobj['STEP 5: BLOT BACK THE MEDIAN IMAGE']['blot'] = False
        adriz_parobj['STEP 6: REMOVE COSMIC RAYS WITH DERIV, DRIZ_CR']['driz_cr'] = False

        astrodrizzle.AstroDrizzle(driz_inputs, configobj=adriz_parobj)

        # Compare results
        outputs = [(outfile, reffile)]
        self.compare_outputs(outputs)

    def test_mef_asn(self):
        """ This test confirms that drizzlepac can correcly process input
        WFPC2 data stored in Multi-extensions FITS(MEF) format.
        """

        # Prepare input files.
        raw_inputs = ['u9yq0703m_c0m.fits', 'u9yq0704m_c0m.fits',
                      'u9yq0707m_c0m.fits', 'u9yq0708m_c0m.fits',
                      'u9yq0703m_c1m.fits', 'u9yq0704m_c1m.fits',
                      'u9yq0707m_c1m.fits', 'u9yq0708m_c1m.fits']
        inputs = [os.path.basename(self.get_input_file('input', i))
                  for i in raw_inputs]

        output = 'wfpc2_mef'
        outfile = '{}_drz.fits'.format(output)
        reffile = 'reference_wfpc2_asn_mef.fits'

        # Update WCS for all inputs
        driz_inputs = updatewcs.updatewcs(inputs[:4], use_db=False)

        # run astrodrizzle now...
        adriz_parobj = teal.load('astrodrizzle', defaults=True)
        adriz_parobj['output'] = output
        adriz_parobj['build'] = True
        adriz_parobj['in_memory'] = True
        adriz_parobj['runfile'] = 'wfpc2_asn_mef.log'
        adriz_parobj['STATE OF INPUT FILES']['preserve'] = False
        adriz_parobj['STATE OF INPUT FILES']['clean'] = True
        adriz_parobj['STEP 2: SKY SUBTRACTION']['skysub'] = False

        astrodrizzle.AstroDrizzle(driz_inputs, configobj=adriz_parobj)

        # Compare results
        outputs = [(outfile, reffile)]
        self.compare_outputs(outputs)
