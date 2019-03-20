import os

from stsci.tools import teal
import drizzlepac
from drizzlepac import astrodrizzle

from ..resources import BaseACS
from ci_watson.hst_helpers import raw_from_asn


class TestAsnRegress(BaseACS):

    def test_hrc_asn(self):
        rootname = 'j8bt06010'
        asn_file = rootname + '_asn.fits'

        # Prepare input files.
        input_file = self.get_data('input', asn_file)

        for raw_file in raw_from_asn(asn_file, suffix='_flt.fits'):
            self.get_input_file('input', raw_file)

        # run astrodrizzle now...
        parObj = teal.load('astrodrizzle', defaults=True)  # get all default values
        parObj['build'] = True
        parObj['runfile'] = 'drizzlepac.run'
        parObj['STATE OF INPUT FILES']['preserve'] = False
        parObj['STATE OF INPUT FILES']['clean'] = True
        parObj['STEP 1: STATIC MASK']['static_sig'] = 3.0
        parObj['STEP 2: SKY SUBTRACTION']['skywidth'] = 0.1
        parObj['STEP 2: SKY SUBTRACTION']['skylower'] = -50.0
        parObj['STEP 2: SKY SUBTRACTION']['skyupper'] = 200.0
        parObj['STEP 3: DRIZZLE SEPARATE IMAGES']['driz_sep_bits'] = 8578
        parObj['STEP 4: CREATE MEDIAN IMAGE']['combine_maskpt'] = 0.7
        parObj['STEP 4: CREATE MEDIAN IMAGE']['combine_nsigma'] = '6 3'
        parObj['STEP 4: CREATE MEDIAN IMAGE']['combine_nhigh'] = 1
        parObj['STEP 6: REMOVE COSMIC RAYS WITH DERIV, DRIZ_CR']['driz_cr_snr'] = '3.0 2.5'
        parObj['STEP 7: DRIZZLE FINAL COMBINED IMAGE']['final_bits'] = 8578
        parObj['STEP 7: DRIZZLE FINAL COMBINED IMAGE']['final_units'] = 'counts'

        astrodrizzle.AstroDrizzle(asn_file, configobj=parObj)

        # Compare results
        outputs = [('j8bt06011_drz.fits', 'reference_asn_regress.fits')]
        self.compare_outputs(outputs)
