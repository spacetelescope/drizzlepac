import os

from stsci.tools import teal
import drizzlepac
from drizzlepac import astrodrizzle
from drizzlepac.helpers.mark import require_bigdata

from ..helpers import BaseACS, raw_from_asn


class TestAsnNarrowband(BaseACS):
    
    def test_acs_narrowband(self):
        rootname = 'j8dw01010'
        asn_file = rootname + '_asn.fits'

        # Prepare input files.
        input_file = self.get_data('input', asn_file)

        for raw_file in raw_from_asn(asn_file, suffix='_flt.fits'):
            self.get_input_file(raw_file)

        # run astrodrizzle now...
        parObj = teal.load('astrodrizzle', defaults=True)  # get all default values
        parObj['build'] = True
        parObj['runfile'] = 'drizzlepac.run'
        parObj['STATE OF INPUT FILES']['preserve'] = False
        parObj['STATE OF INPUT FILES']['clean'] = True
        parObj['STEP 2: SKY SUBTRACTION']['use_static'] = False
        parObj['STEP 2: SKY SUBTRACTION']['sky_bits'] = None
        parObj['STEP 4: CREATE MEDIAN IMAGE']['combine_maskpt'] = 0.7
        parObj['STEP 4: CREATE MEDIAN IMAGE']['combine_nsigma'] = '6 3'
        parObj['STEP 4: CREATE MEDIAN IMAGE']['combine_nhigh'] = 1
        parObj['STEP 6: REMOVE COSMIC RAYS WITH DERIV, DRIZ_CR']['driz_cr_snr'] = '3.0 2.5'
        parObj['STEP 7: DRIZZLE FINAL COMBINED IMAGE']['final_bits'] = 8578
        parObj['STEP 7a: CUSTOM WCS FOR FINAL OUTPUT']['final_wcs'] = True
        parObj['STEP 7a: CUSTOM WCS FOR FINAL OUTPUT']['final_rot'] = 0.0
        
        astrodrizzle.AstroDrizzle(asn_file, configobj=parObj)

        # Compare results
        outputs = [('j8dw01010_drz.fits', 'reference_narrowband.fits')]
        self.compare_outputs(outputs)
