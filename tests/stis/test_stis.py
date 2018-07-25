import os
import shutil

from stsci.tools import teal
from drizzlepac import astrodrizzle
from drizzlepac import pixtopix, pixtosky, skytopix
from stwcs import updatewcs

from ..resources import BaseSTIS


class TestSTIS(BaseSTIS):

    def test_fuv_mama(self):
        """ This test confirms that drizzlepac can correcly apply the
        distortion model for STIS FUV MAMA data and create a combined product.
        """

        # Prepare input files.
        raw_inputs = ['o60q02f3q_flt.fits', 'o60q02f4q_flt.fits',
                      'o60q02f6q_flt.fits', 'o60q02f8q_flt.fits',
                      'o60q02fbq_flt.fits', 'o60q02fcq_flt.fits',
                      'o60q02feq_flt.fits', 'o60q02fgq_flt.fits',
                      'o60q02fiq_flt.fits']
        all_inputs = [self.get_input_file('input', i) for i in raw_inputs]
        inputs = [os.path.basename(i) for i in all_inputs]
        print("[STIS_01] Found inputs: \n{}".format(inputs))

        output = 'final_stis_01'
        outfile = '{}_drz.fits'.format(output)
        reffile = 'reference_stis_01.fits'

        # Update WCS for all inputs
        updatewcs.updatewcs(inputs, use_db=False)

        # run astrodrizzle now...
        adriz_parobj = teal.load('astrodrizzle', defaults=True)
        adriz_parobj['output'] = output
        adriz_parobj['build'] = True
        adriz_parobj['in_memory'] = True
        adriz_parobj['runfile'] = 'stis_01'
        adriz_parobj['STATE OF INPUT FILES']['preserve'] = False
        adriz_parobj['STATE OF INPUT FILES']['clean'] = True
        adriz_parobj['STEP 1: STATIC MASK']['static_sig'] = 3.0
        adriz_parobj['STEP 2: SKY SUBTRACTION']['skystat'] = 'mean'
        adriz_parobj['STEP 2: SKY SUBTRACTION']['skywidth'] = 0.1
        adriz_parobj['STEP 2: SKY SUBTRACTION']['skylower'] = -50.0
        adriz_parobj['STEP 2: SKY SUBTRACTION']['skyupper'] = 200.0
        adriz_parobj['STEP 2: SKY SUBTRACTION']['use_static'] = False
        adriz_parobj['STEP 2: SKY SUBTRACTION']['sky_bits'] = None
        adriz_parobj['STEP 3: DRIZZLE SEPARATE IMAGES']['driz_separate'] = False
        adriz_parobj['STEP 4: CREATE MEDIAN IMAGE']['median'] = False
        adriz_parobj['STEP 5: BLOT BACK THE MEDIAN IMAGE']['blot'] = False
        adriz_parobj['STEP 6: REMOVE COSMIC RAYS WITH DERIV, DRIZ_CR']['driz_cr'] = False

        astrodrizzle.AstroDrizzle(inputs, configobj=adriz_parobj)

        # Compare results
        outputs = [(outfile, reffile)]
        self.compare_outputs(outputs)

    def test_nuv_mama(self):
        """ This test confirms that drizzlepac can correcly apply the
        distortion model for STIS NUV MAMA data and create a combined product.
        """

        # Prepare input files.
        raw_inputs = ['o6bz05pfq_flt.fits', 'o6bz05phq_flt.fits']
        all_inputs = [self.get_input_file('input', i) for i in raw_inputs]
        inputs = [os.path.basename(i) for i in all_inputs]

        output = 'final_stis_03'
        outfile = '{}_drz.fits'.format(output)
        reffile = 'reference_stis_03.fits'

        # Update WCS for all inputs
        updatewcs.updatewcs(inputs, use_db=False)

        # run astrodrizzle now...
        adriz_parobj = teal.load('astrodrizzle', defaults=True)
        adriz_parobj['output'] = output
        adriz_parobj['build'] = True
        adriz_parobj['in_memory'] = True
        adriz_parobj['runfile'] = 'stis_03'
        adriz_parobj['STATE OF INPUT FILES']['preserve'] = False
        adriz_parobj['STATE OF INPUT FILES']['clean'] = True
        adriz_parobj['STEP 1: STATIC MASK']['static_sig'] = 3.0
        adriz_parobj['STEP 2: SKY SUBTRACTION']['skystat'] = 'mean'
        adriz_parobj['STEP 2: SKY SUBTRACTION']['skywidth'] = 0.1
        adriz_parobj['STEP 2: SKY SUBTRACTION']['skylower'] = -50.0
        adriz_parobj['STEP 2: SKY SUBTRACTION']['skyupper'] = 200.0
        adriz_parobj['STEP 2: SKY SUBTRACTION']['use_static'] = False
        adriz_parobj['STEP 2: SKY SUBTRACTION']['sky_bits'] = None
        adriz_parobj['STEP 3: DRIZZLE SEPARATE IMAGES']['driz_separate'] = False
        adriz_parobj['STEP 4: CREATE MEDIAN IMAGE']['median'] = False
        adriz_parobj['STEP 5: BLOT BACK THE MEDIAN IMAGE']['blot'] = False
        adriz_parobj['STEP 6: REMOVE COSMIC RAYS WITH DERIV, DRIZ_CR']['driz_cr'] = False

        astrodrizzle.AstroDrizzle(inputs, configobj=adriz_parobj)

        # Compare results
        outputs = [(outfile, reffile)]
        self.compare_outputs(outputs)

    def test_stis_ccd(self):
        """ This test confirms that drizzlepac can correcly apply the
        distortion model for STIS CCD data and create a combined product.
        """

        # Prepare input files.
        raw_inputs = ['o6cl10arq_flt.fits', 'o6cl10asq_flt.fits',
                      'o6cl10atq_flt.fits', 'o6cl10auq_flt.fits',
                      'o6cl10axq_flt.fits', 'o6cl10ayq_flt.fits',
                      'o6cl10azq_flt.fits', 'o6cl10b0q_flt.fits']

        all_inputs = [self.get_input_file('input', i) for i in raw_inputs]
        inputs = [os.path.basename(i) for i in all_inputs]

        output = 'final_stis_04'
        outfile = '{}_drz.fits'.format(output)
        reffile = 'reference_stis_04.fits'

        # Update WCS for all inputs
        updatewcs.updatewcs(inputs, use_db=False)

        # run astrodrizzle now...
        adriz_parobj = teal.load('astrodrizzle', defaults=True)
        adriz_parobj['output'] = output
        adriz_parobj['build'] = True
        adriz_parobj['in_memory'] = True
        adriz_parobj['runfile'] = 'stis_04'
        adriz_parobj['STATE OF INPUT FILES']['preserve'] = False
        adriz_parobj['STATE OF INPUT FILES']['clean'] = True
        adriz_parobj['STEP 1: STATIC MASK']['static_sig'] = 3.0
        adriz_parobj['STEP 2: SKY SUBTRACTION']['skywidth'] = 0.1
        adriz_parobj['STEP 2: SKY SUBTRACTION']['skylower'] = -50.0
        adriz_parobj['STEP 2: SKY SUBTRACTION']['skyupper'] = 200.0
        adriz_parobj['STEP 2: SKY SUBTRACTION']['use_static'] = False
        adriz_parobj['STEP 2: SKY SUBTRACTION']['sky_bits'] = None
        adriz_parobj['STEP 4: CREATE MEDIAN IMAGE']['combine_maskpt'] = 0.7
        adriz_parobj['STEP 4: CREATE MEDIAN IMAGE']['combine_type'] = 'median'
        adriz_parobj['STEP 4: CREATE MEDIAN IMAGE']['combine_nsigma'] = '6 3'
        adriz_parobj['STEP 4: CREATE MEDIAN IMAGE']['combine_nhigh'] = 1
        adriz_parobj['STEP 7: DRIZZLE FINAL COMBINED IMAGE']['final_bits'] = 528

        astrodrizzle.AstroDrizzle(inputs, configobj=adriz_parobj)

        # Compare results
        outputs = [(outfile, reffile)]
        self.compare_outputs(outputs)

    def test_stis_oiii_ccd(self):
        """ This test confirms that drizzlepac can correcly apply the
        distortion model for STIS F28x50OIII CCD data and create a
        combined product.
        """

        # Prepare input files.
        raw_inputs = ['o47m05f8q_flt.fits', 'o47m05f9q_flt.fits',
                      'o47m05faq_flt.fits', 'o47m05fbq_flt.fits']
        all_inputs = [self.get_input_file('input', i) for i in raw_inputs]
        inputs = [os.path.basename(i) for i in all_inputs]

        output = 'final_stis_09'
        outfile = '{}_drz.fits'.format(output)
        reffile = 'reference_stis_09.fits'

        # Update WCS for all inputs
        updatewcs.updatewcs(inputs, use_db=False)

        # run astrodrizzle now...
        adriz_parobj = teal.load('astrodrizzle', defaults=True)
        adriz_parobj['output'] = output
        adriz_parobj['build'] = True
        adriz_parobj['runfile'] = 'stis_09'
        adriz_parobj['STATE OF INPUT FILES']['preserve'] = False
        adriz_parobj['STATE OF INPUT FILES']['clean'] = False
        adriz_parobj['STEP 1: STATIC MASK']['static_sig'] = 3.0
        adriz_parobj['STEP 2: SKY SUBTRACTION']['skysub'] = False
        adriz_parobj['STEP 4: CREATE MEDIAN IMAGE']['combine_maskpt'] = 0.7
        adriz_parobj['STEP 4: CREATE MEDIAN IMAGE']['combine_type'] = 'median'
        adriz_parobj['STEP 4: CREATE MEDIAN IMAGE']['combine_nsigma'] = '6 3'
        adriz_parobj['STEP 4: CREATE MEDIAN IMAGE']['combine_nhigh'] = 1
        adriz_parobj['STEP 7: DRIZZLE FINAL COMBINED IMAGE']['final_bits'] = 528

        astrodrizzle.AstroDrizzle(inputs, configobj=adriz_parobj)

        # Compare results
        outputs = [(outfile, reffile)]
        self.compare_outputs(outputs)
