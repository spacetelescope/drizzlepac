import os
import shutil
import pytest

from stsci.tools import teal
from stwcs import updatewcs
from drizzlepac import astrodrizzle, tweakreg
from drizzlepac import pixtopix, pixtosky, skytopix

from ..resources import BaseACS


class TestAcsTweak(BaseACS):

    #@pytest.mark.xfail
    def test_tweak(self):
        """ This test confirms that image alignment performed by tweakreg
        works to find the correct alignment solution and that the solution
        can be used to correctly update the headers of the inputs to generate
        a properly aligned, combined mosaic using astrodrizzle.
        """
        input_file1 = 'j94f05bgq_flt.fits'
        input_file2 = 'j9irw1rqq_flt.fits'
        output = 'acs_tweak'

        # Prepare input files.
        input_file1 = self.get_input_file('input', input_file1)
        input_file2 = self.get_input_file('input', input_file2)
        inputs = ','.join([input_file1, input_file2])

        tweak_parObj = teal.load('tweakreg', defaults=True)
        tweak_parObj['input'] = 'j94f05bgq_flt.fits,j9irw1rqq_flt.fits'
        tweak_parObj['writecat'] = False
        tweak_parObj['clean'] = True
        tweak_parObj['OBJECT MATCHING PARAMETERS']['searchrad'] = 3.0
        tweak_parObj['OBJECT MATCHING PARAMETERS']['see2dplot'] = False
        tweak_parObj['OBJECT MATCHING PARAMETERS']['separation'] = 0.0
        tweak_parObj['CATALOG FITTING PARAMETERS']['residplot'] = 'No plot'
        tweak_parObj['CATALOG FITTING PARAMETERS']['nclip'] = 5
        tweak_parObj['CATALOG FITTING PARAMETERS']['sigma'] = 2.8
        imfind_parObj = teal.load('imagefindpars', defaults=True)
        imfind_parObj['computesig'] = False
        imfind_parObj['skysigma'] = 20.0
        imfind_parObj['conv_width'] = 2.0
        imfind_parObj['threshold'] = 500.0
        imfind_parObj['dqbits'] = None
        rfind_parObj = teal.load('refimagefindpars', defaults=True)
        rfind_parObj['computesig'] = False
        rfind_parObj['skysigma'] = 20.0
        rfind_parObj['conv_width'] = 2.0
        rfind_parObj['threshold'] = 500.0
        rfind_parObj['dqbits'] = None

        tweakreg.TweakReg(files=inputs, configobj=tweak_parObj,
                          imagefindcfg=imfind_parObj,
                          refimagefindcfg=rfind_parObj,
                          updatehdr=True, wcsname='TWEAK_regress')

        # run astrodrizzle now...
        adriz_parObj = teal.load('astrodrizzle', defaults=True)  # get all default values
        adriz_parObj['output'] = 'acs_tweak'
        adriz_parObj['build'] = True
        adriz_parObj['in_memory'] = True
        adriz_parObj['runfile'] = 'drizzlepac.run'
        adriz_parObj['STATE OF INPUT FILES']['preserve'] = False
        adriz_parObj['STATE OF INPUT FILES']['clean'] = True
        adriz_parObj['STEP 2: SKY SUBTRACTION']['skywidth'] = 0.3
        adriz_parObj['STEP 2: SKY SUBTRACTION']['use_static'] = False
        adriz_parObj['STEP 2: SKY SUBTRACTION']['sky_bits'] = None
        adriz_parObj['STEP 4: CREATE MEDIAN IMAGE']['combine_maskpt'] = 0.7

        astrodrizzle.AstroDrizzle(inputs, configobj=adriz_parObj)

        # Compare results
        outfile = '{}_drz.fits'.format(output)
        reffile = 'reference_tweak.fits'
        outputs = [(outfile, reffile)]
        for i in range(1,5):
            self.ignore_keywords += ['D00{}DATA'.format(i), 'D00{}MASK'.format(i)]
        self.compare_outputs(outputs)
