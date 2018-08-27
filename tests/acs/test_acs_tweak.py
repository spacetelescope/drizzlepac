import os
import shutil

from stsci.tools import teal
from stwcs import updatewcs
from drizzlepac import astrodrizzle, tweakreg
from drizzlepac import pixtopix, pixtosky, skytopix

from ..resources import BaseACS


class TestAcsTweak(BaseACS):

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

    def test_pixsky1(self):
        """This test verifies that the coordinate transformation tasks
        'pixtopix', 'pixtosky', and 'skytopix' still work as expected.

        This test relies on the truth/comparison file from the `acs_tweak` test
        to define the output reference WCS for these coordinate transformations.
        """

        input_file = self.get_input_file('input', 'j94f05bgq_flt.fits')
        fpath, froot = os.path.split(input_file)

        refwcs_file = self.get_data('truth', 'reference_tweak.fits')
        ref_file1 = 'j94f05bgq_flt_reference.coo'
        ref_file2 = 'j94f05bgq_flt_sky_reference.coo'
        ref_file3 = 'j94f05bgq_flt_sci1_skyxy_reference.coo'
        ref_file4 = 'j94f05bgq_flt_sci2_skyxy_reference.coo'

        # Run test...
        flt_file = updatewcs.updatewcs(froot)
        # create combined results file
        out_refxy_file = froot.replace('_flt.fits', '_flt_refxy.coo')
        out_sky_file = froot.replace('_flt.fits', '_flt_sky.coo')
        inskycat = froot.replace('_flt.fits', '_flt_sky_catalog.coo')
        inskycat = self.get_data('input', inskycat)
        outpixcat_names = []

        outpix = open(out_refxy_file, 'w')
        outsky = open(out_sky_file, 'w')
        for i in [1, 2]:
            outcat = froot.replace('_flt.fits', '_flt_sci%d_refxy.coo' % i)
            incat = froot.replace('_flt.fits', '_flt_sci%d_xy_catalog.coo' % i)
            incat = self.get_data('input', incat)
            pixtopix.tran(froot+'[sci,%d]' % i, '{}[sci,1]'.format(refwcs_file),
                          direction='forward',
                          x=None, y=None, coordfile=incat, colnames='',
                          separator=None, precision=6, output=outcat,
                          verbose=False)

            shutil.copyfileobj(open(outcat, 'r'), outpix)

            outskycat = froot.replace('_flt.fits', '_flt_sci%d_sky.coo' % i)
            pixtosky.xy2rd(froot+'[sci,%d]' % i, coordfile=incat,
                           x=None, y=None, colnames='', separator=None,
                           hms=True, precision=6, output=outskycat,
                           verbose=False)
            shutil.copyfileobj(open(outskycat, 'r'), outsky)

            outpixcat = froot.replace('_flt.fits', '_flt_sci%d_skyxy.coo' % i)
            skytopix.rd2xy(froot+'[sci,%d]' % i, ra=None, dec=None,
                           coordfile=inskycat, colnames=None, precision=6,
                           output=outpixcat, verbose=False)
            outpixcat_names.append(outpixcat)

        # Close combined results files
        outpix.close()
        outsky.close()

        outputs = [(out_refxy_file, ref_file1), (out_sky_file, ref_file2)]
        for oname, rname in zip(outpixcat_names, [ref_file3, ref_file4]):
            outputs.append((oname, rname))

        self.compare_outputs(outputs)
