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
