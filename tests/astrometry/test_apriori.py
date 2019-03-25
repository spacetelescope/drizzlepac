import pytest

import numpy as np

from drizzlepac.hlautils import testutils

from ..resources import BaseACS


class TestAcsApriori(BaseACS):

    @pytest.mark.parameterize('dataset', ['J9I408010'])#, 'J9AV01080', 'JCZGM1010'])
    def test_acs_apriori(self, dataset):
        """This test will perform fits between ALL a priori solutions and GAIA.

        """
        # Perform alignment of all WCS solutions with GAIA
        results_dict = testutils.compare_wcs_alignment(dataset)

        # Now compare results to see whether the a priori solutions succeeded
        # in improving the astrometry compared to default telescope pointing
        # reported by pipeline-defined WCS (IDC_* solution)
        wcsnames = list(results_dict.keys())
        idc_name = None
        for wcs in wcsnames:
            if 'IDC' in wcs and '-' not in wcs:
                idc_name = wcs
                wcsnames.remove(idc_name)
                break
        if idc_name is None:
            raise ValueError

        # Define pipeline-default fit
        pipeline_results = results_dict[idc_name]
        pipeline_offset = np.sqrt(pipeline_results['offset_x']**2 + pipeline_results['offset_y']**2)

        # compare each WCS fit results with pipeline-default fit
        for wcs in wcsnames:
            results = results_dict[wcs]
            # Check that fit for this WCS was successful
            print("Comparing fit for WCS='{}'...".format(wcs))
            assert (results['status'] == 0).all()
            assert (results['fit_qual'] < 2).all()

            offset = np.sqrt(results['offset_x']**2+results['offset_y']**2)
            assert (offset < pipeline_offset).all()
