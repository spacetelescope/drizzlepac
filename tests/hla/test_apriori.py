import pytest

import numpy as np

from drizzlepac.hlautils import testutils

from ..resources import BaseACS


class TestAcsApriori(BaseACS):

    @pytest.mark.parametrize('dataset', ['j95y04010', 'J9I408010',
                                         'ica9t3020', 'icnw34040','ID6Y05010'])
    def test_acs_apriori(self, dataset):
        """This test will perform fits between ALL a priori solutions and GAIA.

        Success criteria:
          * Successful fit to GAIA for dataset
          * Fit was not compromised (fit_qual != 5)
          * radial offset for new WCS < radial offset for IDC_* WCS
          * scale for new WCS within 0.1% of scale for IDC_* WCS
          * rotation for new WCS is within 0.1% of rotation from IDC_* WCS

        ASSUMPTIONS:
          * OPUS-based WCS solutions are ignored
          * Some WCS's may be defined for one image in an ASN but not another,
            in which case, that WCS is ignored (silently).

        """
        # Perform alignment of all WCS solutions with GAIA
        results_dict = testutils.compare_wcs_alignment(dataset)
        limit = 0.001

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

        success = False
        # compare each WCS fit results with pipeline-default fit
        for wcs in wcsnames:
            print("Comparing fit for WCS='{}'...".format(wcs), end=' ')
            results = results_dict[wcs]

            # Check that fit for this WCS was successful
            status = (results['status'] == 0).all()
            # check that fit was not compromised or otherwise invalid
            fit_qual =  (results['fit_qual'] < 5).all()

            # Check radial offset for this WCS compared to radial offset for IDC* WCS
            offset = np.sqrt(results['offset_x']**2+results['offset_y']**2)
            delta = (np.abs(offset - pipeline_offset) < 1).all() and (offset > pipeline_offset).all()
            delta = delta or (offset < pipeline_offset).all()

            # Check that rotation and scale are within
            delta_rot = np.abs(results['rotation']-pipeline_results['rotation'])/pipeline_results['rotation']
            delta_scale = np.abs(results['scale']-pipeline_results['scale'])/pipeline_results['scale']
            rot = (delta_rot < limit).all()
            scale = (delta_scale < limit).all()

            # Determine success/failure of this dataset's fit
            wcs_success = status and fit_qual and delta and rot and scale
            if wcs_success:
                # If at least one WCS succeeds, overall success needs to be set to True
                success = True
                print("SUCCESS")
            else:
                print("FAILED  due to:")
                if not status:
                    print("\t* invalid STATUS")
                if not fit_qual:
                    print("\t* invalid fit quality")
                if not delta:
                    print("\t* increased offset from GAIA.")
                if not rot:
                    print("\t* larger rotation from fit.")
                if not scale:
                    print("\t* larger scale from fit.")

        assert success
