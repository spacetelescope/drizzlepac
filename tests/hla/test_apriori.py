""" This module computes alignment solutions between all "a priori" solutions
    for a dataset and GAIA.

"""
import pytest

import numpy as np

from drizzlepac.hlautils import testutils

from ..resources import BaseACS, BaseWFC3


def compare_apriori(dataset):
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
    for idc_name in wcsnames:
        if 'IDC' in idc_name and '-' not in idc_name:
            wcsnames.remove(idc_name)
            break
    else:
        raise ValueError
    # Define pipeline-default fit
    pipeline_results = results_dict[idc_name]
    pipeline_offset = np.sqrt(pipeline_results['offset_x']**2 +
                              pipeline_results['offset_y']**2)

    success = False
    # compare each WCS fit results with pipeline-default fit
    for wcs in wcsnames:
        print("Comparing fit for WCS='{}'...".format(wcs), end=' ')
        results = results_dict[wcs]

        # Check that fit for this WCS was successful
        status = (results['status'] == 0).all()
        # check that fit was not compromised or otherwise invalid
        fit_qual = (results['fit_qual'] < 5).all()

        # Check radial offset for this WCS compared to radial offset for
        # IDC* WCS
        offset = np.sqrt(results['offset_x']**2 + results['offset_y']**2)
        delta = ((offset < pipeline_offset).all() or
                 np.allclose(offset, pipeline_offset, rtol=0, atol=1))

        # Check that rotation and scale are within
        rot = np.allclose(results['rotation'], pipeline_results['rotation'],
                          rtol=limit, atol=0)
        scale = np.allclose(results['scale'], pipeline_results['scale'],
                            rtol=limit, atol=0)

        # Determine success/failure of this dataset's fit
        if all([status, fit_qual, delta, rot, scale]):
            # If at least one WCS succeeds, overall success
            # needs to be set to True
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


class TestAcsApriori(BaseACS):
    """ Tests which validate whether mosaics can be aligned to an astrometric standard,
        evaluate the quality of the fit, and generate a new WCS.

        * These can only be run using the pytest option:

            pytest --bigdata

        * This test is also compatible with pytest-xdist.
    """

    @pytest.mark.bigdata
    @pytest.mark.parametrize('dataset', ['jb1601020', 'J9I408010'])
    def test_apriori(self, dataset):
        compare_apriori(dataset)


class TestWFC3Apriori(BaseWFC3):
    """ Tests which validate whether mosaics can be aligned to an astrometric
        standard, evaluate the quality of the fit, and generate a new WCS.

        * These can only be run using the pytest option:

            pytest --bigdata

        * This test is also compatible with pytest-xdist.

    """

    @pytest.mark.bigdata
    @pytest.mark.parametrize(
        'dataset', ['ic0g0l010', 'icnw34040']
    )
    def test_apriori(self, dataset):
        compare_apriori(dataset)
