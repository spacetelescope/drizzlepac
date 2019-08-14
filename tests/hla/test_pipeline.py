""" This module tests full pipeline use of the drizzlepac package.

"""
import os
import pytest

from drizzlepac.hlautils import astroquery_utils as aqutils
from drizzlepac import runastrodriz
from astropy.io import fits



def test_astrometric_singleton():
    """ Tests pipeline-style processing of a singleton exposure using runastrodriz.
    """
    dataset_names = ['iaaua1n4q']
    # Get sample data through astroquery
    flcfile = aqutils.retrieve_observation(dataset_names, suffix=['FLC'])[0]
    fltfile = aqutils.retrieve_observation(dataset_names, suffix=['FLT'])[0]
    rawfile = aqutils.retrieve_observation(dataset_names, suffix=['RAW'])[0]

    # Insure environment variables are set for full processing
    os.environ['ASTROMETRY_STEP_CONTROL'] = 'on'
    os.environ['ASTROMETRY_COMPUTE_APOSTERIORI'] = 'on'
    os.environ['ASTROMETRY_APPLY_APRIORI'] = 'on'

    # Run pipeline processing using
    runastrodriz.process(rawfile, force=True, inmemory=True)

    # compare WCSNAMEs from flt and flc files
    flc_wcsname = fits.getval(flcfile, 'wcsname', ext=1)
    flt_wcsname = fits.getval(fltfile, 'wcsname', ext=1)

    # Perform comparisons:
    #   - WCSNAME values should contain '-' from either a priori or a posteriori solution
    #   - WCSNAME value should be the same for FLT and FLC images
    assert('-' in flc_wcsname)
    assert('-' in flt_wcsname)
    assert(flc_wcsname == flt_wcsname)
    
