"""
Computes MMM statistics within photutils apertures.

The functions in this script enable the computation of statistics
within a PhotUtils aperture, which is currently not directly
implemented in PhotUtils itself.  This code is meant to be
imported into other code, and then be usable as a single line to
return all the statistics in a format similar to the
aperture_photometry method in PhotUtils (i.e. an astropy table).

Authors
-------
    - Varun Bajaj, December 2017

Github VCS Info
---------------
This script was borrowed with permission from the following Github repository:

* `spacetelescope/wfc3_photometry/photometry_tools <https://github.com/spacetelescope/wfc3_photometry/tree/master/photometry_tools>`_ at revision 690118f
* copied from https://grit.stsci.edu/HLA/software/blob/master/HLApipeline/HLApipe/scripts/background_median.py

Path
----
HLApipeline/HLApipe/scripts/background_median.py

Dependencies
------------
None.

Inputs
------
None.

Use
---
from background_median import aperture_stats_tbl

stats_tbl = aperture_stats_tbl(data, apertures)

.. note::
    * See the docstring of aperture_stats_tbl for more info.

Classes and Functions
---------------------
"""
import numpy as np

# WAY faster than astropy.stats.sigma_clipped_stats
from scipy.stats import sigmaclip
from astropy.table import Table

def aperture_stats_tbl(data, apertures, method='exact', sigma_clip=True):
    """Computes mean/median/mode/std in Photutils apertures.

    Compute statistics for custom local background methods. This is primarily intended for estimating backgrounds via
    annulus apertures. The intent is that this falls easily into other code to provide background measurements.

    Parameters
    ----------
    data : array
        The data for the image to be measured.

    apertures : photutils PixelAperture object (or subclass)
        The phoutils aperture object to measure the stats in. i.e. the object returned via CirularAperture,
        CircularAnnulus, or RectangularAperture etc.

    method : sting
        The method by which to handle the pixel overlap. Defaults to computing the exact area. NOTE: Currently, this
        will actually fully include a pixel where the aperture has ANY overlap, as a median is also being performed.
        If the method is set to 'center' the pixels will only be included if the pixel's center falls within the
        aperture.

    sigma_clip : Boolean
        Flag to activate sigma clipping of background pixels

    Returns
    -------
    stats_tbl : astropy table
        An astropy Table with the colums X, Y, aperture_mean, aperture_median, aperture_mode, aperture_std,
        aperture_area and a row for each of the positions of the apertures.
    """

    # Get the masks that will be used to identify our desired pixels.
    masks = apertures.to_mask(method=method)

    # Compute the stats of pixels within the masks
    aperture_stats = [calc_aperture_mmm(data, mask, sigma_clip) for mask in masks]

    aperture_stats = np.array(aperture_stats)


    # Place the array of the x y positions alongside the stats
    stacked = np.hstack([apertures.positions, aperture_stats])
    # Name the columns
    names = ['X','Y','aperture_mean','aperture_median','aperture_mode', 'aperture_std', 'aperture_area']
    # Make the table
    stats_tbl = Table(data=stacked, names=names)


    return stats_tbl

def calc_aperture_mmm(data, mask, sigma_clip):
    """Helper function to actually calculate the stats for pixels falling within some Photutils aperture mask on some
    array of data.

    Parameters
    ----------
    data : array
        The data for the image to be measured.

    mask : array
        mask that will be used to identify our desired pixels

    sigma_clip: Boolean
        Flag to activate sigma clipping of background pixels

    Returns
    -------
    mean : flaot
        Mean pixel value

    median : float
        Median pixel value

    mode : float
        mode pixel value

    std : float
        Standard deviation of pixel values

    actual_area : float
        Calculated enclosed area
    """
    cutout = mask.cutout(data, fill_value=np.nan)
    if cutout is None:
        return (np.nan, np.nan, np.nan, np.nan, np.nan)
    else:
        values = cutout * mask.data / mask.data
        values = values[~np.isnan(values)]
        if sigma_clip:
            values, clow, chigh = sigmaclip(values, low=3, high=3)

        mean = np.mean(values)
        median = np.median(values)
        std = np.std(values)

        mode = 3 * median - 2 * mean
        actual_area = (~np.isnan(values)).sum()
        return (mean, median, mode, std, actual_area)
