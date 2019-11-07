"""
Tools for aperture photometry with non native bg/error methods

This function serves to ease the computation of photometric magnitudes
and errors using PhotUtils by replicating DAOPHOT's photometry and
error methods.  The formula for DAOPHOT's error is:

err = sqrt (Poisson_noise / epadu + area * stdev**2 + area**2 * stdev**2 / nsky)

Which gives a magnitude error:

mag_err = 1.0857 * err / flux

Where epadu is electrons per ADU (gain), area is the photometric
aperture area, stdev is the uncertainty in the sky measurement
and nsky is the sky annulus area.  To get the uncertainty in the sky
we must use a custom background tool, which also enables computation of
the mean and median of the sky as well (more robust statistics).
All the stats are sigma clipped.  These are calculated by the
functions in aperture_stats_tbl.

.. note::
    Currently, the background computations will fully include a pixel that has ANY overlap with the background aperture
    (the annulus). This is to simplify the computation of the median, as a weighted median is nontrivial, and slower.
    Copied from https://grit.stsci.edu/HLA/software/blob/master/HLApipeline/HLApipe/scripts/photometry_tools.py
Authors
-------
    - Varun Bajaj, January 2018

Use
---

::

    from photometry_tools import iraf_style_photometry
    phot_aps = CircularAperture((sources['xcentroid'], sources['ycentroid']),r=10.)
    bg_aps = CircularAnnulus((sources['xcentroid'], sources['ycentroid']), r_in=13., r_out=15.)

Simplest call:

::

    photometry_tbl = iraf_style_photometry(phot_aps, bg_aps, data)

Pass in pixelwise error and set background to mean

::

    photometry_tbl = iraf_style_photometry(phot_aps, bg_aps, data, error_array=data_err, bg_method='mean')

Can also set the gain (if image units are DN)

::

    photometry_tbl = iraf_style_photometry(phot_aps, bg_aps, data, epadu=2.5)

Classes and Functions
---------------------
"""
import numpy as np
from astropy.table import Table
from drizzlepac.hlautils.background_median import aperture_stats_tbl
from photutils import aperture_photometry


def iraf_style_photometry(phot_apertures, bg_apertures, data, photflam, photplam, error_array=None,
                          bg_method='mode', epadu=1.0):
    """
    Computes photometry with PhotUtils apertures, with IRAF formulae

    Parameters
    ----------
    phot_apertures : photutils PixelAperture object (or subclass)
        The PhotUtils apertures object to compute the photometry. i.e. the object returned via CirularAperture.

    bg_apertures : photutils PixelAperture object (or subclass)
        The phoutils aperture object to measure the background in. i.e. the object returned via CircularAnnulus.

    data : array
        The data for the image to be measured.

    photflam : float
        inverse sensitivity, in ergs/cm2/angstrom/electron

    photplam : float
        Pivot wavelength, in angstroms

    error_array : array
        (Optional) The array of pixelwise error of the data.  If none, the Poisson noise term in the error computation
        will just be the square root of the flux/epadu. If not none, the aperture_sum_err column output by
        aperture_photometry (divided by epadu) will be used as the Poisson noise term.

    bg_method: string
        {'mean', 'median', 'mode'}, optional. The statistic used to calculate the background. All measurements are
        sigma clipped. Default value is 'mode'. NOTE: From DAOPHOT, mode = 3 * median - 2 * mean.

    epadu : float
        (optional) Gain in electrons per adu (only use if image units aren't e-). Default value is 1.0

    Returns
    -------
        An astropy Table with columns as follows:
        X-Center Y-Center RA DEC ID MagAp1 MagErrAp1 MagAp2 MagErrAp2 MSkyAp2 StdevAp2 FluxAp2 CI Flags
    """
    if bg_method not in ['mean', 'median', 'mode']:
        raise ValueError('Invalid background method, choose either \
                          mean, median, or mode')
    phot = aperture_photometry(data, phot_apertures, error=error_array)
    bg_phot = aperture_stats_tbl(data, bg_apertures, sigma_clip=True)
    names = ['X-Center', 'Y-Center', 'ID']
    x, y = phot_apertures[0].positions.T
    final_stacked = np.stack([x, y, phot["id"].data], axis=1)
    # n_aper = 0
    name_list = 'Flux', 'FluxErr', 'Mag', 'MagErr'
    for aper_string  in ['Ap1', 'Ap2']:
        for col_name in name_list:
            names.append("{}{}".format(col_name,aper_string))

    # for item in list(phot.keys()):
    #     if item.startswith("aperture_sum_") and not item.startswith("aperture_sum_err_"):
    #         aper_size_arcsec = phot_apertures[n_aper].r * platescale
    #         for name in name_list:
    #             names.append("{}_{:.2f}".format(name, aper_size_arcsec))
    #         n_aper += 1
    for aperCtr in range(0, 2):
        ap_area = phot_apertures[aperCtr].area
        bg_method_name = 'aperture_{}'.format(bg_method)

        # NOTE background subtraction below commented out 8/14/19
        flux = phot['aperture_sum_{}'.format(aperCtr)]  # - bg_phot[bg_method_name] * ap_area

        # Need to use variance of the sources
        # for Poisson noise term in error computation.
        #
        # This means error needs to be squared.
        # If no error_array error = flux ** .5

        if error_array is not None:
            flux_error = compute_phot_error(phot['aperture_sum_err_{}'.format(aperCtr)] ** 2.0, bg_phot, bg_method,
                                            ap_area, epadu)
        else:
            flux_error = compute_phot_error(flux, bg_phot, bg_method, ap_area, epadu)

        mag = convert_flux_to_abmag(flux, photflam, photplam)

        # NOTE: Magnitude error calculation comes from computing d(ABMAG)/d(flux).
        # See https://iraf.net/forum/viewtopic.php?showtopic=83932 for details.
        mag_err = 1.0857 * flux_error / flux

        # Build the final data table
        stacked = np.stack([flux, flux_error, mag, mag_err], axis=1)
        final_stacked = np.concatenate([final_stacked, stacked], axis=1)

    # Build final output table
    final_tbl = Table(data=final_stacked, names=names,
                      dtype=[np.float64, np.float64, np.int64, np.float64, np.float64, np.float64, np.float64,
                             np.float64, np.float64, np.float64, np.float64])

    # add sky and std dev columns from background calculation subroutine
    final_tbl.add_column(bg_phot[bg_method_name])
    final_tbl.rename_column(bg_method_name, 'MSkyAp2')
    final_tbl.add_column(bg_phot['aperture_std'])
    final_tbl.rename_column('aperture_std', 'StdevAp2')

    return final_tbl


def compute_phot_error(flux_variance, bg_phot, bg_method, ap_area, epadu=1.0):
    """Computes the flux errors using the DAOPHOT style computation

    Parameters
    ----------
    flux_variance : array
        flux values

    bg_phot : array
        background brightness values.

    bg_method : string
        background method

    ap_area : array
        the area of the aperture in square pixels

    epadu : float
        (optional) Gain in electrons per adu (only use if image units aren't e-). Default value is 1.0

    Returns
    -------
    flux_error : array
        an array of flux errors
    """

    bg_variance_terms = (ap_area * bg_phot['aperture_std'] ** 2.) * (1. + ap_area/bg_phot['aperture_area'])
    variance = flux_variance / epadu + bg_variance_terms
    flux_error = variance ** .5
    return flux_error


def convert_flux_to_abmag(in_flux, photflam, photplam):
    """converts flux (in units of electrons/sec) to ABMAG

    Parameters
    ----------
    in_flux : astropy.table.column.Column object
        flux values to convert to ABMAG, in electrons/second

    photflam : float
        inverse sensitivity, in ergs/cm2/angstrom/electron

    photplam : float
        pivot wavelength, in angstroms

    Returns
    -------
    abmag : astropy.table.column.Column object
        input flux values converted to ABMAG
    """

    # convert flux from units of electrons/second to ergs/cm2/angstrom/second
    f_lambda = in_flux * photflam

    # Convert f_lambda to STMAG
    stmag = -2.5 * np.log10(f_lambda) - 21.10

    # Convert STMAG to ABMAG
    abmag =  stmag - 5.0 * np.log10(photplam) + 18.6921

    return abmag