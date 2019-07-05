import os
import shutil

from scipy import ndimage
from scipy.ndimage import morphology
import numpy as np
import astropy
from astropy.io import fits
from astropy.table import Table
from spherical_geometry.polygon import SphericalPolygon

from stwcs.wcsutil import HSTWCS

from .. import wcs_functions


# Default grid definition file
PCELL_FILENAME = 'allsky_cells.fits'

def get_sky_cells(visit_input, input_path=None, scale=None, cell_size=None):
    """Return all sky cells that overlap the exposures in the input.

    Parameters
    -----------
    visit_input : str or list
        Input specifying the exposures from a single visit; either
        a poller output file or a simple list of exposure filenames.
        Exposures in an input list are assumed to be in the current
        working directory when running the code, unless `input_path`
        has been provided which points to the location of the exposures
        to be processed.

    input_path : str, optional
        Location of input exposures, if provided.  If not provided,
        location will be assumed to be the current working directory.

    scale : float, optional
        User-defined value for the pixel scale, in arcseconds/pixel,
        of the sky cells and projection cells. If `None`, default
        value from grid definition file will be used.

    cell_size : float, optional
        User-specified size, in degrees, for each projection cell.
        If `None`, default value from grid definition file will be
        used.

    Returns
    --------
    sky_cells : list of objects
        List of `SkyCell` objects for all sky cells which overlap the
        exposures provided in `visit_input`.
    """
    # Interpret input
    if isinstance(visit_input, list):
        expnames = visit_input
    else:
        expnames = Table.read(visit_input, format='ascii.fast_no_header')[0]

    # Check that exposures are located in current working directory
    if not os.path.exists(expnames[0]):
        if not input_path:
            msg = "No exposures found in cwd().  Please specify path to files!"
            raise (ValueError, msg)
        bad_files = 0
        for file in expnames:
            fullfile = os.path.join(input_path, file)
            if not os.path.exists(fullfile):
                bad_files.append(fullfile)
                print("Could not find {}".format(fullfile))
                bad_files += 1
                continue
            shutil.copy(fullfile, file)
        if bad_files:
            msg = "Could not find {} specified input files".format(bad_files)
            raise (ValueError, msg)

    # Initialize all sky tessellation object definitions
    # This includes setting the pixel scale.
    sky_grid = AllSky(scale=scale, cell_size=cell_size)

    # build reference wcs for combined footprint of all input exposures
    meta_wcs = wcs_functions.make_mosaic_wcs(expnames, scale=sky_grid.scale)

    # create footprint on the sky (as a tangent plane array) for all input exposures using meta_wcs
    footprint = SkyFootprint(meta_wcs).build(expnames)
    footprint_hdu = footprint.get_footprint_hdu()

    # Use this footprint to identify overlapping sky cells
    sky_cells = sky_grid.find_sky_cells(footprint_hdu)

    return sky_cells


class SkyFootprint(object):

    def __init__(self, meta_wcs):

        self.struct2 = ndimage.generate_binary_structure(2, 2)
        self.meta_wcs = meta_wcs

        self.total_mask = np.zeros(meta_wcs.array_shape, dtype=np.int16)
        self.footprint = None

        self.edges = None
        self.edges_ra = None
        self.edges_dec = None
        self.polygon = None

    def build(self, expnames):
        for exposure in expnames:
            blank = np.zeros(self.meta_wcs.array_shape, dtype=np.int16)
            exp = fits.open(exposure)
            sci_extns = wcs_functions.get_extns(exp)
            for sci in sci_extns:
                wcs = HSTWCS(exp, ext=sci)
                edges = []
                edges += [[0, i] for i in range(wcs.naxis2)]
                edges += [[wcs.naxis1, i] for i in range(wcs.naxis2)]
                edges += [[i, 0] for i in range(wcs.naxis1)]
                edges += [[i, wcs.naxis2] for i in range(wcs.naxis1)]

                sky_edges = wcs.pixel_to_world_values(edges)
                meta_edges = self.meta_wcs.world_to_pixel_values(sky_edges).astype(np.int32)
                blank[meta_edges[:, 0], meta_edges[:, 1]] = 1

                # Fill in outline of each chip
                blank = morphology.binary_dilation(blank, structure=self.struct2)
                blank = morphology.binary_fill_holes(blank)
                blank = morphology.binary_erosion(blank, structure=self.struct2)

            self.total_mask += blank.astype(np.int16)

    # Methods with 'find' compute values
    # Methods with 'get' return values
    def find_footprint(self):
        self.footprint = np.clip(self.total_mask, 0, 1)

    def find_edges(self):
        if not self.footprint:
            self.find_footprint()
        edges = morphology.binary_erosion(self.footprint).astype(np.int16)
        self.edges = self.footprint - edges

    def get_edges_sky(self):
        if not self.edges:
            self.find_edges()
        edges = np.where(self.edges)
        self.edges_ra, self.edges_dec = self.meta_wcs.pixel_to_world_values(edges[0], edges[1])
        return self.edges_ra, self.edges_dec

    def find_polygon(self):
        if not self.edges_ra:
            self.get_edges_sky()
        self.polygon = SphericalPolygon.from_radec(self.edges_ra, self.edges_dec, self.meta_wcs.wcs.crval)

    def _get_fits_hdu(self, data, filename=None, overwrite=True):
        hdulist = fits.HDUList()
        phdu = fits.PrimaryHDU(data=data, header=self.meta_wcs.to_header())
        hdulist.append(phdu)
        if filename:
            hdulist.writeto(filename, overwrite=overwrite)
        return hdulist

    def get_footprint_hdu(self, filename=None, overwrite=True):
        if not self.footprint:
            self.find_footprint()
        return self._get_fits_hdu(self.footprint, filename=filename, overwrite=overwrite)

    def get_edges_hdu(self, filename=None, overwrite=True):
        if not self.edges:
            self.find_edges()
        return self._get_fits_hdu(self.edges, filename=filename, overwrite=overwrite)

    def get_mask_hdu(self, filename=None, overwrite=True):
        return self._get_fits_hdu(self.total_mask, filename=filename, overwrite=overwrite)


class AllSky(object):

    def __init__(self, scale=None, cell_size=None):
        """Setup tesselation based on installed definition of grid."""
        fpath = os.path.abspath(os.path.dirname(__file__))
        froot = os.path.join(os.path.dirname(fpath), 'pars')
        fname = os.path.join(froot, PCELL_FILENAME)

        self.hdu = fits.open(fname)
        self.scale = scale if scale else self.hdu[0].header['PCSCALE']
        self.cell_size = cell_size if cell_size else self.hdu[0].header['PCSIZE']

    def _find_bands(self, footprint):
        """ Select the band or bands which encompass the provided footprint.

        The footprint will be a ndarray mask with pixel value of 1 where the exposures are located,
        and with a fully defined WCS to convert pixel positions to sky coordinates.

        """
        # Interpret footprint to get range of declination in mask
        ra, dec = footprint.get_edges_sky()

        # Select band with projection cell
        rings = self.hdu[1].data
        # find dec zone where rings.dec_min <= dec < rings.dec_max
        bands = np.unique(np.searchsorted(rings.field('dec_max'), dec))

        # special handling at pole where overlap is complicated
        # do extra checks for top 2 rings
        # always start with the ring just below the pole
        nearpole = np.where(bands >= len(rings) - 2)
        bands[nearpole] = len(rings) - 2

        # Identify how may projection cells are in each overlapping band
        nbands = rings[bands].field('nband')

        # Record these values as attributes for use in other methods
        self.bands = bands
        self.nbands = nbands

    def get_sky_cells(self, footprint):

        # Find band[s] that overlap footprint
        self._find_bands(footprint)

        # Define numerical position in band for projection cell
        # self.band_index = self.projection_cell_id - self.band['PROJCELL']
        for band in self.bands:
            # compute band_index, one for each projection cell that overlaps the footprint
            band_indices = []
            pcells = [ProjectionCell(band_index, band) for band_index in band_indices]

        # Find sky cells from identified projection cell(s) that overlap footprint
        sky_cells = []
        for pcell in pcells:
            sky_indices = pcell.find_sky_cells(footprint)
            sky_cells += [SkyCell(sky_x, sky_y, pcell) for sky_x, sky_y in sky_indices]
        return sky_cells


class ProjectionCell(object):

    def __init__(self, index, band):
        """Build projection cell for cell with name `skycell_NNNNN`"""
        self.band_index = index
        self.band = band
        self.cell_id = str(band['PROJCELL'] + index).zfill(5)

        # Generate WCS for projection cell
        self._build_wcs()
        self._build_tangent_plane()


    def _build_wcs(self):
        """Create base WCS definition."""
        crval1 = self.band_index * 360. / self.band['NBAND']
        crval2 = self.band['DEC']
        cd = np.array([[-self.scale / 3600., 0], [0, self.scale / 3600.]], dtype=np.float64)
        self.wcs = astropy.wcs.WCS(naxis=2)
        self.wcs.crpix = [1., 1.]  # 1-based for FITS-compatibility
        self.wcs.crval = [crval1, crval2]
        self.wcs.cd = cd
        self.wcs.ctype = ['RA---TAN', 'DEC--TAN']

    def _build_tangent_plane(self):
        """Create tangent plane defintion with pixels"""
        # define how many pixels across this cell will have based on `self.size`
        # shift CRPIX in WCS to center of projection cell
        pass

    def find_sky_cells(self, footprint):
        """Return the sky cell indices from this projection cell that overlap the input footprint"""
        return 1, 1

class SkyCell(object):

    def __init__(self, x, y, projection_cell):
        self.x_index = x
        self.y_index = y
        self.sky_cell_id = "skycell_{}_x{}y{}".format(projection_cell.cell_id,
                                                      str(x).zfill(3), str(y).zfill(3))
        self._build_wcs(projection_cell.wcs)

    def rescale(self, scale):
        """Return WCS which has a user-defined scale."""
        pass
