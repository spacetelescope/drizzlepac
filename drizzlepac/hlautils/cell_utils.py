import os
import shutil
import sys

from matplotlib import pyplot as plt
from scipy import ndimage
from scipy.ndimage import morphology
import numpy as np
import astropy
from astropy import units as u
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from spherical_geometry.polygon import SphericalPolygon

from PIL import Image, ImageDraw

from stwcs.wcsutil import HSTWCS

from .. import wcs_functions


# Default grid definition file
_fpath = os.path.abspath(os.path.dirname(__file__))
PCELL_PATH = os.path.join(os.path.dirname(_fpath), 'pars')
PCELL_FILENAME = 'allsky_cells.fits'
PCELL_STRLEN = 4

# SkyCell format: "skycell_p0000_x000y000"
SKYCELL_NAME_FMT = f"skycell_p{{:{str(PCELL_STRLEN).zfill(2)}d}}_x{{:03d}}y{{:03d}}"
SKYCELL_NXY = 50
SKYCELL_OVERLAP = 256

NDIMAGE_STRUCT2 = ndimage.generate_binary_structure(2, 2)


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
        expnames = visit_input.copy()
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

    # Check that all exposures have up-to-date WCS solutions
    #  This will weed out exposures which were not processed by the pipeline
    #  such as those with EXPTIME==0
    for filename in expnames:
        with fits.open(filename) as fimg:
            print("Checking {}".format(filename))
            if 'wcsname' not in fimg[1].header:
                expnames.remove(filename)
    if len(expnames) == 0:
        print("No valid exposures to define sky cells")
        return None

    # Initialize all sky tessellation object definitions
    # This includes setting the pixel scale.
    sky_grid = GridDefs(scale=scale, cell_size=cell_size)

    # build reference wcs for combined footprint of all input exposures
    meta_wcs = wcs_functions.make_mosaic_wcs(expnames, scale=sky_grid.scale)

    # create footprint on the sky (as a tangent plane array) for all input exposures using meta_wcs
    footprint = SkyFootprint(meta_wcs)
    footprint.build(expnames)

    # Use this footprint to identify overlapping sky cells
    sky_cells = sky_grid.get_sky_cells(footprint)

    return sky_cells


class SkyFootprint(object):

    def __init__(self, meta_wcs):

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
                edges_x = [0]*wcs.naxis2 + [wcs.naxis1-1]*wcs.naxis2 + list(range(wcs.naxis1)) * 2
                edges_y = list(range(wcs.naxis2)) * 2 + [0]*wcs.naxis1 + [wcs.naxis2-1]*wcs.naxis1

                sky_edges = wcs.pixel_to_world_values(np.vstack([edges_x, edges_y]).T)
                meta_edges = self.meta_wcs.world_to_pixel_values(sky_edges).astype(np.int32)
                # Account for rounding problems with creating meta_wcs
                meta_edges[:,1] = np.clip(meta_edges[:,1], 0, self.meta_wcs.array_shape[0]-1)
                meta_edges[:,0] = np.clip(meta_edges[:,0], 0, self.meta_wcs.array_shape[1]-1)

                # apply meta_edges to blank mask
                # Use PIL to create mask
                parray = np.array(meta_edges.T)
                polygon = list(zip(parray[0], parray[1]))
                nx = self.meta_wcs.array_shape[1]
                ny = self.meta_wcs.array_shape[0]
                img = Image.new('L', (nx, ny) , 0)
                ImageDraw.Draw(img).polygon(polygon, outline=1, fill=1)
                blank = np.array(img)

            self.total_mask += blank.astype(np.int16)

    # Methods with 'find' compute values
    # Methods with 'get' return values
    def find_footprint(self):
        if self.total_mask is None:
            print("Please add exposures before computing footprint...")
        self.footprint = np.clip(self.total_mask, 0, 1)

    def find_edges(self):
        if self.footprint is None:
            self.find_footprint()
        edges = morphology.binary_erosion(self.footprint).astype(np.int16)
        self.edges = self.footprint - edges

    def find_corners(self):
        if self.footprint is None:
            self.find_footprint()

    def get_edges_sky(self):
        if self.edges is None:
            self.find_edges()
        edges = np.where(self.edges)
        self.edges_ra, self.edges_dec = self.meta_wcs.pixel_to_world_values(edges[0], edges[1])
        # close the polygon
        self.edges_ra = np.append(self.edges_ra, [self.edges_ra[0]], axis=0)
        self.edges_dec = np.append(self.edges_dec, [self.edges_dec[0]], axis=0)
        return self.edges_ra, self.edges_dec

    def build_polygon(self):
        if self.edges_ra is None:
            self.get_edges_sky()

        self.polygon = SphericalPolygon.from_radec(self.edges_ra,
                                                   self.edges_dec,
                                                   self.meta_wcs.wcs.crval)

    def _get_fits_hdu(self, data, filename=None, overwrite=True):
        hdulist = fits.HDUList()
        phdu = fits.PrimaryHDU(data=data, header=self.meta_wcs.to_header())
        hdulist.append(phdu)
        if filename:
            hdulist.writeto(filename, overwrite=overwrite)
        return hdulist

    def get_footprint_hdu(self, filename=None, overwrite=True):
        if self.footprint is None:
            self.find_footprint()
        return self._get_fits_hdu(self.footprint, filename=filename, overwrite=overwrite)

    def get_edges_hdu(self, filename=None, overwrite=True):
        if self.edges is None:
            self.find_edges()
        return self._get_fits_hdu(self.edges, filename=filename, overwrite=overwrite)

    def get_mask_hdu(self, filename=None, overwrite=True):
        return self._get_fits_hdu(self.total_mask, filename=filename, overwrite=overwrite)


class GridDefs(object):

    def __init__(self, scale=None, cell_size=None):
        """Setup tesselation based on installed definition of grid."""
        fname = os.path.join(PCELL_PATH, PCELL_FILENAME)

        self.hdu = fits.open(fname)
        self.rings = self.hdu[1].data

        # Extract projection cell defaults
        self.scale = scale if scale else self.hdu[0].header['PC_SCALE']
        self.cell_size = cell_size if cell_size else self.hdu[0].header['PC_SIZE']
        # Extract sky cell defaults
        self.sc_overlap = self.hdu[0].header['SC_OLAP']
        self.sc_nxy = self.hdu[0].header['SC_NXY']

    def find_ring_by_id(self, id):
        return self.rings[np.searchsorted(self.rings['projcell'], id) - 1]

    def get_projection_cells(self, skyfootprint=None, ra=None, dec=None, id=None):
        # Interpret footprint to get range of declination in mask
        if id is None:
            if ra is None:
                ra, dec = skyfootprint.get_edges_sky()
            # Find band[s] that overlap footprint
            self._find_bands(dec)

            self.projection_cells = []
            # Define numerical position in band for projection cell
            # self.band_index = self.projection_cell_id - self.band['PROJCELL']
            for band in self.bands:
                # compute band_index, one for each projection cell that overlaps the footprint
                nra = ra % 360.0
                nband = band['NBAND']
                band_index = np.unique(np.rint(nra * nband / 360.0).astype(int) % nband)
                self.projection_cells += [ProjectionCell(index, band, self.scale) for index in band_index]
        else:
            self.projection_cells = [ProjectionCell(index=i, scale=self.scale) for i in id]

    def get_sky_cells(self, skyfootprint):

        self.get_projection_cells(skyfootprint)

        # Find sky cells from identified projection cell(s) that overlap footprint
        sky_cells = {}
        for pcell in self.projection_cells:
            sky_cells.update(pcell.find_sky_cells(skyfootprint, self.sc_nxy, self.sc_overlap))

        return sky_cells

    def _find_bands(self, dec):
        """ Select the band or bands which encompass the provided footprint.

        The footprint will be a ndarray mask with pixel value of 1 where the exposures are located,
        and with a fully defined WCS to convert pixel positions to sky coordinates.

        """
        if not isinstance(dec, list) and not isinstance(dec, np.ndarray):
            dec = [dec]

        # Select all bands that overlap
        # find dec zones where rings.dec_min <= dec <= rings.dec_max
        maxb = dec[0] <= self.rings.field('dec_max')
        minb = dec[0] >= self.rings.field('dec_min')
        for d in dec:
            maxb = np.bitwise_and(d <= self.rings.field('dec_max'), maxb)
            minb = np.bitwise_and(d >= self.rings.field('dec_min'), minb)

        band_indx = np.where(np.bitwise_and(maxb, minb))[0]
        bands = np.sort(np.unique(band_indx))

        # Record these values as attributes for use in other methods
        self.bands = [self.rings[b] for b in bands]

    def plot(self, projection='aitoff'):
        if not self.projection_cells:
            print("Please run `get_projection_cells()' first...")
            return
        plt.figure()
        plt.subplot(111, projection=projection)
        plt.grid(True)
        for pc in self.projection_cells:
            plt.fill(pc.footprint[:,0], pc.footprint[:,1],
                     facecolor='green', edgecolor='forestgreen',
                     alpha=0.25)
            plt.text(pc.footprint[0,0], pc.footprint[0,1], "{}".format(pc.cell_id),
                     horizontalalignment='right', verticalalignment='bottom')



class ProjectionCell(object):

    def __init__(self, index=None, band=None, scale=None,
                        nxy=None, overlap=None):
        """Build projection cell for cell with name `skycell_NNNNN`

        If `band` is not specified, it will open the grid definitions file to
        obtain the band definition for the cell with the specified `index`.

        """
        if band:
            self.band_index = index
            self.band = band
            self.cell_id = band['PROJCELL'] + index
            self.scale = scale
        else:
            self._from_index(index)

        # Record user-specified sky cell defaults
        if overlap:
            self.sc_overlap = overlap
        if nxy:
            self.sc_nxy = nxy

        # Generate WCS for projection cell
        self._build_wcs()

    def _from_index(self, id):
        grid_defs = GridDefs()
        self.band = grid_defs.find_ring_by_id(id)
        self.band_index = id - self.band['projcell']
        self.cell_id = id
        self.scale = grid_defs.scale
        self.sc_overlap = grid_defs.sc_overlap
        self.sc_nxy = grid_defs.sc_nxy

    def _build_wcs(self):
        """Create base WCS definition."""
        crval1 = self.band_index * 360. / self.band['NBAND']
        crval2 = self.band['DEC']
        cd = np.array([[-self.scale / 3600., 0], [0, self.scale / 3600.]], dtype=np.float64)

        self.wcs = astropy.wcs.WCS(naxis=2)
        self.wcs.wcs.crval = [crval1, crval2]
        self.wcs.wcs.cd = cd
        self.wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']

        # Read in overall size of cell in pixels from grid definitions file
        naxis1 = self.band['XCELL']
        naxis2 = self.band['YCELL']

        # apply new definition to cell WCS
        self.wcs.wcs.crpix = [naxis1 / 2. + 0.5, naxis2 / 2. + 0.5]
        self.wcs.pixel_shape = (naxis1, naxis2)
        self.footprint = self.wcs.calc_footprint()
        # close the polygon
        self.corners = np.append(self.footprint, [self.footprint[0]], axis=0)

    def build_mask(self):
        naxis1, naxis2 = self.wcs.pixel_shape
        edges_x = [0]*naxis2 + [naxis1-1]*naxis2 + list(range(naxis1)) * 2
        edges_y = list(range(naxis2)) * 2 + [0]*naxis1 + [naxis2-1]*naxis1

        polygon = list(zip(edges_x, edges_y))
        img = Image.new("L", (naxis1, naxis2), 0)
        ImageDraw.Draw(img).polygon(polygon, outline=1, fill=1)
        mask = np.array(img)

        self.mask = mask

    def build_polygon(self):

        inner_pix = self.wcs.pixel_to_world_values(2, 2)
        # define polygon on the sky
        self.polygon = SphericalPolygon.from_radec(self.corners[:, 0], self.corners[:, 1], inner_pix)


    def find_sky_cells(self, mosaic, nxy=None, overlap=None):
        """Return the sky cell indices from this projection cell that overlap the input footprint"""
        # Record values for sky cell definitions provided by user, if any.
        if nxy:
            self.sc_nxy = nxy
        if overlap:
            self.sc_overlap = overlap
        print("Looking for sky cells in PROJECTION_CELL {}".format(self.cell_id))

        skycell00 = SkyCell(x=0, y=0, projection_cell=self)
        skycells = {}

        # Get the edges of the mosaic on the sky
        mosaic_ra, mosaic_dec = mosaic.get_edges_sky()
        # Convert edges to positions in projection cell
        mosaic_edges = self.wcs.world_to_pixel_values(mosaic_ra, mosaic_dec)
        # Determine roughly what sky cells overlap this mosaic
        mosaic_edges[0] = (mosaic_edges[0] / skycell00.wcs.pixel_shape[0] + 0.5).astype(np.int32)
        mosaic_edges[1] = (mosaic_edges[1] / skycell00.wcs.pixel_shape[1] + 0.5).astype(np.int32)
        mosaic_xr = [mosaic_edges[0].min() - 1, mosaic_edges[0].max() + 1]
        mosaic_yr = [mosaic_edges[1].min() - 1, mosaic_edges[1].max() + 1]

        print("SkyCell Ranges: {}, {}".format(mosaic_xr, mosaic_yr))
        # for each suspected sky cell or neighbor, look for any pixel by pixel
        #    overlap with input mosaic footprint
        for xi in range(mosaic_xr[0], mosaic_xr[1]):
            for yi in range(mosaic_yr[0], mosaic_yr[1]):
                skycell = SkyCell(x=xi, y=yi, projection_cell=self)
                skycell.build_mask()

                # Translate mosaic edges into SkyCell WCS coordinate frame
                mosaic_xy = skycell.wcs.world_to_pixel_values(mosaic_ra, mosaic_dec)

                # Identify edge pixels which fall outside the sky cell
                #  by comparing to each direction (-X, +X, -Y, +Y) separately
                mosaic_offcell = mosaic_xy[0] < 0
                mosaic_offcell = np.bitwise_or(mosaic_offcell,
                                                mosaic_xy[0] > skycell.wcs.pixel_shape[0])
                mosaic_offcell = np.bitwise_or(mosaic_offcell, mosaic_xy[1] < 0)
                mosaic_offcell = np.bitwise_or(mosaic_offcell,
                                                mosaic_xy[1] > skycell.wcs.pixel_shape[1])

                # With all out of bounds pixels masked out, see if any are left
                sc_overlap = np.any(~mosaic_offcell)

                print("    Checking SkyCell {},{} for overlap: {}".format(xi, yi, sc_overlap))
                if sc_overlap:
                    # We found overlapping pixels from mosaic, so return this SkyCell
                    print("   Found overlap in SkyCell {}".format(skycell.sky_cell_id))
                    skycells[skycell.sky_cell_id] = skycell

        return skycells

    def plot(self, output=None, color='b'):
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection="mollweide")
        ax.fill(self.corners[:, 0], self.corners[:, 1],
                facecolor='green',edgecolor='forestgreen', alpha=0.25)
        if output:
            fig.write(output)
        return ax


class SkyCell(object):

    def __init__(self, name=None, projection_cell=None, x=None, y=None):
        """Define sky cell at position x,y within projection cell.
        X,Y positions need to be 1-based.
        """
        if name:
            self._from_name(name)
        else:
            self.x_index = x
            self.y_index = y
            self.sky_cell_id = SKYCELL_NAME_FMT.format(projection_cell.cell_id, x, y)
            self.projection_cell = projection_cell

        self.overlap = self.projection_cell.sc_overlap  # overlap between sky cells
        self.nxy = self.projection_cell.sc_nxy

        self._build_wcs()

    def _from_name(self, name):
        # parse name into projection cell and sky cell designations
        sc_names = name.split('_')
        pcell_id = int(sc_names[1][1:])

        self.x_index = int(sc_names[2][1:4])
        self.y_index = int(sc_names[2][5:8])
        self.projection_cell = ProjectionCell(index=pcell_id)
        self.sky_cell_id = name

    def __repr__(self):
        return "SkyCell object: {}".format(self.sky_cell_id)

    def rescale(self, scale):
        """Return WCS which has a user-defined scale."""
        pass

    def _build_wcs(self):
        pc_nx = self.projection_cell.wcs.pixel_shape[0]
        pc_ny = self.projection_cell.wcs.pixel_shape[1]
        naxis1 = int((pc_nx - 2 * self.overlap) / self.nxy + 0.5)
        naxis2 = int((pc_ny - 2 * self.overlap) / self.nxy + 0.5)
        crpix1 = self.projection_cell.wcs.wcs.crpix[0] - (((self.x_index) * naxis1) - self.overlap)
        crpix2 = self.projection_cell.wcs.wcs.crpix[1] - (((self.y_index) * naxis2) - self.overlap)

        # apply definitions
        self.wcs = astropy.wcs.WCS(naxis=2)
        self.wcs.wcs.crpix = [crpix1, crpix2]
        self.wcs.wcs.crval = self.projection_cell.wcs.wcs.crval
        self.wcs.wcs.cd = self.projection_cell.wcs.wcs.cd
        self.wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
        self.wcs.pixel_shape = (naxis1, naxis2)
        self.wcs.ltv1 = self.projection_cell.wcs.wcs.crpix[0] - crpix1
        self.wcs.ltv2 = self.projection_cell.wcs.wcs.crpix[1] - crpix2

        self.corners = self.wcs.calc_footprint()
        # close the polygon
        self.corners = np.append(self.corners, [self.corners[0]], axis=0)

    def build_polygon(self):
        inner_pix = self.wcs.pixel_to_world_values(2, 2)
        # define polygon on the sky
        self.polygon = SphericalPolygon.from_radec(self.corners[:, 0],
                                                   self.corners[:, 1],
                                                   inner_pix)
    def build_mask(self):
        naxis1, naxis2 = self.wcs.pixel_shape
        edges_x = [0]*naxis2 + [naxis1-1]*naxis2 + list(range(naxis1)) * 2
        edges_y = list(range(naxis2)) * 2 + [0]*naxis1 + [naxis2-1]*naxis1

        polygon = list(zip(edges_x, edges_y))
        img = Image.new("L", (naxis1, naxis2), 0)
        ImageDraw.Draw(img).polygon(polygon, outline=1, fill=1)
        mask = np.array(img)

        self.mask = mask

#
# Utility functions used in generating or supporting the grid definitions
#
def update_grid_defs(pc_size=5.0, output=None, grid_file=None):
    """Computes updated values for bands and projection cells.

    Parameters
    -----------
    pc_size : float
        Size of each side of the projection cell or width of each band on
        the sky in degrees.  If `None`, the default value will be read in
        from the `PC_SIZE` keyword from the PRIMARY header of the default
        grid definitions file.

    output : str, optional
        Name of output grid definition file.  If `None`, it will write out
        the updated table with the original filename in the current directory
        overwriting any previous file.

    """
    # read in default grid definition file
    if not grid_file:
        grid_file = os.path.join(PCELL_PATH, PCELL_FILENAME)
    grid_defs = fits.open(grid_file)
    grid = grid_defs[1].data
    pc_scale = grid_defs[0].header['PC_SCALE']

    if not pc_size:
        pc_size = grid_defs[0].header['PC_SIZE']

    pos_angle = [0.0 * u.deg, 90.0 * u.deg, 180.0 * u.deg, 270.0 * u.deg]
    pc_edge = pc_size / 2.0 * u.deg

    # Compute size on the sky of the first cell in each band
    # Compute edges using:
    for nband in range(len(grid)):
        c1 = SkyCoord(ra=0. * u.deg, dec=grid[nband]['DEC'] * u.deg, frame='icrs')
        pcell = ProjectionCell(index=0, band=grid[nband], scale=pc_scale)

        # Compute offset to center of each edge, +/- RA and +/- Dec
        c1_edges = [c1.directional_offset_by(p, pc_edge) for p in pos_angle]

        # Convert offset to edge center into distance in pixels from center of cell
        c1_pixels = [np.abs(pcell.wcs.world_to_pixel_values(e.ra, e.dec)) for e in c1_edges]

        # Compute overall size of cell in pixels
        naxis1, naxis2 = (np.array(c1_pixels).sum(axis=0) + 1).astype(np.int)
        # apply new definition to cell WCS
        pcell.wcs.wcs.crpix = [naxis1 / 2. + 0.5, naxis2 / 2. + 0.5]
        pcell.wcs.naxis1 = naxis1
        pcell.wcs.naxis2 = naxis2

        # Determine extent of band
        min_dec, max_dec = compute_band_height(pcell.wcs)
        # Account for wrapping over each pole
        if nband == 0:
            min_dec = grid[nband]['DEC']
        if nband == len(grid) - 1:
            max_dec = grid[nband]['DEC']

        # Update definition for this band in table with newly computed values
        grid[nband]['DEC_MIN'] = min_dec
        grid[nband]['DEC_MAX'] = max_dec
        grid[nband]['XCELL'] = naxis1  # supposed to be sky cell size
        grid[nband]['YCELL'] = naxis2  # supposed to be sky cell size

    # write out updated grid definitions file
    if not output:
        output = PCELL_FILENAME

    # write to path included in 'output', defaulting to current working dir
    grid_defs.writeto(output, overwrite=True)

def compute_band_height(wcs):
    """Compute size in pixels of tangent plane"""
    edges = []
    edges += [[0, i] for i in range(wcs.naxis2)]
    edges += [[wcs.naxis1, i] for i in range(wcs.naxis2)]
    edges += [[i, 0] for i in range(wcs.naxis1)]
    edges += [[i, wcs.naxis2] for i in range(wcs.naxis1)]

    edge_sky = wcs.pixel_to_world_values(edges)
    return min(edge_sky[:, 1]), max(edge_sky[:, 1])
