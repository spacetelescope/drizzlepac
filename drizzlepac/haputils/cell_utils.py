import os
import shutil

from matplotlib import pyplot as plt
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

# SkyCell format: "skycell-p0000x000y000"
SKYCELL_NAME_FMT = f"skycell-p{{:{str(PCELL_STRLEN).zfill(2)}d}}x{{:02d}}y{{:02d}}"
SKYCELL_NXY = 50
SKYCELL_OVERLAP = 256


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

def interpret_scells(sky_cells):
    """Return dict of filenames each with the skycell name they overlap

    Parameters
    ----------
    sky_cells : dict
        Dictionary of sky-cell objects from `get_sky_cells`

    Returns
    -------
    sky_cell_files : dict
        Dictionary of ALL sky-cell IDs as a ';'-delimited string for each
        exposure(sky cell member), with exposure filenames as keys.

    """
    scell_files = {}
    for scell in sky_cells.values():
        for member in scell.members:
            if member not in scell_files:
                scell_files[member] = {}
            scell_files[member][scell.sky_cell_id] = scell

    # convert each entry into a ';'-delimited string instead of a list of IDs
    for member in scell_files:
        scell_files[member]['id'] = ';'.join([id for id in scell_files[member]])

    return scell_files



class SkyFootprint(object):

    def __init__(self, meta_wcs):

        self.meta_wcs = meta_wcs

        # the exp_masks dict records the individual footprints of each exposure
        self.exp_masks = {}
        self.members = []

        self.total_mask = np.zeros(meta_wcs.array_shape, dtype=np.int16)
        self.scaled_mask = np.zeros(meta_wcs.array_shape, dtype=np.float32)
        self.footprint = None

        self.edges = None
        self.edges_ra = None
        self.edges_dec = None
        self.polygon = None

    def build(self, expnames, scale=False, scale_kw='EXPTIME'):
        """ Create mask showing where all input exposures overlap the footprint's WCS

        Notes
        -----
        This method populates the following attributes (initialized as all zeros):
          - total_mask : shows number of chips per pixel
          - scaled_mask : if computed, shows (by default) exposure time per pixel

        Parameters
        -----------
        expnames : list
            List of filenames for all input exposures that overlap the SkyFootprint WCS

        scale : bool, optional
            If specified, scale each chip by the value of the `scale_kw` keyword from the input exposure.

        scale_kw : str, optional
            If `scale` is `True`, get the scaling value from this keyword.  This keyword is assumed to be
            in the PRIMARY header.

        """
        for exposure in expnames:
            if exposure not in self.members:
                self.members.append(exposure)
            self.exp_masks[exposure] = np.zeros(self.meta_wcs.array_shape, dtype=np.int16)
            exp = fits.open(exposure)
            if scale:
                scale_val = fits.getval(exposure, scale_kw)

            sci_extns = wcs_functions.get_extns(exp)
            for sci in sci_extns:
                wcs = HSTWCS(exp, ext=sci)
                edges_x = [0] * wcs.naxis2 + [wcs.naxis1 - 1] * wcs.naxis2 + list(range(wcs.naxis1)) * 2
                edges_y = list(range(wcs.naxis2)) * 2 + [0] * wcs.naxis1 + [wcs.naxis2 - 1] * wcs.naxis1

                sky_edges = wcs.pixel_to_world_values(edges_x, edges_y)
                meta_x, meta_y = self.meta_wcs.world_to_pixel_values(sky_edges[0], sky_edges[1])
                meta_x = meta_x.astype(np.int32)
                meta_y = meta_y.astype(np.int32)

                # Account for rounding problems with creating meta_wcs
                meta_y = np.clip(meta_y, 0, self.meta_wcs.array_shape[0] - 1)
                meta_x = np.clip(meta_x, 0, self.meta_wcs.array_shape[1] - 1)

                # apply meta_edges to blank mask
                # Use PIL to create mask
                # parray = np.array(meta_edges.T)
                parray = (meta_x, meta_y)
                polygon = list(zip(parray[0], parray[1]))
                nx = self.meta_wcs.array_shape[1]
                ny = self.meta_wcs.array_shape[0]
                img = Image.new('L', (nx, ny), 0)
                ImageDraw.Draw(img).polygon(polygon, outline=1, fill=1)
                blank = np.array(img).astype(np.int16)

                if scale:
                    scaled_blank = blank * scale_val
                    self.scaled_mask += scaled_blank
                self.exp_masks[exposure] += blank

            self.total_mask += self.exp_masks[exposure]


    # Methods with 'find' compute values
    # Methods with 'get' return values
    def find_footprint(self, member='total'):
        if self.total_mask is None:
            print("Please add exposures before computing footprint...")
        if member == 'total':
            mask = self.total_mask
        else:
            if member not in self.exp_masks:
                raise ValueError("Member {} not added to footprint".format(member))
            mask = self.exp_masks[member]
        self.footprint = np.clip(mask, 0, 1)
        # self.footprint = np.clip(self.total_mask, 0, 1)

    def find_edges(self, member='total'):
        self.find_footprint(member=member)

        edges = morphology.binary_erosion(self.footprint).astype(np.int16)
        self.edges = self.footprint - edges

    def find_corners(self, member='total'):
        self.find_footprint(member=member)

    def get_edges_sky(self, member='total'):
        self.find_footprint(member=member)
        if self.edges is None:
            self.find_edges()
        edges = np.where(self.edges)
        self.edges_ra, self.edges_dec = self.meta_wcs.pixel_to_world_values(edges[0], edges[1])
        # close the polygon
        self.edges_ra = np.append(self.edges_ra, [self.edges_ra[0]], axis=0)
        self.edges_dec = np.append(self.edges_dec, [self.edges_dec[0]], axis=0)
        return self.edges_ra, self.edges_dec

    def build_polygon(self, member='total'):
        if self.edges_ra is None:
            self.get_edges_sky(member=member)

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

    def get_footprint_hdu(self, filename=None, overwrite=True, mmeber='total'):
        self.find_footprint(member=member)
        return self._get_fits_hdu(self.footprint, filename=filename, overwrite=overwrite)

    def get_edges_hdu(self, filename=None, overwrite=True, member='total'):
        self.find_edges(member=member)
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

    def get_projection_cells(self, skyfootprint=None, member='total',
                             ra=None, dec=None, id=None):
        # Interpret footprint to get range of declination in mask
        if id is None:
            if ra is None:
                ra, dec = skyfootprint.get_edges_sky(member=member)
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

    def get_sky_cells(self, skyfootprint, member='total'):

        self.get_projection_cells(skyfootprint)

        # Find sky cells from identified projection cell(s) that overlap footprint
        sky_cells = {}
        for pcell in self.projection_cells:
            sky_cells.update(pcell.find_sky_cells(skyfootprint,
                                                 nxy=self.sc_nxy,
                                                 overlap=self.sc_overlap))

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
        self.wcs.pscale = self.scale
        self.wcs.orientat = 0.0
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

        member = 'total'

        # Get the edges of the mosaic on the sky
        mosaic_ra, mosaic_dec = mosaic.get_edges_sky(member=member)
        # Convert edges to positions in projection cell
        mosaic_edges_x, mosaic_edges_y = self.wcs.world_to_pixel_values(mosaic_ra, mosaic_dec)

        # Determine roughly what sky cells overlap this mosaic
        mosaic_edges_x = (mosaic_edges_x / skycell00.wcs.pixel_shape[0] + 0.5).astype(np.int32)
        mosaic_edges_y = (mosaic_edges_y / skycell00.wcs.pixel_shape[1] + 0.5).astype(np.int32)
        mosaic_xr = [mosaic_edges_x.min() - 1, mosaic_edges_x.max() + 1]
        mosaic_yr = [mosaic_edges_y.min() - 1, mosaic_edges_y.max() + 1]

        print("SkyCell Ranges: {}, {}".format(mosaic_xr, mosaic_yr))
        # for each suspected sky cell or neighbor, look for any pixel by pixel
        #    overlap with input mosaic footprint
        for xi in range(mosaic_xr[0], mosaic_xr[1]):
            for yi in range(mosaic_yr[0], mosaic_yr[1]):
                skycell = SkyCell(x=xi, y=yi, projection_cell=self)
                skycell.build_mask()

                sc_overlap = self.compute_overlap(skycell, mosaic_ra, mosaic_dec)

                print("    Checking SkyCell {},{} for overlap: {}".format(xi, yi, sc_overlap))
                if sc_overlap:
                    # Within this SkyCell, determine which members of the
                    # mosaic overlap with this SkyCell.
                    for filename in mosaic.members:
                        member_ra, member_dec = mosaic.get_edges_sky(member=filename)
                        member_overlap = self.compute_overlap(skycell,
                                                              member_ra,
                                                              member_dec)
                        if member_overlap:
                            skycell.members.append(filename)
                    # We found overlapping pixels from mosaic, so return this SkyCell
                    print("   Found overlap in SkyCell {}".format(skycell.sky_cell_id))
                    skycells[skycell.sky_cell_id] = skycell

        return skycells

    def compute_overlap(self, skycell, mosaic_ra, mosaic_dec):
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

        return sc_overlap

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

        self.members = []
        self.overlap = self.projection_cell.sc_overlap  # overlap between sky cells
        self.nxy = self.projection_cell.sc_nxy

        self._build_wcs()

    def _from_name(self, name):
        # parse name into projection cell and sky cell designations
        sc_names = name.split('-')
        scell_id = sc_names[1]
        pcell_id = int(scell_id[1:5])

        self.x_index = int(scell_id[6:8])
        self.y_index = int(scell_id[9:11])
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
        self.wcs.orientat = self.projection_cell.wcs.orientat
        self.wcs.pscale = self.projection_cell.wcs.pscale

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


#
#
# CKmeans implementation from github/llimllib/ckmeans
#
#
def ssq(j, i, sum_x, sum_x_sq):
    if (j > 0):
        muji = (sum_x[i] - sum_x[j-1]) / (i - j + 1)
        sji = sum_x_sq[i] - sum_x_sq[j-1] - (i - j + 1) * muji ** 2
    else:
        sji = sum_x_sq[i] - sum_x[i] ** 2 / (i+1)

    return 0 if sji < 0 else sji

def fill_row_k(imin, imax, k, S, J, sum_x, sum_x_sq, N):
    if imin > imax: return

    i = (imin+imax) // 2
    S[k][i] = S[k-1][i-1]
    J[k][i] = i

    jlow = k

    if imin > k:
        jlow = int(max(jlow, J[k][imin-1]))
    jlow = int(max(jlow, J[k-1][i]))

    jhigh = i-1
    if imax < N-1:
        jhigh = int(min(jhigh, J[k][imax+1]))

    for j in range(jhigh, jlow-1, -1):
        sji = ssq(j, i, sum_x, sum_x_sq)

        if sji + S[k-1][jlow-1] >= S[k][i]: break

        # Examine the lower bound of the cluster border
        # compute s(jlow, i)
        sjlowi = ssq(jlow, i, sum_x, sum_x_sq)

        SSQ_jlow = sjlowi + S[k-1][jlow-1]

        if SSQ_jlow < S[k][i]:
            S[k][i] = SSQ_jlow
            J[k][i] = jlow

        jlow += 1

        SSQ_j = sji + S[k-1][j-1]
        if SSQ_j < S[k][i]:
            S[k][i] = SSQ_j
            J[k][i] = j

    fill_row_k(imin, i-1, k, S, J, sum_x, sum_x_sq, N)
    fill_row_k(i+1, imax, k, S, J, sum_x, sum_x_sq, N)

def fill_dp_matrix(data, S, J, K, N):
    sum_x = np.zeros(N, dtype=np.float_)
    sum_x_sq = np.zeros(N, dtype=np.float_)

    # median. used to shift the values of x to improve numerical stability
    shift = data[N//2]

    for i in range(N):
        if i == 0:
            sum_x[0] = data[0] - shift
            sum_x_sq[0] = (data[0] - shift) ** 2
        else:
            sum_x[i] = sum_x[i-1] + data[i] - shift
            sum_x_sq[i] = sum_x_sq[i-1] + (data[i] - shift) ** 2

        S[0][i] = ssq(0, i, sum_x, sum_x_sq)
        J[0][i] = 0

    for k in range(1, K):
        if (k < K-1):
            imin = max(1, k)
        else:
            imin = N-1

        fill_row_k(imin, N-1, k, S, J, sum_x, sum_x_sq, N)

def ckmeans(data, n_clusters):
    if n_clusters <= 0:
        raise ValueError("Cannot classify into 0 or less clusters")
    if n_clusters > len(data):
        raise ValueError("Cannot generate more classes than there are data values")

    # if there's only one value, return it; there's no sensible way to split
    # it. This means that len(ckmeans([data], 2)) may not == 2. Is that OK?
    unique = len(set(data))
    if unique == 1:
        return [data]

    data.sort()
    n = len(data)

    S = np.zeros((n_clusters, n), dtype=np.float_)

    J = np.zeros((n_clusters, n), dtype=np.uint64)

    fill_dp_matrix(data, S, J, n_clusters, n)

    clusters = []
    cluster_right = n-1

    for cluster in range(n_clusters-1, -1, -1):
        cluster_left = int(J[cluster][cluster_right])
        clusters.append(data[cluster_left:cluster_right+1])

        if cluster > 0:
            cluster_right = cluster_left - 1

    return list(reversed(clusters))

##
## HELPER CODE FOR TESTS
##

# partition recipe modified from
# http://wordaligned.org/articles/partitioning-with-python
from itertools import chain, combinations

def sliceable(xs):
    '''Return a sliceable version of the iterable xs.'''
    try:
        xs[:0]
        return xs
    except TypeError:
        return tuple(xs)

def partition_n(iterable, n):
    s = sliceable(iterable)
    l = len(s)
    b, mid, e = [0], list(range(1, l)), [l]
    splits = (d for i in range(l) for d in combinations(mid, n-1))
    return [[s[sl] for sl in map(slice, chain(b, d), chain(d, e))]
            for d in splits]

def squared_distance(part):
    mean = sum(part)/len(part)
    return sum((x-mean)**2 for x in part)

# given a partition, return the sum of the squared distances of each part
def sum_of_squared_distances(partition):
    return sum(squared_distance(part) for part in partition)

# brute force the correct answer by testing every partition.
def min_squared_distance(data, n):
    return min((sum_of_squared_distances(partition), partition)
                for partition in partition_n(data, n))


# if __name__ == "__main__":
def ckmeans_test():
    try:
        ckmeans([], 10)
        1/0
    except ValueError:
        pass

    tests = [
        (([1], 1),                    [[1]]),
        (([0,3,4], 2),                [[0], [3,4]]),
        (([-3,0,4], 2),               [[-3,0], [4]]),
        (([1,1,1,1], 1),              [[1,1,1,1]]),
        (([1,2,3], 3),                [[1], [2], [3]]),
        (([1,2,2,3], 3),              [[1], [2,2], [3]]),
        (([1,2,2,3,3], 3),            [[1], [2,2], [3,3]]),
        (([1,2,3,2,3], 3),            [[1], [2,2], [3,3]]),
        (([3,2,3,2,1], 3),            [[1], [2,2], [3,3]]),
        (([3,2,3,5,2,1], 3),          [[1,2,2], [3,3], [5]]),
        (([0,1,2,100,101,103], 2),    [[0,1,2], [100,101,103]]),
        (([0,1,2,50,100,101,103], 3), [[0,1,2], [50], [100,101,103]]),
        (([-1,2,-1,2,4,5,6,-1,2,-1], 3),
            [[-1, -1, -1, -1], [2, 2, 2], [4, 5, 6]]),
    ]

    for test in tests:
        args, expected = test
        try:
            result = ckmeans(*args)
        except Exception:
            print("✗ {}, {}".format(args[0], args[1], result))
            raise
        errormsg = "✗ ckmeans({}) = {} != {}\n{} > {}".format(
                args, result, expected,
                sum_of_squared_distances(result),
                sum_of_squared_distances(expected))
        assert np.array_equal(result, expected), errormsg
        print("✓ {}".format(result))
