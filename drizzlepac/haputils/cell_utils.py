import os
import shutil
from itertools import chain, combinations

from matplotlib import path
from matplotlib import pyplot as plt
from skimage.feature import corner_peaks, corner_harris
from scipy import ndimage
from scipy.spatial import distance
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

SUPPORTED_SCALES = {'fine': 0.04, 'coarse': 0.12}  # arcseconds/pixel


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
    meta_wcs = wcs_functions.make_mosaic_wcs(expnames, rot=0.0, scale=sky_grid.scale)

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
    """Object for computing footprint of overlapping images.

    This class allows a user to virtually build up a mosaic from a set of
    presumably overlapping exposures.

    Attributes
    -----------
    meta_wcs : `stwcs.wcsutil.HSTWCS`
        WCS of entire mosaic
    members : `list`
        List of filenames for all input exposures that make up the mosaic
    total_mask : `numpy.ndarray`
        Mask of combined footprint of all input exposures scaled by number of exposures at each pixel
    scaled_mask : `numpy.ndarray`
        Mask of combined footprint scaled by each input's exposure value (e.g., EXPTIME)
    footprint : `numpy.ndarray`
        Binary (numpy.int16) mask of combined footprint of all input exposures
    corners : `numpy.ndarray`
        List of vertices of the footprint in sky coordinates
    xy_corners : `numpy.ndarray`
        List of vertices of the footprint in X,Y coordinates
    edge_pixels : `numpy.ndarray`
        List of all pixels along the outside edge of the footprint, in
        counter-clockwise order
    exp_masks : `dict`
        Separate entries for each exposure containing:

        ``"mask"``
            Binary mask showing where this exposure overlaps the mosaic
        ``"xy_corners"``
            positions of the corners of the exposure in mosaic X,Y coordinates
        ``"sky_corners"``
            positions of the corners of the exposure in sky coordinates
            derived using the input exposure's `wcs.calc_footprint()` method
    edges : `numpy.ndarray`
        Binary (numpy.int16) mask containing only the pixels along the outside edge
         of the total footprint

    Methods
    -------
    build(expnames, scale=False, scale_kw='EXPTIME')
        Primary method for building up the mask of the footprint based on the
        provided input exposures.  This optionally allows the user to also compute an
        exposure time mask if desired.
    find_footprint()
        Once all input exposures have been used to build the masks, this method
        generates the total footprint masks.
    find_corners()
        This method computes the vertices of the total mask, something useful for
        creating the region FITS keywords (like S_REGION).
    """
    def __init__(self, meta_wcs):

        self.meta_wcs = meta_wcs

        # the exp_masks dict records the individual footprints of each exposure
        self.exp_masks = {}
        self.members = []
        self.corners = []
        self.sky_corners = []
        self.xy_corners = []
        self.edge_pixels = []

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
            self.exp_masks[exposure] = {'mask': None, 'sky_corners': [], 'xy_corners': []}
            self.exp_masks[exposure]['mask'] = np.zeros(self.meta_wcs.array_shape, dtype=np.int16)
            exp = fits.open(exposure)
            if scale:
                scale_val = fits.getval(exposure, scale_kw)

            sci_extns = wcs_functions.get_extns(exp)
            for sci in sci_extns:
                wcs = HSTWCS(exp, ext=sci)
                # save the footprint for each chip as RA/Dec corner positions
                radec = wcs.calc_footprint().tolist()
                self.exp_masks[exposure]['sky_corners'].append(radec)

                # Also save those corner positions as X,Y positions in the footprint
                xycorners = self.meta_wcs.all_world2pix(radec, 0).astype(np.int32).tolist()
                self.exp_masks[exposure]['xy_corners'].append(xycorners)

                # Now compute RA/Dec of all pixels along each edge
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
                self.exp_masks[exposure]['mask'] += blank

            self.total_mask += self.exp_masks[exposure]['mask']



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
            mask = self.exp_masks[member]['mask']
        self.footprint = np.clip(mask, 0, 1)


    def find_edges(self, member='total'):

        if self.footprint is None:
            self.find_footprint(member=member)

        edges = ndimage.binary_erosion(self.footprint).astype(np.int16)
        self.edges = self.footprint - edges


    def find_corners(self, member='total', sensitivity=0.1):
        """Extract corners from footprint """
        if len(self.members) == 0:
            print("Please add exposures before looking for corners...")
            return

        # Insure footprint has been determined
        if self.footprint is None:
            self.find_footprint(member=member)

        if member == 'total':
            # Use Harris corner detection to identify corners in the
            # total footprint
            mask_corners = corner_peaks(corner_harris(self.footprint),
                                   min_distance=1,
                                   threshold_rel=0.5)
            xy_corners = mask_corners * 0.
            xy_corners[:, 0] = mask_corners[:, 1]
            xy_corners[:, 1] = mask_corners[:, 0]
            self.full_corners = xy_corners.copy()

            # with this Nx2 array of points, we need to order them according
            # to IVOA standards: counter-clockwise when oriented N up starting from
            # northernmost point
            center = self.meta_wcs.wcs.crpix
            dist, phi, deg = cart2pol(xy_corners[:, 0] - center[0], xy_corners[:, 1] - center[1])
            # set '0' to be at -45deg.  This will avoid points at 358 deg being skipped over for a
            # point at 88 degrees, for example, insuring that more points near the top of the image
            # get selected as the starting point.
            deg += 45.0
            deg[deg > 360.] -= 360.
            radial_order = np.argsort(deg)

            # Create a mask from the total footprint consisting solely of the
            # pixels at the outer edge, ordered in clockwise fashion.
            inside = ndimage.binary_erosion(self.footprint).astype(np.int16)
            edge_lists = trace_polygon(inside - ndimage.binary_erosion(inside).astype(np.int16),
                                        sensitivity=sensitivity)
            xy_pixels = None
            for edge_pixels in edge_lists:
                # Start matching the xy_corner positions, one-by-one, to pixels
                # along the edge.  We only want the index of the corner that matches
                # as we travel along the edge in order which will be used to
                # sort the corner positions.
                #
                # start at the position closest to North where `deg` is closest to 0
                corner_dist = distance.cdist(edge_pixels, [xy_corners[radial_order[0]]])
                start_indx = np.where(corner_dist == corner_dist.min())[0][0]
                print(start_indx)
                if start_indx != 0:
                    # re-order edge_pixels so that the list starts at this pixel.
                    ordered_edge = edge_pixels * 0.
                    ordered_edge[:edge_pixels.shape[0] - start_indx] = edge_pixels[start_indx:]
                    ordered_edge[-start_indx:] = edge_pixels[:start_indx]
                else:
                    ordered_edge = edge_pixels

                # Now compute the distances for all the identified corners
                # ordered_dists = distance.cdist(xy_corners, ordered_edge)
                # This results in a list of distances for each xy_corner
                # Now sort by index along ordered edge
                # dist_indx = np.sort([np.where(dist == dist.min())[0][0] for dist in ordered_dists])

                # determine RA/Dec of these positions in the image
                # use this to make sure they are ordered correctly
                if xy_pixels is None:
                    xy_pixels = detect_corners(ordered_edge)
                    xy_pixels = np.concatenate([xy_pixels, [xy_pixels[0]]])  # close the polygon
                    sky_corners = self.meta_wcs.all_pix2world(xy_pixels, 0)
                    ordered_xy = [ordered_edge.copy()]
                else:
                    new_xy = detect_corners(ordered_edge)
                    new_xy = np.concatenate([new_xy, [new_xy[0]]])  # close the polygon
                    xy_pixels = np.concatenate([xy_pixels, new_xy])
                    sky_corners = np.concatenate([sky_corners, self.meta_wcs.all_pix2world(new_xy, 0)])
                    ordered_xy.append(ordered_edge)

            xy_corners = xy_pixels
            corners = sky_corners
            ordered_edge = ordered_xy

            # clean up memory a bit
            del edge_pixels
            del inside
            del corner_dist

        else:
            if member not in self.exp_masks:
                raise ValueError("Member {} not added to footprint".format(member))
            xy_corners = self.exp_masks[member]['xy_corners']
            corners = self.meta_wcs.all_pix2world(xy_corners, 0)

        self.edge_pixels = ordered_edge
        self.xy_corners = xy_corners
        self.corners = corners

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

    def get_footprint_hdu(self, filename=None, overwrite=True, member='total'):
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
            plt.fill(pc.footprint[:, 0], pc.footprint[:, 1],
                     facecolor='green', edgecolor='forestgreen',
                     alpha=0.25)
            plt.text(pc.footprint[0, 0], pc.footprint[0, 1], "{}".format(pc.cell_id),
                     horizontalalignment='right', verticalalignment='bottom')



class ProjectionCell(object):

    def __init__(self, index=None, band=None, scale=None,
                        nxy=None, overlap=None):
        """Build projection cell for cell with name `skycell_NNNNN`

        If `band` is not specified, it will open the grid definitions file to
        obtain the band definition for the cell with the specified `index`.

        """
        self.scale = scale
        if band:
            self.band_index = index
            self.band = band
            self.cell_id = band['PROJCELL'] + index
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
        self.scale = grid_defs.scale if self.scale is None else self.scale
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
        edges_x = [0] * naxis2 + [naxis1 - 1] * naxis2 + list(range(naxis1)) * 2
        edges_y = list(range(naxis2)) * 2 + [0] * naxis1 + [naxis2 - 1] * naxis1

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
                facecolor='green', edgecolor='forestgreen', alpha=0.25)
        if output:
            fig.write(output)
        return ax


class SkyCell(object):

    def __init__(self, projection_cell=None, x=None, y=None, scale="fine"):
        """Define sky cell at position x,y within projection cell.

        Parameters
        ===========
        name : str, optional
            Name of the sky cell in the format 'skycell-p1234x01y01'

        projection_cell : object, optional
            ProjectionCell instance which this SkyCell will be based upon.

        x : int, optional
            X position (1-based) of this SkyCell within the ProjectionCell

        y : int, optional
            Y position (1-based) of this SkyCell within the ProjectionCell

        scale : str or float, optional
            Plate scale to be used to define the SkyCell WCS.  The strings
            'fine' or 'coarse' can be used to refer to default plate scales.
            Default values are specified in the dict `cell_utils.SUPPORTED_SCALES`.
            Alternatively, floating-point values can be provided to specify
            the exact pixel size in arcseconds/pixel should be used.

        """
        # Interpret scale term, if provided
        self.scale = SUPPORTED_SCALES.get(scale, None) if isinstance(scale, str) else scale

        self.x_index = x
        self.y_index = y
        self.sky_cell_id = SKYCELL_NAME_FMT.format(projection_cell.cell_id, x, y)
        self.projection_cell = projection_cell

        self.members = []
        self.overlap = self.projection_cell.sc_overlap  # overlap between sky cells
        self.nxy = self.projection_cell.sc_nxy

        self._build_wcs()

    @classmethod
    def from_name(cls, name, scale="fine"):
        # parse name into projection cell and sky cell designations
        sc_names = name.split('-')
        scell_id = sc_names[1]
        pcell_id = int(scell_id[1:5])
        x = int(scell_id[6:8])
        y = int(scell_id[9:11])
        return cls(projection_cell=ProjectionCell(index=pcell_id), x=x, y=y, scale=scale)

    def __repr__(self):
        return "SkyCell object: {}".format(self.sky_cell_id)

    def rescale(self, scale):
        """Return WCS which has a user-defined scale."""
        pass

    def _build_wcs(self):
        # Determine plate scale ratio between sky cell layer and projection cell
        ratio = self.projection_cell.wcs.pscale / self.scale
        # Define attributes based on projection cell
        pc_nx = self.projection_cell.wcs.pixel_shape[0] * ratio
        pc_ny = self.projection_cell.wcs.pixel_shape[1] * ratio

        naxis1 = int((pc_nx - (2 * self.overlap)) / self.nxy + 0.5)
        naxis2 = int((pc_ny - (2 * self.overlap)) / self.nxy + 0.5)
        crpix1 = (self.projection_cell.wcs.wcs.crpix[0] * ratio) - ((self.x_index * naxis1) - self.overlap)
        crpix2 = (self.projection_cell.wcs.wcs.crpix[1] * ratio) - ((self.y_index * naxis2) - self.overlap)

        # apply definitions
        self.wcs = astropy.wcs.WCS(naxis=2)
        self.wcs.wcs.crpix = [crpix1, crpix2]
        self.wcs.wcs.crval = self.projection_cell.wcs.wcs.crval
        self.wcs.wcs.cd = self.projection_cell.wcs.wcs.cd / ratio
        self.wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
        self.wcs.pixel_shape = (naxis1, naxis2)
        self.wcs.ltv1 = (self.projection_cell.wcs.wcs.crpix[0] * ratio) - crpix1
        self.wcs.ltv2 = (self.projection_cell.wcs.wcs.crpix[1] * ratio) - crpix2
        self.wcs.orientat = self.projection_cell.wcs.orientat
        self.wcs.pscale = self.projection_cell.wcs.pscale / ratio

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
        edges_x = [0] * naxis2 + [naxis1 - 1] * naxis2 + list(range(naxis1)) * 2
        edges_y = list(range(naxis2)) * 2 + [0] * naxis1 + [naxis2 - 1] * naxis1

        polygon = list(zip(edges_x, edges_y))
        img = Image.new("L", (naxis1, naxis2), 0)
        ImageDraw.Draw(img).polygon(polygon, outline=1, fill=1)
        mask = np.array(img)

        self.mask = mask


class SkyCorners(object):
    def __init__(self):
        # attributes keeping track of combined set of inputs
        self.xy_corners = []
        self.sky_corners = []
        self.edges = []
        self.sky_edges = []

        # keep track of how inputs relate to each chip
        self.chips = []

        # attributes related to WCS used to define the output footprint
        self.inside = None
        self.footprint = None
        self.segments = []

    def add_footprint(self, footprint, corners):

        xy_corners = corners.copy()
        xy_corners.append(xy_corners[0])

        sky_corners = footprint.copy()
        sky_corners.append(sky_corners[0])

        edges = []
        sky_edges = []
        # Define edges for each chip based on input corners for the chip
        for i, xy in enumerate(xy_corners[:-1]):
            # Saving edge definitions in 'reverse' order compensates for
            # the clock-wise orientation of the corners from `wcs.calc_footprint()`.
            edge = [xy_corners[i + 1], xy]
            sky_edge = [sky_corners[i + 1], sky_corners[i]]

            chip_edge = {'edge': edge, 'sky_edge': sky_edge}
            edges.append(chip_edge)
            sky_edges.append(sky_edge)
            self.edges.append(chip_edge)

        new_chip = {'sky_corners': sky_corners, 'xy_corners': xy_corners,
                    'edges': edges, 'sky_edges': sky_edges}

        self.chips.append(new_chip)
        self.sky_corners.append(footprint)
        self.xy_corners.append(corners)

    def apply_mask(self, footprint, wcs):
        """Compute segments which make up the outline of the total footprint"""
        # Shrink footprint by 1 pixel, then invert.
        # This results in a mask where pixels inside the footprint (but not
        # the very edge) will be False.
        # This will allow for easy identification of any pixel interior to the footprint
        inside = ndimage.binary_erosion(footprint, iterations=1).astype(np.int16)
        self.inside = trace_polygon(inside - ndimage.binary_erosion(inside).astype(np.int16))


        # Start by identifying what portions of each chip edge makes up
        # the exterior of the total footprint as specified on input
        for edge in self.edges:
            segments, sky_segments = find_segments(edge['edge'][0], edge['edge'][1], self.inside, wcs)
            if len(segments) > 0:
                self.segments.append({'edge': segments, 'sky_edge': sky_segments})

        # Now connect all exterior chip segments into a continuous polygon
        # Start with the corner (start/end of edge) with largest Y value
        start_corner = None
        region = []
        # define what edge will close/finish the polygon when going counter-clockwise
        # This edge defines the last 2 points on the sky for the S_REGION keyword
        last_edge = None
        start_edge = None
        remaining_segments = 0

        for segment in self.segments:
            # The segment may contain 0 or multiple segments
            for start, end in segment['sky_edge']:
                remaining_segments += 1
                # Look for starting corner as corner with largest Dec (closest to +90)
                # which has the smallest RA (if multiple points have identical Dec)
                if (start_corner is not None and ((start[1] > start_corner[1]) or \
                    start[1] == start_corner[1] and start[0] < start_corner[0])) or \
                    start_corner is None:
                    start_corner = start
                    start_edge = [start, end]

                if end[1] > start_corner[1]:
                    start_corner = end
                    last_edge = [start, end]
                    start_edge = None  # Try again to find starting edge

        # After finding starting corner, if last_edge is not None, we know
        # that the starting corner is the 'start' of the segment when
        # going counter-clockwise
        end_corner = None
        if start_edge is None:
            # last_edge will always be defined in this case, so
            # look for edge with start closest to this edge end point
            for edge in self.segments:
                for start, end in self.segments['sky_edge']:
                    if start == last_edge['sky_edge'][1]:
                        start_edge = edge
                        end_corner = end
                        break
        else:
            # start populating the polygon (s_region coordinates)
            end_corner = start_edge['sky_edge'][1]

        # At this point, both the starting edge and end corner are defined.
        # Now start building up region polygon from all exterior segments
        # starting from this edge
        remaining_segments -= 1
        region.append(start_edge)
        region_end = None

        while region_end != end_corner and remaining_segments > 0:
            pass


def cart2pol(x, y, clockwise=False):
    """Converts x,y arrays into radial coordinates of distance and degrees."""
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    deg = 90. - np.rad2deg(phi) if clockwise else np.rad2deg(phi) - 90.
    deg[deg < 0.] += 360.
    return(rho, phi, deg)

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

def detect_corners(edge_pixels):
    """Try to detect pixels based on change of sign of slope"""
    # rotate edge_pixels by a few so starting point is off a corner
    edge = edge_pixels * 0.
    edge[:10] = edge_pixels[-10:]
    edge[10:] = edge_pixels[:-10]

    deltas = edge[2:] - edge[:-2]
    slopes = deltas[:, 1] / deltas[:, 0]

    asign = np.sign(slopes)
    signchange = ((np.roll(asign, 1) - asign) != 0).astype(int)
    corner_indx = np.where(signchange == 1)[0] + 1
    corners = edge[corner_indx]

    return corners

def trace_polygon(input_mask, sensitivity=0.1):
    """Convert mask with only edge pixels into a single contiguous polygon"""
    edge_pixels, mask_residuals = _poly_trace(input_mask)
    edge_lists = [edge_pixels]

    if mask_residuals.sum() > sensitivity * edge_pixels.shape[0]:
        ratio = mask_residuals.sum() / edge_pixels.shape[0]
        # mask_residuals only contains those pixels which were not
        # part of the already traced edge.
        residual_pixels = np.where(mask_residuals == 1)
        residual_pixels = np.array([residual_pixels[1], residual_pixels[0]]).T
        # look for significantly large sets of pixels with consistent distance
        # from the already identified edge_pixels
        residual_dists = distance.cdist(residual_pixels, edge_pixels)
        residual_mins = np.array([resid.min() for resid in residual_dists])
        # count number of edge pixels remaining which are the same distance
        # from pixels in the already identified edge
        num_mins = np.where(np.abs(residual_mins - residual_mins.min()) < 3 * np.sqrt(2))[0].shape[0]
        if num_mins > (sensitivity * ratio * residual_pixels.shape[0]):
            # try tracing the remaining pixels
            edge_pixels, mask_residuals = _poly_trace(mask_residuals)
            if edge_pixels.shape[0] > 0:
                edge_lists.append(edge_pixels)
    return edge_lists

def _poly_trace(input_mask):
    # find a starting point on the mask
    mask = input_mask.copy()
    for x in range(mask.shape[1]):
        pts = np.where(mask[:, x] == 1)[0]
        if len(pts) > 0:
            ystart = pts[0]
            xstart = x
            break
    # Now start tracing looking at a 3x3 set of pixels around that position
    polygon = [[xstart, ystart]]
    # Zero out already identified pixels on the polygon
    mask[ystart, xstart] = 0
    new_x = xstart
    new_y = ystart
    xstart = None
    ystart = None
    new_start = True
    slope = -1
    # determine how many pixels should be in polygon.
    num_pixels = mask.sum()
    while new_start or (num_pixels > 0):
        xstart = new_x
        ystart = new_y
        # Zero out already identified pixels on the polygon
        mask[ystart, xstart] = 0
        box = get_box(mask, xstart, ystart)
        pts = np.where(box == 1)

        if len(pts[0]) == 0:
            dist = distance.cdist([[xstart, ystart]], [polygon[0]])[0][0]
            if dist <= np.sqrt(2):
                # We are back where we started, so quit
                break

        indx = 0
        if len(pts[0]) > 1:
            # Perform some disambiguation to look for
            # pixel which leads to the most pixels going on
            # start with pixels along the same slope that we have been going
            slope_indx = 0 if slope <= 0 else 1
            slope_y = pts[0][slope_indx] + ystart - 1
            slope_x = pts[1][slope_indx] + xstart - 1
            slope_sum = get_box(mask, slope_x, slope_y).sum()
            # Now get sum for the other pixel
            indx2 = 1 if slope < 0 else 0
            y2 = pts[0][indx2] + ystart - 1
            x2 = pts[1][indx2] + xstart - 1
            sum2 = get_box(mask, x2, y2).sum()
            # select point which leads to the largest sum,
            # but favor the previous slope if both directions are equal.
            if slope_sum == sum2:
                indx = slope_indx
            else:
                indx = indx2 if sum2 > slope_sum else slope_indx

        new_y = pts[0][indx] + ystart - 1
        new_x = pts[1][indx] + xstart - 1
        polygon.append([new_x, new_y])
        # reset for next pixel
        num_pixels -= 1
        new_start = False
        slope = (new_y - ystart) / (new_x - xstart)

    # close the polygon
    polygon.append(polygon[0])
    return np.array(polygon, dtype=np.int32), mask

def perp(a):
    b = np.empty_like(a)
    b[0] = -a[1]
    b[1] = a[0]
    return b

# line segment a given by endpoints a1, a2
# line segment b given by endpoints b1, b2
# return
def seg_intersect(start, end):
    a1, a2 = start
    b1, b2 = end
    da = a2 - a1
    db = b2 - b1
    dp = a1 - b1
    dap = perp(da)
    denom = np.dot(dap, db)
    num = np.dot(dap, dp)
    return (num / denom.astype(float)) * db + b1


def get_box(arr, x, y, size=3):
    """subarray extraction with limits checking """
    amin = size // 2
    amax = amin + 1
    ymin = y - amin if y != 0 else 0
    ymax = y + amax if y != arr.shape[0] - 1 else arr.shape[0] - 1
    xmin = x - amin if x != 0 else 0
    xmax = x + amax if x != arr.shape[1] - 1 else arr.shape[1] - 1
    box = arr[ymin: ymax, xmin: xmax]

    return box


def find_segments(start, end, polygon, wcs):
    """Compute the segments of a single edge which make up the footprint

    Parameters
    ----------
    start : tuple
        (X, Y) coordinates for starting corner

    end : tuple
        (X, Y) coordinates for ending corner

    polygon : list or ndarray
        N x 2 list or array of (x, y) pixel positions which trace out the
        outline of the total footprint.  These positions must be in order
        around the polygon.

    Returns
    -------
    segments : list
        List of [(x0,y0), (x1,y1)] pairs representing the start/stop pixels for each exterior
        segment of this edge of the exposure.  This may be empty for an exposure
        completely interior to the overall footprint.

    rd_segments : list
        List of [(RA0,Dec0), (RA1,Dec1)] pairs representing the start/stop pixels for each exterior
        segment of this edge of the exposure.  This may be empty for an exposure
        completely interior to the overall footprint.
    """
    # Convert input polygon pixels into a usable Path object
    poly_path = path.Path(polygon)

    # initialize output values
    segments = []
    sky_segments = []

    # define segment going from starting pixel to ending pixel
    ends = np.array([start, end], dtype=np.int32)


    # Identify all the pixels along the entire length of this segment
    # in the footprint(WCS) specified
    xpts = np.linspace(ends[0, 0], ends[1, 0], abs(ends[1, 0] - ends[0, 0]) + 1)
    slope = (ends[1, 1] - ends[0, 1]) / (ends[1, 0] - ends[0, 0])
    y0 = ends[0, 1] if ends[0, 0] < ends[1, 0] else ends[1, 1]
    x0 = ends[0, 0] if ends[0, 0] < ends[1, 0] else ends[1, 0]
    ypts = y0 + slope * (xpts - x0)

    # Now we have (x,y) pairs for all pixels that comprise this edge in this WCS
    line_pts = np.array([xpts, ypts]).T.astype(np.int32)

    # This limit accounts for aliasing along the edge where an individual pixel
    # may be considered 'inside' when it was merely a matter of integerization of
    # the edge.  This computation also takes into account the orientation of the
    # edge relative to the X,Y pixel grid.
    raw_distances = np.diagonal(distance.cdist(line_pts[1:], line_pts[:-1]))
    dist_limit = raw_distances.min() + raw_distances.max()

    # Apply mask to remove points which are inside the footprint
    masked_line = line_pts[np.invert(poly_path.contains_points(line_pts))]

    # extract segments consisting of starting and ending (x,y) pairs
    # for each contiguous set of pixels remaining (if any)
    if len(masked_line) > 0:
        deltas = masked_line[1:] - masked_line[:-1]
        distances = np.sqrt(np.power(deltas[:, 0], 2) + np.power(deltas[:, 1], 2))
        # distances = np.diagonal(distance.cdist(masked_line[1:], masked_line[:-1]))

        # look for sets of pixels within the 'dist_limit' of each other along the line
        segment_indx = np.where(distances > dist_limit)[0]
        if len(segment_indx) == 0:
            # Only 1 segment
            xy_edge = [masked_line[0], masked_line[-1]]
            rd_edge = wcs.all_pix2world(xy_edge, 0)
            segments.append(xy_edge)
            sky_segments.append(rd_edge)
        else:
            # Multiple segments found
            start_indx = 0
            for i in range(len(segment_indx)):
                xy_edge = [masked_line[start_indx], masked_line[segment_indx[i]]]
                rd_edge = wcs.all_pix2world(xy_edge, 0)
                segments.append(xy_edge)
                sky_segments.append(rd_edge)
                start_indx = segment_indx[i] + 1
            # append last segment
            xy_edge = [masked_line[start_indx], masked_line[-1]]
            rd_edge = wcs.all_pix2world(xy_edge, 0).tolist()
            segments.append(xy_edge)
            sky_segments.append(rd_edge)

    return segments, sky_segments



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
        muji = (sum_x[i] - sum_x[j - 1]) / (i - j + 1)
        sji = sum_x_sq[i] - sum_x_sq[j - 1] - (i - j + 1) * muji ** 2
    else:
        sji = sum_x_sq[i] - sum_x[i] ** 2 / (i + 1)

    return 0 if sji < 0 else sji

def fill_row_k(imin, imax, k, S, J, sum_x, sum_x_sq, N):
    if imin > imax: return

    i = (imin + imax) // 2
    S[k][i] = S[k - 1][i - 1]
    J[k][i] = i

    jlow = k

    if imin > k:
        jlow = int(max(jlow, J[k][imin - 1]))
    jlow = int(max(jlow, J[k - 1][i]))

    jhigh = i - 1
    if imax < N - 1:
        jhigh = int(min(jhigh, J[k][imax + 1]))

    for j in range(jhigh, jlow - 1, -1):
        sji = ssq(j, i, sum_x, sum_x_sq)

        if sji + S[k - 1][jlow - 1] >= S[k][i]: break

        # Examine the lower bound of the cluster border
        # compute s(jlow, i)
        sjlowi = ssq(jlow, i, sum_x, sum_x_sq)

        SSQ_jlow = sjlowi + S[k - 1][jlow - 1]

        if SSQ_jlow < S[k][i]:
            S[k][i] = SSQ_jlow
            J[k][i] = jlow

        jlow += 1

        SSQ_j = sji + S[k - 1][j - 1]
        if SSQ_j < S[k][i]:
            S[k][i] = SSQ_j
            J[k][i] = j

    fill_row_k(imin, i - 1, k, S, J, sum_x, sum_x_sq, N)
    fill_row_k(i + 1, imax, k, S, J, sum_x, sum_x_sq, N)

def fill_dp_matrix(data, S, J, K, N):
    sum_x = np.zeros(N, dtype=np.float_)
    sum_x_sq = np.zeros(N, dtype=np.float_)

    # median. used to shift the values of x to improve numerical stability
    shift = data[N // 2]

    for i in range(N):
        if i == 0:
            sum_x[0] = data[0] - shift
            sum_x_sq[0] = (data[0] - shift) ** 2
        else:
            sum_x[i] = sum_x[i - 1] + data[i] - shift
            sum_x_sq[i] = sum_x_sq[i - 1] + (data[i] - shift) ** 2

        S[0][i] = ssq(0, i, sum_x, sum_x_sq)
        J[0][i] = 0

    for k in range(1, K):
        if (k < K - 1):
            imin = max(1, k)
        else:
            imin = N - 1

        fill_row_k(imin, N - 1, k, S, J, sum_x, sum_x_sq, N)

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
    cluster_right = n - 1

    for cluster in range(n_clusters - 1, -1, -1):
        cluster_left = int(J[cluster][cluster_right])
        clusters.append(data[cluster_left:cluster_right + 1])

        if cluster > 0:
            cluster_right = cluster_left - 1

    return list(reversed(clusters))

##
# # HELPER CODE FOR TESTS
##

# partition recipe modified from
# http://wordaligned.org/articles/partitioning-with-python
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
    splits = (d for i in range(l) for d in combinations(mid, n - 1))
    return [[s[sl] for sl in map(slice, chain(b, d), chain(d, e))]
            for d in splits]

def squared_distance(part):
    mean = sum(part) / len(part)
    return sum((x - mean)**2 for x in part)

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
        1 / 0
    except ValueError:
        pass

    tests = [
        (([1], 1),                          [[1]]),
        (([0, 3, 4], 2),                    [[0], [3, 4]]),
        (([-3, 0, 4], 2),                   [[-3, 0], [4]]),
        (([1, 1, 1, 1], 1),                 [[1, 1, 1, 1]]),
        (([1, 2, 3], 3),                    [[1], [2], [3]]),
        (([1, 2, 2, 3], 3),                 [[1], [2, 2], [3]]),
        (([1, 2, 2, 3, 3], 3),              [[1], [2, 2], [3, 3]]),
        (([1, 2, 3, 2, 3], 3),              [[1], [2, 2], [3, 3]]),
        (([3, 2, 3, 2, 1], 3),              [[1], [2, 2], [3, 3]]),
        (([3, 2, 3, 5, 2, 1], 3),           [[1, 2, 2], [3, 3], [5]]),
        (([0, 1, 2, 100, 101, 103], 2),     [[0, 1, 2], [100, 101, 103]]),
        (([0, 1, 2, 50, 100, 101, 103], 3), [[0, 1, 2], [50], [100, 101, 103]]),
        (([-1, 2, -1, 2, 4, 5, 6, -1, 2, -1], 3),
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
