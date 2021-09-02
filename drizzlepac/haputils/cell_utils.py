import os
import pdb
import shutil
from itertools import chain, combinations

from matplotlib import path
from matplotlib import pyplot as plt
from matplotlib.path import Path
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
        table = Table.read(visit_input, format='ascii.fast_no_header')
        expnames = table['col1'].tolist()

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
        fimg = fits.open(filename)
        print("Validating WCS solutions for {}".format(filename))
        if 'wcsname' not in fimg[1].header:
            expnames.remove(filename)
        fimg.close()
    if len(expnames) == 0:
        print("No valid exposures to define sky cells")
        return None

    # Initialize all sky tessellation object definitions
    # This includes setting the pixel scale.
    sky_grid = GridDefs(scale=scale, cell_size=cell_size)

    # Group input expnames by visit and look for SkyCell overlap
    # on each 'visit' separately, then only update final list of
    # output SkyCells with unique SkyCells merging list of input
    # expnames as needed.  Doing this will minimize the size of
    # the 'meta_wcs' defined for the SkyFootprint used to overlap
    # the SkyCells.  Otherwise, the 'meta_wcs' could become almost
    # arbitrarily large compared to the size of a SkyCell.
    #
    visit_groups = {}
    for exp in expnames:
        exp_visit = extract_visit(exp)
        if exp_visit not in visit_groups:
            visit_groups[exp_visit] = []
        visit_groups[exp_visit].append(exp)
    # at this point, we have a dict of visit IDs with a list of all expnames that go with each visit
    sky_cells = {}
    for visit_id, visit_expnames in visit_groups.items():
        print('Looking for SkyCells that overlap exposures from visit "{}"'.format(visit_id))

        # build reference wcs for combined footprint of all input exposures
        meta_wcs = wcs_functions.make_mosaic_wcs(visit_expnames, rot=0.0, scale=sky_grid.scale)

        print('Visit WCS: \n{}'.format(meta_wcs))

        # create footprint on the sky (as a tangent plane array) for all input exposures using meta_wcs
        footprint = SkyFootprint(meta_wcs)
        footprint.build(visit_expnames)

        # Use this footprint to identify overlapping sky cells
        visit_cells = sky_grid.get_sky_cells(footprint)
        print('Visit "{}" overlapped SkyCells:\n{}'.format(visit_id, visit_cells))

        for scell in visit_cells:
            if scell not in sky_cells:
                sky_cells[scell] = visit_cells[scell]
            else:
                # It overlapped previous exposures, so add these exposures to the SkyCell definition
                sky_cells[scell].members.extend(visit_cells[scell].members)

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
    def __init__(self, meta_wcs, debug=False):

        self.meta_wcs = meta_wcs
        self.debug = debug

        # bounded_wcs corresponds to WCS of bounding box of exposed pixels
        self.bounded_wcs = None
        self.bounding_box = None

        # the exp_masks dict records the individual footprints of each exposure
        self.exp_masks = {}
        self.members = []
        self.corners = []
        self.sky_corners = []
        self.xy_corners = []
        self.edge_pixels = []

        self.total_mask = np.zeros(meta_wcs.array_shape, dtype=np.int16)
        self.scaled_mask = None
        self.footprint = None
        self.footprint_member = None

        self.edges = None
        self.edges_ra = None
        self.edges_dec = None
        self.polygon = None

    def build(self, expnames, scale=False, scale_kw='EXPTIME'):
        """ Create mask showing where all input exposures overlap the SkyFootprint's WCS

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
            self.exp_masks[exposure] = {'sky_corners': [], 'xy_corners': [], 'mask': {}}
            exp = fits.open(exposure)
            if scale:
                # Only assign memory for this array if requested.
                self.scaled_mask = np.zeros(self.meta_wcs.array_shape, dtype=np.float32)
            scale_val = exp[0].header[scale_kw]

            sci_extns = wcs_functions.get_extns(exp)
            if len(sci_extns) == 0 and '_single' in exposure:
                sci_extns = [0]

            for sci in sci_extns:
                wcs = HSTWCS(exp, ext=sci)

                # save the footprint for each chip as RA/Dec corner positions
                # radec = wcs.calc_footprint().tolist()
                radec = calc_wcs_footprint(wcs, offset=1).tolist()
                radec.append(radec[0])  # close the polygon/chip
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

                # check to see whether or not this image falls within meta_wcs at all...
                off_x = (meta_x.max() <= 0) or meta_x.min() > (self.meta_wcs.array_shape[1] - 1)
                off_y = (meta_y.max() <= 0) or meta_y.min() > (self.meta_wcs.array_shape[0] - 1)
                # if this chip falls completely outside SkyCell, skip to next chip
                if off_x or off_y:
                    continue

                # Account for rounding problems with creating meta_wcs
                meta_y = np.clip(meta_y, 0, self.meta_wcs.array_shape[0] - 1)
                meta_x = np.clip(meta_x, 0, self.meta_wcs.array_shape[1] - 1)

                # define subarray spanned by this chip on the SkyCell
                scell_slice = [slice(meta_y.min(), meta_y.max()),
                               slice(meta_x.min(), meta_x.max())]
                scell_ltm = [meta_x.min(), meta_y.min()]
                # Reset range of pixels to be relative to starting pixel position in SkyCell
                meta_x -= scell_ltm[0]
                meta_y -= scell_ltm[1]

                # apply meta_edges to blank mask
                # Use PIL to create mask
                # parray = np.array(meta_edges.T)
                parray = (meta_x, meta_y)
                polygon = list(zip(parray[0], parray[1]))
                nx = self.total_mask[tuple(scell_slice)].shape[1]
                ny = self.total_mask[tuple(scell_slice)].shape[0]
                if nx == 0 or ny == 0:
                    continue
                img = Image.new('L', (nx, ny), 0)
                ImageDraw.Draw(img).polygon(polygon, outline=1, fill=1)
                blank = np.array(img).astype(np.int16)

                # Remember information needed to recreate the mask for this chip
                sci_dict = {}
                sci_dict['scell_slice'] = scell_slice
                sci_dict['scell_ltm'] = scell_ltm
                sci_dict['scale_val'] = scale_val
                sci_dict['polygon'] = polygon
                sci_dict['img_shape'] = (nx, ny)
                self.exp_masks[exposure]['mask'][sci] = sci_dict

                if scale:
                    scaled_blank = blank * scale_val
                    self.scaled_mask[tuple(scell_slice)] += scaled_blank

                self.total_mask[tuple(scell_slice)] += blank
                del blank

            # Only add members which contributed to this footprint
            if exposure not in self.members:
                self.members.append(exposure)

        # Compute the bounded WCS for this mask of exposed pixels
        self.find_bounded_wcs()

    def extract_mask(self, filename):
        """Extract a total_mask from the SCI data directly"""
        # Determine what extension contains the SCI array
        # It could be 0 if 'build=no' for drizzling
        fhdu = fits.open(filename)
        if len(fhdu) > 1:
            sciext = ("SCI", 1)
        else:
            sciext = 0
        # Get the SCI array to use as a mask
        arr = fhdu[sciext].data.copy()
        # Done with image, so close immediately.
        fhdu.close()
        del fhdu

        # If working with drizzled data which has NaN as non-exposed pixel values
        if np.isnan(arr.min()):
            total_mask = (np.isnan(arr) == 0).astype(np.int16)
        else:
            total_mask = (arr != 0).astype(np.int16)
        # Remove 'small' holes in the image due to noise to avoid
        # creating extraneous 'regions' from the image when it is really
        # all one region/chip.
        total_mask = ndimage.binary_fill_holes(total_mask)
        self.total_mask = total_mask

        # clean up as quickly as possible
        del arr

        # populate footprint with this mask
        self.members += [filename]
        self.find_footprint()

        # Now compute the bounded_wcs, if possible.
        self.find_bounded_wcs()

    def find_bounded_wcs(self):
        """Compute the WCS based on the bounding box of exposed pixels(mask) """
        if self.total_mask is None:
            print("Please add exposures before computing bounding box WCS...")

        # start by computing the bounding box for the footprint
        ymin, ymax, xmin, xmax = calc_bounding_box(self.total_mask)
        # make a copy of the full WCS to be revised
        self.bounded_wcs = self.meta_wcs.deepcopy()
        self.bounding_box = [slice(ymin, ymax), slice(xmin, xmax)]

        # Use this box to compute new CRPIX position
        self.bounded_wcs.wcs.crpix -= [xmin, ymin]
        self.bounded_wcs.pixel_shape = [xmax - xmin + 1, ymax - ymin + 1]


    # Methods with 'find' compute values
    # Methods with 'get' return values
    def find_footprint(self, member='total'):
        """Compute a mask for the footprint

        This method converts the full mask built up from all the input
        exposures, or for a single exposure, and converts it to a
        boolean mask.  The mask has a value of 1 where a pixel was
        part of an exposure.

        The resulting mask gets saved as the ``footprint`` attribute.

        Parameters
        ==========
        member : `str`, optional
            Specify what member to compute the footprint for.

        """

        if self.total_mask is None:
            print("Please add exposures before computing footprint...")
        if member == 'total':
            mask = self.total_mask
        else:
            if member not in self.exp_masks:
                raise ValueError("Member {} not added to footprint".format(member))
            # Recompute mask specifically for this member
            #mask = self.exp_masks[member]['mask']
            exp_mask = self.exp_masks[member]['mask']
            mask = np.zeros(self.meta_wcs.array_shape, dtype=np.int16)

            for sci in exp_mask.values():
                img = Image.new('L', sci['img_shape'], 0)
                ImageDraw.Draw(img).polygon(sci['polygon'], outline=1, fill=1)
                blank = np.array(img).astype(np.int16)
                mask[sci['scell_slice']] += blank

        self.footprint = np.clip(mask, 0, 1)
        self.footprint_member = member


    def find_edges(self, member='total'):
        """Computes a mask containing only those pixels along the edge of the footprint."""

        if self.footprint_member != member:
            self.find_footprint(member=member)

        edges = ndimage.binary_erosion(self.footprint).astype(np.int16)
        self.edges = self.footprint - edges


    def find_corners(self, member='total'):
        """Extract corners/vertices from footprint

        This method computes the positions of all the corners that
        make up the footprint.  The corners are computed starting with
        the corner nearest (within 45deg) of vertical as measured from
        the center of the footprint, then proceeds counter-clockwise
        (North to East).  The corners are initially identified using
        the `skimage.corner_harris` function on the footprint mask to
        identify the starting corner which is closest to veritical.  The
        edge pixels are then ordered counter-clockwise, and corners are
        finally confirmed in order where the slope along each edge changes
        sign.

        This results in a list of corner positions
        which can be used to populate the `S_REGION` keyword and traces
        out the outline of the footprint.

        The results are saved as the attributes:

          * ``edge_pixels`` : the complete list of all pixels, as a list
            of `numpy.ndarray` instances, that make up the edge of the
            footprint in counter-clockwise order.  For
            multiple chips in the footprint, there will be a separate
            array of the ordered pixel positions per chip.
          * ``xy_corners`` : `numpy.ndarray` of (X,Y) positions for all
            identified corners from the footprint.
          * ``corners`` : `numpy.ndarray` of (RA, Dec) positions for
            all identified corners from the footprint.

        """
        if len(self.members) == 0:
            print("Please add exposures before looking for corners...")
            return

        # Insure footprint has been determined
        if self.footprint_member != member:
            self.find_footprint(member=member)

        if member == 'total':
            # Use Harris corner detection to identify corners in the
            # total footprint
            # insure footprint has enough signal to detect corners
            fp = np.clip(self.footprint, 0, 1).astype(np.int16)

            # simple trick to remove noise and small regions 3x3 or less.
            scmask = ndimage.binary_dilation(ndimage.binary_erosion(fp, iterations=3), iterations=2)
            # Label each major contiguous region in the mask
            sclabels, nlabels = ndimage.label(scmask)
            slices = ndimage.find_objects(sclabels)

            # For each region, trace the edge, find the Harris corners,
            # then order the Harris corners counter-clockwise around the region
            # using the traced edge pixel positions.
            ordered_xy = []
            ordered_edges = []
            sky_corners = []

            for label, mask_slice in enumerate(slices):
                label += 1
                # get slice with just the region/label of interest
                label_mask = sclabels[mask_slice].copy()
                # make sure no pixels from other regions are present in this mask
                label_mask[label_mask != label] = 0
                # reset label to be a binary mask only
                label_mask[label_mask == label] = 1000
                # insure there is a border all around the region
                # THIS IS CRITICAL to being able to identify corners correctly in slice
                label_mask = ndimage.binary_erosion(label_mask)
                print('extracting corners for region {} in slice {}'.format(label, mask_slice))
                # Perform corner detection on each region/chip separately.
                mask_corners = corner_peaks(corner_harris(label_mask),
                                       min_distance=1,
                                       threshold_rel=0.2)
                xy_corners = mask_corners * 0.
                xy_corners[:, 0] = mask_corners[:, 1]
                xy_corners[:, 1] = mask_corners[:, 0]
                # shift corner positions to full array positions
                xy_corners += (mask_slice[1].start, mask_slice[0].start)

                # Create a mask from the total footprint consisting solely of the
                # pixels at the outer edge, ordered in clockwise fashion.
                #
                # get list of (X,Y) coordinates of all edges from each separate 'region' or chip
                edge_pixels = trace_polygon(label_mask, mask_slice)

                # use the ordering of the traced edge pixels to order the corners in the same way
                cordist = distance.cdist(xy_corners, edge_pixels)  # returns distances for each corner position
                ordered_indices = []
                for distarr, minval in zip(cordist, np.min(cordist, axis=1)):
                    ordered_indices.append(np.where(distarr == minval)[0][0])
                radial_order = np.argsort(ordered_indices)
                ordered_xyc = xy_corners[radial_order].tolist()
                ordered_xyc.append(ordered_xyc[0])  # close polygon

                # save as output values
                ordered_xy.append(np.array(ordered_xyc, dtype=np.float64))
                sky_corners.append(self.meta_wcs.all_pix2world(ordered_xyc, 0))
                ordered_edges.append(edge_pixels)

        else:
            if member not in self.exp_masks:
                raise ValueError("Member {} not added to footprint".format(member))
            xy_corners = self.exp_masks[member]['xy_corners']
            sky_corners = self.meta_wcs.all_pix2world(xy_corners, 0)

        self.edge_pixels = np.array(ordered_edges)
        self.xy_corners = np.array(ordered_xy)
        self.corners = np.array(sky_corners)

    def get_edges_sky(self, member='total'):
        """Returns the sky coordinates of all edge pixels.

        This method uses the WCS of the footprint to convert the (X,Y)
        pixel positions of all the edge pixels and converts them to
        world coordinates.

        The results are saved as the attributes:

          * ``edges_ra``
          * ``edges_dec``

        """
        self.find_footprint(member=member)
        if len(self.corners) == 0:
            self.find_corners()
        self.edges_ra, self.edges_dec = self.meta_wcs.pixel_to_world_values(self.edge_pixels[0][:,0], self.edge_pixels[0][:,1])
        return self.edges_ra, self.edges_dec

    def build_polygon(self, member='total'):
        """Convert the edges into a SphericalPolygon object.

        This method converts the sky coordinates of the footprint edges
        into a `spherical_geometry.SphericalPolygon` instance.

        The results are saved as the ``polygon`` attribute.

        """

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
        """Convert the footprint into a FITS HDUList object

        Parameters
        ----------
        filename : `str`, optional
            If specified, write out the object to the specified FITS file.

        overwrite : `bool`, optional
            Specify whether or not to overwrite a previously written footprint FITS file.

        member : `str`, optional
            Name of member, or 'total', footprint to write out as FITS HDUList object.

        Returns
        --------
        hdu : `fits.PrimaryHDU`
            FITS HDU containing the footprint
        """
        self.find_footprint(member=member)
        return self._get_fits_hdu(self.footprint, filename=filename, overwrite=overwrite)

    def get_edges_hdu(self, filename=None, overwrite=True, member='total'):
        """Convert the edge pixels mask into a FITS HDUList object

        Parameters
        ----------
        filename : `str`, optional
            If specified, write out the object to the specified FITS file.

        overwrite : `bool`, optional
            Specify whether or not to overwrite a previously written mask FITS file.

        member : `str`, optional
            Name of member, or 'total', footprint to write out as FITS HDUList object.

        Returns
        --------
        hdu : `fits.PrimaryHDU`
            FITS HDU containing the edge pixels mask
        """
        self.find_edges(member=member)
        return self._get_fits_hdu(self.edges, filename=filename, overwrite=overwrite)

    def get_mask_hdu(self, filename=None, overwrite=True):
        """Convert the total mask into a FITS HDUList object

        The `total mask` attribute represents the number of exposures per pixel for the mosaic,
        optionally scaled by the exposure time.  This gets written out as a FITS
        PrimaryHDU object by this method.

        Parameters
        ----------
        filename : `str`, optional
            If specified, write out the object to the specified FITS file.

        overwrite : `bool`, optional
            Specify whether or not to overwrite a previously written mask FITS file.

        member : `str`, optional
            Name of member, or 'total', footprint to write out as FITS HDUList object.

        Returns
        --------
        hdu : `fits.PrimaryHDU`
            FITS HDU containing the total mask
        """
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
        ring_id = np.where(self.rings['projcell'] <= id)[0][-1]
        return self.rings[ring_id]

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

        # Get the edges of the mosaic on the sky
        mosaic_ra, mosaic_dec = mosaic.get_edges_sky()
        # Convert edges to positions in projection cell
        mosaic_edges_x, mosaic_edges_y = self.wcs.world_to_pixel_values(mosaic_ra, mosaic_dec)

        # Determine roughly what sky cells overlap this mosaic
        # The 1.5 enforces 1-based indexing for the SkyCell IDs within the Projection Cell
        mosaic_edges_x = (mosaic_edges_x / skycell00.wcs.pixel_shape[0] + 1.5).astype(np.int32)
        mosaic_edges_y = (mosaic_edges_y / skycell00.wcs.pixel_shape[1] + 1.5).astype(np.int32)
        # Define range of SkyCells in X,Y to evaluate for overlap
        # Enforce a range of SkyCell IDs that are 1-based from 1 to self.sc_nxy
        mosaic_xr = [max(1, mosaic_edges_x.min() - 1), min(self.sc_nxy, mosaic_edges_x.max() + 2)]
        mosaic_yr = [max(1, mosaic_edges_y.min() - 1), min(self.sc_nxy, mosaic_edges_y.max() + 2)]

        print("SkyCell Ranges: {}, {}".format(mosaic_xr, mosaic_yr))
        # for each suspected sky cell or neighbor, look for any pixel by pixel
        #    overlap with input mosaic footprint
        for xi in range(mosaic_xr[0], mosaic_xr[1]):
            for yi in range(mosaic_yr[0], mosaic_yr[1]):
                skycell = SkyCell(x=xi, y=yi, projection_cell=self)

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

        # Initialize computed attributes
        self._build_wcs()  # compute .wcs and .corners
        self.mask = None
        self.polygon = None

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
        pc_nx = self.projection_cell.wcs.pixel_shape[0]
        pc_ny = self.projection_cell.wcs.pixel_shape[1]

        # Define size of SkyCells at default/fine plate scale
        # CRPIX of SkyCells should always be at exactly ((pc_nx/self.nxy) * ratio) apart
        # Size of SkyCells needs to self.overlap * 2 larger than the distance between CRPIX values
        sc_nx1 = int(pc_nx / self.nxy + 0.5)
        sc_nx2 = int(pc_ny / self.nxy + 0.5)
        naxis1 = int((sc_nx1 + self.overlap * 2) * ratio)
        naxis2 = int((sc_nx2 + self.overlap * 2) * ratio)
        xindx = self.x_index - 1 if self.x_index > 0 else 0
        yindx = self.y_index - 1 if self.y_index > 0 else 0

        crpix1 = ((self.projection_cell.wcs.wcs.crpix[0] - xindx * sc_nx1) ) * ratio
        crpix2 = ((self.projection_cell.wcs.wcs.crpix[1] - yindx * sc_nx2) ) * ratio

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

def calc_bounding_box(img):
    """Compute bounding box for non-zero section of img """
    rows = np.any(img, axis=1)
    cols = np.any(img, axis=0)
    rmin, rmax = np.where(rows)[0][[0, -1]]
    cmin, cmax = np.where(cols)[0][[0, -1]]

    return rmin, rmax, cmin, cmax


def cart2pol(x, y, clockwise=False):
    """Converts x,y arrays into radial coordinates of distance and degrees."""
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    deg = 90. - np.rad2deg(phi) if clockwise else np.rad2deg(phi) - 90.
    deg[deg < 0.] += 360.
    return(rho, phi, deg)

def pol2cart(rho, phi):
    """Converts radial coordinates of distance and radians into x,y arrays."""
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

def calc_wcs_footprint(wcs, offset=1, sample=4):
    """Calculate the sky coordinates for positions inside each corner."""
    delta_x = wcs.naxis1 // sample
    delta_y = wcs.naxis2 // sample

    new_corners = []
    # along X edge at Y=0
    new_corners = [[offset, offset]]
    for ix in range(0, sample):
        new_corners.append([(delta_x * ix) - offset, offset])
    # along Y edge at X=naxis1
    for iy in range(0, sample):
        new_corners.append([wcs.naxis1 - offset, (delta_y * iy) - offset])
    # along top X edge
    for ix in range(sample, 0, -1):
        new_corners.append([(delta_x * ix) - offset, wcs.naxis2 - offset])
    # along Y edge at X=0
    for iy in range(sample, 0, -1):
        new_corners.append([offset, (delta_y * iy) - offset])

    new_sky = wcs.all_pix2world(new_corners, 0)

    return new_sky

def find_vertices(edge_path, exp_masks):
    """Find vertices the input exposures that overlap this part of the footprint.

    For each input exp_mask, check to see whether it intersects the
    Path of the ordered_edge pixels
    if so, get the list of pixel positions for each edge of the exp_mask
    then determine what points of that line are OUTSIDE (not contained by)
    the Path of the ordered_edge pixels
    the first and last position of all segments outside the Path are the vertices
    only 1 entry for each vertex will be kept in the final list based on
    looking at the cdist the positions relative to each other and deleting
    any duplicate within 1 pixel of another.

    """
    vertices = []
    for image in exp_masks:
        for chip in image:
            chip_path = Path(chip)
            if chip_path.intersects_path(edge_path):
                num_pts = len(chip)
                for indx in range(num_pts - 1):
                    # convert each edge into a list of pixel positions (float values, not int)
                    edge_line = compute_edge(chip[indx], chip[indx + 1], int_pixel=False)
                    edge_bool = edge_path.contains_points(edge_line)
                    edge_pixels = edge_line.astype(np.int32)
                    # Only report integer pixel positions to make it easier to weed out duplicates
                    edge_outside = np.where(~edge_bool)[0]
                    if len(edge_outside) > 0:
                        # start with the first pixel on the line that along the outer edge.
                        vertices.append(edge_pixels[edge_outside[0]].tolist())
                        # Now look for the rest of the segments
                        for i, pos in enumerate(edge_outside[:-1]):
                            if edge_outside[i + 1] - edge_outside[i] == 1:
                                continue
                            else:
                                vertices.append(edge_pixels[edge_outside[i]].tolist())

    # Now, we need to remove duplicate vertices
    all_vertices = set(map(tuple, vertices))
    vertices_arr = np.array(list(all_vertices))

    return vertices_arr


def compute_edge(start, end, int_pixel=True):
    """"Return a list of nb_points equally spaced points between start and end"""
    # If we have 8 intermediate points, we have 8+1=9 spaces
    # between p1 and p2
    nb_points = max(abs(end[0] - start[0]), abs(end[1] - start[1]))
    x_spacing = (end[0] - start[0]) / (nb_points + 1)
    y_spacing = (end[1] - start[1]) / (nb_points + 1)

    xy = [[start[0] + i * x_spacing, start[1] + i * y_spacing]
            for i in range(1, nb_points + 1)]

    if int_pixel:
        xy = np.array(xy, dtype=np.int32)
    else:
        xy = np.array(xy)

    return xy


def trace_polygon(input_mask, mask_slice):
    """Convert mask with only edge pixels into a single contiguous polygon"""

    # now extract the edges to trace
    slice_edges = input_mask.astype(np.int16) - \
                  ndimage.binary_erosion(input_mask).astype(np.int16)

    # trace edge from this region and identify corners of edge
    edge_pixels = _poly_trace(slice_edges)
    # shift the origin of the edge pixels to coincide with the slice from the full mask
    edge_pixels += (mask_slice[1].start, mask_slice[0].start)

    return edge_pixels


def _poly_trace(input_mask, box_size=3):
    # find a starting point on the mask
    # expand input array to include empty border pixels to account for extraction box
    border = (box_size // 2)
    mask = np.zeros((input_mask.shape[0] + (border * 2), input_mask.shape[1] + (border * 2)),
                    dtype=input_mask.dtype)
    mask[border:-border, border:-border] = input_mask.copy()

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
            # try a larger box to see if we can jump this gap
            box = get_box(mask, xstart, ystart, size=5)
            pts = np.where(box == 1)
            if len(pts[0]) == 0:
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

    del mask
    # close the polygon
    polygon.append(polygon[0])
    return np.array(polygon, dtype=np.int32)

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
    ymin = y - amin if y >= 0 else 0
    ymax = y + amax if y <= arr.shape[0] else arr.shape[0]
    xmin = x - amin if x >= 0 else 0
    xmax = x + amax if x <= arr.shape[1] else arr.shape[1]
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
        naxis1, naxis2 = (np.array(c1_pixels).sum(axis=0) + 1).astype(int)
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

def extract_visit(filename):
    """Extract the VISIT ID from the input filename"""
    if filename.startswith('hst_'):
        rootname = filename.split('_')[-2]
        visit_id = rootname[:6]
    else:
        visit_id = filename[:6]

    return visit_id
