"""
Classes to keep track of all WCS and catalog information.
Used by ``TweakReg``.

:Authors: Warren Hack, Mihai Cara

:License: :doc:`LICENSE`

"""
import os
import sys
import copy
import numpy as np

from astropy import wcs as pywcs
import stwcs
from astropy.io import fits
from spherical_geometry.polygon import SphericalPolygon
from stsci.skypac.parseat import FileExtMaskInfo, parse_cs_line
from stsci.skypac import utils as spu

from stwcs.distortion import utils
from stwcs.wcsutil import wcscorr
from stwcs.wcsutil import headerlet
from stwcs.wcsutil import altwcs
from stsci.tools import fileutil as fu
from stsci.stimage import xyxymatch
from stsci.tools import logutil, textutil
try:
    from stsci.tools.bitmask import interpret_bit_flags
except ImportError:
    from stsci.tools.bitmask import interpret_bits_value as interpret_bit_flags

from . import catalogs
from . import linearfit
from . import updatehdr
from . import util
from . import tweakutils
from . import wcs_functions

# DEBUG
IMGCLASSES_DEBUG = False

# use convex hull for images? (this is tighter than chip's bounding box)
IMAGE_USE_CONVEX_HULL = True

log = logutil.create_logger(__name__, level=logutil.logging.NOTSET)

sortKeys = ['minflux', 'maxflux', 'nbright', 'fluxunits']

class Image:
    """ Primary class to keep track of all WCS and catalog information for
        a single input image. This class also performs all matching and fitting.
    """
    def __init__(self,filename,input_catalogs=None,exclusions=None,**kwargs):
        """
        Parameters
        ----------
        filename : str
            Filename for image.

        input_catalogs : list of str or None
            Filename of catalog files for each chip, if specified by user.

        kwargs : dict
            Parameters necessary for processing derived from input configObj object.

        """
        self._im = spu.ImageRef()
        self._dq = spu.ImageRef()
        self.dqbits = interpret_bit_flags(kwargs['dqbits'])

        if 'use_sharp_round' in kwargs:
            self.use_sharp_round = kwargs['use_sharp_round']
        else:
            self.use_sharp_round = False

        self.perform_update = kwargs['updatehdr']
        if self.perform_update:
            self.open_mode = 'update'
        else:
            self.open_mode = 'readonly'

        self.name = filename
        self.filename,self.ext_root = fu.parseFilename(filename)
        self.openFile(openDQ=(self.dqbits is not None))

        if self.ext_root is not None:
            num_sci = 1
        else:
            # try to verify whether or not this image has been updated with
            # a full distortion model
            num_sci = spu.count_extensions(self._im, extname='SCI')

        numwht = spu.count_extensions(self._im, extname='WHT')

        wcsextn = util.findWCSExtn(filename)
        #wnames = altwcs.wcsnames(self._im.hdu, ext=(('SCI',1) if num_sci > 0 else 0))
        # If no WCSNAME keywords were found, raise the possibility that
        # the images have not been updated fully and may result in inaccurate
        # alignment
        # use 'numwht' != 0 to indicate a DRZ file has been specified as input
        #if len(wnames) == 0 and numwht == 0:
        if wcsextn is None:
            print(textutil.textbox('WARNING:\n'
            'Image %s may not have the full correct '%filename+
            'WCS solution in the header as created by stwcs.updatewcs '
            'Image alignment may not be taking into account '
            'the full distortion solution.\n'
            'Turning on the "updatewcs" parameter would insure '
            'that each image uses the full distortion model when '
            'aligning this image.\n', width=60
            ))

        self.rootname = os.path.splitext(os.path.basename(self.filename))[0]
        self.origin = 1
        self.pars = kwargs
        self.exclusions = exclusions
        self.verbose = kwargs['verbose']
        self.interactive = kwargs.get('interactive',True)

        # Record this for use with methods
        if num_sci > 0:
            self.nvers = num_sci
            self.ext_name = 'SCI'
        else:
            self.nvers = 1
            self.ext_name = 'PRIMARY'

        # WCS required, so verify that we can get one
        # Need to count number of SCI extensions
        #  (assume a valid WCS with each SCI extension)
        #TODO: current check for a valid WCS may need a revision to
        # implement a more robust/rigurous check.
        # This verification has already been done when finding 'wcsextn'


        # Need to generate a separate catalog for each chip
        self.chip_catalogs = {}
        xypostypes = 3*[float] + [int] + \
            (3 if self.use_sharp_round else 0)*[float] + [object]
        self.xy_catalog = [np.empty(0, dtype=i) for i in xypostypes]

        self.num_sources = 0

        # Analyze exclusion file list:
        if exclusions is not None:
            nexclusions = len(exclusions)
            if nexclusions >= self.nvers:
                exclusions = exclusions[:self.nvers]
            else:
                exclusions += (self.nvers-nexclusions) * [ None ]
        else:
            exclusions = self.nvers * [ None ]

        # For each SCI extension, generate a catalog and WCS
        print("===  Source finding for image '{}':".format(filename))

        bounding_polygons = []

        chip_filenames = {}
        for sci_extn in range(1,self.nvers+1):
            extnum = fu.findExtname(self._im.hdu, self.ext_name, extver=sci_extn)
            if extnum is None:
                extnum = 0
            chip_filenames[sci_extn] = "{:s}[{:d}]".format(self.filename, extnum)

        for sci_extn in range(1,self.nvers+1):
            chip_filename = chip_filenames[sci_extn]
            wcs = stwcs.wcsutil.HSTWCS(chip_filename)

            if input_catalogs is None:
                # if we already have a set of catalogs provided on input,
                #  we only need the array to get original XY input positions
                source = chip_filename
                catalog_mode='automatic'
            else:
                source = input_catalogs[sci_extn-1]
                catalog_mode='user'

            if exclusions[sci_extn-1] not in [None, 'None', '', ' ', 'INDEF']:
                excludefile = {'region_file': exclusions[sci_extn-1]}
            else:
                excludefile = None

            kwargs['start_id'] = self.num_sources
            catalog = catalogs.generateCatalog(wcs, mode=catalog_mode,
                        catalog=source, src_find_filters=excludefile, **kwargs)

            # creaate DQ mask:
            if self.dqbits is None:
                mask = None
            else:
                # make sure DQ data are available:
                dq_available = False
                if not self._dq.closed and self._dqext is not None and \
                   len(self._dqext) > 0:
                    try:
                        data = self._dq.hdu[self._dqext[sci_extn-1]].data
                        dq_available = True
                    except KeyError:
                        pass

                if dq_available:
                    mask = np.logical_not(np.bitwise_and(data, ~self.dqbits)
                                          ).astype(np.uint8)
                    del data
                else:
                    mask = None
                    print(textutil.textbox(
                            "WARNING: User requested to use "
                            "DQ mask for source finding, but DQ data are "
                            "not available.\n"
                            "DQ mask WILL NOT be used for source finding.",
                            indent = 5), file=sys.stderr)

            # read in and convert all catalog positions to RA/Dec
            catalog.buildCatalogs(exclusions=None, mask=mask)

            self.num_sources += catalog.num_objects
            self.chip_catalogs[sci_extn] = {'catalog':catalog,'wcs':wcs}

            # Merge input X,Y positions from all chips into a single catalog
            # This will be used for writing output matched source catalogs
            nxypos = 0 if catalog.xypos is None else min(len(self.xy_catalog),len(catalog.xypos))
            if nxypos > 0:
                nsrc = catalog.xypos[0].shape[0]
                for i in range(nxypos):
                    self.xy_catalog[i] = np.append(self.xy_catalog[i], catalog.xypos[i])

                # add "source origin":
                self.xy_catalog[-1] = np.append(self.xy_catalog[-1],
                                                np.asarray(nsrc*[source]))

                # Compute bounding convex hull for the reference catalog:
                if IMAGE_USE_CONVEX_HULL and self.xy_catalog is not None:
                    # if catalog.xypos[0].shape[0] < 3:
                    xy_vertices = np.asarray(convex_hull(
                        list(map(tuple,np.asarray([catalog.xypos[0],catalog.xypos[1]]).transpose()))),
                                                 dtype=np.float64)
                    if xy_vertices.shape[0] > 2:
                        rdv = wcs.all_pix2world(xy_vertices, 1)
                        bounding_polygons.append(SphericalPolygon.from_radec(rdv[:,0], rdv[:,1]))
                    else:
                        bounding_polygons.append(SphericalPolygon.from_wcs(wcs))

                    if IMGCLASSES_DEBUG:
                        all_ra, all_dec = wcs.all_pix2world(
                            catalog.xypos[0], catalog.xypos[1], 1)
                        _debug_write_region_fk5('dbg_'+self.rootname+'_bounding_polygon.reg',
                                                list(zip(*rdv.transpose())),
                                                list(zip(*[catalog.radec[0], catalog.radec[1]])),
                                                append=sci_extn > 1)
                else:
                    bounding_polygons.append(SphericalPolygon.from_wcs(wcs))

        npoly = len(bounding_polygons)
        if npoly > 1:
            self.skyline = SphericalPolygon.multi_union(bounding_polygons)
        elif npoly == 1:
            self.skyline = bounding_polygons[0]
        else:
            self.skyline = SphericalPolygon([])

        self.catalog_names = {}
        # Build full list of all sky positions from all chips
        self.buildSkyCatalog()
        tsrc = 0 if self.all_radec is None else len(self.all_radec[0])
        print("===  FINAL number of objects in image '{:s}': {:d}"
              .format(filename, tsrc))

        if self.pars['writecat']:
            catname = self.rootname+"_sky_catalog.coo"
            self.catalog_names['match'] = self.rootname+"_xy_catalog.match"
            self.write_skycatalog(catname)
            self.catalog_names['sky'] = catname # Keep track of catalogs being written out
            for nsci in range(1,self.nvers+1):
                catname = "%s_sci%d_xy_catalog.coo"%(self.rootname,nsci)
                self.chip_catalogs[nsci]['catalog'].writeXYCatalog(catname)
                # Keep track of catalogs being written out
                if 'input_xy' not in self.catalog_names:
                    self.catalog_names['input_xy'] = []
                self.catalog_names['input_xy'].append(catname)
            self.catalog_names['fitmatch'] = self.rootname+"_catalog_fit.match"

        print()

        # Set up products which need to be computed by methods of this class
        self.outxy = None
        self.refWCS = None # reference WCS assigned for the final fit
        self.matches = {'image':None,'ref':None} # stores matched list of coordinates for fitting
        self.fit = None # stores result of fit
        self.match_pars = None
        self.fit_pars = None
        self.identityfit = False # set to True when matching/fitting to itself
        self.goodmatch = True # keep track of whether enough matches were found for a fit
        self.figure_id = 1

        self.next_key = ' '

        self.quit_immediately = False
        # close file handle... for now.
        self.close()

    def close(self):
        """ Close any open file handles and flush updates to disk
        """
        self._im.release()
        self._dq.release()

    def openFile(self, openDQ=False):
        """ Open file and set up filehandle for image file
        """

        if self._im.closed:
            if not self._dq.closed:
                self._dq.release()
            assert(self._dq.closed)
            fi = FileExtMaskInfo(clobber=False,
                                 doNotOpenDQ=not openDQ,
                                 im_fmode=self.open_mode)
            fi.image = self.name
            self._im = fi.image
            fi.append_ext(spu.get_ext_list(self._im, extname='SCI'))
            fi.finalize()
            self._im = fi.image
            self._dq = fi.DQimage
            self._imext = fi.fext
            self._dqext = fi.dqext

    def get_wcs(self):
        """ Helper method to return a list of all the input WCS objects associated
            with this image.
        """
        wcslist = []
        for chip in self.chip_catalogs:
            wcslist.append(self.chip_catalogs[chip]['wcs'])
        return wcslist

    def buildSkyCatalog(self):
        """ Convert sky catalog for all chips into a single catalog for
            the entire field-of-view of this image.
        """
        self.all_radec = None
        self.all_radec_orig = None
        ralist = []
        declist = []
        fluxlist = []
        idlist = []
        for scichip in self.chip_catalogs:
            skycat = self.chip_catalogs[scichip]['catalog'].radec
            xycat = self.chip_catalogs[scichip]['catalog'].xypos
            if skycat is not None:
                ralist.append(skycat[0])
                declist.append(skycat[1])
                if xycat is not None and len(xycat) > 2:
                    fluxlist.append(xycat[2])
                    idlist.append(xycat[3])
                elif len(skycat) > 2:
                    fluxlist.append(skycat[2])
                    idlist.append(skycat[3])
                else:
                    fluxlist.append([999.0]*len(skycat[0]))
                    idlist.append(np.arange(len(skycat[0])))

                self.all_radec = [np.concatenate(ralist),np.concatenate(declist),
                        np.concatenate(fluxlist),np.concatenate(idlist)]
                self.all_radec_orig = copy.deepcopy(self.all_radec)

    def buildDefaultRefWCS(self):
        """ Generate a default reference WCS for this image. """
        self.default_refWCS = None
        if self.use_wcs:
            wcslist = []
            for scichip in self.chip_catalogs:
                wcslist.append(self.chip_catalogs[scichip]['wcs'])
            self.default_refWCS = utils.output_wcs(wcslist)

    def transformToRef(self,ref_wcs,force=False):
        """ Transform sky coords from ALL chips into X,Y coords in reference WCS.
        """
        if not isinstance(ref_wcs, pywcs.WCS):
            print(textutil.textbox('Reference WCS not a valid HSTWCS object'),
                  file=sys.stderr)
            raise ValueError
        # Need to concatenate catalogs from each input
        if self.outxy is None or force:
            outxy = ref_wcs.wcs_world2pix(self.all_radec[0],self.all_radec[1],self.origin)
            # convert outxy list to a Nx2 array
            self.outxy = np.column_stack([outxy[0][:,np.newaxis],outxy[1][:,np.newaxis]])
            if self.pars['writecat']:
                catname = self.rootname+"_refxy_catalog.coo"
                self.write_outxy(catname)
                self.catalog_names['ref_xy'] = catname

    def sortSkyCatalog(self):
        """ Sort and clip the source catalog based on the flux range specified
        by the user. It keeps a copy of the original full list in order to
        support iteration.

        """
        if len(self.all_radec_orig[2].nonzero()[0]) == 0:
            warn_str = "Source catalog NOT trimmed by flux/mag. No fluxes read in for sources!"
            print('\nWARNING: ',warn_str,'\n')
            log.warning(warn_str)
            return
        clip_catalog = False
        clip_prefix = ''
        for k in sortKeys:
            for p in self.pars.keys():
                pindx = p.find(k)
                if pindx >= 0 and self.pars[p] is not None:
                    log.info('found a match for %s to %s'%(
                                str(p),str(self.pars[p])))
                    # find prefix (if any)
                    clip_prefix = p[:pindx].strip()
                    #Only clip the catalog if one of the keys is specified
                    # in the catalog parameters, not the source finding pars
                    if clip_prefix and 'units' not in p:
                        clip_catalog = True
                        break
            if clip_catalog:
                break
        all_radec = None

        if clip_catalog:

            # Start by clipping by any specified flux range
            if self.pars[clip_prefix+'maxflux'] is not None or \
                    self.pars[clip_prefix+'minflux'] is not None:
                clip_catalog = True
                if self.pars[clip_prefix+'minflux'] is not None:
                    fluxmin = self.pars[clip_prefix+'minflux']
                else:
                    fluxmin = self.all_radec[2].min()

                if self.pars[clip_prefix+'maxflux'] is not None:
                    fluxmax = self.pars[clip_prefix+'maxflux']
                else:
                    fluxmax = self.all_radec[2].max()

                # apply flux limit clipping
                minindx = self.all_radec_orig[2] >= fluxmin
                maxindx = self.all_radec_orig[2] <= fluxmax
                flux_indx = np.bitwise_and(minindx,maxindx)
                all_radec = []
                all_radec.append(self.all_radec_orig[0][flux_indx])
                all_radec.append(self.all_radec_orig[1][flux_indx])
                all_radec.append(self.all_radec_orig[2][flux_indx])
                all_radec.append(np.arange(len(self.all_radec_orig[0][flux_indx])))

            if clip_prefix+'nbright' in self.pars and \
                    self.pars[clip_prefix+'nbright'] is not None:
                clip_catalog = True

                nbright = self.pars[clip_prefix+'nbright']

                # pick out only the brightest 'nbright' sources
                if self.pars[clip_prefix+'fluxunits'] == 'mag':
                    nbslice = slice(None,nbright)
                else:
                    nbslice = slice(nbright,None)

                if all_radec is None:
                    # work on copy of all original data
                    all_radec = copy.deepcopy(self.all_radec_orig)
                # find indices of brightest
                nbright_indx = np.argsort(all_radec[2])[nbslice]
                self.all_radec[0] = all_radec[0][nbright_indx]
                self.all_radec[1] = all_radec[1][nbright_indx]
                self.all_radec[2] = all_radec[2][nbright_indx]
                self.all_radec[3] = np.arange(len(all_radec[0][nbright_indx]))

            else:
                if all_radec is not None:
                    self.all_radec = copy.deepcopy(all_radec)

    def match(self,refimage, quiet_identity, **kwargs):
        """ Uses xyxymatch to cross-match sources between this catalog and
            a reference catalog (refCatalog).
        """
        ref_outxy = refimage.outxy
        refWCS = refimage.wcs
        refname = refimage.name
        ref_inxy = refimage.xy_catalog

        cat_src_type = kwargs['cat_src_type']
        del kwargs['cat_src_type']

        if not quiet_identity:
            print("Matching sources from \'{}\' with sources from "
                  "reference {} \'{}\'"
                  .format(self.name, cat_src_type, refname))
        #self.sortSkyCatalog() # apply any catalog sorting specified by the user
        self.transformToRef(refWCS)
        self.refWCS = refWCS
        # extract xyxymatch parameters from input parameters
        matchpars = kwargs.copy()
        self.match_pars = matchpars
        minobj = matchpars['minobj'] # needed for later
        del matchpars['minobj'] # not needed in xyxymatch

        self.goodmatch = True

        # Check to see whether or not it is being matched to itself
        if refname.strip() == self.name.strip():
            self.identityfit = True
            if not quiet_identity:
                log.info('NO fit performed for reference image: %s\n'%self.name)
        else:
            # convert tolerance from units of arcseconds to pixels, as needed
            radius = matchpars['searchrad']
            if matchpars['searchunits'] == 'arcseconds':
                radius /= refWCS.pscale

            # Determine xyoff (X,Y offset) and tolerance to be used with xyxymatch
            if matchpars['use2dhist']:
                xsh, ysh, maxval, flux, zpmat, qual = _estimate_2dhist_shift(
                    self.outxy,
                    ref_outxy,
                    searchrad=radius
                )
                xyoff = (xsh, ysh)

                if matchpars['see2dplot']:
                    zpstd = max(10, flux // 5) if qual else 10
                    title_str = ("Histogram of offsets: Peak has {:d} matches "
                                 "at ({:0.4g}, {:0.4g})"
                                 .format(maxval, xsh, ysh))
                    hist_name = None
                    if not self.interactive:
                        hist_name = 'hist2d_{0}.png'.format(self.rootname)

                    plot_pars = {
                        'data': zpmat,
                        'figure_id': self.figure_id,
                        'vmax': zpstd,
                        'xp': xsh,
                        'yp': ysh,
                        'searchrad': radius,
                        'title_str': title_str,
                        'plotname': hist_name,
                        'interactive': self.interactive
                    }

                    tweakutils.plot_zeropoint(plot_pars)

                if matchpars['see2dplot'] and ('residplot' in matchpars and
                                               'No' in matchpars['residplot']):
                    if self.interactive:
                        prompt = ("Press ENTER for next image, \n" +
                                  "     'n' to continue without updating header or \n" +
                                  "     'q' to quit immediately...\n")

                        if sys.version_info[0] >= 3:
                            a = input(prompt)
                        else:
                            a = raw_input(prompt)

                    else:
                        a = ' '

                    if 'n' in a.lower():
                        self.perform_update = False

                    if 'q' in a.lower():
                        self.quit_immediately = True

                if matchpars['see2dplot']:
                    self.figure_id += 1

            else:
                xoff = 0.
                yoff = 0.
                if not util.is_blank(matchpars['xoffset']):
                    xoff = matchpars['xoffset']
                if not util.is_blank(matchpars['yoffset']):
                    yoff = matchpars['yoffset']
                xyoff = (xoff, yoff)

            matches = xyxymatch(self.outxy, ref_outxy, origin=xyoff,
                                tolerance=matchpars['tolerance'],
                                separation=matchpars['separation'])

            if len(matches) > minobj:
                self.matches['image'] = np.column_stack([matches['input_x'][:,
                                np.newaxis],matches['input_y'][:,np.newaxis]])
                self.matches['ref'] = np.column_stack([matches['ref_x'][:,
                                np.newaxis],matches['ref_y'][:,np.newaxis]])
                self.matches['ref_idx'] = matches['ref_idx']
                self.matches['img_idx'] = self.all_radec[3][matches['input_idx']]
                self.matches['input_idx'] = matches['input_idx']
                self.matches['img_RA'] = self.all_radec[0][matches['input_idx']]
                self.matches['img_DEC'] = self.all_radec[1][matches['input_idx']]

                self.matches['ref_orig_xy'] = np.column_stack([
                                    np.array(ref_inxy[0])[matches['ref_idx']][:,np.newaxis],
                                    np.array(ref_inxy[1])[matches['ref_idx']][:,np.newaxis]])

                self.matches['img_orig_xy'] = np.column_stack([
                    np.array(self.xy_catalog[0])[matches['input_idx']][:,np.newaxis],
                    np.array(self.xy_catalog[1])[matches['input_idx']][:,np.newaxis]])
                self.matches['src_origin'] = ref_inxy[-1][matches['ref_idx']]
                print('Found %d matches for %s...'%(len(matches),self.name))

                if self.pars['writecat']:
                    matchfile = open(self.catalog_names['match'],mode='w+')
                    matchfile.write('#Reference: %s\n'%refname)
                    matchfile.write('#Input: %s\n'%self.name)
                    title = '#Ref_X        Ref_Y        '
                    title += 'Input_X    Input_Y        '
                    title += 'Ref_X0    Ref_Y0        '
                    title += 'Input_X0    Input_Y0        '
                    title += 'Ref_ID    Input_ID        '
                    title += 'Ref_Source\n'
                    fmtstr = '%0.6f    %0.6f        '*4
                    fmtstr += '%d    %d    %s\n'
                    matchfile.write(title)
                    for i in range(len(matches['input_x'])):
                        linestr = fmtstr%\
                            (matches['ref_x'][i],matches['ref_y'][i],\
                             matches['input_x'][i],matches['input_y'][i],
                            self.matches['ref_orig_xy'][:,0][i],
                            self.matches['ref_orig_xy'][:,1][i],
                            self.matches['img_orig_xy'][:,0][i],
                            self.matches['img_orig_xy'][:,1][i],
                            matches['ref_idx'][i],matches['input_idx'][i],
                            self.matches['src_origin'][i])
                        matchfile.write(linestr)
                    matchfile.close()
            else:
                warnstr = textutil.textbox('WARNING: \n'+
                    'Not enough matches (< %d) found for input image: %s'%(minobj,self.name))
                for line in warnstr.split('\n'):
                    log.warning(line)
                print(warnstr)
                self.goodmatch = False

    def performFit(self,**kwargs):
        """ Perform a fit between the matched sources.

            Parameters
            ----------
            kwargs : dict
                Parameter necessary to perform the fit; namely, *fitgeometry*.

            Notes
            -----
            This task still needs to implement (eventually) interactive iteration of
                   the fit to remove outliers.
        """
        assert(self.refWCS is not None)
        pars = kwargs.copy()
        self.fit_pars = pars

        self.fit = {'offset':[0.0,0.0],'rot':0.0,'scale':[1.0],'rms':[0.0,0.0],
                    'rms_keys':{'RMS_RA':0.0,'RMS_DEC':0.0,'NMATCH':0},
                    'fit_matrix':[[1.0,0.0],[0.0,1.0]], 'src_origin':[None]}

        if not self.identityfit:
            if self.matches is not None and self.goodmatch:
                self.fit = linearfit.iter_fit_all(
                    self.matches['image'],self.matches['ref'],
                    self.matches['img_idx'],self.matches['ref_idx'],
                    xyorig=self.matches['img_orig_xy'],
                    uvorig=self.matches['ref_orig_xy'],
                    mode=pars['fitgeometry'],nclip=pars['nclip'],
                    sigma=pars['sigma'],minobj=pars['minobj'],
                    center=self.refWCS.wcs.crpix,
                    verbose=self.verbose)

                self.fit['rms_keys'] = self.compute_fit_rms()
                radec_fit = self.refWCS.all_pix2world(self.fit['fit_xy'],1)
                self.fit['fit_RA'] = radec_fit[:,0]
                self.fit['fit_DEC'] = radec_fit[:,1]
                self.fit['src_origin'] = self.matches['src_origin']

                print('Computed ',pars['fitgeometry'],' fit for ',self.name,': ')
                if pars['fitgeometry'] == 'shift':
                    print("XSH: {:.4f}  YSH: {:.4f}"
                          .format(self.fit['offset'][0],
                                  self.fit['offset'][1]))
                elif pars['fitgeometry'] == 'rscale' and self.fit['proper']:
                    print("XSH: {:.4f}  YSH: {:.4f}    ROT: {:.10g}    "
                          "SCALE: {:.6f}".format(
                              self.fit['offset'][0],
                              self.fit['offset'][1],
                              self.fit['rot'],
                              self.fit['scale'][0]))
                elif pars['fitgeometry'] == 'general' or \
                     (pars['fitgeometry'] == 'rscale' and not self.fit['proper']):
                    print("XSH: {:.4f}  YSH: {:.4f}    PROPER ROT: {:.10g}    "
                          "".format(
                              self.fit['offset'][0],
                              self.fit['offset'][1],
                              self.fit['rot']))
                    print("<ROT>: {:.10g}  SKEW: {:.10g}    ROT_X: {:.10g}  "
                          "ROT_Y: {:.10g}".format(
                              self.fit['rotxy'][2],
                              self.fit['skew'],
                              self.fit['rotxy'][0],
                              self.fit['rotxy'][1]))
                    print("<SCALE>: {:.10g}  SCALE_X: {:.10g}  "
                          "SCALE_Y: {:.10g}".format(
                              self.fit['scale'][0],
                              self.fit['scale'][1],
                              self.fit['scale'][2]))
                else:
                    assert(False)

                print('FIT XRMS: {:<7.2g}    FIT YRMS: {:<7.2g}'
                      .format(*self.fit['rms']))
                print('FIT RMSE: {:<7.2g}    FIT MAE: {:<7.2g}\n'
                      .format(self.fit['rmse'], self.fit['mae']))
                print('RMS_RA: %.2g (deg)   RMS_DEC: %.2g (deg)\n'%(
                        self.fit['rms_keys']['RMS_RA'],
                        self.fit['rms_keys']['RMS_DEC']))
                print('Final solution based on ',self.fit['rms_keys']['NMATCH'],' objects.')

                self.write_fit_catalog()

                # Plot residuals, if requested by the user
                if 'residplot' in pars and "No" not in pars['residplot']:
                    xy = self.fit['img_coords']
                    resids = self.fit['resids']
                    xy_fit = xy + resids
                    title_str = 'Residuals\ for\ {0}\ using\ {1:6d}\ sources'.format(
                        self.name.replace('_','\_'),self.fit['rms_keys']['NMATCH'])
                    if not self.interactive:
                        resid_name = 'residuals_{0}.png'.format(self.rootname)
                        vector_name = resid_name.replace('residuals','vector')
                    else:
                        resid_name = None
                        vector_name = None
                    if pars['residplot'] == 'both':
                        tweakutils.make_vector_plot(None,
                            data=[xy[:,0],xy[:,1],xy_fit[:,0],xy_fit[:,1]],
                            figure_id=self.figure_id, vector=True,
                            labelsize=pars['labelsize'],
                            plotname=vector_name, title=title_str)
                        ptype=False # Setup
                        self.figure_id += 1
                    elif pars['residplot'] == 'vector':
                        ptype = True
                    else:
                        ptype = False

                    # Generate new plot
                    tweakutils.make_vector_plot(None,
                        data=[xy[:,0],xy[:,1],xy_fit[:,0],xy_fit[:,1]],
                        figure_id=self.figure_id, vector=ptype,
                        ylimit=pars['ylimit'], labelsize=pars['labelsize'],
                        plotname=resid_name, title=title_str)
                    if self.interactive:
                        prompt = ("Press ENTER for next image, \n" +
                                  "      'n' to continue without updating header or \n" +
                                  "      'q' to quit immediately...\n")
                        if sys.version_info[0] >= 3:
                            a = input(prompt)
                        else:
                            a = raw_input(prompt)
                    else:
                        a = ' '
                    if 'n' in a.lower():
                        self.perform_update = False
                    if 'q' in a.lower():
                        self.quit_immediately = True
            else:
                self.fit['offset'] = [np.nan,np.nan]
                self.fit['rot'] = np.nan
                self.fit['scale'] = [np.nan]

    def compute_fit_rms(self):
        # start by interpreting the fit to get the RMS values
        if not self.identityfit and self.goodmatch:
            crpix = self.refWCS.wcs.crpix + self.fit['rms']
            crval_rms = self.refWCS.wcs_pix2world([crpix],1)[0]
            rms_ra,rms_dec = np.abs(crval_rms - self.refWCS.wcs.crval)
            nmatch = self.fit['resids'].shape[0]
        else:
            rms_ra = 0.0
            rms_dec = 0.0
            nmatch = 0
        return {'RMS_RA':rms_ra,'RMS_DEC':rms_dec,'NMATCH':nmatch}

    def updateHeader(self, wcsname=None, reusename=False):
        """ Update header of image with shifts computed by *perform_fit()*.
        """
        # Insure filehandle is open and available...
        self.openFile()

        verbose_level = 1
        if not self.perform_update:
            verbose_level = 0
        # Create WCSCORR table to keep track of WCS revisions anyway
        if self.perform_update:
            wcscorr.init_wcscorr(self._im.hdu)

        extlist = []
        wcscorr_extname = self.ext_name
        if self.ext_name == "PRIMARY":
            extlist = [0]
        else:
            for ext in range(1,self.nvers+1):
                extlist.append((self.ext_name,ext))
                # add WCSNAME to SCI headers, if not provided (such as for
                # drizzled images directly obtained from the archive pre-AD)
                if ('wcsname' not in self._im.hdu[self.ext_name,ext].header and
                    self._im.hdu.fileinfo(0)['filemode'] == 'update'):
                    self._im.hdu[self.ext_name,ext].header['wcsname'] = 'Default'

        if not self.identityfit and self.goodmatch and \
                self.fit['offset'][0] != np.nan:
            updatehdr.updatewcs_with_shift(self._im.hdu, self.refWCS,
                wcsname=wcsname, reusename=reusename,
                fitgeom=self.fit_pars['fitgeometry'],
                xsh=self.fit['offset'][0],ysh=self.fit['offset'][1],
                rot=self.fit['rot'],scale=self.fit['scale'][0],
                fit=self.fit['fit_matrix'], verbose=verbose_level,
                xrms=self.fit['rms_keys']['RMS_RA'],
                yrms=self.fit['rms_keys']['RMS_DEC'])

            wnames = altwcs.wcsnames(self._im.hdu,ext=extlist[0])

            altkeys = []
            for k in wnames:
                if wnames[k] == wcsname:
                    altkeys.append(k)
            if len(altkeys) > 1 and ' ' in altkeys:
                altkeys.remove(' ')
            if len(altkeys) == 0:
                next_key = ' '
            else:
                next_key = altkeys[-1]
            if self.perform_update:
                log.info('    Writing out new WCS to alternate WCS: "%s"'%next_key)

            self.next_key = next_key
        else: #if self.identityfit or not self.goodmatch:
            if reusename:
                # Look for key of WCS with this name
                next_key = altwcs.getKeyFromName(self._im.hdu[extlist[0]].header,wcsname)
                # This wcsname is new, so start fresh
                if next_key is None:
                    next_key = altwcs.next_wcskey(self._im.hdu[extlist[0]].header)
            else:
                # Find key for next WCS and save again to replicate an updated solution
                next_key = altwcs.next_wcskey(self._im.hdu[extlist[0]].header)

            if self.perform_update:
                # archive current WCS as alternate WCS with specified WCSNAME
                # Start by archiving original PRIMARY WCS
                wnames = altwcs.wcsnames(self._im.hdu,ext=extlist[0])

                # Define a default WCSNAME in the case that the file to be
                # updated did not have the WCSNAME keyword defined already
                # (as will happen when updating images that have not been
                #  updated using updatewcs).
                if len(wnames) == 0:
                    pri_wcsname = None
                else:
                    # Safeguard against headers not having WCSNAME defined
                    # This would occur if they were written out by something
                    # other than stwcs.updatewcs v
                    if ' ' not in wnames:
                        self._im.hdu[extlist[0]].header['wscname'] = ''
                        wnames[' '] = ''
                    pri_wcsname = wnames[' ']

                next_pkey = altwcs.getKeyFromName(fits.getheader(self.name, extlist[0], memmap=False),pri_wcsname)
                log.info('    Saving Primary WCS to alternate WCS: "%s"'%next_pkey)

                altwcs.archiveWCS(self._im.hdu, extlist,
                                    wcskey=next_pkey, wcsname=pri_wcsname,
                                    reusekey=True)
                if reusename:
                    # Look for key of WCS with this name
                    next_key = altwcs.getKeyFromName(self._im.hdu[extlist[0]].header,wcsname)
                    # This wcsname is new, so start fresh
                    if next_key is None:
                        next_key = altwcs.next_wcskey(self._im.hdu[extlist[0]].header)
                else:
                    # Find key for next WCS and save again to replicate an updated solution
                    next_key = altwcs.next_wcskey(self._im.hdu[extlist[0]].header)
                    # update WCSNAME to be the new name
                    for ext in extlist:
                        self._im.hdu[ext].header['WCSNAME'] = wcsname

                # save again using new WCSNAME
                altwcs.archiveWCS(self._im.hdu, extlist,
                    wcskey=next_key,wcsname=wcsname, reusekey=reusename)
            self.next_key = ' '

        # add FIT values to image's PRIMARY header
        fimg = self._im.hdu

        if wcsname in ['',' ',None,"INDEF"]:
            wcsname = 'TWEAK'
        # Record values for the fit with both the PRIMARY WCS being updated
        # and the alternate WCS which will be created.
        assert(not self._im.closed)

        for ext in extlist:
            self._im.hdu[ext].header['FITNAME'+next_key] = wcsname
            for kw in self.fit['rms_keys']:
                self._im.hdu[ext].header.set(kw+next_key,
                                     self.fit['rms_keys'][kw],
                                     after='FITNAME'+next_key)

        if self.perform_update:
            log.info('Updating WCSCORR table with new WCS solution "%s"'%wcsname)
            wcscorr.update_wcscorr(self._im.hdu, wcs_id=wcsname,
                                   extname=self.ext_name)

    def writeHeaderlet(self,**kwargs):
        """ Write and/or attach a headerlet based on update to PRIMARY WCS
        """
        # Insure filehandle is open and available...
        self.openFile()

        pars = kwargs.copy()
        rms_pars = self.fit['rms_keys']

        str_kw = ['descrip','history','author','hdrfile']
        for kw in str_kw:
            if pars[kw] == '': pars[kw] = None

        # Call function with properly interpreted input parameters
        # Syntax: write_headerlet(filename, hdrname, output, sciext='SCI',
        #                    wcsname=None, wcskey=None, destim=None,
        #                    sipname=None, npolfile=None, d2imfile=None,
        #                    author=None, descrip=None, history=None,
        #                    rms_ra=None, rms_dec=None, nmatch=None, catalog=None,
        #                    attach=True, clobber=False):
        headerlet.write_headerlet(self._im.hdu, pars['hdrname'],
                output=pars['hdrfile'],
                wcsname=None, wcskey=self.next_key, destim=None,
                sipname=None, npolfile=None, d2imfile=None,
                author=pars['author'], descrip=pars['descrip'],
                history=pars['history'],
                nmatch=rms_pars['NMATCH'],catalog=pars['catalog'],
                attach=pars['attach'], clobber=pars['clobber']
            )

    def write_skycatalog(self,filename):
        """ Write out the all_radec catalog for this image to a file.
        """
        if self.all_radec is None:
            return
        ralist = self.all_radec[0]#.tolist()
        declist = self.all_radec[1]#.tolist()
        f = open(filename,'w')
        f.write("#Sky positions for: "+self.name+'\n')
        f.write("#RA        Dec\n")
        f.write("#(deg)     (deg)\n")
        for i in range(len(ralist)):
            f.write('%0.12f  %0.12f\n'%(ralist[i],declist[i]))
        f.close()

    def get_xy_catnames(self):
        """ Return a string with the names of input_xy catalog names
        """
        catstr = self.name+'  '
        if 'input_xy' in self.catalog_names:
            for xycat in self.catalog_names['input_xy']:
                catstr += '  '+xycat
        return catstr + '\n'

    def write_fit_catalog(self):
        """ Write out the catalog of all sources and resids used in the final fit.
        """
        if self.pars['writecat']:
            log.info('Creating catalog for the fit: {:s}'.format(self.catalog_names['fitmatch']))
            f = open(self.catalog_names['fitmatch'],'w')
            f.write('# Input image: {:s}\n'.format(self.filename))
            f.write('# Coordinate mapping parameters: \n')
            f.write('#    X and Y rms: {:.2g}  {:.2g}\n'
                    .format(self.fit['rms'][0], self.fit['rms'][1]))
            f.write('#    X and Y shift: {:.4f}  {:.4f}\n'
                    .format(self.fit['offset'][0], self.fit['offset'][1]))
            f.write('#    X and Y scale: {:.10g}  {:.10g}\n'
                    .format(self.fit['scale'][1], self.fit['scale'][2]))
            f.write('#    X and Y rotation: {:.10g}  {:.10g}\n'
                    .format(self.fit['rotxy'][0], self.fit['rotxy'][1]))
            f.write('#    <rotation> and skew: {:.10g}  {:.10g}\n'
                    .format(self.fit['rotxy'][2], self.fit['skew']))

            f.write('# \n# Input Coordinate Listing\n')
            f.write('#     Column 1: X (reference)\n')
            f.write('#     Column 2: Y (reference)\n')
            f.write('#     Column 3: X (input)\n')
            f.write('#     Column 4: Y (input)\n')
            f.write('#     Column 5: X (fit)\n')
            f.write('#     Column 6: Y (fit)\n')
            f.write('#     Column 7: X (residual)\n')
            f.write('#     Column 8: Y (residual)\n')
            f.write('#     Column 9: Original X (reference)\n')
            f.write('#     Column 10: Original Y (reference)\n')
            f.write('#     Column 11: Original X (input)\n')
            f.write('#     Column 12: Original Y (input)\n')
            f.write('#     Column 13: Ref ID\n')
            f.write('#     Column 14: Input ID\n')
            f.write('#     Column 15: Input EXTVER ID \n')
            f.write('#     Column 16: RA (fit)\n')
            f.write('#     Column 17: Dec (fit)\n')
            f.write('#     Column 18: Ref source provenience\n')
            #
            # Need to add chip ID for each matched source to the fitmatch file
            # The chip information can be extracted from the following source:
            #
            #     self.chip_catalogs[sci_extn] = {'catalog':catalog,'wcs':wcs}
            #     xypos = catalog.xypos
            img_chip_id = self.fit['img_indx'].copy()
            for sci_extn in range(1,self.nvers+1):
                catalog = self.chip_catalogs[sci_extn]['catalog']
                if catalog.xypos is not None:
                    img_indx_orig = self.chip_catalogs[sci_extn]['catalog'].xypos[3]
                    chip_min = img_indx_orig.min()
                    chip_max = img_indx_orig.max()
                    cid = np.logical_and((img_chip_id >= chip_min),(img_chip_id <= chip_max))
                    img_chip_id[cid] = sci_extn
            #
            f.write('#\n')
            f.close()

            xydata = [[self.fit['ref_coords'][:,0],self.fit['ref_coords'][:,1],
                      self.fit['img_coords'][:,0],self.fit['img_coords'][:,1],
                      self.fit['fit_xy'][:,0],self.fit['fit_xy'][:,1],
                      self.fit['resids'][:,0],self.fit['resids'][:,1],
                      self.fit['ref_orig_xy'][:,0],
                      self.fit['ref_orig_xy'][:,1],
                      self.fit['img_orig_xy'][:,0],
                      self.fit['img_orig_xy'][:,1]],
                      [self.fit['ref_indx'],self.fit['img_indx'],img_chip_id],
                      [self.fit['fit_RA'],self.fit['fit_DEC']],
                      [self.fit['src_origin']]
                    ]

            tweakutils.write_xy_file(self.catalog_names['fitmatch'],xydata,
                append=True,format=["%15.6f","%8d","%20.12f","   %s"])

    def write_outxy(self,filename):
        """ Write out the output(transformed) XY catalog for this image to a file.
        """
        f = open(filename,'w')
        f.write("#Pixel positions for: "+self.name+'\n')
        f.write("#X           Y\n")
        f.write("#(pix)       (pix)\n")
        for i in range(self.all_radec[0].shape[0]):
            f.write('%f  %f\n'%(self.outxy[i,0],self.outxy[i,1]))
        f.close()

    def get_shiftfile_row(self):
        """ Return the information for a shiftfile for this image to provide
            compatability with the IRAF-based MultiDrizzle.
        """
        if self.fit is not None:
            rowstr = '%s    %0.6f  %0.6f    %0.6f     %0.6f   %0.6f  %0.6f\n'%(
                    self.name,self.fit['offset'][0],self.fit['offset'][1],
                    self.fit['rot'],self.fit['scale'][0],
                    self.fit['rms'][0],self.fit['rms'][1])
        else:
            rowstr = None
        return rowstr

    def clean(self):
        """ Remove intermediate files created.
        """
        #TODO: add cleaning of mask files, *if* created ...
        for f in self.catalog_names:
            if 'match' in f:
                if os.path.exists(self.catalog_names[f]):
                    log.info('Deleting intermediate match file: %s'%
                                self.catalog_names[f])
                    os.remove(self.catalog_names[f])
            else:
                for extn in f:
                    if os.path.exists(extn):
                        log.info('Deleting intermediate catalog: %d'%extn)
                        os.remove(extn)


class RefImage:
    """ This class provides all the information needed by to define a reference
    tangent plane and list of source positions on the sky.

    .. warning::
        When ``wcs_list`` is a Python list of ``WCS`` objects,
        each element must be an instance of `stwcs.wcsutil.HSTWCS`.

    """
    def __init__(self, wcs_list, catalog, xycatalog=None, cat_origin=None, **kwargs):
        assert(isinstance(xycatalog, list) if xycatalog is not None else True)
        hdulist = None

        self.pars = kwargs

        if isinstance(wcs_list, str):
            # Input was a filename for the reference image
            #froot, fextn = fu.parseFilename(wcs_list)
            fi = parse_cs_line(wcs_list, fnamesOnly=False, doNotOpenDQ=True,
                               im_fmode='readonly')
            fi[0].release_all_images()
            if len(fi) != 1:
                ValueError('Reference image file name must contain a single '
                           'file name specification.')
            froot = fi[0].image
            if len(fi[0].fext) == 0:
                # there are no 'SCI' extensions in the input file -> use
                # primary HDU:
                fextn = 0
            else:
                fextn = fi[0].fext[0] # <- we assume that only one
                                      # extension is specified or that the
                                      # first one should be used

            hdulist = fits.open(froot, mode='readonly', memmap=False)
            try:
                self.wcs = stwcs.wcsutil.HSTWCS(hdulist, fextn)
                if _is_wcs_distorted(self.wcs):
                    if self.wcs.instrument == 'DEFAULT':
                        raise ValueError("Distorted non-HST reference images "
                                         "are not supported.")
                    log.warn("\nReference image contains a distorted WCS.\n"
                             "Using the undistorted version of this WCS.\n")
                    self.wcs = utils.output_wcs([self.wcs], undistort=True)
            except KeyError as e:
                if not e.args[0].startswith('Unsupported instrument'):
                    raise e
                # try using astropy.wcs:
                self.wcs = pywcs.WCS(hdulist[fextn].header, hdulist)
                if _is_wcs_distorted(self.wcs):
                    raise ValueError("Distorted non-HST reference images "
                                     "are not supported.")
            finally:
                if hdulist is not None:
                    hdulist.close()

            self.wcs.filename = froot

        elif isinstance(wcs_list, list):
            # generate a reference tangent plane from a list of STWCS objects
            undistort = _is_wcs_distorted(wcs_list[0])
            if undistort and wcs_list[0].instrument == 'DEFAULT':
                raise ValueError("Distorted non-HST reference images "
                                 "are not supported.")

            self.wcs = utils.output_wcs(wcs_list, undistort=undistort)
            self.wcs.cd = wcs_functions.make_perfect_cd(self.wcs)

            if 'ref_wcs_name' in kwargs:
                self.wcs.filename = kwargs['ref_wcs_name']
            else:
                self.wcs.filename = wcs_list[0].filename

        else:
            # only a single WCS provided, so use that as the definition
            if isinstance(wcs_list, stwcs.wcsutil.HSTWCS):
                 # User provided full HSTWCS object
                self.wcs = wcs_list
                froot = self.wcs.filename
                if _is_wcs_distorted(self.wcs):
                    if self.wcs.instrument == 'DEFAULT':
                        raise ValueError("Distorted non-HST reference images "
                                         "are not supported.")
                    log.warn("\nReference image contains a distorted WCS.\n"
                             "Using the undistorted version of this WCS.\n")
                    self.wcs = utils.output_wcs([self.wcs], undistort=True)
                    self.wcs.filename = froot

            else:
                raise TypeError("Unsupported 'wcs_list' type.")

        # assert(not _is_wcs_distorted(self.wcs))
        if _is_wcs_distorted(self.wcs):
            warnstr = textutil.textbox(
                "WARNING: Reference image appears to have a distorted WCS. \n"
                "This is likely due to a deviation of WCS' CD (or PC) "
                "matrix from orthogonality. \n"
                "Residual distortions in the reference WCS will result "
                "in sub-optimal image alignment.\n"
                "Check orthogonality of the CD (or PC) matrix!"
            )
            for line in warnstr.split('\n'):
                log.warning(line)
            print(warnstr)

        self.dirty = False

        if 'use_sharp_round' in kwargs:
            self.use_sharp_round = kwargs['use_sharp_round']
        else:
            self.use_sharp_round = False

        self.name = self.wcs.filename
        self.refWCS = None

        # See if computation of the bounding polygon for
        # the reference catalog was requested:
        find_bounding_polygon = kwargs.pop('find_bounding_polygon', True)

        # Interpret the provided catalog
        self.catalog = catalogs.RefCatalog(None, catalog,**kwargs)
        self.catalog.buildCatalogs()

        self.all_radec = self.catalog.radec
        nradec = self.all_radec[0].shape[0]
        nradec_col = len(self.all_radec)
        if nradec_col == 2:
            self.all_radec.append(np.zeros(nradec, dtype=float))
            self.all_radec.append(np.arange(nradec, dtype=int))
            nradec_col += 2
        elif nradec_col == 3:
            self.all_radec.append(np.arange(nradec, dtype=int))
            nradec_col += 1

        # add input source positions to RefCatalog for reporting in final
        # matched output catalog files
        if xycatalog is not None:
            xypostypes = 3*[float] + [int] + \
                (3 if self.use_sharp_round else 0)*[float] + [object]
            self.xy_catalog = [np.asarray(cat, dtype=dt) for dt,cat in \
                               zip(xypostypes,xycatalog)]
            ncol = len(self.xy_catalog)
            nobj = self.xy_catalog[0].shape[0]
            assert(nradec == nobj)

            assert(ncol >= 2)
            if ncol == 2: # account for flux column
                self.xy_catalog.append(np.zeros(nobj, dtype=float))
                ncol = 3

            if ncol == 3: # add source ID column
                self.xy_catalog.append(np.arange(nobj)+self.start_id)
                ncol = 4

            if self.use_sharp_round:
                # if necessary, add sharpness and roundness columns:
                for i in range(ncol, 7):
                    self.xy_catalog.append(np.zeros(nobj, dtype=float))
                ncol = len(self.xy_catalog)

        else:
            # Convert RA/Dec positions of source from refimage into
            # X,Y positions based on WCS of refimage
            self.xy_catalog = []
            if 'refxyunits' in self.pars and self.pars['refxyunits'] == 'pixels':
                xypos = [self.all_radec[0], self.all_radec[1]]
            else:
                xypos = self.wcs.wcs_world2pix(self.all_radec[0], self.all_radec[1],1)
            nobj = xypos[0].shape[0]
            assert(nradec == nobj)
            self.xy_catalog.append(xypos[0])
            self.xy_catalog.append(xypos[1])
            self.xy_catalog.append(copy.copy(self.all_radec[2]))
            self.xy_catalog.append(copy.copy(self.all_radec[3]))
            if self.use_sharp_round:
                for i in [5,6,7]:
                    self.xy_catalog.append(np.zeros(nobj, dtype=float))
                ncol = 7
            else:
                ncol = 4

        # add source provenience (origin)
        if (self.use_sharp_round and ncol == 7) or \
           (not self.use_sharp_round and ncol == 4):
            if cat_origin is None:
                self.xy_catalog.append(np.asarray(nobj*[self.name],
                                                  dtype=object))
            else:
                if isinstance(cat_origin, str):
                    self.xy_catalog.append(np.asarray(nobj*[cat_origin],
                                                      dtype=object))
                elif isinstance(cat_origin, list) or \
                     isinstance(cat_origin, np.ndarray):
                    self.xy_catalog.append(np.asarray(cat_origin,
                                                      dtype=object))
                    orig_len = len(cat_origin)
                    if orig_len < nobj:
                        self.xy_catalog[-1] = np.append(self.xy_catalog[-1],
                            np.asarray(nobj*[''], dtype=object))
                else:
                    raise TypeError("Parameter 'cat_origin' must be "
                        "a string, a list, or a numpy.ndarray")

        self.outxy = None
        self.origin = 1
        if self.all_radec is not None:
            # convert sky positions to X,Y positions on reference tangent plane
            self.transformToRef()

        # Compute bounding convex hull for the reference catalog:
        if (find_bounding_polygon or IMAGE_USE_CONVEX_HULL) and \
           self.outxy is not None:
            xy_vertices = np.asarray(convex_hull(list(map(tuple,self.outxy))),
                                     dtype=np.float64)
            if xy_vertices.shape[0] > 2:
                rdv = self.wcs.wcs_pix2world(xy_vertices, 1)
                self.skyline = SphericalPolygon.from_radec(rdv[:,0], rdv[:,1])
            else:
                self.skyline = SphericalPolygon([])

            if IMGCLASSES_DEBUG:
                _debug_write_region_fk5('dbg_tweakreg_refcat_bounding_polygon.reg',
                    list(zip(*rdv.transpose())), list(zip(*self.all_radec)),
                    self.xy_catalog[-1])
        else:
            self.skyline = SphericalPolygon([])

    def clear_dirty_flag(self):
        self.dirty = False

    def set_dirty(self):
        self.dirty = True

    def write_skycatalog(self, filename, show_flux=False, show_id=False):
        """ Write out the all_radec catalog for this image to a file.
        """
        f = open(filename,'w')
        f.write("#Sky positions for cumulative reference catalog. Initial catalog from: "+self.name+'\n')
        header1 = "#RA        Dec"
        header2 = "#(deg)     (deg)"
        if show_flux:
            header1 += "        Flux"
            header2 += "       (counts)"
        if show_id:
            header1 += "        ID        Origin"
            header2 += ""
        header1 += "\n"
        header2 += "\n"

        f.write(header1)
        f.write(header2)

        show_details = show_flux or show_id
        flux_end_char = ''
        if show_details and show_flux:
            if show_id:
                flux_end_char = '\t'

        for i in range(self.all_radec[0].shape[0]):
            src_line = "{:.7f}  {:.7f}" \
                .format(self.all_radec[0][i], self.all_radec[1][i])
            if show_details:
                #src_line += " #"
                if show_flux:
                    #src_line += " flux: {:.5g}{:s}" \
                        #.format(self.xy_catalog[2][i], flux_end_char)
                    src_line += "  {:.5g}".format(self.xy_catalog[2][i])
                if show_id:
                    #src_line += " ID: {:d}\torigin: '{:s}'" \
                        #.format(self.xy_catalog[3][i], self.xy_catalog[-1][i])
                    src_line += "  {:d}  {:s}".format(self.xy_catalog[3][i],
                                                      self.xy_catalog[-1][i])
            f.write(src_line + '\n')
        f.close()

    def transformToRef(self):
        """ Transform reference catalog sky positions (self.all_radec)
        to reference tangent plane (self.wcs) to create output X,Y positions.
        """
        if 'refxyunits' in self.pars and self.pars['refxyunits'] == 'pixels':
            log.info('Creating RA/Dec positions for reference sources...')
            self.outxy = np.column_stack([self.all_radec[0][:,np.newaxis],self.all_radec[1][:,np.newaxis]])
            skypos = self.wcs.wcs_pix2world(self.all_radec[0],self.all_radec[1],self.origin)
            self.all_radec[0] = skypos[0]
            self.all_radec[1] = skypos[1]
        else:
            log.info('Converting RA/Dec positions of reference sources from "%s" to '%self.name+
                        'X,Y positions in reference WCS...')
            self.refWCS = self.wcs
            outxy = self.wcs.wcs_world2pix(self.all_radec[0],self.all_radec[1],self.origin)
            # convert outxy list to a Nx2 array
            self.outxy = np.column_stack([outxy[0][:,np.newaxis],outxy[1][:,np.newaxis]])

    def append_not_matched_sources(self, image):
        assert(hasattr(image, 'fit') and hasattr(image, 'matches'))
        if not image.goodmatch or image.identityfit:
            return
        matched_idx = image.matches['input_idx']

        # append new sources to the reference catalog:
        old_ref_nsources = self.outxy.shape[0]
        adding_nsources = image.outxy.shape[0] - matched_idx.shape[0]
        if adding_nsources < 1:
            return
        print("Adding {} new sources to the reference catalog "
              "for a total of {} sources."
              .format(adding_nsources, old_ref_nsources + adding_nsources))
        not_matched_mask = np.ones(image.outxy.shape[0], dtype=np.bool)
        not_matched_mask[matched_idx] = False
        not_matched_outxy = image.outxy[not_matched_mask]

        # apply corrections based on the fit:
        new_outxy = np.dot(not_matched_outxy-image.fit['offset']-self.wcs.wcs.crpix,
                           image.fit['fit_matrix'].transpose())+self.wcs.wcs.crpix

        # convert to RA & DEC:
        new_radec = self.wcs.wcs_pix2world(new_outxy, 1)

        self.outxy = np.append(self.outxy, new_outxy, 0)
        id1 = self.all_radec[3][-1] + 1
        id2 = id1 + len(self.all_radec[3])
        self.all_radec = [np.append(self.all_radec[0], new_radec[:,0], 0),
                          np.append(self.all_radec[1], new_radec[:,1], 0),
                          np.append(self.all_radec[2], image.all_radec[2][not_matched_mask], 0),
                          np.append(self.all_radec[3], np.arange(id1,id2), 0)]

        # Append original image coordinates:
        self.xy_catalog[0] = np.append(self.xy_catalog[0],
            np.asarray(image.xy_catalog[0])[not_matched_mask], 0)
        self.xy_catalog[1] = np.append(self.xy_catalog[1],
            np.asarray(image.xy_catalog[1])[not_matched_mask], 0)
        ncol = len(self.xy_catalog)
        self.use_sharp_round
        for i,col in enumerate(image.xy_catalog[2:-1]):
            i2 = i + 2
            if i2 < ncol-1:
                self.xy_catalog[i2] = np.append(self.xy_catalog[i2], col, 0)
        self.xy_catalog[-1] = np.append(self.xy_catalog[-1],
                                        np.asarray(image.xy_catalog[-1])[not_matched_mask], 0)

        #self.skyline = self.skyline.union(skyline)
        xy_vertices = np.asarray(convex_hull(list(map(tuple,self.outxy))),
                                 dtype=np.float64)
        rdv = self.wcs.wcs_pix2world(xy_vertices, 1)
        self.skyline = SphericalPolygon.from_radec(rdv[:,0], rdv[:,1])
        if IMGCLASSES_DEBUG:
            _debug_write_region_fk5('dbg_tweakreg_refcat_bounding_polygon.reg',
                                    list(zip(*rdv.transpose())),
                                    list(zip(*self.all_radec)),
                                    self.xy_catalog[-1])

        self.set_dirty()

    def get_shiftfile_row(self):
        """ Return the information for a shiftfile for this image to provide
            compatability with the IRAF-based MultiDrizzle.
        """
        rowstr = '%s    0.0  0.0    0.0     1.0    0.0 0.0\n'%(self.name)
        return rowstr

    def clean(self):
        """ Remove intermediate files created
        """
        if not util.is_blank(self.catalog.catname) and os.path.exists(self.catalog.catname):
            os.remove(self.catalog.catname)

    def close(self):
        pass


def build_referenceWCS(catalog_list):
    """ Compute default reference WCS from list of Catalog objects.
    """
    wcslist = []
    for catalog in catalog_list:
        for scichip in catalog.catalogs:
            wcslist.append(catalog.catalogs[scichip]['wcs'])
    return utils.output_wcs(wcslist)


def verify_hdrname(**pars):
    if not pars['headerlet']:
        return

    # Insure a valid value for hdrname
    if util.is_blank(pars['hdrname']) and \
        not util.is_blank(pars['wcsname']):
            pars['hdrname'] = pars['wcsname']

    # interpret input strings
    if pars['hdrname'].strip() == '':
        warnstr = 'WARNING:\n'
        warnstr += '    Headerlet generation requested, but \n'
        warnstr += '    no valid "hdrname" parameter value was provided!\n'
        warnstr += '    The requested headerlets can NOT be generated! \n'
        warnstr += '\n'
        warnstr += '    Please provide valid "hdrname" or "wcsname" value \n'
        warnstr += '        in your next run to write out headerlets.\n'
        print(textutil.textbox(warnstr))
        return


def _debug_write_region_fk5(fname, ivert, points, src_origin=None, append=False):
    if append:
        fh = open(fname, 'a')
    else:
        fh = open(fname, 'w')
        fh.write('# Region file format: DS9 version 4.1\n')
        fh.write('global color=green dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        fh.write('fk5\n')
    fh.write('polygon({}) # color=red width=2\n'
             .format(','.join(map(str,[x for e in ivert for x in e]))))
    if src_origin is None:
        for p in points:
            fh.write('point({}) # point=x\n'
                     .format(','.join(map(str,[e for e in p[:2]]))))
    else:
        for p,orig in zip(points,src_origin):
            fh.write('point({}) # point=x text={{{}}}\n'
                     .format(','.join(map(str,[e for e in p[:2]])), orig))
    fh.close()


def convex_hull(points):
    """Computes the convex hull of a set of 2D points.

    Implements `Andrew's monotone chain algorithm <http://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain>`_.
    The algorithm has O(n log n) complexity.

    Credit: `<http://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain>`_

    Parameters
    ----------

    points : list of tuples
        An iterable sequence of (x, y) pairs representing the points.

    Returns
    -------
    Output : list
        A list of vertices of the convex hull in counter-clockwise order,
        starting from the vertex with the lexicographically smallest
        coordinates.
    """

    # Sort the points lexicographically (tuples are compared lexicographically).
    # Remove duplicates to detect the case we have just one unique point.
    points = sorted(set(points))

    # Boring case: no points or a single point, possibly repeated multiple times.
    if len(points) <= 1:
        return points

    # 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
    # Returns a positive value, if OAB makes a counter-clockwise turn,
    # negative for clockwise turn, and zero if the points are collinear.
    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    # Build lower hull
    lower = []
    for p in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    # Build upper hull
    upper = []
    for p in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    # Concatenation of the lower and upper hulls gives the convex hull.
    # Last point of each list is omitted because it is repeated at the beginning of the other list.
    return lower[:-1] + upper


def _is_wcs_distorted(wcs):
    return not (_is_cd_unitary(wcs) and wcs.sip is None and \
                wcs.cpdis1 is None and wcs.cpdis2 is None and \
                wcs.det2im1 is None and wcs.det2im2 is None)


def _is_cd_unitary(wcs):
    # set maximum acceptable deviation of matrix elements from zero -
    # a measure of closeness to matrix being unitary:
    maxerr = 50.0 * np.finfo(np.float32).eps

    cd = None
    if hasattr(wcs.wcs, 'cd'):
        cd = wcs.wcs.cd
    elif hasattr(wcs.wcs, 'pc'):
        if not hasattr(wcs.wcs, 'cdelt'):
            raise ValueError("Inconsistent WCS specification (missing CDELT "
                             "with a valid PC).")
        cd = np.dot(np.diag(wcs.wcs.cdelt), wcs.wcs.pc)

    if cd is None:
        return False

    shape = cd.shape
    assert(len(shape) == 2 and shape[0] == shape[1])

    pixarea = np.abs(np.linalg.det(cd))
    if pixarea <= 0.0:
        raise ValueError("Singular CD (or PC & CDELT) matrix.")

    # NOTE: Technically, below we should use np.dot(cd, np.conjugate(cd.T))
    # However, I am not aware of complex CD/PC matrices...
    I = np.dot(cd, cd.T) / pixarea
    cd_unitary_err = np.amax(np.abs(I - np.eye(shape[0])))
    # or, use RMS instead of max error:
    #cd_unitary_err = np.sqrt(np.mean(np.abs(I-np.eye(shape[0]))**2))

    return (cd_unitary_err < maxerr)


def _xy_2dhist(imgxy, refxy, r):
    # This code replaces the C version (arrxyzero) from carrutils.c
    # It is about 5-8 times slower than the C version.
    dx = np.subtract.outer(imgxy[:, 0], refxy[:, 0]).ravel()
    dy = np.subtract.outer(imgxy[:, 1], refxy[:, 1]).ravel()
    idx = np.where((dx < r + 0.5) & (dx >= -r - 0.5) &
                   (dy < r + 0.5) & (dy >= -r - 0.5))
    h = np.histogram2d(dx[idx], dy[idx], 2 * r + 1,
                       [[-r - 0.5, r + 0.5], [-r - 0.5, r + 0.5]])
    return h[0].T


def _estimate_2dhist_shift(imgxy, refxy, searchrad=3.0):
    """ Create a 2D matrix-histogram which contains the delta between each
        XY position and each UV position. Then estimate initial offset
        between catalogs.
    """
    print("Computing initial guess for X and Y shifts...")

    # create ZP matrix
    zpmat = _xy_2dhist(imgxy, refxy, r=searchrad)

    nonzeros = np.count_nonzero(zpmat)
    if nonzeros == 0:
        # no matches within search radius. Return (0, 0):
        print("WARNING: No matches found within a search radius of {:g} "
              "pixels.".format(searchrad))
        return 0.0, 0.0, 0, 0, zpmat, False

    elif nonzeros == 1:
        # only one non-zero bin:
        yp, xp = np.unravel_index(np.argmax(zpmat), zpmat.shape)
        maxval = int(np.ceil(zpmat[yp, xp]))
        xp -= searchrad
        yp -= searchrad
        print("Found initial X and Y shifts of {:.4g}, {:.4g} "
              "based on a single non-zero bin and {} matches"
              .format(xp, yp, maxval))
        return xp, yp, maxval, maxval, zpmat, True

    (xp, yp), fit_status, fit_sl = _find_peak(zpmat, peak_fit_box=5,
                                              mask=zpmat > 0)

    if fit_status.startswith('ERROR'):
        print("WARNING: No valid shift found within a search radius of {:g} "
              "pixels.".format(searchrad))
        maxval = int(np.ceil(zpmat.max()))
        return 0.0, 0.0, maxval, maxval, zpmat, False

    xp -= searchrad
    yp -= searchrad

    if fit_status == 'WARNING:EDGE':
        print(
            "WARNING: Found peak in the 2D histogram lies at the edge of "
            "the histogram. Try increasing 'searchrad' for improved results."
        )

    # Attempt to estimate "significance of detection":
    maxval = zpmat.max()
    zpmat_mask = (zpmat > 0) & (zpmat < maxval)

    if np.any(zpmat_mask):
        bkg = zpmat[zpmat_mask].mean()
        sig = maxval / np.sqrt(bkg)

    flux = int(zpmat[fit_sl].sum())
    print("Found initial X and Y shifts of {:.4g}, {:.4g} "
          "with significance of {:.4g} and {:d} matches"
          .format(xp, yp, sig, flux))

    return xp, yp, int(np.ceil(maxval)), flux, zpmat, True


def _find_peak(data, peak_fit_box=5, mask=None):
    """
    Find location of the peak in an array. This is done by fitting a second
    degree 2D polynomial to the data within a `peak_fit_box` and computing the
    location of its maximum. An initial
    estimate of the position of the maximum will be performed by searching
    for the location of the pixel/array element with the maximum value.

    Parameters
    ----------
    data : numpy.ndarray
        2D data.

    peak_fit_box : int, optional
        Size (in pixels) of the box around the initial estimate of the maximum
        to be used for quadratic fitting from which peak location is computed.
        It is assumed that fitting box is a square with sides of length
        given by ``peak_fit_box``.

    mask : numpy.ndarray, optional
        A boolean type `~numpy.ndarray` indicating "good" pixels in image data
        (`True`) and "bad" pixels (`False`). If not provided all pixels
        in `image_data` will be used for fitting.

    Returns
    -------
    coord : tuple of float
        A pair of coordinates of the peak.

    fit_status : str
        Status of the peak search. Currently the following values can be
        returned:

        - ``'SUCCESS'``: Fit was successful and peak is not on the edge of
          the input array;
        - ``'ERROR:NODATA'``: Not enough valid data to perform the fit; The
          returned coordinate is the center of input array;
        - ``'WARNING:EDGE'``: Peak lies on the edge of the input array.
          Returned coordinates are the result of a discreet search;
        - ``'WARNING:BADFIT'``: Performed fid did not find a maximum. Returned
          coordinates are the result of a discreet search;
        - ``'WARNING:CENTER-OF-MASS'``: Returned coordinates are the result
          of a center-of-mass estimate instead of a polynomial fit. This is
          either due to too few points to perform a fit or due to a
          failure of the polynomial fit.

    fit_box : a tuple of two tuples
        A tuple of two tuples of the form ``((x1, x2), (y1, y2))`` that
        indicates pixel ranges used for fitting (these indices can be used
        directly for slicing input data)

    """
    # check arguments:
    data = np.asarray(data, dtype=np.float64)
    ny, nx = data.shape

    # find index of the pixel having maximum value:
    if mask is None:
        jmax, imax = np.unravel_index(np.argmax(data), data.shape)
        coord = (float(imax), float(jmax))

    else:
        j, i = np.indices(data.shape)
        i = i[mask]
        j = j[mask]

        if i.size == 0:
            # no valid data:
            coord = ((nx - 1.0) / 2.0, (ny - 1.0) / 2.0)
            return coord, 'ERROR:NODATA', np.s_[0:ny-1, 0:nx-1]

        ind = np.argmax(data[mask])
        imax = i[ind]
        jmax = j[ind]
        coord = (float(imax), float(jmax))

    if data[jmax, imax] < 1:
        # no valid data: we need some counts in the histogram bins
        coord = ((nx - 1.0) / 2.0, (ny - 1.0) / 2.0)
        return coord, 'ERROR:NODATA', np.s_[0:ny-1, 0:nx-1]

    # choose a box around maxval pixel:
    x1 = max(0, imax - peak_fit_box // 2)
    x2 = min(nx, x1 + peak_fit_box)
    y1 = max(0, jmax - peak_fit_box // 2)
    y2 = min(ny, y1 + peak_fit_box)

    # if peak is at the edge of the box, return integer indices of the max:
    if imax == x1 or imax == x2 or jmax == y1 or jmax == y2:
        return (float(imax), float(jmax)), 'WARNING:EDGE', np.s_[y1:y2, x1:x2]

    # expand the box if needed:
    if (x2 - x1) < peak_fit_box:
        if x1 == 0:
            x2 = min(nx, x1 + peak_fit_box)
        if x2 == nx:
            x1 = max(0, x2 - peak_fit_box)

    if (y2 - y1) < peak_fit_box:
        if y1 == 0:
            y2 = min(ny, y1 + peak_fit_box)
        if y2 == ny:
            y1 = max(0, y2 - peak_fit_box)

    if x2 - x1 == 0 or y2 - y1 == 0:
        # not enough data:
        coord = ((nx - 1.0) / 2.0, (ny - 1.0) / 2.0)
        return coord, 'ERROR:NODATA', np.s_[y1:y2, x1:x2]

    # fit a 2D 2nd degree polynomial to data:
    xi = np.arange(x1, x2)
    yi = np.arange(y1, y2)
    x, y = np.meshgrid(xi, yi)
    x = x.ravel()
    y = y.ravel()
    v = np.vstack((np.ones_like(x), x, y, x*y, x*x, y*y)).T
    d = data[y1:y2, x1:x2].ravel()
    if mask is not None:
        m = mask[y1:y2, x1:x2].ravel()
        v = v[m]
        d = d[m]

    if d.size == 0 or np.max(d) <= 0:
        coord = ((nx - 1.0) / 2.0, (ny - 1.0) / 2.0)
        return coord, 'ERROR:NODATA', np.s_[y1:y2, x1:x2]

    if d.size < 6:
        # we need at least 6 points to fit a 2D quadratic polynomial
        # attempt center-of-mass instead:
        dt = d.sum()
        xc = np.dot(v[:, 1], d) / dt
        yc = np.dot(v[:, 2], d) / dt
        return (xc, yc), 'WARNING:CENTER-OF-MASS', np.s_[y1:y2, x1:x2]

    try:
        c = np.linalg.lstsq(v, d, rcond=None)[0]
    except np.linalg.LinAlgError:
        print("WARNING: Least squares failed!\n{}".format(c))

        # attempt center-of-mass instead:
        dt = d.sum()
        xc = np.dot(v[:, 1], d) / dt
        yc = np.dot(v[:, 2], d) / dt
        return (xc, yc), 'WARNING:CENTER-OF-MASS', np.s_[y1:y2, x1:x2]

    # find maximum of the polynomial:
    _, c10, c01, c11, c20, c02 = c
    det = 4 * c02 * c20 - c11**2
    if det <= 0 or ((c20 > 0.0 and c02 >= 0.0) or (c20 >= 0.0 and c02 > 0.0)):
        # polynomial does not have max. return maximum value in the data:
        return coord, 'WARNING:BADFIT', np.s_[y1:y2, x1:x2]

    xm = (c01 * c11 - 2.0 * c02 * c10) / det
    ym = (c10 * c11 - 2.0 * c01 * c20) / det

    if 0.0 <= xm <= (nx - 1.0) and 0.0 <= ym <= (ny - 1.0):
        coord = (xm, ym)
        fit_status = 'SUCCESS'

    else:
        xm = 0.0 if xm < 0.0 else min(xm, nx - 1.0)
        ym = 0.0 if ym < 0.0 else min(ym, ny - 1.0)
        fit_status = 'WARNING:EDGE'

    return coord, fit_status, np.s_[y1:y2, x1:x2]
