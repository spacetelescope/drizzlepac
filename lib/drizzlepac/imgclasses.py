import copy,os
import numpy as np

import pywcs
import stwcs
import pyfits

from stwcs import distortion
from stwcs.distortion import utils
from stwcs.wcsutil import wcscorr
from stwcs.wcsutil import headerlet
from stwcs.wcsutil import altwcs
from stsci.tools import fileutil as fu
from stsci.stimage import xyxymatch
from stsci.tools import logutil,textutil

import catalogs
import linearfit
import updatehdr
import util
import tweakutils

log = logutil.create_logger(__name__)

sortKeys = ['fluxmax','fluxmin','nbright','fluxunits']

class Image(object):
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
        self.perform_update = kwargs['updatehdr']
        if self.perform_update:
            self.open_mode = 'update'
        else:
            self.open_mode = 'readonly'
        self.name = filename
        self.openFile()

        # try to verify whether or not this image has been updated with
        # a full distortion model
        numsci,sciname = util.count_sci_extensions(filename)
        numwht = fu.countExtn(filename,extname='WHT')

        if sciname == 'PRIMARY': sciext = 0
        else: sciext = 1
        wnames = altwcs.wcsnames(self.hdulist,ext=sciext)
        # If no WCSNAME keywords were found, raise the possibility that
        # the images have not been updated fully and may result in inaccurate
        # alignment
        # use 'numwht' != 0 to indicate a DRZ file has been specified as input
        if len(wnames) == 0 and numwht == 0:
            print textutil.textbox('WARNING:\n'
            'Image %s may not have the full correct '%filename+
            'WCS solution in the header as created by stwcs.updatewcs '
            'Image alignment may not be taking into account the full '
            'the full distortion solution. \n'
            'Turning on the "updatewcs" parameter would insure '
            'that each image uses the full distortion model when '
            'aligning this image.\n', width=60
            )

        self.rootname = filename[:filename.find('.')]
        self.origin = 1
        self.pars = kwargs
        self.exclusions = exclusions
        self.verbose = kwargs['verbose']
        self.interactive = kwargs.get('interactive',True)

        print 'Defining source catalogs for: ',filename

        if input_catalogs is not None and kwargs['xyunits'] == 'degrees':
            # Input was a catalog of sky positions, so no WCS or image needed
            use_wcs = False
            num_sci = 0
        else:
            # WCS required, so verify that we can get one
            # Need to count number of SCI extensions
            #  (assume a valid WCS with each SCI extension)
            num_sci,extname = util.count_sci_extensions(filename)
            if num_sci < 1:
                print >> sys.stderr,textutil.textbox(
                        'ERROR: No Valid WCS available for %s'%filename)
                raise InputError
            use_wcs = True
        # Record this for use with methods
        self.use_wcs = use_wcs
        self.num_sci = num_sci
        self.ext_name = extname

        # Need to generate a separate catalog for each chip
        self.chip_catalogs = {}
        self.xy_catalog = [[],[],[],[]]
        num_sources = 0
        # For each SCI extension, generate a catalog and WCS
        for sci_extn in range(1,num_sci+1):
            extnum = fu.findExtname(pyfits.open(filename),extname,extver=sci_extn)
            if extnum is None: extnum = 0
            chip_filename = filename+'[%d]'%(extnum)
            if use_wcs:
                wcs = stwcs.wcsutil.HSTWCS(chip_filename)
            if input_catalogs is None:
                # if we already have a set of catalogs provided on input,
                #  we only need the array to get original XY input positions
                source = chip_filename
                catalog_mode='automatic'
            else:
                source = input_catalogs[sci_extn-1]
                catalog_mode='user'

            if exclusions is not None and len(exclusions) >= sci_extn:
                if exclusions[sci_extn-1] not in ['None',' ','INDEF']:
                    excludefile=exclusions[sci_extn-1]
                else:
                    excludefile = None
            else:
                excludefile = None
            kwargs['start_id'] = num_sources
            catalog = catalogs.generateCatalog(wcs,mode=catalog_mode,catalog=source,**kwargs)
            # read in and convert all catalog positions to RA/Dec
            catalog.buildCatalogs(exclusions=excludefile)
            num_sources += catalog.num_objects
            self.chip_catalogs[sci_extn] = {'catalog':catalog,'wcs':wcs}
            # Merge input X,Y positions from all chips into a single catalog
            # This will be used for writing output matched source catalogs
            for i,col in enumerate(self.xy_catalog):
                if catalog.xypos is not None and len(catalog.xypos) > 0 and i < len(catalog.xypos):
                    col = col.extend(catalog.xypos[i])

        self.catalog_names = {}
        # Build full list of all sky positions from all chips
        self.buildSkyCatalog()
        if self.pars['writecat']:
            catname = self.rootname+"_sky_catalog.coo"
            self.catalog_names['match'] = self.rootname+"_xy_catalog.match"
            self.write_skycatalog(catname)
            self.catalog_names['sky'] = catname # Keep track of catalogs being written out
            for nsci in range(1,num_sci+1):
                catname = "%s_sci%d_xy_catalog.coo"%(self.rootname,nsci)
                self.chip_catalogs[nsci]['catalog'].writeXYCatalog(catname)
                # Keep track of catalogs being written out
                if 'input_xy' not in self.catalog_names:
                    self.catalog_names['input_xy'] = []
                self.catalog_names['input_xy'].append(catname)
            self.catalog_names['fitmatch'] = self.rootname+"_catalog_fit.match"

        print '\n'
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
        if self.hdulist:
            self.hdulist.close()
            self.hdulist = None

    def openFile(self):
        """ Open file and set up filehandle for image file
        """
        if not hasattr(self, 'hdulist') or self.hdulist is None:
            self.hdulist = fu.openImage(self.name,mode=self.open_mode)

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
                if len(xycat) > 2:
                    fluxlist.append(xycat[2])
                    idlist.append(xycat[3])
                else:
                    fluxlist.append([999.0]*len(skycat[0]))
                    idlist.append(np.arange(len(skycat[0])))


                self.all_radec = [np.concatenate(ralist),np.concatenate(declist),
                        np.concatenate(fluxlist),np.concatenate(idlist)]
                self.all_radec_orig = copy.deepcopy(self.all_radec)

    def buildDefaultRefWCS(self):
        """ Generate a default reference WCS for this image.
        """
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
            print >> sys.stderr, textutil.textbox(
                            'Reference WCS not a valid HSTWCS object')
            raise ValueError
        # Need to concatenate catalogs from each input
        if self.outxy is None or force:
            outxy = ref_wcs.wcs_sky2pix(self.all_radec[0],self.all_radec[1],self.origin)
            # convert outxy list to a Nx2 array
            self.outxy = np.column_stack([outxy[0][:,np.newaxis],outxy[1][:,np.newaxis]])
            if self.pars['writecat']:
                catname = self.rootname+"_refxy_catalog.coo"
                self.write_outxy(catname)
                self.catalog_names['ref_xy'] = catname

    def sortSkyCatalog(self):
        """ Sort and clip the source catalog based on the flux range specified
            by the user.
            It keeps a copy of the original full list in order to support iteration.
        """
        if len(self.all_radec_orig[2].nonzero()[0]) == 0:
            warn_str = "Source catalog NOT trimmed by flux/mag. No fluxes read in for sources!"
            print '\nWARNING: ',warn_str,'\n'
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
            if self.pars[clip_prefix+'fluxmax'] is not None or \
                    self.pars[clip_prefix+'fluxmin'] is not None:
                clip_catalog = True
                if self.pars[clip_prefix+'fluxmin'] is not None:
                    fluxmin = self.pars[clip_prefix+'fluxmin']
                else:
                    fluxmin = self.all_radec[2].min()

                if self.pars[clip_prefix+'fluxmax'] is not None:
                    fluxmax = self.pars[clip_prefix+'fluxmax']
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

    def match(self,refimage, **kwargs):
        """ Uses xyxymatch to cross-match sources between this catalog and
            a reference catalog (refCatalog).
        """
        ref_outxy = refimage.outxy
        refWCS = refimage.wcs
        refname = refimage.name
        ref_inxy = refimage.xy_catalog

        print 'Matching sources from ',self.name,' with sources from reference image',refname
        self.sortSkyCatalog() # apply any catalog sorting specified by the user
        self.transformToRef(refWCS)
        self.refWCS = refWCS
        # extract xyxymatch parameters from input parameters
        matchpars = kwargs.copy()
        self.match_pars = matchpars
        minobj = matchpars['minobj'] # needed for later
        del matchpars['minobj'] # not needed in xyxymatch

        # Check to see whether or not it is being matched to itself
        if (refname.strip() == self.name.strip()):# or (
#                ref_outxy.shape == self.outxy.shape) and (
#                ref_outxy == self.outxy).all():
            self.identityfit = True
            log.info('NO fit performed for reference image: %s\n'%self.name)
        else:
            # convert tolerance from units of arcseconds to pixels, as needed
            radius = matchpars['searchrad']
            if matchpars['searchunits'] == 'arcseconds':
                radius /= refWCS.pscale

            # Determine xyoff (X,Y offset) and tolerance to be used with xyxymatch
            use2d = matchpars['use2dhist']
            xyxysep = matchpars['separation']
            if use2d:
                hist_name = None
                if not self.interactive:
                    hist_name = 'hist2d_{0}.png'.format(self.rootname)
                zpxoff,zpyoff,flux,zpqual = tweakutils.build_xy_zeropoint(self.outxy,
                                    ref_outxy,searchrad=radius,
                                    histplot=matchpars['see2dplot'],
                                    figure_id = self.figure_id, plotname=hist_name)
                if matchpars['see2dplot'] and ('residplot' in matchpars and
                                               'No' in matchpars['residplot']):
                    if self.interactive:
                        a = raw_input("Press ENTER for next image, \n     'n' to continue without updating header or \n     'q' to quit immediately...\n")
                    else:
                        a = ' '
                    if 'n' in a.lower():
                        self.perform_update = False
                    if 'q' in a.lower():
                        self.quit_immediately = True

                if matchpars['see2dplot']:
                    self.figure_id += 1
                if zpqual is not None:
                    xyoff = (zpxoff,zpyoff)
                    # set tolerance as well
                    # This value allows initial guess to be off by 1 in both and
                    # still pick up the identified matches
                    xyxytolerance = 1.5
                    #xyxysep = 0.
                else:
                    use2d = False
            if not use2d:
                xoff = 0.
                yoff = 0.
                if not util.is_blank(matchpars['xoffset']):
                    xoff = matchpars['xoffset']
                if not util.is_blank(matchpars['yoffset']):
                    yoff = matchpars['yoffset']
                xyoff = (xoff,yoff)
                # set tolerance based on user-specified input
                xyxytolerance = matchpars['tolerance']

            matches = xyxymatch(self.outxy,ref_outxy,origin=xyoff,
                                tolerance=xyxytolerance,separation=xyxysep)

            if len(matches) > minobj:
                self.matches['image'] = np.column_stack([matches['input_x'][:,
                                np.newaxis],matches['input_y'][:,np.newaxis]])
                self.matches['ref'] = np.column_stack([matches['ref_x'][:,
                                np.newaxis],matches['ref_y'][:,np.newaxis]])
                self.matches['ref_idx'] = matches['ref_idx']
                self.matches['img_idx'] = self.all_radec[3][matches['input_idx']]
                self.matches['ref_orig_xy'] = np.column_stack([
                                    np.array(ref_inxy[0])[matches['ref_idx']][:,np.newaxis],
                                    np.array(ref_inxy[1])[matches['ref_idx']][:,np.newaxis]])
                self.matches['img_orig_xy'] = np.column_stack([
                    np.array(self.xy_catalog[0])[matches['input_idx']][:,np.newaxis],
                    np.array(self.xy_catalog[1])[matches['input_idx']][:,np.newaxis]])
                print 'Found %d matches for %s...'%(len(matches),self.name)

                if self.pars['writecat']:
                    matchfile = open(self.catalog_names['match'],mode='w+')
                    matchfile.write('#Reference: %s\n'%refname)
                    matchfile.write('#Input: %s\n'%self.name)
                    title = '#Ref_X        Ref_Y        '
                    title += 'Input_X    Input_Y        '
                    title += 'Ref_X0    Ref_Y0        '
                    title += 'Input_X0    Input_Y0        '
                    title += 'Ref_ID    Input_ID\n'
                    fmtstr = '%0.6f    %0.6f        '*4
                    fmtstr += '%d    %d\n'
                    matchfile.write(title)
                    for i in xrange(len(matches['input_x'])):
                        linestr = fmtstr%\
                            (matches['ref_x'][i],matches['ref_y'][i],\
                             matches['input_x'][i],matches['input_y'][i],
                            self.matches['ref_orig_xy'][:,0][i],
                            self.matches['ref_orig_xy'][:,1][i],
                            self.matches['img_orig_xy'][:,0][i],
                            self.matches['img_orig_xy'][:,1][i],
                            matches['ref_idx'][i],matches['input_idx'][i])
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
        pars = kwargs.copy()
        self.fit_pars = pars

        self.fit = {'offset':[0.0,0.0],'rot':0.0,'scale':[1.0],'rms':[0.0,0.0],
                    'rms_keys':{'RMS_RA':0.0,'RMS_DEC':0.0,'NMATCH':0},
                    'fit_matrix':[[1.0,0.0],[0.0,1.0]]}

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
                if pars['fitgeometry'] != 'general':
                    self.fit['fit_matrix'] = None

                print 'Computed ',pars['fitgeometry'],' fit for ',self.name,': '
                print 'XSH: %0.6g  YSH: %0.6g    ROT: %0.6g    SCALE: %0.6g'%(
                    self.fit['offset'][0],self.fit['offset'][1],
                    self.fit['rot'],self.fit['scale'][0])
                print 'XRMS: %0.6g    YRMS: %0.6g\n'%(
                        self.fit['rms'][0],self.fit['rms'][1])
                print 'RMS_RA: %g (deg)   RMS_DEC: %g (deg)\n'%(
                        self.fit['rms_keys']['RMS_RA'],
                        self.fit['rms_keys']['RMS_DEC'])
                print 'Final solution based on ',self.fit['rms_keys']['NMATCH'],' objects.'

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
                        plotname=resid_name, title=title_str)
                    if self.interactive:
                        a = raw_input("Press ENTER for next image, \n     'n' to continue without updating header or \n     'q' to quit immediately...\n")
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
            crval_rms = self.refWCS.wcs_pix2sky([crpix],1)[0]
            rms_ra,rms_dec = np.abs(crval_rms - self.refWCS.wcs.crval)
            nmatch = self.fit['resids'].shape[0]
        else:
            rms_ra = 0.0
            rms_dec = 0.0
            nmatch = 0
        return {'RMS_RA':rms_ra,'RMS_DEC':rms_dec,'NMATCH':nmatch}

    def updateHeader(self,wcsname=None):
        """ Update header of image with shifts computed by *perform_fit()*.
        """
        # Insure filehandle is open and available...
        self.openFile()

        verbose_level = 1
        if not self.perform_update:
            verbose_level = 0
        # Create WCSCORR table to keep track of WCS revisions anyway
        if self.perform_update:
            wcscorr.init_wcscorr(self.hdulist)

        extlist = []
        wcscorr_extname = self.ext_name
        if self.num_sci == 1 and self.ext_name == "PRIMARY":
            extlist = [0]
        else:
            for ext in range(1,self.num_sci+1):
                extlist.append((self.ext_name,ext))
                # add WCSNAME to SCI headers, if not provided (such as for
                # drizzled images directly obtained from the archive pre-AD)
                if ('wcsname' not in self.hdulist[self.ext_name,ext].header and
                    self.hdulist.fileinfo(0)['filemode'] == 'update'):
                    self.hdulist[self.ext_name,ext].header['wcsname'] = 'Default'

        next_key = altwcs.next_wcskey(pyfits.getheader(self.name,extlist[0]))

        if not self.identityfit and self.goodmatch and \
                self.fit['offset'][0] != np.nan:
            updatehdr.updatewcs_with_shift(self.hdulist,self.refWCS,wcsname=wcsname,
                xsh=self.fit['offset'][0],ysh=self.fit['offset'][1],
                rot=self.fit['rot'],scale=self.fit['scale'][0],
                fit=self.fit['fit_matrix'], verbose=verbose_level,
                xrms=self.fit['rms_keys']['RMS_RA'],yrms=self.fit['rms_keys']['RMS_DEC'])

            wnames = altwcs.wcsnames(self.hdulist,ext=extlist[0])
            altkeys = []
            for k in wnames:
                if wnames[k] == wcsname:
                    altkeys.append(k)
            if len(altkeys) > 1 and ' ' in altkeys:
                altkeys.remove(' ')
            next_key = altkeys[-1]
            if self.perform_update:
                log.info('    Writing out new WCS to alternate WCS: "%s"'%next_key)

            self.next_key = next_key
        else: #if self.identityfit or not self.goodmatch:
            if self.perform_update:
                log.info('    Saving Primary WCS to alternate WCS: "%s"'%next_key)
                # archive current WCS as alternate WCS with specified WCSNAME
                # Start by archiving original PRIMARY WCS
                wnames = altwcs.wcsnames(self.hdulist,ext=extlist[0])
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
                        self.hdulist[extlist[0]].header['wscname'] = ''
                        wnames[' '] = ''
                    pri_wcsname = wnames[' ']
                altwcs.archiveWCS(self.hdulist,extlist,
                                    wcskey=next_key,wcsname=pri_wcsname,
                                    reusekey=True)
                # Find key for next WCS and save again to replicate an updated solution
                next_key = altwcs.next_wcskey(self.hdulist[extlist[0]].header)
                # save again using new WCSNAME
                altwcs.archiveWCS(self.hdulist,extlist,
                    wcskey=next_key,wcsname=wcsname)
                # update WCSNAME to be the new name
                for ext in extlist:
                    self.hdulist[ext].header['WCSNAME'] = wcsname
            self.next_key = ' '

        # add FIT values to image's PRIMARY header
        fimg = self.hdulist

        if wcsname in ['',' ',None,"INDEF"]:
            wcsname = 'TWEAK'
        # Record values for the fit with both the PRIMARY WCS being updated
        # and the alternate WCS which will be created.
        for ext in extlist:
            fimg[ext].header.update('FITNAME'+next_key,wcsname)
            for kw in self.fit['rms_keys']:
                fimg[ext].header.update(kw+next_key,
                        self.fit['rms_keys'][kw],after='FITNAME'+next_key)

        if self.perform_update:
            log.info('Updating WCSCORR table with new WCS solution "%s"'%wcsname)
            wcscorr.update_wcscorr(fimg,wcs_id=wcsname,extname=self.ext_name)
        #fimg.close()
        self.hdulist = fimg

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
        headerlet.write_headerlet(self.hdulist, pars['hdrname'],
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
        for i in xrange(len(ralist)):
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
            log.info('Creating catalog for the fit: %s'%self.catalog_names['fitmatch'])
            f = open(self.catalog_names['fitmatch'],'w')
            f.write('# Input image: %s\n'%self.rootname)
            f.write('# Coordinate mapping parameters: \n')
            f.write('#    X and Y rms: %15.6g  %15.6g\n'%(self.fit['rms'][0],
                                                        self.fit['rms'][1]))
            f.write('#    X and Y shift: %15.6g  %15.6g\n'%(self.fit['offset'][0],
                                                            self.fit['offset'][1]))
            f.write('#    X and Y scale: %15.6g  %15.6g\n'%(self.fit['scale'][1],
                                                            self.fit['scale'][2]))
            f.write('#    X and Y rotation: %15.6g \n'%(self.fit['rot']))
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
            #
            # Need to add chip ID for each matched source to the fitmatch file
            # The chip information can be extracted from the following source:
            #
            #     self.chip_catalogs[sci_extn] = {'catalog':catalog,'wcs':wcs}
            #     xypos = catalog.xypos
            img_chip_id = self.fit['img_indx'].copy()
            for sci_extn in range(1,self.num_sci+1):
                if self.chip_catalogs[sci_extn]['catalog'].xypos is not None:
                    img_indx_orig = self.chip_catalogs[sci_extn]['catalog'].xypos[-1]
                    chip_min = img_indx_orig.min()
                    chip_max = img_indx_orig.max()
                    cid = np.bitwise_and((img_chip_id >= chip_min),(img_chip_id <= chip_max))
                    img_chip_id = np.where(cid, sci_extn, img_chip_id)
            #
            f.write('#\n')
            f.close()
            fitvals = self.fit['img_coords']+self.fit['resids']
            xydata = [[self.fit['ref_coords'][:,0],self.fit['ref_coords'][:,1],
                      self.fit['img_coords'][:,0],self.fit['img_coords'][:,1],
                      fitvals[:,0],fitvals[:,1],
                      self.fit['resids'][:,0],self.fit['resids'][:,1],
                      self.fit['ref_orig_xy'][:,0],
                      self.fit['ref_orig_xy'][:,1],
                      self.fit['img_orig_xy'][:,0],
                      self.fit['img_orig_xy'][:,1]],
                      [self.fit['ref_indx'],self.fit['img_indx'],img_chip_id]
                    ]
            tweakutils.write_xy_file(self.catalog_names['fitmatch'],xydata,
                                        append=True,format=["%15.6f","%8d"])

    def write_outxy(self,filename):
        """ Write out the output(transformed) XY catalog for this image to a file.
        """
        f = open(filename,'w')
        f.write("#Pixel positions for: "+self.name+'\n')
        f.write("#X           Y\n")
        f.write("#(pix)       (pix)\n")
        for i in xrange(self.all_radec[0].shape[0]):
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

class RefImage(object):
    """ This class provides all the information needed by to define a reference
        tangent plane and list of source positions on the sky.
    """
    def __init__(self,wcs_list,catalog,xycatalog=None,**kwargs):
        if isinstance(wcs_list,str):
            # Input was a filename for the reference image
            froot,fextn = fu.parseFilename(wcs_list)
            if fextn is None:
                num_sci,extname = util.count_sci_extensions(froot)
                if num_sci < 1:
                    fextn='[0]'
                else:
                    fextn = '[%s,1]'%extname
                fname = froot+fextn
            else:
                fname = wcs_list
            self.wcs = stwcs.wcsutil.HSTWCS(fname)
        elif isinstance(wcs_list,list):
            # generate a reference tangent plane from a list of STWCS objects
            undistort = True
            if wcs_list[0].sip is None:
                undistort=False
            self.wcs = utils.output_wcs(wcs_list,undistort=undistort)
            self.wcs.filename = wcs_list[0].filename
        else:
            # only a single WCS provided, so use that as the definition
            if not isinstance(wcs_list,stwcs.wcsutil.HSTWCS): # User only provided a filename
                self.wcs = stwcs.wcsutil.HSTWCS(wcs_list)
            else: # User provided full HSTWCS object
                self.wcs = wcs_list

        self.name = self.wcs.filename
        self.refWCS = None
        # Interpret the provided catalog
        self.catalog = catalogs.RefCatalog(None,catalog,**kwargs)
        self.catalog.buildCatalogs()
        self.all_radec = self.catalog.radec

        # add input source positions to RefCatalog for reporting in final
        # matched output catalog files
        if xycatalog is not None:
            self.xy_catalog = xycatalog
        else:
            # Convert RA/Dec positions of source from refimage into
            # X,Y positions based on WCS of refimage
            self.xy_catalog = [[],[],[],[]]
            xypos = self.wcs.all_sky2pix(self.all_radec[0],self.all_radec[1],1)
            self.xy_catalog[0] = xypos[0]
            self.xy_catalog[1] = xypos[1]
            self.xy_catalog[2] = np.zeros(xypos[0].shape[0],dtype=np.float32)
            self.xy_catalog[3] = np.arange(xypos[0].shape[0],dtype=np.int32)

        self.outxy = None
        self.origin = 1
        self.pars = kwargs
        if self.all_radec is not None:
            # convert sky positions to X,Y positions on reference tangent plane
            self.transformToRef()


    def write_skycatalog(self,filename):
        """ Write out the all_radec catalog for this image to a file.
        """
        f = open(filename,'w')
        f.write("#Sky positions for: "+self.name+'\n')
        f.write("#RA        Dec\n")
        f.write("#(deg)     (deg)\n")
        for i in xrange(self.all_radec[0].shape[0]):
            f.write('%g  %g\n'%(self.all_radec[0][i],self.all_radec[1][i]))
        f.close()

    def transformToRef(self):
        """ Transform reference catalog sky positions (self.all_radec)
        to reference tangent plane (self.wcs) to create output X,Y positions.
        """
        if 'refxyunits' in self.pars and self.pars['refxyunits'] == 'pixels':
            log.info('Creating RA/Dec positions for reference sources...')
            self.outxy = np.column_stack([self.all_radec[0][:,np.newaxis],self.all_radec[1][:,np.newaxis]])
            skypos = self.wcs.wcs_pix2sky(self.all_radec[0],self.all_radec[1],self.origin)
            self.all_radec = np.column_stack([skypos[0][:,np.newaxis],skypos[1][:,np.newaxis]])
        else:
            log.info('Converting RA/Dec positions of reference sources from "%s" to '%self.name+
                        'X,Y positions in reference WCS...')
            self.refWCS = self.wcs
            outxy = self.wcs.wcs_sky2pix(self.all_radec[0],self.all_radec[1],self.origin)
            # convert outxy list to a Nx2 array
            self.outxy = np.column_stack([outxy[0][:,np.newaxis],outxy[1][:,np.newaxis]])


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
        print textutil.textbox(warnstr)
        return
