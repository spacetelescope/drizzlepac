import copy,os
import numpy as np

import pywcs
import stwcs
import pyfits

from stwcs import distortion
from stwcs.distortion import utils
from stwcs.wcsutil import wcscorr
from pytools import fileutil as fu
from stimage import xyxymatch

import catalogs
import linearfit
import updatehdr
import util

class Image(object):
    """ Primary class to keep track of all WCS and catalog information for
        a single input image. This class also performs all matching and fitting.
    """
    def __init__(self,filename,input_catalogs=None,**kwargs):
        """
        Parameters
        ----------
        filename: str
            filename for image

        input_catalogs: list of str or None
            filename of catalog files for each chip, if specified by user

        kwargs: dict
            parameters necessary for processing derived from input configObj object

        """
        self.name = filename
        self.rootname = filename[:filename.find('.')]
        self.origin = 1
        self.pars = kwargs
        if input_catalogs is not None and kwargs['xyunits'] == 'degrees':
            # Input was a catalog of sky positions, so no WCS or image needed
            use_wcs = False
            num_sci = 0
        else:
            # WCS required, so verify that we can get one
            # Need to count number of SCI extensions
            #  (assume a valid WCS with each SCI extension)
            num_sci,extname = count_sci_extensions(filename)
            if num_sci < 1:
                print 'ERROR: No Valid WCS available for %s',filename
                raise InputError
            use_wcs = True
        # Record this for use with methods
        self.use_wcs = use_wcs

        # Need to generate a separate catalog for each chip
        self.chip_catalogs = {}
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
            catalog = catalogs.generateCatalog(wcs,mode=catalog_mode,catalog=source,**kwargs)
            catalog.buildCatalogs() # read in and convert all catalog positions to RA/Dec
            self.chip_catalogs[sci_extn] = {'catalog':catalog,'wcs':wcs}

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

        # Set up products which need to be computed by methods of this class
        self.outxy = None
        self.refWCS = None # reference WCS assigned for the final fit
        self.matches = {'image':None,'ref':None} # stores matched list of coordinates for fitting
        self.fit = None # stores result of fit
        self.match_pars = None
        self.fit_pars = None
        self.identityfit = False # set to True when matching/fitting to itself
        self.goodmatch = True # keep track of whether enough matches were found for a fit
        
    def get_wcs(self):
        """ Helper method to return a list of all the input WCS objects associated
            with this image
        """
        wcslist = []
        for chip in self.chip_catalogs:
            wcslist.append(self.chip_catalogs[chip]['wcs'])
        return wcslist

    def buildSkyCatalog(self):
        """ Convert sky catalog for all chips into a single catalog for
            the entire field-of-view of this image
        """
        ralist = []
        declist = []
        fluxlist = []
        for scichip in self.chip_catalogs:
            skycat = self.chip_catalogs[scichip]['catalog'].radec
            xycat = self.chip_catalogs[scichip]['catalog'].xypos
            if skycat is not None:
                ralist.append(skycat[0])
                declist.append(skycat[1])
                if len(xycat) > 2:
                    fluxlist.append(xycat[2])
                else:
                    fluxlist.append([999.0]*len(skycat[0]))

        self.all_radec = [np.concatenate(ralist),np.concatenate(declist),np.concatenate(fluxlist)]
        self.all_radec_orig = copy.deepcopy(self.all_radec)
        

    def buildDefaultRefWCS(self):
        """ Generate a default reference WCS for this image
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
            print 'Reference WCS not a valid HSTWCS object'
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
            by the user
            It keeps a copy of the original full list in order to support iteration
        """
        _sortKeys = ['fluxmax','fluxmin','nbright']
        clip_catalog = False
        clip_prefix = ''
        for k in _sortKeys:
            for p in self.pars.keys():
                pindx = p.find(k)
                if pindx >= 0 and self.pars[p] is not None:
                    clip_catalog = True
                    print 'found a match for ',p,self.pars[p]
                    # find prefix (if any)
                    clip_prefix = p[:pindx].strip()

        if clip_catalog:
            # Start by clipping by any specified flux range
            if self.pars[clip_prefix+'fluxmax'] is not None or self.pars[clip_prefix+'fluxmin'] is not None:
                clip_catalog = True
                if self.pars[clip_prefix+'fluxmin'] is not None:
                    fluxmin = self.pars[clip_prefix+'fluxmin']
                else:
                    fluxmin = self.all_radec[2].min()

                if self.pars[clip_prefix+'fluxmax'] is not None:
                    fluxmin = self.pars[clip_prefix+'fluxmax']
                else:
                    fluxmin = self.all_radec[2].max()

            if self.pars[clip_prefix+'nbright'] is not None:
                clip_catalog = True
                # pick out only the brightest 'nbright' sources
                if self.pars[clip_prefix+'fluxunits'] == 'mag':
                    nbslice = slice(None,nbright)
                else:
                    nbslice = slice(nbright,None)

            all_radec = copy.deepcopy(self.all_radec_orig) # work on copy of all original data
            nbright_indx = np.argsort(all_radec[2])[nbslice] # find indices of brightest
            self.all_radec[0] = all_radec[0][nbright_indx]
            self.all_radec[1] = all_radec[1][nbright_indx]
            self.all_radec[2] = all_radec[2][nbright_indx]


    def match(self,ref_outxy, refWCS, refname, **kwargs):
        """ Uses xyxymatch to cross-match sources between this catalog and
            a reference catalog (refCatalog).
        """
        self.sortSkyCatalog() # apply any catalog sorting specified by the user
        self.transformToRef(refWCS)
        self.refWCS = refWCS
        # extract xyxymatch parameters from input parameters
        matchpars = kwargs.copy()
        self.match_pars = matchpars
        minobj = matchpars['minobj'] # needed for later
        del matchpars['minobj'] # not needed in xyxymatch

        # Check to see whether or not it is being matched to itself
        if (ref_outxy.shape == self.outxy.shape) and (ref_outxy == self.outxy).all():
            self.identityfit = True
            print 'NO fit performed for reference image: ',self.name,'\n'
        else:
            xoff = 0.
            yoff = 0.
            if matchpars['xoffset'] is not None:
                xoff = matchpars['xoffset']
            if matchpars['yoffset'] is not None:
                yoff = matchpars['yoffset']
            
            # convert tolerance from units of arcseconds to pixels, as needed
            radius = matchpars['searchrad']
            if matchpars['searchunits'] == 'arcseconds':
                radius /= refWCS.pscale

            xyoff = (xoff,yoff)
            matches = xyxymatch(self.outxy,ref_outxy,origin=xyoff,tolerance=radius,separation=matchpars['separation'])
            if len(matches) > minobj:
                self.matches['image'] = np.column_stack([matches['input_x'][:,np.newaxis],matches['input_y'][:,np.newaxis]])
                self.matches['ref'] = np.column_stack([matches['ref_x'][:,np.newaxis],matches['ref_y'][:,np.newaxis]])
                print 'Found %d matches for %s...'%(len(matches),self.name)

                if self.pars['writecat']:
                    matchfile = open(self.catalog_names['match'],mode='w+')
                    matchfile.write('#Reference: %s\n'%refname)
                    matchfile.write('#Input: %s\n'%self.name)
                    matchfile.write('#Ref_X        Ref_Y            Input_X        Input_Y         Ref_ID    Input_ID\n')
                    for i in xrange(len(matches['input_x'])):
                        linestr = "%0.6f    %0.6f        %0.6f    %0.6f        %d    %d\n"%\
                            (matches['ref_x'][i],matches['ref_y'][i],\
                             matches['input_x'][i],matches['input_y'][i],
                            matches['ref_idx'][i],matches['input_idx'][i])
                        matchfile.write(linestr)
                    matchfile.close()
            else:
                print 'WARNING: Not enough matches found for input image: ',self.name
                self.goodmatch = False


    def performFit(self,**kwargs):
        """ Perform a fit between the matched sources

            Parameters
            ----------
            kwargs: dict
                Parameter necessary to perform the fit; namely, *fitgeometry*

            Notes
            -----
            This task still needs to implement (eventually) interactive iteration of
                   the fit to remove outliers
        """
        pars = kwargs.copy()
        self.fit_pars = pars

        if not self.identityfit:
            if self.matches is not None and self.goodmatch:
                if pars['fitgeometry'] in ['rscale','general']:
                    self.fit = linearfit.fit_arrays(self.matches['image'],self.matches['ref'])
                else:
                    self.fit = linearfit.fit_shifts(self.matches['image'],self.matches['ref'])
                print 'Computed fit for ',self.name,': '
                print 'XSH: %0.6g  YSH: %0.6g    ROT: %0.6g    SCALE: %0.6g\n'%(
                    self.fit['offset'][0],self.fit['offset'][1], 
                    self.fit['rot'],self.fit['scale'][0])
        else:
            self.fit = {'offset':[0.0,0.0],'rot':0.0,'scale':[1.0]}

    def updateHeader(self,wcsname=None):
        """ Update header of image with shifts computed by *perform_fit()*
        """

        if not self.identityfit and self.goodmatch:
            updatehdr.updatewcs_with_shift(self.name,self.refWCS,wcsname=wcsname,
                xsh=self.fit['offset'][0],ysh=self.fit['offset'][1],rot=self.fit['rot'],scale=self.fit['scale'][0])
        if self.identityfit:
            # Update header using 'updatewcs'
            stwcs.updatewcs.updatewcs(self.name,wcsname=wcsname)
            # Create WCSCORR table to keep track of WCS revisions anyway
            wcscorr.init_wcscorr(self.name)

    def write_skycatalog(self,filename):
        """ Write out the all_radec catalog for this image to a file
        """
        ralist = self.all_radec[0].tolist()
        declist = self.all_radec[1].tolist()
        f = open(filename,'w')
        f.write("#Sky positions for: "+self.name+'\n')
        f.write("#RA        Dec\n")
        f.write("#(deg)     (deg)\n")
        for i in xrange(len(ralist)):
            f.write('%0.8g  %0.8g\n'%(ralist[i],declist[i]))
        f.close()

    def write_outxy(self,filename):
        """ Write out the output(transformed) XY catalog for this image to a file
        """
        f = open(filename,'w')
        f.write("#Pixel positions for: "+self.name+'\n')
        f.write("#X           Y\n")
        f.write("#(pix)       (pix)\n")
        for i in xrange(self.all_radec[0].shape[0]):
            f.write('%g  %g\n'%(self.outxy[i,0],self.outxy[i,1]))
        f.close()

    def get_shiftfile_row(self):
        """ Return the information for a shiftfile for this image to provide
            compatability with the IRAF-based MultiDrizzle
        """
        if self.fit is not None:
            rowstr = '%s    %g  %g    %g     %g\n'%(self.name,self.fit['offset'][0],self.fit['offset'][1],self.fit['rot'],self.fit['scale'][0])
        else:
            rowstr = None
        return rowstr


    def clean(self):
        """ Remove intermediate files created
        """            
        for f in self.catalog_names:
            if f == 'match':
                if os.path.exists(self.catalog_names[f]): 
                    print 'Deleting intermediate match file: ',self.catalog_names[f]
                    os.remove(self.catalog_names[f])
            else:
                for extn in f:
                    if os.path.exists(extn): 
                        print 'Deleting intermediate catalog: ',extn
                        os.remove(extn)

class RefImage(object):
    """ This class provides all the information needed by to define a reference
        tangent plane and list of source positions on the sky.
    """
    def __init__(self,wcs_list,catalog,**kwargs):
        if isinstance(wcs_list,list):
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
        self.origin = 1
        self.pars = kwargs
        
        # convert sky positions to X,Y positions on reference tangent plane
        self.transformToRef()

    def write_skycatalog(self,filename):
        """ Write out the all_radec catalog for this image to a file
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
        to reference tangent plane (self.wcs) to create output X,Y positions
        """
        if self.pars.has_key('refxyunits') and self.pars['refxyunits'] == 'pixels':
            print 'Creating RA/Dec positions for reference sources...'
            self.outxy = np.column_stack([self.all_radec[0][:,np.newaxis],self.all_radec[1][:,np.newaxis]])
            skypos = self.wcs.wcs_pix2sky(self.all_radec[0],self.all_radec[1],self.origin)
            self.all_radec = np.column_stack([skypos[0][:,np.newaxis],skypos[1][:,np.newaxis]])
        else:
            print 'Converting RA/Dec positions of reference sources to X,Y positions in reference WCS...'
            self.refWCS = self.wcs
            outxy = self.wcs.wcs_sky2pix(self.all_radec[0],self.all_radec[1],self.origin)
            # convert outxy list to a Nx2 array
            self.outxy = np.column_stack([outxy[0][:,np.newaxis],outxy[1][:,np.newaxis]])
            

    def get_shiftfile_row(self):
        """ Return the information for a shiftfile for this image to provide
            compatability with the IRAF-based MultiDrizzle
        """
        rowstr = '%s    0.0  0.0    0.0     1.0\n'%(self.name)
        return rowstr

    def clean(self):
        """ Remove intermediate files created
        """
        if not util.is_blank(self.catalog.catname) and os.path.exists(self.catalog.catname):
            os.remove(self.catalog.catname)

def build_referenceWCS(catalog_list):
    """ Compute default reference WCS from list of Catalog objects
    """
    wcslist = []
    for catalog in catalog_list:
        for scichip in catalog.catalogs:
            wcslist.append(catalog.catalogs[scichip]['wcs'])
    return utils.output_wcs(wcslist)

def count_sci_extensions(filename):
    """ Return the number of SCI extensions and the EXTNAME from a input MEF file
    """
    num_sci = 0
    extname = 'SCI'
    for extn in fu.openImage(filename):
        if extn.header.has_key('extname') and extn.header['extname'] == extname:
            num_sci += 1
    if num_sci == 0:
        extname = 'PRIMARY'
        num_sci = 1

    return num_sci,extname
