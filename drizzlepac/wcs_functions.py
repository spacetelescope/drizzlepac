"""

:Authors: Warren Hack

:License: :doc:`LICENSE`

"""
from astropy.io import fits as pyfits
import copy
import numpy as np
from numpy import linalg

from stsci.tools import fileutil, logutil
from . import util

from astropy import wcs
from stwcs import wcsutil
from stwcs.distortion import coeff_converter, utils
from stwcs.wcsutil import altwcs

DEFAULT_WCS_PARS = {'ra': None, 'dec': None, 'scale': None, 'rot': None,
                    'outnx': None, 'outny': None,
                    'crpix1': None, 'crpix2': None}

log = logutil.create_logger(__name__, level=logutil.logging.NOTSET)


# Default mapping function based on astropy.wcs
class WCSMap:
    """ Sample class to demonstrate how to define a coordinate transformation
    """
    def __init__(self, input, output, origin=1):
        # Verify that we have valid WCS input objects
        self.checkWCS(input, 'Input')
        self.checkWCS(output, 'Output')

        self.input = input
        self.output = copy.deepcopy(output)

        self.origin = origin
        self.shift = None
        self.rot = None
        self.scale = None

    def checkWCS(self, obj, name):
        try:
            assert isinstance(obj, wcs.WCS)
        except AssertionError:
            print(name + ' object needs to be an instance or subclass of a astropy.wcs.WCS object.')
            raise

    def forward(self, pixx, pixy):
        """ Transform the input pixx,pixy positions in the input frame
            to pixel positions in the output frame.

            This method gets passed to the drizzle algorithm.
        """
        # This matches WTRAXY results to better than 1e-4 pixels.
        skyx, skyy = self.input.all_pix2world(pixx, pixy, self.origin)
        result = self.output.wcs_world2pix(skyx, skyy, self.origin)
        return result

    def backward(self, pixx, pixy):
        """ Transform pixx,pixy positions from the output frame back onto their
            original positions in the input frame.
        """
        skyx, skyy = self.output.wcs_pix2world(pixx, pixy, self.origin)
        result = self.input.all_world2pix(skyx, skyy, self.origin)
        return result

    def get_pix_ratio(self):
        """ Return the ratio of plate scales between the input and output WCS.
            This is used to properly distribute the flux in each pixel in 'tdriz'.
        """
        return self.output.pscale / self.input.pscale

    def xy2rd(self, wcs, pixx, pixy):
        """ Transform input pixel positions into sky positions in the WCS provided.
        """
        return wcs.all_pix2world(pixx, pixy, 1)
    def rd2xy(self, wcs, ra, dec):
        """ Transform input sky positions into pixel positions in the WCS provided.
        """
        return wcs.wcs_world2pix(ra, dec, 1)

def get_pix_ratio_from_WCS(input, output):
    """ [Functional form of .get_pix_ratio() method of WCSMap]"""
    return output.pscale / input.pscale
##
#
# ### Default no-op transformation
#
##
class IdentityMap:
    def __init__(self, input, output):
        print('Applying identity transformation...')
        self.input = input
        self.output = output

    def forward(self, pixx, pixy):
        return pixx, pixy
##
#
# ### Linear transformation mapper
#
##
class LinearMap:
    def __init__(self, xsh=0.0, ysh=0.0, rot=0.0, scale=1.0):
        # Define rotation matrix
        _theta = np.deg2rad(rot)
        _mrot = np.zeros(shape=(2, 2), dtype=np.float64)
        _mrot[0] = (np.cos(_theta), np.sin(_theta))
        _mrot[1] = (-np.sin(_theta), np.cos(_theta))
        # apply scaling factor to rotation matrix
        self.transform = _mrot * scale
        # define offset to be applied to positions after rotation
        self.offset = [[xsh], [ysh]]

    def forward(self, pixx, pixy):
        return np.dot(self.transform, [pixx, pixy]) + self.offset


# Stand-alone functions for WCS handling
def get_hstwcs(filename, hdulist, extnum):
    """ Return the HSTWCS object for a given chip. """
    hdrwcs = wcsutil.HSTWCS(hdulist, ext=extnum)
    hdrwcs.filename = filename
    hdrwcs.expname = hdulist[extnum].header['expname']
    hdrwcs.extver = hdulist[extnum].header['extver']

    return hdrwcs


def build_hstwcs(crval1, crval2, crpix1, crpix2, naxis1, naxis2, pscale, orientat):
    """ Create an HSTWCS object for a default instrument without distortion
        based on user provided parameter values.
    """
    wcsout = wcsutil.HSTWCS()
    wcsout.wcs.crval = np.array([crval1, crval2])
    wcsout.wcs.crpix = np.array([crpix1, crpix2])
    wcsout.naxis1 = naxis1
    wcsout.naxis2 = naxis2
    wcsout.wcs.cd = fileutil.buildRotMatrix(orientat) * [-1, 1] * pscale / 3600.0
    # Synchronize updates with astropy.wcs/WCSLIB objects
    wcsout.wcs.set()
    wcsout.setPscale()
    wcsout.setOrient()
    wcsout.wcs.ctype = ['RA---TAN', 'DEC--TAN']

    return wcsout


def update_linCD(cdmat, delta_rot=0.0, delta_scale=1.0, cx=[0.0, 1.0], cy=[1.0, 0.0]):
    """ Modify an existing linear CD matrix with rotation and/or scale changes
        and return a new CD matrix.  If 'cx' and 'cy' are specified, it will
        return a distorted CD matrix.

        Only those terms which are varying need to be specified on input.
    """
    rotmat = fileutil.buildRotMatrix(delta_rot) * delta_scale
    new_lincd = np.dot(cdmat, rotmat)

    cxymat = np.array([[cx[1], cx[0]], [cy[1], cy[0]]])
    new_cd = np.dot(new_lincd, cxymat)

    return new_cd


def create_CD(orient, scale, cx=None, cy=None):
    """ Create a (un?)distorted CD matrix from the basic inputs.

    The 'cx' and 'cy' parameters, if given, provide the X and Y coefficients of
    the distortion as returned by reading the IDCTAB.  Only the first 2 elements
    are used and should correspond to the 'OC[X/Y]10' and 'OC[X/Y]11' terms in that
    order as read from the expanded SIP headers.

    The units of 'scale' should be 'arcseconds/pixel' of the reference pixel.
    The value of 'orient' should be the absolute orientation on the sky of the
    reference pixel.

    """
    cxymat = np.array([[cx[1], cx[0]], [cy[1], cy[0]]])
    rotmat = fileutil.buildRotMatrix(orient) * scale/3600.
    new_cd = np.dot(rotmat, cxymat)
    return new_cd

def ddtohms(xsky, ysky, verbose=False, precision=6):
    """ Convert sky position(s) from decimal degrees to HMS format. """
    xskyh = xsky / 15.
    xskym = (xskyh - np.floor(xskyh)) * 60.
    xskys = (xskym - np.floor(xskym)) * 60.
    yskym = (np.abs(ysky) - np.floor(np.abs(ysky))) * 60.
    yskys = (yskym - np.floor(yskym)) * 60.

    fmt = "%." + repr(precision) + "f"
    if isinstance(xskyh, np.ndarray):
        rah, dech = [], []
        for i in range(len(xskyh)):
            rastr = repr(int(xskyh[i])) + ':' + repr(int(xskym[i])) + ':' + fmt % (xskys[i])
            decstr = repr(int(ysky[i])) + ':' + repr(int(yskym[i])) + ':' + fmt % (yskys[i])
            rah.append(rastr)
            dech.append(decstr)
            if verbose:
                print('RA = ', rastr, ', Dec = ', decstr)
    else:
        rastr = repr(int(xskyh)) + ':' + repr(int(xskym)) + ':' + fmt % (xskys)
        decstr = repr(int(ysky)) + ':' + repr(int(yskym)) + ':' + fmt % (yskys)
        rah = rastr
        dech = decstr
        if verbose:
            print('RA = ', rastr, ', Dec = ', decstr)

    return rah, dech

# Functions for applying pixel-based transformations
#
def build_pixel_transform(chip, output_wcs):
    driz_pars = {}
    driz_pars['pxg'] = np.zeros([2, 2], dtype=np.float32)
    driz_pars['pyg'] = np.zeros([2, 2], dtype=np.float32)
    # Use default C mapping function
    _inwcs = np.zeros([8], dtype=np.float64)
    driz_pars['inwcs'] = convertWCS(output_wcs.wcs, _inwcs)

    indx = chip.outputNames['data'].find('.fits')
    driz_pars['coeffs_name'] = chip.outputNames['data'][:indx] + '_coeffs' + str(chip.detnum) + '.dat'
    #
    # Need to compute and write out coeffs files for each chip as well.
    #
    xcoeffs, ycoeffs = coeff_converter.sip2idc(chip.wcs)
    # account for the case where no IDCSCALE has been set, due to a
    # lack of IDCTAB or to 'coeffs=False'.
    idcscale = chip.wcs.idcscale
    if idcscale is None: idcscale = chip.wcs.pscale
    xcoeffs /= idcscale
    ycoeffs /= idcscale
    driz_pars['coeffs'] = [xcoeffs, ycoeffs]

    abxt, cdyt = wcsfit(chip.wcs, output_wcs)
    driz_pars['abxt'] = abxt
    driz_pars['cdyt'] = cdyt

    driz_pars['delta_rot'] = np.rad2deg(np.arctan2(abxt[1], cdyt[0]))
    # Compute scale from fit to allow WFPC2 (and similar) data to be handled correctly
    driz_pars['scale'] = 1. / np.sqrt(abxt[0]**2 + abxt[1]**2)
    driz_pars['tddalpha'] = chip.header['tddalpha']
    driz_pars['tddbeta'] = chip.header['tddbeta']

    return driz_pars

#
# Possibly need to generate a stand-alone interface for this function.
#
# ### Primary interface for creating the output WCS from a list of HSTWCS objects
def make_outputwcs(imageObjectList, output, configObj=None, perfect=False):
    """ Computes the full output WCS based on the set of input imageObjects
        provided as input, along with the pre-determined output name from
        process_input.  The user specified output parameters are then used to
        modify the default WCS to produce the final desired output frame.
        The input imageObjectList has the outputValues dictionary
        updated with the information from the computed output WCS.
        It then returns this WCS as a WCSObject(imageObject)
        instance.

    """
    if not isinstance(imageObjectList, list):
        imageObjectList = [imageObjectList]

    # Compute default output WCS, replace later if user specifies a refimage
    hstwcs_list = []
    undistort = True
    for img in imageObjectList:
        chip_wcs = copy.deepcopy(img.getKeywordList('wcs'))
        # IF the user turned off use of coeffs (coeffs==False)
        if not configObj['coeffs']:
            for cw in chip_wcs:
                # Turn off distortion model for each input
                cw.sip = None
                cw.cpdis1 = None
                cw.cpdis2 = None
                cw.det2im = None
            undistort = False
        hstwcs_list += chip_wcs
    if not undistort and len(hstwcs_list) == 1:
        default_wcs = hstwcs_list[0].deepcopy()
    else:
        default_wcs = utils.output_wcs(hstwcs_list, undistort=undistort)

    if perfect:
        default_wcs.wcs.cd = make_perfect_cd(default_wcs)

    # Turn WCS instances into WCSObject instances
    outwcs = createWCSObject(output, default_wcs, imageObjectList)

    # Merge in user-specified attributes for the output WCS
    # as recorded in the input configObj object.
    final_pars = DEFAULT_WCS_PARS.copy()

    # More interpretation of the configObj needs to be done here to translate
    # the input parameter names to those understood by 'mergeWCS' as defined
    # by the DEFAULT_WCS_PARS dictionary.
    single_step = configObj[util.getSectionName(configObj, 3)]
    singleParDict = configObj[util.getSectionName(configObj, '3a')].copy()
    if single_step['driz_separate'] and singleParDict['driz_sep_wcs']:
        single_pars = DEFAULT_WCS_PARS.copy()
        del singleParDict['driz_sep_wcs']
        keyname = 'driz_sep_'
        for key in singleParDict:
            k = key[len(keyname):]
            if k != 'refimage':
                single_pars[k] = singleParDict[key]

        # Now, account for any user-specified reference image
        def_wcs = default_wcs.deepcopy()
        single_ref = singleParDict[keyname + 'refimage']
        if single_ref:
            if isinstance(single_ref, wcs.WCS):
                default_wcs = single_ref
            else:
                default_wcs = wcsutil.HSTWCS(singleParDict[keyname + 'refimage'])

        # ## Create single_wcs instance based on user parameters
        outwcs.single_wcs = mergeWCS(default_wcs, single_pars)
        # restore global default WCS to original value so single_drizzle WCS does not
        # influence final_drizzle WCS
        default_wcs = def_wcs.deepcopy()

    final_step = configObj[util.getSectionName(configObj, 7)]
    finalParDict = configObj[util.getSectionName(configObj, '7a')].copy()
    if final_step['driz_combine'] and finalParDict['final_wcs']:
        del finalParDict['final_wcs']
        keyname = 'final_'
        for key in finalParDict:
            k = key[len(keyname):]
            if k != 'refimage':
                final_pars[k] = finalParDict[key]

        # Now, account for any user-specified reference image
        final_ref = finalParDict[keyname + 'refimage']
        if final_ref:
            if isinstance(final_ref, wcs.WCS):
                default_wcs = final_ref
                if hasattr(final_ref, 'filename'):
                    rootname = final_ref.filename
                else:
                    rootname = ""
                print('Creating OUTPUT WCS from WCS object based on {}'.format(rootname))
            else:
                rootname, extnum = fileutil.parseFilename(finalParDict[keyname + 'refimage'])
                extnum = util.findWCSExtn(finalParDict[keyname + 'refimage'])
                print('Creating OUTPUT WCS from {}[{}]'.format(rootname, extnum))
                default_wcs = wcsutil.HSTWCS('{}[{}]'.format(rootname, extnum))

        # ## Create single_wcs instance based on user parameters
        outwcs.final_wcs = mergeWCS(default_wcs, final_pars)
        outwcs.wcs = outwcs.final_wcs.copy()

    # Apply user settings to create custom output_wcs instances
    # for each drizzle step
    updateImageWCS(imageObjectList, outwcs)

    return outwcs


def make_perfect_cd(wcs):
    """ Create a perfect (square, orthogonal, undistorted) CD matrix from the
        input WCS.
    """
    def_scale = (wcs.pscale) / 3600.
    def_orientat = np.deg2rad(wcs.orientat)
    perfect_cd = def_scale * np.array(
        [[-np.cos(def_orientat),np.sin(def_orientat)],
         [np.sin(def_orientat),np.cos(def_orientat)]]
    )
    return perfect_cd


def calcNewEdges(wcs, shape):
    """
    This method will compute sky coordinates for all the pixels around
    the edge of an image AFTER applying the geometry model.

    Parameters
    ----------
    wcs : obj
        HSTWCS object for image

    shape : tuple
        numpy shape tuple for size of image

    Returns
    -------
    border : arr
        array which contains the new positions for
        all pixels around the border of the edges in alpha,dec

    """

    naxis1 = shape[1]
    naxis2 = shape[0]
    # build up arrays for pixel positions for the edges
    # These arrays need to be: array([(x,y),(x1,y1),...])
    numpix = naxis1*2 + naxis2*2
    border = np.zeros(shape=(numpix,2),dtype=np.float64)

    # Now determine the appropriate values for this array
    # We also need to account for any subarray offsets
    xmin = 1.
    xmax = naxis1
    ymin = 1.
    ymax = naxis2

    # Build range of pixel values for each side
    # Add 1 to make them consistent with pixel numbering in IRAF
    # Also include the LTV offsets to represent position in full chip
    #   since the model works relative to full chip positions.
    xside = np.arange(naxis1) + xmin
    yside = np.arange(naxis2) + ymin

    #Now apply them to the array to generate the appropriate tuples
    #bottom
    _range0 = 0
    _range1 = naxis1
    border[_range0:_range1,0] = xside
    border[_range0:_range1,1] = ymin
    #top
    _range0 = _range1
    _range1 = _range0 + naxis1
    border[_range0:_range1,0] = xside
    border[_range0:_range1,1] = ymax
    #left
    _range0 = _range1
    _range1 = _range0 + naxis2
    border[_range0:_range1,0] = xmin
    border[_range0:_range1,1] = yside
    #right
    _range0 = _range1
    _range1 = _range0 + naxis2
    border[_range0:_range1,0] = xmax
    border[_range0:_range1,1] = yside

    edges = wcs.all_pix2world(border[:,0],border[:,1],1)
    return edges


def computeEdgesCenter(edges):
    alpha = np.deg2rad(edges[0])
    dec = np.deg2rad(edges[1])

    xmean = np.mean(np.cos(dec)*np.cos(alpha))
    ymean = np.mean(np.cos(dec)*np.sin(alpha))
    zmean = np.mean(np.sin(dec))

    crval1 = np.rad2deg(np.arctan2(ymean,xmean))%360.0
    crval2 = np.rad2deg(np.arctan2(zmean,np.sqrt(xmean*xmean+ymean*ymean)))

    return crval1,crval2


#### Utility functions for working with WCSObjects
def createWCSObject(output,default_wcs,imageObjectList):
    """Converts a astropy.wcs WCS object into a WCSObject(baseImageObject) instance."""
    from . import imageObject
    outwcs = imageObject.WCSObject(output)
    outwcs.default_wcs = default_wcs
    outwcs.wcs = default_wcs.copy()
    outwcs.final_wcs = default_wcs.copy()
    outwcs.single_wcs = default_wcs.copy()

    outwcs.updateContextImage(imageObjectList[0].createContext)

    #
    # Add exptime information for use with drizzle
    #
    outwcs._exptime,outwcs._expstart,outwcs._expend = util.compute_texptime(imageObjectList)

    outwcs.nimages = util.countImages(imageObjectList)

    return outwcs

def removeAllAltWCS(hdulist,extlist):
    """
    Removes all alternate WCS solutions from the header
    """
    original_logging_level = log.level
    log.setLevel(logutil.logging.WARNING)

    try:
        hdr = hdulist[extlist[0]].header
        wkeys = altwcs.wcskeys(hdr)
        if ' ' in wkeys:
            wkeys.remove(' ')
        for extn in extlist:
            for wkey in wkeys:
                if wkey == 'O':
                    continue
                altwcs.deleteWCS(hdulist,extn,wkey)

            # Forcibly remove OPUS WCS Keywords, since deleteWCS will not do it
            hwcs = readAltWCS(hdulist,extn,wcskey='O')

            if hwcs is None:
                continue

            for k in hwcs.keys():
                if k not in ['DATE-OBS','MJD-OBS'] and k in hdr:
                    try:
                        del hdr[k]
                    except KeyError:
                        pass
    except:
        raise

    finally:
        log.setLevel(original_logging_level) # restore original logging level


def updateImageWCS(imageObjectList, output_wcs):

     # Update input imageObjects with output WCS information
    for img in imageObjectList:
        img.updateOutputValues(output_wcs)


def restoreDefaultWCS(imageObjectList, output_wcs):
    """ Restore WCS information to default values, and update imageObject
        accordingly.
    """
    if not isinstance(imageObjectList,list):
        imageObjectList = [imageObjectList]

    output_wcs.restoreWCS()

    updateImageWCS(imageObjectList, output_wcs)


def _rotateCD(cd, theta):
    # This is the old (pre-astropy PR #5189 -
    # see https://github.com/astropy/astropy/pull/5189) as this the version
    # of the rotateCD expected in 'mergeWCS' below
    theta = np.deg2rad(theta)
    cth = np.cos(theta)
    sth = np.sin(theta)
    return np.dot(cd, [[cth, sth], [-sth, cth]])


def _py2round(x):
    """
    This function returns a rounded up value of the argument, similar
    to Python 2.
    """
    if hasattr(x, '__iter__'):
        rx = np.empty_like(x)
        m = x >= 0.0
        rx[m] = np.floor(x[m] + 0.5)
        m = np.logical_not(m)
        rx[m] = np.ceil(x[m] - 0.5)
        return rx

    else:
        if x >= 0.0:
            return np.floor(x + 0.5)
        else:
            return np.ceil(x - 0.5)


def _check_custom_WCS_pars(par1_name, par2_name, user_pars):
    par1 = par1_name in user_pars and user_pars[par1_name] is not None
    par2 = par2_name in user_pars and user_pars[par2_name] is not None

    if par1 != par2:
        if par1:
            par1n = par1_name
            par2n = par2_name
        else:
            par1n = par2_name
            par2n = par1_name
        raise ValueError("When WCS parameter '{}' is specified, '{}' must "
                         "be specified as well.".format(par1n, par2n))

    return par1


def _check_close_scale(scale, ref):
    rtol = 10.0 * np.finfo(float).eps
    atol = 100.0 * np.finfo(float).tiny
    return np.isclose(scale, ref, atol=atol, rtol=rtol)


def mergeWCS(default_wcs, user_pars):
    """ Merges the user specified WCS values given as dictionary derived from
        the input configObj object with the output astropy.wcs object computed
        using distortion.output_wcs().

        The user_pars dictionary needs to have the following set of keys::

            user_pars = {'ra':None,'dec':None,'scale':None,'rot':None,
                         'outnx':None,'outny':None,'crpix1':None,'crpix2':None}
    """
    #
    # Start by making a copy of the input WCS...
    #
    outwcs = default_wcs.deepcopy()

    # If there are no user set parameters, just return a copy of
    # the original WCS:
    if all([upar is None for upar in user_pars.values()]):
        return outwcs

    if _check_custom_WCS_pars('ra', 'dec', user_pars):
        _crval = (user_pars['ra'], user_pars['dec'])
    else:
        _crval = None

    if ('scale' in user_pars and user_pars['scale'] is not None and
        not _check_close_scale(user_pars['scale'], outwcs.pscale)):
        _scale = user_pars['scale']
        _ratio = outwcs.pscale / _scale
    else:
        _ratio = None
        _scale = None

    if ('rot' not in user_pars) or user_pars['rot'] is None:
        _delta_rot = None
    else:
        _delta_rot = outwcs.orientat - user_pars['rot']
        if _delta_rot == 0.0:
            _delta_rot = None

    if _check_custom_WCS_pars('crpix1', 'crpix2', user_pars):
        _crpix = (user_pars['crpix1'], user_pars['crpix2'])
    else:
        _crpix = None

    shape = None

    if _check_custom_WCS_pars('outnx', 'outny', user_pars):
        shape = (
            int(_py2round(user_pars['outnx'])),
            int(_py2round(user_pars['outny']))
        )

        if shape[0] < 1 or shape[1] < 1:
            raise ValueError("Custom WCS output image size smaller than 1")

        if _crpix is None:
            # make sure new image is centered on the CRPIX of the old WCS:
            _crpix = ((shape[0] + 1.0) / 2.0, (shape[1] + 1.0) / 2.0)

    else:
        naxis1, naxis2 = outwcs.pixel_shape
        if _delta_rot is None:
            # no rotation is involved

            if _ratio is not None:
                # apply scale only:

                # compute output image shape:
                shape = (
                    max(1, int(np.ceil(_ratio * naxis1))),
                    max(1, int(np.ceil(_ratio * naxis2)))
                )

                # update CRPIX:
                if _crpix is None:
                    _crpix = 1.0 + _ratio * (outwcs.wcs.crpix - 1.0)

        else:
            _corners = np.array(
                [[0.5, 0.5],
                 [naxis1 + 0.5, 0.5],
                 [0.5, naxis2 + 0.5],
                 [naxis1 + 0.5, naxis2 + 0.5]]
            ) - outwcs.wcs.crpix

            if _ratio is not None:
                # scale corners:
                _corners *= _ratio

            # rotate corners and find new image range:
            ((_xmin, _xmax), (_ymin, _ymax)) = util.getRotatedSize(_corners,
                                                                   _delta_rot)

            # compute output image shape:
            # NOTE: _py2round may be replaced with np.ceil
            shape = (
                max(1, int(np.ceil(_xmax - _xmin))),
                max(1, int(np.ceil(_ymax - _ymin)))
            )

            if _crpix is None:
                # update CRPIX:
                _crpix = (-_xmin + 0.5, -_ymin + 0.5)


    # Set up the new WCS based on values from old one:

    if _ratio is not None:
        # Update plate scale
        outwcs.wcs.cd = outwcs.wcs.cd / _ratio
        outwcs.pscale = _scale

    # update orientation
    if _delta_rot is not None:
        outwcs.wcs.cd = _rotateCD(outwcs.wcs.cd, _delta_rot)
        outwcs.orientat -= _delta_rot

    if shape is not None:
        # update size:
        outwcs.pixel_shape = shape

    # update reference position
    if _crpix is not None:
        outwcs.wcs.crpix = np.array(_crpix, dtype=np.float64)

    if _crval is not None:
        outwcs.wcs.crval = np.array(_crval, dtype=np.float64)

    return outwcs


def convertWCS(inwcs,drizwcs):
    """ Copy WCSObject WCS into Drizzle compatible array."""
    drizwcs[0] = inwcs.crpix[0]
    drizwcs[1] = inwcs.crval[0]
    drizwcs[2] = inwcs.crpix[1]
    drizwcs[3] = inwcs.crval[1]
    drizwcs[4] = inwcs.cd[0][0]
    drizwcs[5] = inwcs.cd[1][0]
    drizwcs[6] = inwcs.cd[0][1]
    drizwcs[7] = inwcs.cd[1][1]

    return drizwcs

def updateWCS(drizwcs,inwcs):
    """ Copy output WCS array from Drizzle into WCSObject."""
    crpix = np.array([drizwcs[0],drizwcs[2]], dtype=np.float64)
    crval = np.array([drizwcs[1],drizwcs[3]], dtype=np.float64)
    cd = np.array([[drizwcs[4],drizwcs[6]],[drizwcs[5],drizwcs[7]]], dtype=np.float64)
    inwcs.cd = cd
    inwcs.crval = crval
    inwc.crpix = crpix
    inwcs.pscale = N.sqrt(N.power(inwcs.cd[0][0],2)+N.power(inwcs.cd[1][0],2)) * 3600.
    inwcs.orient = N.arctan2(inwcs.cd[0][1],inwcs.cd[1][1]) * 180./N.pi


def wcsfit(img_wcs, ref_wcs):
    """
    Perform a linear fit between 2 WCS for shift, rotation and scale.
    Based on the WCSLIN function from 'drutil.f'(Drizzle V2.9) and modified to
    allow for differences in reference positions assumed by PyDrizzle's
    distortion model and the coeffs used by 'drizzle'.

    Parameters
    ----------
    img  : obj
        ObsGeometry instance for input image

    ref_wcs : obj
        Undistorted WCSObject instance for output frame

    """
    # Define objects that we need to use for the fit...
    #in_refpix = img_geom.model.refpix
    wmap = WCSMap(img_wcs,ref_wcs)
    cx, cy = coeff_converter.sip2idc(img_wcs)
    # Convert the RA/Dec positions back to X/Y in output product image
    #_cpix_xyref = np.zeros((4,2),dtype=np.float64)

    # Start by setting up an array of points +/-0.5 pixels around CRVAL1,2
    # However, we must shift these positions by 1.0pix to match what
    # drizzle will use as its reference position for 'align=center'.
    _cpix = (img_wcs.wcs.crpix[0],img_wcs.wcs.crpix[1])
    _cpix_arr = np.array([_cpix,(_cpix[0],_cpix[1]+1.),
                       (_cpix[0]+1.,_cpix[1]+1.),(_cpix[0]+1.,_cpix[1])], dtype=np.float64)
    # Convert these positions to RA/Dec
    _cpix_rd = wmap.xy2rd(img_wcs,_cpix_arr[:,0],_cpix_arr[:,1])
    #for pix in xrange(len(_cpix_rd[0])):
    _cpix_xref,_cpix_yref = wmap.rd2xy(ref_wcs,_cpix_rd[0],_cpix_rd[1])
    _cpix_xyref = np.zeros((4,2),dtype=np.float64)
    _cpix_xyref[:,0] = _cpix_xref
    _cpix_xyref[:,1] = _cpix_yref

    """
    # needed to handle correctly subarrays and wfpc2 data
    if img_wcs.delta_refx == 0.0 and img_wcs.delta_refy == 0.0:
        offx, offy = (0.0,0.0)
    else:
        offx, offy = (1.0, 1.0)
    """
    offx, offy = (0.0,0.0)

    # Now, apply distortion model to input image XY positions
    #_cpix_xyc = np.zeros((4,2),dtype=np.float64)
    _cpix_xyc = utils.apply_idc(_cpix_arr, cx, cy, img_wcs.wcs.crpix, img_wcs.pscale, order=1)

    # Need to get the XDELTA,YDELTA values included here in order to get this
    # to work with MDTng.
    #if in_refpix:
    #    _cpix_xyc += (in_refpix['XDELTA'], in_refpix['YDELTA'])

    # Perform a fit between:
    #       - undistorted, input positions: _cpix_xyc
    #       - X/Y positions in reference frame: _cpix_xyref
    abxt,cdyt = fitlin(_cpix_xyc,_cpix_xyref)

    # This correction affects the final fit when you are fitting
    # a WCS to itself (no distortion coeffs), so it needs to be
    # taken out in the coeffs file by modifying the zero-point value.
    #  WJH 17-Mar-2005
    abxt[2] -= ref_wcs.wcs.crpix[0] + offx
    cdyt[2] -= ref_wcs.wcs.crpix[1] + offy

    return abxt,cdyt


def fitlin(imgarr,refarr):
    """ Compute the least-squares fit between two arrays.
        A Python translation of 'FITLIN' from 'drutil.f' (Drizzle V2.9).
    """
    # Initialize variables
    _mat = np.zeros((3,3),dtype=np.float64)
    _xorg = imgarr[0][0]
    _yorg = imgarr[0][1]
    _xoorg = refarr[0][0]
    _yoorg = refarr[0][1]
    _sigxox = 0.
    _sigxoy = 0.
    _sigxo = 0.
    _sigyox = 0.
    _sigyoy = 0.
    _sigyo = 0.

    _npos = len(imgarr)
    # Populate matrices
    for i in range(_npos):
        _mat[0][0] += np.power((imgarr[i][0] - _xorg),2)
        _mat[0][1] += (imgarr[i][0] - _xorg) * (imgarr[i][1] - _yorg)
        _mat[0][2] += (imgarr[i][0] - _xorg)
        _mat[1][1] += np.power((imgarr[i][1] - _yorg),2)
        _mat[1][2] += imgarr[i][1] - _yorg

        _sigxox += (refarr[i][0] - _xoorg)*(imgarr[i][0] - _xorg)
        _sigxoy += (refarr[i][0] - _xoorg)*(imgarr[i][1] - _yorg)
        _sigxo += refarr[i][0] - _xoorg
        _sigyox += (refarr[i][1] - _yoorg)*(imgarr[i][0] -_xorg)
        _sigyoy += (refarr[i][1] - _yoorg)*(imgarr[i][1] - _yorg)
        _sigyo += refarr[i][1] - _yoorg

    _mat[2][2] = _npos
    _mat[1][0] = _mat[0][1]
    _mat[2][0] = _mat[0][2]
    _mat[2][1] = _mat[1][2]

    # Now invert this matrix
    _mat = linalg.inv(_mat)

    _a  = _sigxox*_mat[0][0]+_sigxoy*_mat[0][1]+_sigxo*_mat[0][2]
    _b  = -1*(_sigxox*_mat[1][0]+_sigxoy*_mat[1][1]+_sigxo*_mat[1][2])
    #_x0 = _sigxox*_mat[2][0]+_sigxoy*_mat[2][1]+_sigxo*_mat[2][2]

    _c  = _sigyox*_mat[1][0]+_sigyoy*_mat[1][1]+_sigyo*_mat[1][2]
    _d  = _sigyox*_mat[0][0]+_sigyoy*_mat[0][1]+_sigyo*_mat[0][2]
    #_y0 = _sigyox*_mat[2][0]+_sigyoy*_mat[2][1]+_sigyo*_mat[2][2]

    _xt = _xoorg - _a*_xorg+_b*_yorg
    _yt = _yoorg - _d*_xorg-_c*_yorg

    return [_a,_b,_xt],[_c,_d,_yt]


def fitlin_rscale(xy,uv,verbose=False):
    """ Performs a linear, orthogonal fit between matched
        lists of positions 'xy' (input) and 'uv' (output).

        Output: (same as for fit_arrays_general)
    """
    mu = uv[:,0].mean()
    mv = uv[:,1].mean()
    mx = xy[:,0].mean()
    my = xy[:,1].mean()

    u = uv[:,0] - mu
    v = uv[:,1] - mv
    x = xy[:,0] - mx
    y = xy[:,1] - my

    Sxx = np.dot(x,x)
    Syy = np.dot(y,y)
    Sux = np.dot(u,x)
    Suy = np.dot(u,y)
    Svx = np.dot(v,x)
    Svy = np.dot(v,y)

    # implement parity check
    if (np.dot(Sux,Svy) > 0):
        p = 1
    else:
        p = -1

    XX = p*Sux + Svy
    YY = Suy - p*Svx

    # derive output values
    theta_deg = np.rad2deg(np.arctan2(YY,XX))% 360.0
    scale = np.sqrt(XX**2 + YY**2) / (Sxx+Syy)
    shift = (mu-mx,mv-my)
    if verbose:
        print('Linear RSCALE fit: rotation = ',theta_deg,'  scale = ',scale,'  offset = ',shift)
    coeffs = scale * fileutil.buildRotMatrix(-theta_deg)

    P = [coeffs[0,0],coeffs[0,1],shift[0]]
    Q = [coeffs[1,1],coeffs[1,0],shift[1]]
    return P,Q

def fitlin_clipped(xy,uv,verbose=False,mode='rscale',nclip=3,reject=3):
    """ Perform a clipped fit based on the number of iterations and rejection limit
        (in sigma) specified by the user. This will more closely replicate the results
        obtained by 'geomap' using 'maxiter' and 'reject' parameters.
    """
    fitting_funcs = {'rscale':fitlin_rscale,'general':fitlin}
    # Get the fitting function to be used
    fit_func = fitting_funcs[mode.lower()]
    # Perform the initial fit
    P,Q = fit_func(xy,uv)
    xyc = apply_fitlin(xy,P,Q)

    # compute residuals from fit for input positions
    dx = uv[:,0] - xyc[0]
    dy = uv[:,1] - xyc[1]
    fit_rms = [dx.std(),dy.std()]

    if nclip > 0:
        data = xy.copy()
        outdata = uv.copy()

    numclipped = 0
    for i in range(nclip):
        iterclipped = 0
        xyc = apply_fitlin(data,P,Q)

        # compute residuals from fit for input positions
        dx = outdata[:,0] - xyc[0]
        dy = outdata[:,1] - xyc[1]

        # find indices of outliers in x and y
        xout = np.where(np.abs(dx - dx.mean()) > reject*dx.std())
        yout = np.where(np.abs(dy - dy.mean()) > reject*dy.std())
        # concatenate those indices and sort them
        outliers_indx = xout[0].tolist()+yout[0].tolist()
        outliers_indx.sort()
        # define the full range of indices for the data points left
        full_indx = list(range(data.shape[0]))
        # remove all unique indices specified in outliers from full range
        for o in outliers_indx:
            # only remove if it has not been removed already
            # accounts for the same point being an outlier in both x and y
            if full_indx.count(o) > 0:
                full_indx.remove(o)
            iterclipped += 1

        if iterclipped == 0:
            break

        numclipped += iterclipped
        if verbose:
            print('Removed a total of ',numclipped,' points through iteration ',i+1)
        # create clipped data
        data_iter = np.zeros([len(full_indx),2],dtype=data.dtype)
        if verbose:
            print('Iter #',i+1,' data:',data.shape,data_iter.shape,len(full_indx))

        data_iter[:,0] = data[:,0][full_indx]
        data_iter[:,1] = data[:,1][full_indx]
        outdata_iter = np.zeros([len(full_indx),2],dtype=data.dtype)
        outdata_iter[:,0] = outdata[:,0][full_indx]
        outdata_iter[:,1] = outdata[:,1][full_indx]

        # perform the fit again with the clipped data and go to the next iteration
        data = data_iter
        outdata = outdata_iter
        P,Q = fit_func(data,outdata)
        # compute residuals from fit for input positions
        xyc = apply_fitlin(data,P,Q)
        dx = outdata[:,0] - xyc[0]
        dy = outdata[:,1] - xyc[1]
        fit_rms = [dx.std(),dy.std()]

    if verbose:
        print('Fit clipped ',numclipped,' points over ',nclip,' iterations.')
    return P,Q,fit_rms

def apply_fitlin(data,P,Q):
    # start by parsing coefficients into affine matrix and offset
    fit = np.array([[P[0],P[1]],[Q[1],Q[0]]])
    xsh = P[2]
    ysh = Q[2]

    if fit is not None:
        xy1 = np.dot(data,fit)
        xy1x = xy1[:,0] + xsh
        xy1y = xy1[:,1] + ysh
    return xy1x,xy1y

def readAltWCS(fobj, ext, wcskey=' ', verbose=False):
    """
    Reads in alternate primary WCS from specified extension.

    Parameters
    ----------
    fobj : str, `astropy.io.fits.HDUList`
        fits filename or fits file object
        containing alternate/primary WCS(s) to be converted
    wcskey : str
        [" ",A-Z]
        alternate/primary WCS key that will be replaced by the new key
    ext : int
        fits extension number
    Returns
    -------
    hdr: fits.Header
        header object with ONLY the keywords for specified alternate WCS
    """
    if isinstance(fobj, str):
        fobj = fits.open(fobj, memmap=False)

    hdr = altwcs._getheader(fobj, ext)
    try:
        original_logging_level = log.level
        log.setLevel(logutil.logging.WARNING)
        nwcs = wcs.WCS(hdr, fobj=fobj, key=wcskey)

    except KeyError:
        if verbose:
            print('readAltWCS: Could not read WCS with key %s' % wcskey)
            print('            Skipping %s[%s]' % (fobj.filename(), str(ext)))
        return None

    finally:
        log.setLevel(original_logging_level) # restore original logging level

    hwcs = nwcs.to_header()

    if nwcs.wcs.has_cd():
        hwcs = altwcs.pc2cd(hwcs, key=wcskey)
    return hwcs

def make_mosaic_wcs(filenames, rot=None, scale=None):
    """Combine WCSs from all input files into a single meta-WCS

    Parameters
    ----------
    filenames : list
        list of filenames to process

    rot : float
        desired rotation of meta-WCS. Default value is None.

    scale : float
        desired scale of meta-WCS. Default value is None.

    Returns
    --------
    mosaic_wcs : HSTWCS object
        the merged composite WCS
    """
    if not isinstance(filenames, list):
        filenames = [filenames]

    # Compile list of WCSs for all chips from all input filenames
    hstwcs_list = []
    for f in filenames:
        with pyfits.open(f) as hdulist:
            hstwcs_list.extend([get_hstwcs(f, hdulist, extnum) for extnum in get_extns(f)])
    # Combine them into a single mosaic WCS
    output_wcs = utils.output_wcs(hstwcs_list, undistort=True)
    output_wcs.wcs.cd = make_perfect_cd(output_wcs)

    mosaic_wcs = mergeWCS(output_wcs, {'rot': rot, 'scale': scale})
    return mosaic_wcs

def create_mosaic_pars(mosaic_wcs):
    """Return dict of values to use with AstroDrizzle to specify output mosaic

    Parameters
    -----------
    mosaic_wcs : object
        WCS as generated by make_mosaic_wcs

    Returns
    --------
    mosaic_pars : dict
        This dictionary can be used as input to astrodrizzle.AstroDrizzle with
        the syntax 'AstroDrizzle(filenames, **mosaic_pars)'

    """
    mosaic_pars = dict(
        driz_sep_rot=mosaic_wcs.orientat,
        driz_sep_scale=mosaic_wcs.pscale,
        driz_sep_outnx=mosaic_wcs.array_shape[1],
        driz_sep_outny=mosaic_wcs.array_shape[0],
        driz_sep_ra=mosaic_wcs.wcs.crval[0],
        driz_sep_dec=mosaic_wcs.wcs.crval[1],
        driz_sep_crpix1=mosaic_wcs.wcs.crpix[0],
        driz_sep_crpix2=mosaic_wcs.wcs.crpix[1],
        final_rot=mosaic_wcs.orientat,
        final_scale=mosaic_wcs.pscale,
        final_outnx=mosaic_wcs.array_shape[1],
        final_outny=mosaic_wcs.array_shape[0],
        final_ra=mosaic_wcs.wcs.crval[0],
        final_dec=mosaic_wcs.wcs.crval[1],
        final_crpix1=mosaic_wcs.wcs.crpix[0],
        final_crpix2=mosaic_wcs.wcs.crpix[1]
        )
    return mosaic_pars

def get_extns(fimg, extname='SCI'):
    """
    Examines the specified fits image and returns the extension number(s) whose
     name(s) match the specified by input value.

    Parameters
    -----------
    fimg : string
        Name of the image to probe
    extname : string
        Extension name to search for in specified image. Default value is 'SCI'.

    Returns
    --------
    extns : list
        List of matching extension numbers

    """
    extns = []
    close_img = False
    if isinstance(fimg, str):
        fimg = pyfits.open(fimg)
        close_img = True

    for i,e in enumerate(fimg):
        if 'extname' in e.header and e.header['extname'] == extname:
            extns.append(i)

    if close_img:
        fimg.close()

    return extns
