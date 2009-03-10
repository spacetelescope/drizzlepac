import numpy as np
from numpy import linalg

from updatewcs import wcsutil
from updatewcs.distortion import utils

from pytools import fileutil
import util
import imageObject
import updatewcs
from updatewcs import pywcs
from updatewcs import distortion
from updatewcs.distortion import coeff_converter,utils


DEFAULT_WCS_PARS = {'ra':None,'dec':None,'psize':None,'orient':None,
                     'outnx':None,'outny':None,'crpix1':None,'crpix2':None,
                     'crval1':None,'crval2':None}


# Default mapping function based on PyWCS 
class WCSMap:
    def __init__(self,input,output):
        # Verify that we have valid WCS input objects
        self.checkWCS(input,'Input')
        self.checkWCS(output,'Output')

        self.input = input
        self.output = output
        

    def checkWCS(self,obj,name):
        try:
            assert isinstance(obj, pywcs.WCS)
        except AssertionError:
            print name +' object needs to be an instance or subclass of a PyWCS object.'
            raise
    def forward(self,pixx,pixy,origin=1):
        # This matches WTRAXY results to better than 1e-4 pixels.
        skyx,skyy = self.input.all_pix2sky(pixx,pixy,origin)
        return self.output.wcs_sky2pix(skyx,skyy,origin) 

    def xy2rd(self,wcs,pixx,pixy):
        return wcs.wcs_pix2sky(pixx,pixy,1)
    def rd2xy(self,wcs,ra,dec):
        return wcs.wcs_sky2pix(ra,dec,1)
    
    def backward(self,pixx,pixy):
        # Not needed for drizzle or blot.
        skyx,skyy = self.output.all_pix2sky(pixx,pixy,1)
        return self.input.wcs_sky2pix(skyx,skyy,1)
    
                    
def get_hstwcs(filename,hdulist,extnum):
    ''' Return the HSTWCS object for a given chip.
    
    '''
    hdrwcs = wcsutil.HSTWCS(hdulist,ext=extnum)
    hdrwcs.filename = filename
    hdrwcs.expname = hdulist[extnum].header['expname']
    hdrwcs.extver = hdulist[extnum].header['extver']
    
    return hdrwcs

#
# Possibly need to generate a stand-alone interface for this function.
#
# Primary interface for creating the output WCS from a list of HSTWCS objects
def make_outputwcs(imageObjectList,output,configObj=None):
    """ Computes the full output WCS based on the set of input imageObjects
        provided as input, along with the pre-determined output name from
        process_input.  The user specified output parameters are then used to
        modify the default WCS to produce the final desired output frame.
        The input imageObjectList has the outputValues dictionary
        updated with the information from the computed output WCS. 
        It then returns this WCS as a WCSObject(imageObject) 
        instance.
    """
    if not isinstance(imageObjectList,list): 
        imageObjectList = [imageObjectList]
        
    if configObj['refimage'].strip() in ['',None]:        
        # Compute default output WCS, if no refimage specified
        hstwcs_list = []
        for img in imageObjectList:
            hstwcs_list += img.getKeywordList('wcs')
        default_wcs = utils.output_wcs(hstwcs_list)
    else:
        # Otherwise, simply use the reference image specified by the user
        default_wcs = wcsutil.HSTWCS(configObj['refimage'])

    # Turn WCS instances into WCSObject instances
    outwcs = createWCSObject(output,default_wcs,default_wcs,imageObjectList)
    
    # Merge in user-specified attributes for the output WCS
    # as recorded in the input configObj object.
    final_pars = DEFAULT_WCS_PARS.copy()
         
    # More interpretation of the configObj needs to be done here to translate
    # the input parameter names to those understood by 'mergeWCS' as defined
    # by the DEFAULT_WCS_PARS dictionary.
    single_step = util.getSectionName(configObj,3)
    if single_step and configObj[single_step]['driz_separate']: 
        single_pars = DEFAULT_WCS_PARS.copy()
        single_pars['ra'] = configObj['ra']
        single_pars['dec'] = configObj['dec']
        #single_pars.update(configObj['STEP 3: DRIZZLE SEPARATE IMAGES'])
        single_keys = {'outnx':'driz_sep_outnx','outny':'driz_sep_outny',
                        'rot':'driz_sep_rot', 'scale':'driz_sep_scale'}
        for key in single_keys.keys():
            single_pars[key] = configObj['STEP 3: DRIZZLE SEPARATE IMAGES'][single_keys[key]]
        ### Create single_wcs instance based on user parameters
        outwcs.single_wcs = mergeWCS(default_wcs,single_pars)
        

    final_step = util.getSectionName(configObj,7)
    if final_step and configObj[final_step]['driz_combine']: 
        final_pars = DEFAULT_WCS_PARS.copy()
        final_pars['ra'] = configObj['ra']
        final_pars['dec'] = configObj['dec']
        final_keys = {'outnx':'final_outnx','outny':'final_outny','rot':'final_rot', 'scale':'final_scale'}
        #final_pars.update(configObj['STEP 7: DRIZZLE FINAL COMBINED IMAGE'])
        for key in final_keys.keys():
            final_pars[key] = configObj['STEP 7: DRIZZLE FINAL COMBINED IMAGE'][final_keys[key]]
        ### Create single_wcs instance based on user parameters
        outwcs.final_wcs = mergeWCS(default_wcs,final_pars)
        outwcs.wcs = outwcs.final_wcs.copy()

    # Apply user settings to create custom output_wcs instances 
    # for each drizzle step
    updateImageWCS(imageObjectList,outwcs)
    
    return outwcs


def createWCSObject(output,default_wcs,final_wcs,imageObjectList):
    """Converts a PyWCS WCS object into a WCSObject(baseImageObject) instance."""
    outwcs = imageObject.WCSObject(output)
    outwcs.default_wcs = default_wcs
    outwcs.wcs = final_wcs

    #
    # Add exptime information for use with drizzle
    #
    outwcs._exptime,outwcs._expstart,outwcs._expend = util.compute_texptime(imageObjectList)
        
    outwcs.nimages = util.countImages(imageObjectList)
     
    return outwcs

def updateImageWCS(imageObjectList,output_wcs):
    
     # Update input imageObjects with output WCS information
    for img in imageObjectList:
        img.updateOutputValues(output_wcs)
   
def restoreDefaultWCS(imageObjectList,output_wcs):
    """ Restore WCS information to default values, and update imageObject
        accordingly.
    """
    if not isinstance(imageObjectList,list): 
        imageObjectList = [imageObjectList]

    output_wcs.restoreWCS()
    
    updateImageWCS(imageObjectList,output_wcs)

def mergeWCS(default_wcs,user_pars):
    """ Merges the user specified WCS values given as dictionary derived from 
        the input configObj object with the output PyWCS object computed 
        using distortion.output_wcs().
        
        The user_pars dictionary needs to have the following set of keys:
        user_pars = {'ra':None,'dec':None,'psize':None,'orient':None,
                     'outnx':None,'outny':None,'crpix1':None,'crpix2':None,
                     'crval1':None,'crval2':None}
    """
    #
    # Start by making a copy of the input WCS...
    #    
    outwcs = default_wcs.copy()    

    # If there are no user set parameters, just return a copy of the original WCS
    if user_pars == DEFAULT_WCS_PARS:
        return outwcs

    if (not user_pars.has_key('ra')) or user_pars['ra'] == None:
        _crval = None
    else:
        _crval = (user_pars['ra'],user_pars['dec'])

    if (not user_pars.has_key('psize')) or user_pars['psize'] == None:
        _ratio = 1.0
        _psize = None
        # Need to resize the WCS for any changes in pscale
    else:
        _ratio = outwcs.pscale / user_pars['psize']
        _psize = user_pars['psize']

    if (not user_pars.has_key('orient')) or user_pars['orient'] == None:
        _orient = None
        _delta_rot = 0.
    else:
        _orient = user_pars['orient']
        _delta_rot = outwcs.orientat - user_pars['orient']

    _mrot = fileutil.buildRotMatrix(_delta_rot)

    if (not user_pars.has_key('outnx')) or user_pars['outnx'] == None:
        _corners = np.array([[0.,0.],[outwcs.naxis1,0.],[0.,outwcs.naxis2],[outwcs.naxis1,outwcs.naxis2]])
        _corners -= (outwcs.naxis1/2.,outwcs.naxis2/2.)
        _range = util.getRotatedSize(_corners,_delta_rot)
        shape = ((_range[0][1] - _range[0][0])*_ratio,(_range[1][1]-_range[1][0])*_ratio)
        old_shape = (outwcs.naxis1*_ratio,outwcs.naxis2*_ratio)

        _crpix = (shape[0]/2., shape[1]/2.)

    else:
        shape = [user_pars['outnx'],user_pars['outny']]
        if user_pars['crpix1'] == None:
            _crpix = (shape[0]/2.,shape[1]/2.)
        else:
            _crpix = [user_pars['crpix1'],user_pars['crpix2']]

    # Set up the new WCS based on values from old one.
    # Update plate scale
    outwcs.wcs.cd *= _ratio
    outwcs.pscale /= _ratio
    #Update orientation
    outwcs.rotateCD(_delta_rot)
    outwcs.orientat += -_delta_rot
    # Update size
    outwcs.naxis1 =  int(shape[0])
    outwcs.naxis2 =  int(shape[1])
    # Update reference position
    outwcs.wcs.crpix =_crpix
    if _crval is not None:
        outwcs.wcs.crval = _crval

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
    inwcs.crpix[0]    = drizwcs[0]
    inwcs.crval[0]   = drizwcs[1]
    inwcs.crpix[1]   = drizwcs[2]
    inwcs.crval[1]   = drizwcs[3]
    inwcs.cd[0][0]     = drizwcs[4]
    inwcs.cd[1][0]     = drizwcs[5]
    inwcs.cd[0][1]     = drizwcs[6]
    inwcs.cd[1][1]     = drizwcs[7]
    inwcs.pscale = N.sqrt(N.power(inwcs.cd[0][0],2)+N.power(inwcs.cd[1][0],2)) * 3600.
    inwcs.orient = N.arctan2(inwcs.cd[0][1],inwcs.cd[1][1]) * 180./N.pi


def wcsfit(img_wcs, ref_wcs):
    """
    Perform a linear fit between 2 WCS for shift, rotation and scale.
    Based on 'WCSLIN' from 'drutil.f'(Drizzle V2.9) and modified to
    allow for differences in reference positions assumed by PyDrizzle's
    distortion model and the coeffs used by 'drizzle'.

    Parameters:
        img      - ObsGeometry instance for input image
        ref_wcs  - Undistorted WCSObject instance for output frame
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
    for i in xrange(_npos):
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