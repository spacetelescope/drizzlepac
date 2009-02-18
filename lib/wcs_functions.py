import numpy as np
import copy 

from updatewcs import wcsutil
from updatewcs.distortion import utils

from pytools import fileutil
import util
import imageObject
import updatewcs
from updatewcs import pywcs

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
    def forward(self,pixx,pixy):
        skyx,skyy = self.input.all_pix2sky(pixx,pixy,1)
        return self.output.wcs_sky2pix(skyx,skyy,1)
    def backward(self,pixx,pixy):
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
    
    hstwcs_list = []
    for img in imageObjectList:
        hstwcs_list += img.getKeywordList('wcs')

    # Compute default output WCS
    default_wcs = utils.output_wcs(hstwcs_list)
           
    # Merge in user-specified attributes for the output WCS
    # as recorded in the input configObj object.
    user_pars = DEFAULT_WCS_PARS.copy()
    
    # More interpretation of the configObj needs to be done here to translate
    # the input parameter names to those understood by 'mergeWCS' as defined
    # by the DEFAULT_WCS_PARS dictionary.
    if configObj is not None:
        user_pars.update(configObj)

    # Apply user settings to output_wcs
    final_wcs = mergeWCS(default_wcs,user_pars)

    # Turn WCS instances into WCSObject instances
    outwcs = createWCSObject(output,default_wcs,final_wcs,imageObjectList)
    
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
    
    outwcs = copy.copy(default_wcs)
    outwcs.wcs = copy.copy(default_wcs.wcs)
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
