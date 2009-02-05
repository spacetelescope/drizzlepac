import numpy as np

from updatewcs import wcsutil
from updatewcs.distortion import utils

from pytools import fileutil
import util

DEFAULT_WCS_PARS = {'ra':None,'dec':None,'psize':None,'orient':None,
                     'outnx':None,'outny':None,'crpix1':None,'crpix2':None,
                     'crval1':None,'crval2':None}
                    
def get_hstwcs(filename,hdulist,hdr,primary_hdr):
    ''' Return the HSTWCS object for a given chip.
    
    '''
    hdrwcs = wcsutil.HSTWCS(primary_hdr,hdr,hdulist)
    hdrwcs.filename = filename
    hdrwcs.expname = hdr['expname']
    hdrwcs.extver = hdr['extver']
    
    return hdrwcs

#
# Possibly need to generate a stand-alone interface for this function.
#
# Primary interface for creating the output WCS from a list of HSTWCS objects
def make_outputwcs(imageObjectList,configObj=None):
    hstwcs_list = []
    for img in imageObjectList:
        hstwcs_list += img.getKeywordList('wcs')
        
    output_wcs = utils.output_wcs(hstwcs_list)
    
    # Merge in user-specified attributes for the output WCS
    # as recorded in the input configObj object.
    user_pars = DEFAULT_WCS_PARS.copy()
    
    # More interpretation of the configObj needs to be done here to translate
    # the input parameter names to those understood by 'mergeWCS' as defined
    # by the DEFAULT_WCS_PARS dictionary.
    if configObj is not None:
        user_pars.update(configObj)

    # Apply user settings to output_wcs
    mergeWCS(output_wcs,user_pars)
    
    return output_wcs

def mergeWCS(outwcs,user_pars):
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
