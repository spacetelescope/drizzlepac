"""
Provides function for manipulating WCS in images.

:Authors: Warren Hack

:License: :doc:`LICENSE`

"""
import sys, types, os, copy
from . import util
import numpy as np
from astropy.io import fits

from stsci.tools import fileutil,teal
from . import outputimage, wcs_functions, processInput,util
import stwcs
from stwcs import distortion, wcsutil
from stwcs.wcsutil import headerlet

__taskname__ = 'buildwcs'

# This is specifically NOT intended to match the package-wide version information.
__version__ = '0.1.0'
__version_date__ = '22-June-2011'

# These default parameter values have the same keys as the parameters from
# the configObj interface
default_user_wcs = {'pscale':None,'orientat':None,'crval1':None,'crval2':None,
                    'crpix1':None,'crpix2':None,'naxis1':None,'naxis2':None}

user_hstwcs_pars = {'outscale':'pscale','orient':'orientat',
                    'raref':'crval1','decref':'crval2',
                    'xrefpix':'crpix1','yrefpix':'crpix2',
                    'outnx':'naxis1','outny':'naxis2'}
model_attrs = ['cpdis1','cpdis2','det2im','det2im1','det2im2',
                    'ocx10','ocx11','ocy10','ocy11','sip']

#
#### User level interface run from TEAL
#

def buildwcs(outwcs, configObj=None,editpars=False,**input_dict):
    if input_dict is None:
        input_dict = {}
    input_dict['outwcs'] = outwcs

    # If called from interactive user-interface, configObj will not be
    # defined yet, so get defaults using EPAR/TEAL.
    #
    # Also insure that the input_dict (user-specified values) are folded in
    # with a fully populated configObj instance.
    configObj = util.getDefaultConfigObj(__taskname__,configObj,input_dict,loadOnly=(not editpars))
    if configObj is None:
        return

    if not editpars:
        run(configObj,wcsmap=wcsmap)

def run(configObj,wcsmap=None):
    """ Interpret parameters from TEAL/configObj interface as set interactively
        by the user and build the new WCS instance
    """

    distortion_pars = configObj['Distortion Model']

    outwcs = build(configObj['outwcs'], configObj['wcsname'],
            configObj['refimage'], undistort = configObj['undistort'],
            usecoeffs=distortion_pars['applycoeffs'], coeffsfile=distortion_pars['coeffsfile'],
            **configObj['User WCS Parameters'])


#### Low-level interface using Python objects only
#

def build(outname, wcsname, refimage, undistort=False,
          applycoeffs=False, coeffsfile=None, **wcspars):
    """ Core functionality to create a WCS instance from a reference image WCS,
        user supplied parameters or user adjusted reference WCS.
        The distortion information can either be read in as part of the reference
        image WCS or given in 'coeffsfile'.

        Parameters
        ----------
        outname   : string
            filename of output WCS
        wcsname   : string
            WCSNAME ID for generated WCS
        refimage  : string
            filename of image with source WCS used as basis for output WCS
        undistort : bool
            Create an undistorted WCS?
        applycoeffs : bool
            Apply coefficients from refimage to generate undistorted WCS?
        coeffsfile  : string
            If specified, read distortion coeffs from separate file


    """

    # Insure that the User WCS parameters have values for all the parameters,
    # even if that value is 'None'
    user_wcs_pars = convert_user_pars(wcspars)
    userwcs = wcspars['userwcs']

    """
    Use cases to document the logic required to interpret the parameters

    WCS generation based on refimage/userwcs parameters
    -------------------------------------------------------------
    refimage == None, userwcs == False:
        *NO WCS specified*
        => print a WARNING message and return without doing anything
    refimage == None, userwcs == True:
        => Create WCS without a distortion model entirely from user parameters*
    refimage != None, userwcs == False:
        => No user WCS parameters specified
        => Simply use refimage WCS as specified
    refimage != None, userwcs == True:
        => Update refimage WCS with user specified values*

    Apply distortion and generate final headerlet using processed WCS
    -----------------------------------------------------------------
    refimage == None, userwcs == True:
        *Output WCS generated entirely from user supplied parameters*
        Case 1: applycoeffs == False, undistort == True/False (ignored)
            => no distortion model to interpret
            => generate undistorted headerlet with no distortion model
        Case 2: applycoeffs == True/False, undistort == True
            => ignore any user specified distortion model
            => generate undistorted headerlet with no distortion model
        Case 3: applycoeffs == True, undistort == False
            => WCS from scratch combined with distortion model from another image
            => generate headerlet with distortion model

    refimage != None, userwcs == True/False:
        *Output WCS generated from reference image possibly modified by user parameters*
        Case 4: applycoeffs == False, undistort == True
            => If refimage has distortion, remove it
            => generate undistorted headerlet with no distortion model
        Case 5: applycoeffs == False, undistort == False
            => Leave refimage distortion model (if any) unmodified
            => generate a headerlet using same distortion model (if any) as refimage
        Case 6: applycoeffs == True, undistort == False
            => Update refimage with distortion model with user-specified model
            => generate a headerlet with a distortion model
        Case 7: applycoeffs == True, undistort == True
            => ignore user specified distortion model and undistort WCS
            => generate a headerlet without a distortion model
    """
    ### Build WCS from refimage and/or user pars
    if util.is_blank(refimage) and not userwcs:
        print('WARNING: No WCS specified... No WCS created!')
        return
    customwcs = None
    if util.is_blank(refimage) and userwcs:
        # create HSTWCS object from user parameters
        complete_wcs = True
        for key in user_wcs_pars:
            if util.is_blank(user_wcs_pars[key]):
                complete_wcs = False
                break
        if complete_wcs:
            customwcs = wcs_functions.build_hstwcs(user_wcs_pars['crval1'],user_wcs_pars['crval2'],
                user_wcs_pars['crpix1'],user_wcs_pars['crpix2'],
                user_wcs_pars['naxis1'],user_wcs_pars['naxis2'],
                user_wcs_pars['pscale'],user_wcs_pars['orientat'])
        else:
            print('WARNING: Not enough WCS information provided by user!')
            raise ValueError

    if not util.is_blank(refimage):
        refwcs = stwcs.wcsutil.HSTWCS(refimage)
    else:
        refwcs = customwcs

    ### Apply distortion model (if any) to update WCS
    if applycoeffs and not util.is_blank(coeffsfile):
        if not util.is_blank(refimage):
            replace_model(refwcs, coeffsfile)
        else:
            if not undistort:
                add_model(refwcs,coeffsfile)
                # Only working with custom WCS from user, no distortion
                # so apply model to WCS, including modifying the CD matrix
                apply_model(refwcs)

    ### Create undistorted WCS, if requested
    if undistort:
        outwcs = undistortWCS(refwcs)
    else:
        outwcs = refwcs

    if userwcs:
        # replace (some/all?) WCS values from refimage with user WCS values
        # by running 'updatewcs' functions on input WCS
        outwcs = mergewcs(outwcs,customwcs,user_wcs_pars)


    ### Create the final headerlet and write it out, if specified
    if not util.is_blank(refimage):
        template = refimage
    elif not util.is_blank(coeffsfile):
        template = coeffsfile
    else:
        template = None


    # create default WCSNAME if None was given
    wcsname = create_WCSname(wcsname)
    print('Creating final headerlet with name ',wcsname,' using template ',template)
    outhdr = generate_headerlet(outwcs,template,wcsname,outname=outname)

    # synchronize this new WCS with the rest of the chips in the image
    for ext in outhdr:
        if 'extname' in ext.header and ext.header['extname'] == 'SIPWCS':
            ext_wcs = wcsutil.HSTWCS(ext)
            stwcs.updatewcs.makewcs.MakeWCS.updateWCS(ext_wcs,outwcs)

    return outwcs

def create_WCSname(wcsname):
    """ Verify that a valid WCSNAME has been provided, and if not, create a
        default WCSNAME based on current date.
    """
    if util.is_blank(wcsname):
        ptime = fileutil.getDate()
        wcsname = "User_"+ptime

    return wcsname

def convert_user_pars(wcspars):
    """ Convert the parameters provided by the configObj into the corresponding
        parameters from an HSTWCS object
    """
    default_pars = default_user_wcs.copy()
    for kw in user_hstwcs_pars:
        default_pars[user_hstwcs_pars[kw]] = wcspars[kw]
    return default_pars

def mergewcs(outwcs, customwcs, wcspars):
    """ Merge the WCS keywords from user specified values into a full HSTWCS object
        This function will essentially follow the same algorithm as used by
        updatehdr only it will use direct calls to updatewcs.Makewcs methods
        instead of using 'updatewcs' as a whole
    """
    # start by working on a copy of the refwcs
    if outwcs.sip is not None:
        wcslin = stwcs.distortion.utils.undistortWCS(outwcs)
        outwcs.wcs.cd = wcslin.wcs.cd
        outwcs.wcs.set()
        outwcs.setOrient()
        outwcs.setPscale()
    else:
        wcslin = outwcs
    if customwcs is None:
        # update valid pars from wcspars
        if wcspars['crval1'] is not None:
            outwcs.wcs.crval = np.array([wcspars['crval1'],wcspars['crval2']])
        if wcspars['crpix1'] is not None:
            outwcs.wcs.crpix = np.array([wcspars['crpix1'],wcspars['crpix2']])
        if wcspars['naxis1'] is not None:
            outwcs.pixel_shape = (wcspars['naxis1'], wcspars['naxis2'])
            outwcs.wcs.crpix = np.array(outwcs.pixel_shape) / 2.0

        pscale = wcspars['pscale']
        orient = wcspars['orientat']
        if pscale is not None or orient is not None:
            if pscale is None: pscale = wcslin.pscale
            if orient is None: orient = wcslin.orientat
            pix_ratio = pscale/wcslin.pscale
            delta_rot = wcslin.orientat - orient
            delta_rot_mat = fileutil.buildRotMatrix(delta_rot)
            outwcs.wcs.cd = np.dot(outwcs.wcs.cd,delta_rot_mat)*pix_ratio
            # apply model to new linear CD matrix
            apply_model(outwcs)
    else:
        # A new fully described WCS was provided in customwcs
        outwcs.wcs.cd = customwcs.wcs.cd
        outwcs.wcs.crval = customwcs.wcs.crval
        outwcs.wcs.crpix = customwcs.wcs.crpix
        outwcs.pixel_shape = customwcs.pixel_shape
    return outwcs

def add_model(refwcs, newcoeffs):
    """ Add (new?) distortion model to existing HSTWCS object
    """
    # Update refwcs with distortion model
    for kw in model_attrs:
        if newcoeffs.__dict__[key] is not None:
            refwcs.__dict__[key] = newcoeffs.__dict__[key]

def apply_model(refwcs):
    """ Apply distortion model to WCS, including modifying
        CD with linear distortion terms

    """
    # apply distortion model to CD matrix
    if 'ocx10' in refwcs.__dict__ and refwcs.ocx10 is not None:
        linmat = np.array([[refwcs.ocx11,refwcs.ocx10],[refwcs.ocy11,refwcs.ocy10]])/refwcs.idcscale
        refwcs.wcs.cd = np.dot(refwcs.wcs.cd,linmat)

        refwcs.wcs.set()
        refwcs.setOrient()
        refwcs.setPscale()

def replace_model(refwcs, newcoeffs):
    """ Replace the distortion model in a current WCS with a new model
        Start by creating linear WCS, then run
    """
    print('WARNING:')
    print('    Replacing existing distortion model with one')
    print('    not necessarily matched to the observation!')
    # create linear version of WCS to be updated by new model
    wcslin = stwcs.distortion.utils.undistortWCS(refwcs)
    outwcs = refwcs.deepcopy()
    outwcs.wcs.cd = wcslin.wcs.cd
    outwcs.wcs.set()
    outwcs.setOrient()
    outwcs.setPscale()
    # add new model to updated WCS object
    add_model(outwcs,newcoeffs)
    # Update CD matrix with new model
    apply_model(outwcs)

    # replace original input WCS with newly updated WCS
    refwcs = outwcs.deepcopy()

def undistortWCS(refwcs):
    """ Generate an undistorted HSTWCS from an HSTWCS object with a distortion model
    """
    wcslin = stwcs.distortion.utils.output_wcs([refwcs])

    outwcs = stwcs.wcsutil.HSTWCS()
    outwcs.wcs = wcslin.wcs
    outwcs.wcs.set()
    outwcs.setPscale()
    outwcs.setOrient()
    outwcs.sip = None

    # Update instrument specific keywords
    outwcs.inst_kw = refwcs.inst_kw
    for kw in refwcs.inst_kw:
        outwcs.__dict__[kw] = refwcs.__dict__[kw]
    outwcs.pixel_shape = wcslin.pixel_shape

    return outwcs

def generate_headerlet(outwcs,template,wcsname,outname=None):
    """ Create a headerlet based on the updated HSTWCS object

        This function uses 'template' as the basis for the headerlet.
        This file can either be the original wcspars['refimage'] or
        wcspars['coeffsfile'], in this order of preference.

        If 'template' is None, then a simple Headerlet will be
        generated with a single SIPWCS extension and no distortion
    """
    # Create header object from HSTWCS object
    siphdr = True
    if outwcs.sip is None:
        siphdr = False
    outwcs_hdr = outwcs.wcs2header(sip2hdr=siphdr)
    outwcs_hdr['NPIX1'] = outwcs.pixel_shape[0]
    outwcs_hdr['NPIX2'] = outwcs.pixel_shape[1]

    # create headerlet object in memory; either from a file or from scratch
    if template is not None and siphdr:
        print('Creating headerlet from template...')
        fname,extn = fileutil.parseFilename(template)
        extnum = fileutil.parseExtn(extn)
        extname = ('sipwcs',extnum[1])
        hdrlet = headerlet.createHeaderlet(fname,wcsname)
        # update hdrlet with header values from outwcs
        for kw in outwcs_hdr.items():
            hdrlet[extname].header[kw[0]] = kw[1]
        hdrlet[extname].header['WCSNAME'] = wcsname
    else:
        print('Creating headerlet from scratch...')
        hdrlet = fits.HDUList()
        hdrlet.append(fits.PrimaryHDU())
        siphdr = fits.ImageHDU(header=outwcs_hdr)
        siphdr.header['EXTNAME'] = 'SIPWCS'
        siphdr.header['WCSNAME'] = wcsname
        hdrlet.append(siphdr)

    # Write out header to a file as the final product
    if outname is not None:
        if outname.find('_hdr.fits') < 0:
            outname += '_hdr.fits'
        if os.path.exists(outname):
            print('Overwrite existing file "%s"'%outname)
            os.remove(outname)
        hdrlet.writeto(outname)
        print('Wrote out headerlet :',outname)


def help(file=None):
    """
    Print out syntax help for running a drizzlepac task.

    Parameters
    ----------
    file : str (Default = None)
        If given, write out help to the filename specified by this parameter
        Any previously existing file with this name will be deleted before
        writing out the help.

    """
    helpstr = getHelpAsString(docstring=True, show_ver = True)
    if file is None:
        print(helpstr)
    else:
        if os.path.exists(file): os.remove(file)
        f = open(file, mode = 'w')
        f.write(helpstr)
        f.close()


def getHelpAsString(docstring = False, show_ver = True):
    """
    return useful help from a file in the script directory called
    __taskname__.help

    """
    install_dir = os.path.dirname(__file__)
    taskname = util.base_taskname(__taskname__, '')
    htmlfile = os.path.join(install_dir, 'htmlhelp', taskname + '.html')
    helpfile = os.path.join(install_dir, taskname + '.help')

    if docstring or (not docstring and not os.path.exists(htmlfile)):
        if show_ver:
            helpString = os.linesep + \
                ' '.join([__taskname__, 'Version', __version__,
                ' updated on ', __version_date__]) + 2*os.linesep
        else:
            helpString = ''
        if os.path.exists(helpfile):
            helpString += teal.getHelpFileAsString(taskname, __file__)
        else:
            if __doc__ is not None:
                helpString += __doc__ + os.linesep
    else:
        helpString = 'file://' + htmlfile

    return helpString
