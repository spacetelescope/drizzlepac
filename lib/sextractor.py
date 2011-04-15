# ======================================================================
#
# Cosmograil: cosmograil.tools.sextractor
#
# sextractor module.
#
# Author: Laurent Le Guillou <laurentl@ster.kuleuven.ac.be>
#
# $Id: sextractor.py,v 1.2 2005/07/06 21:40:43 hack Exp $
#
# ======================================================================
#
# "sextractor": wrapper around SExtractor.
#
# ======================================================================
#
# $Log: sextractor.py,v $
# Revision 1.2  2005/07/06 21:40:43  hack
# Tweakshifts version 0.5.0 (WJH):
#   - added support for SExtractor PSET and user-supplied SExtractor config file
#   - added 'nbright' parameter for selecting only 'nbright' objects for matching
#   - redefined 'ascend' to 'fluxunits' of 'counts/cps/mag'
#   - fixed bug in countSExtractorObjects()reported by Andy
#   - turned off overwriting of output WCS file
#
# Revision 1.15  2005/06/29 13:07:41  hack
# Added Python interface to SExtractor to STSDAS$Python for use with 'tweakshifts'. WJH
# Added 3 more parameters to config
#
# Revision 1.14  2005/02/14 19:27:31  laurentl
# Added write facilities to rdb module.
#
# Revision 1.13  2005/02/14 17:47:02  laurentl
# Added iterator interface
#
# Revision 1.12  2005/02/14 17:16:30  laurentl
# clean now removes the NNW config file too.
#
# Revision 1.2  2005/02/14 17:13:49  laurentl
# *** empty log message ***
#
# Revision 1.1  2005/02/14 11:34:10  laurentl
# quality monitor now uses SExtractor wrapper.
#
# Revision 1.10  2005/02/11 14:40:35  laurentl
# minor changes
#
# Revision 1.9  2005/02/11 14:32:44  laurentl
# Fixed bugs in setup()
#
# Revision 1.8  2005/02/11 13:50:08  laurentl
# Fixed bugs in setup()
#
# Revision 1.7  2005/02/10 20:15:14  laurentl
# Improved SExtractor wrapper.
#
# Revision 1.6  2005/02/10 17:46:35  laurentl
# Greatly improved the SExtractor wrapper.
#
# Revision 1.5  2005/02/09 23:32:50  laurentl
# Implemented SExtractor wrapper
#
# Revision 1.4  2005/02/04 05:00:09  laurentl
# *** empty log message ***
#
# Revision 1.3  2005/01/06 13:37:11  laurentl
# *** empty log message ***
# 
#
# ======================================================================


"""
A wrapper for SExtractor

A wrapper for SExtractor, the Source Extractor.
by Laurent Le Guillou
version: 1.15 - last modified: 2005-07-06

This wrapper allows you to configure SExtractor, run it and get
back its outputs without the need of editing SExtractor
configuration files. by default, configuration files are created
on-the-fly, and SExtractor is run silently via python.

Tested on SExtractor versions 2.2.1 and 2.3.2.


Example of use:

-----------------------------------------------------------------

    import sextractor

    # Create a SExtractor instance
    sex = sextractor.SExtractor()

    # Modify the SExtractor configuration
    sex.config['GAIN'] = 0.938
    sex.config['PIXEL_SCALE'] = .19
    sex.config['VERBOSE_TYPE'] = "FULL"
    sex.config['CHECKIMAGE_TYPE'] = "BACKGROUND"

    # Add a parameter to the parameter list
    sex.config['PARAMETERS_LIST'].append('FLUX_BEST')

    # Lauch SExtractor on a FITS file
    sex.run("nf260002.fits")

    # Read the resulting catalog [first method, whole catalog at once]
    catalog = sex.catalog()
    for star in catalog:
        print star['FLUX_BEST'], star['FLAGS']
        if (star['FLAGS'] & sextractor.BLENDED):
            print "This star is BLENDED"

    # Read the resulting catalog [second method, whole catalog at once]
    catalog_name = sex.config['CATALOG_NAME']
    catalog_f = sextractor.open(catalog_name)
    catalog = catalog_f.readlines()
    for star in catalog:
        print star['FLUX_BEST'], star['FLAGS']
        if (star['FLAGS'] & sextractor.BLENDED):
            print "This star is BLENDED"
    catalog_f.close()

    # Read the resulting catalog [third method, star by star]
    catalog_name = sex.config['CATALOG_NAME']
    catalog_f = sextractor.open(catalog_name)
    star = catalog_f.readline()
    while star:
        print star['FLUX_BEST'], star['FLAGS']
        if (star['FLAGS'] & sextractor.BLENDED):
            print "This star is BLENDED"
        star = catalog_f.readline()
    catalog_f.close()

    # Removing the configuration files, the catalog and
    # the check image
    sex.clean(config=True, catalog=True, check=True)

-----------------------------------------------------------------


"""

# ======================================================================
from __future__ import division # confidence high

import __builtin__

import sys
import os
import popen2
import exceptions
import re
import copy

from sexcatalog import *


# ======================================================================

__version__ = "1.15.0 (2005-07-06)"

# ======================================================================

class SExtractorException(Exception):
    pass

# ======================================================================

nnw_config = \
"""NNW
# Neural Network Weights for the SExtractor star/galaxy classifier (V1.3)
# inputs:	9 for profile parameters + 1 for seeing.
# outputs:	``Stellarity index'' (0.0 to 1.0)
# Seeing FWHM range: from 0.025 to 5.5'' (images must have 1.5 < FWHM < 5 pixels)
# Optimized for Moffat profiles with 2<= beta <= 4.

 3 10 10  1

-1.56604e+00 -2.48265e+00 -1.44564e+00 -1.24675e+00 -9.44913e-01 -5.22453e-01  4.61342e-02  8.31957e-01  2.15505e+00  2.64769e-01
 3.03477e+00  2.69561e+00  3.16188e+00  3.34497e+00  3.51885e+00  3.65570e+00  3.74856e+00  3.84541e+00  4.22811e+00  3.27734e+00

-3.22480e-01 -2.12804e+00  6.50750e-01 -1.11242e+00 -1.40683e+00 -1.55944e+00 -1.84558e+00 -1.18946e-01  5.52395e-01 -4.36564e-01 -5.30052e+00
 4.62594e-01 -3.29127e+00  1.10950e+00 -6.01857e-01  1.29492e-01  1.42290e+00  2.90741e+00  2.44058e+00 -9.19118e-01  8.42851e-01 -4.69824e+00
-2.57424e+00  8.96469e-01  8.34775e-01  2.18845e+00  2.46526e+00  8.60878e-02 -6.88080e-01 -1.33623e-02  9.30403e-02  1.64942e+00 -1.01231e+00
 4.81041e+00  1.53747e+00 -1.12216e+00 -3.16008e+00 -1.67404e+00 -1.75767e+00 -1.29310e+00  5.59549e-01  8.08468e-01 -1.01592e-02 -7.54052e+00
 1.01933e+01 -2.09484e+01 -1.07426e+00  9.87912e-01  6.05210e-01 -6.04535e-02 -5.87826e-01 -7.94117e-01 -4.89190e-01 -8.12710e-02 -2.07067e+01
-5.31793e+00  7.94240e+00 -4.64165e+00 -4.37436e+00 -1.55417e+00  7.54368e-01  1.09608e+00  1.45967e+00  1.62946e+00 -1.01301e+00  1.13514e-01
 2.20336e-01  1.70056e+00 -5.20105e-01 -4.28330e-01  1.57258e-03 -3.36502e-01 -8.18568e-02 -7.16163e+00  8.23195e+00 -1.71561e-02 -1.13749e+01
 3.75075e+00  7.25399e+00 -1.75325e+00 -2.68814e+00 -3.71128e+00 -4.62933e+00 -2.13747e+00 -1.89186e-01  1.29122e+00 -7.49380e-01  6.71712e-01
-8.41923e-01  4.64997e+00  5.65808e-01 -3.08277e-01 -1.01687e+00  1.73127e-01 -8.92130e-01  1.89044e+00 -2.75543e-01 -7.72828e-01  5.36745e-01
-3.65598e+00  7.56997e+00 -3.76373e+00 -1.74542e+00 -1.37540e-01 -5.55400e-01 -1.59195e-01  1.27910e-01  1.91906e+00  1.42119e+00 -4.35502e+00

-1.70059e+00 -3.65695e+00  1.22367e+00 -5.74367e-01 -3.29571e+00  2.46316e+00  5.22353e+00  2.42038e+00  1.22919e+00 -9.22250e-01 -2.32028e+00


 0.00000e+00 
 1.00000e+00 
"""

# ======================================================================

class SExtractor:
    """
    A wrapper class to transparently use SExtractor.

    """

    _SE_config = { 
        "CATALOG_NAME":
        {"comment": "name of the output catalog",
         "value": "py-sextractor.cat"},
        
        "CATALOG_TYPE":
        {"comment":
         '"NONE","ASCII_HEAD","ASCII","FITS_1.0" or "FITS_LDAC"',
         "value": "ASCII_HEAD"},
        
        "PARAMETERS_NAME":
        {"comment": "name of the file containing catalog contents",
         "value": "py-sextractor.param"},
        
        "DETECT_TYPE":
        {"comment": '"CCD" or "PHOTO"',
         "value": "CCD"},
        
        "FLAG_IMAGE":
        {"comment": "filename for an input FLAG-image",
         "value": "flag.fits"},
        
        "DETECT_MINAREA":
        {"comment": "minimum number of pixels above threshold",
         "value": 5},
        
        "DETECT_THRESH":
        {"comment": "<sigmas> or <threshold>,<ZP> in mag.arcsec-2",
         "value": 1.5},
        
        "ANALYSIS_THRESH":
        {"comment": "<sigmas> or <threshold>,<ZP> in mag.arcsec-2",
         "value": 1.5},
        
        "FILTER":
        {"comment": 'apply filter for detection ("Y" or "N")',
         "value": 'Y'},
        
        "FILTER_NAME":
        {"comment": "name of the file containing the filter",
         "value": "py-sextractor.conv"},
        
        "DEBLEND_NTHRESH":
        {"comment": "Number of deblending sub-thresholds",
         "value": 32},
        
        "DEBLEND_MINCONT":
        {"comment": "Minimum contrast parameter for deblending",
         "value": 0.005},
        
        "CLEAN":
        {"comment": "Clean spurious detections (Y or N)",
         "value": 'Y'},
        
        "CLEAN_PARAM":
        {"comment": "Cleaning efficiency",
         "value": 1.0},
        
        "MASK_TYPE":
        {"comment": 'type of detection MASKing: can be one of "NONE", "BLANK" or "CORRECT"',
         "value": "CORRECT"},
        
        "PHOT_APERTURES":
        {"comment": "MAG_APER aperture diameter(s) in pixels",
         "value": 5},
        
        "PHOT_AUTOPARAMS":
        {"comment": 'MAG_AUTO parameters: <Kron_fact>,<min_radius>',
         "value": [2.5, 3.5]},
        
        "SATUR_LEVEL":
        {"comment": "level (in ADUs) at which arises saturation",
         "value": 50000.0},
        
        "MAG_ZEROPOINT":
        {"comment": "magnitude zero-point",
         "value": 0.0},
        
        "MAG_GAMMA":
        {"comment": "gamma of emulsion (for photographic scans)",
         "value": 4.0},
        
        "GAIN":
        {"comment": "detector gain in e-/ADU",
         "value": 0.0},
        
        "PIXEL_SCALE":
        {"comment": "size of pixel in arcsec (0=use FITS WCS info)",
         "value": 1.0},
        
        "SEEING_FWHM":
        {"comment": "stellar FWHM in arcsec",
         "value": 1.2},
        
        "STARNNW_NAME":
        {"comment": "Neural-Network_Weight table filename",
         "value": "py-sextractor.nnw"},
        
        "BACK_SIZE":
        {"comment": "Background mesh: <size> or <width>,<height>",
         "value": 64},

        "BACK_TYPE":
        {"comment": "Type of background to subtract: MANUAL or AUTO generated",
         "value": 'MANUAL'},

        "BACK_VALUE":
        {"comment": "User-supplied constant value to be subtracted as sky",
         "value": "0.0,0.0"},
        
        "BACK_FILTERSIZE":
        {"comment": "Background filter: <size> or <width>,<height>",
         "value": 3},
        
        "BACKPHOTO_TYPE":
        {"comment": 'can be "GLOBAL" or "LOCAL"',
         "value": "GLOBAL"},

        "BACKPHOTO_THICK":
        {"comment": "Thickness in pixels of the background local annulus",
         "value": 24},
        
        "CHECKIMAGE_TYPE":
        {"comment": 'can be one of "NONE", "BACKGROUND", "MINIBACKGROUND", "-BACKGROUND", "OBJECTS", "-OBJECTS", "SEGMENTATION", "APERTURES", or "FILTERED"',
         "value": "NONE"},
        
        "CHECKIMAGE_NAME":
        {"comment": "Filename for the check-image",
         "value": "check.fits"},
        
        "MEMORY_OBJSTACK":
        {"comment": "number of objects in stack",
         "value": 3000},
        
        "MEMORY_PIXSTACK":
        {"comment": "number of pixels in stack",
         "value": 300000},
        
        "MEMORY_BUFSIZE":
        {"comment": "number of lines in buffer",
         "value": 1024},
        
        "VERBOSE_TYPE":
        {"comment": 'can be "QUIET", "NORMAL" or "FULL"',
         "value": "QUIET"},

        # -- Extra-keys (will not be saved in the main configuration file

        "PARAMETERS_LIST":
        {"comment": '[Extra key] catalog contents (to put in PARAMETERS_NAME)',
         "value": ["NUMBER", "FLUX_BEST", "FLUXERR_BEST", 
                   "X_IMAGE", "Y_IMAGE", "FLAGS", "FWHM_IMAGE"]},

        "CONFIG_FILE":
        {"comment": '[Extra key] name of the main configuration file',
         "value": "py-sextractor.sex"},

        "FILTER_MASK":
        {"comment": 'Array to put in the FILTER_MASK file',
         "value": [[1, 2, 1],
                   [2, 4, 2],
                   [1, 2, 1]]}
        }

    
    # -- Special config. keys that should not go into the config. file.

    _SE_config_special_keys = ["PARAMETERS_LIST", "CONFIG_FILE", "FILTER_MASK"]


    # -- Dictionary of all possible parameters (from sexcatalog.py module)

    _SE_parameters = SExtractorfile._SE_keys



    def __init__(self):
        """
        SExtractor class constructor.
        """

        self.config = (
            dict([(k, copy.deepcopy(SExtractor._SE_config[k]["value"]))\
                  for k in SExtractor._SE_config.keys()]))

        # print self.config

        self.program = None
        self.version = None


    def setup(self, path=None):
        """
        Look for SExtractor program ('sextractor', or 'sex').
        If a full path is provided, only this path is checked.
        Raise a SExtractorException if it failed.
        Return program and version if it succeed.
        """

        # -- Finding sextractor program and its version
        # first look for 'sextractor', then 'sex'

        candidates = ['sextractor', 'sex']

        if (path):
            candidates = [path]
        
        selected=None
        for candidate in candidates:
            try:
                (_out_err, _in) = popen2.popen4(candidate)
                versionline = _out_err.read()
                if (versionline.find("SExtractor") != -1):
                    selected=candidate
                    break
            except IOError:
                continue
                
        if not(selected):
            raise SExtractorException, \
                  """
                  Cannot find SExtractor program. Check your PATH,
                  or provide the SExtractor program path in the constructor.
                  """

        _program = selected

        # print versionline
        _version_match = re.search("[Vv]ersion ([0-9\.])+", versionline)
        if not _version_match:
            raise SExtractorException, \
                  "Cannot determine SExtractor version."

        _version = _version_match.group()[8:]
        if not _version:
            raise SExtractorException, \
                  "Cannot determine SExtractor version."

        # print "Use " + self.program + " [" + self.version + "]"

        return _program, _version



    def update_config(self):
        """
        Update the configuration files according to the current
        in-memory SExtractor configuration.
        """

        # -- Write filter configuration file

        # First check the filter itself

        filter = self.config['FILTER_MASK']
        rows = len(filter)
        cols = len(filter[0])   # May raise ValueError, OK

        filter_f = __builtin__.open(self.config['FILTER_NAME'], 'w')
        filter_f.write("CONV NORM\n")
        filter_f.write("# %dx%d Generated from sextractor.py module.\n" %
                       (rows, cols))
        for row in filter:
            filter_f.write(" ".join(map(repr, row)))
            filter_f.write("\n")
            
        filter_f.close()

        # -- Write parameter list file

        parameters_f = __builtin__.open(self.config['PARAMETERS_NAME'], 'w')
        for parameter in self.config['PARAMETERS_LIST']:
            print >>parameters_f, parameter

        parameters_f.close()

        # -- Write NNW configuration file

        nnw_f = __builtin__.open(self.config['STARNNW_NAME'], 'w')
        nnw_f.write(nnw_config)
        nnw_f.close()


        # -- Write main configuration file

        main_f = __builtin__.open(self.config['CONFIG_FILE'], 'w')

        for key in self.config.keys():
            if (key in SExtractor._SE_config_special_keys):
                continue

            if (key == "PHOT_AUTOPARAMS"): # tuple instead of a single value
                value = " ".join(map(str, self.config[key]))
            else:
                value = str(self.config[key])
            
            
            print >>main_f, ("%-16s       %-16s # %s" %
                             (key, value, SExtractor._SE_config[key]['comment']))

        main_f.close()


    def run(self, file, updateconfig=True, clean=False, path=None):
        """
        Run SExtractor.

        If updateconfig is True (default), the configuration
        files will be updated before running SExtractor.

        If clean is True (default: False), configuration files 
        (if any) will be deleted after SExtractor terminates.

        """

        if updateconfig:
            self.update_config()

        # Try to find SExtractor program
        # This will raise an exception if it failed

        self.program, self.version = self.setup(path)

        commandline = (
            self.program + " -c " + self.config['CONFIG_FILE'] + " " + file)
        # print commandline

        rcode = os.system(commandline)

        if (rcode):
            raise SExtractorException, \
                  "SExtractor command [%s] failed." % commandline
            
        if clean:
            self.clean()



    def catalog(self):
        """
        Read the output catalog produced by the last SExtractor run.
        Output is a list of dictionaries, with a dictionary for
        each star: {'param1': value, 'param2': value, ...}.
        """

        output_f = SExtractorfile(self.config['CATALOG_NAME'], 'r')
        c = output_f.read()
        output_f.close()

        return c


    def clean(self, config=True, catalog=False, check=False):
        """
        Remove the generated SExtractor files (if any).
        If config is True, remove generated configuration files.
        If catalog is True, remove the output catalog.
        If check is True, remove output check image.
        """

        try:
            if (config):
                os.unlink(self.config['FILTER_NAME'])
                os.unlink(self.config['PARAMETERS_NAME'])
                os.unlink(self.config['STARNNW_NAME'])
                os.unlink(self.config['CONFIG_FILE'])
            if (catalog):
                os.unlink(self.config['CATALOG_NAME'])
            if (check):
                os.unlink(self.config['CHECKIMAGE_NAME'])
                
        except OSError:
            pass



# ======================================================================
