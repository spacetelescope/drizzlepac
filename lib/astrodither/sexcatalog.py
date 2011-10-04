# ======================================================================
#
# Cosmograil: cosmograil.tools.sexcatalog
#
# sexcatalog module.
#
# Author: Laurent Le Guillou <laurentl@ster.kuleuven.ac.be>
#
# $Id: sexcatalog.py,v 1.1 2005/06/29 13:07:41 hack Exp $
#
# ======================================================================
#
# "sexcatalog": python module to read and parse SExtractor catalogs
# A simple interface to read SExtractor text catalogs
#
# ======================================================================
#
# $Log: sexcatalog.py,v $
# Revision 1.1  2005/06/29 13:07:41  hack
# Added Python interface to SExtractor to STSDAS$Python for use with 'tweakshifts'. WJH
#
# Revision 1.9  2005/02/14 19:27:31  laurentl
# Added write facilities to rdb module.
#
# Revision 1.8  2005/02/14 17:47:02  laurentl
# Added iterator interface
#
# Revision 1.7  2005/02/14 17:16:30  laurentl
# clean now removes the NNW config file too.
#
# Revision 1.2  2005/02/14 17:13:49  laurentl
# *** empty log message ***
#
# Revision 1.1  2005/02/14 11:34:10  laurentl
# quality monitor now uses SExtractor wrapper.
#
# Revision 1.5  2005/02/11 14:40:35  laurentl
# minor changes
#
# Revision 1.4  2005/02/10 20:15:14  laurentl
# Improved SExtractor wrapper.
#
# Revision 1.2  2005/02/09 23:32:50  laurentl
# Implemented SExtractor wrapper
#
# Revision 1.1  2005/01/06 12:29:25  laurentl
# Added a SExtractor wrapper module. Renamed sextractor.py sexcatalog.py.
#
# Revision 1.1  2004/12/09 03:06:23  laurentl
# Changed tree structure
#
# Revision 1.5  2004/11/26 18:26:59  laurentl
# Added a module to manage the data tree.
#
# Revision 1.4  2004/11/24 15:11:31  laurentl
# Fixed a lot of bugs in sexcatalog module.
#
# Revision 1.2  2004/11/23 22:38:23  laurentl
# Added sexcatalog module.
#
#
# ======================================================================


"""
A simple interface to manipulate SExtractor ASCII catalogs

A simple interface to manipulate SExtractor ASCII catalogs
through a file-like API (open, read, readline, etc.).
For the moment only reading ('r' mode) is supported.
by Laurent Le Guillou
version: 0.1.5 - last modified: 2005-02-14

Future: implement a 'w' mode to be able to save catalogs
in SExtractor format.

Examples:

-----------------------------------------------------------------

    # Through sexcatalog module
    import sexcatalog

    # Read a SExtractor ASCII catalog

    # First method: read the whole catalog at once
    catalog_f = sexcatalog.open(catalog_name)
    catalog = catalog_f.readlines()
    for star in catalog:
        print star['FLUX_BEST'], star['FLAGS']
        if (star['FLAGS'] & sexcatalog.BLENDED):
            print "This star is BLENDED"
    catalog_f.close()

    # Second method: read the catalog star by star
    catalog_f = sexcatalog.open(catalog_name)
    for star in catalog_f:
        print star['FLUX_BEST'], star['FLAGS']
        if (star['FLAGS'] & sexcatalog.BLENDED):
            print "This star is BLENDED"
    catalog_f.close()

    # -------------

    # Through sextractor module
    import sextractor

    # Read a SExtractor ASCII catalog

    # First method: read the whole catalog at once
    catalog_f = sextractor.open(catalog_name)
    catalog = catalog_f.readlines()
    for star in catalog:
        print star['FLUX_BEST'], star['FLAGS']
        if (star['FLAGS'] & sextractor.BLENDED):
            print "This star is BLENDED"
    catalog_f.close()

    # Second method: read the catalog star by star
    catalog_f = sextractor.open(catalog_name)
    star = catalog_f.readline()
    while star:
        print star['FLUX_BEST'], star['FLAGS']
        if (star['FLAGS'] & sextractor.BLENDED):
            print "This star is BLENDED"
        star = catalog_f.readline()
    catalog_f.close()

-----------------------------------------------------------------


"""

# ======================================================================
from __future__ import division # confidence high

import __builtin__

import sys
import exceptions

# ======================================================================

__version__ = "0.1.5 (2005-02-14)"

# ======================================================================

# -- FLAGS meaning

NEIGHBOURS         =   1
BLENDED            =   2
SATURATED          =   4
TRUNCATED          =   8
CORRUPTED_APER     =  16
CORRUPTED_ISO      =  32
OVERFLOW_DEBLEND   =  64
OVERFLOW_EXTRACT   = 128


class WrongSExtractorfileException(Exception):
    pass

class SExtractorfile:
    """
    A class to manipulate SExtractor ASCII catalogs.

    For the moment only reading ('r' mode) is supported.
      
    """

    _SE_keys = \
             {"NUMBER"         : {"comment": "Running object number",
                                  "infunc": int,
                                  "format": "%10d",
                                  "unit": ""},
              
              "FLAGS"          : {"comment": "Extraction flags",
                                  "infunc": int,
                                  "format": "%3d",
                                  "unit": ""},
              
              "FLUX_ISO"       : {"comment": "Isophotal flux",
                                  "infunc": float,
                                  "format": "%12g",
                                  "unit": "count"},

              "FLUXERR_ISO"    : {"comment": "RMS error for isophotal flux",
                                  "infunc": float,
                                  "format": "%12g",
                                  "unit": "count"},

              "MAG_ISO"        : {"comment": "Isophotal magnitude",
                                  "infunc": float,
                                  "format": "%8.4f",
                                  "unit": "mag"},

              "MAGERR_ISO"     : {"comment":
                                  "RMS error for isophotal magnitude",
                                  "infunc": float,
                                  "format": "%8.4f",
                                  "unit": "mag"},

              "FLUX_ISOCOR"    : {"comment": "Corrected isophotal flux",
                                  "infunc": float,
                                  "format": "%12g",
                                  "unit": "count"},

              "FLUXERR_ISOCOR" : {"comment":
                                  "RMS error for corrected isophotal flux", 
                                  "infunc": float,
                                  "format": "%12g",
                                  "unit": "count"},

              "MAG_ISOCOR"     : {"comment": "Corrected isophotal magnitude",
                                  "infunc": float,
                                  "format": "%8.4f",
                                  "unit": "mag"},

              "MAGERR_ISOCOR"  : {"comment":
                                  "RMS error for corrected isophotal magnitude",
                                  "infunc": float,
                                  "format": "%8.4f",
                                  "unit": "mag"},

              "FLUX_AUTO"      : {"comment":
                                  "Flux within a Kron-like elliptical aperture",
                                  "infunc": float,
                                  "format": "%12g",
                                  "unit": "count"},

              "FLUXERR_AUTO"   : {"comment": "RMS error for AUTO flux",
                                  "infunc": float,
                                  "format": "%12g",
                                  "unit": "count"},
              
              "MAG_AUTO"       : {"comment":
                                  "Kron-like elliptical aperture magnitude", 
                                  "infunc": float,
                                  "format": "%8.4f",
                                  "unit": "mag"},
              
              "MAGERR_AUTO"    : {"comment": "RMS error for AUTO magnitude",
                                  "infunc": float,
                                  "format": "%8.4f",
                                  "unit": "mag"},

              "FLUX_BEST"      : {"comment":
                                  "Best of FLUX_AUTO and FLUX_ISOCOR",
                                  "infunc": float,
                                  "format": "%12g",
                                  "unit": "count"},

              "FLUXERR_BEST"   : {"comment": "RMS error for BEST flux",
                                  "infunc": float,
                                  "format": "%12g",
                                  "unit": "count"},
              
              "MAG_BEST"       : {"comment": "Best of MAG_AUTO and MAG_ISOCOR",
                                  "infunc": float,
                                  "format": "%8.4f",
                                  "unit": "mag"},

              "MAGERR_BEST"    : {"comment": "RMS error for MAG_BEST",
                                  "infunc": float,
                                  "format": "%8.4f",
                                  "unit": "mag"},

              "KRON_RADIUS"    : {"comment":
                                  "Kron apertures in units of A or B",
                                  "infunc": float,
                                  "format": "%5.2f",
                                  "unit": ""},
              "BACKGROUND"     : {"comment": "Background at centroid position",
                                  "infunc": float,
                                  "format": "%12g",
                                  "unit": "count"},

              "THRESHOLD"      : {"comment":
                                  "Detection threshold above background",
                                  "infunc": float,
                                  "format": "%12g",
                                  "unit": "count"},
              
              "MU_THRESHOLD"   : {"comment":
                                  "Detection threshold above background", 
                                  "infunc": float,
                                  "format": "%8.4f",
                                  "unit": "mag * arcsec**(-2)"},

              "FLUX_MAX"       : {"comment": "Peak flux above background",
                                  "infunc": float,
                                  "format": "%12g",
                                  "unit": "count"},
              
              "MU_MAX"         : {"comment":
                                  "Peak surface brightness above background", 
                                  "infunc": float,
                                  "format": "%8.4f",
                                  "unit": "mag * arcsec**(-2)"},  

              "ISOAREA_WORLD"  : {"comment":
                                  "Isophotal area above Analysis threshold", 
                                  "infunc": float,
                                  "format": "%12g",
                                  "unit": "deg**2"},

              "XMIN_IMAGE"     : {"comment":
                                  "Minimum x-coordinate among detected pixels",
                                  "infunc": int,
                                  "format": "%10d",
                                  "unit": "pixel"},
              
              "YMIN_IMAGE"     : {"comment":
                                  "Minimum y-coordinate among detected pixels",
                                  "infunc": int,
                                  "format": "%10d",
                                  "unit": "pixel"},
              
              "XMAX_IMAGE"     : {"comment":
                                  "Maximum x-coordinate among detected pixels",
                                  "infunc": int,
                                  "format": "%10d",
                                  "unit": "pixel"},
              
              "YMAX_IMAGE"     : {"comment":
                                  "Maximum y-coordinate among detected pixels",
                                  "infunc": int,
                                  "format": "%10d",
                                  "unit": "pixel"},
              
              "X_IMAGE"        : {"comment": "Object position along x",
                                  "infunc": float,
                                  "format": "%10.3f",
                                  "unit": "pixel"},
              
              "Y_IMAGE"        : {"comment": "Object position along y",
                                  "infunc": float,
                                  "format": "%10.3f",
                                  "unit": "pixel"},
              
              "X_WORLD"        : {"comment":
                                  "Barycenter position along world x axis",
                                  "infunc": float,
                                  "format": "%15e",
                                  "unit": "deg"},
              
              "Y_WORLD"        : {"comment":
                                  "Barycenter position along world y axis",
                                  "infunc": float,
                                  "format": "%15e",
                                  "unit": "deg"},
              
              "ALPHA_SKY"      : {"comment":
                                  "Right ascension of barycenter (native)",
                                  "infunc": float,
                                  "format": "%11.7f",
                                  "unit": "deg"},
              
              "DELTA_SKY"      : {"comment":
                                  "Declination of barycenter (native)",
                                  "infunc": float,
                                  "format": "%+11.7f",
                                  "unit": "deg"},

              "ALPHA_J2000"    : {"comment":
                                  "Right ascension of barycenter (J2000)",
                                  "infunc": float,
                                  "format": "%11.7f",
                                  "unit": "deg"},
              
              "DELTA_J2000"    : {"comment":
                                  "Declination of barycenter (J2000)",
                                  "infunc": float,
                                  "format": "%+11.7f",
                                  "unit": "deg"},

              "ALPHA_B1950"    : {"comment":
                                  "Right ascension of barycenter (B1950)",
                                  "infunc": float,
                                  "format": "%11.7f",
                                  "unit": "deg"},
              
              "DELTA_B1950"    : {"comment":
                                  "Declination of barycenter (B1950)",
                                  "infunc": float,
                                  "format": "%+11.7f",
                                  "unit": "deg"},

              "X2_IMAGE"       : {"comment": "Variance along x",
                                  "infunc": float,
                                  "format": "%15e",
                                  "unit": "pixel**2"},
              
              "Y2_IMAGE"       : {"comment": "Variance along y",
                                  "infunc": float,
                                  "format": "%15e",
                                  "unit": "pixel**2"},
              
              "XY_IMAGE"       : {"comment": "Covariance between x and y",
                                  "infunc": float,
                                  "format": "%15e",
                                  "unit": "pixel**2"},
              
              "CXX_IMAGE"      : {"comment": "Cxx object ellipse parameter",
                                  "infunc": float,
                                  "format": "%12e",
                                  "unit": "pixel**(-2)"},
              
              "CYY_IMAGE"      : {"comment": "Cyy object ellipse parameter",
                                  "infunc": float,
                                  "format": "%12e",
                                  "unit": "pixel**(-2)"},
              
              "CXY_IMAGE"      : {"comment": "Cxy object ellipse parameter",
                                  "infunc": float,
                                  "format": "%12e",
                                  "unit": "pixel**(-2)"},
              
              "A_IMAGE"        : {"comment": "Profile RMS along major axis",
                                  "infunc": float,
                                  "format": "%9.3f",
                                  "unit": "pixel"},
              
              "B_IMAGE"        : {"comment": "Profile RMS along minor axis",
                                  "infunc": float,
                                  "format": "%9.3f",
                                  "unit": "pixel"},
              
              "THETA_IMAGE"    : {"comment": "Position angle (CCW/x)",
                                  "infunc": float,
                                  "format": "%5.1f",
                                  "unit": "deg"},
              
              "ELONGATION"     : {"comment": "A_IMAGE/B_IMAGE",
                                  "infunc": float,
                                  "format": "%8.3f",
                                  "unit": ""},
              
              "ELLIPTICITY"    : {"comment": "1 - B_IMAGE/A_IMAGE",
                                  "infunc": float,
                                  "format": "%8.3f",
                                  "unit": ""},

              "ERRX2_IMAGE"    : {"comment": "Variance of position along x",
                                  "infunc": float,
                                  "format": "%15e",
                                  "unit": "pixel**2"},
              
              "ERRY2_IMAGE"    : {"comment": "Variance of position along y",
                                  "infunc": float,
                                  "format": "%15e",
                                  "unit": "pixel**2"},
              
              "ERRXY_IMAGE"    : {"comment":
                                  "Covariance of position between x and y",
                                  "infunc": float,
                                  "format": "%15e",
                                  "unit": "pixel**2"},
              
              "ERRCXX_IMAGE"   : {"comment": "Cxx error ellipse parameter",
                                  "infunc": float,
                                  "format": "%12g",
                                  "unit": "pixel**(-2)"},
              
              "ERRCYY_IMAGE"   : {"comment": "Cyy error ellipse parameter",
                                  "infunc": float,
                                  "format": "%12g",
                                  "unit": "pixel**(-2)"},
              
              "ERRCXY_IMAGE"   : {"comment": "Cxy error ellipse parameter",
                                  "infunc": float,
                                  "format": "%12g",
                                  "unit": "pixel**(-2)"},
              
              "ERRA_IMAGE"     : {"comment":
                                  "RMS position error along major axis",
                                  "infunc": float,
                                  "format": "%8.4f",
                                  "unit": "pixel"},

              "ERRB_IMAGE"     : {"comment":
                                  "RMS position error along minor axis",
                                  "infunc": float,
                                  "format": "%8.4f",
                                  "unit": "pixel"},
              
              "ERRTHETA_IMAGE" : {"comment":
                                  "Error ellipse position angle (CCW/x)",
                                  "infunc": float,
                                  "format": "%5.1f",
                                  "unit": "deg"},
              
              "FWHM_IMAGE"     : {"comment": "FWHM assuming a gaussian core",
                                  "infunc": float,
                                  "format": "%8.2f",
                                  "unit": "pixel"},
              
              "X2_WORLD"       : {"comment": "Variance along X-WORLD (alpha)",
                                  "infunc": float,
                                  "format": "%15e",
                                  "unit": "deg**2"},

              "Y2_WORLD"       : {"comment": "Variance along Y-WORLD (delta)",
                                  "infunc": float,
                                  "format": "%15e",
                                  "unit": "deg**2"},
              
              "XY_WORLD"       : {"comment":
                                  "Covariance between X-WORLD and Y-WORLD",
                                  "infunc": float,
                                  "format": "%15e",
                                  "unit": "deg**2"},

              "CXX_WORLD"      : {"comment":
                                  "Cxx object ellipse parameter (WORLD units)",
                                  "infunc": float,
                                  "format": "%12e",
                                  "unit": "deg**(-2)"},

              "CYY_WORLD"      : {"comment":
                                  "Cyy object ellipse parameter (WORLD units)",
                                  "infunc": float,
                                  "format": "%12e",
                                  "unit": "deg**(-2)"},
              
              "CXY_WORLD"      : {"comment":
                                  "Cxy object ellipse parameter (WORLD units)",
                                  "infunc": float,
                                  "format": "%12e",
                                  "unit": "deg**(-2)"},
              
              "A_WORLD"        : {"comment":
                                  "Profile RMS along major axis (world units)",
                                  "infunc": float,
                                  "format": "%12g",
                                  "unit": "deg"},
              
              "B_WORLD"        : {"comment":
                                  "Profile RMS along minor axis (world units)",
                                  "infunc": float,
                                  "format": "%12g",
                                  "unit": "deg"},
              
              "THETA_WORLD"    : {"comment": "Position angle (CCW/world-x)",
                                  "infunc": float,
                                  "format": "%5.1f",
                                  "unit": "deg"},
              
              "THETA_SKY"      : {"comment":
                                  "Position angle (east of north) (native)",
                                  "infunc": float,
                                  "format": "%+6.2f",
                                  "unit": "deg"},
              
              "THETA_J2000"    : {"comment":
                                  "Position angle (east of north) (J2000)",
                                  "infunc": float,
                                  "format": "%+6.2f",
                                  "unit": "deg"},
              
              "THETA_B1950"    : {"comment":
                                  "Position angle (east of north) (B1950)",
                                  "infunc": float,
                                  "format": "%+6.2f",
                                  "unit": "deg"},
              
              "ERRX2_WORLD"    : {"comment":
                                  "Variance of position along X-WORLD (alpha)",
                                  "infunc": float,
                                  "format": "%15e",
                                  "unit": "deg**2"},

              "ERRY2_WORLD"    : {"comment":
                                  "Variance of position along Y-WORLD (delta)",
                                  "infunc": float,
                                  "format": "%15e",
                                  "unit": "deg**2"},
              
              "ERRXY_WORLD"    : {"comment":
                                  "Covariance of position X-WORLD/Y-WORLD",
                                  "infunc": float,
                                  "format": "%15e",
                                  "unit": "deg**2"},
              
              "ERRCXX_WORLD"   : {"comment":
                                  "Cxx error ellipse parameter (WORLD units)", 
                                  "infunc": float,
                                  "format": "%12g",
                                  "unit": "deg**(-2)"},
              
              "ERRCYY_WORLD"   : {"comment":
                                  "Cyy error ellipse parameter (WORLD units)",
                                  "infunc": float,
                                  "format": "%12g",
                                  "unit": "deg**(-2)"},
              
              "ERRCXY_WORLD"   : {"comment":
                                  "Cxy error ellipse parameter (WORLD units)",
                                  "infunc": float,
                                  "format": "%12g",
                                  "unit": "deg**(-2)"},
              
              "ERRA_WORLD"     : {"comment":
                                  "World RMS position error along major axis",
                                  "infunc": float,
                                  "format": "%12g",
                                  "unit": "pixel"},

              "ERRB_WORLD"     : {"comment":
                                  "World RMS position error along minor axis",
                                  "infunc": float,
                                  "format": "%12g",
                                  "unit": "pixel"},

              "ERRTHETA_WORLD" : {"comment":
                                  "Error ellipse pos. angle (CCW/world-x)",
                                  "infunc": float,
                                  "format": "%5.1f",
                                  "unit": "deg"},
              
              "ERRTHETA_SKY"   : {"comment":
                                  "Native error ellipse pos." + \
                                  "angle (east of north)",
                                  "infunc": float,
                                  "format": "%5.1f",
                                  "unit": "deg"},

              "ERRTHETA_J2000" : {"comment":
                                  "J2000 error ellipse pos." + \
                                  "angle (east of north)",
                                  "infunc": float,
                                  "format": "%5.1f",
                                  "unit": "deg"},
              
              "ERRTHETA_B1950" : {"comment":
                                  "B1950 error ellipse pos." + \
                                  "angle (east of north)",
                                  "infunc": float,
                                  "format": "%5.1f",
                                  "unit": "deg"},

              "FWHM_WORLD"     : {"comment": "FWHM assuming a gaussian core",
                                  "infunc": float,
                                  "format": "%12g",
                                  "unit": "deg"},

              "CLASS_STAR"     : {"comment": "S/G classifier output",
                                  "infunc": float,
                                  "format": "%5.2f",
                                  "unit": ""}
              }
    

    def __init__(self, name, mode='r'):
        self.name = name
        self.mode = mode
        self.closed = True
        
        self._file = None
        self._keys = list()
        self._keys_positions = {}
        self._output = None
        self._firstline = True
        
        if self.mode != 'r':
            raise ValueError, \
                  'only read-only access is now implemented.'
        
        self._file = __builtin__.open(self.name, self.mode)
        self.closed = False
        
        # Reading header

        self._line = self._file.readline()
        if not(self._line):
            raise WrongSExtractorfileException, \
                  'not a SExtractor text catalog (empty file)'

        while (self._line):
            __ll = (self._line).replace('\n', '')
            if __ll[0] == '#':   # Still in header
                columns = __ll.split()
                if len(columns) < 3:
                    raise WrongSExtractorfileException, \
                          'not a SExtractor text catalog (invalid header)'
                name=columns[2]
                if not(name in SExtractorfile._SE_keys.keys()):
                    raise WrongSExtractorfileException, \
                          'not a SExtractor text catalog (unknown keyword %s)'\
                          % name
                self._keys_positions[name]=int(columns[1])-1
                self._keys.append(name)
            else:
                break
            self._line = self._file.readline()


        if not(self._keys):
            raise WrongSExtractorfileException, \
                  'not a SExtractor text catalog (empty header)'
            
        self._outdict = dict([(k, None) for k in self._keys])
        self._firstline = True


    def __del__(self):
        self.close()
        

    def __iter__(self):
        return self


    def next(self):
        rr = self.readline()
        if not(rr):
            raise StopIteration
        return rr


    def __nonzero__(self):
        return self._file


    def keys(self):
        "Return the list of available parameters."
        return self._keys


    def getcolumns(self):
        "Return the list of available parameters."
        return self.keys()


    def readline(self):
        """
        Read and analyse the next line of the SExtractor catalog
        and return a dictionary {'param1': value, 'param2': value, ...}.
        """
        if not(self._firstline):
            self._line = self._file.readline()

        self._firstline = False

        if not(self._line):
            return None
        __ll = (self._line).replace('\n', '')
        __values = __ll.split()
        self._outdict.update(dict(zip(self._keys, __values)))
        for i in self._keys:
            self._outdict[i] = (
                SExtractorfile._SE_keys[i]["infunc"](self._outdict[i]))
        
        return self._outdict.copy()



    def read(self):
        """
        Read the file until EOF and return a list of dictionaries.
        """
        __result = []

        __ll = self.readline()
        while __ll:
            __result.append(__ll)
            __ll = self.readline()
            
        return list(__result)


    def readlines(self):
        return self.read()


    def close(self):
        """
        Close the SExtractor file.
        """
        if self._file:
            if not(self._file.closed):
                self._file.close()
        self.closed = True


# ======================================================================

def open(name, mode='r'):
    """
    Factory function.
    Open a SExtractor file and return a SExtractor file object.
    """ 
    return SExtractorfile(name, mode)

# ======================================================================
