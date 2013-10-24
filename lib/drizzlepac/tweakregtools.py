import string,os

import numpy as np
import stsci.ndimage as ndimage

from stsci.tools import asnutil, irafglob, parseinput, fileutil
import pyfits
import astrolib.coords as coords


import stsci.imagestats as imagestats

from . import findobj
from . import cdriz

def generate_reg_files(refimg, wcs_reg, images, chip_reg, append=False,
                       outpath='./regions'):
