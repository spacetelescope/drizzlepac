import os,string
from stsci.tools import teal
import tweakreg

__taskname__ = 'imagefindpars'

# This is specifically NOT intended to match the package-wide version information.
__version__ = tweakreg.__version__
__vdate__ = tweakreg.__vdate__

def getHelpAsString(docstring=False):
    """ Teal help function to describe the parameters being set
        Descriptions of parameters come from .help file
    """
    install_dir = os.path.dirname(__file__)
    htmlfile = os.path.join(install_dir,'htmlhelp',__taskname__+'.html')
    helpfile = os.path.join(install_dir,__taskname__+'.help')
    if docstring or (not docstring and not os.path.exists(htmlfile)):
        helpString = __taskname__+' Version '+__version__+' updated on '+__vdate__+'\n\n'
        if os.path.exists(helpfile):
            helpString += teal.getHelpFileAsString(__taskname__,__file__)
    else:
        helpString = 'file://'+htmlfile

    return helpString


def help():
    print getHelpAsString(docstring=True)

__doc__ = getHelpAsString(docstring=True)
