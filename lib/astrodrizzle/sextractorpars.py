import os, string
from stsci.tools import teal

__taskname__ = 'sextractorpars'


def getHelpAsString():
    """ Teal help function to describe the parameters being set
        Descriptions of parameters come from .help file
    """
    helpString = teal.getHelpFileAsString(__taskname__,__file__)

    return helpString

def run(configobj=None):
    """ Sextractor parameters which can be set by `tweakreg`_.
    """
    pass

run.__doc__ += getHelpAsString()