import os,string
from stsci.tools import teal

__taskname__ = 'imagefindpars'


def getHelpAsString():
    """ Teal help function to describe the parameters being set
        Descriptions of parameters come from .help file
    """
    helpString = teal.getHelpFileAsString(__taskname__,__file__)

    return helpString

def run(configobj=None):
    """ Imagefind parameters to control operation of built-in source extraction algorithm
    """
    pass

run.__doc__ += getHelpAsString()