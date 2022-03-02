"""
:Authors: Warren Hack, Mihai Cara

:License: :doc:`LICENSE`

"""
import os, string
from stsci.tools import teal
from . import util
from . import __version__

__taskname__ = 'imagefindpars'


__doc__ = util._def_help_functions(
    locals(), module_file=__file__, task_name=__taskname__, module_doc=__doc__
)

