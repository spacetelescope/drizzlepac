"""

:Authors: Warren Hack, Mihai Cara

:License: :doc:`../LICENSE`

"""
from stsci.tools import teal
from drizzlepac import util
from drizzlepac import __version__

__taskname__ = 'refimagefindpars'


__doc__ = util._def_help_functions(
    locals(), module_file=__file__, task_name=__taskname__, module_doc=__doc__
)
