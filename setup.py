
# This setup.py does not follow our usual pattern, because it is
# the setup for pytools -- it can't use pytools to install itself..

import sys

# We can't import this because we are not installed yet, so
# exec it instead.  We only want it to initialize itself,
# so we don't need to keep the symbol table.

syms = { }
f=open("lib/stsci_distutils_hack.py","r")
exec f in syms
f.close()

syms['__set_svn_version__']()
syms['__set_setup_date__']()

if "version" in sys.argv :
    sys.exit(0)

from distutils.core import setup

from defsetup import setupargs, pkg

setup(
    name =              pkg,
    packages =          [ pkg ],
    package_dir =       { pkg : 'lib' },
    **setupargs
    )
