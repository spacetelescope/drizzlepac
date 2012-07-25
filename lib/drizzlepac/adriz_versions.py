try:
    from . import svn_version
except:
    svn_version = None

if svn_version:
    sversion = 'dev'+svn_version.__svn_version__
else:
    sversion = ''

__version__ = '1.0.4'
__full_version__ = __version__+sversion
__vdate__ = '25-July-2012'
