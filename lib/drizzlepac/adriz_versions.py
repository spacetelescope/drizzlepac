try:
    from . import svn_version
except:
    svn_version = None

if svn_version:
    sversion = 'dev'+svn_version.__svn_version__
else:
    sversion = ''

__version__ = '1.0.6'
__full_version__ = __version__+sversion
__vdate__ = '13-Aug-2012'

def main():
    print '%s(%s)'%(__version__,__vdate__)

if __name__ == "__main__":
    main()
