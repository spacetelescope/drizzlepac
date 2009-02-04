# python version tagger
#
# This file (version.py) is copied into every package that wants
# to store the subversion identifiers.  It lives at the top level of
# stsci_python.  If you are looking at this in some other directory
# (either lower in the hierarchy or in another package), you are looking
# at a copy.
#
# If you need to edit this file, edit the original in stsci_python, and then 
# copy it to these places:
#   stsci_python/multidrizzle/version.py 
#   stsci_python/pydrizzle/version.py    
#   stsci_python/pytools/version.py
#
# If you need to _use_ this file in your package, edit the list above to
# add your package, commit the changes, then copy the committed file into 
# your package.
#
# To use this file in your package, place this code in your setup.py file:
#
#   # gather our subversion information and save it; quit if that is all we need
#   import version
#   version.__set_svn_version__()
#   if "versiononly" in sys.argv[:2] :
#       sys.exit(0)
#
# After that, you can
#   # find out what subversion information applies to yourpackage
#   import yourpackage.svn_version
#   print yourpackage.svn_version.__svn_version__
#   print yourpackage.svn_version.__full_svn_info__
#
import os.path
import re

#
# This is the entry point.  All you need to do is call this function from your
# setup.py according to the example above.  It will create a file called
# lib/svn_version.py

def __set_svn_version__(path="./", fname='svn_version.py', fullInfo=False):
    #
    # path is the package where the version information will be stored.  Default
    # is "this package", but from a higher level package, you can specify a directory
    # of a package to process
    #
    # fname is the name of the file to store the version information in.  Never change
    # this.
    #
    # fullInfo=False suppresses the "svn info" data.  I don't know why you would
    # ever want this.
    #

    info = None
    rev = __get_svn_rev__(path)
    version_file = os.path.join(path,'lib',fname)

    # if we are unable to determine the revision, we default to leaving the
    # revision file unchanged.  Otherwise, we fill it in with whatever
    # we have

    if rev is None:
        if os.path.exists(version_file) :
            return
        revision = 'Unable to determine SVN revision'
    else:
        if ( rev == 'exported' or rev == 'unknown' ) and os.path.exists(version_file):
            return
        revision = str(rev)

    #
    if fullInfo:
        info = __get_full_info__(path)
    
    # now we can write the version information

    f = open(version_file,'w')
    f.write("\n__svn_version__ = %s\n" % repr(revision))

    # info will be a multi-line string.  We are not using repr(info)
    # for readability; the output of "svn info" can't not contain '''
    # unless you are doing something bad.
    f.write("\n__full_svn_info__ = '''\n%s'''\n" % info)
    f.close()

    
def __get_svn_rev__(path):
    m = None
    try:
        # with popen3,  stderr goes into a pipe where we ignore it, 
        # This means the user does not see errors.
        (sin, sout, serr) = os.popen3('svnversion')

        # pick up the first line of output
        m=sout.read().strip()

        # if it looks like valid svnversion output, return it
        if m == 'exported' :
            return m
        if re.match('^[0-9][0-9:]*[A-Z]*$',m) :
            return m

        # if we get here, it was not valid - that probably means
        # an error of some kind.
    except:
        pass

    # If we get this far, we did not get anything useful from svnversion.
    # This code rummages through some kinds of .svn directories to do
    # what svnversion would do.
    revision = None
    entries = os.path.join(path,'.svn','entries')
    if os.path.isfile(entries):
        f = open(entries)
        fstr = f.read()
        f.close()
        if fstr[:5] == '<?xml':  # pre 1.4
            m = re.search(r'revision="(?P<revision>\d+)"',fstr)
            if m:
                revision = int(m.group('revision'))
        else:  # non-xml entries file --- check to be sure that
            m = re.search(r'dir[\n\r]+(?P<revision>\d+)', fstr)
            if m:
                revision = int(m.group('revision'))
    return revision

def __get_full_info__(path):
    info = None
    try:
        # with popen3,  stderr goes into a pipe where we ignore it, 
        # This means the user does not see errors.
        (sin, sout, serr) = os.popen3('svn info %s' % path)

        # pick up all the lines of output
        info = [l.strip() for l in sout.readlines()]

        # if no output, there was an error and we don't know anything
        if len(info) == 0 :
            return "unknown"

        # there was output, so join it all together
        return '\n'.join(info)

    except:
        pass
    return "unknown"

