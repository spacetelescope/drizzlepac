"""
This module supports the interpretation of the ``MDRIZTAB`` for
processing as used in the pipeline.

:Authors: Warren Hack, Ivo Busko, Christopher Hanley

:License: :doc:`LICENSE`

"""
import string, os

from astropy.io import fits
import numpy as np

from stsci.tools import fileutil

def getMdriztabParameters(files):
    """ Gets entry in MDRIZTAB where task parameters live.
        This method returns a record array mapping the selected
        row.
    """

    # Get the MDRIZTAB table file name from the primary header.
    # It is gotten from the first file in the input list. No
    # consistency checks are performed.
    _fileName = files[0]
    _header = fileutil.getHeader(_fileName)
    if 'MDRIZTAB' in _header:
        _tableName = _header['MDRIZTAB']
    else:
        raise KeyError("No MDRIZTAB found in file " + _fileName)

    _tableName = fileutil.osfn(_tableName)

    # Now get the filters from the primary header.
    _filters = fileutil.getFilterNames(_header)

    # Specifically check to see whether the MDRIZTAB file can be found
    mtab_path = os.path.split(_tableName)[0] # protect against no path given for _tableName
    if mtab_path and not os.path.exists(mtab_path): # check path first, if given
        raise IOError("Directory for MDRIZTAB '%s' could not be accessed!"%mtab_path)
    if not os.path.exists(_tableName): # then check for the table itself
        raise IOError("MDRIZTAB table '%s' could not be found!"%_tableName)

    # Open MDRIZTAB file.
    try:
        _mdriztab = fits.open(_tableName, memmap=False)
    except:
        raise IOError("MDRIZTAB table '%s' not valid!" % _tableName)

    # Look for matching rows based on filter name. If no
    # match, pick up rows for the default filter.
    _rows = _getRowsByFilter(_mdriztab, _filters)
    if _rows == []:
        _rows = _getRowsByFilter(_mdriztab, 'ANY')

    # Now look for the row that matches the number of images.
    # The logic below assumes that rows for a given filter
    # are arranged in ascending order of the 'numimage' field.
    _nimages = len(files)
    _row = 0
    for i in _rows:
        _numimages = _mdriztab[1].data.field('numimages')[i]
        if _nimages >= _numimages:
            _row = i
    print('- MDRIZTAB: AstroDrizzle parameters read from row %s.'%(_row+1))

    mpars = _mdriztab[1].data[_row]
    _mdriztab.close()

    interpreted = _interpretMdriztabPars(mpars)

    if "staticfile" in interpreted:
        interpreted.pop("staticfile")

    return interpreted

def _getRowsByFilter(table, filters):
    rows = []
    for i in range(table[1].data.shape[0]):
        _tfilters = table[1].data.field('filter')[i]
        if _tfilters == filters:
            rows.append(i)
    return rows

def _interpretMdriztabPars(rec):
    """
    Collect task parameters from the MDRIZTAB record and
    update the master parameters list with those values

    Note that parameters read from the MDRIZTAB record must
    be cleaned up in a similar way that parameters read
    from the user interface are.
    """
    tabdict = {}
    # for each entry in the record...
    for indx in range(len(rec.array.names)):
        # ... get the name, format, and value.
        _name = rec.array.names[indx]
        _format = rec.array.formats[indx]
        _value = rec.field(_name)

        # Translate names from MDRIZTAB columns names to
        # input parameter names found in IRAF par file.
        #
        #if _name.find('final') > -1: _name = 'driz_'+_name
        if _name in ['shiftfile','mdriztab']:
            continue
        drizstep_names = ['driz_sep_','final_']
        if _name in ['refimage','bits']:
            for dnames in drizstep_names:
                tabdict[dnames+_name] = _value
            continue
        if _name in ['driz_sep_bits','final_bits']:
            tabdict[_name] = str(_value)
            continue
        if _name == 'coeffs':
            _val = True
            if _value in ['INDEF',None,"None",'',' ']: _val = False
            tabdict[_name] = _val
            continue

        par_table = {'subsky':'skysub','crbitval':'crbit','readnoise':'rdnoise'}
        if _name in par_table:
            _name = par_table[_name]

        # We do not care about the first two columns at this point
        # as they are only used for selecting the rows
        if _name != 'filter' and _name != 'numimages':
            # start by determining the format type of the parameter
            _fmt = findFormat(_format)

            # Based on format type, apply proper conversion/cleaning
            if (_fmt == 'a') or (_fmt == 'A'):
                _val = cleanBlank(_value)
                if _val is None:
                    _val = ''

            elif (_format == 'i1') or (_format=='1L'):
                _val = toBoolean(_value)

            elif (_format == 'i4') or (_format == '1J'):
                _val = cleanInt(_value)

            elif ('E' in _format) or (_format == 'f4') :
                _val = cleanNaN(_value)

            else:
                print('MDRIZTAB column ',_name,' has unrecognized format',_format)
                raise ValueError
            if _name in ['ra','dec']:
                for dnames in drizstep_names:
                    tabdict[dnames+_name] = _val

            else:
                tabdict[_name] = _val

    return tabdict

def toBoolean(flag):
    if (flag == 1):
        return True
    return False

def cleanNaN(value):
    a = np.array(value)
#    b = np.where(np.isnan(a))
    if np.any(np.isnan(a)): return None
    return float(value)

def cleanInt(value):
    # THIS MAY BE MACHINE-DEPENDENT !!!
    if value == -2147483647:
    # Try to use 'sys.maxint' as value (WJH)
    #if value == -sys.maxint:
        return None
    return int(value)

def cleanBlank(value):
    if value.strip() == '':
        return None
    return value

def findFormat(format):
    # Parses record array format string for type
    _fmt = None
    for ltr in string.ascii_letters:
        if format.find(ltr) > -1:
            _fmt = ltr
            break
    return _fmt
