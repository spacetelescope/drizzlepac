#
#   Authors: Warren Hack, Ivo Busko, Christopher Hanley
#   Program: mdzhandler.py
#   Purpose: Module that handles the MDRIZTAB reference file.
import string

import pyfits
import numpy as np

from pytools import fileutil


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
    if _header.has_key('MDRIZTAB'):
        _tableName = _header['MDRIZTAB']
    else:
        raise KeyError, "No MDRIZTAB found in file " + _fileName

    _tableName = fileutil.osfn(_tableName)

    # Now get the filters from the primary header.
    _filters = fileutil.getFilterNames(_header)

    # Open MDRIZTAB file.
    try:
        _mdriztab = pyfits.open(_tableName)
    except:
        raise IOError,"MDRIZTAB table '%s' not valid!" % _tableName

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
    print '- MDRIZTAB: MultiDrizzle parameters read from row %s.'%(_row+1)

    mpars = _mdriztab[1].data[_row]
    _mdriztab.close() 

    return _interpretMdriztabPars(mpars)

def _getRowsByFilter(table, filters):
    rows = []
    for i in xrange(table[1].data.shape[0]):
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
    for indx in xrange(len(rec.array.names)):
        # ... get the name, format, and value.
        _name = rec.array.names[indx]
        _format = rec.array.formats[indx]
        _value = rec.field(_name)
         
        # Translate names from MDRIZTAB columns names to 
        # input parameter names found in IRAF par file.
        #
        if _name.find('final') > -1: _name = 'driz_'+_name
        elif _name == 'subsky': _name = 'skysub'
        elif _name == 'crbitval': _name = 'crbit'
        elif _name == 'readnoise': _name = 'rdnoise'
                    
        # We do not care about the first two columns at this point
        # as they are only used for selecting the rows
        if _name != 'filter' and _name != 'numimages':
            # start by determining the format type of the parameter
            _fmt = findFormat(_format)
            
            # Based on format type, apply proper conversion/cleaning
            if (_fmt == 'a') or (_fmt == 'A'):
                _val = cleanBlank(_value)
                if _val == None: _val = ''
            elif (_format == 'i1') or (_format=='1L'):
                _val = toBoolean(_value)
            elif (_format == 'i4') or (_format == '1J'):
                _val = cleanInt(_value)
            elif (_format == 'f4') or (_format == '1E'):
                _val = cleanNaN(_value)
            else:
                print 'MDRIZTAB column ',_name,' has unrecognized format',_format 
                raise ValueError
            if _name.find('fillval') > -1 and _val == None:
                _val = 'INDEF'
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
    return value

def cleanInt(value):
    # THIS MAY BE MACHINE-DEPENDENT !!!
    if value == -2147483647:
    # Try to use 'sys.maxint' as value (WJH)
    #if value == -sys.maxint:
        return None
    return value

def cleanBlank(value):
    if value.strip() == '':
        return None
    return value

def findFormat(format):
    # Parses record array format string for type
    _fmt = None
    for ltr in string.letters:
        if format.find(ltr) > -1:
            _fmt = ltr
            break
    return _fmt
