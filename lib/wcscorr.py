import os,copy
import pyfits
import numpy as np

from pytools import fileutil
from stwcs.wcsutil import altwcs


###
### WCSEXT table related keyword archive functions
### 
def init_wcscorr(input,force=False):
    """ This function will initialize the WCSCORR table if it is not already present,
        and look for WCS keywords with a prefix of 'O' as the original OPUS
        generated WCS as the initial row for the table or use the current WCS
        keywords as initial row if no 'O' prefix keywords are found. 
        
        This function will NOT overwrite any rows already present.
        
        This function works on all SCI extensions at one time.
    """

    if isinstance(input,str):
        # input must be a filename, so open as PyFITS object
        fimg = pyfits.open(input,mode='update')
        need_to_close = True
    else:
        fimg = input
        need_to_close = False

    # Verify that a WCSCORR extension does not already exist...
    for extn in fimg:
        if extn.header.has_key('extname') and extn.header['extname'] == 'WCSCORR':
            if not force:
                return
            else:
                del fimg['WCSCORR']
    # define the primary columns of the WCSEXT table with initial rows for each
    # SCI extension for the original OPUS solution
    numsci = count_ext(fimg)
    # create new table with more rows than needed initially to make it easier to
    # add new rows later
    wcsext = create_wcscorr(numrows=numsci,padding=numsci*4)
    # Assign the correct EXTNAME value to this table extension
    wcsext.header.update('TROWS',numsci*2,comment='Number of updated rows in table')
    wcsext.header.update('EXTNAME','WCSCORR',comment='Table with WCS Update history')
    
    # Now copy original OPUS values into table
    for extver in xrange(1,numsci+1):
        rowind = find_wcscorr_row(wcsext.data,{'WCS_ID':'OPUS','EXTVER':extver,'WCS_key':'O'})
        # There should only EVER be a single row for each extension with OPUS values
        rownum = np.where(rowind)[0][0]
        #print 'Archiving OPUS WCS in row number ',rownum,' in WCSCORR table for SCI,',extver

        hdr = fimg['SCI',extver].header

        # Look for old-style archived keywords for OPUS-generated values
        prefix = 'O'
        if not hdr.has_key('O'+wcs_keys[0]):
            # No such keywords, so turn off this prefix and use primary WCS keywords
            prefix = ''
        if wcsext.data.field('CRVAL1')[rownum] != 0:
            # If we find values for these keywords already in the table, do not
            # overwrite them again
            print 'WCS keywords already updated...'
            break   
        for key in wcs_keys:
            wcsext.data.field(key)[rownum] = hdr[(prefix+key)[:8]]
        # Now get any keywords from PRIMARY header needed for WCS updates
        for key in prihdr_keys:
            wcsext.data.field(key)[rownum] = fimg[0].header[key]
            
    
    # Now copy current WCS values into table
    for extver in xrange(1,numsci+1):
        hdr = fimg['SCI',extver].header

        # identify next empty row
        rowind = find_wcscorr_row(wcsext.data,selections={'wcs_id':' '})
        rows = np.where(rowind)
        if len(rows[0]) > 0:
            rownum = np.where(rowind)[0][0]
        else:
            print 'No available rows found for updating. '
        #print 'Archiving current WCS row number ',rownum,' in WCSCORR table for SCI,',extver
        
        # Update selection columns for this row with relevant values
        idctab =fileutil.osfn(fimg[0].header['idctab'])
        idcname = os.path.split(idctab)[-1]
        idcname = idcname[:idcname.find('_')]
        wcsext.data.field('WCS_ID')[rownum] = 'IDC_'+idcname
        wcsext.data.field('EXTVER')[rownum] = extver

        # Look for standard WCS keyword values
        for key in wcs_keys:            
            wcsext.data.field(key)[rownum] = hdr[key]
        # Now get any keywords from PRIMARY header needed for WCS updates
        for key in prihdr_keys:
            wcsext.data.field(key)[rownum] = fimg[0].header[key]
                
    # Append this table to the image FITS file
    fimg.append(wcsext)
    # force an update now
    # set the verify flag to 'warn' so that it will always succeed, but still
    # tell the user if PyFITS detects any problems with the file as a whole
    fimg.flush('warn')
    
    if need_to_close:
        fimg.close()

def find_wcscorr_row(wcstab, selections):
    """ Return an array of indices from the table (NOT HDU) 'wcstab' that matches the 
        selections specified by the user. 
        
        The row selection criteria must be specified as a dictionary with 
        column name as key and value(s) representing the valid desired row values.
        For example, {'wcs_id':'OPUS','extver':2}.
    """
    mask = None
    for i in selections:
        bmask = (wcstab.field(i) == selections[i])
        if mask is None:
            mask = bmask.copy()
        else:
            mask = np.logical_and(mask,bmask)
        del bmask
    return mask

def archive_wcs_file(image):
    """ Update WCSCORR table with rows for each SCI extension to record the 
        newly updated WCS keyword values.
    """
    fimg = pyfits.open(image,mode='update')
    numsci = count_ext(fimg)
    for extn in range(1,numsci+1):
        update_wcscorr(fimg,fimg['sci',extn].header)
    fimg.close()

def update_wcscorr(fimg,hdr,selections=None):
    """ Update WCSCORR table with a new row for this extension header. It 
    copies the current set of WCS keywords as a new row of the table based on 
    keyed WCSs as per Paper I Multiple WCS standard). 
    """
    # Now update the table...
    if selections is None:
        # define the WCS ID for this update
        wcs_key = altwcs.wcskeys(hdr)[-1]
        selections = {'WCS_ID':'TWEAK_'+fileutil.getDate(),'EXTVER':hdr['extver'],'WCS_key':wcs_key}
        
    # create new table for hdr and populate it with the newly updated values
    new_table = create_wcscorr()
    # ===> NOTE: need to add checks to insure that these IDs are unique
    # Update selection column values
    for key in selections:
        new_table.data.field(key)[0] = selections[key]

    for key in wcs_keys:
        new_table.data.field(key)[0] = hdr[key]

    for key in prihdr_keys:
        new_table.data.field(key)[0] = fimg[0].header[key]
            
    # Now, we need to merge this into the existing table
    old_table = fimg['WCSCORR']
    rowind = find_wcscorr_row(old_table.data,{'wcs_id':' '})
    old_nrows = np.where(rowind)[0][0]

    # check to see if there is room for the new row
    if (old_nrows + 1) > old_table.data.shape[0]:
        # if not, create a new table with 'pad_rows' new empty rows
        upd_table = pyfits.new_table(old_table.columns,nrows=old_table.data.shape[0]+pad_rows)
    else:
        upd_table = old_table
    # Now, add 
    for name in old_table.columns.names:
        upd_table.data.field(name)[old_nrows:old_nrows+1] = new_table.data.field(name)
    upd_table.header.update('TROWS',old_nrows+1)

    # replace old extension with newly updated table extension
    fimg['WCSCORR'] = upd_table
    
def restore_file_from_wcscorr(image,id='OPUS',wcskey=''):
    """ Copies the values of the WCS from the WCSCORR based on ID specified by user.
    The default will be to restore the original OPUS-derived values to the Primary WCS.
    If wcskey is specified, the WCS with that key will be updated instead. 
    """
    fimg = pyfits.open(image,mode='update')
    numsci = count_ext(fimg)
    wcs_table = fimg['WCSCORR']
    orig_rows = (wcs_table.data.field('WCS_ID') == 'OPUS')
    for extn in range(1,numsci+1):
        # find corresponding row from table
        ext_rows = (wcs_table.data.field('EXTVER') == extn)
        erow = np.where(np.logical_and(ext_rows,orig_rows))[0][0]
        for key in wcs_keys:
            tkey = key
            if 'orient' in key.lower():
                key = 'ORIENT'
            if wcskey == '':
                skey = key
            else:
                skey = key[:7]+wcskey
            fimg['sci',extn].header[skey] = wcs_table.data.field(tkey)[erow]
        for key in prihdr_keys:
            if wcskey == '':
                pkey = key
            else:
                pkey = key[:7]+wcskey
            fimg[0].header.update(pkey,wcs_table.data.field(key)[erow])
                
    # close the image now that the update has been completed.
    fimg.close()
    
def create_wcscorr(descrip=False, numrows=1, padding=0):
    """ Return the basic definitions for a WCSCORR table.
    The dtype definitions for the string columns are set to the maximum allowed so 
    that all new elements will have the same max size which will be automatically
    truncated to this limit upon updating (if needed).

    The table is initialized with rows corresponding to the OPUS solution
    for all the 'SCI' extensions.
    """ 
    trows = numrows + padding
    # define initialized arrays as placeholders for column data
    def_float64_zeros = np.array([0.0]*trows,dtype=np.float64)
    def_float64_ones = def_float64_zeros + 1.0
    def_range_1nrows = np.array(range(1,numrows+1),dtype=np.int16)
    def_float_col = {'format':'D','array':def_float64_zeros.copy()}
    def_float1_col = {'format':'D','array':def_float64_ones.copy()}
    def_str40_col = {'format':'40A','array':np.array(['']*trows,dtype="S40")}
    def_int32_col = {'format':'J','array':np.array([0]*trows,dtype=np.int32)}
    
    # If more columns are needed, simply add their definitions to this list
    col_names = [('CRVAL1',def_float_col),('CRVAL2',def_float_col),
                ('CRPIX1',def_float_col),('CRPIX2',def_float_col),
                ('CD1_1',def_float_col),('CD1_2',def_float_col),
                ('CD2_1',def_float_col),('CD2_2',def_float_col),
                ('ORIENTAT',def_float_col),('PA_V3',def_float_col),
                ('Delta_RA',def_float_col),('Delta_Dec',def_float_col),
                ('RMS_RA',def_float_col),('RMS_Dec',def_float_col),
                ('Delta_Orientat',def_float_col),('Delta_Scale',def_float1_col),
                ('NMatch',def_int32_col),('Catalog',def_str40_col)]
            
    # Define selector columns
    id_col = pyfits.Column(name='WCS_ID',format='24A',array=np.array(['OPUS']*numrows+['']*padding,dtype="S24"))
    extver_col = pyfits.Column(name='EXTVER',format='I',array=np.array(range(1,numrows+1),dtype=np.int16))
    wcskey_col = pyfits.Column(name='WCS_key',format='A',array=np.array(['O']*numrows+['']*padding,dtype="S"))
    # create list of remaining columns to be added to table
    col_list = [id_col,extver_col,wcskey_col] # start with selector columns
    for c in col_names:
        cdef = copy.deepcopy(c[1])
        col_list.append(pyfits.Column(name=c[0],format=cdef['format'],array=cdef['array']))

    if descrip:
        col_list.append(pyfits.Column(name='Descrip',format='128A',array=np.array(['Original WCS computed by OPUS']*numrows,dtype="S128")))

    # Now create the new table from the column definitions 
    newtab = pyfits.new_table(pyfits.ColDefs(col_list),nrows=trows)

    return newtab

def count_ext(fimg,extname='SCI'):
    """ Return the number of 'extname' extensions, defaulting to counting the
    number of SCI extensions.
    """
    n = 0
    for e in fimg:
        if e.header.has_key('extname') and e.header['extname'] == extname:
            n += 1
    return n
