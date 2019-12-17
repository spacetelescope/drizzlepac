#!/usr/bin/env python
import os
import pdb
import sys
from astropy.table import Table
from drizzlepac.devutils.comparison_tools import compare_sourcelists
from drizzlepac.devutils.comparison_tools.read_hla import read_hla_catalog


def correct_hla_classic_ra_dec(orig_hla_classic_sl_name,cattype):
    """
    This subroutine runs Rick White's read_hla_catalog script to convert the HLA Classic sourcelist RA and DEC values
    into same reference frame used by the HAP sourcelists. A new version of the input file with the converted RA and DEC
    values is written to the current working directory named <INPUT SOURCELIST NAME>_corrected.txt.

    Parameters
    ----------
    orig_hla_classic_sl_name : string
        name of the HLA Classic sourcelist whose RA and DEC values will be converted.

    cattype : string
        HLA Classic catalog type. Either 'sex' (source extractor) or 'dao' (DAOphot).

    Returns
    -------
    mod_sl_name : string
        Name of the new version of the input file with the converted RA and DEC values
    """
    mod_sl_name = os.path.basename(orig_hla_classic_sl_name)

    # Execute read_hla_catalog.read_hla_catalog() to convert RA and Dec values
    dataset = mod_sl_name.replace("_{}phot.txt".format(cattype),"")
    modcat = read_hla_catalog.read_hla_catalog(dataset, cattype=cattype, applyomega=True, multiwave=False, verbose=False, trim=False)


    # Identify RA and Dec column names in the new catalog table object
    for ra_col_title in ["ra", "RA", "ALPHA_J2000", "alpha_j2000"]:
        if ra_col_title in modcat.colnames:
            true_ra_col_title = ra_col_title
            print("RA Col_name: {}".format(true_ra_col_title))
            break
    for dec_col_title in ["dec", "DEC", "Dec", "DELTA_J2000", "delta_j2000"]:
        if dec_col_title in modcat.colnames:
            true_dec_col_title = dec_col_title
            print("DEC Col_name: {}".format(true_dec_col_title))
            break

    # get HLA Classic sourcelist data, replace existing RA and Dec column data with the converted RA and Dec column data
    cat = Table.read(orig_hla_classic_sl_name, format='ascii')
    cat['RA'] = modcat[true_ra_col_title]
    cat['DEC'] = modcat[true_dec_col_title]

    # Write updated version of HLA Classic sourcelist to current working directory
    mod_sl_name = mod_sl_name.replace(".txt","_corrected.txt")
    print("RA/DEC corrected version of HLA Classic file {} written to {}.".format(orig_hla_classic_sl_name,mod_sl_name))
    cat.write(mod_sl_name,format="ascii.csv")

    return mod_sl_name



if __name__ == "__main__":
    orig_hla_classic_sl_name = "orig/hst_11708_02_wfc3_ir_f125w_sexphot.txt"
    imgname = "orig/hst_11708_02_wfc3_ir_f125w_drz.fits"
    mod_sl_name = correct_hla_classic_ra_dec(orig_hla_classic_sl_name,'sex')

    return_status = compare_sourcelists.comparesourcelists([orig_hla_classic_sl_name, mod_sl_name], [imgname, imgname],
                                                           plotGen="screen", diffMode="absolute", verbose=True)