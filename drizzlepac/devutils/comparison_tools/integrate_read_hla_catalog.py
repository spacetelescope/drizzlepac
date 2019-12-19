#!/usr/bin/env python
import os
import pdb
import sys
from astropy.table import Table
from drizzlepac.devutils.comparison_tools import compare_sourcelists
from drizzlepac.devutils.comparison_tools.read_hla import read_hla_catalog


def correct_hla_classic_ra_dec(orig_hla_classic_sl_name, cattype, log_level):
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

    log_level : int
        The desired level of verboseness in the log statements displayed on the screen and written to the .log file.

    Returns
    -------
    mod_sl_name : string
        Name of the new version of the input file with the converted RA and DEC values
    """

    mod_sl_name = os.path.basename(orig_hla_classic_sl_name)

    # Execute read_hla_catalog.read_hla_catalog() to convert RA and Dec values
    dataset = mod_sl_name.replace("_{}phot.txt".format(cattype),"")
    modcat = read_hla_catalog.read_hla_catalog(dataset, cattype=cattype, applyomega=True, multiwave=False,
                                               verbose=True, trim=False, log_level=log_level)


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
    hla_sourcelist_name = "hla_classic/hst_11665_06_wfc3_uvis_f555w_daophot.txt"
    hla_imgname = "hla_classic/hst_11665_06_wfc3_uvis_f555w_drz.fits"
    hap_sourcelist_name = "hst_11665_06_wfc3_uvis_f555w_ib4606_point-cat.ecsv"
    hap_imgname = "hst_11665_06_wfc3_uvis_f555w_ib4606_drc.fits"
    updated_hla_sourcelist_name = correct_hla_classic_ra_dec(hla_sourcelist_name,'dao',20)

    return_status = compare_sourcelists.comparesourcelists([updated_hla_sourcelist_name,hap_sourcelist_name], [hla_imgname, hap_imgname],
                                                           plotGen="screen", diffMode="absolute", verbose=True)