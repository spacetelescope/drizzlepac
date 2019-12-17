#!/usr/bin/env python
import os
import pdb
import sys
from astropy.table import Table
from drizzlepac.devutils.comparison_tools import compare_sourcelists
from drizzlepac.devutils.comparison_tools.read_hla import read_hla_catalog


dataset = "hst_11708_02_wfc3_ir_f125w"
cattype = "sex"
orig_hla_classic_sl_name = "orig/{}_{}phot.txt".format(dataset,cattype)
imgname = "orig/hst_11708_02_wfc3_ir_f125w_drz.fits"
modcat = read_hla_catalog.read_hla_catalog(dataset, cattype=cattype, applyomega=True, multiwave=False, verbose=False, trim=False)
cat = Table.read(orig_hla_classic_sl_name, format='ascii')
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
cat['RA'] = modcat[true_ra_col_title]
cat['DEC'] = modcat[true_dec_col_title]
mod_sl_name = os.path.basename(orig_hla_classic_sl_name)
cat.write(mod_sl_name,format="ascii.csv")

return_status = compare_sourcelists.comparesourcelists([orig_hla_classic_sl_name,mod_sl_name], [imgname, imgname],plotGen="screen",diffMode="absolute",verbose=True)
