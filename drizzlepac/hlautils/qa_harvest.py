import os
import json
import pandas as pd
import glob


def make_dataset_df(dirname, pattern='*_photometry.json'):
    """Convert dir full of JSON files into a DataFrame"""

    jpatt = os.path.join(dirname, pattern)
    hdr = None

    pdtabs = []
    for jfilename in sorted(glob.glob(jpatt)):
        print(jfilename)
        with open(jfilename) as jfile:
            resids = json.load(jfile)
        pdindx = None
        if hdr is None:
            hdr = resids['header']
        rootname = hdr['FILENAME'].replace('.fits', '')
        k = resids['data'].keys()
        for key in k:
            dat = resids['data'][key]['data']

            det = dat['detector']
            filtname = dat['filter_name']
            del dat['detector']
            del dat['filter_name']
            if pdindx is None:
                pdindx = '-'.join([rootname, det, filtname])
            for dk in dat.keys():
                for di in dat[dk].keys():
                    hdr.update(dict([('-'.join([dk.split(" ")[2], di]), dat[dk][di])]))

        pdtabs.append(pd.DataFrame(hdr, index=[pdindx]))
    if len(pdtabs) == 0:
        allpd = None
    else:
        allpd = pd.concat(pdtabs)
    return allpd


def make_master_df(dirname, pattern='*.json', num=None):
    dirs = sorted(glob.glob(os.path.join(dirname, '*')))
    allpd = None
    for d in dirs[:num]:
        print(d)
        pdtab = make_dataset_df(d, pattern=pattern)
        if pdtab is not None:
            if allpd is None:
                allpd = pdtab
            else:
                allpd = allpd.append(pdtab)

    return allpd