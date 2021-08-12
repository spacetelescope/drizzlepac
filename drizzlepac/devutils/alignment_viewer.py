"""Viewer for aligned datasets

This code can be used to summarize a set of results after running runastrodriz.
The results will be a multi-page PDF document with a title page that summarizes
the number of tests found along with stats on the number of various types of
WCS solutions in the final drizzle products.  A separate page will then be generated
for each dataset with either a DRZ or DRC product.  That page will contain:
    - a full-frame view of the final product image
    - a view of the centeral 256x256 region from the final product image
    - positions of any GAIA catalog sources overplotted on both images
    - a summary of the product; including:
        * WCSNAME
        * WCSTYPE value
        * Target name
        * Instrument configuration
        * total exposure time for the final product image
        * number of GAIA sources in the field
        * If an aposteriori fit was found, statistics on the fit such as shift, rot, scale and RMS

If you are in a directory with 2 sub-directories 'iaal01hxq' and 'iacs01t4q'
containing the standard-pipeline processing results for those 2 datasets, then
the report can be generated in the Python shell (like ipython) using:

>>> from drizzlepac.devutils import alignment_viewer
>>> d = alignment_viewer.Datasets('.', num_levels=1)
>>> d.create(pdfname='pipeline_results.pdf')

The "num_levels" parameter provides a mechanism for looking for processed products in that many directories below the "parent_dir" (starting directory).

Alternatively, a single dataset can be summarized by simply using:

>>> from drizzlepac.devutils import alignment_viewer
>>> figure = create_product_page("iaal01hxq/iaal01hxq_drc.fits")
>>> figure.savefig("iaal01hxq_summmary.pdf")  # to write out a PDF file


"""

import glob
import os
import json

from datetime import datetime

from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import markers

from astropy.io import fits
import numpy as np

from stwcs import wcsutil

from drizzlepac.haputils import astrometric_utils as amutils

plt.ioff()

class Datasets:
    def __init__(self, parent_dir, num_levels=1):
        # determine how many datasets are available to review
        levels = os.path.sep.join('*' * num_levels)
        dirnames = os.path.join(parent_dir, levels)
        self.parent_dir = os.path.split(os.path.abspath(parent_dir))[1]
        self.dirnames = dirnames

        # Determine list of products to be reviewed
        drznames = glob.glob(os.path.join(dirnames, '*drz.fits'))
        prodnames = glob.glob(os.path.join(dirnames, '*drc.fits'))
        # Now, add any DRZ-only products to the list of DRC products
        for drz in drznames:
            drzc = drz.replace('drz.fits', 'drc.fits')
            if drzc not in prodnames:
                prodnames.append(drz)

        # Insure that pytest-temp directories are not included either
        self.prodnames = [d for d in prodnames if 'current' not in d]

        # Extract key information from products for statistical use
        #
        # Define values to be included in summary
        self.wcsnames = []
        self.instruments = []
        self.detectors = []
        self.targnames = []
        self.texptimes = []
        print("Processing {} products".format(len(self.prodnames)))
        for p in self.prodnames:
            with fits.open(p) as prod:
                phdu = prod[0].header
                sci = prod[1].header if len(prod) > 1 else phdu
                self.wcsnames.append(sci['wcsname'])

    def create_summary(self):
        font_size = 12
        first_page = plt.figure(figsize=(8.5, 11))
        first_page.clf()
        txt = 'Summary of Alignment Results from {}'.format(self.parent_dir)
        first_page.text(0.5, 0.95, txt, transform=first_page.transFigure, size=font_size, ha="center")
        txt = 'Total Datasets : {}'.format(len(self.prodnames))
        first_page.text(0.1, 0.85, txt, transform=first_page.transFigure, size=font_size, ha="left")

        # compute basic stats on final WCS solutions
        apost = len([w for w in self.wcsnames if 'FIT' in w])
        defwcs = len([w for w in self.wcsnames if '-' not in w or 'None' in w or w.strip()==''])
        apri = len(self.wcsnames) - apost - defwcs
        relonly = len([w for w in self.wcsnames if 'NONE' in w ])

        txt = 'Datasets with a posteriori WCS : {}'.format(apost)
        first_page.text(0.1, 0.8, txt, transform=first_page.transFigure, size=font_size, ha="left")
        txt = 'Datasets with a posteriori relative alignment : {}'.format(relonly)
        first_page.text(0.1, 0.75, txt, transform=first_page.transFigure, size=font_size, ha="left")
        txt = 'Datasets with a priori WCS : {}'.format(apri)
        first_page.text(0.1, 0.7, txt, transform=first_page.transFigure, size=font_size, ha="left")
        txt = 'Datasets with pipeline-default WCS : {}'.format(defwcs)
        first_page.text(0.1, 0.65, txt, transform=first_page.transFigure, size=font_size, ha="left")

        now = datetime.now()
        rtime = datetime.strftime(now, "%b %d %Y  %H:%M:%S")
        txt = 'Report created at {}'.format(rtime)
        first_page.text(0.1, 0.1, txt, transform=first_page.transFigure, size=font_size, ha="left")

        return first_page

    def create(self, pdfname='multipage_pdf.pdf', num_datasets=None):
        if num_datasets is not None:
            prodnames = self.prodnames[:num_datasets]
            wcsnames = self.wcsnames[:num_datasets]
        else:
            prodnames = self.prodnames
            wcsnames = self.wcsnames

        json_summary = {}
        with PdfPages(pdfname) as pdf:
            first_page = self.create_summary()
            pdf.savefig(first_page)
            plt.close()

            plt.ioff()
            # Now generate a separate page for each dataset
            for p, w in zip(prodnames, wcsnames):
                result, summary = create_product_page(p, wcsname=w)
                if result is not None:
                    pdf.savefig(result)
                    plt.close()
                    json_summary[os.path.basename(p)] = summary
            plt.ion()

        with open(pdfname.replace('.pdf', '_summary.json'), 'w') as jsonfile:
            json.dump(json_summary, jsonfile)


def create_product_page(prodname, zoom_size=128, wcsname="",
                        gcolor='magenta', fsize=8):
    """Create a matplotlib Figure() object which summarizes this product FITS file."""

    # obtain image data to display
    with fits.open(prodname) as prod:
        data = prod[1].data
        phdr = prod[0].header
        if 'texptime' not in phdr:
            return None, None
        targname = phdr['targname']
        inst = phdr['instrume']
        det = phdr['detector']
        texptime = phdr['texptime']
        inexp = phdr['d001data'].split('[')[0]
        wcstype = prod[1].header['wcstype']
        wcs = wcsutil.HSTWCS(prod, ext=1)
        hdrtab = prod['hdrtab'].data
        filters = ';'.join([phdr[f] for f in phdr['filter*']])
        dateobs = phdr['date-obs']  # human-readable date
        expstart = phdr['expstart']  # MJD float value
        asnid = phdr.get('asn_id', '')

    center = (data.shape[0] // 2, data.shape[1] // 2)
    prod_path = os.path.split(prodname)[0]

    data = np.nan_to_num(data, nan=0.0)

    # Get GAIA catalog
    refcat = amutils.create_astrometric_catalog(prodname, existing_wcs=wcs,
                                                output=None)
    refx, refy = wcs.all_world2pix(refcat['RA'], refcat['DEC'], 0)
    # Remove points outside the full-size image area
    rx = []
    ry = []
    zx = []
    zy = []
    for x, y in zip(refx, refy):
        if 0 < x < data.shape[1] and 0 < y < data.shape[0]:
            rx.append(x)
            ry.append(y)
        if -zoom_size < x - center[1] < zoom_size and -zoom_size < y - center[0] < zoom_size:
            zx.append(x - center[1] + zoom_size)
            zy.append(y - center[0] + zoom_size)

    # Define subplot regions on page
    fig = plt.figure(constrained_layout=True, figsize=(8.5, 11))
    gs = fig.add_gridspec(ncols=4, nrows=5)

    # title plots
    rootname = os.path.basename(prodname)
    img_title = "{} image of {} with WCSNAME={}".format(rootname, targname, wcsname)
    fig.suptitle(img_title, ha='center', va='top', fontsize=fsize)

    # Define image display
    fig_img = fig.add_subplot(gs[:3, :])
    fig_zoom = fig.add_subplot(gs[3:, 2:])
    fig_summary = fig.add_subplot(gs[3:, :2])

    # Compute display range
    dmax = (data.max() // 10) if data.max() <= 1000. else 100
    dscaled = np.log10(np.clip(data, -0.1, dmax) + 0.10001)
    # identify zoom region around center of data
    zoom = dscaled[center[0] - zoom_size:center[0] + zoom_size,
                   center[1] - zoom_size:center[1] + zoom_size]

    # display full image
    fig_img.imshow(dscaled, cmap='gray', origin='lower')
    # display zoomed section
    fig_zoom.imshow(zoom, cmap='gray', origin='lower')

    # define markerstyle
    mstyle = markers.MarkerStyle(marker='o')
    mstyle.set_fillstyle('none')
    # plot GAIA sources onto full size image
    fig_img.scatter(rx, ry, marker=mstyle, alpha=0.35, c=gcolor, s=3)
    fig_zoom.scatter(zx, zy, marker=mstyle, alpha=0.35, c=gcolor)

    # Print summary info
    pname = os.path.split(prodname)[1]
    fig_summary.text(0.01, 0.95, "Summary for {}".format(pname), fontsize=fsize)
    fig_summary.text(0.01, 0.9, "WCSNAME: {}".format(wcsname), fontsize=fsize)
    fig_summary.text(0.01, 0.85, "TARGET: {}".format(targname), fontsize=fsize)
    fig_summary.text(0.01, 0.8, "Instrument: {}/{}".format(inst, det), fontsize=fsize)
    fig_summary.text(0.01, 0.75, "Filters: {}".format(filters), fontsize=fsize)
    fig_summary.text(0.01, 0.7, "Total Exptime: {}".format(texptime), fontsize=fsize)
    fig_summary.text(0.01, 0.65, "WCSTYPE: {}".format(wcstype), fontsize=fsize)
    fig_summary.text(0.01, 0.5, "Total # of GAIA sources: {}".format(len(refx)), fontsize=fsize)
    fig_summary.text(0.01, 0.45, "# of GAIA matches: {}".format(len(rx)), fontsize=fsize)

    # Get extended information about observation
    hdrtab_cols = hdrtab.columns.names
    mtflag = get_col_val(hdrtab, 'mtflag', default="")
    gyromode = get_col_val(hdrtab, 'gyromode', default='N/A')


    # populate JSON summary info
    summary = dict(wcsname=wcsname, targname=targname, asnid=asnid,
                    dateobs=dateobs, expstart=expstart,
                    instrument=(inst, det), exptime=texptime,
                    wcstype=wcstype, num_gaia=len(refx), filters=filters,
                    rms_ra=-1, rms_dec=-1, nmatch=-1, catalog="")
    obs_kws = ['gyromode', 'fgslock', 'aperture', 'mtflag', 'subarray',
                'obstype', 'obsmode', 'scan_typ', 'photmode']
    for kw in obs_kws:
        summary[kw] = get_col_val(hdrtab, kw, default="")


    if 'FIT' in wcsname:
        # Look for FIT RMS and other stats from headerlet
        exp = fits.open(os.path.join(prod_path, inexp))
        for ext in exp:
            if 'extname' in ext.header and ext.header['extname'] == 'HDRLET' \
                and ext.header['wcsname'] == wcsname:
                hdrlet = ext.headerlet
                rms_ra = hdrlet[0].header['rms_ra']
                rms_dec = hdrlet[0].header['rms_dec']
                nmatch = hdrlet[0].header['nmatch']
                catalog = hdrlet[0].header['catalog']

                fit_vals = dict(rms_ra=rms_ra, rms_dec=rms_dec, nmatch=nmatch, catalog=catalog)
                summary.update(fit_vals)
                break
        exp.close()
        try:
            fig_summary.text(0.01, 0.4, "RMS: RA={:0.3f}mas, DEC={:0.3f}mas".format(rms_ra, rms_dec), fontsize=fsize)
            fig_summary.text(0.01, 0.35, "# matches: {}".format(nmatch), fontsize=fsize)
            fig_summary.text(0.01, 0.3, "Matched to {} catalog".format(catalog), fontsize=fsize)
        except:
            fig_summary.text(0.01, 0.35, "No MATCH to GAIA")
            print("Data without a match to GAIA: {},{}".format(inexp, wcsname))


    return fig, summary

def get_col_val(hdrtab, keyword, default=None):
    val = hdrtab[0][keyword.upper()] if keyword.upper() in hdrtab.columns.names else default
    if isinstance(val, (bool, np.bool_)):
        val = str(val)
    return val
