"""Viewer for aligned datasets"""

import glob
import os

from datetime import datetime

from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import markers

from astropy.io import fits
import numpy as np

from stwcs import wcsutil

from drizzlepac.hlautils import astrometric_utils as amutils

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
        drcnames = glob.glob(os.path.join(dirnames, '*drc.fits'))
        # Now, remove duplicates (keep drc, del drz)
        prodnames = [d if d.replace('drz.fits', 'drc.fits') not in drcnames else d.replace('drz.fits', 'drc.fits') for d in drznames]
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
        first_page = plt.figure(figsize=(8.2, 11.69))
        first_page.clf()
        txt = 'Summary of Alignment Results from {}'.format(self.parent_dir)
        first_page.text(0.05, 0.95, txt, transform=first_page.transFigure, size=font_size, ha="center")
        txt = 'Total Datasets : {}'.format(len(self.prodnames))
        first_page.text(0.1, 0.85, txt, transform=first_page.transFigure, size=font_size, ha="left")

        # compute basic stats on final WCS solutions
        apost = len([w for w in self.wcsnames if 'FIT' in w])
        defwcs = len([w for w in self.wcsnames if '-' not in w])
        apri = len(self.wcsnames) - apost - defwcs
        relonly = len([w for w in self.wcsnames if 'NONE' in w])

        txt = 'Datasets with a posteriori WCS : {}'.format(apost)
        first_page.text(0.1, 0.8, txt, transform=first_page.transFigure, size=font_size, ha="left")
        txt = 'Datasets with a posteriori relative alignment : {}'.format(relonly)
        first_page.text(0.1, 0.75, txt, transform=first_page.transFigure, size=font_size, ha="left")
        txt = 'Datasets with a priori WCS : {}'.format(apri)
        first_page.text(0.1, 0.7, txt, transform=first_page.transFigure, size=font_size, ha="left")
        txt = 'Datasets with pipeline-default WCS : {}'.format(defwcs)
        first_page.text(0.1, 0.65, txt, transform=first_page.transFigure, size=font_size, ha="left")

        now = datetime.now()
        rtime = datetime.strftime(now, "b %d %Y  %H:%M:%S")
        txt = 'Report created at {}'.format(rtime)
        first_page.text(0.1, 0.1, txt, transform=first_page.transFigure, size=font_size, ha="center")

        return first_page

    def create(self, pdfname='multipage_pdf.pdf'):

        with PdfPages(pdfname) as pdf:
            first_page = self.create_summary()
            pdf.savefig(first_page)
            plt.close()

            # Now generate a separate page for each dataset
            for p, w in zip(self.prodnames, self.wcsnames):
                result = create_product_page(p, wcsname=w)
                pdf.savefig(result)
                plt.close()


def create_product_page(prodname, zoom_size=128, wcsname=""):

    # obtain image data to display
    with fits.open(prodname) as prod:
        data = prod[1].data
        phdr = prod[0].header
        targname = phdr['targname']
        inst = phdr['instrume']
        det = phdr['detector']
        texptime = phdr['texptime']
        inexp = phdr['d001data'].split('[')[0]
        wcstype = prod[1].header['wcstype']
        wcs = wcsutil.HSTWCS(prod, ext=1)
    center = (data.shape[0] // 2, data.shape[1] // 2)
    prod_path = os.path.split(prodname)[0]

    data = np.nan_to_num(data, 0.0)

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
    fig = plt.figure(constrained_layout=True, figsize=(7, 10))
    gs = fig.add_gridspec(ncols=4, nrows=5)

    # title plots
    img_title = "{} image of {} with WCSNAME={}".format(prodname, targname, wcsname)
    plt.title(img_title, loc='center', fontsize=8)

    # Define image display
    fig_img = fig.add_subplot(gs[:3, :])
    fig_zoom = fig.add_subplot(gs[3:, 2:])
    fig_summary = fig.add_subplot(gs[3:, :2])

    # Compute display range
    dmax = (data.max() // 10)
    dscaled = np.log10(np.clip(data, -0.9, dmax) + 1)
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
    fig_img.scatter(rx, ry, marker=mstyle, alpha=0.25, c='cyan', s=3)
    fig_zoom.scatter(zx, zy, marker=mstyle, alpha=0.25, c='cyan')

    # Print summary info
    fsize = 8
    pname = os.path.split(prodname)[1]
    fig_summary.text(0.01, 0.95, "Summary for {}".format(pname), fontsize=fsize)
    fig_summary.text(0.01, 0.9, "WCSNAME: {}".format(wcsname), fontsize=fsize)
    fig_summary.text(0.01, 0.85, "TARGET: {}".format(targname), fontsize=fsize)
    fig_summary.text(0.01, 0.8, "Instrument: {}/{}".format(inst, det), fontsize=fsize)
    fig_summary.text(0.01, 0.7, "Total Exptime: {}".format(texptime), fontsize=fsize)
    fig_summary.text(0.01, 0.65, "WCSTYPE: {}".format(wcstype), fontsize=fsize)
    fig_summary.text(0.01, 0.5, "# of GAIA sources: {}".format(len(rx)), fontsize=fsize)

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
                break
        exp.close()
        fig_summary.text(0.01, 0.4, "RMS: RA={:0.3}mas, DEC={:0.3}mas".format(rms_ra, rms_dec), fontsize=fsize)
        fig_summary.text(0.01, 0.35, "# matches: {}".format(nmatch), fontsize=fsize)
        fig_summary.text(0.01, 0.3, "Matched to {} catalog".format(catalog), fontsize=fsize)

    return fig
