# write image and catalogs with and without WCS correction

from drizzlepac.devutils.comparison_tools.read_hla.read_hla_catalog import read_hla_catalog
from drizzlepac.devutils.comparison_tools.read_hla.read_hla_image import write_hla_image

if __name__ == "__main__":
    dataset = 'hst_11708_02_wfc3_ir_f125w'
    cattype = "sex"

    # without astrometry correction

    catfile = "{}_{}phot.cat".format(dataset,cattype)
    xyfile = "{}_{}phot.xy".format(dataset,cattype)
    imgfile = "{}_drz.fits".format(dataset)
    cat = read_hla_catalog(dataset, cattype=cattype, applyomega=False, multiwave=True)
    cat.write(catfile, format="ascii", overwrite=True)
    cat["ra","dec"].write(xyfile, format="ascii", overwrite=True)
    write_hla_image(dataset, imgfile, applyomega=False, overwrite=True)
    print("Wrote {} {} {}".format(catfile,imgfile,xyfile))

    # with astrometry correction

    catfile = "{}_{}phot_corr.cat".format(dataset,cattype)
    xyfile = "{}_{}phot_corr.xy".format(dataset,cattype)
    imgfile = "{}_corr.fits".format(dataset)
    cat = read_hla_catalog(dataset, cattype=cattype, multiwave=True)
    cat.write(catfile, format="ascii", overwrite=True)
    cat["ra","dec"].write(xyfile, format="ascii", overwrite=True)
    write_hla_image(dataset, imgfile, overwrite=True)
    print("Wrote {} {} {}".format(catfile,imgfile,xyfile))
