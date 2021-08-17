
"""The module is a high-level interface to astrocut for use with HAP SVM and MVM files."""

from astrocut import fits_cut
from astropy import units as u
from astropy.units.quantity import Quantity
from astropy.coordinates import SkyCoord
from drizzlepac.haputils import cell_utils as cu
from astroquery.mast import Observations
from astropy.table import Table, vstack
from stsci.tools import logutil
import numpy as np
import math
import os
import shutil
import sys

__taskname__ = 'hapcut_utils'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout, 
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)


def mvm_id_filenames(sky_coord, cutout_size, verbose=False):
    """
    This function retrieves a table of MVM drizzled image filenames with additional
    information from the archive.  The user can then further cull the table to use as
    input to obtain a list of files from the archive.

    Parameters
    ----------
    sky_coord : str or `~astropy.coordinates.SkyCoord` object
        The position around which to cutout. It may be specified as a string ("ra dec" in degrees)
        or as the appropriate `~astropy.coordinates.SkyCoord` object.

    cutout_size : int, array-like, `~astropy.units.Quantity`
        The size of the cutout array. If ``cutout_size`` is a scalar number or a scalar
        `~astropy.units.Quantity`, then a square cutout of ``cutout_size`` will be created.
        If ``cutout_size`` has two elements, they should be in ``(ny, nx)`` order.  Scalar numbers
        in ``cutout_size`` are assumed to be in units of arcseconds. `~astropy.units.Quantity` objects
        must be in angular units.

    verbose : bool
        Default False.  If True, additional intermediate information is printed.

    Returns
    -------
    filtered_dp_table : `~astropy.table.Table` object

    """

    # If the cutout_size is not an astropy.units.Quantity object, the scalar(s)
    # are assumed to be arcseconds.  The variable must be cast as a Quantity.
    if not isinstance(cutout_size, Quantity):
        cutout_size *= u.arcsec
        cutout_size = np.atleast_1d(cutout_size)
        if len(cutout_size) == 1:
            cutout_size = np.repeat(cutout_size, 2)

    if not isinstance(sky_coord, SkyCoord):
        sky_coord = SkyCoord(sky_coord, unit="deg")

    # From HST data, Search for the list of images based upon: coordinates, search region, data
    # product type, and the instrument name (with wildcard), project (HAP), and observation
    # collection (HST).  Use the wildcard to get all the detectors for the instrument.  Multiple
    # instruments cannot be searched at the same time.  Use the diagonal of the cutout to define
    # the search radius for the archive.  Images which fall outside the desired cutout need to
    # be filtered from the solution later.
    radius = math.ceil(math.sqrt(math.pow(cutout_size.value[0], 2) + math.pow(cutout_size.value[1], 2)) / 2.0)

    # Careful - the radius must be a str or Quantity
    radius *= u.arcsec
    if verbose:
        print("Performing query for ACS images. Radius: {}.".format(radius))

    # Note: calib_level does not seem to work
    acs_query_table = Observations.query_criteria(coordinates=sky_coord,
                                                  radius=radius,
                                                  dataproduct_type="IMAGE",
                                                  instrument_name="ACS*",
                                                  project="HAP", 
                                                  calib_level=3,
                                                  obs_collection="HST") 

    if verbose:
        print("Performing query for WFC3 images.")
    wfc3_query_table = Observations.query_criteria(coordinates=sky_coord,
                                                   radius=radius,
                                                   dataproduct_type="IMAGE",
                                                   instrument_name="WFC3*",
                                                   project="HAP", 
                                                   calib_level=3,
                                                   obs_collection="HST") 

    query_table = vstack([acs_query_table, wfc3_query_table])

    # Compute the limits of the cutout region
    deg_cutout_size = cutout_size.to(u.deg)
    ra_min = sky_coord.ra.value - deg_cutout_size.value[0]
    ra_max = sky_coord.ra.value + deg_cutout_size.value[0]
    dec_min = sky_coord.dec.value - deg_cutout_size.value[1]
    dec_max = sky_coord.dec.value + deg_cutout_size.value[1]

    # Filter the output as necessary to include only MVM filenames (MVM prefix: hst_skycell).
    # Also, filter out images which are not actually in the requested cutout region as the
    # archive search had to be done using a radius.
    good_rows = []
    updated_table = None
    for i, old_row in enumerate(query_table):
        #if old_row["obs_id"].startswith("hst_skycell"):
        if old_row["obs_id"].startswith("hst_"):
            if old_row["s_ra"] >= ra_min and old_row["s_ra"] <= ra_max and \
               old_row["s_dec"] >= dec_min and old_row["s_dec"] <= dec_max:
                good_rows.append(old_row)
    updated_query_table = Table(rows=good_rows, names=query_table.colnames)

    # Catch the case where no files are found which satisfy the criteria
    if not updated_query_table:
        print("WARNING: Query for objects within {} of {} returned NO RESULTS!".format(radius, sky_coord))
        return updated_query_table

    # Get the data product list associated with the elements of the table
    if verbose:
        print("Get the product list for all entries in the query table.")
    dp_table = Observations.get_product_list(updated_query_table)

    # Filter on MVM drizzled product
    suffix = ["DRZ", "DRC"]
    if verbose:
        print("Filter the product list table for only {} filenames.".format(suffix))
    filtered_dp_table = Observations.filter_products(dp_table,
                                                     productSubGroupDescription=suffix,
                                                     extension="fits")

    if not filtered_dp_table:
        print("WARNING: No drizzle product files (DRZ/DRC) found.")
        return filtered_dp_table

    # Write the filtered data product table out to a file.  This allows for further
    # manipulation of the information before a list of filenames is distilled from
    # the table.
    # Output filename in the form: mvm_query_ra.dddd_sdec.dddd_radius_cutout.ecsv.
    #                              mvm_query_84.9208_s69.1483_71_cutout.ecsv
    ns = "s" if sky_coord.dec.degree < 0.0 else "n"
    query_filename = "mvm_query_" + str(sky_coord.ra.degree) + "_" + ns + str(abs(sky_coord.dec.degree)) + \
                     "_{:.0f}_cutout".format(radius.value) + ".ecsv" 

    if verbose:
        print("Writing out the filter product list table to {}.".format(query_filename))
    filtered_dp_table.write(query_filename, format="ascii.ecsv")

    return filtered_dp_table


def mvm_retrieve_files(products, archive=False, clobber=False, verbose=False):
    """
    This function retrieves a table of image filenames with additional information 
    from the archive.  The function returns a list of filenames available on disk.

    Parameters
    ----------
    products : str, list or Table
        Either a filename string, a list of filename strings or a Table of products
        (as returned by mvm_id_filenames). 

    archive : Boolean, optional
        Retain copies of the downloaded files in the astroquery created
        sub-directories? Default is "False".

    clobber : Boolean, optional
        Download and Overwrite existing files? Default is "False".

    verbose : bool
        Default False.  If True, additional intermediate information is printed.

    Returns
    -------
    local_files : list
        List of filenames

    Note: Code here cribbed from retrieve_obsevation in astroquery_utils module.
    """

    # FIXME BEG
    if type(products) == str:
        products = [products]

    if not isinstance(products, Table):
        new_table = Table()
        new_table['productFilename'] = products
        products = new_table.copy()
    # FIXME END

    all_images = []
    all_images = products['productFilename'].tolist()
    if not clobber:
        rows_to_remove = []
        for row_idx, row in enumerate(products):
            fname = row['productFilename']
            if os.path.isfile(fname):
                print(fname + " already exists. File download skipped.")
                rows_to_remove.append(row_idx)
        products.remove_rows(rows_to_remove)

    # Actual download of products
    manifest = Observations.download_products(products, mrp_only=False)

    if not clobber:
        for rownum in rows_to_remove[::-1]:
            if manifest:
                manifest.insert_row(rownum,
                                    vals=[all_images[rownum], "LOCAL", "None", "None"])
            else:
                return all_images

    download_dir = None
    local_files = []
    for file, file_status in zip(manifest['Local Path'], manifest['Status']):
        if file_status != "LOCAL":
            # Identify what sub-directory was created by astroquery for the
            # download
            if download_dir is None:
                download_dir = os.path.dirname(os.path.abspath(file))
            # Move or copy downloaded file to current directory
            local_file = os.path.abspath(os.path.basename(file))
            if archive:
                shutil.copy(file, local_file)
            else:
                shutil.move(file, local_file)
            # Record what files were downloaded and their current location
            local_files.append(os.path.basename(local_file))
        else:
            local_files.append(file)
    if not archive:
        # Remove astroquery created sub-directories
        shutil.rmtree('mastDownload')

    return(local_files)


def make_the_cut(input_files, sky_coord, cutout_size, verbose=False):
    """
    This function makes the actual cut in the input images. As such it isa a high-level
    interface for the `˜astrocut.cutout.fits_cut` functionality.

    Parameters
    ----------
    input_files : list
        List of fits image filenames from which to create cutouts. The SCI image is assumed to be in the
        first extension with the weight image in the second extension.

    sky_coord : str or `~astropy.coordinates.SkyCoord` object
        The position around which to cutout. It may be specified as a string ("ra dec" in degrees)
        or as the appropriate `~astropy.coordinates.SkyCoord` object.

    cutout_size : int, array-like, `~astropy.units.Quantity`
        The size of the cutout array. If ``cutout_size`` is a scalar number or a scalar
        `~astropy.units.Quantity`, then a square cutout of ``cutout_size`` will be created.
        If ``cutout_size`` has two elements, they should be in ``(ny, nx)`` order.  Scalar numbers
        in ``cutout_size`` are assumed to be in units of arcseconds. `~astropy.units.Quantity` objects
        must be in angular units.

    single_outfile : bool 
        Default True. If True, return all cutouts in a single fits file with one cutout per extension,
        if False return cutouts in individual fits files. If returing a single file the filename will 
        have the form: <cutout_prefix>_<ra>_<dec>_<size x>_<size y>.fits. If returning multiple files
        each will be named: <original filemame base>_<ra>_<dec>_<size x>_<size y>.fits.

    output_dir : str
        Default value '.'. The directory to save the cutout file(s) to.

    verbose : bool
        Default False. If True, additional intermediate information is printed.

    Returns
    -------
    response : str or list
        If single_outfile is True returns the single output filepath. Otherwise returns a list of all 
        the output filepaths.
        If memory_only is True a list of `~astropy.io.fit.HDUList` objects is returned instead of
        file name(s).

    Note: Output filename of the form: hst_skycell_cutout-p<pppp>-ra<##>d<####>-dec<n|s><##>d<####>
    """

    # Set the values for fits_cut that we are not allowing the user to modify
    CORRECT_WCS = False
    EXTENSION = [1, 2]  # SCI and WHT
    OUTPUT_PREFIX = "hst_skycell_cutout-"
    MEMORY_ONLY = True # This code will modify the output before it is written.

    # Making sure we have an array of images
    if type(input_files) == str:
        input_files = [input_files]

    # If the cutout_size is not an astropy.units.Quantity object, the scalar(s)
    # are assumed to be arcseconds.  The variable must be cast as a Quantity.
    if not isinstance(cutout_size, Quantity):
        cutout_size *= u.arcsec

    if not isinstance(sky_coord, SkyCoord):
        sky_coord = SkyCoord(sky_coord, unit="deg")

    # Call the cutout workhorse
    out_HDUList = fits_cut(input_files, sky_coord, cutout_size, correct_wcs=CORRECT_WCS,
                           extension=EXTENSION, single_outfile=True, cutout_prefix=OUTPUT_PREFIX,
                           output_dir=".", memory_only=MEMORY_ONLY, verbose=False)
     
    # Fix up the EXTNAME to be more illustrative of the actual data as
    # all of the EXTNAMEs are "CUTOUT" and all the EXTVERs are "1".
    # sci_detector_filter
    extlist = out_HDUList[1:]
    for index in range(len(extlist)):
        extlist[index].header["EXTNAME"] = extlist[index].header["O_EXT_NM"] + \
                                           "_CUTOUT_" + \
                                           extlist[index].header["ORIG_FLE"].split("_")[6].upper()

    # Finally, construct an MVM-style output filename
    #output_filename = OUTPUT_PREFIX + 
        

    #return out_HDUList
