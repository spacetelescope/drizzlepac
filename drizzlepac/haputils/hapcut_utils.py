"""The module is a high-level interface to astrocut for use with HAP SVM and MVM files."""

from astrocut import fits_cut
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table, vstack, unique
from astropy.units.quantity import Quantity
from astroquery.mast import Observations
from drizzlepac.haputils import cell_utils as cu
from pprint import pprint
from stsci.tools import logutil

import astrocut
import copy
import glob
import math
import numpy as np
import os
import shutil
import sys

__taskname__ = 'hapcut_utils'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger("hapcut", level=logutil.logging.NOTSET, stream=sys.stdout, 
                            filename="hapcut_utility.log", format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)


def mvm_id_filenames(sky_coord, cutout_size, log_level=logutil.logging.INFO):
    """
    This function retrieves a table of MVM drizzled image filenames with additional
    information from the archive.  The user can then further cull the table to use as
    input to obtain a list of files from the archive.  This function will return filter-level
    products. At this time, both ACS and WFC3 are searched by default.

    Parameters
    ----------
    sky_coord : str or `astropy.coordinates.SkyCoord <https://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html>`_ object
        The position around which to cutout. It may be specified as a string ("ra dec" in degrees)
        or as the appropriate
        `astropy.coordinates.SkyCoord <https://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html>`_
        object.

    cutout_size : int, array-like, `astropy.units.Quantity <https://docs.astropy.org/en/stable/api/astropy.units.quantity.Quantity.html>`_
        The size of the cutout array. If ``cutout_size`` is a scalar number or a scalar
        `astropy.units.Quantity <https://docs.astropy.org/en/stable/api/astropy.units.quantity.Quantity.html>`_,
        then a square cutout of ``cutout_size`` will be created.
        If ``cutout_size`` has two elements, they should be in ``(ny, nx)`` order.  Scalar numbers
        in ``cutout_size`` are assumed to be in units of arcseconds.
        `astropy.units.Quantity <https://docs.astropy.org/en/stable/api/astropy.units.quantity.Quantity.html>`_ objects
        must be in angular units.

    log_level : int, optional
        The desired level of verbosity in the log statements displayed on the screen and written to the
        .log file. Default value is 20, or 'info'.

    Returns
    -------
    final_table : `astropy.table.Table <https://docs.astropy.org/en/stable/api/astropy.table.Table.html>`_ object

    This utility also writes an output ECSV file version of the in-memory filtered data product table,
    final_table.  The output filename is in the form: 
    mvm_query-ra<###>d<####>-dec<n|s><##>d<####>_<radius>_cutout.ecsv 
    (e.g., mvm_query-ra84d8208-decs69d8516_354_cutout.ecsv, where radius has been computed from the
    cutout dimensions.
    """

    # set logging level to user-specified level
    log.setLevel(log_level)

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
    log.info("Radius for query: {}.".format(radius))
    log.info("Performing query for ACS images.")

    acs_query_table = Observations.query_criteria(coordinates=sky_coord,
                                                  radius=radius,
                                                  dataproduct_type="IMAGE",
                                                  instrument_name="ACS*",
                                                  project="HAP", 
                                                  obs_collection="HST") 

    log.info("Performing query for WFC3 images.")
    wfc3_query_table = Observations.query_criteria(coordinates=sky_coord,
                                                   radius=radius,
                                                   dataproduct_type="IMAGE",
                                                   instrument_name="WFC3*",
                                                   project="HAP", 
                                                   obs_collection="HST") 

    query_table = vstack([acs_query_table, wfc3_query_table])
    del acs_query_table
    del wfc3_query_table

    # Catch the case where no files are found which satisfied the Query
    if not query_table:
        log.warning("Query for objects within {} of {} returned NO RESULTS!".format(radius, (str_ra, str_dec)))
        return query_table

    # Compute the limits of the cutout region
    deg_cutout_size = cutout_size.to(u.deg)
    ra_min = sky_coord.ra.degree - deg_cutout_size.value[0]
    ra_max = sky_coord.ra.degree + deg_cutout_size.value[0]
    dec_min = sky_coord.dec.degree - deg_cutout_size.value[1]
    dec_max = sky_coord.dec.degree + deg_cutout_size.value[1]
    str_ra = "{:.4f}".format(sky_coord.ra.degree)
    str_dec = "{:.4f}".format(sky_coord.dec.degree)

    # Filter the output as necessary to include only MVM filenames (MVM prefix: hst_skycell).
    # Also, filter out images which are not actually in the requested cutout region as the
    # archive search had to be done using a radius.
    good_rows = []
    updated_query_table = None
    for old_row in query_table:
        if old_row["obs_id"].startswith("hst_skycell"):
            if old_row["s_ra"] >= ra_min and old_row["s_ra"] <= ra_max and \
               old_row["s_dec"] >= dec_min and old_row["s_dec"] <= dec_max:
                good_rows.append(old_row)
    
    # Catch the case where no files are found which satisfy the clean up criteria
    if len(good_rows) == 0:
        log.warning("Query for objects within cutout {} of {} returned NO RESULTS!".format(cutout_size, (str_ra, str_dec)))
        return updated_query_table
    
    # Make the cleaned up table
    updated_query_table = Table(rows=good_rows, names=query_table.colnames)
    del query_table

    # Get the data product list associated with the elements of the table
    log.info("Get the product list for all entries in the query table.")
    dp_table = Observations.get_product_list(updated_query_table)
    del updated_query_table

    # Filter on MVM drizzled products only
    suffix = ["DRZ", "DRC"]
    log.info("Filter the product list table for only {} filenames.".format(suffix))
    filtered_dp_table = Observations.filter_products(dp_table,
                                                     productSubGroupDescription=suffix,
                                                     extension="fits")

    if not filtered_dp_table:
        log.warning("No MVM drizzle product datasets (DRZ/DRC) found within {} of {}.".format(radius, (str_ra, str_dec)))
        return filtered_dp_table
    del dp_table

    # Need to filter out any non-hst-skycell entries AGAIN which may have
    # crept back into the list via the get_product_list() function.
    good_rows = []
    output_table = None
    for old_row in filtered_dp_table:
        if old_row["obs_id"].startswith("hst_skycell"):
            good_rows.append(old_row)
    
    # Catch the case where no files are found which satisfy the criteria
    if len(good_rows) == 0:
        log.warning("After filtering datasets there are NO RESULTS within {} of {}!".format(radius, (str_ra, str_dec)))
        return output_table
    
    # Make the output table
    output_table = Table(rows=good_rows, names=filtered_dp_table.colnames)
    del filtered_dp_table

    # Finally, make sure the entries are unique
    final_table = None
    final_table = unique(output_table, keys="productFilename")
    del output_table

    # Write the table to a file.  This allows for further manipulation of
    # the information before a list of filenames is distilled from the table.
    # Output filename in the form: mvm_query-ra<###>d<####>-dec<n|s><##>d<####>_<radius>_cutout.ecsv
    # (e.g., mvm_query-ra84d9208-decs69d1483_71_cutout.ecsv), where radius has been computed from the
    # cutout dimensions.
    #
    # Get the whole number and fractional components of the RA and Dec
    ns = "s" if sky_coord.dec.degree < 0.0 else "n"
    ra_whole = int(sky_coord.ra.value)
    ra_frac  = str(sky_coord.ra.value).split(".")[1][0:4]
    dec_whole = abs(int(sky_coord.dec.value))
    dec_frac = str(sky_coord.dec.value).split(".")[1][0:4]
    log.info("coords2. {} {} {}".format(sky_coord.ra.value, sky_coord.dec.value, dec_frac))

    query_filename = "mvm_query-ra" + str(ra_whole) + "d" + ra_frac + "-dec" + ns + \
                     str(dec_whole) + "d" + dec_frac + "_{:.0f}".format(radius.value) + "_cutout.ecsv"

    log.info("Writing out the MVM product list table to {}.".format(query_filename))
    log.info("Number of entries in table: {}.".format(len(final_table)))
    final_table.write(query_filename, format="ascii.ecsv")

    return final_table


def mvm_retrieve_files(products, archive=False, clobber=False, log_level=logutil.logging.INFO):
    """
    This function retrieves specified files from the archive - unless the file is found
    to be locally resident on disk.  Upon completion, The function returns a list of 
    filenames available on disk. 

    Parameters
    ----------
    products : `astropy.table.Table <https://docs.astropy.org/en/stable/api/astropy.table.Table.html>`_ object
        A Table of products as returned by the mvm_id_filenames function. 

    archive : Boolean, optional
        Retain copies of the downloaded files in the astroquery created
        sub-directories? Default is "False".

    clobber : Boolean, optional
        Download and Overwrite existing files? Default is "False".

    log_level : int, optional
        The desired level of verbosity in the log statements displayed on the screen and written to the
        .log file. Default value is 20, or 'info'.

    Returns
    -------
    local_files : list
        List of filenames

    Note: Code here cribbed from retrieve_obsevation in astroquery_utils module.
    """

    # set logging level to user-specified level
    log.setLevel(log_level)

    # Determine if the files of interest are already on the local disk. If so,
    # remove the filename from the download list.
    all_images = []
    all_images = products['productFilename'].tolist()
    if not clobber:
        rows_to_remove = []
        for row_idx, row in enumerate(products):
            fname = row['productFilename']
            if os.path.isfile(fname):
                log.info(fname + " already exists. File download skipped.")
                rows_to_remove.append(row_idx)
        products.remove_rows(rows_to_remove)

    # Only download files as necessary
    if products:
        # Actual download of products
        log.info("Downloading files now...")
        manifest = Observations.download_products(products, mrp_only=False)
    else:
        log.info("There are no files to download as they are all resident on disk.")

    # Manifest has the following columns: "Local Path", "Status", "Message", and "URL"
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

    return local_files


def make_the_cut(input_files, sky_coord, cutout_size, output_dir=".", log_level=logutil.logging.INFO, verbose=False):
    """
    This function makes the actual cut in the input MVM drizzled filter- and exposure-level FITS
    files. As such it is a high-level interface for the
    `astrocut.cutouts.fits_cut <https://astrocut.readthedocs.io/en/latest/astrocut/index.html>`_ functionality.

    Parameters
    ----------
    input_files : list
        List of fits image filenames from which to create cutouts.

    sky_coord : str or `astropy.coordinates.SkyCoord <https://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html>`_ object
        The position around which to cutout. It may be specified as a string ("ra dec" in degrees)
        or as the appropriate `astropy.coordinates.SkyCoord
        <https://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html>`_ object.

    cutout_size : int, array-like, `astropy.units.Quantity <https://docs.astropy.org/en/stable/api/astropy.units.quantity.Quantity.html>`_
        The size of the cutout array. If ``cutout_size`` is a scalar number or a scalar
        `astropy.units.Quantity <https://docs.astropy.org/en/stable/api/astropy.units.quantity.Quantity.html>`_,
        then a square cutout of ``cutout_size`` will be created.
        If ``cutout_size`` has two elements, they should be in ``(ny, nx)`` order.  Scalar numbers
        in ``cutout_size`` are assumed to be in units of arcseconds.
        `astropy.units.Quantity <https://docs.astropy.org/en/stable/api/astropy.units.quantity.Quantity.html>`_ objects
        must be in angular units.

    output_dir : str
        Default value '.'. The directory where the cutout file(s) will be saved.

    log_level : int, optional
        The desired level of verbosity in the log statements displayed on the screen and written to the
        .log file. Default value is 20, or 'info'.

    verbose : bool
        Default False. If True, additional intermediate information is printed for the underlying 
        `spacetelescope.astrocut <https://astrocut.readthedocs.io/en/latest/>`_ utilities.

    Returns
    -------
    response : list
        Returns a list of all the output filenames.


    Note: For each input file designated for a cutout, there will be a corresponding output file.
        Since both the SCI and WHT extensions of the input files are actually cut, individual fits files
        will contain two image extensions, a SCI followed by the WHT.

        While the standard pipeline processing does not produce an MVM exposure-level drizzled
        product, it is possible for a user to turn on this capability in the pipeline while performing
        custom processing.  As such this routine will perform cutouts of the exposure-level drizzled files.

        Each filter-level output filename will be of the form:
          hst_cutout_skycell-p<pppp>-ra<##>d<####>-dec<n|s><##>d<####>_instrument_detector_filter[_platescale].fits
        Each exposure-level filename will be of the form:
          hst_cutout_skycell-p<pppp>-ra<##>d<####>-dec<n|s><##>d<####>_instrument_detector_filter[_platescale]-ipppssoo.fits

        where platescale has the value of "coarse" representing 0.12"/pixel for WFC3/IR, or there
        is no platescale value present which is the default and represents a "fine" platescale of 0.04"/pixel.

    """

    # set logging level to user-specified level
    log.setLevel(log_level)

    # Set the values for fits_cut that we are not allowing the user to modify
    CORRECT_WCS = False
    EXTENSION = [1, 2]  # SCI and WHT
    OUTPUT_PREFIX = "hst_cutout_skycell-"
    MEMORY_ONLY = True # This code will modify the output before it is written.
    SINGLE_OUTFILE = False

    # Making sure we have an array of images
    if type(input_files) == str:
        input_files = [input_files]

    # Retain the fully qualified pathname for each file
    input_paths = {}
    for files in input_files:
        input_paths[os.path.basename(files)] = os.path.dirname(files)

    # If the cutout_size is not an astropy.units.Quantity object, the scalar(s)
    # are assumed to be arcseconds.  The variable must be cast as a Quantity.
    if not isinstance(cutout_size, Quantity):
        cutout_size *= u.arcsec

    if not isinstance(sky_coord, SkyCoord):
        sky_coord = SkyCoord(sky_coord, unit="deg")

    # Call the cutout workhorse
    # SINGLE_OUTFILE = FALSE ==> MULTIPLE FILES: For each file cutout, there is an HDUList
    # comprised of a PHDU and one or more EHDUs.  The out_HDUList is then a list of HDULists.
    #
    # Loop over the input list so if there is an exception with a file, the
    # remaining files can still be used to generate cutout images.
    tmp_HDUList = []
    out_HDUList = []
    for infile in input_files:
        try:
            tmp_HDUList = fits_cut(infile, sky_coord, cutout_size, correct_wcs=CORRECT_WCS,
                                   extension=EXTENSION, single_outfile=SINGLE_OUTFILE, cutout_prefix=OUTPUT_PREFIX,
                                   output_dir=output_dir, memory_only=MEMORY_ONLY, verbose=True)
   
            # Copy and append the first (and it turns out the only) entry/list in the list
            out_HDUList.append(copy.deepcopy(tmp_HDUList[0]))
        except Exception as x_cept:
            log.error("")
            log.error("Exception encountered during the cutout process: {}".format(x_cept))
            log.error("No cutout files were created for file: {}.".format(infile))

    # hst_cutout_skycell-p<pppp>-ra<##>d<####>-dec<n|s><##>d<####>_instrument_detector_filter[_platescale][-ipppssoo].fits
    # Get the whole number and fractional components of the RA and Dec
    ra_whole = int(sky_coord.ra.value)
    ra_frac  = str(sky_coord.ra.value).split(".")[1][0:4]
    dec_whole = abs(int(sky_coord.dec.value))
    dec_frac = str(sky_coord.dec.value).split(".")[1][0:4]
    ns = "s" if sky_coord.dec.degree < 0.0 else "n"

    filename_list = []
    for HDU in out_HDUList:
   
        # Update only the image extensions
        extlist = HDU[1:]
        
        # Update the EXTNAME for all of the EHDUs
        for index in range(len(extlist)):
            input_filename = extlist[index].header["ORIG_FLE"]
            tokens = input_filename.split("_")
            skycell = tokens[1].split("-")[1]
            instr = tokens[2]
            detector = tokens[3]
            filter = tokens[4]
            label_plus = tokens[5]
            old_extname= extlist[index].header["O_EXT_NM"].strip().upper()
            extlist[index].header["EXTNAME"] = old_extname

            # Determine if the file is WFC3/IR which has both a "fine" (default) and
            # "coarse" platescale.
            plate_scale = "_coarse" if label_plus.upper().find("COARSE") != -1 else ""

            # Since the multiple output cutout files can also be input to the CutoutsCombiner,
            # there is some additional keyword manipulation done in the header.
            #
            # SCI extensions are followed by WHT extensions - when the WHT extension
            # has been updated, it is time to write out the file.
            if old_extname == "WHT":

                # Construct an MVM-style output filename with detector and filter
                output_filename = OUTPUT_PREFIX + skycell + "-ra" + str(ra_whole) + \
                                  "d" + ra_frac + "-dec" + ns + str(dec_whole) + "d" + \
                                  dec_frac + "_" + instr + "_" + detector + "_" + filter + plate_scale + ".fits"

                # Determine if the original file were a filter-level or exposure-level MVM product
                # ORIG_FLE filter-level: hst_skycell-p1253x05y09_acs_wfc_f658n_all_drc.fits
                # ORIG_FLE filter-level: hst_skycell-p0081x14y15_wfc3_ir_f128n_coarse-all_drz.fits
                # ORIG_FLE filter-level: hst_skycell-p0081x14y15_wfc3_ir_f128n_all_drz.fits (fine scale)
                # ORIG_FLE exposure-level: hst_skycell-p0081x14y15_wfc3_ir_f128n_coarse-all-ibp505mf_drz.fits
                # NOTE: Be careful of the WFC3/IR filenames which can include "coarse".
                ef_discriminant = label_plus.split("-")[-1]
                if ef_discriminant.upper() != "ALL":
                    product_type="EXPOSURE"
                    output_filename = output_filename.replace(".fits", "-" + ef_discriminant + ".fits")
                else:
                    product_type="FILTER"

                # Examples of output cutout filenames:
                # hst_cutout_skycell-p0081x14y15-ra84d9207-decs69d8516_wfc3_uvis_f275w.fits
                # hst_cutout_skycell-p0081x14y15-ra84d9207-decs69d8516_acs_wfc_f814w-jbp505jg.fits
                # hst_cutout_skycell-p0081x14y15-ra84d9207-decs69d8516_wfc3_ir_f128n_coarse.fits
                # hst_cutout_skycell-p0081x14y15-ra84d9207-decs69d8516_wfc3_ir_f128n_coarse-ibp505mf.fits
                cutout_path = os.path.join(output_dir, output_filename)

                log.info("Cutout FITS filename: {}".format(cutout_path))

                # Retain some keywords written in the PHDU of the cutout file
                # by the astrocut software
                ra_obj = HDU[0].header["RA_OBJ"]
                dec_obj = HDU[0].header["DEC_OBJ"]

                # Replace the minimal primary header written by the astrocut
                # software with the primary header from the corresponding input file,
                # so we can retain a lot of information from the observation
                HDU[0].header = fits.getheader(os.path.join(input_paths[input_filename], input_filename))

                # Put the new RA/DEC_OBJ keywords back
                HDU[0].header["RA_OBJ"] = (ra_obj, "[deg] right ascension")
                HDU[0].header["DEC_OBJ"] = (dec_obj, "[deg] declination")

                # Update PHDU FILENAME keyword with the new filename
                HDU[0].header['FILENAME'] = output_filename

                # Insert the new keyword, ORIG_FLE, in the PHDU which is the
                # *input* filename.  This keyword is also in the EHDUs.
                HDU[0].header["ORIG_FLE"] = input_filename

                output_HDUs = fits.HDUList(HDU)
                output_HDUs.writeto(cutout_path, overwrite=True)

                filename_list.append(output_filename)

    # Clean up any files left by `˜astrocut.cutouts.fits_cut`
    try:
        cruft_filenames = glob.glob(output_dir + "/hst_skycell*_astrocut.fits")
        if cruft_filenames:
            for cf in cruft_filenames:
               os.remove(cf)
    except Exception as x_cept:
        log.warning("")
        log.warning("Exception encountered: {}.".format(x_cept))
        log.warning("The following residual files could not be deleted from disk. " \
                    "Please delete these files to avoid confusion at your earliest convenience:")
        pprint(cruft_filenames)

    return filename_list


def mvm_combine(cutout_files, output_dir=".", log_level=logutil.logging.INFO):
    """
    This function combines multiple MVM skycell cutout images from the same detector/filter combination
    to create a single view of the requested data.  All of the functions in this module are designed to
    work in conjunction with one another, so the cutout images should be on the user's local disk.  This
    task is a high-level wrapper for the
    `astrocut.cutout_processing.combine
    <https://astrocut.readthedocs.io/en/latest/astrocut/index.html#combining-cutouts>`_ functionality.

    Specifically, this routine will combine filter-level cutouts from multiple skycells, all sharing 
    the same detector and filter.  This routine will also combine exposure-level cutouts from
    multiple skycells, all sharing the same detector, filter, and ipppssoo.  Images which do not
    share a detector and filter with any other image will be ignored. Individual exposures from
    a single skycell will also be ignored.

    Parameters
    ----------
    cutout_files : list
        List of fits image cutout filenames where the cutouts are presumed to have been created
        with `~drizzlepac.haputils.hapcut_utils.make_the_cut`.

    output_dir : str
        Default value '.' - The directory where the output combined files will be saved.

    log_level : int, optional
        The desired level of verbosity in the log statements displayed on the screen and written to the
        .log file. Default value is 20, or 'info'.

    """
    # def mvm_combine(cutout_files, img_combiner=None, output_dir=".", log_level=logutil.logging.INFO):
    img_combiner = None

    # set logging level to user-specified level
    log.setLevel(log_level)

    # Make sure the cutout_files are really a list of MULTIPLE filenames
    if type(cutout_files) == str or type(cutout_files) == list and len(cutout_files) < 2:
        log.error("The 'mvm_combine' function requires a list of MULTIPLE cutout filenames where" \
                  " the files were generated by 'make_the_cut'.")

    # Sort the cutout filenames by detector (primary) and filter (secondary)
    cutout_files.sort(key = lambda x: (x.split("_")[3], x.split("_")[4]))

    # Report the cutout files submitted for the combination process
    log.info("Input cutout files:")
    for cf in cutout_files:
        log.info("File: {}".format(cf))

    # Examples of input cutout filenames
    # Filter-level
    # hst_cutout_skycell-p0081x14y15-ra84d9207-decs69d8516_wfc3_uvis_f275w.fits
    # hst_cutout_skycell-p0081x14y15-ra84d9207-decs69d8516_wfc3_ir_f128n_coarse.fits
    # Exposure-level
    # hst_cutout_skycell-p0081x14y15-ra84d9207-decs69d8516_acs_wfc_f814w-jbp505jg.fits
    # hst_cutout_skycell-p0081x14y15-ra84d9207-decs69d8516_wfc3_ir_f128n_coarse-ibp505mf.fits
    #
    # Combined filter-level files will be generated for each detector/filter combination
    # Combined exposure-level files will be generated for each detector/filter combination
    #     where the ipppssoo is the same
    #
    # Walk the sorted input list and create filter-level and exposure-level dictionaries
    filter_dict = {}
    exposure_dict = {}
    for cfile in cutout_files:

        # Since the filename could be modified, open the file and read the FILENAME keyword
        hdu0 = fits.getheader(cfile, ext=0)
        cf = hdu0["FILENAME"].replace(".fits", "")

        # Parse to get the important information
        tokens = cf.split("_")
        instr = tokens[3]
        detector = tokens[4]
        filter = tokens[5].split("-")[0]
        str_tmp = tokens[-1].split("-")
        ipppssoo = ""
        if len(str_tmp) > 1:
            ipppssoo = str_tmp[1]

        # Based upon type of input file, filter-level or exposure-level, populate
        # the appropriate dictionary
        instr_det_filt_ippp = ""
        instr_det_filt = ""
        if ipppssoo:
            instr_det_filt_ippp = instr + "_" + detector + "_" + filter + "_" + ipppssoo
            if instr_det_filt_ippp not in exposure_dict:
                exposure_dict[instr_det_filt_ippp] = [cfile]
            else:
                exposure_dict[instr_det_filt_ippp].append(cfile)
        else:
            instr_det_filt = instr + "_" + detector + "_" + filter
            if instr_det_filt not in filter_dict:
                filter_dict[instr_det_filt] = [cfile]
            else:
                filter_dict[instr_det_filt].append(cfile)

    # FILTER-LEVEL COMBINATION
    # For each detector/filter, generate the output filename and perform the combine
    log.info("")
    log.info("=== Combining filter-level files ===")
    __combine_cutouts(filter_dict, type="FILTER", img_combiner=img_combiner, output_dir=output_dir, log_level=log_level)

    # EXPOSURE-LEVEL COMBINATION
    log.info("")
    log.info("=== Combining exposure-level files ===")
    __combine_cutouts(exposure_dict, type="EXPOSURE", img_combiner=img_combiner, output_dir=output_dir, log_level=log_level)

    log.info("Cutout combination is done.")


def __combine_cutouts(input_dict, type="FILTER", img_combiner=None, output_dir=".", log_level=logutil.logging.INFO):
    """
    This private function performs the actual combine of the multiple MVM skycell cutout images.

    Parameters
    ----------
    input_dict : dictionary 
        A dictionary where the key is the instr_detector_filter or instr_detector_filter_ipppssoo string and
        the corresponding value is a list of filenames corresponding to the key.

    type : string
        A string to indicate whether the input_dict variable is for a filter-level or exposure-level
        dictionary

    img_combiner : func
        The function to be used to combine the images

    output_dir : str
        Default value '.'. The directory to save the cutout file(s) to.

    log_level : int, optional
        The desired level of verbosity in the log statements displayed on the screen and written to the
        .log file. Default value is 20, or 'info'.

    """

    # set logging level to user-specified level
    log.setLevel(log_level)

    # Output prefix
    OUTPUT_PREFIX = "hst_combined_skycells-"

    for key, file_list in input_dict.items():

        # If there are multiple files to combine, then do it
        if len(file_list) > 1:

            # Construct the combined filename based on the first file in the list
            # Example: hst_combined_skycells-ra84d9207-decs69d8516_wfc3_uvis_f275w.fits
            filename = fits.getheader(file_list[0], ext=0)['FILENAME']
            fname = filename.replace(".fits", "")
            sky_tokens = fname.split("_")[2].split("-")
            skycell = sky_tokens[1][1:5]
            ra = sky_tokens[2]
            dec = sky_tokens[3]

            instr = key.split("_")[0]
            detector = key.split("_")[1]
            filter = key.split("_")[2]
            if type.upper() == "EXPOSURE":
                exposure = key.split("_")[3]
                output_filename = os.path.join(output_dir, OUTPUT_PREFIX + ra + "-" + dec + "_" + \
                                               instr + "_" + detector + "_" + filter + "_" + exposure + ".fits")
            else:
                output_filename = os.path.join(output_dir, OUTPUT_PREFIX + ra + "-" + dec + "_" + \
                                               instr + "_" + detector + "_" + filter + ".fits")
    

            # Combine the SCI and then the WHT extensions in the specified files
            log.info("Combining the SCI and then the WHT extensions of the input cutout files.")
            try:
                combined_cutout = astrocut.CutoutsCombiner(file_list, img_combiner=img_combiner).combine(output_file=output_filename, \
                                                           memory_only=True)
            except Exception as x_cept:
                log.warning("The cutout combine was not successful for files, {}, due to {}.".format(file_list, x_cept))
                log.warning("Processing continuuing on next possible set of data.")
                continue

            log.info("The combined output filename is {}.\n".format(output_filename))

            # Add the FILENAME keyword to the PHDU of the in-memory output
            if output_filename.startswith("./"):
                output_filename = output_filename.replace("./", "")
            combined_cutout[0].header["FILENAME"] = output_filename
 
            # Update the EXTNAMEs of the EHDUs
            combined_cutout[1].header["EXTNAME"] = "SCI"
            combined_cutout[2].header["EXTNAME"] = "WHT"

            # Write out the file
            combined_cutout.writeto(output_filename, overwrite=True)

        # Only a single file 
        else:
            log.warning("There is only one file for this detector/filter[/ipppssoo] combination, so there" \
                        " is nothing to combine.")
            log.warning("File {} will be ignored for combination purposes.\n".format(file_list))

