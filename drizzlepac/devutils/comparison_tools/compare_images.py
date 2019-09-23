#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4 ai :

"""
This script compares images with the same name in two user-specified paths and displays various measures of their
differences.

The following items are inspected:

    - File structure (Number of FITS extensions, extension names, image/table dimensions)
    - Header values
    - Pixel values in image extensions

Path
----
HLApipeline/regression_testing/compare_images.py

Dependencies
------------
None.

Inputs
------
* Required inputs
    1: *pathNames*
        * A space-separated pair of paths where images to compare can be found.

* Optional inputs:
    1: -c *comparisonType*
        * Compare just image pixel data, just image header data, or both image pixel and image header data?
        * Input choices: "pixel data", "header data", or "both"
        * Default value: "both"

    2: -d *drzOnly*
        * Restrict comparisons to only drz.fits files?
        * Input choices: "True" or "False"
        * Default value: "False"

    3: -v *verbose*
        * Print verbose output to screen?
        * Input choices: "True" or "False"
        * Default value: "True"

Example
-------

Inspect common drz.fits images found in both 'foo/' and 'bar/' subdirectories; Inspect file structure and pixel data but not header information with verbose output turned off::

    $HLApipeline/regression_testing/compare_images.py foo bar -c "pixel data" -d True -v False

Classes and Functions
---------------------
"""
import argparse
import os
import pdb
import sys
#import time
from astropy.io import fits
import numpy as np
#-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
def compareFileStructure(imgHDU1,imgHDU2,verbose):
    """
    Compare the structure of the two input fits files. Performs the following:

        - compare number of extensions
        - compare extension names
        - compare image extension sizes
        - compare table extension sizes

    :param imgHDU1: name of first fits file to compare
    :param imgHDU2: name of second fits file to compare
    :param verbose: Print verbose output to screen?
    :type imgHDU1: astropy.io.fits.hdu.hdulist.HDUList class
    :type imgHDU2: astropy.io.fits.hdu.hdulist.HDUList class
    :type verbose: Boolean
    :return: Logical 'True' if everything is ok, and logical 'False' if dissimilarities are found.
    """
    if verbose: print("> > > > > COMPARE FILE STRUCTURE < < < < <")

    returnValue = False

    if len(imgHDU1) == len(imgHDU2):
        returnValue = True
        for ext_num in range(0, len(imgHDU1)):
            img1Info = get_ext_info(imgHDU1[ext_num])
            img2Info = get_ext_info(imgHDU2[ext_num])

            for dictKey in list(img1Info.keys()):
                if img1Info[dictKey] != img2Info[dictKey]:
                    returnValue = False
            if  not returnValue:
                if verbose: print((ext_num, img1Info['name'],img2Info['name'],"     ",img1Info['type'],img2Info['type'],"     ",img1Info['size'],img2Info['size']))
    if returnValue:
        if verbose: print("OK")
    if verbose: print()
    return(returnValue)
# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
def compareHeaderValues(imgHDU1,imgHDU2,verbose):
    """
    Compare the headers of the two input fits files.

    :param imgHDU1: name of first fits file to compare
    :param imgHDU2: name of second fits file to compare
    :param verbose: Print verbose output to screen?
    :type imgHDU1: astropy.io.fits.hdu.hdulist.HDUList class
    :type imgHDU2: astropy.io.fits.hdu.hdulist.HDUList class
    :type verbose: Boolean
    :return: Logical 'True' if everything is ok, and logical 'False' if dissimilarities are found.
    """
    if verbose: print("> > > > > COMPARE HEADERS < < < < <")

    out_list = []
    exclude_list = ['DATE', 'IRAF-TLM']  # list of header keys to ignore from search.
    for ext_num in range(0, len(imgHDU1)):
        hdr_list1 = list(imgHDU1[ext_num].header.keys())
        hdr_list2 = list(imgHDU2[ext_num].header.keys())
        for hdr_title in hdr_list1:
            if ((hdr_title != 'HISTORY') and (exclude_list.count(hdr_title) == 0) and (len(hdr_title) != 0)):
                try:
                    hdr_value1 = imgHDU1[ext_num].header[hdr_title]
                except:
                    hdr_value1 = "--MISSING--"
                try:
                    hdr_value2 = imgHDU2[ext_num].header[hdr_title]
                except:
                    hdr_value2 = "--MISSING--"
                if hdr_value1 != hdr_value2:
                    out_list.append("[%s] %s: %s %s" % (str(ext_num), hdr_title, hdr_value1, hdr_value2))
    returnValue = True
    if len(out_list) > 0:
        returnValue = False
        if verbose:
            for item in out_list: print(item)
    if returnValue:
        if verbose: print("OK")
    if verbose: print()
    return (returnValue)
# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
def comparePixelValues(imgHDU1,imgHDU2,verbose):
    """
    Compare the pixel values in all extensions the two input fits files.

    :param imgHDU1: name of first fits file to compare
    :param imgHDU2: name of second fits file to compare
    :param verbose: Print verbose output to screen?
    :type imgHDU1: astropy.io.fits.hdu.hdulist.HDUList class
    :type imgHDU2: astropy.io.fits.hdu.hdulist.HDUList class
    :type verbose: Boolean
    :return: Logical 'True' if everything is ok, and logical 'False' if dissimilarities are found.
    """
    returnValue = True
    if verbose: print("> > > > > COMPARE PIXEL VALUES < < < < <")
    imgExtenList=[]
    #Determine which extensions are image extension
    for ext_num in range(0, len(imgHDU1)):
        try:
            extTest=imgHDU1[ext_num].shape
        except:
            extTest=()
        if extTest: imgExtenList.append(ext_num)
    #compute image statistics
    statMeasureList=['mean','stdev','median','max','min']
    for ext_num in imgExtenList:
        img1Stats = get_ext_stats(imgHDU1[ext_num].data)
        img2Stats = get_ext_stats(imgHDU2[ext_num].data)
        for item in statMeasureList:
            if img1Stats[item] != img2Stats[item]:
                pct_diff=((img2Stats[item]-img1Stats[item])/np.abs(img1Stats[item]))*100.0
                if verbose: print(("[{}] {}:{} {} {} ({}%)".format(ext_num, item," "*(6-len(item)), img1Stats[item], img2Stats[item],pct_diff)))
                returnValue = False
    if returnValue:
        if verbose: print("OK")
    if verbose: print()
    return (returnValue)
# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
def findFiles(rootPath,drzOnly):
    """
    Recursively search though specified path and any subdirectories therein for fits files (any filename ending in ".fits")

    :param rootPath: top level path to recursively search through
    :param drzOnly: Restrict comparisons to only drz.fits files?
    :type rootPath: string
    :type drzOnly: Boolean
    :return: list of fits files
    """
    if not os.path.exists(rootPath): sys.exit("Error! Directory %s not found. Exiting..."%(rootPath))
    file_list=[]
    #excludeList=['hlet.fits'] #filetypes to ignore
    excludeList = []
    if drzOnly:
        endsWithString="drz.fits"
    else:
        endsWithString = ".fits"
    for root, dirs, files in os.walk(rootPath):
        for name in files:
            fname = os.path.join(root, name)
            if fname.endswith(endsWithString):
                fname=fname.replace(rootPath, "")
                if fname.split("_")[-1] not in excludeList:
                    file_list.append(fname)
    return(file_list)
# -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
def get_ext_info(hdu):
    """
    Returns name, type and size of a specified HDU extension

    :param hdu: a single extension of a multiextension HDU item
    :type hdu: astropy.io.fits.hdu.image class
    :return: 3- element dictionary with text values for extension name, extension type, and extension size.
    """
    outDict={}
    outDict["name"]=hdu.name
    outDict["type"] = str(hdu).split(" object")[0].split(".")[-1]
    if ((outDict["type"] == 'ImageHDU') and (hdu.shape)):
        outDict["size"] = "{} x {}".format(hdu.shape[0],hdu.shape[1])
    elif "table" in outDict["type"].lower():
        nRows = len(hdu.data)
        nCols = len(hdu.data[0])
        outDict["size"] = "{}R x {}C".format(nRows,nCols)
    else:
        outDict["size"] = "()"
    return(outDict)
#-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
def get_ext_stats(imgRA):
    """
    Returns dictionary of mean, standard deviation, median, mode, max and min values for the specified image array.

    Ignores pixels with np.inf and np.nan values

    :param imgRA: array of pixel data to statistically examine
    :type imgRA: numpy.ndarray
    :return: five-element dictionary with values for mean, standard deviation, median, mode, max and min
    """
    statsDict={}
    if imgRA.shape:
        try:
            infIdx = np.where(imgRA == np.inf)
            imgRA[infIdx] = np.nan
        except:
            pass
        statsDict['mean'] = np.nanmean(imgRA)
        statsDict['stdev'] = np.nanstd(imgRA)
        statsDict['median'] = np.nanmedian(imgRA)
        statsDict['max'] = np.nanmax(imgRA)
        statsDict['min'] = np.nanmin(imgRA)
    else:
        statsDict['mean'] = "n/a"
        statsDict['stdev'] = "n/a"
        statsDict['median'] = "n/a"
        statsDict['max'] = "n/a"
        statsDict['min'] = "n/a"
    return(statsDict)
#-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
def make_file_list(pathNames,drzOnly):
    """
    Returns a list of fits files that exist in both paths specified in pathNames.

    :param pathNames: paths to search for fits flies
    :param drzOnly: Restrict comparisons to only drz.fits files?
    :type pathNames: two-element list
    :type drzOnly: Boolean
    :return: a list a fits files that exist in both paths specified in pathNames
    """
    file_list1 = findFiles(pathNames[0],drzOnly)
    file_list2 = findFiles(pathNames[1],drzOnly)
    if not file_list1: sys.exit("Error! No fits files found in directory %s. Exiting..." % (pathNames[0]))
    if not file_list2: sys.exit("Error! No fits files found in directory %s. Exiting..." % (pathNames[1]))
    unique_sorted_filelist = sorted(list(set(file_list1+file_list2)))
    outputList=[]
    for imgName in unique_sorted_filelist:
        imgName0=pathNames[0]+imgName
        imgName1=pathNames[1]+imgName
        #print imgName,imgName0,": ",os.path.exists(imgName0),imgName1,": ",os.path.exists(imgName1)
        if os.path.exists(imgName0) and os.path.exists(imgName1): outputList.append(imgName)
    if not outputList: sys.exit("Error! No fits files were found to exist in both %s and %s. Exiting..."%(pathNames[0],pathNames[1]))
    return outputList
#-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
def run(pathNames,compType,drzOnly,verbose):
    """
    Main calling subroutine.

    :param pathNames: A pair of paths where images to compare can be found.
    :param compType: Compare just image pixel data, just image header data, or both image pixel and image header data?
    :param drzOnly: Restrict comparisons to only drz.fits files?
    :param verbose: Print verbose output to screen?
    :type pathNames: list
    :type compType: string
    :type drzOnly: Boolean
    :type verbose: Boolean
    :return: "OK" if all tests were passed, or "FAILURE" if inconsistencies were found.
    """
    #1: get list of images common to both path1 and path2
    overallStatus = "OK"
    maxNameSepLength = 0
    for ctr in range(0,2):
        if not pathNames[ctr].endswith("/"): pathNames[ctr]=pathNames[ctr]+"/"
    print(pathNames)
    uniqueFileList=make_file_list(pathNames,drzOnly)
    check_status=["Skipped","Skipped","Skipped"]
    imgCtr=1
    numImgs=len(uniqueFileList)
    for imgName in uniqueFileList:
        imgName1 = pathNames[0] + imgName
        imgName2 = pathNames[1] + imgName
        imgHDU1 = fits.open(imgName1)
        imgHDU2 = fits.open(imgName2)

        name_seperator = "{}> {}/{}: {} <{}".format("-"*58,imgCtr,numImgs,imgName,"-"*58)
        if len(name_seperator) >maxNameSepLength: maxNameSepLength = len(name_seperator)
        print(name_seperator)
        check_status[0]=compareFileStructure(imgHDU1,imgHDU2,verbose)
        if compType in ["header data","both"]: check_status[1]=compareHeaderValues(imgHDU1,imgHDU2,verbose)
        if compType in ["pixel data", "both"]:check_status[2] = comparePixelValues(imgHDU1, imgHDU2,verbose)

        print()
        print(("File Structure... {}".format(check_status[0])))
        print(("Header Values.... {}".format(check_status[1])))
        print(("Pixel Data....... {}".format(check_status[2])))
        print()
        print()

        for item in check_status:
            if item == False:
                overallStatus = "FAILED"
        imgCtr +=1
    print(("="*maxNameSepLength))
    print(("Overall Status... {}".format(overallStatus)))
    return(overallStatus)
#=======================================================================================================================
if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description='Compare all images common to the user-specified pair of paths')
    # required positional input arguments
    PARSER.add_argument('pathNames', nargs=2, help='A space-separated pair of paths where images to compare can be found.')
    # optional input arguments
    PARSER.add_argument('-c', '--comparisonType', required=False, choices=["pixel data", "header data","both"], default="both",help='Compare just image pixel data, just image header data, or both image pixel and image header data?')
    PARSER.add_argument('-d', '--drzOnly', required=False, choices=["True",'False'], default="False",help='Restrict comparisons to only *drz.fits files (True/False)? Default value is "False')
    PARSER.add_argument('-v', '--verbose', required=False, choices=["True",'False'], default="True",help='Print verbose output to screen (True/False)? Default value is "True".')
    ARGS = PARSER.parse_args()

    if ARGS.drzOnly == "True": ARGS.drzOnly = True
    if ARGS.drzOnly == "False": ARGS.drzOnly = False

    if ARGS.verbose == "True": ARGS.verbose = True
    if ARGS.verbose == "False": ARGS.verbose = False

    #startTime = time.time()
    runStatus=run(ARGS.pathNames,ARGS.comparisonType,ARGS.drzOnly,ARGS.verbose)
    #print "Processing time = %s seconds" % (time.time() - startTime)


