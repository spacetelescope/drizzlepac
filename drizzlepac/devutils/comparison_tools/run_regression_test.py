#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4 ai :

"""
This script automates the execution of regression_testing/compare_images.py and
regression_testing/compare_sourcelists.py.

* compare_images compares all images common to both user-specified input paths.
* compare_sourcelists compares all sourcelists common to both user-specified input paths' logs/ subdirectory.

Path
----
HLApipeline/regression_testing/run_regression_test.py

Dependencies
------------
* drizzlepac/drizzlepac/devutils/comparison_tools/compare_images.py
* drizzlepac/drizzlepac/devutils/comparison_tools/compare_sourcelists.py

Inputs
------
* Required inputs
    1: *pathNames*
        * A space-separated pair of paths where images to compare can be found.

* Optional general inputs:
    1: -v *verbose*
        * Print verbose output to screen?
        * Input choices: "True" or "False"
        * Default value: "True"

* Optional inputs for compare_images:
    1: -c *comparisonType*
        * Compare just image pixel data, just image header data, or both image pixel and image header data?
        * Input choices: "pixel data", "header data", or "both"
        * Default value: "both"

    2: -d *drzOnly*
        * Restrict comparisons to only drz.fits files?
        * Input choices: "True" or "False"
        * Default value: "False"

* Optional inputs for compare_sourcelists:
    1: -m *diffMode*
        * How should the comp-ref difference be calculated? "absolute" is simply the straight comp-ref difference. "peman" is the mean percent difference ((C-R)/avg(R)) x 100. "pdynamic" is the dynamic percent difference ((C-R)/R) x 100
        * Input choices: "absolute", "pmean" or "pdynamic"
        * Default value: "pmean"

    2: -p *plotGen*
        * Generate plots?
        * Input choices: "screen", "file", or "none'
        * Default value: "none'

Example
-------

Inspect common drz.fits images found in both 'foo/' and 'bar/' subdirectories; Inspect file structure and pixel data
but not header information with verbose output turned off; Use absolute differences when comparing sourcelists, and
generate plots::

    $HLApipeline/regression_testing/compare_images.py foo bar -c "pixel data" -d True -v False -p True -m absolute

Classes and Functions
---------------------
"""
import argparse
from drizzlepac.devutils.comparison_tools import compare_images
from drizzlepac.devutils.comparison_tools import compare_sourcelists
import glob


def convert_string_tf_to_logical_tf(in_tf):
    """
    simple function to convert string 'True' or 'False' value to a Boolean
    *True* or Boolean *False* value.

    :param in_tf: 'True'/'False' value
    :type in_tf: string
    :return: **in_tf** value as a Boolean.
    """
    in_tf = in_tf.lower()
    if in_tf == "true": rv = True
    if in_tf == "false": rv = False
    return (rv)
#-----------------------------------------------------------------------------------------------------------------------
def find_matched_sourcelists(ref_path,comp_path):
    """
    Generate a list of sourcelists found in both ref_path/logs and comp_path/logs.

    :param ref_path: input reference path
    :param comp_path: input comparision path
    :type ref_path: string
    :type comp_path: sting
    :return: list of sourcelists found in both ref_path/logs and comp_path/logs.
    """
    refList=[]
    compList=[]
    for item in glob.glob(ref_path+"/logs/*phot.txt"):
        refList.append(item.split("logs/")[-1])
    for item in glob.glob(ref_path+"*.ecsv"):
        refList.append(item)
    for item in glob.glob(comp_path+"/logs/*phot.txt"):
        compList.append(item.split("logs/")[-1])
    for item in glob.glob(comp_path+"*.ecsv"):
        compList.append(item)
    outList=list(set(refList).intersection(set(compList)))
    return(outList)
#=======================================================================================================================
if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description='Compare Images and Sourcelists')
    # required positional input arguments
    PARSER.add_argument('pathNames', nargs=2, help='A space-separated pair of paths where images to compare can be found.')

    # optional sourcelist comparision input arguments
    PARSER.add_argument('-m', '--diffMode', required=False, choices=["absolute", "pmean","pdynamic"], default="pmean",
                        help='How should the comp-ref difference be calculated? "absolute" is simply the stright comp-ref difference. "peman" is the mean percent difference ((C-R)/avg(R)) x 100. "pdynamic" is the dynamic percent difference ((C-R)/R) x 100. Default value is "pmean".')
    PARSER.add_argument('-p', '--plotGen', required=False, choices=["screen", "file", "none"], default="none",
                        help='Generate Plots? "screen" displays plots on-screen. "file" saves them to a .pdf file, and "none" skips all plot generation.')

    # optional image comparision input arguments
    PARSER.add_argument('-c', '--comparisonType', required=False, choices=["pixel data", "header data", "both"],
                        default="both",
                        help='Compare just image pixel data, just image header data, or both image pixel and image header data? Default value is "both".')
    PARSER.add_argument('-d', '--drzOnly', required=False, choices=["True", 'False'], default="False",
                        help='Restrict comparisons to only *drz.fits files (True/False)? Default value is "False".')

    # optional general input arguments
    PARSER.add_argument('-v', '--verbose', required=False, choices=["True", 'False'], default="True",
                        help='Print verbose output to screen (True/False)? Default value is "True".')
    ARGS = PARSER.parse_args()

    #Convert "True" / "False" input values from type string to type Boolean
    ARGS.drzOnly = convert_string_tf_to_logical_tf(ARGS.drzOnly)
    ARGS.verbose = convert_string_tf_to_logical_tf(ARGS.verbose)

    #Remove any trailing "/" characters from the input paths
    for ctr in range(0,2):
        if ARGS.pathNames[ctr].endswith("/"):
            ARGS.pathNames[ctr]=ARGS.pathNames[ctr][:-1]
    padding=68
    print(("{}>{}<{}".format("=" * padding, " COMPARE SOURCELISTS ", "=" * padding)))
    overallStatus="OK"
    slList = find_matched_sourcelists(ARGS.pathNames[0], ARGS.pathNames[1])
    slCtr =1

    for slName in slList:
        name_seperator = "{}> {}/{}: {} <{}".format("-" * 56, slCtr, len(slList), slName, "-" * 56)
        print(name_seperator)
        slPair = [ARGS.pathNames[0] + "/logs/" + slName, ARGS.pathNames[1] + "/logs/" + slName]
        if slName.endswith("daophot.txt"): imgname=slName.replace("daophot.txt","drz.fits")
        if slName.endswith("sexphot.txt"): imgname = slName.replace("sexphot.txt", "drz.fits")
        imgPair = [ARGS.pathNames[0] + "/" + imgname, ARGS.pathNames[1] + "/" + imgname]

        runStatus= compare_sourcelists.run(slPair, imgPair, ARGS.plotGen, ARGS.diffMode, ARGS.verbose)
        if runStatus != "OK": overallStatus = "FAILURE"
        slCtr+=1

    print(("{}>{}<{}".format("=" * padding, " COMPARE  IMAGES ", "=" * padding)))
    runStatus = compare_images.run(ARGS.pathNames, ARGS.comparisonType, ARGS.drzOnly, ARGS.verbose)
    if runStatus != "OK": overallStatus = "FAILURE"

    print()
    print()
    print()
    print(("OVERALL REGRESSION TEST RESULT: ",overallStatus))
