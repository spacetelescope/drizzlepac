# Pinning Versions of Dependencies 12/17/24

## Context

To ensure that our HAP pipeline products are consistent with each build, and to ensure that the drizzlepac code is working, sometimes we have to force the code to only install a certain version of a dependency package. 

## Decision

We aim to avoid pinning the versions of packages on the main branch (where we develop), but each build does pin all dependency version.

## Consequences 

This will avoid sudden changes to the build, while new bugs may pop up where we are developing on main. This leaves us in a better place to fix them and update the code. 


# Updating WCS in Headers 11/14/23

## Context

In updating WCS information there was a worry that information about datasets used in research would be lost.

## Decision

In order to ensure reproducibility and the ability to share consistent datasets "headerlets" were created. These hold the wcs inforamtion without any data.  

## Consequences

Headerlets can now be used to easily update and select WCS for datasets. Headerlets are files with the _hlet.fits suffix


# Updating APERTURE keyword through poller 11/14/23

## Context

Two types of poller files exist. Full poller files, that include all of the necessary information, and simple poller files that include only the filenames; additional information is taken from the header. We needed a way to update a header keyword ("APERTURE") and it was found that the easiest way to do this was to pass the keyword through the second column of the poller file.

## Decision

The second column of a poller file is now reserved (for WFPC2) for passing the aperture keyword through the code to update the header keyword. For other instruments, the code treats a two-column poller file as a simple (filename-only) poller file. 

## Consequences

Caution must be taken is using variations of the poller file while processing WFPC2 data.


# Primary Drizzlepac Module run functions 02/29/24

## Context

The primary modules for drizzlepac have main functions named "run". This comes from the TEAL interface's need to identify modules. There are also a number of cases when the primary function of a module is unclear (e.g. ablot.py)

## Decision

The primary function in the drizzlepac modules will be renamed to 'main()'. Any user interface will be renamed to 'user_main()'. Exceptions are commonly used standalone tools like Tweakback, Tweakreg, and Astrodrizzle. 

## Consequences

Calls of just "main" may be confusing if the function is loaded from any module. Calls to "main" should be left in their associated namespaces e.g. createMedian.main() not main(). 
=======
# The Use of the TEAL Interface 11/14/23

## Context

The code can be run using the interactive GUI TEAL, however, it is hard to maintain. TEAL is also used in the code to load the json parameter files. 

## Decision

In order to make the code more easily maintainable, we will no longer support the use of the GUI for running tasks, however, TEAL will still be used in the background to load the parameter files as there is no current alternative for parsing that data. 

## Consequences

TEAL will need to be included in drizzlepac until a replacement for parsing the json files can be found. 
