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

The primary function in the drizzlepac modules will be renamed to 'main()'. Any user interface will be renamed to 'user_main()'. 

## Consequences

Calls of just "main" may be confusing if the function is loaded from any module. Calls to "main" should be left in their associated namespaces e.g. createMedian.main() not main(). 